#include <map>
#include <vector>
#include <chrono>
#include <string>
#include <cmath>
#include "math.h"
#include <complex>
#include <iostream>

#include "ChebTools/ChebTools.h"

/// Raise a CheyshevExpansion to a power
template <typename NType>
auto pow(const ChebTools::ChebyshevExpansion &ce, NType n){
    std::function<Eigen::ArrayXd(const Eigen::ArrayXd &)> f = [n](const Eigen::ArrayXd& y) {return y.pow(n); };
    return ce.apply(f);
}

/// Raise an EigenArray to a power
auto pow(const Eigen::ArrayXd& v, std::size_t n) {
    return v.pow(static_cast<int>(n));
}

const static double k_Boltzmann = 1.380649e-23;
const static double N_AV = 8.314462618/k_Boltzmann; 
const static double PI = 3.141592654;

/// Coefficients for one fluid
struct SAFTCoeffs{
    std::string name; ///< Name of fluid
    double m, ///< number of segments
        sigma_Angstrom, ///< [A] segment diameter
        epsilon_over_k; ///< [K] depth of pair potential divided by Boltzman constant
    std::string BibTeXKey; ///< The BibTeXKey for the reference for these coefficients
};

/// Manager class for PCSAFT coefficients
class PCSAFTLibrary{
    std::map<std::string, SAFTCoeffs> coeffs;
public:
    PCSAFTLibrary(){
        insert_normal_fluid("Methane", 1.0000, 3.7039, 150.03, "Gross-IECR-2001");
        insert_normal_fluid("Ethane", 1.6069, 3.5206, 191.42, "Gross-IECR-2001");
    }
    void insert_normal_fluid(const std::string &name, double m, const double sigma_Angstrom, const double epsilon_over_k, const std::string &BibTeXKey){
        SAFTCoeffs coeff;
        coeff.name = name;
        coeff.m = m;
        coeff.sigma_Angstrom = sigma_Angstrom;
        coeff.epsilon_over_k = epsilon_over_k;
        coeff.BibTeXKey = BibTeXKey;
        coeffs.insert(std::pair<std::string, SAFTCoeffs>(name, coeff));
    }
    const auto & get_normal_fluid(const std::string &name){
        auto it = coeffs.find(name);
        if (it != coeffs.end()){
            return it->second;
        }
        else{
        }
    }
};

/// Eqn. A.11
/// Erratum: should actually be 1/RHS of equation A.11 according to sample
/// FORTRAN code
template <typename Eta, typename Mbar>
auto C1(const Eta &eta, Mbar mbar){
    return 1.0/(1.0
        + mbar*(8.0*eta-2.0*eta*eta)/pow(1.0-eta, 4)
        + (1.0-mbar)*(20.0*eta - 27.0*eta*eta + 12.0*pow(eta, 3) - 2.0*pow(eta, 4))/pow((1.0 - eta)*(2.0 - eta), 2));
}
/// Eqn. A.31
template <typename Eta, typename Mbar>
auto C2(const Eta &eta, Mbar mbar){
    return -pow(C1(eta, mbar), 2)*(
        mbar*(-4.0*eta*eta+20.0*eta+8.0)/pow(1.0-eta, 5)
        +(1.0-mbar)*(2.0*eta*eta*eta+12.0*eta*eta-48.0*eta+40.0)/pow((1.0-eta)*(2.0-eta), 3)
    );
}
/// Eqn. A.18
template<typename TYPE>
Eigen::ArrayXd get_a(TYPE mbar){
    static Eigen::ArrayXd a_0 = (Eigen::ArrayXd(7) << 0.9105631445, 0.6361281449, 2.6861347891, -26.547362491, 97.759208784, -159.59154087, 91.297774084).finished();
    static Eigen::ArrayXd a_1 = (Eigen::ArrayXd(7) <<  -0.3084016918, 0.1860531159, -2.5030047259, 21.419793629, -65.255885330, 83.318680481, -33.746922930).finished();
    static Eigen::ArrayXd a_2 = (Eigen::ArrayXd(7) <<  -0.0906148351, 0.4527842806, 0.5962700728, -1.7241829131, -4.1302112531, 13.776631870, -8.6728470368).finished();
    return a_0+(mbar-1.0)/mbar*a_1 + (mbar-1.0)/mbar*(mbar-2.0)/mbar*a_2;
}
/// Eqn. A.19
template<typename TYPE>
Eigen::ArrayXd get_b(TYPE mbar){
    // See https://stackoverflow.com/a/35170514/1360263
    static Eigen::ArrayXd b_0 = (Eigen::ArrayXd(7) << 0.7240946941, 2.2382791861, -4.0025849485, -21.003576815, 26.855641363, 206.55133841, -355.60235612).finished() ;
    static Eigen::ArrayXd b_1 = (Eigen::ArrayXd(7) << -0.5755498075, 0.6995095521, 3.8925673390, -17.215471648, 192.67226447, -161.82646165, -165.20769346).finished() ;
    static Eigen::ArrayXd b_2 = (Eigen::ArrayXd(7) << 0.0976883116, -0.2557574982, -9.1558561530, 20.642075974, -38.804430052, 93.626774077, -29.666905585).finished() ;
    return b_0+(mbar-1.0)/mbar*b_1 + (mbar-1.0)/mbar*(mbar-2.0)/mbar*b_2;
}
/// Residual contribution from hard-sphere (Eqn. A.26)
template<typename VecType>
auto Z_hs(const VecType &zeta){
    auto Upsilon = 1.0-zeta[3];
    return zeta[3]/Upsilon
           + 3.0*zeta[1]*zeta[2]/(zeta[0]*pow(Upsilon, 2))
           + (3.0*pow(zeta[2], 3) - zeta[3]*pow(zeta[2], 3))/(zeta[0]*pow(Upsilon, 3));
}
/// Derivative term from Eqn. A.27
template<typename zVecType, typename dVecType>
auto rho_A3_dgij_HS_drhoA3(const zVecType &zeta, const dVecType &d, 
                           std::size_t i, std::size_t j){
    auto Upsilon = 1.0-zeta[3];
    return zeta[3]/pow(Upsilon, 2)
           + d[i]*d[j]/(d[i]+d[j])*(3.0*zeta[2]/pow(Upsilon, 2) + 6.0*zeta[2]*zeta[3]/pow(Upsilon, 3))
           + pow(d[i]*d[j]/(d[i]+d[j]), 2)*(4.0*pow(zeta[2], 2)/pow(Upsilon, 3) + 6.0*pow(zeta[2], 2)*zeta[3]/pow(Upsilon, 4));
}
/// Term from Eqn. A.7
template<typename zVecType, typename dVecType>
auto gij_HS(const zVecType &zeta, const dVecType &d, 
            std::size_t i, std::size_t j){
    auto Upsilon = 1.0-zeta[3];
    return 1.0/(Upsilon) + d[i]*d[j]/(d[i]+d[j])*3.0*zeta[2]/pow(Upsilon, 2)
        + pow(d[i]*d[j]/(d[i]+d[j]), 2)*2.0*pow(zeta[2], 2)/pow(Upsilon, 3);
}
/// Eqn. A.16, Eqn. A.29
template <typename Eta>
auto get_I1(const Eta &eta, double mbar){
    auto avec = get_a(mbar);
    Eta summer_I1 = 0.0, summer_etadI1deta = 0.0;
    for (std::size_t i = 0; i < 7; ++i){
        auto increment = avec(i)*pow(eta, static_cast<int>(i));
        summer_I1 += increment;
        summer_etadI1deta += increment*(i+1.0);
    }
    return std::make_tuple(summer_I1, summer_etadI1deta);
}
/// Eqn. A.17, Eqn. A.30
template <typename Eta>
auto get_I2(const Eta &eta, double mbar){
    auto bvec = get_b(mbar);
    Eta summer_I2 = 0.0*eta, summer_etadI2deta = 0.0*eta;
    for (std::size_t i = 0; i < 7; ++i){
        auto increment = bvec(i) * pow(eta, static_cast<int>(i));
        summer_I2 += increment;
        summer_etadI2deta += increment*(i+1.0);
    }
    return std::make_tuple(summer_I2, summer_etadI2deta);
}

PCSAFTLibrary library;


/**
Sum up two array-like objects that can each have different container types and value types
*/
template<typename VecType1, typename VecType2>
auto sumproduct(const VecType1 &v1, const VecType2 &v2){
    using ResultType = decltype(v1[0]*v2[0]);
    ResultType summer = 0.0;
    for (auto i = 0; i < v1.size(); ++i){
        summer += v1[i]*v2[i];
    }
    return summer;
}

/**
Sum up three array-like objects that can each have different container types and value types
*/
template<typename VecType1, typename VecType2, typename VecType3>
auto sumproduct(const VecType1 &v1, const VecType2 &v2, const VecType3 &v3){
    using ResultType = decltype(v1[0]*v2[0]*v3[0]);
    ResultType summer = 0.0;
    for (auto i = 0; i < v1.size(); ++i){
        summer += v1[i]*v2[i]*v3[i];
    }
    return summer;
}

/**
Sum up three array-like objects that can each have different container types and value types
*/
template<typename VecType1, typename NType>
auto powvec(const VecType1 &v1, NType n){
    auto o = v1;
    for (auto i = 0; i < v1.size(); ++i){
        o[i] = pow(v1[i], n);
    }
    return o;
}

/// Parameters for model evaluation
class SAFTCalc {
public:
    // Just temperature dependent things
    std::vector<double> d;

    // These things also have composition dependence
    double mbar,
           m2_epsilon_sigma3_bar, ///< Eq. A. 12
           m2_epsilon2_sigma3_bar; ///< Eq. A. 13
};

/// A class used to evaluate mixtures using PC-SAFT model
class PCSAFTMixture{
private:
    std::vector<double> m, ///< number of segments
        sigma_Angstrom, ///< 
        epsilon_over_k; ///< depth of pair potential divided by Boltzman constant
    std::vector<std::string> names;
    double k_ij; ///< binary interaction parameter
public:
    PCSAFTMixture(const std::vector<std::string> names) : names(names)
    {
        m.resize(names.size()); 
        sigma_Angstrom.resize(names.size()); 
        epsilon_over_k.resize(names.size());
        auto i = 0;
        for (auto name : names){
            const SAFTCoeffs & coeff = library.get_normal_fluid(name);
            m[i] = coeff.m;
            sigma_Angstrom[i] = coeff.sigma_Angstrom;
            epsilon_over_k[i] = coeff.epsilon_over_k;
            i++;
        }
        k_ij = 0;
    };
    void print_info(){
        std::cout << "i m sigma / A e/k / K \n  ++++++++++++++" << std::endl;
        for (auto i = 0; i < m.size(); ++i){
            std::cout << i <<  " " << m[i] << " " << sigma_Angstrom[i] << " " << epsilon_over_k[i] << std::endl;
        }
    }
    template<typename VecType1>
    auto max_rhoN(double T, const VecType1 &mole_fractions){
    	auto N = mole_fractions.size();
    	std::vector<decltype(T)> d(N);
        for (auto i = 0; i < N; ++i){
            d[i] = sigma_Angstrom[i]*(1.0-0.12*exp(-3.0*epsilon_over_k[i]/T));
        }
        return 6*0.74/PI/sumproduct(mole_fractions,m,powvec(d,3));
    }
    template<typename RhoType, typename TTYPE, typename VecType>
    auto calc_Z(const RhoType &rhomolar, TTYPE T, const VecType &mole_fractions) const{

        using std::pow;
        std::size_t N = m.size();

        SAFTCalc c;
        c.d.resize(N);
        c.m2_epsilon_sigma3_bar = 0;
        c.m2_epsilon2_sigma3_bar = 0;
        for (std::size_t i = 0; i < N; ++i){
            c.d[i] = sigma_Angstrom[i]*(1.0-0.12*exp(-3.0*epsilon_over_k[i]/T)); // [A]
            for (std::size_t j = 0; j < N; ++j){
                // Eq. A.5
                auto sigma_ij = 0.5*sigma_Angstrom[i] + 0.5*sigma_Angstrom[j];
                auto eij_over_k = sqrt(epsilon_over_k[i]*epsilon_over_k[j])*(1.0-k_ij);
                c.m2_epsilon_sigma3_bar += mole_fractions[i]*mole_fractions[j]*m[i]*m[j]*eij_over_k/T*pow(sigma_ij, 3);
                c.m2_epsilon2_sigma3_bar += mole_fractions[i]*mole_fractions[j]*m[i]*m[j]*pow(eij_over_k/T, 2)*pow(sigma_ij, 3);
            }
        }
        c.mbar = sumproduct(mole_fractions, m);

        /// Convert from molar density to number density in molecules/Angstrom^3
        auto rho_A3 = rhomolar*N_AV*1e-30; //[molecules (not moles)/A^3]
        
        /// Evaluate the components of zeta
        std::vector<RhoType> zeta(4);
        for (std::size_t n = 0; n < 4; ++n){
            // Eqn A.8
            zeta[n] = (PI/6.0)*rho_A3*sumproduct(mole_fractions, m, powvec(c.d, static_cast<double>(n)));
        }

        /// Packing fraction is the 4-th value in zeta
        auto eta = zeta[3];

        auto [I1, etadI1deta] = get_I1(eta, c.mbar);
        auto [I2, etadI2deta] = get_I2(eta, c.mbar);

        RhoType summer = 0.0*rhomolar;
        for (std::size_t i = 0; i < N; ++i){
            summer -= mole_fractions[i]*(m[i]-1.0)/gij_HS(zeta, c.d, i, i)*rho_A3_dgij_HS_drhoA3(zeta, c.d, i, i);
        }
        auto Z_hs_ = Z_hs(zeta);
        auto Z_hc = c.mbar*Z_hs_ + summer;
        auto Z_disp = -2*PI*rho_A3*etadI1deta*c.m2_epsilon_sigma3_bar
                      -PI*rho_A3*c.mbar*(C1(eta,c.mbar)*etadI2deta + C2(eta,c.mbar)*eta*I2)*c.m2_epsilon2_sigma3_bar;
        auto Z = 1.0 + Z_hc + Z_disp; //[-]
        return Z;
    }
    template<typename Rho, typename TTYPE, typename VecType>
    auto calc_p(const Rho &rhomolar, TTYPE T, const VecType& mole_fractions) const{
        /// Convert from molar density to number density
        auto p = calc_Z(rhomolar, T, mole_fractions)*k_Boltzmann*T*rhomolar*N_AV; //[Pa]
        return p;
    }
};

class ChebyshevSAFT {
private:
    PCSAFTMixture& mix;
    std::vector<ChebTools::ChebyshevExpansion> expansions, //< p(v)
                                               derivatives, //< \f$ dpdv|T \f$
                                               integrals; //< \f$\int p dv \f$
    double m_T = 0;
public:
    ChebyshevSAFT(PCSAFTMixture& mix) : mix(mix) {}; 
    auto T() { return m_T; }

    template<typename VecType>
    void make_pv_expansions(double T, const VecType& mole_fractions, double vmolarmin, double vmolarmax, int Ndegree, int Ndivisions) {
        auto factor = pow(vmolarmax / vmolarmin, 1.0/(Ndivisions-1.0));
        
        for (auto idivision = 0; idivision < Ndivisions-1; ++idivision) {
            auto xmin = vmolarmin * pow(factor, idivision), xmax = xmin * factor;

            expansions.emplace_back(ChebTools::ChebyshevExpansion::factory(
                Ndegree, 
                [this, T, &mole_fractions](const double& v) {return this->mix.calc_p(1 / v, T, mole_fractions); }, 
                xmin, xmax)
            );

            //auto func = [this, T, &mole_fractions](const Eigen::ArrayXd& v) {return this->mix.calc_p(1 / v, T, mole_fractions); };
            //auto x_nodes_n11 = ChebTools::get_CLnodes(Ndegree);
            //Eigen::ArrayXd x_k = ((xmax - xmin) * x_nodes_n11.array() + (xmax + xmin)) / 2.0;
            //Eigen::ArrayXd f = func(x_k);

            // Build expansion from nodal values in one fell swoop
            //expansions.emplace_back( ChebTools::ChebyshevExpansion::factoryf(Ndegree, f, xmin, xmax) );
        }
        for (auto& expansion : expansions) {
            derivatives.emplace_back(expansion.deriv(1));
            integrals.emplace_back(expansion.integrate(1));
        }
        m_T = T;
    }
    auto get_pv_roots(double p) {
        std::vector<double> out;
        std::size_t iexpansion = 0;
        for (auto &expan : expansions) {
            auto roots = (expan-p).real_roots(true);
            for (auto& root : roots) {
                if (derivatives[iexpansion].y_Clenshaw(root) < 0) {
                    out.push_back(root);
                }
            }
            iexpansion++;
        }
        return out;
    }
    auto get_extrema() {
        std::vector<double> roots, values;
        std::size_t iexpansion = 0;
        for (auto& deriv : derivatives) {
            auto rhos = deriv.real_roots(true);
            for (auto& rho : rhos) {
                auto val = expansions[iexpansion].y_Clenshaw(rho);
                roots.push_back(rho); values.push_back(val);
            }
            iexpansion++;
        }
        return std::make_tuple(roots, values);
    }
    auto get_Maxwell_condition(double v1, double v2) {
        double output = 0.0;
        
        // Find left and right intersections
        std::size_t iexpansion = 0, ileft=0, iright=0;
        for (auto& deriv : derivatives) {
            if (v1 >= deriv.xmin() && v1 <= deriv.xmax()) {
                ileft = iexpansion;
            }
            if (v2 >= deriv.xmin() && v2 <= deriv.xmax()) {
                iright = iexpansion;
            }
            iexpansion++;
        }
        // Add contribution from the volume root val to right edge of expansion
        auto& il = integrals[ileft];
        output += il.y_Clenshaw(il.xmax()) - il.y_Clenshaw(v1);
        // Add all integrals until...
        for (auto iex = ileft + 1; iex < iright; ++iex) {
            auto& i = integrals[iex];
            output += i.y_Clenshaw(i.xmax()) - i.y_Clenshaw(i.xmin());
        }
        // Add left edge to val
        auto& ir = integrals[iright];
        output += ir.y_Clenshaw(v2) - ir.y_Clenshaw(ir.xmin());
        return output;
    }
    auto get_psat() {
        auto [vroots, vals] = get_extrema();
        auto pmax = *std::max_element(vals.begin(), vals.end());
        if (vals.size() == 2) {
            auto presid = [this](double lnpval) -> double {
                auto pval = exp(lnpval);
                auto volumeroots = get_pv_roots(pval);
                auto rho0 = 1 / volumeroots[0], rho1 = 1 / volumeroots[1];
                auto pdeltaV = (volumeroots[1] - volumeroots[0]) * pval; // must be positive
                auto maxwell = get_Maxwell_condition(volumeroots[0], volumeroots[1]) - pdeltaV;
                return maxwell;
            };
            auto ce = ChebTools::ChebyshevExpansion::factory(20, presid, log(1e-6), log(pmax*0.999));
            auto lnproots = ce.real_roots(true);
            auto p = exp(lnproots[0]);
            return p;
        }
        else {
            throw std::invalid_argument("");
        }
    }
};

template <typename TYPE>
void do_calc(){
    std::vector<std::string> names = {"Methane"};
    for (double z0 = 0.1; z0 < 1.0; z0 += 101){
        std::vector<double> z = {1.0};
        PCSAFTMixture mix(names);
        mix.print_info();
        double T = 100;
        auto max_rhoN = mix.max_rhoN(T, z)*1e30; // [molecule/m^3]
        auto max_rhomolar = max_rhoN/N_AV;
        std::cout << "max rhomolar: " << max_rhomolar << std::endl; 
        TYPE rhomolar = 25591.48146033746;
        double h = 1e-100;
        if constexpr (std::is_same<TYPE, std::complex<double>>::value){
            rhomolar += TYPE(0.0, 1e-100);
        }
        auto val = mix.calc_p(rhomolar, T, z);
        std::cout << z0 <<  " " << val << std::endl;
        
        if constexpr (std::is_same<TYPE, std::complex<double>>::value){
            std::cout << "dZdrho: " << std::imag(val)/h << std::endl;
        }

        if constexpr (std::is_same<TYPE, double>::value){
            {
                auto startTime = std::chrono::high_resolution_clock::now();
                auto N = 1000;
                auto rhomolar = 3800.0;
                double pcheb = 0;
                for (auto i = 0; i < N; ++i) {
                    pcheb += mix.calc_p(rhomolar, T, z);
                }
                auto endTime = std::chrono::high_resolution_clock::now();
                double elap = std::chrono::duration<double>(endTime - startTime).count();
                std::cout << "elapsed to call p(rho): " << elap/N*1e6 << " us/call" << std::endl;
            }

            auto tic = std::chrono::high_resolution_clock::now();
            double rootsum = 0.0;
            auto N = 100;
            for (auto i = 0; i < N; ++i) {
                ChebyshevSAFT chebSAFT(mix);
                chebSAFT.make_pv_expansions(T, z, 1/(0.9*max_rhomolar), 1/1e-14, 8, 100);
                rootsum += chebSAFT.get_psat();
            }
            auto toc = std::chrono::high_resolution_clock::now();
            double elapp = std::chrono::duration<double>(toc - tic).count()/N;
            std::cout << "elapsed (total): " << elapp * 1e6 << " us; result: " <<  rootsum/N << std::endl;
        }
    }
}
int main(){
    do_calc<double>();
    //do_calc<std::complex<double> >();
}
