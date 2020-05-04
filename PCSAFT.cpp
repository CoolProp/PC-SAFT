
#include <map>
#include <vector>
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

const static double k_Boltzmann = 1.3806622169047228e-23;
const static double PI = 3.141592654;
const static double N_AV = 6.022e23;

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
auto C1(Eta eta, Mbar mbar){
    return 1.0/(1.0
        + mbar*(8.0*eta-2.0*eta*eta)/pow(1.0-eta, 4)
        + (1.0-mbar)*(20.0*eta - 27.0*eta*eta + 12.0*pow(eta, 3) - 2.0*pow(eta, 4))/pow((1.0 - eta)*(2.0 - eta), 2));
}
/// Eqn. A.31
template <typename Eta, typename Mbar>
auto C2(Eta eta, Mbar mbar){
    return -pow(C1(eta, mbar), 2)*(
        mbar*(-4.0*eta*eta+20.0*eta+8.0)/pow(1.0-eta, 5)
        +(1.0-mbar)*(2.0*eta*eta*eta+12.0*eta*eta-48.0*eta+40.0)/pow((1.0-eta)*(2.0-eta), 3)
    );
}
/// Eqn. A.18
template<typename TYPE>
auto a_i(std::size_t i, TYPE mbar){
    static TYPE a_0[7] = { 0.9105631445, 0.6361281449, 2.6861347891, -26.547362491, 97.759208784, -159.59154087, 91.297774084 };
    static TYPE a_1[7] = { -0.3084016918, 0.1860531159, -2.5030047259, 21.419793629, -65.255885330, 83.318680481, -33.746922930 };
    static TYPE a_2[7] = { -0.0906148351, 0.4527842806, 0.5962700728, -1.7241829131, -4.1302112531, 13.776631870, -8.6728470368 };
    return a_0[i]+(mbar-1.0)/mbar*a_1[i] + (mbar-1.0)/mbar*(mbar-2.0)/mbar*a_2[i];
}
/// Eqn. A.19
template<typename TYPE>
auto b_i(std::size_t i, TYPE mbar){
    static TYPE b_0[7] = { 0.724094694, 2.2382791861, -4.0025849485, -21.003576815, 26.855641363, 206.55133841, -355.60235612 };
    static TYPE b_1[7] = { -0.5755498075, 0.6995095521, 3.8925673390, -17.215471648, 192.67226447, -161.82646165, -165.20769346 };
    static TYPE b_2[7] = { 0.0976883116, -0.2557574982, -9.1558561530, 20.642075974, -38.804430052, 93.626774077, -29.666905585 };
    return b_0[i]+(mbar-1.0)/mbar*b_1[i] + (mbar-1.0)/mbar*(mbar-2.0)/mbar*b_2[i];
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
/// Eqn. A.16
template <typename Eta>
auto I1(Eta eta, double mbar){
    Eta summer = 0.0;
    for (std::size_t i = 0; i < 7; ++i){
        summer += a_i(i, mbar)*pow(eta, i);
    }
    return summer;
}
/// Eqn. A.29
template <typename Eta>
auto d_etaI1_deta(Eta eta, double mbar){
    Eta summer = 0.0*eta;
    for (std::size_t j = 0; j < 7; ++j){
        summer += a_i(j, mbar)*(j+1.0)*pow(eta, static_cast<int>(j));
    }
    return summer;
}
/// Eqn. A.17
template <typename Eta>
auto I2(Eta eta, double mbar){
    Eta summer = 0.0*eta;
    for (std::size_t i = 0; i < 7; ++i){
        summer += b_i(i, mbar)*pow(eta, static_cast<int>(i));
    }
    return summer;
}
/// Eqn. A.30
template <typename Eta>
auto d_etaI2_deta(Eta eta, double mbar){
    Eta summer = 0.0*eta;
    for (std::size_t j = 0; j < 7; ++j){
        summer += b_i(j, mbar)*(j+1.0)*pow(eta, static_cast<int>(j));
    }
    return summer;
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

/// A class used to evaluate mixtures using PC-SAFT model
class PCSAFTMixture{
private:
    std::vector<double> m, ///< number of segments
                        sigma_Angstrom, ///< 
                        epsilon_over_k, ///< depth of pair potential divided by Boltzman constant
                        mole_fractions;
    std::vector<std::string> names;

    double mbar, ///< mean segment number in the system
        k_ij, ///< binary interaction parameter
        m2_epsilon_sigma3_bar, ///< Eqn. A. 12
        m2_epsilon2_sigma3_bar; ///< Eqn. A. 13
public:
    PCSAFTMixture(const std::vector<std::string> names, 
                  const std::vector<double> &mole_fractions) 
        :mole_fractions(mole_fractions), names(names)
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
        mbar = sumproduct(mole_fractions, m);
    };
    void print_info(){
        for (auto i = 0; i < m.size(); ++i){
            std::cout << i <<  " " << m[i] << " " << sigma_Angstrom[i] << " " << epsilon_over_k[i] << " " << mole_fractions[i] << std::endl;
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
    template<typename RhoType, typename TTYPE>
    auto calc_Z(RhoType rhomolar, TTYPE T){

        using std::pow;

        std::size_t N = m.size();
        m2_epsilon_sigma3_bar = 0;
        m2_epsilon2_sigma3_bar = 0;
        std::vector<decltype(T)> d(N);
        for (std::size_t i = 0; i < N; ++i){
            d[i] = sigma_Angstrom[i]*(1.0-0.12*exp(-3.0*epsilon_over_k[i]/T));
            for (std::size_t j = 0; j < N; ++j){
                // Eqn A.5
                auto sigma_ij = 0.5*sigma_Angstrom[i] + 0.5*sigma_Angstrom[j];
                auto eij_over_k = sqrt(epsilon_over_k[i]*epsilon_over_k[j])*(1.0-k_ij);
                m2_epsilon_sigma3_bar += mole_fractions[i]*mole_fractions[j]*m[i]*m[j]*eij_over_k/T*pow(sigma_ij, 3);
                m2_epsilon2_sigma3_bar += mole_fractions[i]*mole_fractions[j]*m[i]*m[j]*pow(eij_over_k/T, 2)*pow(sigma_ij, 3);
            }   
        }

        /// Convert from molar density to number density in molecules/Angstrom^3
        auto rho_A3 = rhomolar*N_AV*1e-30; //[molecules (not moles)/A^3]
        auto eta = (PI/6.0)*rho_A3*sumproduct(mole_fractions,m,powvec(d,3));

        /// Evaluate the components of zeta
        std::vector<RhoType> zeta;
        for (std::size_t n = 0; n < 4; ++n){
            // Eqn A.8
            zeta.push_back(
                (PI/6.0)*rho_A3*sumproduct(mole_fractions, m, powvec(d, static_cast<double>(n)))
            );
        }

        RhoType summer = 0.0*rhomolar;
        for (std::size_t i = 0; i < N; ++i){
            summer -= mole_fractions[i]*(m[i]-1.0)/gij_HS(zeta, d, i, i)*rho_A3_dgij_HS_drhoA3(zeta, d, i, i);
        }
        auto Z_hs_ = Z_hs(zeta);
        auto Z_hc = mbar*Z_hs_ + summer;
        if constexpr (std::is_same<RhoType, ChebTools::ChebyshevExpansion>::value){
            std::cout << (1/zeta[0]).y_Clenshaw(3000.0) << std::endl;
            std::cout << "zeta[0]: " << zeta[0].coef() << std::endl;
        }
        else{
            std::cout << (1/zeta[0]) << std::endl;
        }
        
        auto Z_disp = -2*PI*rho_A3*d_etaI1_deta(eta, mbar)*m2_epsilon_sigma3_bar
                      -PI*rho_A3*mbar*(C1(eta,mbar)*d_etaI2_deta(eta, mbar) + C2(eta,mbar)*eta*I2(eta, mbar))*m2_epsilon2_sigma3_bar;
        auto Z = 1.0 + Z_hc + Z_disp; //[-]
        return Z;
    }
    template<typename Rho, typename TTYPE>
    auto calc_p(Rho rhomolar, TTYPE T){
        /// Convert from molar density to number density of molecules in molecules/Angstroms^3
        auto rho_A3 = rhomolar*N_AV*1e-30; //[molecules (not moles)/A^3]
        auto p = calc_Z(rhomolar, T)*k_Boltzmann*T*rho_A3*1e30; //[Pa]
        return p;
    }
};

template <typename TYPE>
void do_calc(){
    std::vector<std::string> names = {"Methane","Ethane"};
    for (double z0 = 0.1; z0 < 1.0; z0 += 101){
        std::vector<double> z = {z0, 1.0-z0};
        PCSAFTMixture mix(names, z);
        //mix.print_info();
        double T = 200;
        TYPE rhomolar = 3000.0;
        double h = 1e-100;
        if constexpr (std::is_same<TYPE, std::complex<double>>::value){
            rhomolar += TYPE(0.0, 1e-100);
        }
        auto val = mix.calc_p(rhomolar, T);
        std::cout << z0 <<  " " << val << std::endl;
        
        if constexpr (std::is_same<TYPE, std::complex<double>>::value){
            std::cout << "dZdrho: " << std::imag(val)/h << std::endl;
        }

        if constexpr (std::is_same<TYPE, double>::value){
            auto rhomolar_ = ChebTools::ChebyshevExpansion::factory(11, [](double x){return x;}, 0.001, 10000);
            // std::cout << "||" << rhomolar_.coef() << "||" << std::endl;
            auto pcheb = mix.calc_p(rhomolar_, T);
            auto diff = (pcheb-val);
            //std::cout << diff.y_Clenshaw(rhomolar) << std::endl;
            // std::cout << diff.get_node_function_values() << std::endl;
            // std::cout << diff.coef().size() << std::endl;
            // std::cout << diff.coef() << std::endl;
            // std::cout << diff.companion_matrix(diff.coef());
            // auto roots = (pcheb-p).real_roots(true);
            // for (auto root: roots){
            //     std::cout << root << std::endl;
            // }
        }
    }
}
int main(){
    do_calc<double>();
    //do_calc<std::complex<double> >();
}
