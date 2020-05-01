
#include <map>
#include <vector>
#include <string>
#include <cmath>
#include "math.h"
#include <complex>
#include <iostream>

#include "ChebTools/ChebTools.h"

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
    const SAFTCoeffs & get_normal_fluid(const std::string &name){
        std::map<std::string, SAFTCoeffs>::iterator it = coeffs.find(name);
        if (it != coeffs.end()){
            return it->second;
        }
        else{
        }
    }
};

PCSAFTLibrary library;

/// A class used to evaluate mixtures using PC-SAFT model
template <typename TYPE>
class PCSAFTMixture{
private:
    std::vector<TYPE> m, ///< number of segments
                        sigma_Angstrom, ///< 
                        epsilon_over_k, ///< depth of pair potential divided by Boltzman constant
                        d, ///< temperature-dependent segment diameter
                        zeta,
                        mole_fractions;
    std::vector<std::string> names;

    TYPE mbar, ///< mean segment number in the system
        k_ij, ///< binary interaction parameter
        m2_epsilon_sigma3_bar, ///< Eqn. A. 12
        m2_epsilon2_sigma3_bar; ///< Eqn. A. 13
public:
    PCSAFTMixture(const std::vector<std::string> names, const std::vector<TYPE> mole_fractions) : mole_fractions(mole_fractions), names(names)
    {
        m.clear(); sigma_Angstrom.clear(); epsilon_over_k.clear();
        for (std::vector<std::string>::const_iterator it = names.begin(); it != names.end(); ++it){
            const SAFTCoeffs & coeff = library.get_normal_fluid(*it);
            m.push_back(coeff.m);
            sigma_Angstrom.push_back(coeff.sigma_Angstrom);
            epsilon_over_k.push_back(coeff.epsilon_over_k);
        }
    };
    void init(){
        
        k_ij = 0;

        mbar = 0;
        std::size_t N = m.size();
        for (std::size_t i = 0; i < N; ++i){
            // Eqn A.5
            mbar += mole_fractions[i]*m[i];
        }
    }
    TYPE calc_Z(TYPE rhomolar, TYPE T){

        std::size_t N = m.size();
        m2_epsilon_sigma3_bar = 0;
        m2_epsilon2_sigma3_bar = 0;
        d.resize(N);
        for (std::size_t i = 0; i < N; ++i){
            d[i] = sigma_Angstrom[i]*(1.0-0.12*exp(-3.0*epsilon_over_k[i]/T));
            for (std::size_t j = 0; j < N; ++j){
                // Eqn A.5
                TYPE sigma_ij = 0.5*sigma_Angstrom[i] + 0.5*sigma_Angstrom[j];
                TYPE eij_over_k = sqrt(epsilon_over_k[i]*epsilon_over_k[j])*(1.0-k_ij);
                m2_epsilon_sigma3_bar += mole_fractions[i]*mole_fractions[j]*m[i]*m[j]*eij_over_k/T*pow(sigma_ij, 3);
                m2_epsilon2_sigma3_bar += mole_fractions[i]*mole_fractions[j]*m[i]*m[j]*pow(eij_over_k/T, 2)*pow(sigma_ij, 3);
            }   
        }

        /// Convert from molar density to total number density of molecules in mol/Angstroms^3
        TYPE rho_A3 = rhomolar*N_AV/10e30; //[molecules (not moles)/A^3]

        TYPE summer = 0;
        for (std::size_t i = 0; i < N; ++i){
            summer += mole_fractions[i]*m[i]*pow(d[i], 3);
        }
        TYPE eta = PI*rho_A3/6.0*summer;

        /// Cache evaluations of the components of zeta
        zeta.resize(4);
        for (std::size_t n = 0; n < 4; ++n){   
            TYPE summer = 0;
            for (std::size_t i = 0; i < N; ++i){
                // Eqn A.8
                summer += mole_fractions[i]*m[i]*pow(d[i], static_cast<int>(n));
            }
            zeta[n] = PI*rho_A3/6.0*summer;
        }

        summer = 0;
        for (std::size_t i = 0; i < N; ++i){
            summer -= mole_fractions[i]*(m[i]-1.0)/gij_HS(i, i)*rho_A3_dgij_HS_drhoA3(i, i);
        }
        TYPE Z_hc = mbar*Z_hs() + summer;
        
        TYPE Z_disp = -2*PI*rho_A3*d_etaI1_deta(eta, mbar)*m2_epsilon_sigma3_bar
                        -PI*rho_A3*mbar*(C1(eta)*d_etaI2_deta(eta, mbar) + C2(eta)*eta*I2(eta, mbar))*m2_epsilon2_sigma3_bar;
        TYPE Z = 1.0 + Z_hc + Z_disp; //[-]
        return Z;
    }
    TYPE calc_p(TYPE rhomolar, TYPE T){
        
        /// Convert from molar density to total number density of molecules in mol/Angstroms^3
        TYPE rho_A3 = rhomolar*N_AV/10e30; //[molecules (not moles)/A^3]
        
        TYPE p = calc_Z(rhomolar, T)*k_Boltzmann*T*rho_A3*1e30; //[Pa]
        return p;
    }
    /// Eqn. A.18
    TYPE a_i(std::size_t i, TYPE mbar){
        static TYPE a_0[7] = { 0.9105631445, 0.6361281449, 2.6861347891, -26.547362491, 97.759208784, -159.59154087, 91.297774084 };
        static TYPE a_1[7] = { -0.3084016918, 0.1860531159, -2.5030047259, 21.419793629, -65.255885330, 83.318680481, -33.746922930 };
        static TYPE a_2[7] = { -0.0906148351, 0.4527842806, 0.5962700728, -1.7241829131, -4.1302112531, 13.776631870, -8.6728470368 };
        return a_0[i]+(mbar-1.0)/mbar*a_1[i] + (mbar-1.0)/mbar*(mbar-2.0)/mbar*a_2[i];
    }
    /// Eqn. A.19
    TYPE b_i(std::size_t i, TYPE mbar){
        static TYPE b_0[7] = { 0.724094694, 2.2382791861, -4.0025849485, -21.003576815, 26.855641363, 206.55133841, -355.60235612 };
        static TYPE b_1[7] = { -0.5755498075, 0.6995095521, 3.8925673390, -17.215471648, 192.67226447, -161.82646165, -165.20769346 };
        static TYPE b_2[7] = { 0.0976883116, -0.2557574982, -9.1558561530, 20.642075974, -38.804430052, 93.626774077, -29.666905585 };
        return b_0[i]+(mbar-1.0)/mbar*b_1[i] + (mbar-1.0)/mbar*(mbar-2.0)/mbar*b_2[i];
    }
    /// Residual contribution from hard-sphere (Eqn. A.26)
    TYPE Z_hs(){
        return zeta[3]/(1.0-zeta[3]) 
               + 3.0*zeta[1]*zeta[2]/(zeta[0]*pow(1.0-zeta[3], 2)) 
               + (3.0*pow(zeta[2], 3) - zeta[3]*pow(zeta[2], 3))/(zeta[0]*pow(1.0-zeta[3], 3));
    }
    /// Derivative term from Eqn. A.27
    TYPE rho_A3_dgij_HS_drhoA3(std::size_t i, std::size_t j){
        return zeta[3]/pow(1.0-zeta[3], 2)
               + d[i]*d[j]/(d[i]+d[j])*(3.0*zeta[2]/pow(1.0-zeta[3], 2) + 6.0*zeta[2]*zeta[3]/pow(1.0-zeta[3], 3))
               + pow(d[i]*d[j]/(d[i]+d[j]), 2)*(4.0*pow(zeta[2], 2)/pow(1.0-zeta[3], 3) + 6.0*pow(zeta[2], 2)*zeta[3]/pow(1.0-zeta[3], 4));
    }
    /// Term from Eqn. A.7
    TYPE gij_HS(std::size_t i, std::size_t j){
        return 1.0/(1.0-zeta[3]) + d[i]*d[j]/(d[i]+d[j])*3.0*zeta[2]/pow(1.0-zeta[3], 2)
            + pow(d[i]*d[j]/(d[i]+d[j]), 2)*2.0*pow(zeta[2], 2)/pow(1.0-zeta[3], 3);
    }
    /// Eqn. A.16
    TYPE I1(TYPE eta, TYPE mbar){
        TYPE summer = 0;
        for (std::size_t i = 0; i < 7; ++i){
            summer += a_i(i, mbar)*pow(eta, i);
        }
        return summer;
    }
    /// Eqn. A.29
    TYPE d_etaI1_deta(TYPE eta, TYPE mbar){
        TYPE summer = 0;
        for (std::size_t j = 0; j < 7; ++j){
            summer += a_i(j, mbar)*(j+1.0)*pow(eta, j);
        }
        return summer;
    }
    /// Eqn. A.17
    TYPE I2(TYPE eta, TYPE mbar){
        TYPE summer = 0;
        for (std::size_t i = 0; i < 7; ++i){
            summer += b_i(i, mbar)*pow(eta, i);
        }
        return summer;
    }
    /// Eqn. A.30
    TYPE d_etaI2_deta(TYPE eta, TYPE mbar){
        TYPE summer = 0;
        for (std::size_t j = 0; j < 7; ++j){
            summer += b_i(j, mbar)*(j+1.0)*pow(eta, j);
        }
        return summer;
    }
    /// Eqn. A.11
    /// Erratum: should actually be 1/RHS of equation A.11 according to sample
    /// FORTRAN code
    TYPE C1(TYPE eta){
        return 1.0/(1.0
            + mbar*(8.0*eta-2.0*eta*eta)/pow(1.0-eta, 4)
            + (1.0-mbar)*(20.0*eta - 27.0*eta*eta + 12.0*pow(eta, 3) - 2.0*pow(eta, 4))/pow((1.0 - eta)*(2.0 - eta), 2));
    }
    /// Eqn. A.31
    TYPE C2(TYPE eta){
        return -pow(C1(eta), 2)*(
            mbar*(-4.0*eta*eta+20.0*eta+8.0)/pow(1.0-eta, 5)
            +(1.0-mbar)*(2.0*eta*eta*eta+12.0*eta*eta-48.0*eta+40.0)/pow((1.0-eta)*(2.0-eta), 3)
        );
    }
};

template <typename TYPE>
void do_calc(){
    std::vector<std::string> names(2, "Methane"); names[1] = "Ethane";
    std::vector<TYPE> z(2, 0);
    for (double z0 = 0.001; z0 < 1.0; z0 += 0.01){
        z[0] = z0; z[1] = 1.0-z0;
        PCSAFTMixture<TYPE> mix(names, z);
        mix.init();
        TYPE rhomolar = 3000.0, T = 200.0;
        double h = 1e-100;
        if constexpr (std::is_same<TYPE, std::complex<double>>::value){
            T += TYPE(0.0, 1e-100);
        }
        std::cout << z0 <<  " " << mix.calc_p(rhomolar, T) << std::endl;
        auto val = mix.calc_p(rhomolar, T);
        if constexpr (std::is_same<TYPE, std::complex<double>>::value){
            std::cout << "dZdT: " << std::imag(val)/h << std::endl;
        }
    }
}
int main(){
    do_calc<double>();
    do_calc<std::complex<double> >();
    auto x = ChebTools::ChebyshevExpansion::from_powxn(1, 0, 4000);
}