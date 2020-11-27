#ifndef _BETHEBLOCH_H_SEEN_
#define _BETHEBLOCH_H_SEEN_

#include <cmath>

// Almost entirely from the PDG
// and pdg.lbl.gov/AtomicNuclearProperties
namespace BetheBloch_Utils {

  // Fine structure
  const double a = 1/137.035999139;
  // Fine structure squared
  const double a_2 = a*a;
  // Electron mass in MeV
  const double Me = 0.5109989461;
  // Squared electron mass in MeV2
  const double Me_2 = Me*Me;
  // Muon mass in MeV
  const double Mm = 105.658389;
  // Muon mass in MeV2
  const double Mm_2 = Mm*Mm;
  // Compton wavelength of electron (1E11 cm);
  const double Le = 3.8616;
  // Compton wavelength of electron2 (1E22 cm2);
  const double Le_2 = Le*Le;
  // good old pi
  const double pi = M_PI;
  // Avogadros number in 10E-23
  const double Na = 6.0222140857;

  // Sternheimer factor
  const double Stern2ln10 = 2*log(10);

  // coefficient for dE/dx in /mol cm2
  const double K = 0.307075;

  inline double RelativisticBeta(double m, double E) {
    if (E < m) return 0;

    double E2 = E*E;
    double m2 = m*m;

    return sqrt(E2-m2)/E;
  }

  inline double RelativisticGamma(double m, double E) {
    if (E < m) return 0;

    return E/m;
  }

  // Convert energy to momentum
  inline double EnergyToMomentum(double m, double E) {
    if (E < m) return 0;
    return sqrt( E*E - m*m );
  }

  // Maximum energy transfer to electron from incoming muon (causing ionisation), in MeV
  inline double MaximumEnergyTransfer(double E) {
    double me = Me;
    double me_2 = Me_2;
    double mm = Mm;
    double mm_2 = Mm_2;
    double p = EnergyToMomentum(mm, E);
    double p_2 = p*p;

    return 2 * me * p_2 / ( me_2 + mm_2 + 2*me*E );
  }

};


// Generic material class
// Only support iron and polystyrene for now
// Taken right from LBL and PDG
class Material {

  public:

    enum MaterialType { kPolyStyrene, kIron, kGraphite, kGArgon, kLArgon, kWater, kUnknown };

    std::string MaterialName() {

      switch (fMaterialType) {
        case kPolyStyrene:
          return "Polystyrene";
          break;

        case kIron:
          return "Iron";
          break;

        case kGraphite:
          return "Graphite";
          break;

        case kGArgon:
          return "GAr";
          break;

        case kLArgon:
          return "LAr";
          break;

        case kWater:
          return "H_{2}O";
          break;


        default:
          return "unknown";
          break;
      }

      return "unknown";
    }

    Material(MaterialType type) {
      fMaterialType = type;

      switch (type) {
        // Polystyrene: https://pdg.lbl.gov/2020/AtomicNuclearProperties/MUE/muE_polystyrene.pdf
        case kPolyStyrene:
          Z_A = 0.53768;
          rho = 1.060;
          I = 68.7;
          a = 0.16454;
          m = 3.2224;
          x0 = 0.1647;
          x1 = 2.5031;
          Cbar = 3.2999;
          d0 = 0.00;
          break;

          // Iron: https://pdg.lbl.gov/2020/AtomicNuclearProperties/MUE/muE_iron_Fe.pdf
        case kIron:
          Z_A = 26/55.845;
          rho = 7.874;
          I = 286.0;
          a = 0.14680;
          m = 2.9632;
          x0 = -0.0012;
          x1 = 3.1531;
          Cbar = 4.2911;
          d0 = 0.12;
          break;

          // Graphite
          // https://pdg.lbl.gov/2020/AtomicNuclearProperties/MUE/muE_carbon_graphite_C.pdf
        case kGraphite:
          Z_A = 6/12.0107;
          rho = 2.210;
          I = 78.0;
          a = 0.20762;
          m = 2.9532;
          x0 = -0.0090;
          x1 = 2.4817;
          Cbar = 2.8926;
          d0 = 0.14;
          break;

          // Gaseous Argon
          // https://pdg.lbl.gov/2020/AtomicNuclearProperties/MUE/muE_argon_gas_Ar.pdf
        case kGArgon:
          Z_A = 18/39.948;
          rho = 1.662E-3;
          I = 188;
          a = 0.19714;
          m = 2.9618;
          x0 = 1.7635;
          x1 = 4.4855;
          Cbar = 11.9480;
          d0 = 0.00;
          break;

          // Liquid Argon
          // https://pdg.lbl.gov/2020/AtomicNuclearProperties/MUE/muE_liquid_argon.pdf
        case kLArgon:
          Z_A = 18/39.948;
          rho = 1.396;
          I = 188;
          a = 0.19559;
          m = 3.0000;
          x0 = 0.2000;
          x1 = 3.0000;
          Cbar = 5.2146;
          d0 = 0.00;
          break;

          // Liquid water
          // https://pdg.lbl.gov/2020/AtomicNuclearProperties/MUE/muE_water_liquid.pdf
        case kWater:
          Z_A = 0.55509;
          rho = 1.000;
          I = 79.7;
          a = 0.09116;
          m = 3.4773;
          x0 = 0.2400;
          x1 = 2.8004;
          Cbar = 3.5017;
          d0 = 0.00;
          break;


        default:
          std::cerr << "Material not supported" << std::endl;
          throw;
      }
    }

    // Z/A
    double Z_A;
    // Ionisation
    double I; // eV
    // density
    double rho; // g/cm3

    // Numbers for density correction factors to Bethe
    double x0;
    double x1;
    double a;
    double m;
    double Cbar;
    double d0;

    MaterialType fMaterialType;


};

class BetheBloch_Calculator {

  public:

    BetheBloch_Calculator() = delete;

    BetheBloch_Calculator(Material::MaterialType type) :
      fMaterial(type) {
      };

    // Density correction factor a la PDG (Sternheimer) MeV
    // eq 34.7 pdg
    inline double DensityCorrectionFactor(double E) {
      // Update from muons in Iron from pdg
      const double x0 = fMaterial.x0;
      const double x1 = fMaterial.x1;
      const double a = fMaterial.a;
      const double m = fMaterial.m;
      const double Cbar = fMaterial.Cbar;
      const double d0 = fMaterial.d0;

      const double beta = BetheBloch_Utils::RelativisticBeta(BetheBloch_Utils::Mm, E);
      const double gamma = BetheBloch_Utils::RelativisticGamma(BetheBloch_Utils::Mm, E);
      const double X = log10( beta*gamma );

      if (x0 < X && X < x1) {
        return ( BetheBloch_Utils::Stern2ln10 * X - Cbar + a * pow(x1-X,m) );
      }
      if (X > x1) {
        return ( BetheBloch_Utils::Stern2ln10 * X - Cbar );
      }

      // Conductor
      if (X < x0) {
        return d0*pow(10,2*(X-x0));
      }

      return 0;
    }

    // Calculate the ionsiation for muons in Bethe Bloch in MeV/(gr/cm2)
    inline double Calc_dEdx(double E) {
      // z over A
      double Z_A = fMaterial.Z_A;
      double beta = BetheBloch_Utils::RelativisticBeta(BetheBloch_Utils::Mm, E);
      double beta_2 = beta*beta;
      double gamma = BetheBloch_Utils::RelativisticGamma(BetheBloch_Utils::Mm, E);
      double gamma_2 = gamma*gamma;
      // Ionisation from PDG
      double I = fMaterial.I * 1E-6; //286E-6; // 286 eV -> 286E-6 in MeV
      double I_2 = I*I;
      double Em = BetheBloch_Utils::MaximumEnergyTransfer(E);
      double Em_2 = Em*Em;
      double E_2 = E*E;
      double d = DensityCorrectionFactor(E);

      // K = 4pi NA re2 me c2
      double de_dx = BetheBloch_Utils::K * Z_A * (1/beta_2) *
        (0.5*log( 2*BetheBloch_Utils::Me*beta_2*gamma_2*Em/I_2 ) - beta_2 - 0.5*d );

      // Correction term in Groom not present in PDG
      //double de_dx2 = 10* a_2 * 2*pi * Na * Le_2 * Z_A * (Me / beta_2) *
      //( log( 2*Me*beta_2*gamma_2*Em/I_2 ) - 2*beta_2 + 0.25*(Em_2/E_2) - d );

      return de_dx; // in MeV/(gr/cm^2)
    }

    // Calculate the ionsiation for muons in Bethe Bloch in MeV/(gr/cm2)
    // Gaussian CLT
    inline double Calc_dEdx_Straggling(double E) {
      // z over A
      double Z_A = fMaterial.Z_A;
      double beta = BetheBloch_Utils::RelativisticBeta(BetheBloch_Utils::Mm, E);
      double beta_2 = beta*beta;
      double Em = BetheBloch_Utils::MaximumEnergyTransfer(E);
      double E_2 = E*E;
      double thick = 1.0;
      double thickness_den = fMaterial.rho * thick;

      double Xi = (BetheBloch_Utils::K/2*Z_A*1/beta_2)*thickness_den;
      double kappa = Xi/Em;
      double sigma = Xi*sqrt(1-beta_2)/2*kappa;

      return sigma;
    }


    // Calculate the ionsiation for muons in Bethe Bloch MeV/(gr/cm2)
    // eq 34.12 pdg
    inline double Calc_dEdx_mostprob(double E) {
      // z over A
      double Z_A = fMaterial.Z_A;
      double beta = BetheBloch_Utils::RelativisticBeta(BetheBloch_Utils::Mm, E);
      double beta_2 = beta*beta;
      double gamma = BetheBloch_Utils::RelativisticGamma(BetheBloch_Utils::Mm, E);
      double gamma_2 = gamma*gamma;
      // Ionisation from PDG
      double I = fMaterial.I*1E-6; //286E-6; // 286 eV -> 286E-6 in MeV
      double E_2 = E*E;
      double d = DensityCorrectionFactor(E);

      // Set some thickness
      //const double thick = 1./7.85; // thickness in g/cm2 -> try 1cm with 7.85 g/cm3 density
      //const double thick = 7.85;
      const double thick = fMaterial.rho * 1; // 1 cm of material
      const double j = 0.200; // from pdg

      double eps = (BetheBloch_Utils::K/2)*Z_A*thick/beta_2; // convenient

      double de_dx = eps *
        (log( 2*BetheBloch_Utils::Me*beta_2*gamma_2/I ) + log(eps/I) + j - beta_2 - d );

      //double de_dx2 = 10* a_2 * 2*pi * Na * Le_2 * Z_A * (Me / beta_2) *
      //( log( 2*Me*beta_2*gamma_2*Em/I_2 ) - 2*beta_2 + 0.25*(Em_2/E_2) - d );

      return de_dx; // in MeV/(gr/cm^2)
    }

    Material fMaterial;
};

#endif
