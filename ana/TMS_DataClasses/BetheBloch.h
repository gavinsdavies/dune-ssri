#include <iostream>
#include <cmath>
#include <assert.h>
#include "TGraph.h"
#include "TCanvas.h"

// Almost entirely from the PDG
// and pdg.lbl.gov/AtomicNuclearProperties
namespace BetheBloch {

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
  double MaximumEnergyTransfer(double E) {
    double me = Me;
    double me_2 = Me_2;
    double mm = Mm;
    double mm_2 = Mm_2;
    double p = EnergyToMomentum(mm, E);
    double p_2 = p*p;

    return 2 * me * p_2 / ( me_2 + mm_2 + 2*me*E );
  }

  // Density correction factor a la PDG (Sternheimer)
  // eq 34.7 pdg
  // MeV
  double DensityCorrectionFactor(double E) {
    /*
       const double X0 = 0.10;
       const double X1 = 3;
       const double a = 0.127;
       const double m = 3.29;
       const double C = -4.62;
       */
    // Update from muons in Iron from pdg
    const double X0 = -0.0012;
    const double X1 = 3.1531;
    const double a = 0.14680;
    const double m = 3.29;
    const double Cbar = 4.2911;
    const double d0 = 0.12;

    const double beta = RelativisticBeta(Mm, E);
    const double gamma = RelativisticGamma(Mm, E);
    const double X = log10( beta*gamma );

    if (X0 < X && X < X1) {
      return ( Stern2ln10 * X - Cbar + a * pow(X1-X,m) );
    }
    if (X > X1) {
      return ( Stern2ln10 * X - Cbar );
    }

    // Conductor
    if (X < X0) {
      return d0*pow(10,2*(X-X0));
    }

    return 0;
  }

  // Calculate the ionsiation for muons in Bethe Bloch
  // MeV/(gr/cm2)
  double Calc_dEdx(double E) {
    // z over A
    double Z_A = 26./56;
    double beta = RelativisticBeta(Mm, E);
    double beta_2 = beta*beta;
    double gamma = RelativisticGamma(Mm, E);
    double gamma_2 = gamma*gamma;
    // Ionisation from PDG
    double I = 286E-6; // 286 eV -> 286E-6 in MeV
    double I_2 = I*I;
    double Em = MaximumEnergyTransfer(E);
    double Em_2 = Em*Em;
    double E_2 = E*E;
    double d = DensityCorrectionFactor(E);

    // K = 4pi NA re2 me c2
    double de_dx = K * Z_A * (1/beta_2) *
      (0.5*log( 2*Me*beta_2*gamma_2*Em/I_2 ) - beta_2 - 0.5*d );

    //double de_dx2 = 10* a_2 * 2*pi * Na * Le_2 * Z_A * (Me / beta_2) *
    //( log( 2*Me*beta_2*gamma_2*Em/I_2 ) - 2*beta_2 + 0.25*(Em_2/E_2) - d );

    return de_dx; // in MeV/(gr/cm^2)
  }

  // Calculate the ionsiation for muons in Bethe Bloch
  // 34.12 pdg
  // MeV/(gr/cm2)
  double Calc_dEdx_mostprob(double E) {
    // z over A
    double Z_A = 26./56;
    double beta = RelativisticBeta(Mm, E);
    double beta_2 = beta*beta;
    double gamma = RelativisticGamma(Mm, E);
    double gamma_2 = gamma*gamma;
    // Ionisation from PDG
    double I = 286E-6; // 286 eV -> 286E-6 in MeV
    double E_2 = E*E;
    double d = DensityCorrectionFactor(E);
    //const double thick = 1./7.85; // thickness in g/cm2 -> try 1cm with 7.85 g/cm3 density
    //const double thick = 7.85;
    const double thick = 1.0;
    const double j = 0.200; // from pdg

    double eps = (K/2)*Z_A*thick/beta_2; // convenient

    // K = 4pi NA re2 me c2
    double de_dx = eps *
      (log( 2*Me*beta_2*gamma_2/I ) + log(eps/I) + j - beta_2 - d );

    //double de_dx2 = 10* a_2 * 2*pi * Na * Le_2 * Z_A * (Me / beta_2) *
    //( log( 2*Me*beta_2*gamma_2*Em/I_2 ) - 2*beta_2 + 0.25*(Em_2/E_2) - d );

    return de_dx; // in MeV/(gr/cm^2)
  }

}
