#include "TMS_Kalman.h"

TMS_Kalman::TMS_Kalman() : 
  Bethe(Material::kPolyStyrene),
  MSC(Material::kPolyStyrene),
  ForwardFitting(false) {
}

// Take a collection of hits
TMS_Kalman::TMS_Kalman(std::vector<TMS_Hit> &Candidates) : 
  Bethe(Material::kPolyStyrene), 
  MSC(Material::kPolyStyrene),
  ForwardFitting(false) {

  int nCand = Candidates.size();
  // Start at 20 MeV
  mom = 20;
  // And muon mass
  mass = BetheBloch_Utils::Mm;
  // Variance in energy
  total_en_var = 0;
  en = sqrt(mom*mom+mass*mass);
  for (int i = 0; i < nCand; ++i) {
    Predict(Candidates[i], Candidates[i+1]);
  }

  std::cout << "Total momentum: " << mom << std::endl;
  std::cout << "Total energy: " << en << std::endl;
  std::cout << "Total energy variation: " << total_en_var << std::endl;
}

// Update to the next step
void TMS_Kalman::Update() {

}

// Predict the next step
void TMS_Kalman::Predict(TMS_Hit &hit1, TMS_Hit &hit2) {

  // Right out of MINOS
  // Can probably do this better; e.g. calculate noise matrix at the same time as updating others (avoid duplicate calculations)
  /*
  BuildNoiseMatrix();
  ExtrapolateCovariance();
  KalmanGain();
  UpdateStateVector();
  UpdateCovariance();
  MoveArrays();
  */


  // Update the energy
  //double mom = CurrentState.qp;

  // Cheat for now
  double xval = hit1.GetTrueHit().GetX();
  double yval = hit1.GetTrueHit().GetY();
  double zval = hit1.GetTrueHit().GetZ();

  double xval2 = hit2.GetTrueHit().GetX();
  double yval2 = hit2.GetTrueHit().GetY();
  double zval2 = hit2.GetTrueHit().GetZ();

  // Make TVector3s of the two points
  TVector3 start(xval,yval,zval); // Start
  TVector3 stop(xval2,yval2,zval2); // Stop

  std::cout << "Going from " << start.X() << " " << start.Y() << " " << start.Z() << std::endl;
  std::cout << "To " << stop.X() << " " << stop.Y() << " " << stop.Z() << std::endl;

  // Get the materials between the two points
  std::vector<std::pair<TGeoMaterial*, double> > Materials = TMS_Geom::GetInstance().GetMaterials(start, stop);

  std::cout << "Looping over " << Materials.size() << " materials" << std::endl;

  // Loop over the materials between the two projection points
  int counter = 0;
  for (auto material : Materials) {

    // Read these directly from a TGeoManager
    double density = material.first->GetDensity()/(CLHEP::g/CLHEP::cm3); // now in g/cm3 (edep-sim geometry is in CLHEP units)
    double thickness = material.second/10.; // in cm (was in mm in geometry)

    std::cout << "Material " << counter << " = " << material.first->GetName() << std::endl;
    std::cout << "  density: " << density << std::endl;
    std::cout << "  thickness: " << thickness << std::endl;
    std::cout << "  thickness*density = " << density*thickness << std::endl;

    // Skip if density or thickness is small
    if (density*thickness < 0.1) {
      std::cout << "  Skipping material, to little path length to bother" << std::endl;
      continue;
    }

    // Update the Bethe Bloch calculator to use this material
    Material matter(density);

    Bethe.fMaterial = matter;

    std::cout << "energy before: " << en << std::endl;
    // Subtract off the energy loss for this material
    if (ForwardFitting) en -= Bethe.Calc_dEdx(en)*density*thickness;
    else                en += Bethe.Calc_dEdx(en)*density*thickness;
    std::cout << "energy after: " << en << std::endl;

    // Updated momentum squared based on the energy loss
    //double momup2 = en*en-mass*mass;

    //double momup = momup2 > 0 ? sqrt(momup2) : 0;

    // Update the state vector
    //FutureState.qp = 1./momup;

    // Variance assuming Gaussian straggling
    double en_var = Bethe.Calc_dEdx_Straggling(en)*density*thickness;
    total_en_var += en_var*en_var;
    std::cout << "energy var: " << en_var*en_var << std::endl;

    // Set the material for the multiple scattering
    MSC.fMaterial = matter;
    // Calculate this before or after the energy subtraction/addition?
    MSC.Calc_MS(en, thickness*density);

    counter++;
  }

  // Updated momentum^2
  double p_2_up = en*en-BetheBloch_Utils::Mm*BetheBloch_Utils::Mm;
  double p_up;
  if (p_2_up > 0) p_up = sqrt(p_2_up);
  else {
    std::cerr << "negative momentum squared, setting momentum to 1 MeV" << std::endl;
    p_up = 1;
  }

  std::cout << "total energy variation: " << total_en_var << std::endl;
  // Momentum variance
  double p_var = (2*en/p_up)*(2*en/p_up) * total_en_var;
  double qp_var = 1./(p_up*p_up*p_up*p_up) * p_var;
  NoiseMatrix(4,4) = qp_var;

  // Finally when we've iterated through the materials we can get the combined effect from the multiple-scattering
  double ms = MSC.Calc_MS_Sigma();

  std::cout << "Multiple scattering: " << ms << std::endl;

  // Now proceed with Wolin and Ho (Nucl Inst A329 1993 493-500)
  // covariance for multiple scattering
  // Also see MINOS note on Kalman filter (John Marshall, Nov 15 2005)
  /*
  double ax2 = state.ax*state.ax;
  double ay2 = state.ay*state.ay;
  */

}


