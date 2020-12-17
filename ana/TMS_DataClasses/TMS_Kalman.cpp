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
  ForwardFitting(false),
  Talk(false)
  //Talk(true)
{

  // Empty the KalmanStates
  KalmanNodes.empty();

  // Save the number of initial candidates
  int nCand = Candidates.size();
  // And muon mass
  mass = BetheBloch_Utils::Mm;

  // Make a new Kalman state for each hit
  KalmanNodes.reserve(nCand);
  double PreviousZ = -999;
  for (int i = 0; i < nCand; ++i) {
    TMS_Hit hit = Candidates[i];
    //if (hit.GetBar().GetBarType() != TMS_Bar::kYBar) continue;
    double x = hit.GetTrueHit().GetX();
    double y = hit.GetTrueHit().GetY();
    double z = hit.GetTrueHit().GetZ();

    double DeltaZ = z-PreviousZ;

    // This also initialises the state vectors in each of the nodes
    TMS_KalmanNode Node(x,y,z,DeltaZ);
    KalmanNodes.emplace_back(std::move(Node));

    PreviousZ = z;
  }

  RunKalman();

}

// Update to the next step
void TMS_Kalman::RunKalman() {

  int nCand = KalmanNodes.size();
  for (int i = 0; i < nCand; ++i) {
    Predict(KalmanNodes[i]);
  }
  std::cout << "Momentum of last hit: " << 1./KalmanNodes.back().CurrentState.qp << std::endl;
}

// Predict the next step
// Here we use the previous state to predict the current state
// We also calculate the updated noise matrix from multiple scattering and energy loss
// Account for energy loss, multiple scattering, and bending due to the magnetic field
void TMS_Kalman::Predict(TMS_KalmanNode &Node) {

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

  TMS_KalmanState &PreviousState = Node.PreviousState;
  TMS_KalmanState &CurrentState = Node.CurrentState;

  // Propagate the current state
  TMatrixD &Transfer = Node.TransferMatrix;
  TVectorD PreviousVec(5);
  PreviousVec[0] = PreviousState.x;
  PreviousVec[1] = PreviousState.y;
  PreviousVec[2] = PreviousState.dxdz;
  PreviousVec[3] = PreviousState.dydz;
  PreviousVec[4] = PreviousState.qp;
  // Just the propagator matrix influence
  TVectorD UpdateVec = Transfer*(PreviousVec);
  // Now construct the current state (z of CurrentState is already set to be z+dz)
  CurrentState.x = UpdateVec[0];
  CurrentState.y = UpdateVec[1];
  CurrentState.dxdz = UpdateVec[2];
  CurrentState.dydz = UpdateVec[3];
  // Don't update q/p until later (when we've done the energy loss calculation)

  // Update the energy
  double mom = 1./PreviousState.qp;
  // Initial energy before energy loss
  double en_initial = sqrt(mom*mom+mass*mass);
  // The energy we'll be changing
  double en = en_initial;

  if (Talk) std::cout << "mom: " << mom << std::endl;

  // Cheat for now
  double xval = PreviousState.x;
  double yval = PreviousState.y;
  double zval = PreviousState.z;

  double xval2 = CurrentState.x;
  double yval2 = CurrentState.y;
  double zval2 = CurrentState.z;

  // Make TVector3s of the two points
  TVector3 start(xval,yval,zval); // Start
  TVector3 stop(xval2,yval2,zval2); // Stop

  if (Talk) {
    std::cout << "Going from " << start.X() << " " << start.Y() << " " << start.Z() << std::endl;
    std::cout << "To " << stop.X() << " " << stop.Y() << " " << stop.Z() << std::endl;
  }

  // Get the materials between the two points
  std::vector<std::pair<TGeoMaterial*, double> > Materials = TMS_Geom::GetInstance().GetMaterials(start, stop);

  if (Talk) std::cout << "Looping over " << Materials.size() << " materials" << std::endl;
  double TotalPathLength = 0;
  double TotalLength = 0;

  // Loop over the materials between the two projection points
  int counter = 0;
  double total_en_var = 0;
  if (Talk) std::cout << "Energy before: " << en_initial << std::endl;
  for (auto material : Materials) {

    // Read these directly from a TGeoManager
    double density = material.first->GetDensity()/(CLHEP::g/CLHEP::cm3); // now in g/cm3 (edep-sim geometry is in CLHEP units)
    double thickness = material.second/10.; // in cm (was in mm in geometry)

    if (Talk) {
      std::cout << "Material " << counter << " = " << material.first->GetName() << std::endl;
      std::cout << "  density: " << density << std::endl;
      std::cout << "  thickness: " << thickness << std::endl;
      std::cout << "  thickness*density = " << density*thickness << std::endl;
    }

    // Skip if density or thickness is small
    if (density*thickness < 0.1) {
      if (Talk) std::cout << "  Skipping material, to little path length to bother" << std::endl;
      continue;
    }

    TotalPathLength += density*thickness;
    TotalLength += thickness;

    // Update the Bethe Bloch calculator to use this material
    Material matter(density);

    Bethe.fMaterial = matter;

    if (Talk) {
      std::cout << "Mean: " << Bethe.Calc_dEdx(en)*density*thickness << std::endl;
      std::cout << "Most prob: " << Bethe.Calc_dEdx_mostprob(en)*density*thickness << std::endl;
      std::cout << "energy before: " << en << std::endl;
    }

    // Subtract off the energy loss for this material
    if (ForwardFitting) en -= Bethe.Calc_dEdx(en)*density*thickness;
    else                en += Bethe.Calc_dEdx(en)*density*thickness;
    if (Talk) std::cout << "energy after: " << en << std::endl;


    // Variance assuming Gaussian straggling
    double en_var = Bethe.Calc_dEdx_Straggling(en)*density*thickness;
    total_en_var += en_var*en_var;
    if (Talk) std::cout << "energy var: " << en_var*en_var << std::endl;

    // Set the material for the multiple scattering
    MSC.fMaterial = matter;
    // Calculate this before or after the energy subtraction/addition?
    MSC.Calc_MS(en, thickness*density);

    counter++;
  }
  if (Talk) {
    std::cout << "Total path length (g/cm2): " << TotalPathLength << std::endl;
    std::cout << "Total length (cm): " << TotalLength << std::endl;
  }

  // Updated momentum^2
  double p_2_up = en*en-BetheBloch_Utils::Mm*BetheBloch_Utils::Mm;
  double p_up;
  if (p_2_up > 0) p_up = sqrt(p_2_up);
  else {
    std::cerr << "negative momentum squared, setting momentum to 1 MeV" << std::endl;
    p_up = 1;
  }

  // Update the state's q/p
  CurrentState.qp = 1./p_up;

  if (Talk) {
    std::cout << "total energy variation: " << total_en_var << std::endl;
    std::cout << "energy after: " << en << std::endl;
    std::cout << "Delta(Energy): " << en-en_initial << std::endl;
  }

  double p_var = (2*en/p_up)*(2*en/p_up) * total_en_var;
  double qp_var = 1./(p_up*p_up*p_up*p_up) * p_var;

  if (Talk) {
    std::cout << "momentum after: " << p_up << std::endl;
    std::cout << "total momentum variation: " << p_var << std::endl;
  }
  // Momentum variance

  // Set the noise matrix
  TMatrixD &NoiseMatrix = Node.NoiseMatrix;

  NoiseMatrix(4,4) = qp_var;

  // Finally when we've iterated through the materials we can get the combined effect from the multiple-scattering
  double ms = MSC.Calc_MS_Sigma();

  if (Talk) std::cout << "Multiple scattering: " << ms << std::endl;

  // Now proceed with Wolin and Ho (Nucl Inst A329 1993 493-500)
  // covariance for multiple scattering
  // Also see MINOS note on Kalman filter (John Marshall, Nov 15 2005)
  double ax2 = CurrentState.dxdz; // get the unit vectors
  double ay2 = CurrentState.dydz; // get the unit vectors

  double norm = 1+ax2+ay2; // 1+P3^2+P4^2 in eq 16, 17, 18 in Wolin and Ho
  double covAxAx = norm*ms*(1+ax2);// eq 16 Wolin and Ho
  double covAyAy = norm*ms*(1+ay2);// eq 17 Wolin and Ho
  double covAxAy = norm*ms*ax2*ay2;// eq 18 Wolin and Ho

  // Check the z0 variable!
  if (!ForwardFitting) TotalPathLength = -1*TotalPathLength;

  // Build the covariance matrix in Wolin and Ho after eq 18
  // Equation 15 in Robert Harr Calculation of Track and Vertex Errors for Detector Design Studies
  // Also equation 48; kappa in Harr is "norm" here
  double TotalPathLengthSq = TotalPathLength*TotalPathLength;
  NoiseMatrix(0,0) = covAxAx * TotalPathLengthSq / 3.;
  NoiseMatrix(1,1) = covAyAy * TotalPathLengthSq / 3.;
  NoiseMatrix(2,2) = covAxAx;
  NoiseMatrix(3,3) = covAyAy;

  // Not sure about the negative signs...? Wolin and Ho have them in, Harr doesn't
  NoiseMatrix(1,0) = NoiseMatrix(0,1) =    covAxAy * TotalPathLengthSq/3.;
  NoiseMatrix(2,0) = NoiseMatrix(0,2) = -1*covAxAx * TotalPathLength/2.;
  NoiseMatrix(3,0) = NoiseMatrix(0,3) = -1*covAxAy * TotalPathLength/2.;

  NoiseMatrix(2,1) = NoiseMatrix(1,2) = -1*covAxAy * TotalPathLength/2.;
  NoiseMatrix(3,1) = NoiseMatrix(3,1) = -1*covAyAy * TotalPathLength/2.;

  NoiseMatrix(3,2) = NoiseMatrix(2,3) =    covAxAy;

  // I think that's it
  // Other than the B-field...
}


