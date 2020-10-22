#include "TMS_Hit.h"

// The constructor for a hit in the TMS, from a edep-sim hit
/*
TMS_Hit::TMS_Hit(double x, double y, double z, double t, double E) {

  // Save the true hit
  TrueHit = TMS_TrueHit(x, y, z, t, E);

  // Save the true particle
  //TMS_TrueParticle TrueParticle;

  // Save the bar
  Bar = TMS_Bar(x,y,z);

  // Save the energy deposit
  SetE(E);

  SetT(t);
}
*/

// Set the bar and truehit
TMS_Hit::TMS_Hit(TG4HitSegment &edep_seg) : 
  TrueHit(edep_seg), 
  Bar(edep_seg),
  EnergyDeposit(edep_seg.GetEnergyDeposit()),
  // Define time as the average between start and stop of hit
  Time((edep_seg.GetStop().T()+edep_seg.GetStart().T())/2)
{


  // The true particle
  //TrueParticle = TMS_TrueParticle(edep_seg);

}

void TMS_Hit::Print() {
  std::cout << "Printing TMS hit" << std::endl;
  std::cout << "EnergyDeposit: " << EnergyDeposit << std::endl;
  std::cout << "Time: " << Time << std::endl;
  std::cout << "Bar: " << std::endl;
  Bar.Print();

  std::cout << "TrueHit: " << std::endl;
  TrueHit.Print();
}
