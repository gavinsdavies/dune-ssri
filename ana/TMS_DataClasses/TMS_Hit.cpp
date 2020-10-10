#include "TMS_Hit.h"

// The constructor for a hit in the TMS, from a edep-sim hit
TMS_Hit::TMS_Hit(double x, double y, double z, double t, double E) {

  // Save the true hit
  TrueHit = TMS_TrueHit(x, y, z, t, E);

  // Save the true particle
  //TMS_TrueParticle TrueParticle;

  // Save the bar
  Bar = FindBar(x,y,z);

  // Save the energy deposit
  SetE(E);
}

TMS_Hit::TMS_Hit(TG4HitSegment &edep_seg) {
  // Save the true hit
  //TrueHit = TMS_TrueHit(edep_seg);

  // Save the energy deposit
  EnergyDeposit = edep_seg.GetEnergyDeposit();

  // The true particle
  TrueParticle = TMS_TrueParticle(edep_seg);

  // Get the bar
  Bar = TMS_Bar(edep_seg);
}

// Find which bar a given x,y,z position corresponds to
// Maybe this function should be moved to the singleton instead
TMS_Bar TMS_Hit::FindBar(double x, double y, double z) {

  // Use the ROOT geometry to figure it out if available
  TGeoManager *geom = TMS_Geom::GetInstance().GetGeometry();

  geom->FindNode(x,y,z);
  TGeoNavigator *nav = geom->GetCurrentNavigator();
  std::string NodeName = std::string(nav->GetCurrentNode()->GetName());
  // cd up in the geometry to find the right name
  while (NodeName.find("modeulelayervol_PV") == std::string::npos && NodeName.find("volWorld") == std::string::npos) {
    nav->CdUp();
    NodeName = std::string(nav->GetCurrentNode()->GetName());
  }

  // The bar that we return
  TMS_Bar Bar;

  // If we've reached the world volume we don't have a scintillator hit -> return some mad bad value
  if (NodeName.find("volWorld") != std::string::npos) {
    // Since the bar has already been created as a "error" in the above empty constructor we can just return
    return Bar;
  }

  // Get the node
  //TGeoNode *node = nav->GetCurrentNode();
  // Get the volume
  //TGeoVolume *vol = node->GetVolume();
  //TGeoShape *shape = vol->GetShape();

  // Assume we now have the bar information in the navigator
  //x = ;

  //int BarNumber = nav->GetCurrentNode()->GetNumber();

  // If ROOT geometry isn't available, use this hard-coded mess
  return Bar;
}
