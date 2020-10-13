#include "TMS_Bar.h"

TMS_Bar::TMS_Bar() {
  // Default to an error
  BarOrient = kError;
  x = -999.9999;
  y = -999.9999;
  z = -999.9999;
  xw = -999.9999;
  yw = -999.9999;
  zw = -999.9999;
  BarNumber = -999;
}

TMS_Bar::TMS_Bar(double xpos, double ypos, double zpos) {
  x = xpos;
  y = ypos;
  z = zpos;
  // Figure out what to do here... It neds to find itself?
  BarNumber = -999;
  PlaneNumber = FindBar(x,y,z);
}

TMS_Bar::TMS_Bar(double xpos, double ypos, double zpos, int BarNo) {
  x = xpos;
  y = ypos;
  z = zpos;
  BarNumber = BarNo;
  PlaneNumber = FindBar(x,y,z);
}

// Construct a bar from a hit
TMS_Bar::TMS_Bar(TG4HitSegment &edep_seg) {
  double avgx = ( edep_seg.GetStart().X() + edep_seg.GetStop().X() ) / 2;
  double avgy = ( edep_seg.GetStart().Y() + edep_seg.GetStop().Y() ) / 2;
  double avgz = ( edep_seg.GetStart().Z() + edep_seg.GetStop().Z() ) / 2;

  PlaneNumber = FindBar(avgx, avgy, avgz);
  x = avgx;
  y = avgy;
  z = avgz;

  // Now set the width and orientation of the bar
  if (PlaneNumber % 2 == 0) BarOrient = kXBar;
  else BarOrient = kYBar;

  BarNumber = -1;
  GlobalBarNumber = -1;

  xw = -1;
  yw = -1;
  zw = -1;
}

// Find which bar a given x,y,z position corresponds to
// Maybe this function should be moved to the singleton instead
int TMS_Bar::FindBar(double x, double y, double z) {

  // Use the ROOT geometry to figure it out if available
  TGeoManager *geom = TMS_Geom::GetInstance().GetGeometry();

  // Find which node this position is equivalent too
  geom->FindNode(x,y,z);
  TGeoNavigator *nav = geom->GetCurrentNavigator();
  std::string NodeName = std::string(nav->GetCurrentNode()->GetName());
  // cd up in the geometry to find the right name
  while (NodeName.find("modulelayervol_PV") == std::string::npos && NodeName.find("volWorld") == std::string::npos) {
    nav->CdUp();
    NodeName = std::string(nav->GetCurrentNode()->GetName());
  }

  // If we've reached the world volume we don't have a scintillator hit -> return some mad bad value
  if (NodeName.find("volWorld") != std::string::npos) {
    // Since the bar has already been created as a "error" in the above empty constructor we can just return
    std::cout << "Bar position not found in TMS_Bar::FindBar!" << std::endl;
    return -1;
  }

  // Get the node
  //TGeoNode *node = nav->GetCurrentNode();
  // Get the volume
  //TGeoVolume *vol = node->GetVolume();
  //TGeoShape *shape = vol->GetShape();

  // Assume we now have the bar information in the navigator
  //x = ;

  int Number = nav->GetCurrentNode()->GetNumber();
  std::cout << Number << std::endl;

  // If ROOT geometry isn't available, use this hard-coded mess
  return Number;
}

std::string TMS_Bar::BarType_ToString(BarType bar) {
  if (bar == kXBar) return std::string("X-bar");
  else if (bar == kYBar) return std::string("Y-bar");
  else if (bar == kUBar) return std::string("U-bar");
  else if (bar == kVBar) return std::string("V-bar");
  return std::string("ERROR");
}

void TMS_Bar::Print() {
  std::cout << "Printing TMS bar: " << std::endl;
  std::cout << "(x, y, z) = " << "(" << x << ", " << y << ", " << z << ")" << std::endl;
  std::cout << "(xw, yw, zw) = " << "(" << xw << ", " << yw << ", " << zw << ")" << std::endl;
  std::cout << "BarOrient: " << BarType_ToString(BarOrient) << " (" << BarOrient << ")" << std::endl;
  std::cout << "PlaneNumber: " << PlaneNumber << std::endl;
  std::cout << "BarNumber: " << BarNumber << std::endl;
  std::cout << "GlobalBarNumber: " << GlobalBarNumber << std::endl;
}
