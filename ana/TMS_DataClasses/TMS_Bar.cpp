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
  //PlaneNumber = FindBar(x,y,z);
}

// Construct a bar from a hit
TMS_Bar::TMS_Bar(TG4HitSegment &edep_seg) {

  BarNumber = -1;
  GlobalBarNumber = -1;
  PlaneNumber = -1;

  xw = -1;
  yw = -1;
  zw = -1;

  // Get the average distance
  double avgx = ( edep_seg.GetStart().X() + edep_seg.GetStop().X() ) / 2;
  double avgy = ( edep_seg.GetStart().Y() + edep_seg.GetStop().Y() ) / 2;
  double avgz = ( edep_seg.GetStart().Z() + edep_seg.GetStop().Z() ) / 2;

  x = avgx;
  y = avgy;
  z = avgz;

  // Find the bar in the geometry
  FindModules(avgx, avgy, avgz);

  // Then correct the position of the bar?

  // Now set the width and orientation of the bar
  if (PlaneNumber % 2 == 0) BarOrient = kXBar;
  else BarOrient = kYBar;

}

// Find which bar a given x,y,z position corresponds to
// Maybe this function should be moved to the singleton instead
bool TMS_Bar::FindModules(double xval, double yval, double zval) {

  // Use the ROOT geometry to figure it out if available
  TGeoManager *geom = TMS_Geom::GetInstance().GetGeometry();

  // Find which node this position is equivalent too
  std::string NodeName = std::string(geom->FindNode(xval,yval,zval)->GetName());

  // The position of the hit bar
  double Translation[] = {0., 0., 0.};
  for (int i = 0; i < 3; ++i) Translation[i] = geom->GetCurrentMatrix()->GetTranslation()[i];

  // cd up in the geometry to find the right name
  while (NodeName.find(TMS_Const::TMS_TopLayerName) == std::string::npos) {

    // We've found the plane number
    if (NodeName.find(TMS_Const::TMS_ModuleLayerName) != std::string::npos) {
      PlaneNumber = geom->GetCurrentNode()->GetNumber();
    }

    // This is the furthest down hit we have
    else if (NodeName.find(TMS_Const::TMS_ScintLayerName) != std::string::npos) {
      BarNumber = geom->GetCurrentNode()->GetNumber();
    }

    else if (NodeName.find(TMS_Const::TMS_ModuleName) != std::string::npos) {
      GlobalBarNumber = geom->GetCurrentNode()->GetNumber();
    }

    geom->CdUp();
    NodeName = std::string(geom->GetCurrentNode()->GetName());
  }

  // Update the hit value to be the bar
  x = Translation[0];
  // y-bar information doesn't seem to work?
  //y = Translation[1];
  z = Translation[2];

  // If we've reached the world volume we don't have a scintillator hit -> return some mad bad value
  if (BarNumber == -1 || PlaneNumber == -1 || GlobalBarNumber == -1) {
    // Since the bar has already been created as a "error" in the above empty constructor we can just return
    //std::cout << "Bar number or plane number not found in Geometry" << std::endl;
    return false;
  }

  return true;
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
  while (NodeName.find(TMS_Const::TMS_ModuleLayerName) == std::string::npos && 
      NodeName.find(TMS_Const::TMS_TopLayerName) == std::string::npos) {
    nav->CdUp();
    NodeName = std::string(nav->GetCurrentNode()->GetName());
  }

  // If we've reached the world volume we don't have a scintillator hit -> return some mad bad value
  if (NodeName.find(TMS_Const::TMS_TopLayerName) != std::string::npos) {
    // Since the bar has already been created as a "error" in the above empty constructor we can just return
    std::cout << "Bar position not found in TMS_Bar::FindBar!" << std::endl;
    return -1;
  }

  int Number = nav->GetCurrentNode()->GetNumber();

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
