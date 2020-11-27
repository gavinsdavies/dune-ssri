#include "TMS_Bar.h"

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
  // Updates the x,y,z,xw,yw,zw
  FindModules(x, y, z);
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

  // cd up in the geometry to find the right name
  while (NodeName.find(TMS_Const::TMS_TopLayerName) == std::string::npos) {

    // We've found the plane number
    if (NodeName.find(TMS_Const::TMS_ModuleLayerName) != std::string::npos) {
      PlaneNumber = geom->GetCurrentNode()->GetNumber();
    }

    // This is the furthest down hit we have: scintillator bar
    else if (NodeName.find(TMS_Const::TMS_ScintLayerName) != std::string::npos) {
      BarNumber = geom->GetCurrentNode()->GetNumber();

      // Get the width
      TGeoBBox *box = dynamic_cast<TGeoBBox*>(geom->GetCurrentVolume()->GetShape());
      // ROOT saves half of the width of a shape
      xw = 2*box->GetDX();
      yw = 2*box->GetDY();
      zw = 2*box->GetDZ();

      // Do a sanity check (CHEATING!)
      // Know the bars are 1cm in z and 4cm in x
      if (zw != 10 || xw != 40) {
        std::cerr << "width of " << NodeName << " not as expected!" << std::endl;
        std::cerr << "xwidth: " << xw << std::endl;
        std::cerr << "zwidth: " << zw << std::endl;
        throw;
      }

      for (int i = 0; i < 3; ++i) Translation[i] = geom->GetCurrentMatrix()->GetTranslation()[i];
    }

    else if (NodeName.find(TMS_Const::TMS_ModuleName) != std::string::npos) {
      GlobalBarNumber = geom->GetCurrentNode()->GetNumber();
    }

    geom->CdUp();
    NodeName = std::string(geom->GetCurrentNode()->GetName());
  }

  // Update the hit value to be the bar, not the exact hit position
  x = Translation[0];
  // y-bar information doesn't seem to work?
  //y = Translation[1];
  z = Translation[2];

  // For the y-bar hack around it
  // Know y start and end points of detector, so split up into 40mm slices and see in which slice this falls in
  y = FindYbar(y); 

  // If we've reached the world volume we don't have a scintillator hit -> return some mad bad value
  if (BarNumber == -1 || PlaneNumber == -1 || GlobalBarNumber == -1) {
    // Since the bar has already been created as a "error" in the above empty constructor we can just return
    return false;
  }

  // Set the width and orientation of the bar
  if (PlaneNumber % 2 == 0) BarOrient = kXBar;
  else BarOrient = kYBar;

  // If this is a y-bar, remove the y coordinate
  if (BarOrient == kXBar) {
    x = -99999000;
    // Flip the widths
    //std::cout << "flipping width "  << yw << " to be " << xw << std::endl;
    double tempyw = yw;
    yw = xw;
    xw = tempyw; 
    //std::cout << "flipped width "  << yw << " to be " << xw << std::endl;
  } else if (BarOrient == kYBar) {
    //std::cout << "for ybar: " << yw << " " << xw << std::endl;
    y = -99999000;
    // Don't need to flip the widths because they're already correct (yw = large, xw = 4cm)
  } else {
    x = -99999000;
    y = -99999000;
    z = -99999000;
  }

  // Reset the geom navigator node level in case it's used again
  geom->FindNode(xval,yval,zval);

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

double TMS_Bar::FindYbar(double yval) {
  const double ymin = (-250+TMS_Const::TMS_Det_Offset[1])*10;
  //const double ymax = (100+TMS_Const::TMS_Det_Offset[1])*10;

  // The total range of y
  //double yrange = ymax-ymin;
  // Splits into how many 40mm slices
  //int nSlices = yrange/40;
  int bin = (yval-ymin)/40;
  // Return the center of the bin
  double val = ymin+bin*40+20;
  //std::cout << yval << " corrected to " << val << std::endl;

  return val;
}
