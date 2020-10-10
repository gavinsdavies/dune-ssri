#include "TMS_Bar.h"

// To get the geometry
#include "TMS_Geom.h"

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
}

TMS_Bar::TMS_Bar(double xpos, double ypos, double zpos, int BarNo) {
  x = xpos;
  y = ypos;
  z = zpos;
  BarNumber = BarNo;
}

TMS_Bar::TMS_Bar(TG4HitSegment &edep_seg) {
}

