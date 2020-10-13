#ifndef _TMS_BAR_H_SEEN_
#define _TMS_BAR_H_SEEN_

// Include the constants
#include "TMS_Constants.h"
// To get the geometry
#include "TMS_Geom.h"

#include "EDepSim/TG4HitSegment.h"

class TMS_Bar {
  public:

    TMS_Bar();
    TMS_Bar(double xpos, double ypos, double zpos);
    TMS_Bar(double xpos, double ypos, double zpos, int BarNo);
    TMS_Bar(TG4HitSegment &edep_seg);
    //~TMS_Bar();

    // Enum for the x, y, U, V bar orientation
    enum BarType { kXBar, kYBar, kUBar, kVBar, kError };
    std::string BarType_ToString(BarType bar);

    // Getter functions
    int GetBarNumber() { return BarNumber; };
    BarType GetBarType() { return BarOrient; };
    double GetX() { return x; };
    double GetY() { return y; };
    double GetZ() { return z; };

    double GetXw() { return xw; };
    double GetYw() { return yw; };
    double GetZw() { return zw; };

    void Print();

    int FindBar(double x, double y, double z);

    bool FindModules(double x, double y, double z);

  private:
    // Plane that the bar belongs in
    int PlaneNumber;
    // The bar number in this plane
    int BarNumber;
    // The global bar number (0-100) 
    int GlobalBarNumber;
    // The bar start positions
    double x;
    double y;
    double z;
    // The bar widths
    double xw;
    double yw;
    double zw;
    // Which type of bar
    BarType BarOrient;
};

#endif
