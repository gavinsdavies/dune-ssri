#ifndef _TMS_BAR_H_SEEN_
#define _TMS_BAR_H_SEEN_

// Include the constants
#include "TMS_Constants.h"
// To get the geometry
#include "TMS_Geom.h"

#include "EDepSim/TG4HitSegment.h"

#include "TGeoBBox.h"

class TMS_Bar {
  public:

    /*
    TMS_Bar();
    TMS_Bar(double xpos, double ypos, double zpos);
    TMS_Bar(double xpos, double ypos, double zpos, int BarNo);
    ~TMS_Bar();
    */
    TMS_Bar(TG4HitSegment &edep_seg);

    // Enum for the x, y, U, V bar orientation
    enum BarType { kXBar, kYBar, kUBar, kVBar, kError };
    std::string BarType_ToString(BarType bar);

    // Getter functions
    int GetBarNumber() { return BarNumber; };
    int GetPlaneNumber() { return PlaneNumber; };
    int GetGlobalBarNumber() { return GlobalBarNumber; };
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

    // Find if a 2D point is inside the bar
    // x here denotes the other view than z
    // can be both x and y views (depending on bar type)!
    bool Contains(double x, double z) {

      //std::cout << "BarOrient: " << BarType_ToString(BarOrient) << std::endl;
      //std::cout << "Getting (x,z)=(" << x << "," << z << ")" << std::endl;

      double zmin = GetZ()-GetZw()/2;
      double zmax = GetZ()+GetZw()/2;
      //std::cout << "zmin,zmax: " << zmin << "," << zmax << std::endl;
      if (z > zmax || z < zmin) return false;

      double xmin = -9999999999999;
      double xmax = 9999999999999;
      if (BarOrient == kXBar) {
        xmin = GetY()-GetYw()/2;
        xmax = GetY()+GetYw()/2;
      } else if (BarOrient == kYBar) {
        xmin = GetX()-GetXw()/2;
        xmax = GetX()+GetXw()/2;
      }
      //std::cout << "xmin,xmax: " << xmin << "," << xmax << std::endl;

      if (x > xmax || x < xmin) return false;

      return true;
    }

  private:
    // Plane that the bar belongs in
    int PlaneNumber;
    // The bar number in this plane
    int BarNumber;
    // The global bar number (0-100) 
    int GlobalBarNumber;
    // All in mm units!
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
