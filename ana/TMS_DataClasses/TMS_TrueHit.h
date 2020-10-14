#ifndef _TMS_TRUEHIT_H_
#define _TMS_TRUEHIT_H_

#include <vector>
#include <iostream>

// Include the constants
#include "TMS_Constants.h"

#include "EDepSim/TG4HitSegment.h"

// Essentially a copy of the edep-sim THit
class TMS_TrueHit {
  public:
    TMS_TrueHit(TG4HitSegment &edep_seg);

    double GetX() {return x;};
    double GetY() {return y;};
    double GetZ() {return z;};
    double GetT() {return t;};
    double GetE() {return EnergyDeposit; };

    void SetX(double pos) {x = pos;};
    void SetY(double pos) {y = pos;};
    void SetZ(double pos) {z = pos;};
    void SetT(double pos) {t = pos;};
    void SetE(double E) {EnergyDeposit = E;};

    void Print();

  private:
    double x;
    double y;
    double z;
    double t;
    double EnergyDeposit;
    int PrimaryId;
};

#endif
