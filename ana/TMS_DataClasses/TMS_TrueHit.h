#ifndef _TMS_TRUEHIT_H_
#define _TMS_TRUEHIT_H_

#include <vector>

// Include the constants
#include "TMS_Constants.h"

#include "EDepSim/TG4HitSegment.h"

// Essentially a copy of the edep-sim THit
class TMS_TrueHit {
  public:
    TMS_TrueHit(TG4HitSegment &edep_seg);
    TMS_TrueHit(double x, double y, double z, double t, double E);

    TMS_TrueHit();
    ~TMS_TrueHit();

    double GetX() {return x;};
    double GetY() {return y;};
    double GetZ() {return z;};
    double GetT() {return t;};
    double GetE() {return EnergyDeposit; };
    const std::vector<int> &GetContributors() {return Contributors;};

    void SetX(double pos) {x = pos;};
    void SetY(double pos) {y = pos;};
    void SetZ(double pos) {z = pos;};
    void SetT(double pos) {t = pos;};
    void SetE(double E) {EnergyDeposit = E;};
    void AddContributor(int Con) { Contributors.push_back(Con); };

  private:
    double x;
    double y;
    double z;
    double t;
    double EnergyDeposit;
    std::vector<int> Contributors;
};

#endif
