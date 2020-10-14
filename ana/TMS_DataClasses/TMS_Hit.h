#ifndef _TMS_HIT_H_SEEN_
#define _TMS_HIT_H_SEEN_

#include <string>

// Include the constants
#include "TMS_Constants.h"
#include "TMS_Bar.h"
#include "TMS_TrueParticle.h"
#include "TMS_TrueHit.h"
// To get the geometry
#include "TMS_Geom.h"

// Include the edep-sim
#include "EDepSim/TG4HitSegment.h"

// Maybe should be in singleton?
#include "TGeoManager.h"
#include "TGeoNavigator.h"

// A low-level hit
class TMS_Hit {

  public:
    void Print();
    // The constructor for the TMS hit
    TMS_Hit(TG4HitSegment &edep_seg);
    //TMS_Hit(double x, double y, double z, double t, double E);
    //TMS_Hit();
    //~TMS_Hit();

    const TMS_Bar &GetBar() { return Bar; };
    void SetBar(TMS_Bar bar) { Bar = bar; };

    // The true particle that created this hit
    const TMS_TrueParticle &GetTrueParticle();

    // The true hit
    const TMS_TrueHit &GetTrueHit();

    // Over-riders (maybe delete in future)
    void SetTrueParticle(TMS_TrueParticle part) {TrueParticle = part;};
    void SetTrueHit(TMS_TrueHit hit) {TrueHit = hit;};

    void SetE(double E) {EnergyDeposit = E;};
    void SetT(double t) {Time = t;};

    double GetE() {return EnergyDeposit;};
    double GetT() {return Time;};

  private:
    // The true hit (x,y,z,t) --- does not quantise hit into bars
    TMS_TrueHit TrueHit;
    // The true particle that created this hit
    TMS_TrueParticle TrueParticle;
    // The bar that registered the hit
    TMS_Bar Bar;
    // The energy deposited
    double EnergyDeposit;
    // The timing of the hit
    double Time;
};

#endif
