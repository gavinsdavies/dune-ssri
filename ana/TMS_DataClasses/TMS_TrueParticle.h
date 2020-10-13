#ifndef _TMS_TRUEPARTICLE_H_SEEN_
#define _TMS_TRUEPARTICLE_H_SEEN_

// edep-sim classes for TG4PrimaryParticle
#include "EDepSim/TG4PrimaryVertex.h"
#include "EDepSim/TG4HitSegment.h"

#include "TLorentzVector.h"

#include <iostream>

class TMS_TrueParticle {
  public:

    TMS_TrueParticle();
    // Construct directly from edep-sim
    TMS_TrueParticle(const TG4PrimaryParticle &particle);
    //TMS_TrueParticle(const TG4HitSegment &hit);
    //~TMS_TrueParticle();

    // Give four vector setters
    void SetFourVector(TLorentzVector Vec) { FourVector = Vec; };
    void SetPx(double px) { FourVector.SetPx(px); };
    void SetPy(double py) { FourVector.SetPy(py); };
    void SetPz(double pz) { FourVector.SetPz(pz); };
    void SetE(double E) { FourVector.SetE(E); };

    void SetParent(int num) { Parent = num; };
    void SetPDG(int num) { PDG = num; };
    void Print();

  private:
    int PDG;
    TLorentzVector FourVector;
    int Parent;
    int TrackId;
};

#endif
