#ifndef _TMS_TRUEPARTICLE_H_SEEN_
#define _TMS_TRUEPARTICLE_H_SEEN_

// edep-sim classes for TG4PrimaryParticle
#include "EDepSim/TG4PrimaryVertex.h"
#include "EDepSim/TG4HitSegment.h"

#include "TLorentzVector.h"

#include <iostream>

class TMS_TrueParticle {
  public:

    TMS_TrueParticle() :
      Parent(-999), TrackId(-999), PDG(-999), TMS_Momentum(TVector3(-999, -999, -999)) {
    }

    // Construct directly from edep-sim
    TMS_TrueParticle(int ParentVal, int TrackVal, int PDGVal, TVector3 Momentum, TLorentzVector Position) : 
      Parent(ParentVal), TrackId(TrackVal), PDG(PDGVal), TMS_Momentum(Momentum), TMS_EntryPoint(Position)
  {
    }

    void AddPoint(TLorentzVector &Position, TVector3 &Momentum) {
      PositionPoints.push_back(Position);
      MomentumPoints.push_back(Momentum);
    }

    void AddPoint(TLorentzVector &Position) {
      PositionPoints.push_back(Position);
    }

    void SetParent(int num) { Parent = num; };
    void SetTrackId(int num) { TrackId = num; };
    void SetPDG(int num) { PDG = num; };
    void Print();

    TVector3 GetInitialTMSMomentum() { return TMS_Momentum; };
    TLorentzVector GetInitialTMSPoint() { return TMS_EntryPoint; };
    int GetPDG() { return PDG; };
    int GetParent() { return Parent; };
    int GetTrackId() { return TrackId; };

    std::vector<TLorentzVector> &GetPositionPoints() { return PositionPoints; };
    std::vector<TVector3> &GetMomentumPoints() { return MomentumPoints; };


  private:
    int Parent;
    int TrackId;
    int PDG;
    TVector3 TMS_Momentum;
    TLorentzVector TMS_EntryPoint;
    std::vector<TLorentzVector> PositionPoints;
    std::vector<TVector3> MomentumPoints;
};

#endif
