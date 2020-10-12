#ifndef _TMS_EVENT_H_SEEN_
#define _TMS_EVENT_H_SEEN_

#include <string>
#include <iostream>

// Include the constants
#include "TMS_Constants.h"
#include "TMS_Hit.h"
#include "TMS_TrueParticle.h"
//#include "TMS_Tracks.h"

// The edep-sim event class
#include "EDepSim/TG4Event.h"

// The general event class
class TMS_Event {
  public:
    TMS_Event(TG4Event &event);
    ~TMS_Event();

    // The getters once the class is completed
    std::vector<TMS_Hit> GetHits() {return TMS_Hits;};
    // Reconstructed tracks
    //std::vector<TMS_Track> GetTracks() {return TMS_Tracks;};
    // The true particles
    std::vector<TMS_TrueParticle> GetParticles() { return TMS_TrueParticles; };

  private:
    // Hits
    std::vector<TMS_Hit> TMS_Hits;
    // Reconstructed tracks
    //std::vector<TMS_Track> TMS_Tracks;
    // True particles
    std::vector<TMS_TrueParticle> TMS_TrueParticles;

    // Spill number (can have many events in a spill)?
    int SpillNumber;
};

#endif
