#ifndef _TMS_EVENT_H_SEEN_
#define _TMS_EVENT_H_SEEN_

// Include the constants
#include "TMS_Constants.h"
#include "TMS_Hit.h"
#include "TMS_Particle.h"

#include "TG4Event.h"

#include <string>
#include <iostream>

// The general event class
class TMS_Event {
  public:
    TMS_Event(TG4Event &event);
    ~TMS_Event();

    // The getters once the class is completed
    std::vector<TMS_Hit> GetHits() {return TMS_Hits;};
    // Reconstructed tracks
    std::vector<TMS_Track> GetTracks() {return TMS_Tracks;};
    // The true particles
    std::vector<TMS_TrueParticle> GetParticles() { return TMS_TrueParticles; };

  private:
    // Hits
    std::vector<TMS_Hit> TMS_Hits;
    // Reconstructed tracks
    std::vector<TMS_Track> TMS_Tracks;
    // True particles
    std::vector<TMS_TrueParticle> TMS_TrueParticles;

    // Spill number (can have many events in a spill)?
    int SpillNumber;
}

#endif
