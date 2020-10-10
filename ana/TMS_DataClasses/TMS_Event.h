#ifndef _TMS_EVENT_H_SEEN_
#define _TMS_EVENT_H_SEEN_

// Include the constants
#include "TMS_Constants.h"

// The general event class
class TMS_Event {
  public:
    TMS_Event();
    ~TMS_Event();

    std::vector<TMS_Hit> GetHits() {return TMS_Hits;};
    std::vector<TMS_Track> GetTracks() {return TMS_Tracks;};
  private:
    // Hits
    std::vector<TMS_Hit> TMS_Hits;
    // Reconstructed tracks
    std::vector<TMS_Track> TMS_Tracks;
    // Spill number (can have many events in a spill)?
    int SpillNumber;
}

#endif
