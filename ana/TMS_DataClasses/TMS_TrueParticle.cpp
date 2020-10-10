#include "TMS_TrueParticle.h"

TMS_TrueParticle(TG4PrimaryParticle &particle) {
  FourVector = particle.GetMomentum();
  PDG = particle.GetPDGCode();
  TrackId = particle.GetTrackId();
}

TMS_TrueParticle() :
  PDG(-999);
  FourVector(-999.99, -999.99, -999.99, -999.99);
  Parent(-999);
  TrackId(-999);
}

// Set the true particle from a segment
// Can't set from TG4HitSegment because it does not know about the trajectories
// Needs to be set from the event?
/*
TMS_TrueParticle(TG4HitSegment &segment) {

  // Get the primary particle from the segment
  int PrimaryId = segment.Contrib[0];

  FourVector = particle.GetMomentum();
  PDG = particle.GetPDGCode();
  TrackId = particle.GetTrackId();
}
*/
