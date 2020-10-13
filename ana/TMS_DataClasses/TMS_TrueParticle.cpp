#include "TMS_TrueParticle.h"

TMS_TrueParticle::TMS_TrueParticle(const TG4PrimaryParticle &particle) {
  FourVector = particle.GetMomentum();
  PDG = particle.GetPDGCode();
  TrackId = particle.GetTrackId();
  Parent = -999;
}

TMS_TrueParticle::TMS_TrueParticle() :
  PDG(-999),
  FourVector(-999.99, -999.99, -999.99, -999.99),
  Parent(-999),
  TrackId(-999)
{
}

// Set the true particle from a segment
// Can't set from TG4HitSegment because it does not know about the trajectories
// Needs to be set from the event?
/*
TMS_TrueParticle::TMS_TrueParticle(const TG4HitSegment &segment) {

  // Get the primary particle from the segment
  int PrimaryId = segment.Contrib[0];

  FourVector = particle.GetMomentum();
  PDG = particle.GetPDGCode();
  TrackId = particle.GetTrackId();
}
*/

void TMS_TrueParticle::Print() {
  std::cout << "Printing TMS_TrueParticle class: " << std::endl;
  std::cout << "  PDG = " << PDG << std::endl;
  std::cout << "  (Px, Py, Pz, E) = " << " (" << FourVector.Px() << ", " << FourVector.Py() << ", " << FourVector.Pz() << ", " << FourVector.E() << ")" << std::endl;
  std::cout << "  Parent: " << Parent << std::endl;
  std::cout << "  TrackId: " << TrackId << std::endl;
}
