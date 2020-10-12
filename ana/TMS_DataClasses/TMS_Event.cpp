#include "TMS_Event.h"

// Start the relatively tedious process of converting into TMS products!
TMS_Event::TMS_Event(TG4Event &event) {

  // Check the integrity of the event
  //CheckIntegrity();

  // Loop over the primary vertices
  for (TG4PrimaryVertexContainer::iterator it = event.Primaries.begin(); it != event.Primaries.end(); ++it) {

    TG4PrimaryVertex vtx = *it;
    // Get the particle trajectories
    std::vector<TG4PrimaryParticle> particles = vtx.Particles;
    // Loop over the particles in the vertex
    for (TG4PrimaryVertex::PrimaryParticles::iterator jt = particles.begin(); jt != particles.end(); ++jt) {
      TG4PrimaryParticle particle = *jt;
      TMS_TrueParticle truepart = TMS_TrueParticle(particle);
      TMS_TrueParticles.push_back(truepart);
    }

    // Do we need the tracjectory points?
    // Maybe?
    for (TG4TrajectoryContainer::iterator jt = event.Trajectories.begin(); jt != event.Trajectories.end(); ++jt) {
      TG4Trajectory traj = *jt;

      // traj.GetTrackId()
      // traj.GetParentId()
      // traj.GetName()
      // traj.GetPDGCode()
      // traj.GetInitialMomentum
      // traj.Points
    // Then loop over the points
      for (std::vector<TG4TrajectoryPoint>::iterator kt = traj.Points.begin(); kt != traj.Points.end(); kt++) {
        TG4TrajectoryPoint pt = *kt;
        // pt.GetProcess()
        // pt.GetSubprocess()
        // pt.GetMomentum()
        // pt.GetPosition();
      }
    }

    // Loop over each hit
    for (TG4HitSegmentDetectors::iterator jt = event.SegmentDetectors.begin(); jt != event.SegmentDetectors.end(); ++jt) {
      // Only look at TMS hits
      std::string DetString = (*jt).first;
      if (DetString != "rmmsvol") continue;

      TG4HitSegmentContainer tms_hits = (*jt).second;
      for (TG4HitSegmentContainer::iterator kt = tms_hits.begin(); kt != tms_hits.end(); ++kt) {
        TG4HitSegment edep_hit = *kt;
        TMS_Hit hit = TMS_Hit(edep_hit);
        TMS_Hits.push_back(hit);
      }
    }

  }
}
