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
      const TG4PrimaryParticle particle = *jt;
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
      if (DetString != TMS_Const::TMS_VolumeName) continue;

      TG4HitSegmentContainer tms_hits = (*jt).second;
      for (TG4HitSegmentContainer::iterator kt = tms_hits.begin(); kt != tms_hits.end(); ++kt) {
        TG4HitSegment edep_hit = *kt;
        TMS_Hit hit = TMS_Hit(edep_hit);
        TMS_Hits.push_back(hit);
      }
    }

  }
}

void TMS_Event::Print() {
  std::cout << std::endl;
  std::cout << "*** " << std::endl;
  std::cout << "Printing TMS_Event class from "  << __FILE__ << std::endl;
  std::cout << "  Using geometry: " << TMS_Geom::GetInstance().GetGeometry()->GetName() << ", " << TMS_Geom::GetInstance().GetGeometry()->GetTitle() << std::endl;
  std::cout << "  From: " << TMS_Geom::GetInstance().GetFileName() << std::endl;
  std::cout << "  N Truth particles: " << TMS_TrueParticles.size() << std::endl;
  std::cout << "  N Hits: " << TMS_Hits.size() << std::endl;
  std::cout << "  IsEmpty: " << IsEmpty() << std::endl;

  int PartCount = 0;
  std::cout << "Printing particle stack: " << std::endl;
  for (std::vector<TMS_TrueParticle>::iterator it = TMS_TrueParticles.begin(); it != TMS_TrueParticles.end(); ++it, ++PartCount) {
    std::cout << "Particle " << PartCount << std::endl;
    (*it).Print();
  }

  int HitCount = 0;
  std::cout << "Printing hit stack: " << std::endl;
  for (std::vector<TMS_Hit>::iterator it = TMS_Hits.begin(); it != TMS_Hits.end(); ++it, ++HitCount) {
    std::cout << "Hit "  << HitCount << std::endl;
    (*it).Print();
  }
}

