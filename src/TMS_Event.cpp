#include "TMS_Event.h"

// Initialise the event counter to 0
int TMS_Event::EventNumber = 0;

// Start the relatively tedious process of converting into TMS products!
TMS_Event::TMS_Event(TG4Event &event) {
  bool OnlyMuon = false;
  bool LArOnly = true;
  bool TMSLArOnly = true;

  // Check the integrity of the event
  //CheckIntegrity();

  // Loop over the primary vertices
  for (TG4PrimaryVertexContainer::iterator it = event.Primaries.begin(); it != event.Primaries.end(); ++it) {

    TG4PrimaryVertex vtx = *it;
    /*
    std::vector<TG4PrimaryParticle> particles = vtx.Particles;

    // Loop over the particles in the vertex
    std::cout << particles.size() << " particles" << std::endl;
    for (TG4PrimaryVertex::PrimaryParticles::iterator jt = particles.begin(); jt != particles.end(); ++jt) {
      const TG4PrimaryParticle particle = *jt;
      std::cout << particle.GetName() << std::endl;
      std::cout << particle.GetPDGCode() << std::endl;
      std::cout << particle.GetTrackId() << std::endl;
      TMS_TrueParticle truepart = TMS_TrueParticle(particle);
      TMS_TrueParticles.emplace_back(truepart);
    }
    */

    bool skip = true;

    for (TG4TrajectoryContainer::iterator jt = event.Trajectories.begin(); jt != event.Trajectories.end(); ++jt) {
      TG4Trajectory traj = *jt;

      // Only bother to save the muon for now
      int PDGcode = traj.GetPDGCode();
      if (OnlyMuon && abs(PDGcode) != 13) continue;

      // Only from fundamental vertex
      int ParentId = traj.GetParentId();
      if (ParentId != -1) continue;
      int TrackId = traj.GetTrackId();

      // This could be upstream in the LAr or anywhere the muon is born
      // The starting point of the trajectory
      //TVector3 Momentum = traj.Points[0].GetMomentum(); 
      //TLorentzVector Position = traj.Points[0].GetPosition();

      // Let's identify the detector the first point is in; this way can isolate on LAr and TMS
      TG4TrajectoryPoint start = traj.Points[0];
      TGeoNode *vol = TMS_Geom::GetInstance().GetGeometry()->FindNode(start.GetPosition().X(), start.GetPosition().Y(), start.GetPosition().Z());

      std::string VolumeName = vol->GetName();
      bool InLAr = false;
      bool InTMS = false;
      if (VolumeName.find(TMS_Const::LAr_ActiveName) != std::string::npos) InLAr = true;
      if (VolumeName.find(TMS_Const::TMS_EDepSim_VolumeName) != std::string::npos) InTMS = true;
      if (LArOnly && !InLAr) continue;
      if (TMSLArOnly && !InLAr && !InTMS) continue;
      skip = false;

      // Save down the trajectory points of the true particles
      for (std::vector<TG4TrajectoryPoint>::iterator kt = traj.Points.begin(); kt != traj.Points.end(); kt++) {
        TG4TrajectoryPoint pt = *kt;
        // Check the point against the geometry
        TGeoNode *vol = TMS_Geom::GetInstance().GetGeometry()->FindNode(pt.GetPosition().X(), pt.GetPosition().Y(), pt.GetPosition().Z());
        VolumeName = vol->GetName();
        // Only look at TMS hits
        if (VolumeName.find(TMS_Const::TMS_VolumeName) == std::string::npos && 
            VolumeName.find(TMS_Const::TMS_ModuleLayerName) == std::string::npos) continue;

        TLorentzVector Position = pt.GetPosition();
        TVector3 Momentum = pt.GetMomentum();

        // Make the true particle that created this trajectory
        TMS_TrueParticle muon(ParentId, TrackId, PDGcode, Momentum, Position);

        TMS_TrueParticles.push_back(muon);
        // We have the TMS starting point, now exit the loop
        break;
      }
    }

    if (!skip) {
      // Loop over each hit
      for (TG4HitSegmentDetectors::iterator jt = event.SegmentDetectors.begin(); jt != event.SegmentDetectors.end(); ++jt) {
        // Only look at TMS hits
        std::string DetString = (*jt).first;
        if (DetString != TMS_Const::TMS_EDepSim_VolumeName) continue;

        TG4HitSegmentContainer tms_hits = (*jt).second;
        for (TG4HitSegmentContainer::iterator kt = tms_hits.begin(); kt != tms_hits.end(); ++kt) {
          TG4HitSegment edep_hit = *kt;
          TMS_Hit hit = TMS_Hit(edep_hit);
          TMS_Hits.push_back(hit);

          // Now associate the hits with the muon
          int PrimaryId = edep_hit.GetPrimaryId();
          // Loop through the True Particle list and associate
          for (auto &TrueParticle : TMS_TrueParticles) {
            // Check the primary ID
            if (TrueParticle.GetTrackId() != PrimaryId) continue;
            TLorentzVector Position = (edep_hit.GetStop()+edep_hit.GetStart());
            Position *= 0.5;
            TrueParticle.AddPoint(Position);
          }
        }
      }
    }
  }

  EventNumber++;
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

