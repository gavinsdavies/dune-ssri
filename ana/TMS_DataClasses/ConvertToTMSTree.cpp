// python code was very slow...
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TGeoManager.h"
#include "TStopwatch.h"

// EDepSim includes
#include "EDepSim/TG4Event.h"
#include "EDepSim/TG4PrimaryVertex.h"

// TMS includes
#include "TMS_Geom.h"
#include "TMS_Event.h"
/*
#include "TMS_Constants.h"
#include "TMS_Hit.h"
#include "TMS_Bar.h"
#include "TMS_TrueHit.h"
#include "TMS_TrueParticle.h"
*/

bool dumpSSRITree(std::string filename, std::string output_filename) {
  std::cout << "Got " << filename << ", writing to " << output_filename << std::endl;

  TFile *input = new TFile(filename.c_str());
  // The EDepSim events
  TTree *events = (TTree*)input->Get("EDepSimEvents");
  // The generator pass-through information
  TTree *gRoo = (TTree*)input->Get("DetSimPassThru/gRooTracker");
  // Get the detector geometry
  TGeoManager *geom = (TGeoManager*)input->Get("EDepSimGeometry");

  // Get the event
  TG4Event *event = NULL;
  events->SetBranchAddress("Event", &event);

  // Load up the geometry
  TMS_Geom::GetInstance().SetGeometry(geom);

  int N_entries = events->GetEntries();

  std::cout << "Starting loop over " << N_entries << " entries..." << std::endl;
  TStopwatch Timer;
  Timer.Start();

  for (int i = 0; i < N_entries; ++i) {
    events->GetEntry(i);
    gRoo->GetEntry(i);

    if (i % (N_entries/10) == 0) {
      std::cout << "Processed " << i << "/" << N_entries << " (" << double(i)*100./N_entries << "%)" << std::endl;
    }

    TMS_Event tms_event = TMS_Event(*event);

    /*
    // Loop over the primary vertices
    for (TG4PrimaryVertexContainer::iterator it = event->Primaries.begin(); it != event->Primaries.end(); ++it) {

      // The vertex
      TG4PrimaryVertex vtx = *it;

      // Loop over the particles in the vertex
      std::vector<TG4PrimaryParticle> particles = vtx.Particles;
      for (TG4PrimaryVertex::PrimaryParticles::iterator jt = particles.begin(); jt != particles.end(); ++jt) {
        TG4PrimaryParticle particle = *jt;
      }

      TG4Trajectory LepTraj = event->Trajectories[LepIndex];
      // Loop over the points that the lepton trajectory left
      std::vector<TG4TrajectoryPoint> LepPoints = LepTraj.Points;
      for (TG4Trajectory::TrajectoryPoints::iterator jt = LepPoints.begin(); jt != LepPoints.end(); ++jt) {
        TG4TrajectoryPoint pt = *jt;
        TGeoNode *node = geom->FindNode(Point.X(), Point.Y(), Point.Z());
        if (node == NULL) {
          std::cerr << "Couldn't find node for this point" << std::endl;
          continue;
        }
        std::string VolumeName = node->GetName();
        bool active = false;
      }

      for (TG4HitSegmentDetectors::iterator jt = event->SegmentDetectors.begin(); jt != event->SegmentDetectors.end(); ++jt) {

        // Get the name of the active detector
        std::string DetString = (*jt).first;
//
        // Remember what this segment is in
        bool ThisSegInLAr = false;
        bool ThisSegInTMS = false;

        if (DetString == "rmmsvol") {
          ThisSegInTMS = true;
        }

        // The container for all of the segments
        TG4HitSegmentContainer hitsegments = (*jt).second;

        for (TG4HitSegmentContainer::iterator kt = hitsegments.begin(); kt != hitsegments.end(); ++kt) {
          // The current segment
          TG4HitSegment seg = (*kt);

          // Loop over the contributors to this track segment
          for (TG4HitSegment::Contributors::iterator ct = seg.Contrib.begin(); ct != seg.Contrib.end(); ++ct) {
            int Contributor = *ct;
            int Contributor = 0;
             Get the trajectory that produced this hit
            TG4Trajectory Trajectory = event->Trajectories[Contributor];
            int ParentId = Trajectory.GetParentId();
            int pdg = Trajectory.GetPDGCode();
            int TrackId = Trajectory.GetTrackId();
          }
        }

      } 
    } // End loop over the vertices in the event
    */
  } // End loop over all the events

  Timer.Stop();
  std::cout << "Event loop took " << Timer.RealTime() << "s for " << N_entries << " (" << Timer.RealTime()/N_entries << " s/event)" << std::endl;

  return true;
}

int main(int argc, char **argv) {
  if (argc != 2 && argc !=3) {
    std::cerr << "Need one argument: [EDepSim output file] [Output filename]" << std::endl;
    return -1;
  }

  std::string EDepSimFile = std::string(argv[1]);

  std::string OutputFile;
  // If only two arguments are given, 
  if (argc == 2) {
    std::string filename = std::string(argv[1]);
    OutputFile = filename.substr(0, filename.find(".root"));
    OutputFile += "_output.root";
  } else {
    OutputFile = std::string(argv[2]);
  }

  bool ok = dumpSSRITree(EDepSimFile, OutputFile);
  if (ok) return 0;
  else return -1;
}
