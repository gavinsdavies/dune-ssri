// python code was very slow...
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TGeoManager.h"
#include "TStopwatch.h"

#include "EDepSim/TG4Event.h"
#include "EDepSim/TG4PrimaryVertex.h"

const double mass_mu = 105.6583755;
const double extra_trk_length = 24.35;

bool dumpSSRITree(std::string filename) {
  std::cout << "Got " << filename << std::endl;

  // Offset offsets (all?) the vertex position in {x,y,z}, don't know why
  const double offset[] = { 0., 5.5, 411. };
  // Fiducial volume cuts? Never used
  const double fvLo[] = { -300., -100., 50. };
  const double fvHi[] = { 300., 100., 450. };
  // Never used
  const double collarLo[] = { -320., -120., 30. };
  const double collarHi[] = { 320., 120., 470. };

  const int nplanes_cut = 8;

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
  // For the true event we just care about the neutrino
  // GENIE saves max 100 particles in gRooTracker format
  const int kNPmax = 100;
  int True_NuPDG[kNPmax];
  double True_NuE[kNPmax][4];
  gRoo->SetBranchAddress("StdHepPdg", True_NuPDG);
  gRoo->SetBranchAddress("StdHepP4", True_NuE);

  int N_entries = events->GetEntries();
  int ient = 0;

  //////////////////////////////////////////////////
  // Setup the output tree
  std::string output_filename = filename.substr(0, filename.find(".root"));
  output_filename += "_output.root";
  TFile *output = new TFile(output_filename.c_str(), "RECREATE" );
  TTree *tree_out = new TTree("tree", "tree");

  std::string ReactionCode;
  int ievt;
  float E_nu;
  int PDG_nu;
  float Lepton[3];
  float Vertex[3];
  float LeptonDeath[3];
  float MuonDeath[3];
  float MuonBirth[3];
  float LepRMMS[3];
  int LepPDG;
  float LepKE;
  float RMMS_KE;
  float LepE;
  float MuonExitPt[3];
  float MuonExitMom[3];
  float MuonExitKE;
  int MuonReco;
  float MuonScintLen;
  float MuonArLen;
  float MuonScintLen_Old;
  float MuonScintEnergy;
  int nFinalState;
  const int MAX_PARTICLES = 100;
  int PDGFinalState[MAX_PARTICLES];
  float PxFinalState[MAX_PARTICLES];
  float PyFinalState[MAX_PARTICLES];
  float PzFinalState[MAX_PARTICLES];
  float EFinalState[MAX_PARTICLES];
  std::vector<float> xPos;
  std::vector<float> yPos;
  std::vector<float> zPos;

  tree_out->Branch("reac", &ReactionCode);
  tree_out->Branch("ievt", &ievt, "ievt/I");
  tree_out->Branch("Ev", &E_nu, "Ev/F");
  tree_out->Branch("PDGv", &PDG_nu, "PDGv/I");
  tree_out->Branch("p3lep", Lepton, "p3lep[3]/F");
  tree_out->Branch("vtx", Vertex, "vtx[3]/F");
  tree_out->Branch("lepDeath", LeptonDeath, "lepDeath[3]/F");
  tree_out->Branch("muonDeath", MuonDeath, "muonDeath[3]/F");
  tree_out->Branch("muonBirth", MuonBirth, "muonBirth[3]/F");
  tree_out->Branch("lepRMMS", LepRMMS, "lepRMMS[3]/F");
  tree_out->Branch("lepPdg", &LepPDG, "lepPDG/I");
  tree_out->Branch("lepKE", &LepKE, "lepKE/F");
  tree_out->Branch("rmmsKE", &RMMS_KE, "rmmsKE/F");
  tree_out->Branch("lepE", &LepE, "lepE/F");
  tree_out->Branch("muonExitPt", MuonExitPt, "muonExitPt[3]/F");
  tree_out->Branch("muonExitMom", MuonExitMom, "muonExitMom[3]/F");
  tree_out->Branch("muonExitKE", &MuonExitKE, "muonExitKE/F");
  tree_out->Branch("muonReco", &MuonReco, "muonReco/I");
  tree_out->Branch("muScintLen", &MuonScintLen, "muScintLen/F");
  tree_out->Branch("muArLen", &MuonArLen, "muArLen/F");
  tree_out->Branch("muScintLenOld", &MuonScintLen_Old, "muScintLenOld/F");
  tree_out->Branch("muScintEnergy", &MuonScintEnergy, "muScintEnergy/F");
  tree_out->Branch("nFS", &nFinalState, "nFS/I");
  tree_out->Branch("fsPdg", PDGFinalState, "fsPdg[nFS]/I");
  tree_out->Branch("fsPx", PxFinalState, "fsPx[nFS]/F");
  tree_out->Branch("fsPy", PyFinalState, "fsPy[nFS]/F");
  tree_out->Branch("fsPz", PzFinalState, "fsPz[nFS]/F");
  tree_out->Branch("fsE", EFinalState, "fsE[nFS]/F");
  tree_out->Branch("xpt", &xPos);
  tree_out->Branch("ypt", &yPos);
  tree_out->Branch("zpt", &zPos);
  // Finished the output branches
  //////////////////////////////////////////////////

  // Now loop
  std::cout << "Starting loop over " << N_entries << " entries..." << std::endl;
  TStopwatch Timer;
  Timer.Start();
  for (int i = 0; i < N_entries; ++i) {
    events->GetEntry(i);
    gRoo->GetEntry(i);

    E_nu = True_NuE[0][3];
    PDG_nu = True_NuPDG[0];

    if (i % (N_entries/10) == 0) {
      std::cout << "Processed " << i << "/" << N_entries << " (" << double(i)*100./N_entries << "%)" << std::endl;
      std::cout << "  E_nu: " << E_nu << " PDG: " << PDG_nu << std::endl;
    }

    // Loop over the primary vertices
    for (TG4PrimaryVertexContainer::iterator it = event->Primaries.begin(); it != event->Primaries.end(); ++it) {
      // The vertex
      TG4PrimaryVertex vtx = *it;

      // Now reset the variables
      ievt = LepPDG = MuonReco = -999;
      LepKE = RMMS_KE = LepE = MuonExitKE = MuonScintLen = MuonArLen = MuonScintLen_Old = MuonScintEnergy;
      nFinalState = 0;
      for (int j = 0; j < 3; ++j) {
        Lepton[j] = Vertex[j] = LeptonDeath[j] = MuonDeath[j] = MuonBirth[j] = LepRMMS[j] = MuonExitPt[j] = MuonExitMom[j] = -999;
      }
      for (int j = 0; j < MAX_PARTICLES; ++j) {
        PDGFinalState[j] = PxFinalState[j] = PyFinalState[j] = PzFinalState[j] = EFinalState[j] = -999;
      }
      xPos.clear();
      yPos.clear();
      zPos.clear();
      ReactionCode.clear();
      // Done resetting variables

      ReactionCode = vtx.GetReaction();
      ievt = i;
      // Technicall a TLorentzVector
      for (int j = 0; j < 3; ++j) Vertex[j] = (vtx.GetPosition()[j])/10. - offset[j]; // in cm

      int LepIndex = -1;
      std::vector<TG4PrimaryParticle> particles = vtx.Particles;
      // Loop over the particles in the vertex
      for (TG4PrimaryVertex::PrimaryParticles::iterator jt = particles.begin(); jt != particles.end(); ++jt) {
        TG4PrimaryParticle particle = *jt;
        // Get some fundamentals about the particle
        double E = particle.GetMomentum()[3];
        double p = particle.GetMomentum().Vect().Mag();
        double m = particle.GetMomentum().Mag();
        int pdg = particle.GetPDGCode();

        PDGFinalState[nFinalState] = pdg;
        PxFinalState[nFinalState] = particle.GetMomentum().Px();
        PyFinalState[nFinalState] = particle.GetMomentum().Py();
        PzFinalState[nFinalState] = particle.GetMomentum().Pz();
        EFinalState[nFinalState] = particle.GetMomentum().E();

        // Find the lepton with matching PDG to incoming neutrino
        // If NC, PDG matches
        if (PDG_nu == pdg) LepIndex = nFinalState;
        // If neutrino, look for charged lepton PDG of PDG_nu-1
        else if (PDG_nu > 0 && PDG_nu == pdg+1) LepIndex = nFinalState;
        // If anti-neutrino, look for charged lepton PDG of PDG_nu+1
        else if (PDG_nu < 0 && PDG_nu == pdg-1) LepIndex = nFinalState;

        // Increment final state particle counter
        nFinalState++;
      }

      if (LepIndex == -1) {
        std::cerr << "No charged lepton found in event " << ievt << std::endl;
      }

      // Now extract the Lepton kinematics
      TG4PrimaryParticle lepton = particles[LepIndex];
      for (int j = 0; j < 3; ++j) Lepton[j] = lepton.GetMomentum()[j];
      LepPDG = lepton.GetPDGCode();
      LepKE = lepton.GetMomentum().E()-lepton.GetMomentum().M();

      // If there is a muon, determine how to reconstruct its momentum and charge
      bool Exit = false;
      bool InTMS = false;

      // Get the trajectory of that lepton trajector
      TG4Trajectory LepTraj = event->Trajectories[LepIndex];
      // Loop over the points that the lepton trajectory left
      std::vector<TG4TrajectoryPoint> LepPoints = LepTraj.Points;
      // Save the previous point
      TLorentzVector PreviousPoint(-999, -999, -999, -999);
      for (TG4Trajectory::TrajectoryPoints::iterator jt = LepPoints.begin(); jt != LepPoints.end(); ++jt) {
        TG4TrajectoryPoint pt = *jt;
        TLorentzVector Point = pt.GetPosition();
        TGeoNode *node = geom->FindNode(Point.X(), Point.Y(), Point.Z());
        if (node == NULL) {
          std::cerr << "Couldn't find node for this point" << std::endl;
          continue;
        }
        std::string VolumeName = node->GetName();
        bool active = false;
        // Get the position of this point
        TVector3 Pos(Point.X()/10.-offset[0], Point.Y()/10.-offset[1], Point.Z()/10.-offset[2]);

        // Easiest to just do a string comparison, unfortunately
        // If in active LAr volume or pixels, update the exit points
        if (VolumeName.find("LAr") != std::string::npos || 
            VolumeName.find("PixelPlane") != std::string::npos || 
            VolumeName.find("sPlane") != std::string::npos) {
          MuonExitPt[0] = Pos.X();
          MuonExitPt[1] = Pos.Y();
          MuonExitPt[2] = Pos.Z();
          // What momentum did it have when it exited the LAr
          MuonExitMom[0] = pt.GetMomentum().Px();
          MuonExitMom[1] = pt.GetMomentum().Py();
          MuonExitMom[2] = pt.GetMomentum().Pz();
          // Assume it's a muon to get a KE estimate
          MuonExitKE = sqrt(pt.GetMomentum().Mag2()+mass_mu*mass_mu) - mass_mu;
        // If not in tracked LAr
        } else {

          // Update that trajectory has exited the LAr
          if (!Exit) {
            MuonExitPt[0] = Pos.X();
            MuonExitPt[1] = Pos.Y();
            MuonExitPt[2] = Pos.Z();
            Exit = true;
          }

          // Check if it's inside the TMS
          // Read what the kinetic energy is when entering the TMS
          if (!InTMS && (VolumeName.find("RMMS") != std::string::npos || 
                          VolumeName.find("modulelayer") != std::string::npos)) {
            RMMS_KE = sqrt(pt.GetMomentum().Mag2()+mass_mu*mass_mu) - mass_mu;
            InTMS = true;
          }
        }
        // Update the previous point
        PreviousPoint = Point;
      }

      // Get the last point to see where the trajectory went
      TLorentzVector EndPoint = LepPoints.back().GetPosition();
      TGeoNode *EndNode = geom->FindNode(EndPoint.X(), EndPoint.Y(), EndPoint.Z());
      if (EndNode == NULL) {
        std::cerr << "Found no end node for this event" << std::endl;
      } else {
        LeptonDeath[0] = EndPoint.X()/10.-offset[0];
        LeptonDeath[1] = EndPoint.Y()/10.-offset[1];
        LeptonDeath[2] = EndPoint.Z()/10.-offset[2];
        std::string EndVolName = EndNode->GetName();
        if (EndVolName.find("LArActive") != std::string::npos) MuonReco = 1;
        else if (EndVolName.find("RMMS") != std::string::npos) MuonReco = 2;
        else MuonReco = 0;
      }

      tree_out->Fill();
    }
  }
  Timer.Stop();
  std::cout << "Event loop took " << Timer.RealTime() << "s for " << N_entries << " (" << Timer.RealTime()/N_entries << " s/event)" << std::endl;

  std::cout << "Writing to " << output->GetName() << "..." << std::endl;
  output->cd();
  tree_out->Write();
  output->Close();

  return true;
}

int main(int argc, char **argv) {
  if (argc != 2) {
    std::cerr << "Need one argument: the EDepSim output file" << std::endl;
    return -1;
  }
  bool ok = dumpSSRITree(std::string(argv[1]));
  if (ok) return 0;
  else return -1;
}
