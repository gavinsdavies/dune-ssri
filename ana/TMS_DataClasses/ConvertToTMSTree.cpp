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
// Geometry singleton
#include "TMS_Geom.h"
// Event class
#include "TMS_Event.h"
// Event viewer singleton
#include "TMS_EventViewer.h"

/*
#include "TMS_Constants.h"
#include "TMS_Hit.h"
#include "TMS_Bar.h"
#include "TMS_TrueHit.h"
#include "TMS_TrueParticle.h"
*/

bool dumpSSRITree(std::string filename, std::string output_filename) {
  std::cout << "Got " << filename << ", writing to " << output_filename << std::endl;

  // The input file
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
  TMS_Geom::GetInstance().SetFileName(filename);

  int N_entries = events->GetEntries();

  std::cout << "Starting loop over " << N_entries << " entries..." << std::endl;
  TStopwatch Timer;
  Timer.Start();

  for (int i = 0; i < N_entries; ++i) {
    events->GetEntry(i);
    gRoo->GetEntry(i);

    if (i > 100) break;
    if (i % (N_entries/10) == 0) {
      std::cout << "Processed " << i << "/" << N_entries << " (" << double(i)*100./N_entries << "%)" << std::endl;
    }

    // Make a TMS event
    TMS_Event tms_event = TMS_Event(*event);
    // Dump information
    //tms_event.Print();

    // View it
    TMS_EventViewer::GetViewer().Draw(tms_event);

  } // End loop over all the events

  Timer.Stop();
  std::cout << "Event loop took " << Timer.RealTime() << "s for " << N_entries << " entries (" << Timer.RealTime()/N_entries << " s/entries)" << std::endl;

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
