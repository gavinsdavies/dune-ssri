#ifndef _TMSCONSTANTS_H_SEEN_
#define _TMSCONSTANTS_H_SEEN_

#include <string>

// Constants
namespace TMS_KinConst {
  const double mass_mu = 105.6583755; // Muon mass in MeV/c2
}

// General constants for the TMS
namespace TMS_Const {

  // Dead region (area between LAr and TMS) track length contribution, in g/cm2
  const double Dead_Region_Track_Length = 24.35;

  // Material densities, in g/cm3
  const double LAr_density = 1.3954;
  const double TMS_Steel_density = 7.85;
  const double TMS_Scint_density = 1.05;

  // Z positions of the scintillator bars (cm)
  const double TMS_Thin_Start = 730.8;
  const double TMS_Trans_Start = 939.8;
  const double TMS_Thick_Start = 949.3;

  // Approximate starting and end positions of TMS detector in geometry for plotting hits, in {x,y,z}
  // in mm by default!
  const double TMS_Start[] = {-4000, -3500, 10000};
  const double TMS_End[] = {4000, 500, 20000};

  // Gap for TMS region that is thin iron layer (cm)
  const double TMS_Thin_gap = 5.5;
  // Gap for TMS region that is between thin and thick regions (cm)
  const double TMS_Transition_gap = 9.5;
  // Gap for TMS region that is thick iron layer (cm)
  const double TMS_Thick_gap = 8.0;

  // TMS scintillator width (1 cm)
  const double TMS_Scint_Width = 1;
  //TMS steel width in thin region (1.5 cm);
  const double TMS_Thin_Steel_Width = 1.5;
  //TMS steel width in thick region (4.0 cm);
  const double TMS_Thick_Steel_Width = 4.0;

  // Offsets to put the TMS in the middle
  const double TMS_Det_Offset[] = { 0., 0., 0. };
  //const double TMS_Det_Offset[] = { 0., 5.5, 411. };
  //const double TMS_Det_Offset[] = { 0., 0.0, 411. };

  // Needs translating by the TMS_Const::TMS_Det_Offset array
  // Start and end of top dead region in x
  const double TMS_Dead_Top[] = {171.7, 180.4};
  // Start and end of central dead region in x
  const double TMS_Dead_Center[] = {-3.3, 3.3};
  // Start and end of bottom dead region in x
  const double TMS_Dead_Bottom[] = {-180.4, -171.7};

  // Some distance after the end of the TMS
  const double TMS_End_z = (1500+TMS_Const::TMS_Det_Offset[2])*10;

  // Start and end of top dead region in x
  const double TMS_Dead_Top_T[] = {(TMS_Dead_Top[0]+TMS_Det_Offset[1])*10, (TMS_Dead_Top[1]+TMS_Det_Offset[1])*10};
  // Start and end of central dead region in x
  const double TMS_Dead_Center_T[] = {(TMS_Dead_Center[0]+TMS_Det_Offset[1])*10, (TMS_Dead_Center[0]+TMS_Det_Offset[1])*10};
  // Start and end of bottom dead region in x
  const double TMS_Dead_Bottom_T[] = {(TMS_Dead_Bottom[0]+TMS_Det_Offset[1])*10, (TMS_Dead_Bottom[1]+TMS_Det_Offset[1])*10};

  // Volume name of TMS related hits
  const std::string TMS_VolumeName = "TMS";
  // Volume name for edep-sim SegmentDetectors
  const std::string TMS_EDepSim_VolumeName = "volTMS";
  // To find in z
  const std::string TMS_ModuleLayerName = "modulelayervol_PV";
  // To find scintillator "box"
  const std::string TMS_ModuleName = "ModuleBoxvol_PV";
  // To find scintillator "box"
  const std::string TMS_ScintLayerName = "scinBoxlvTMS_PV";
  // The top layer name
  const std::string TMS_TopLayerName = "volWorld";

  const std::string LAr_ActiveName = "volTPCActive";
}

#endif
