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
  const double offset[] = { 0., 5.5, 411. };

  // Volume name of TMS related hits
  const std::string TMS_VolumeName = "rmmsvol";
  // To find in z
  const std::string TMS_ModuleLayerName = "modulelayervol_PV";
  // To find scintillator "box"
  const std::string TMS_ModuleName = "ModuleBoxvol_PV";
  // To find scintillator "box"
  const std::string TMS_ScintLayerName = "scinBoxlvRMMS_PV";
  // The top layer name
  const std::string TMS_TopLayerName = "volWorld";
}

#endif
