#include "TMS_TrueHit.h"

TMS_TrueHit::TMS_TrueHit() :
  x(-999.99),
  y(-999.99),
  z(-999.99),
  t(-999.99),
  EnergyDeposit(-999.99)
{};

TMS_TrueHit::TMS_TrueHit(double x, double y, double z, double t, double E) {
  SetX(x);
  SetY(y);
  SetZ(z);
  SetT(t);
  SetE(E);
}

TMS_TrueHit::TMS_TrueHit(TG4HitSegment &edep_seg) {
  // Set the energy
  SetE(edep_seg.GetEnergyDeposit());
  // Set the primary contributor

  TLorentzVector start = edep_seg.GetStart();
  SetX(start.X());
  SetY(start.Y());
  SetZ(start.Z());
  SetT(start.T());
}
