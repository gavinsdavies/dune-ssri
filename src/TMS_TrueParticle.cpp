#include "TMS_TrueParticle.h"

void TMS_TrueParticle::Print() {
  std::cout << "Printing TMS_TrueParticle class: " << std::endl;
  std::cout << "  Parent: " << Parent << std::endl;
  std::cout << "  TrackId: " << TrackId << std::endl;
  std::cout << "  PDG: " << PDG << std::endl;
  std::cout << "  TMS momentum: " << TMS_Momentum.Mag() << std::endl;
  TMS_Momentum.Print();
  std::cout << "  TMS entry point: " << std::endl;
  TMS_EntryPoint.Print();
  std::cout << "  Size of points vector: " << PositionPoints.size() << std::endl;
  std::cout << "  Position and momentum at points: " << std::endl;
  for (size_t i = 0; i < PositionPoints.size(); ++i) {
    std::cout << "  Point " << i << "/" << PositionPoints.size() << std::endl;
    PositionPoints[i].Print();
    MomentumPoints[i].Print();
  }
}
