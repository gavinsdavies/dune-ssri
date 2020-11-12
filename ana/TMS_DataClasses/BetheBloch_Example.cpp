#include <iostream>

#include "TGraph.h"
#include "TCanvas.h"

#include "BetheBloch.h"

// Simple example showing Bethe-Bloch ionisation for iron and polystyrene
// Add your favourite element in the header under fMaterial
int main() {
  // Number of calculations we wish to make
  const int npoints = 10000;
  // Graphs that we draw to
  TGraph *graph = new TGraph(npoints);
  TGraph *graph2 = new TGraph(npoints);

  // The calculator for Iron
  BetheBloch_Calculator iron(Material::kIron);
  // The calculator for Polystyrene
  BetheBloch_Calculator poly(Material::kPolyStyrene);

  // What energy we want to start drawing at
  double minimum = BetheBloch_Utils::Mm*1.2;
  for (int i = 0; i < npoints; ++i) {
    double energy = minimum+(i*0.5);
    // Calculate the dEdx for polystyrene
    double dedx = poly.Calc_dEdx(energy);
    // Calculate the dEdx for iron
    double dedx2 = iron.Calc_dEdx(energy);
    graph->SetPoint(i, energy, dedx);
    graph2->SetPoint(i, energy, dedx2);
  }
  TCanvas *canv = new TCanvas("canv", "canv", 1024, 1024);
  canv->Print("dedx.pdf[");

  // Draw
  graph->Draw("AL");
  graph->SetLineColor(kRed);
  canv->Print("dedx.pdf");

  graph2->Draw("AL");
  graph2->SetLineColor(kGreen);
  canv->Print("dedx.pdf");

  graph->Draw();
  graph2->Draw("same");
  canv->Print("dedx.pdf");

  canv->Print("dedx.pdf]");
}
