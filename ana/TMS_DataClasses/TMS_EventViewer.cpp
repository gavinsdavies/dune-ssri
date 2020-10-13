#include "TMS_EventViewer.h"

TMS_EventViewer::TMS_EventViewer() {
  gStyle->SetOptStat(0);
  gStyle->SetNumberContours(255);

  const double zmin = 700;
  const double zmax = 1450;
  const double xmin = -400;
  const double xmax = 400;
  const double ymin = -250;
  const double ymax = 100;
  // Scint bars are 4 by 1 cm
  const int nbinsz = (zmax-zmin)/4;
  const int nbinsx = (xmax-xmin)/4;
  const int nbinsy = (ymax-ymin)/4;

  // The 2D views
  xz_view = new TH2D("TMS_Viewer_xz", "TMS xz;z (cm); x (cm)", nbinsz, zmin, zmax, nbinsx, xmin, xmax);
  yz_view = new TH2D("TMS_Viewer_yz", "TMS yz;z (cm); y (cm)", nbinsz, zmin, zmax, nbinsy, ymin, ymax);

  xz_view->SetMinimum(-0.01);
  yz_view->SetMinimum(-0.01);
  xz_view->SetMaximum(4);
  yz_view->SetMaximum(4);

  // The canvas
  Canvas = new TCanvas("TMS_EventViewer", "TMS_EventViewer", 1024, 1024);
  Canvas->SetLeftMargin(Canvas->GetLeftMargin()*0.8);

  // Full view from inspecting all hits
  xz_box_Full = new TBox(730, -348.5, 1415, 348.5);
  xz_box_Full->SetLineColor(kGreen);
  xz_box_Full->SetFillStyle(0);

  yz_box_Full = new TBox(730, -234, 1415, 87);
  yz_box_Full->SetLineColor(kGreen);
  yz_box_Full->SetFillStyle(0);

  // FV just taking 50 cm in from the full
  xz_box_FV = new TBox(730, -300, 1365, 300);
  xz_box_FV->SetLineColor(kRed);
  xz_box_FV->SetLineStyle(kDashed);
  xz_box_FV->SetFillStyle(0);

  yz_box_FV = new TBox(730, -184, 1365, 37);
  yz_box_FV->SetLineColor(kRed);
  yz_box_FV->SetLineStyle(kDashed);
  yz_box_FV->SetFillStyle(0);
}

// Draw the finished event
void TMS_EventViewer::Draw(TMS_Event &event) {

  xz_view->Reset();
  yz_view->Reset();

  std::vector<TMS_Hit> TMS_Hits = event.GetHits();

  // Check that there are hits
  if (TMS_Hits.size() < 50) {
    std::cerr << "Trying to draw an event that has no hits in the TMS, returning..." << std::endl;
    return;
  }


  // Loop over the hits and add them
  for (std::vector<TMS_Hit>::iterator it = TMS_Hits.begin(); it != TMS_Hits.end(); ++it) {
    TMS_Bar bar = (*it).GetBar();  
    // Transform into the new basis
    double x = bar.GetX()/10-TMS_Const::TMS_Det_Offset[0];
    double y = bar.GetY()/10-TMS_Const::TMS_Det_Offset[1];
    double z = bar.GetZ()/10-TMS_Const::TMS_Det_Offset[2];
    int BarType = bar.GetBarType();
    double e = (*it).GetE();
    // Bar along x
    if (BarType == TMS_Bar::BarType::kXBar) {
      xz_view->Fill(z, x, e);
    }
    // Bar along y
    else if (BarType == TMS_Bar::BarType::kYBar) {
      yz_view->Fill(z, y, e);
    }
  }

  Canvas->Clear();
  Canvas->Divide(2);
  int EventNumber = event.GetEventNumber();
  TString Canvasname = Form("test_Event_%i", EventNumber);

  Canvas->cd(1); 
  xz_view->Draw("colz");
  xz_box_FV->Draw("same");
  xz_box_Full->Draw("same");

  Canvas->cd(2); 
  yz_view->Draw("colz");
  yz_box_FV->Draw("same");
  yz_box_Full->Draw("same");

  Canvas->Print(Canvasname+".pdf");
}
