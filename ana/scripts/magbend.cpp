void magbend(std::string filename) {
  gStyle->SetPalette(55);

  TFile *file = new TFile(filename.c_str());
  TTree *tree = (TTree*)file->Get("tree");
  //TChain *tree = new TChain("tree");
  //TString fileset = filename.c_str();
  //fileset = fileset(0,fileset.Last('/')+1);

  // Get directory, add all root files
  //tree->Add(fileset+"*.root");

  int nEntries = tree->GetEntries();
  std::vector<float> *xptLep = NULL;
  std::vector<float> *yptLep = NULL;
  std::vector<float> *zptLep = NULL;
  float vtx[3];
  float muonDeath[3];
  float muonBirth[3];
  float muonExitPt[3];
  float muonExitKE;
  int lepPdg;

  //tree->Print();
  //return;

  tree->SetBranchStatus("*", false);

  tree->SetBranchStatus("xptLep", true);
  tree->SetBranchAddress("xptLep", &xptLep);
  tree->SetBranchStatus("yptLep", true);
  tree->SetBranchAddress("yptLep", &yptLep);
  tree->SetBranchStatus("zptLep", true);
  tree->SetBranchAddress("zptLep", &zptLep);

  tree->SetBranchStatus("vtx", true);
  tree->SetBranchAddress("vtx", vtx);

  tree->SetBranchStatus("muonDeath", true);
  tree->SetBranchAddress("muonDeath", muonDeath);
  tree->SetBranchStatus("muonBirth", true);
  tree->SetBranchAddress("muonBirth", muonBirth);
  tree->SetBranchStatus("muonExitPt", true);
  tree->SetBranchAddress("muonExitPt", muonExitPt);

  tree->SetBranchStatus("muonExitKE", true);
  tree->SetBranchAddress("muonExitKE", &muonExitKE);

  tree->SetBranchStatus("lepPdg", true);
  tree->SetBranchAddress("lepPdg", &lepPdg);

  TH1D *LepAll = new TH1D("Lep_All", "Lep_All", 100, -200, 200);
  TH1D *LepCorrect = new TH1D("Lep_Correct", "Lep_Correct",  100, -200, 200);
  LepCorrect->SetLineColor(kBlue);
  TH1D *LepWrong = new TH1D("Lep_Wrong", "Lep_Wrong",  100, -200, 200);
  LepWrong->SetLineColor(kRed);

  for (int i = 0; i < nEntries; ++i) {
    tree->GetEntry(i);
    if (abs(vtx[0]) > 300. || abs(vtx[1]) > 100. || vtx[2] < 50. || vtx[2] > 350.) continue;

    if (muonBirth[2] > 735. || muonDeath[2] > 1365. || 
        muonDeath[1] > 90-50. || muonDeath[1] < -235+50. ||
        abs(muonDeath[0]) > 300. ||
        muonExitKE <= 0 ||
        abs(muonDeath[0]) > 165. || abs(muonDeath[0]) < 10. ) continue;

    double x1 = muonExitPt[0];
    double z1 = muonExitPt[2];

    double x2 = muonBirth[0];
    double z2 = muonBirth[2];

    double x3 = muonDeath[0];
    double z3 = muonDeath[2];

    double signed_dist = (-1*(z2-z1)*x3 + (x2-x1)*z3 + x1*z2 - z1*x2)/sqrt((x2-x1)*(x2-x1)+(z2-z1)*(z2-z1));

    LepAll->Fill(signed_dist);
    if ((lepPdg > 0 && signed_dist > 0) || (lepPdg < 0 && signed_dist < 0)) LepCorrect->Fill(signed_dist);
    else if ((lepPdg < 0 && signed_dist > 0) || (lepPdg > 0 && signed_dist < 0)) LepWrong->Fill(signed_dist);

  }

  TCanvas *canv = new TCanvas("canv", "canv", 1024, 1024);
  TString canvname = filename.c_str();
  canvname.ReplaceAll(".root", "_signdist");
  std::cout << canvname << std::endl;
  canv->Print(canvname+".pdf[");

  LepAll->Draw();
  canv->Print(canvname+".pdf");

  LepCorrect->Draw();
  canv->Print(canvname+".pdf");

  LepWrong->Draw();
  canv->Print(canvname+".pdf");

  LepAll->Draw();
  LepCorrect->Draw("same");
  LepWrong->Draw("same");

  LepCorrect->SetTitle(Form("%2.2f/%2.2f=%2.2f", LepCorrect->Integral(), LepAll->Integral(), double(LepCorrect->Integral())/LepAll->Integral()));
  LepWrong->SetTitle(Form("%2.2f/%2.2f=%2.2f", LepWrong->Integral(), LepAll->Integral(), double(LepWrong->Integral())/LepAll->Integral()));
  canv->BuildLegend();

  canv->Print(canvname+".pdf");
  canv->Print(canvname+".pdf]");

  TString outputname = canvname+".root";
  TFile *output = new TFile(outputname, "recreate");
  LepAll->Write();
  LepCorrect->Write();
  LepWrong->Write();
  output->Close();
}
