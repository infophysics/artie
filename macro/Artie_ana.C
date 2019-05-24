void Artie_ana() {
  gROOT->Reset();
  
  // Draw histos filled by Geant4 simulation 
  //   
  TFile f("Artie.root");      
  TCanvas* c1 = new TCanvas("c1", "  ");
  //c1->SetLogy(0);
  //c1->SetLogx(1);
  c1->cd();
  c1->Update();
  
  //TH1D* hist1 = (TH1D*)f.Get("hTofDetected");
  //hist1->Draw("HIST");

  TH1D* hist2 = (TH1D*)f.Get("hPrimEnergyDetector");
  //hist2->Draw("HIST");

  //TH2D* hist2 = (TH2D*)f.Get("hTof_vs_energy");
  //hist2->Draw("COLZ");

  TH1* h1 = new TH1D("h1", "Cross Section", 1000, 0.0, 0.1);

  Double_t CroSec = 0;
  Double_t TransCoeff = 0;
  Double_t BinEntry = 0;

  for(int i=1; i<=1000; i++){
    BinEntry = hist2->GetBinContent(i);
    //produced 10,000,000 neutrons with uniform random distribution
    //There are 1000 bins
    //So there are 10000 neutrons in each energy bin
    TransCoeff = BinEntry/10000;
    if (TransCoeff == 0)
    {
      CroSec = 0; //this not right. what should this be?
    } else {
      CroSec = (-1/4.2)*TMath::Log(TransCoeff);
    }
    h1->SetBinContent(i,CroSec);
  }

  h1->Draw("HIST");
  
}