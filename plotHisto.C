{
  gROOT->Reset();
  
  // Draw histos filled by Geant4 simulation 
  //   
  TFile f("Artie.root");      
  TCanvas* c1 = new TCanvas("c1", "  ");
  //c1->SetLogy(1);
  //c1->SetLogx(1);
  c1->cd();
  c1->Update();
  
  //TH1D* hist1 = (TH1D*)f.Get("hTofDetected");
  //hist1->Draw("HIST");

  TH1D* hist2 = (TH1D*)f.Get("hPrimEnergyDetector");
  //hist2->Draw("HIST");

  //TH2D* hist2 = (TH2D*)f.Get("hTof_vs_energy");
  //hist2->Draw("COLZ");

  TH1* h1 = new TH1D("h1", "Cross Section vs Energy", 1000, 0.04, 0.09);
  h1->GetXaxis()->SetTitle("MeV");
  h1->GetYaxis()->SetTitle("#sigma(E) (in barns)");

  TH1* h2 = new TH1D("h1", "Trans Coeff vs Energy", 1000, 0.04, 0.09);
  h2->GetXaxis()->SetTitle("MeV");
  h2->GetYaxis()->SetTitle("Trans Coeff");

  Double_t CroSec = 0;
  Double_t TransCoeff = 0;
  Double_t BinEntry = 0;

  for(int i=1; i<=1000; i++){
    BinEntry = hist2->GetBinContent(i);
    //produced 10,000,000 (or 1,000,000) neutrons with uniform random distribution
    //There are 1000 bins
    //So there are 10000 (or 1000) neutrons in each energy bin
    TransCoeff = BinEntry/1000000; // BinEntry/10000;
    h2->SetBinContent(i,TransCoeff);
    if (TransCoeff == 0)
    {
      CroSec = 0; //this not right. what should this be?
    } else {
      //cout << BinEntry << endl;
      CroSec = (-1/4.2)*TMath::Log(TransCoeff);
    }
    h1->SetBinContent(i,CroSec);
  }

  //h1->Draw("HIST");
  h2->Draw("HIST");



  
//  TH1D* hist2 = (TH1D*)f.Get("2");
//  hist2->Draw("HIST");
//  
//  TH1D* hist3 = (TH1D*)f.Get("3");
//  hist3->Draw("HIST");
//  
//  TH1D* hist4 = (TH1D*)f.Get("4");
//  hist4->Draw("HIST");
//  
//  TH1D* hist5 = (TH1D*)f.Get("5");
//  hist5->Draw("HIST");
//  
//  TH1D* hist6 = (TH1D*)f.Get("6");
//  hist6->Draw("HIST");
//  
//  TH1D* hist7 = (TH1D*)f.Get("7");
//  hist7->Draw("HIST");      
}  
