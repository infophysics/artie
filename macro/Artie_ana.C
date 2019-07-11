{ //void Artie_ana() 
  gROOT->Reset();
  
  // Draw histos filled by Geant4 simulation 
  //   
  TFile f1("Artie_argon_fill.root");      
  TFile f2("Artie_argon_fill_AlWin.root"); 
  TFile f3("Artie_argon_fill_TeflonWin.root"); 
  TFile f4("Artie_argon_fill_MylarWin.root"); 
  TFile f5("Artie_vacuum_fill.root");
  TFile f6("Artie_vacuum_fill_AlWin.root"); 
  TFile f7("Artie_vacuum_fill_TeflonWin.root"); 
  TFile f8("Artie_vacuum_fill_MylarWin.root");
  TCanvas* c1 = new TCanvas("c1", "  ");
  TCanvas* c2 = new TCanvas("c2", "  ");
  c1->SetLogy(1);
  //c1->SetLogx(1);
  c1->cd();
  c1->Update();

  //TH1D* hist1 = (TH1D*)f1.Get("hPrimEnergyBeam");
  
  TH1D* hist1 = (TH1D*)f1.Get("hPrimEnergyDetector");  
  //hist1->Draw("HIST");

  TH1D* hist2 = (TH1D*)f5.Get("hPrimEnergyDetector");
  //hist2->SetTitle("Energy Spectrum of transmitted Neutrons getting into the detector");
  //hist2->Draw("HIST");

  TH1D* hist3 = (TH1D*)f2.Get("hPrimEnergyDetector");  
  //hist1->Draw("HIST");

  TH1D* hist4 = (TH1D*)f6.Get("hPrimEnergyDetector");

  TH1D* hist5 = (TH1D*)f3.Get("hPrimEnergyDetector");  
  //hist1->Draw("HIST");

  TH1D* hist6 = (TH1D*)f7.Get("hPrimEnergyDetector");

  TH1D* hist7 = (TH1D*)f4.Get("hPrimEnergyDetector");  
  //hist1->Draw("HIST");

  TH1D* hist8 = (TH1D*)f8.Get("hPrimEnergyDetector");

/*
  TH1D* hist3 = (TH1D*)f.Get("hTofDetected");
  //hist3->Draw("HIST");

  TH2D* hist4 = (TH2D*)f.Get("hTof_vs_energy");
  //hist4->Draw("COLZ");

  TH1D* hist5 = (TH1D*)f.Get("hEnergyScatteringOut");
  hist5->SetTitle("Energy of the Neutrons scattering out of the LAr Target in to the LAr container (for Stainless Steel)");
  //hist5->Draw("HIST");

  TH1D* hist6 = (TH1D*)f.Get("hEnergyScatteringIn");
  hist6->SetTitle("Energy of the Neutrons scattering in to the LAr Target from the LAr container (For Stainless Steel)");
  //hist6->Draw("HIST");
*/
  
  TH1* h1 = new TH1D("h1", "Cross Section vs Energy (for different windows)", 1000, 0.04, 0.09);
  h1->GetXaxis()->SetTitle("MeV");
  h1->GetYaxis()->SetTitle("#sigma(E) (in barns)");

  TH1* h2 = new TH1D("h1", "Trans Coeff vs Energy (for different windows)", 1000, 0.04, 0.09);
  h2->GetXaxis()->SetTitle("MeV");
  h2->GetYaxis()->SetTitle("Trans Coeff");

  Double_t CroSec[4] = {0,0,0,0};
  Double_t TransCoeff[4] = {0,0,0,0};
  Double_t BinEntry1[4] = {0,0,0,0};
  Double_t BinEntry2[4] = {0,0,0,0};

  for(int i=1; i<=1000; i++){
    BinEntry1[0] = hist1->GetBinContent(i);
    BinEntry2[0] = hist2->GetBinContent(i);
    //produced 1,000,000 neutrons with uniform random distribution
    //There are 1000 bins
    //So there are 1000 neutrons in each energy bin
    TransCoeff[0] = BinEntry1[0]/BinEntry2[0];
    if (TransCoeff[0] == 0)
    {
      CroSec[0] = 0; //this is not right. what should this be?
    } else {
      CroSec[0] = (-1/4.2)*TMath::Log(TransCoeff[0]);
    }
    h1->SetBinContent(i,CroSec[0]);
    h2->SetBinContent(i,TransCoeff[0]);
  }

  h1->Draw("HIST");
  c2->cd();
  h2->Draw("HIST");

  //////////////////////// Aluminium

  TH1* h3 = new TH1D("h3", "Cross Section vs Energy", 1000, 0.04, 0.09);
  //h3->GetXaxis()->SetTitle("MeV");
  //h3->GetYaxis()->SetTitle("#sigma(E) (in barns)");

  TH1* h4 = new TH1D("h4", "Trans Coeff vs Energy", 1000, 0.04, 0.09);
  //h2->GetXaxis()->SetTitle("MeV");
  //h2->GetYaxis()->SetTitle("Trans Coeff");

  for(int i=1; i<=1000; i++){
    BinEntry1[1] = hist3->GetBinContent(i);
    BinEntry2[1] = hist4->GetBinContent(i);
    //produced 1,000,000 neutrons with uniform random distribution
    //There are 1000 bins
    //So there are 1000 neutrons in each energy bin
    TransCoeff[1] = BinEntry1[1]/BinEntry2[1];
    if (TransCoeff[1] == 0)
    {
      CroSec[1] = 0; //this is not right. what should this be?
    } else {
      CroSec[1] = (-1/4.2)*TMath::Log(TransCoeff[1]);
    }
    h3->SetBinContent(i,CroSec[1]);
    h4->SetBinContent(i,TransCoeff[1]);
  }

  h3->SetLineColor(kRed);
  h4->SetLineColor(kRed);

  c1->cd();
  h3->Draw("SAME");
  c2->cd();
  h4->Draw("SAME");

  //////////////////////// Teflon

  TH1* h5 = new TH1D("h1", "Cross Section vs Energy", 1000, 0.04, 0.09);
  //h1->GetXaxis()->SetTitle("MeV");
  //h1->GetYaxis()->SetTitle("#sigma(E) (in barns)");

  TH1* h6 = new TH1D("h1", "Trans Coeff vs Energy", 1000, 0.04, 0.09);
  //h2->GetXaxis()->SetTitle("MeV");
  //h2->GetYaxis()->SetTitle("Trans Coeff");

  for(int i=1; i<=1000; i++){
    BinEntry1[2] = hist5->GetBinContent(i);
    BinEntry2[2] = hist6->GetBinContent(i);
    //produced 1,000,000 neutrons with uniform random distribution
    //There are 1000 bins
    //So there are 1000 neutrons in each energy bin
    TransCoeff[2] = BinEntry1[2]/BinEntry2[2];
    if (TransCoeff[2] == 0)
    {
      CroSec[2] = 0; //this is not right. what should this be?
    } else {
      CroSec[2] = (-1/4.2)*TMath::Log(TransCoeff[2]);
    }
    h5->SetBinContent(i,CroSec[2]);
    h6->SetBinContent(i,TransCoeff[2]);
  }

  h5->SetLineColor(kGreen);
  h6->SetLineColor(kGreen);

  c1->cd();
  h5->Draw("SAME");
  c2->cd();
  h6->Draw("SAME");

  /////////////////////// Mylar

  TH1* h7 = new TH1D("h1", "Cross Section vs Energy", 1000, 0.04, 0.09);
  //h1->GetXaxis()->SetTitle("MeV");
  //h1->GetYaxis()->SetTitle("#sigma(E) (in barns)");

  TH1* h8 = new TH1D("h1", "Trans Coeff vs Energy", 1000, 0.04, 0.09);
  //h2->GetXaxis()->SetTitle("MeV");
  //h2->GetYaxis()->SetTitle("Trans Coeff");

  for(int i=1; i<=1000; i++){
    BinEntry1[3] = hist7->GetBinContent(i);
    BinEntry2[3] = hist8->GetBinContent(i);
    //produced 1,000,000 neutrons with uniform random distribution
    //There are 1000 bins
    //So there are 1000 neutrons in each energy bin
    TransCoeff[3] = BinEntry1[3]/BinEntry2[3];
    if (TransCoeff[3] == 0)
    {
      CroSec[3] = 0; //this is not right. what should this be?
    } else {
      CroSec[3] = (-1/4.2)*TMath::Log(TransCoeff[3]);
    }
    h7->SetBinContent(i,CroSec[3]);
    h8->SetBinContent(i,TransCoeff[3]);
  }

  h7->SetLineColor(kOrange);
  h8->SetLineColor(kOrange);

  c1->cd();
  h7->Draw("SAME");
  c2->cd();
  h8->Draw("SAME");

  c1->cd();
  auto legend1 = new TLegend(0.7,0.7,0.9,0.9);
   legend1->SetHeader("Legend","C"); // option "C" allows to center the header
   legend1->AddEntry(h1,"Kapton (0.00762 cm)","l");
   legend1->AddEntry(h3,"Aluminium (1 mm)","l");
   legend1->AddEntry(h5,"Teflon (1 mm)","l");
   legend1->AddEntry(h7,"Mylar (1 mm)","l");   
   legend1->Draw();

  c2->cd();
  auto legend2 = new TLegend(0.7,0.7,0.9,0.9);
   legend2->SetHeader("Legend","C"); // option "C" allows to center the header
   legend2->AddEntry(h2,"Kapton (0.00762 cm)","l");
   legend2->AddEntry(h4,"Aluminium (1 mm)","l");
   legend2->AddEntry(h6,"Teflon (1 mm)","l");
   legend2->AddEntry(h8,"Mylar (1 mm)","l");   
   legend2->Draw();

}