{
	gROOT->Reset();

	// Draw histos filled by Geant4 simulation 
     
  TFile f1("Artie_argon_fill.root");      
/*  TFile f2("Artie_argon_fill_AlWin.root"); 
  TFile f3("Artie_argon_fill_TeflonWin.root"); 
  TFile f4("Artie_argon_fill_MylarWin.root"); 
*/  TFile f5("Artie_vacuum_fill.root");
/*  TFile f6("Artie_vacuum_fill_AlWin.root"); 
  TFile f7("Artie_vacuum_fill_TeflonWin.root"); 
  TFile f8("Artie_vacuum_fill_MylarWin.root");
*/
  TFile f9("Artie_argon_fill_Alcont.root");
  TFile f10("Artie_vacuum_fill_Alcont.root");
  
  TCanvas* c1 = new TCanvas("c1", "  ");
  TCanvas* c2 = new TCanvas("c2", "  ");
  //TCanvas* c3 = new TCanvas("c3", "  ");
  c1->SetLogy(1);
  //c1->SetLogx(1);
  //c1->cd();
  //c1->Update();

  //TH1D* hist1 = (TH1D*)f1.Get("hPrimEnergyBeam");
  
  TH1D* hist1 = (TH1D*)f1.Get("hPrimEnergyDetector");
  //c3->cd();
  //hist1->Draw("HIST");
  //hist1->Draw("HIST");

  // TH1D* hist11 = (TH1D*)f1.Get("hEnergyStrayNeutrons");
  // 

  TH1D* hist2 = (TH1D*)f5.Get("hPrimEnergyDetector");
  TH1D* hist3 = (TH1D*)f1.Get("hPrimEnergyBeam");



  //hist2->SetTitle("Energy Spectrum of transmitted Neutrons getting into the detector");
  //hist2->Draw("HIST");
/*
  TH1D* hist3 = (TH1D*)f2.Get("hPrimEnergyDetector");  
  //hist1->Draw("HIST");

  TH1D* hist4 = (TH1D*)f6.Get("hPrimEnergyDetector");

  TH1D* hist5 = (TH1D*)f3.Get("hPrimEnergyDetector");  
  //hist1->Draw("HIST");

  TH1D* hist6 = (TH1D*)f7.Get("hPrimEnergyDetector");

  TH1D* hist7 = (TH1D*)f4.Get("hPrimEnergyDetector");  
  //hist1->Draw("HIST");

  TH1D* hist8 = (TH1D*)f8.Get("hPrimEnergyDetector");
*/
  TH1D* hist9 = (TH1D*)f9.Get("hPrimEnergyDetector");
  TH1D* hist10 = (TH1D*)f10.Get("hPrimEnergyDetector");
  TH1D* hist11 = (TH1D*)f9.Get("hPrimEnergyBeam");



  TH1* h1 = new TH1D("h1", "Cross Section vs Energy", 100, 0.04, 0.09);
  h1->GetXaxis()->SetTitle("MeV");
  h1->GetYaxis()->SetTitle("#sigma(E) (in barns)");

  TH1* h2 = new TH1D("h1", "Trans Coeff vs Energy", 100, 0.04, 0.09);
  h2->GetXaxis()->SetTitle("MeV");
  h2->GetYaxis()->SetTitle("Trans Coeff");

  Double_t CroSec[4] = {0,0,0,0};
  Double_t TransCoeff[4] = {0,0,0,0};
  Double_t BinEntry1[4] = {0,0,0,0};
  Double_t BinEntry2[4] = {0,0,0,0};
  Double_t BinEntry3[4] = {0,0,0,0};

  for(int i=1; i<=100; i++){
    BinEntry1[0] = hist1->GetBinContent(i); //argon fill
    BinEntry2[0] = hist2->GetBinContent(i); //vacuum fill
    BinEntry3[0] = hist3->GetBinContent(i); //beam neutrons
    //produced 1,000,000 neutrons with uniform random distribution
    //There are 1000 bins
    //So there are 1000 neutrons in each energy bin
    TransCoeff[0] = 1 + (BinEntry1[0]-BinEntry2[0])/BinEntry3[0];
    Double_t TransError = 0;
    Double_t ArgonError = 0;
    Double_t VacuumError = 0;
    Double_t BeamError = 0;
    Double_t CrossError = 0;
    if (TransCoeff[0] == 0)
    {
      CroSec[0] = 0; //this is not right. what should this be?
    } else {
      CroSec[0] = (-1/4.2)*TMath::Log(TransCoeff[0]);

      ArgonError = 1/TMath::Sqrt(BinEntry1[0]);
      VacuumError = 1/TMath::Sqrt(BinEntry2[0]);
      BeamError = 1/TMath::Sqrt(BinEntry3[0]);

      TransError = TransCoeff[0] * ( (BeamError + ArgonError + VacuumError)/(BinEntry1[0] + BinEntry2[0] - BinEntry3[0]) + BeamError/BinEntry3[0] );
      CrossError = (1/4.2)*(TransError/TransCoeff[0]);
    }
    h1->SetBinContent(i,CroSec[0]);
    h2->SetBinContent(i,TransCoeff[0]);
    h2->SetBinError(i, TransError);
    h1->SetBinError(i,CrossError);
        }

  cout << h1->GetBinContent(34) << ", " << h1->GetBinError(34) << endl;

  c1->cd();
  h1->Draw("HIST");
  c2->cd();
  h2->Draw("HIST");

  ////////////////// Al

  TH1* h3 = new TH1D("h3", "Cross Section vs Energy", 100, 0.04, 0.09);
  //h3->GetXaxis()->SetTitle("MeV");
  //h3->GetYaxis()->SetTitle("#sigma(E) (in barns)");

  TH1* h4 = new TH1D("h4", "Trans Coeff vs Energy", 100, 0.04, 0.09);
  //h2->GetXaxis()->SetTitle("MeV");
  //h2->GetYaxis()->SetTitle("Trans Coeff");

  for(int i=1; i<=100; i++){
    BinEntry1[1] = hist9->GetBinContent(i);
    BinEntry2[1] = hist10->GetBinContent(i);
    BinEntry3[1] = hist11->GetBinContent(i);
    Double_t TransError = 0;
    Double_t ArgonError = 0;
    Double_t VacuumError = 0;
    Double_t BeamError = 0;
    Double_t CrossError = 0;
    //produced 1,000,000 neutrons with uniform random distribution
    //There are 1000 bins
    //So there are 1000 neutrons in each energy bin
    TransCoeff[1] = 1 + (BinEntry1[1]-BinEntry2[1])/BinEntry3[1];
    if (TransCoeff[1] == 0)
    {
      CroSec[1] = 0; //this is not right. what should this be?
    } else {
      CroSec[1] = (-1/4.2)*TMath::Log(TransCoeff[1]);

      ArgonError = 1/TMath::Sqrt(BinEntry1[1]);
      VacuumError = 1/TMath::Sqrt(BinEntry2[1]);
      BeamError = 1/TMath::Sqrt(BinEntry3[1]);

      TransError = TransCoeff[1] * ( (BeamError + ArgonError + VacuumError)/(BinEntry1[1] + BinEntry2[1] - BinEntry3[1]) + BeamError/BinEntry3[1] );
      CrossError = (1/4.2)*(TransError/TransCoeff[1]);
    }
    h3->SetBinContent(i,CroSec[1]);
    h4->SetBinContent(i,TransCoeff[1]);
    h4->SetBinError(i, TransError);
    h3->SetBinError(i,CrossError);
  }

  cout << h3->GetBinContent(34) << ", " << h3->GetBinError(34) << endl;

  h3->SetLineColor(kRed);
  h4->SetLineColor(kRed);

  c1->cd();
  h3->Draw("SAME");
  c2->cd();
  h4->Draw("SAME");

  c1->cd();
  auto legend1 = new TLegend(0.7,0.7,0.9,0.9);
   legend1->SetHeader("Legend","C"); // option "C" allows to center the header
   legend1->AddEntry(h3,"Aluminium","l");
   legend1->AddEntry(h1,"Stainless Steel","l");  
   legend1->Draw();

  c2->cd();
  auto legend2 = new TLegend(0.7,0.7,0.9,0.9);
   legend2->SetHeader("Legend","C"); // option "C" allows to center the header
   legend2->AddEntry(h4,"Aluminium","l");
   legend2->AddEntry(h2,"Stainless Steel","l");  
   legend2->Draw();
}