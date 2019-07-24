{
	gROOT->Reset();

	// Draw histos filled by Geant4 simulation 
     
  TFile f1("Artie_argon_fill_Alcont.root");      
/*  TFile f2("Artie_argon_fill_AlWin.root"); 
  TFile f3("Artie_argon_fill_TeflonWin.root"); 
  TFile f4("Artie_argon_fill_MylarWin.root"); 
*/  TFile f5("Artie_vacuum_fill_Alcont.root");
/*  TFile f6("Artie_vacuum_fill_AlWin.root"); 
  TFile f7("Artie_vacuum_fill_TeflonWin.root"); 
  TFile f8("Artie_vacuum_fill_MylarWin.root");
*/
  // TFile f9("Artie_argon_fill_Alcont.root");
  // TFile f10("Artie_vacuum_fill_Alcont.root");
  
  TCanvas* c1 = new TCanvas("c1", "  ");
  TCanvas* c2 = new TCanvas("c2", "  ");
  TCanvas* c3 = new TCanvas("c3", "  ");
  TCanvas* c4 = new TCanvas("c4", "  ");

  TCanvas* c5 = new TCanvas("c5", "  ");
  TCanvas* c6 = new TCanvas("c6", "  ");
  c1->SetLogy(1);
  //c5->SetLogy(1);
  //c1->SetLogx(1);
  //c1->cd();
  //c1->Update();

  TH1D* hist3 = (TH1D*)f1.Get("hTofDetected");
  hist3->SetTitle("Time of flight of the neutrons entering the detector (Argon Fill)");
  c3->cd();
  hist3->Draw("HIST");

  TH1D* hist4 = (TH1D*)f5.Get("hTofDetected");
  hist4->SetTitle("Time of flight of the neutrons entering the detector (Vacuum Fill)");
  c4->cd();
  hist4->Draw("HIST");



  TH1* h1 = new TH1D("h1", "Cross Section vs Energy", 100, 0, 0.1);
  h1->SetTitle("Cross Section vs Energy (Standard Method)");
  h1->GetXaxis()->SetTitle("MeV");
  h1->GetYaxis()->SetTitle("#sigma(E) (in barns)");

  TH1* h2 = new TH1D("h1", "Trans Coeff vs Energy", 100, 0, 0.1);
  h2->SetTitle("Trans Coeff vs Energy (Standard Method)");
  h2->GetXaxis()->SetTitle("MeV");
  h2->GetYaxis()->SetTitle("Trans Coeff");

  TH1* h3 = new TH1D("h3", "Energy of the neutrons entering the detector (Argon Fill)", 100, 0, 0.1);
  h3->GetXaxis()->SetTitle("MeV");
  TH1* h4 = new TH1D("h4", "Energy of the neutrons entering the detector (Vacuum Fill)", 100, 0, 0.1);
  h4->GetXaxis()->SetTitle("MeV");

  Double_t CroSec[4] = {0,0,0,0};
  Double_t TransCoeff[4] = {0,0,0,0};
  Double_t BinEntry1[4] = {0,0,0,0};
  Double_t BinEntry2[4] = {0,0,0,0};
  Double_t BinEntry3[4] = {0,0,0,0};

  Double_t tofArgon = 0;
//  Double_t tofVacuum = 0;

  Double_t tofArgonCount = 0;
  Double_t tofVacuumCount = 0;

  Double_t denom = 0;
 // Double_t someC = 0;


  Double_t mN = 939.565; //Neutron Mass
  Double_t c = 299.792458; //Speed of Light
  Double_t L = 71.0; //Flight Path length (from neutron source to the detector)

  Double_t Eargon = 0;
//  Double_t Evacuum = 0;

  for (int i = 1; i <= 1000; ++i)
  {
    tofArgonCount = hist3->GetBinContent(i);
    tofVacuumCount = hist4->GetBinContent(i);

    tofArgon = 0.1*i;

    denom = TMath::Sqrt(1 - (L*L)/(tofArgon*tofArgon*c*c));
    //cout << denom << endl;

    Eargon = mN * (1/denom - 1);
//    Evacuum = mN * (1/TMath::Sqrt(1-(L*L/tofVacuum*tofVacuum*c*c)) - 1);

//    cout << Eargon << endl;
    //cout << Evacuum << endl;

    h3->Fill(Eargon, tofArgonCount);
    h4->Fill(Eargon, tofVacuumCount);
  }

  c5->cd();
  h3->Draw("HIST");
  c6->cd();
  h4->Draw("HIST");

  for (int i = 1; i <=100; ++i)
  {
    BinEntry1[0] = h3->GetBinContent(i); //argon fill
    BinEntry2[0] = h4->GetBinContent(i); //vacuum fill

    if (BinEntry2[0] == 0)
    {
      TransCoeff[0] = 0;
    } else {
      TransCoeff[0] = BinEntry1[0]/BinEntry2[0];
    }

    if (TransCoeff[0] == 0)
    {
      CroSec[0] = 0; //this is not right. what should this be?
    } else {
      CroSec[0] = (-1/4.2)*TMath::Log(TransCoeff[0]);
    } 

    h1->SetBinContent(i,CroSec[0]);
    h2->SetBinContent(i,TransCoeff[0]); 
  }

  c1->cd();
  h1->Draw("HIST");
  c2->cd();
  h2->Draw("HIST");

  
}