

#include "ArtieStyle.C"
//#include "ArtieTree.h"
#include "ArtieTree.C"
#include <TH2.h>
#include <TH1.h>
#include <TCanvas.h>   

//////////////////////////////////// Parameters ////////////////////////////////////
bool _UseTrueEnergy = false;
bool _UseLogScale = true;
double _T0Res = 125.0; //ns
double _TDetRes = 1.0; //ns
double _FlightPath = 71.9; //m
double _Energy_min = 40.0; // keV
double _Energy_max = 70.0; // keV
TString _TimeSmearOption = "MC"; // "MC" or "Gaus"
TRandom3 _randgen;
// will create a configuration file to store the parameters

//////////////////////////////////// Histograms ////////////////////////////////////
// beam energy 
TH1D *hBeam_Energy = 0; 
// time of flight
TH1D *hTof_argon = 0;
TH1D *hTof_vacuum = 0;
// neutron energy
TH1D *hEnergy_argon = 0;
TH1D *hEnergy_vacuum = 0; 
TH1D *hTransmission = 0;
TH1D *hCrossSection = 0;
// moderator response function
TFile *fModeratorFunction = 0;
TH2D *hModeratorFunction = 0;

TCanvas *c[20];

//////////////////////////////////// Functiion declaration ////////////////////////////////////
void Init();
void Plot();
double GetEnergyFromTOF(double tof);
double GetTofFromEnergy(double e);   
TH2D* GetModeratorFunction(TString rootfile, TString histname);
double GetTimeSmear(double e, TH2D* hist);

////////////////////////////////////  Main function ////////////////////////////////////
void ArtieAnalysis(TString argonfile, TString vacuumfile){    
	
	// Set plot style
	SetArtieStyle();
	
	// Initialize Histograms
	Init();
	
	// Get Moderator function
  TFile *fModeratorFunc = new TFile("resolution13a.root");
  hModeratorFunction = (TH2D*)fModeratorFunc->Get("tally15");
	
  // Load Root files
  ArtieTree targ(argonfile);
  ArtieTree tvac(vacuumfile);
  if (targ.fChain == 0 || tvac.fChain == 0) return;
  Long64_t nentries_arg = targ.fChain->GetEntriesFast();
  Long64_t nentries_vac = tvac.fChain->GetEntriesFast();
  Long64_t nentries = TMath::Min(nentries_arg, nentries_vac);
  //nentries = 1000;

  // Loop over entries and fill histograms
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
  	 if(jentry%10000==0) cout<<"Event "<<jentry<<" starts"<<endl;
     targ.LoadTree(jentry);
     tvac.LoadTree(jentry);
     targ.fChain->GetEntry(jentry);
     tvac.fChain->GetEntry(jentry);
     hBeam_Energy->Fill(targ.gen_energy*1.e6); // convert to eV
     if(_UseTrueEnergy) { 
       hEnergy_argon->Fill(targ.arrival_e*1000); // convert to keV
       hEnergy_vacuum->Fill(tvac.arrival_e*1000); // convert to keV  	
     }
     else {
     	 // get true time-of-flight
     	 double tof_arg = targ.arrival_time/1000; // ns -> us
     	 double tof_vac = tvac.arrival_time/1000; // ns -> us
     	 if(_TimeSmearOption == "MC") {
     	 	 // use true energy and true tof to determine the relative moderator smearing
     	   double tof_arg_smear = GetTimeSmear(targ.arrival_e, hModeratorFunction);
     	   double tof_vac_smear = GetTimeSmear(tvac.arrival_e, hModeratorFunction);
     	 	 tof_arg = tof_arg*(1+tof_arg_smear);
         tof_vac = tof_vac*(1+tof_vac_smear);
     	 }
       else if(_TimeSmearOption == "Gaus") {
       	 // use Gaussian moderator response 
       	 double tof_arg_smear = _randgen.Gaus(0, _T0Res/1000/2.35); // FWHM -> sigma
     	   double tof_vac_smear = _randgen.Gaus(0, _T0Res/1000/2.35); // FWHM -> sigma
     	   tof_arg += tof_arg_smear;
         tof_vac += tof_vac_smear;
       }
     	 
     	 // add detector time resolution
     	 tof_arg += _randgen.Gaus(0, _TDetRes/1000);
     	 tof_vac += _randgen.Gaus(0, _TDetRes/1000);
     	 
     	 hTof_argon->Fill(tof_arg*1000); // ns
     	 hTof_vacuum->Fill(tof_vac*1000); // ns
     	 
     	 // get energy from TOF
     	 double energy_arg = GetEnergyFromTOF(tof_arg); // us -> keV
     	 double energy_vac = GetEnergyFromTOF(tof_vac); // us -> keV
     	 
       hEnergy_argon->Fill(energy_arg); 
       hEnergy_vacuum->Fill(energy_vac);
     }
  }
  
  // Background subtraction (to be implemented)
  
  // Loop over energy bins
  for (int i=0; i<hEnergy_argon->GetNbinsX(); i++) {
  	int ArgCount = 0, VacCount = 0;
    double ArgError = 0., VacError = 0.;
    double TransCoeff = 0., TransCoeffError = 0.;
    double Xsection = 0., XsectionError = 0.;
    if(hEnergy_vacuum->GetBinContent(i+1)!=0) { 
    	// Entries in each energy bin
    	ArgCount = hEnergy_argon->GetBinContent(i+1);
    	VacCount = hEnergy_vacuum->GetBinContent(i+1);
    	if(VacCount!=0 && ArgCount!=0) {
    	  // Transmission coefficient from energy 
    	  TransCoeff = ArgCount*1.0/VacCount;
    	  // Cross section from energy
        Xsection = (-1./4.2)*TMath::Log(TransCoeff);
    	  // Uncertainty estimate
    	  ArgError = TMath::Sqrt(hEnergy_argon->GetBinContent(i+1))/ArgCount;
        VacError = TMath::Sqrt(hEnergy_vacuum->GetBinContent(i+1))/VacCount;
    	  TransCoeffError = sqrt(pow(ArgError, 2) + pow(VacError, 2)); 
        XsectionError = (1/4.2)*TransCoeffError;
      }
    }
    // Set Bin Content
    hTransmission->SetBinContent(i+1, TransCoeff);
    hTransmission->SetBinError(i+1, TransCoeffError);
    hCrossSection->SetBinContent(i+1, Xsection);
    hCrossSection->SetBinError(i+1, XsectionError);  
  } 
  
  // Make plots
  Plot();
}

//////////////////////////////////// Get energy from Time-of-Flight  ////////////////////////////////////
double GetEnergyFromTOF(double tof){
  double MN = 939.*1000.0;  // keV 
  double C  = 3.e2;         // m / us
  //double L = 71.03;        //  m
  double L = _FlightPath;
  double beta = (L / tof) / C;
  double gamma = sqrt(1.0 / (1 - beta*beta));
  return (gamma - 1.0)*MN;    
}

// Get Time-of-Flight from Energy
double GetTofFromEnergy(double e){
  double MN = 939.*1000.0;  // keV 
  double C  = 3.e2;         // m / us
  //double L = 71.03;        //  m
  double L = _FlightPath;
  double gamma = sqrt(1-1/(1+e/MN)/(1+e/MN));
  return L/C/gamma;
}

//////////////////////////////////// Initialize  ////////////////////////////////////
void Init() {
	
// Make variable-bin histogram for beam energy
  const int nlogbins=500;        
  double xmin = 1.e-3; //eV
  double xmax = 1.e6; //eV
  double *xbins    = new double[nlogbins+1];
  double xlogmin = TMath::Log10(xmin);
  double xlogmax = TMath::Log10(xmax);
  double dlogx   = (xlogmax-xlogmin)/((double)nlogbins);
  for (int i=0;i<=nlogbins;i++) { 
     double xlog = xlogmin+ i*dlogx;
     xbins[i] = TMath::Exp( TMath::Log(10) * xlog ); 
  }
  hBeam_Energy = new TH1D("hBeam_Energy", "", nlogbins, xbins);
	
// Make fix-bin histograms
	int emin  = _Energy_min;
  int emax  = _Energy_max;
  int nbins = emax-emin;
  
  int tofmin = TMath::Floor(GetTofFromEnergy(emax)*1000); // keV -> ns
  int tofmax = TMath::Ceil(GetTofFromEnergy(emin)*1000); // keV -> ns
  int ntofbins = (tofmax - tofmin)/100; // 10 ns per bin
  cout<<"tofmin, tofmax = "<<tofmin<<", "<<tofmax<<endl;
  
  hTof_argon = new TH1D("hTof_argon", "Argon-fill: Time-of-Flight", ntofbins, tofmin, tofmax);
  hTof_vacuum = new TH1D("hTof_vacuum", "Vacuum-fill: Time-of-Flight", ntofbins, tofmin, tofmax);
  hEnergy_argon = new TH1D("hEnergy_argon", "Argon-fill: energy spectrum", nbins, emin, emax);
  hEnergy_vacuum = new TH1D("hEnergy_vacuum", "Vacuum-fill: energy spectrum", nbins, emin, emax);
  hTransmission = new TH1D("hTransmission", "Transmission Coefficient", nbins, emin, emax);
  hCrossSection = new TH1D("hCrossSection", "Total cross-section", nbins, emin, emax);
  
  ArtieSimulation(); // need a fix!
}

//////////////////////////////////// Moderator resonse ////////////////////////////////////
TH2D* GetModeratorFunction(TString rootfile, TString histname){
  TFile *f = new TFile(rootfile);
  TH2D* hist = (TH2D*)f->Get(histname);
  return hist;
} 

//////////////////////////////////// Energy smearing ////////////////////////////////////
double GetTimeSmear(double e, TH2D *hist2d) {
  TAxis *xaxis = hist2d->GetXaxis();
  int binx = xaxis->FindBin(e);
  TH1D* hist_py = hist2d->ProjectionY("", binx, binx+1);
  return hist_py->GetRandom()/7; // scaled to 70 m
}

//////////////////////////////////// Make Plots ////////////////////////////////////
void Plot() {
	
	TFile *f = new TFile("Artie_histograms.root","recreate");
	
	int ir= 0;
	c[ir] = new TCanvas(Form("c%d", ir), Form("c%d", ir), 20*ir, 20*ir, 800, 600);
	c[ir]->cd();
	c[ir]->Update();
	c[ir]->SetLogx();
	c[ir]->SetLogy();
	hBeam_Energy->GetXaxis()->SetTitle("Energy [eV]");
  hBeam_Energy->GetYaxis()->SetTitle("Count [/bin]");
  hBeam_Energy->Draw("HIST");
  hBeam_Energy->Write();
  
  ir++;
  c[ir] = new TCanvas(Form("c%d", ir), Form("c%d", ir), 20*ir, 20*ir, 800, 600);
  c[ir]->cd();	          
  hTof_vacuum->GetXaxis()->SetTitle("Time-of-Flight [ns]");
  hTof_vacuum->GetYaxis()->SetTitle("Count [/10ns]");
  hTof_vacuum->GetYaxis()->SetRangeUser(0,hEnergy_vacuum->GetMaximum()*1.6);
  hTof_vacuum->SetLineColor(4);
  hTof_vacuum->SetMarkerStyle(20);
  hTof_vacuum->SetMarkerColor(4);    
  hTof_vacuum->SetMinimum(0);
  hTof_vacuum->Draw("PE");
  hTof_argon->GetXaxis()->SetTitle("Time-of-Flight [ns]");
  hTof_argon->GetYaxis()->SetTitle("Count [/10ns]");
  hTof_argon->SetLineColor(2);
  hTof_argon->SetMarkerStyle(20);
  hTof_argon->SetMarkerColor(2);
  hTof_argon->Draw("PEsame");
  hTof_vacuum->Write();
  hTof_argon->Write();
  
  TLegend *legend1 = new TLegend(0.55,0.7,0.85,0.85);
  legend1->AddEntry(hEnergy_vacuum,"Vacuum-fill","l");
  legend1->AddEntry(hEnergy_argon,"Argon-fill","l");
  legend1->Draw("same");
	
  ir++;
  c[ir] = new TCanvas(Form("c%d", ir), Form("c%d", ir), 20*ir, 20*ir, 800, 600);
  c[ir]->cd();	          
  hEnergy_vacuum->GetXaxis()->SetTitle("Energy [keV]");
  hEnergy_vacuum->GetYaxis()->SetTitle("Count [/keV]");
  hEnergy_vacuum->GetYaxis()->SetRangeUser(0,hEnergy_vacuum->GetMaximum()*1.6);
  hEnergy_vacuum->SetLineColor(4);
  hEnergy_vacuum->SetMarkerStyle(20);
  hEnergy_vacuum->SetMarkerColor(4);    
  hEnergy_vacuum->SetMinimum(0);
  hEnergy_vacuum->Draw("PE");
  hEnergy_argon->GetXaxis()->SetTitle("Energy [keV]");
  hEnergy_argon->GetYaxis()->SetTitle("Count [/keV]");
  hEnergy_argon->SetLineColor(2);
  hEnergy_argon->SetMarkerStyle(20);
  hEnergy_argon->SetMarkerColor(2);
  hEnergy_argon->Draw("PEsame");
  hEnergy_vacuum->Write();
  hEnergy_argon->Write();
  
  TLegend *legend2 = new TLegend(0.55,0.7,0.85,0.85);
  legend2->AddEntry(hEnergy_vacuum,"Vacuum-fill","l");
  legend2->AddEntry(hEnergy_argon,"Argon-fill","l");
  legend2->Draw("same");
  
  ir++;
  c[ir] = new TCanvas(Form("c%d", ir), Form("c%d", ir), 20*ir, 20*ir, 800, 600);
  c[ir]->cd();
  hTransmission->GetXaxis()->SetTitle("Kinetic Energy (keV)");
  hTransmission->GetYaxis()->SetTitle("Transmission coefficient (/keV)");
  hTransmission->SetLineColor(4);
  hTransmission->SetMarkerStyle(20);
  hTransmission->SetMarkerColor(4);
  hTransmission->Draw("PE");
  hTransmission->Write();
  
  ir++;
  c[ir] = new TCanvas(Form("c%d", ir), Form("c%d", ir), 20*ir, 20*ir, 800, 600);
  c[ir]->cd();
  c[ir]->SetLogy();   
  hCrossSection->GetXaxis()->SetTitle("Kinetic Energy (keV)");
  hCrossSection->GetYaxis()->SetTitle("#sigma_{tot} (b)");
  hCrossSection->SetLineColor(4);
  hCrossSection->SetMarkerStyle(20);
  hCrossSection->SetMarkerColor(4);
  hCrossSection->Draw("PE");
  hCrossSection->Write();
  
  ir++;
  c[ir] = new TCanvas(Form("c%d", ir), Form("c%d", ir), 20*ir, 20*ir, 800, 600);
  c[ir]->cd();
  c[ir]->SetLogx(); 
  c[ir]->SetLogy(); 
  hModeratorFunction->GetYaxis()->SetTitle("#DeltaT/T_{true}");
  hModeratorFunction->Draw("colz");
  
//  ir++;
//  c[ir] = new TCanvas(Form("c%d", ir), Form("c%d", ir), 20*ir, 20*ir, 800, 600);
//  c[ir]->cd();
//  TAxis *xaxis = hModeratorFunction->GetXaxis();
//  int binx = xaxis->FindBin(57./1000);
//  TH1D *hist_py = hModeratorFunction->ProjectionY("", binx, binx+1);
//  hist_py->SetName("hsmear");
//  hist_py->Draw("HIST");
//  hist_py->Write();
  
  f->Close();
}
