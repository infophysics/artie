#include "ArtieStyle.C"
//#include "ArtieTree.h"
#include "ArtieTree.C"
#include <TH2.h>
#include <TH1.h>
#include <TCanvas.h>   

//////////////////////////////////// Parameters ////////////////////////////////////
bool UseTrueEnergy = false;
double T0Res = 125.0; //ns
double TDetRes = 1.0; //ns
TRandom3 randgen;

//////////////////////////////////// Histograms ////////////////////////////////////
// beam energy 
TH1D *hBeam_Energy = 0; 
// True energy from MC
TH1D *hEnergy_argon = 0;
TH1D *hEnergy_vacuum = 0; 
TH1D *hTransmission = 0;
TH1D *hCrossSectiion = 0;

TCanvas *c1 = 0, *c2 = 0, *c3 = 0, *c4 = 0;

//////////////////////////////////// Functiions ////////////////////////////////////
void Init();
void Plot();
double GetEnergyFromTOF(double tof);
double GetTOFFromEnergy(double e);    

////////////////////////////////////  Main function ////////////////////////////////////
void ArtieAnalysis(TString argonfile, TString vacuumfile){    
	
	// Set plot style
	SetArtieStyle();
	
	// Initialize Histograms
	Init();
	
  // Load Root files
  ArtieTree targ(argonfile);
  ArtieTree tvac(vacuumfile);
  if (targ.fChain == 0 || tvac.fChain == 0) return;
  Long64_t nentries_arg = targ.fChain->GetEntriesFast();
  Long64_t nentries_vac = tvac.fChain->GetEntriesFast();
  Long64_t nentries = TMath::Min(nentries_arg, nentries_vac);

  // Loop over entries and fill histograms
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
     targ.LoadTree(jentry);
     tvac.LoadTree(jentry);
     targ.fChain->GetEntry(jentry);
     tvac.fChain->GetEntry(jentry);
     hBeam_Energy->Fill(targ.gen_energy*1.e6); // convert to eV
     if(UseTrueEnergy) { 
       hEnergy_argon->Fill(targ.arrival_e*1000); // convert to keV
       hEnergy_vacuum->Fill(tvac.arrival_e*1000); // convert to keV  	
     }  
     else {
     	 double tof_arg = targ.arrival_time/1000 + randgen.Uniform(-0.5*T0Res/1000, 0.5*T0Res/1000) + randgen.Gaus(0, TDetRes/1000);
     	 double tof_vac = tvac.arrival_time/1000 + randgen.Uniform(-0.5*T0Res/1000, 0.5*T0Res/1000) + randgen.Gaus(0, TDetRes/1000);;
     	 double energy_arg = GetEnergyFromTOF(tof_arg); // ns -> keV
     	 double energy_vac = GetEnergyFromTOF(tof_vac); // ns -> keV
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
    // Set Bin Content
    hTransmission->SetBinContent(i+1, TransCoeff);
    hTransmission->SetBinError(i+1, TransCoeffError);
    hCrossSectiion->SetBinContent(i+1, Xsection);
    hCrossSectiion->SetBinError(i+1, XsectionError);  
  }
  // Get Cross-section
  
  
  // Make plots
  Plot();
}

//////////////////////////////////// Get energy from Time-of-Flight  ////////////////////////////////////
double GetEnergyFromTOF(double tof){
  double MN = 939.*1000.0;  // keV 
  double C  = 3.e2;         // m / mus
  double L = 71.03;        //  m
  double beta = (L / tof) / C;
  double gamma = sqrt(1.0 / (1 - beta*beta));
  return (gamma - 1.0)*MN;    
}

// Get Time-of-Flight from Energy
double GetTOFFromEnergy(double e){
  double MN = 939.*1000.0;  // keV 
  double C  = 3.e2;         // m / mus
  double L = 71.03;        //  m
  double gamma = sqrt(1-1/(1+e/MN)/(1+e/MN));
  return L/C/gamma;
}

//////////////////////////////////// Initialize  ////////////////////////////////////
void Init() {
	
// Make variable-bin histogram for beam energy
  const int nlogbins=30;        
  double xmin = 1.; //eV
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
	double emin  = 40.;
  double emax  = 70.;
  double nbins = 30;
  
  hEnergy_argon = new TH1D("hEnergy_argon", "Argon-fill: energy spectrum", nbins, emin, emax);
  hEnergy_vacuum = new TH1D("hEnergy_vacuum", "Vacuum-fill: energy spectrum", nbins, emin, emax);
  hTransmission = new TH1D("hTransmission", "Transmission Coefficient", nbins, emin, emax);
  hCrossSectiion = new TH1D("hCrossSectiion", "Total cross-section", nbins, emin, emax);
  
  c1 = new TCanvas("c1", "c1");
  c2 = new TCanvas("c2", "c2");
  c3 = new TCanvas("c3", "c3");
  c4 = new TCanvas("c4", "c4");
  ArtieSimulation();
//  
 	
}

//////////////////////////////////// Make Plots ////////////////////////////////////
void Plot() {
	
	c1->cd();
	c1->Update();
	c1->SetLogx();
	c1->SetLogy();
	hBeam_Energy->GetXaxis()->SetTitle("Energy [eV]");
  hBeam_Energy->GetYaxis()->SetTitle("Count [/bin]");
  hBeam_Energy->Draw("HIST");
	
  c2->cd();	          
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
  
  TLegend *legend2 = new TLegend(0.55,0.7,0.85,0.85);
  legend2->AddEntry(hEnergy_vacuum,"Vacuum-fill","l");
  legend2->AddEntry(hEnergy_argon,"Argon-fill","l");
  legend2->Draw("same");
  
  c3->cd();
  hTransmission->GetXaxis()->SetTitle("Kinetic Energy (keV)");
  hTransmission->GetYaxis()->SetTitle("Transmission coefficient (/keV)");
  hTransmission->SetLineColor(4);
  hTransmission->SetMarkerStyle(20);
  hTransmission->SetMarkerColor(4);
  hTransmission->Draw("PE");
  
  c4->cd();
  c4->SetLogy();   
  hCrossSectiion->GetXaxis()->SetTitle("Kinetic Energy (keV)");
  hCrossSectiion->GetYaxis()->SetTitle("#sigma_{tot} (b)");
  hCrossSectiion->SetLineColor(4);
  hCrossSectiion->SetMarkerStyle(20);
  hCrossSectiion->SetMarkerColor(4);
  hCrossSectiion->Draw("PE");
  
  
}


