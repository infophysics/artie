// To set this as default, you need a .rootrc file in your home directory,
// containing the following line:
// Rint.Logon: /full/path/to/rootlogon.C

#include "TColor.h"
#include "TH1.h"
#include "TLatex.h"
#include "TROOT.h"
#include "TStyle.h"  
#include "ArtieStyle.h" 

void SetArtieStyle ()
{
  static TStyle* artieStyle = 0;
  std::cout << "\nApplying ARTIE style settings...\n" << std::endl ;
  if ( artieStyle==0 ) artieStyle = ArtieStyle();
  gROOT->SetStyle("artieStyle");
  gROOT->ForceStyle();
}

TStyle* ArtieStyle()
{
  printf("Using ARTIE default plot style \n");

  // Defaults to classic style, but that's OK, we can fix it
  TStyle* artieStyle = new TStyle("artieStyle", "ARTIE Style");

  // Centre title
  artieStyle->SetTitleAlign(22);
  artieStyle->SetTitleX(.5);
  artieStyle->SetTitleY(.95);
  artieStyle->SetTitleBorderSize(0);

  // No info box
  artieStyle->SetOptStat(0);

  //set the background color to white
  artieStyle->SetFillColor(10);
  artieStyle->SetFrameFillColor(10);
  artieStyle->SetCanvasColor(10);
  artieStyle->SetPadColor(10);
  artieStyle->SetTitleFillColor(0);
  artieStyle->SetStatColor(10);

  // set canvas options 
  artieStyle->SetCanvasDefW(600);
  artieStyle->SetCanvasDefH(500); 
  artieStyle->SetCanvasColor(0); // canvas...
  artieStyle->SetCanvasBorderMode(0);
  artieStyle->SetCanvasBorderSize(0);   
  artieStyle->SetPadBorderMode(0); 
  artieStyle->SetPadBottomMargin(0.18); //margins...
  artieStyle->SetPadTopMargin(0.12);
  artieStyle->SetPadLeftMargin(0.18);
  artieStyle->SetPadRightMargin(0.09);
  artieStyle->SetPadGridX(0); // grids, tickmarks
  artieStyle->SetPadGridY(0);
  artieStyle->SetPadTickX(1);
  artieStyle->SetPadTickY(1);
  artieStyle->SetFrameBorderMode(0);
  artieStyle->SetPaperSize(20,24); // US letter size 

  // Set the default line color for a fit function to be red
  artieStyle->SetFuncColor(kRed);

  // Marker settings
  artieStyle->SetMarkerStyle(kFullCircle);

  // No border on legends
  artieStyle->SetLegendBorderSize(0);

  // Disabled for violating ARTIE style guidelines
  // Scientific notation on axes
  //  TGaxis::SetMaxDigits(3);

  // Axis titles
  artieStyle->SetTitleSize(.06, "xyz");
  artieStyle->SetTitleOffset(1.2, "xyz");
  // More space for y-axis to avoid clashing with big numbers
  artieStyle->SetTitleOffset(1.3, "y");
  // This applies the same settings to the overall plot title
  artieStyle->SetTitleSize(.06, "");
  artieStyle->SetTitleOffset(.9, "");
  // Axis labels (numbering)
  artieStyle->SetLabelSize(.06, "xyz");
  artieStyle->SetLabelOffset(.005, "xyz");

  // Prevent ROOT from occasionally automatically zero-suppressing
  artieStyle->SetHistMinimumZero();

  // Thicker lines
  artieStyle->SetHistLineWidth(2);
  artieStyle->SetFrameLineWidth(2);
  artieStyle->SetFuncWidth(2);

  // Set the number of tick marks to show
  artieStyle->SetNdivisions(506, "xyz");

  // Set the tick mark style
  artieStyle->SetPadTickX(1);
  artieStyle->SetPadTickY(1);

  // Fonts
  const int kArtieFont = 42;
  artieStyle->SetStatFont(kArtieFont);
  artieStyle->SetLabelFont(kArtieFont, "xyz");
  artieStyle->SetTitleFont(kArtieFont, "xyz");
  artieStyle->SetTitleFont(kArtieFont, ""); // Apply same setting to plot titles
  artieStyle->SetTextFont(kArtieFont);
  artieStyle->SetLegendFont(kArtieFont);

  // Get moodier colours for colz
  const Int_t NRGBs = 5;
  const Int_t NCont = 255;
  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  artieStyle->SetNumberContours(NCont);
  
  return artieStyle;
}

// Put a "ARTIE Preliminary" tag in the corner
void ArtiePreliminary(TString data="")
{
  TLatex* prelim = new TLatex(.9, .95, "ARTIE "+data+" Preliminary");
  prelim->SetTextColor(kBlue);
  prelim->SetNDC();
  prelim->SetTextSize(2/30.);
  prelim->SetTextAlign(32);
  prelim->Draw();
}

// Put a "ARTIE Preliminary" tag on the right
void ArtiePreliminarySide(TString data="")
{
  TLatex* prelim = new TLatex(.93, .9, "ARTIE "+data+" Preliminary");
  prelim->SetTextColor(kBlue);
  prelim->SetNDC();
  prelim->SetTextSize(2/30.);
  prelim->SetTextAngle(270);
  prelim->SetTextAlign(12);
  prelim->Draw();
}

// Put a "ARTIE Simulation" tag in the corner
void ArtieSimulation()
{
  TLatex* prelim = new TLatex(.9, .95, "ARTIE Simulation");
  prelim->SetTextColor(kGray+1);
  prelim->SetNDC();
  prelim->SetTextSize(2/30.);
  prelim->SetTextAlign(32);
  prelim->Draw();
}

// Put a "ARTIE Simulation" tag on the right
void ArtieSimulationSide()
{
  TLatex* prelim = new TLatex(.93, .9, "ARTIE Simulation");
  prelim->SetTextColor(kGray+1);
  prelim->SetNDC();
  prelim->SetTextSize(2/30.);
  prelim->SetTextAngle(270);
  prelim->SetTextAlign(12);
  prelim->Draw();
}

void ArtieCenterTitles(TH1* histo)
{
  histo->GetXaxis()->CenterTitle();
  histo->GetYaxis()->CenterTitle();
  histo->GetZaxis()->CenterTitle();
}
