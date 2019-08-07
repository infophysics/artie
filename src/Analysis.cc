#include <iostream>
#include <algorithm>

#include "TFile.h"
#include "TTree.h"

#include "G4UIdirectory.hh"
#include "G4Step.hh"
#include "G4Event.hh"
#include "G4Run.hh"
#include "G4UIdirectory.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"

#include "Analysis.hh"
#include "PrimaryGenerator.hh"
#include "Detector.hh"

using namespace std;

Analysis * Analysis::instance_ = NULL;

Analysis * Analysis::instance(){
  static Analysis a;
  return instance_;
}

Analysis::Analysis()
: 
  G4UImessenger(),
  analysisDir_(0),
  ntupleFilenameCmd_(0),
  dummyIntCmd_(0),
  dummyDoubleCmd_(0),
  file_(0),
  ntuple_(0)
{ 
  // set default values:  
  ntuple_filename = "analysis.root";
  dummy_int     = 0;
  dummy_double  = 0.0;

  analysisDir_ = new G4UIdirectory("/analysis/");
  analysisDir_->SetGuidance("Analysis commands");

  ntupleFilenameCmd_ = new G4UIcmdWithAString("/analysis/ntuple_filename", this);
  ntupleFilenameCmd_->SetParameterName ("ntuple_filename", true);
  ntupleFilenameCmd_->SetDefaultValue (ntuple_filename);  

  dummyIntCmd_ = new G4UIcmdWithAnInteger("/analysis/dummy_int", this);
  dummyIntCmd_->SetParameterName ("dummy_int", true);
  dummyIntCmd_->SetDefaultValue (dummy_int);  

  dummyDoubleCmd_ = new G4UIcmdWithADouble("/analysis/dummy_double", this);
  dummyDoubleCmd_->SetParameterName ("dummy_double", true);
  dummyDoubleCmd_->SetDefaultValue (dummy_double);  

  instance_ = this;
}

Analysis::~Analysis()
{
  instance_ = NULL;
  delete ntupleFilenameCmd_;
  delete dummyDoubleCmd_;
  delete dummyIntCmd_;
  delete analysisDir_;
  if (file_) delete file_;
}

void Analysis::SetNewValue(G4UIcommand * command,G4String arg){
  std::istringstream is((const char *) arg);

  if (command == ntupleFilenameCmd_){
    is >> ntuple_filename;
    cout << "INFO:  Set ntuple filename argument " << arg << " parsed as " << ntuple_filename << "\n";
  }

  if (command == dummyIntCmd_){
    is >> dummy_int;
    cout << "INFO:  Set dummy int argument " << arg << " parsed as " << dummy_int << "\n";
  }

  if (command == dummyDoubleCmd_){
    is >> dummy_double;
    cout << "INFO:  Set dummy int argument " << arg << " parsed as " << dummy_double << "\n";
  }
}

void Analysis::Book()
{ 
  // Creating a tree container to handle histograms and ntuples.
  // This tree is associated to an output file.

 
  file_ = new TFile(ntuple_filename,"RECREATE");
  if (! file_) {
    G4cout << " Analysis::Book :" 
           << " problem creating the ROOT TFile "
           << G4endl;
    return;
  }
  
  // create ntuple
  ntuple_ = new TTree("artie", "");
  ntuple_->Branch("gen_energy",    &gen_energy,     "gen_energy/D");
  ntuple_->Branch("arrival_time",  &arrival_time,   "arrival_time/D");
  ntuple_->Branch("arrival_e",     &arrival_e,      "arrival_e/D");
  ntuple_->Branch("num_elastic",   &num_elastic,    "num_elastic/I");
  ntuple_->Branch("num_inelastic", &num_inelastic,  "num_elastic/I");
  ntuple_->Branch("num_ncapture",  &num_ncapture,   "num_elastic/I");
  ntuple_->Branch("num_fission",   &num_fission,    "num_elastic/I");
  ntuple_->Branch("num_scatter",   &num_scatter,    "num_scatter/I");
  ntuple_->Branch("num_scatout",   &num_scatout,    "num_scatout/I");
  ntuple_->Branch("z_scatter",     &z_scatter,      "z_scatter/I");
  ntuple_->Branch("first_gas",     &first_gas,      "first_gas/I");
  ntuple_->Branch("max_dphi",      &max_dphi,       "max_dphi/D");
  ntuple_->Branch("max_dp",        &max_dp,         "max_dp/D");
  ntuple_->Branch("max_de",        &max_de,         "max_de/D");

  // example variable length ntuple variable:
  //// save reduced grid:
  //ntuple_->Branch("cell_max", &cell_max, "cell_max/I");
  //ntuple_->Branch("cell_edep",cell_edep,"cell_edep[cell_max]/F");

  G4cout << "\n----> Output file is open in " << ntuple_filename << G4endl;

  // Reset run parameters:
  num_events       = 0;
  events_detected  = 0;
  events_elastic   = 0; 
  events_inelastic = 0; 
  events_ncapture  = 0; 
  events_fission   = 0; 
  events_scatter   = 0;
  events_scatdet   = 0;
}

void Analysis::Save()
{ 
  if (! file_) return;
  
  file_->Write();
  file_->Close();
  
  G4cout << "\n----> Histograms and ntuples are saved\n" << G4endl;
}

void Analysis::FillNtuple(){
  num_events++;

  const PrimaryGenerator* generator
   = static_cast<const PrimaryGenerator*>
     (G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
  const G4ParticleGun* particleGun = generator->GetParticleGun();
  gen_energy = particleGun->GetParticleEnergy();


  if (!ntuple_) return;
  ntuple_->Fill();
}

void Analysis::Step(const G4Step* step)
{
  //obtain the detector (needed for volumes)
  const Detector* detector
   = static_cast<const Detector*>
   (G4RunManager::GetRunManager()->GetUserDetectorConstruction()); 

  //basic step information
  const G4StepPoint* pre = step->GetPreStepPoint();   
  const G4StepPoint* post = step->GetPostStepPoint();
  const G4VPhysicalVolume* preStepPhysical = pre->GetPhysicalVolume();
  const G4VPhysicalVolume* postStepPhysical = post->GetPhysicalVolume();
  // Sanity Checks: (crashes without these...)
  if(preStepPhysical == 0 || postStepPhysical == 0) return;
  if(preStepPhysical->GetCopyNo() == -1 && postStepPhysical->GetCopyNo() == -1) return;
  
  // Only consider the initial neutron
  //if (step->GetTrack()->GetTrackID() != 1) return;

  // Only consider neutrons:
  if (step->GetTrack()->GetDefinition()->GetParticleName() != "neutron") return;

  // remaining analysis is only for before striking target:
  if (arrival_time > 0){
    return;
  }
  
  // If we have just reached the detector, record the time and energy
  G4LogicalVolume* volume = post->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
  if (step->IsFirstStepInVolume()){    
    if (volume == detector->GetLogicDetector()){
      arrival_e = post->GetKineticEnergy();
      arrival_time   = step->GetTrack()->GetLocalTime();
    }  
  }

  // Keep track of how often each process occurred:
  const G4VProcess* process   = post->GetProcessDefinedStep();
  G4String procName = process->GetProcessName();
  if (procName == "hadElastic"){
    num_elastic++;
  }
  if (procName == "neutronInelastic"){
    num_inelastic++;
  }
  if (procName == "nCapture"){
    num_ncapture++;
  }
  if (procName == "nFission"){
    num_fission++;
  }
  // running tallies for all processes encountered:
  std::map<G4String,G4int>::iterator it = process_counts.find(procName);
  if ( it == process_counts.end()) {
    process_counts[procName] = 1;
  }
  else {
    process_counts[procName]++; 
  }

  // Quantify Scattering:

  // calculate delta phi, delta energy, and delta momentum.
  G4ThreeVector pa = pre->GetMomentumDirection();
  G4ThreeVector pb = post->GetMomentumDirection();
  double dphi = pa.angle(pb);
  double dp   = (pa - pb).mag();
  double de   = fabs(post->GetKineticEnergy() - pre->GetKineticEnergy());

  // The current position 
  G4ThreeVector pos = step->GetTrack()->GetPosition();

  // check for scattering (|dp| > 0):
  //  -> in target gas
  //  -> outside target gas before detection
  //  -> note if first scatter was in target gas
  //  -> note the z position of the first scatter
  if (dp > 0){
    if (volume == detector->GetLogicTarget()){
      num_scatter++;
      if (num_scatter == 1){
	if (num_scatout == 0){
	  z_scatter = pos.z();
	  first_gas = 1;
	}
      }
    } else {
      num_scatout++;    
      if ((num_scatout == 1) && (num_scatter == 0)){
	z_scatter = pos.z();
      }
    }    
  }

  // maximum angle/momentum/energy change prior to arrival in detector
  if (volume != detector->GetLogicDetector()){  
    if (dp > max_dp){
      max_dp = dp;
    }
    if (de > max_de){
      max_de = de;
    }
    if (dphi > max_dphi){
      max_dphi    = dphi;
    }    
  }
}

void Analysis::BeginEvent(const G4Event*)
{ 
  // reset event parameters:
  gen_energy    = 0;
  arrival_time  = 0;
  arrival_e     = 0;
  num_elastic   = 0;
  num_inelastic = 0;
  num_ncapture  = 0;
  num_fission   = 0;
  num_scatter   = 0;
  num_scatout   = 0;
  first_gas     = 0;
  z_scatter     = 0;
  max_dphi      = 0;
  max_dp        = 0;
  max_de        = 0;
}

void Analysis::EndEvent(const G4Event*)
{   
  if (num_elastic > 0)
    events_elastic++;
  if (num_inelastic > 0)
    events_inelastic++;
  if (num_ncapture > 0)
    events_ncapture++;
  if (num_fission > 0)
    events_fission++;
  if (first_gas>0){
    events_scatter++;
  } else if (num_scatout>0){
    events_scatout++;
  }
  if (arrival_time > 0)
    events_detected++;
  if (((num_scatter > 0) || (num_scatout > 0)) && (arrival_time > 0))
    events_scatdet++;

  FillNtuple();
}

void Analysis::BeginRun(const G4Run*)
{ 
  //inform the runManager whether to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);
  // Prepare the ntuple
  Book();
}

void Analysis::EndRun(const G4Run* run)
{
  G4int run_events = run->GetNumberOfEvent();
  if (run_events == 0) return;
  Save();       

  G4cout << "************************************************" << G4endl;
  G4cout << " Process Counts:  " << G4endl;
  std::map<G4String,G4int>::iterator it;    
  for (it = process_counts.begin(); it != process_counts.end(); it++) {
    G4String procName = it->first;
    G4int    count    = it->second;
    G4cout << "\t" << procName << "= " << count;
    G4cout << G4endl;
  }
  G4cout << G4endl;
  G4cout << "************************************************" << G4endl;
  G4cout << " Run Summary:  " << G4endl;
  G4cout << " - Events in Run:      " << run_events      << G4endl;
  G4cout << " - Events in Ntuple:   " << num_events      << G4endl;
  G4cout << " - Reaching Detector:  " << events_detected << G4endl;
  G4cout << " - Events w/ Process Elastic Scatter:       " << events_elastic  << G4endl;
  G4cout << " - Events w/ Process Inelastic Scatter:     " << events_inelastic << G4endl;
  G4cout << " - Events w/ Process Neutron Capture:       " << events_ncapture  << G4endl;
  G4cout << " - Events w/ Process Fission:               " << events_fission   << G4endl;
  G4cout << " - Events w/ First Scatter in Target Gas:   " << events_scatter   << G4endl;
  G4cout << " - Events w/ First Scatter Outside Target:  " << events_scatout   << G4endl;
  G4cout << " - Events Detected after Scatter:           " << events_scatdet   << G4endl;
  G4cout << "************************************************" << G4endl;
}


// G4 zealously deletes our classes for us, so need this bit of misdirection for a clean exits:

class AnalysisRunAction :  public G4UserRunAction {
public:
  virtual void   BeginOfRunAction(const G4Run*r) {Analysis::instance()->BeginRun(r); }
  virtual void   EndOfRunAction(const G4Run*r) {Analysis::instance()->EndRun(r); }  
};

class AnalysisStepAction :  public G4UserSteppingAction{
  virtual void   UserSteppingAction(const G4Step* step) {Analysis::instance()->Step(step); }
};
													   
class AnalysisEventAction :   public G4UserEventAction{
  virtual void   BeginOfEventAction (const G4Event* evt) {Analysis::instance()->BeginEvent(evt); }
  virtual void   EndOfEventAction   (const G4Event* evt) {Analysis::instance()->EndEvent(evt);   }
};

void Analysis::Init::BuildForMaster() const
{
  SetUserAction(new AnalysisRunAction());
}

void Analysis::Init::Build() const
{
  SetUserAction(new AnalysisRunAction);
  SetUserAction(new AnalysisEventAction);
  SetUserAction(new AnalysisStepAction);
}  

