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

#include "Analysis.hh"
#include "PrimaryGenerator.hh"

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
  dummyIntCmd_(0),
  dummyDoubleCmd_(0),
  file_(0),
  ntuple_(0)
{ 
  // set default values:  
  dummy_int     = 0;
  dummy_double  = 0.0;

  analysisDir_ = new G4UIdirectory("/analysis/");
  analysisDir_->SetGuidance("Analysis commands");

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
  delete dummyDoubleCmd_;
  delete dummyIntCmd_;
  delete analysisDir_;
  if (file_) delete file_;
}

void Analysis::SetNewValue(G4UIcommand * command,G4String arg){
  std::istringstream is((const char *) arg);

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

  G4String fileName = "pdout.root";
  file_ = new TFile(fileName,"RECREATE");
  if (! file_) {
    G4cout << " Analysis::Book :" 
           << " problem creating the ROOT TFile "
           << G4endl;
    return;
  }
  
  // create ntuple
  ntuple_ = new TTree("artie", "");
  ntuple_->Branch("gen_energy",   &gen_energy,   "gen_energy/D");
  ntuple_->Branch("num_scatter",  &num_scatter,  "num_scatter/I");
  ntuple_->Branch("num_target",   &num_target,   "num_target/I");
  ntuple_->Branch("arrival_time", &arrival_time, "arrival_time/D");
  ntuple_->Branch("scatter_z", &arrival_time, "arrival_time/D");

  // example variable length ntuple variable:
  //// save reduced grid:
  //ntuple_->Branch("cell_max", &cell_max, "cell_max/I");
  //ntuple_->Branch("cell_edep",cell_edep,"cell_edep[cell_max]/F");

  G4cout << "\n----> Output file is open in " << fileName << G4endl;

  // Reset run parameters:
  num_events      = 0;
  events_scatter  = 0;
  events_target   = 0;
  events_detector = 0;
  events_scatdet  = 0;
  events_dblscat  = 0;
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
  (void)step;
}

void Analysis::BeginEvent(const G4Event*)
{ 
  // reset event parameters:
  gen_energy  = 0;
  num_scatter = 0;
  num_target  = 0;
}

void Analysis::EndEvent(const G4Event*)
{   
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
  
  // Run conditions
  //  note: There is no primary generator action object for "master"
  //        run manager for multi-threaded mode.
  //const PrimaryGeneratorAction* generatorAction
  // = static_cast<const PrimaryGeneratorAction*>
  //   (G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
  //G4String runCondition;
  //if (generatorAction)
  // {

  //  runCondition += G4BestUnit(particleEnergy,"Energy");
  //}
         
  G4cout << "***************************************" << G4endl;
  G4cout << " Run Summary:  " << G4endl;
  G4cout << " - Events in Run:  " << run_events << G4endl;
  G4cout << " - Our count:      " << num_events << G4endl;
  G4cout << "***************************************" << G4endl;
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

