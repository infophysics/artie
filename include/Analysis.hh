#ifndef Analysis_h
#define Analysis_h 1

#include "globals.hh"
#include "G4UImessenger.hh"
#include "G4UserRunAction.hh"
#include "G4UserEventAction.hh"
#include "G4UserSteppingAction.hh"
#include "G4VUserActionInitialization.hh"

class G4UIcmdWithADouble;
class G4UIcmdWithAnInteger;
class G4Step;
class G4Event;
class G4Run;
class G4UIdirectory;
class TFile;
class TTree;
class TH1D;
class G4LogicalVolume;

class Analysis: 
  public G4UImessenger
{
  public:
  class Init :  public G4VUserActionInitialization {
  public:
    Init() : G4VUserActionInitialization() {}
    virtual ~Init() {}  
    virtual void   BuildForMaster() const;
    virtual void   Build() const;          
  };

  Analysis();
  virtual ~Analysis();  
  
  static Analysis * instance();

  void   Step(const G4Step*);
  void   BeginRun(const G4Run*);
  void   EndRun(const G4Run*);
  void   BeginEvent(const G4Event*);
  void   EndEvent(const G4Event*);
  void   Book();
  void   Save();
  void   FillNtuple();      

  // config parameters:
  int    dummy_int;
  double dummy_double;  

  // set config parameters:
  virtual void SetNewValue(G4UIcommand * command,G4String arg);

private:
  G4UIdirectory*        analysisDir_; 
  //G4UIcmdWithAString*   ntupleFileCmd_; 
  G4UIcmdWithAnInteger* dummyIntCmd_; 
  G4UIcmdWithADouble*   dummyDoubleCmd_;

  // run monitoring variables:
  int num_events;      // total number of events
  int events_scatter;  // events scatter somewhere before detector
  int events_target;   // events sattering in target
  int events_detector; // events hitting detector
  int events_scatdet;  // events that scattered then reached detector
  int events_dblscat;  // events that scattered two or more times

  // event variables:
  double gen_energy;   // energy of generated neutron
  int    num_scatter;  // number of scatters
  int    num_target;   // number of target scatters
  double arrival_time; // true arrival time at detector
  double scatter_z;    // z location of first scatter

  static Analysis * instance_;
  TFile * file_;
  TTree * ntuple_;
};


#endif

    
