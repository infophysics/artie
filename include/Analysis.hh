#ifndef Analysis_h
#define Analysis_h 1

#include <map>

#include "globals.hh"
#include "G4UImessenger.hh"
#include "G4UserRunAction.hh"
#include "G4UserEventAction.hh"
#include "G4UserSteppingAction.hh"
#include "G4VUserActionInitialization.hh"

class G4UIdirectory;
class G4UIcmdWithADouble;
class G4UIcmdWithAnInteger;
class G4UIcmdWithAString;
class G4Step;
class G4Event;
class G4Run;
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
    Init() : G4VUserActionInitialization() { Analysis::instance();}
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
  G4String ntuple_filename;
  int    dummy_int;
  double dummy_double;  

  // set config parameters:
  virtual void SetNewValue(G4UIcommand * command,G4String arg);

private:
  G4UIdirectory*        analysisDir_; 
  G4UIcmdWithAString*   ntupleFilenameCmd_; 
  G4UIcmdWithAnInteger* dummyIntCmd_; 
  G4UIcmdWithADouble*   dummyDoubleCmd_;

  // run monitoring variables:
  int num_events;        // total number of events
  int events_detected;   // events reaching detector
  int events_elastic;    // events with elastic scatter
  int events_inelastic;  // events with inelastic scatter
  int events_ncapture;   // events with neutron capture
  int events_fission;    // events with neutron capture
  int events_scatter;    // events which scattered first in gas
  int events_scatout;    // events which scattered first outside gas
  int events_scatdet;    // events with scatter that still reached detector 
  std::map<G4String,G4int> process_counts;        

  // event variables:
  double gen_energy;   // energy of generated neutron
  double arrival_time; // true arrival time at detector
  double arrival_e;    // true arrival energy at detector
  int num_elastic;     // elastic scatters
  int num_inelastic;   // inelastic scatters
  int num_ncapture;    // neutron capture
  int num_fission;     // fission
  int num_scatter;     // the total number of scatters in the target gas
  int num_scatout;     // the number of scatters outside the target gas
  int first_gas;       // did neutron scatter in the target gas first?
  int z_scatter;       // the z position of the first scatter
  double max_dphi;     // maximum direction change in any step
  double max_dp;       // maximum momentum change in any step
  double max_de;       // maximum energy change in any step

  static Analysis * instance_;
  TFile * file_;
  TTree * ntuple_;
};


#endif

    
