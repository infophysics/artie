#ifndef PrimaryGenerator_h
#define PrimaryGenerator_h 1

#include "globals.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4VUserActionInitialization.hh"
#include "G4ParticleGun.hh"
#include "G4GenericMessenger.hh"
#include "TFile.h"
#include "TH1D.h"
#include "TGraph.h"
#include "TMath.h"


class G4Event;

class PrimaryGenerator : public G4VUserPrimaryGeneratorAction
{
public:
  class Init : public G4VUserActionInitialization
  {
  public:
    Init();
    virtual ~Init();
    virtual void BuildForMaster() const;
    virtual void Build() const;
  };

  PrimaryGenerator();    
  ~PrimaryGenerator();
  void UseUserDefinedEnergy(TString inputfile, TString histname);
  void SetEnergyRange(G4double elow, G4double ehigh);
  virtual void GeneratePrimaries(G4Event*);
  const G4ParticleGun* GetParticleGun() const {return fParticleGun;};
  

private:
  G4ParticleGun*  fParticleGun;        //pointer a to G4 service class
  TFile* fEnergyFile = 0; 
  TH1D* fEnergy = 0;
  G4GenericMessenger* fMessenger;      //pinter to generic messanger class
  G4bool fUseUserDefinedEnergy = false;  // energy option
  G4double fEnergyLowCut = 0; // in keV
  G4double fEnergyHighCut = 0; // in keV
};

#endif


