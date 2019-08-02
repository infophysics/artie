#include "Randomize.hh"
#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "PrimaryGenerator.hh"
#include "Detector.hh"
#include "G4RunManager.hh"

PrimaryGenerator::PrimaryGenerator()
: G4VUserPrimaryGeneratorAction(),fParticleGun(0)
{
  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);
  
  // default particle kinematic

  G4ParticleDefinition* particle
           = G4ParticleTable::GetParticleTable()->FindParticle("neutron");
  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));
  fParticleGun->SetParticleEnergy(57*keV);
  fParticleGun->SetParticlePosition(G4ThreeVector(0.*cm,0.*cm,0.*cm));

}

PrimaryGenerator::~PrimaryGenerator()
{
  delete fParticleGun;
}

void PrimaryGenerator::GeneratePrimaries(G4Event* anEvent)
{
  //obtain the detector (needed for tzero location)
  const Detector* detector
   = static_cast<const Detector*>
    (G4RunManager::GetRunManager()->GetUserDetectorConstruction()); 
  double tz = detector->TzeroLocation();

  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));

  // tzero position should be taken from geometry, hard-coded for now:
  fParticleGun->SetParticlePosition(G4ThreeVector(0.*cm,0.*cm,tz));
//  fParticleGun->SetParticleEnergy(57*keV);
  //G4double e = 57*keV;//100*G4UniformRand()*keV; //
  G4double e = (40+50*G4UniformRand())*keV; //
  fParticleGun->SetParticleEnergy(e);
  fParticleGun->GeneratePrimaryVertex(anEvent);  
}

PrimaryGenerator::Init::Init()
 : G4VUserActionInitialization()
{}

PrimaryGenerator::Init::~Init()
{}

void PrimaryGenerator::Init::BuildForMaster() const
{}

void PrimaryGenerator::Init::Build() const
{
  PrimaryGenerator* primary = new PrimaryGenerator();
  SetUserAction(primary);    
}  
