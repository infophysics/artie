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
: G4VUserPrimaryGeneratorAction(),fParticleGun(0), fUseUserDefinedEnergy(false), fEnergyLowCut(0), fEnergyHighCut(1000)
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
  
  // 
  fMessenger = new G4GenericMessenger(this,"/primaryGenerator/", "...doc...");
  fMessenger->DeclareMethod("UseUserDefinedEnergy", &PrimaryGenerator::UseUserDefinedEnergy, "...doc...");
  fMessenger->DeclareMethod("SetEnergyRange", &PrimaryGenerator::SetEnergyRange, "...doc...");
}

PrimaryGenerator::~PrimaryGenerator()
{
  delete fParticleGun;
  if(fUseUserDefinedEnergy) {
    fEnergyFile->Close();
    delete fEnergyFile;
//    if(fEnergy!=0) {
//    	delete fEnergy; 
//    	fEnergy = 0;
//    }
  } 
}

void PrimaryGenerator::UseUserDefinedEnergy(TString inputfile, TString objectname)
{ 
	std::cout<<"PrimaryGenerator::UseUserDefinedEnergy(): Load energy distribution from a ROOT file"<<std::endl;
	// Open the energy file
	fUseUserDefinedEnergy = true;
  fEnergyFile = new TFile(inputfile);
  //fEnergy = (TH1D*)fEnergyFile->Get(objectname);
  TGraph *gr = (TGraph*)fEnergyFile->Get(objectname);
  
  // Make variable-bin histogram for beam energy
  const int nlogbins=500;        
  double xmin = 1.e-3; //eV
  double xmax = 1.e7; //eV
  double *xbins    = new double[nlogbins+1];
  double xlogmin = TMath::Log10(xmin);
  double xlogmax = TMath::Log10(xmax);
  double dlogx   = (xlogmax-xlogmin)/((double)nlogbins);
  for (int i=0;i<=nlogbins;i++) {
     double xlog = xlogmin+ i*dlogx;
     xbins[i] = TMath::Exp( TMath::Log(10) * xlog ); 
  }
  fEnergy = new TH1D("hBeam_Energy", "", nlogbins, xbins);
  auto nPoints = gr->GetN(); // number of points 
  for(int i=0; i < nPoints; ++i) {
     double x,y;
     gr->GetPoint(i, x, y); //eV
     if(x/1000>fEnergyLowCut && x/1000<fEnergyHighCut) fEnergy->Fill(x,y);
  }
}

void PrimaryGenerator::SetEnergyRange(G4double elow, G4double ehigh)
{ 
	std::cout<<"PrimaryGenerator::SetEnergyRange(): Set energy range for incident particles"<<std::endl;
	fEnergyLowCut = elow;
	fEnergyHighCut = ehigh;
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
  G4double e = 0;
  
  if(fUseUserDefinedEnergy) {
  	e = fEnergy->GetRandom()/1000*keV;	
  }
  else e = (fEnergyLowCut+(fEnergyHighCut-fEnergyLowCut)*G4UniformRand())*keV;

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
