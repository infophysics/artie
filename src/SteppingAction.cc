//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file SteppingAction.cc
/// \brief Implementation of the SteppingAction class
//
// $Id: SteppingAction.cc 71404 2013-06-14 16:56:38Z maire $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"
#include "Run.hh"
#include "TrackingAction.hh"
#include "HistoManager.hh"
#include "Randomize.hh"
#include "DetectorConstruction.hh"
#include "EventAction.hh"

#include "G4RunManager.hh"
                           
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(DetectorConstruction* det, TrackingAction* TrAct, EventAction* event)
: G4UserSteppingAction(), fDetector(det), fTrackingAction(TrAct), fEventAction(event)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  //track informations
  const G4StepPoint* prePoint = aStep->GetPreStepPoint();   
  const G4StepPoint* postPoint = aStep->GetPostStepPoint();
  const G4ParticleDefinition* particle = aStep->GetTrack()->GetDefinition(); 
  const G4VPhysicalVolume* preStepPhysical = prePoint->GetPhysicalVolume();
  const G4VPhysicalVolume* postStepPhysical = postPoint->GetPhysicalVolume();

  //RNG
  G4double mu = 0.*CLHEP::ns;
  G4double sigma = 120.0*CLHEP::ns;
  
  // The track does not exist
  if(preStepPhysical == 0 || postStepPhysical == 0) return;
  // Both steps are in the World
  if(preStepPhysical->GetCopyNo() == -1 && postStepPhysical->GetCopyNo() == -1) return;
  
  // count processes
  const G4VProcess* process   = postPoint->GetProcessDefinedStep();
  Run* run = static_cast<Run*>(
        G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  run->CountProcesses(process);
  
  // get step volumes
  G4LogicalVolume* preVolume = prePoint->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
  G4LogicalVolume* endVolume = postPoint->GetTouchableHandle()->GetVolume()->GetLogicalVolume();  

  // incident neutron
  //
  if (aStep->GetTrack()->GetTrackID() == 1) { 
  	G4double ekin  = postPoint->GetKineticEnergy();
    G4double trackl = aStep->GetTrack()->GetTrackLength();
    G4double time   = aStep->GetTrack()->GetLocalTime();           
    fTrackingAction->UpdateTrackInfo(ekin,trackl,time);
    G4AnalysisManager::Instance()->FillH1(7,ekin);	  
  }
/*
  if (aStep->GetTrack()->GetTrackID() == 1
  && aStep->GetTrack()->GetDefinition()->GetParticleName() == "neutron"
  && postPoint->GetStepStatus() == fGeomBoundary // step crossing bondary
  && preVolume == fDetector->GetLogicTarget() // step prevolume is LAr target
  && endVolume == fDetector->GetLogicLArcontainer()) { // step postvolume is LAr container
    G4double KinEOut = postPoint->GetKineticEnergy();
    G4AnalysisManager::Instance()->FillH1(10,KinEOut);  
  }

  if (aStep->GetTrack()->GetTrackID() == 1
  && aStep->GetTrack()->GetDefinition()->GetParticleName() == "neutron"
  && postPoint->GetStepStatus() == fGeomBoundary // step crossing bondary
  && preVolume == fDetector->GetLogicLArcontainer() // step prevolume is LAr container
  && endVolume == fDetector->GetLogicTarget()) { // step postvolume is LAr target 
    G4double KinEIn = postPoint->GetKineticEnergy();
    G4AnalysisManager::Instance()->FillH1(11,KinEIn);  
  }
*/ 

  if (aStep->GetTrack()->GetDefinition()->GetParticleName() == "neutron"){ 

    G4double trackid = aStep->GetTrack()->GetTrackID();
    G4AnalysisManager::Instance()->FillH1(13,trackid);

  }

  if (aStep->GetTrack()->GetTrackID() == 1
  && aStep->GetTrack()->GetDefinition()->GetParticleName() == "neutron"
  && postPoint->GetStepStatus() == fGeomBoundary // step crossing bondary
  && preVolume == fDetector->GetLogicTarget() // step prevolume is LAr target
  && endVolume == fDetector->GetLogicLArcontainer()) { // step postvolume is LAr container
    fEventAction->AddStrayNeutronCounter();  // Adding to the counter
  }
  
  
  // Neutron Detected
  if (aStep->GetTrack()->GetTrackID() == 1
  && aStep->GetTrack()->GetDefinition()->GetParticleName() == "neutron"
  && postPoint->GetStepStatus() == fGeomBoundary // step crossing bondary
  && preVolume == fDetector->GetLogicBeamLineV() // step prevolume is beam line volume
  //&& preVolume == fDetector->GetLogicWorld()
  && endVolume == fDetector->GetLogicDetector()
  ){ // step postvolume is neutron detector) 
    
  	   G4double RandGausNum = G4RandGauss::shoot(mu,sigma);
  	   G4double ekin  = postPoint->GetKineticEnergy();
       G4double trackl = aStep->GetTrack()->GetTrackLength();
       G4double time   = aStep->GetTrack()->GetLocalTime() + RandGausNum;           
       //fTrackingAction->MarkTrackInfo(ekin,trackl,time);
    if(fEventAction->GetStrayNeutronCounter() == 0) { //including only those neutrons which didn't scatter out of the target
       G4AnalysisManager::Instance()->FillH1(6,time);	  
       G4AnalysisManager::Instance()->FillH1(8,ekin);	 
       G4AnalysisManager::Instance()->FillH2(0,ekin,time);
    } else {
      G4AnalysisManager::Instance()->FillH1(12,ekin);
    }	 
  }

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


