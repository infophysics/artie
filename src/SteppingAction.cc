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

#include "G4RunManager.hh"
                           
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(DetectorConstruction* det, TrackingAction* TrAct)
: G4UserSteppingAction(), fDetector(det), fTrackingAction(TrAct)
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
  
  if (aStep->GetTrack()->GetTrackID() == 1
  && aStep->GetTrack()->GetDefinition()->GetParticleName() == "neutron"
  && postPoint->GetStepStatus() == fGeomBoundary // step crossing bondary
  && preVolume == fDetector->GetLogicWorld() // step prevolume is air
  && endVolume == fDetector->GetLogicDetector()) { // step postvolume is neutron detector) 
  	   G4double ekin  = postPoint->GetKineticEnergy();
       G4double trackl = aStep->GetTrack()->GetTrackLength();
       G4double time   = aStep->GetTrack()->GetLocalTime();           
       //fTrackingAction->MarkTrackInfo(ekin,trackl,time);
       G4AnalysisManager::Instance()->FillH1(6,time);	  
       G4AnalysisManager::Instance()->FillH1(8,ekin);	  
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


