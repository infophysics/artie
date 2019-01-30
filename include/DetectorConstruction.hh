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
/// \file DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class
//
// $Id: DetectorConstruction.hh 66586 2012-12-21 10:48:39Z ihrivnac $
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4LogicalVolume;
class G4Material;
class DetectorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
  
    DetectorConstruction();
   ~DetectorConstruction();

  public:
  
    virtual G4VPhysicalVolume* Construct();

    G4Material* 
    MaterialWithSingleIsotope(G4String, G4String, G4double, G4int, G4int);
         
    void SetSize     (G4double);              
    void SetMaterial (G4String);            

  public:
  
     const G4VPhysicalVolume* GetWorld()      {return fPWorld;};           
                    
     G4double           GetSize()       {return fWorldSizeX;}; //usless function (J.Wang)     
     G4Material*        GetMaterial()   {return fWorldMater;};
     
     const G4LogicalVolume* GetLogicWorld()           {return fLWorld;};   
     const G4LogicalVolume* GetLogicDetector()        {return fLogicDetector;};   
     
     void               PrintParameters();
                       
  private:
  
     // world 
     G4VPhysicalVolume* fPWorld;
     G4LogicalVolume*   fLWorld;
     G4double           fWorldSizeX;
     G4double           fWorldSizeY;
     G4double           fWorldSizeZ;
     G4Material*        fWorldMater;   
    
     
     // gas container
     G4double           fGascontainerLength; 
     G4double           fGascontainerRadius; 
     G4LogicalVolume*   fLogicGascontainer;
     G4VPhysicalVolume* fPhysiGascontainer;
     G4Material*        fGascontainerMater; 
     
     // gas insulator
     G4double           fGasinsulatorLength; 
     G4double           fGasinsulatorRadius; 
     G4LogicalVolume*   fLogicGasinsulator;
     G4VPhysicalVolume* fPhysiGasinsulator;
     G4Material*        fGasinsulatorMater; 
     
     // liquid argon container   
     G4double           fLArcontainerLength; 
     G4double           fLArcontainerRadius; 
     G4LogicalVolume*   fLogicLArcontainer;
     G4VPhysicalVolume* fPhysiLArcontainer;
     G4Material*        fLArcontainerMater; 
     
     // liquid argon target
     G4double           fTargetLength; 
     G4double           fTargetRadius; 
     G4LogicalVolume*   fLogicTarget;
     G4VPhysicalVolume* fPhysiTarget;
     G4Material*        fTargetMater; 
     
     // kapton window to be added
     
     // neutron collimator
     G4double           fCollimatorShieldThickness;
     G4double           fCollimatorSolidLength;
     G4double           fCollimatorSolidRadius;
     G4double           fCollimatorHollowLength;
     G4double           fCollimatorHollowRadius;
     G4LogicalVolume*   fLogicCollimator;
     G4VPhysicalVolume* fPhysiCollimator;
     G4Material*        fCollimatorMater;
     
     // toy neutron detector
     G4double						fDetectorLength;
     G4double           fDetectorRadius;
     G4LogicalVolume*   fLogicDetector;
     G4VPhysicalVolume* fPhysiDetector;
     G4Material*        fDetectorMater; 
     G4double           fDetectorPositionZ;
     
     
     DetectorMessenger* fDetectorMessenger;

  private:
    
     void               DefineMaterials();
     G4VPhysicalVolume* ConstructVolumes();    
      
     
     
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#endif

