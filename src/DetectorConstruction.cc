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
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
//
// $Id: DetectorConstruction.cc 70755 2013-06-05 12:17:48Z ihrivnac $
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SubtractionSolid.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4RunManager.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"

#include "G4RotationMatrix.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
:G4VUserDetectorConstruction(),
 fPWorld(0), fLWorld(0), fWorldMater(0), fDetectorMessenger(0)
{
  fWorldSizeX = 4*m;
  fWorldSizeY = 4*m;
  fWorldSizeZ = 200*m;
  DefineMaterials();

  //Room
  fRoomSizeX = 2*m;
  fRoomSizeY = 2*m;
  fRoomSizeZ = 180*m;

  fRoomThickness = 1*m;
  
  // liquid argon target
  fTargetLength = 200.*cm; 
  fTargetRadius = 2.5*cm/2;  // DN25 (ID = 25 mm), Flange size 2-1/8" (OD)
  
  // liquid argon container
  fLArcontainerLength = 200.*cm; 
  fLArcontainerInnerRadius = fTargetRadius;
  fLArcontainerOuterRadius = (1.+3./8)*fInch2cm/2*cm; // OD = 1-3/8"
  
  // kapton window
  fKaptonThickness = 0.00762*cm; //0.1*cm; //
  
  // thermal insulator
  fInsulatorLength = fLArcontainerLength+2*fKaptonThickness;
  fInsulatorInnerRadius = fLArcontainerOuterRadius;
  fInsulatorOuterRadius = 1.5*fInch2cm/2*cm; //1.5" OD
  
  // Insulator container
  fInsulatorContainerLength = fLArcontainerLength  + 2*fKaptonThickness;
  fInsulatorContainerInnerRadius = fLArcontainerOuterRadius;
  fInsulatorContainerOuterRadius = 2.75*fInch2cm/2*cm; 
  
  
  //Fill and vent lines
  fPipeRadius = 2.*cm;
  fPipeLength = 10.*cm;
 
  // neutron collimator to be defined
  fCollimatorShieldThickness = 10.*cm;
  fCollimatorHollowLength = 90.*cm;   
  fCollimatorHollowRadius = fTargetRadius;
  fCollimatorSolidLength = fCollimatorHollowLength + fCollimatorShieldThickness;    
  fCollimatorSolidRadius = fCollimatorHollowRadius + fCollimatorShieldThickness;    

  
  // toy neutron detector
  fDetectorLength = 20.*cm;
  fDetectorRadius = 2.*cm; // not known
  fDetectorPositionZ = 70.*m;
  fDetectorMessenger = new DetectorMessenger(this);
  
  // Buffer volume
  fBufferLength = 5*cm;    
  fBufferInnerRadius = fLArcontainerInnerRadius;
  fBufferOuterRadius = fLArcontainerOuterRadius;

  //Beam Line
  fBeamLineLength = fDetectorPositionZ - 0.5*fDetectorLength - 0.5*fInsulatorContainerLength- fBufferLength;
  fBeamLineRadiusOUT = 20.*cm;
  fBeamLineRadiusIN = 18.*cm;
  
  


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ delete fDetectorMessenger;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  return ConstructVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
  // specific element name for thermal neutronHP
  // (see G4ParticleHPThermalScatteringNames.cc)

  G4int Z, A, a, ncomponents, natoms;
  G4double fractionmass, abundance;

  // pressurized water
  G4Element* H  = new G4Element("TS_H_of_Water" ,"H" , 1., 1.0079*g/mole);
  G4Element* O  = new G4Element("Oxygen"        ,"O" , 8., 16.00*g/mole);
  G4Material* H2O = 
  new G4Material("Water_ts", 1.000*g/cm3, ncomponents=2,
                         kStateLiquid, 593*kelvin, 150*bar);
  H2O->AddElement(H, natoms=2);
  H2O->AddElement(O, natoms=1);
  H2O->GetIonisation()->SetMeanExcitationEnergy(78.0*eV);
  
  // heavy water
  G4Isotope* H2 = new G4Isotope("H2",1,2);
  G4Element* D  = new G4Element("TS_D_of_Heavy_Water", "D", 1);
  D->AddIsotope(H2, 100*perCent);  
  G4Material* D2O = new G4Material("HeavyWater", 1.11*g/cm3, ncomponents=2,
                        kStateLiquid, 293.15*kelvin, 1*atmosphere);
  D2O->AddElement(D, natoms=2);
  D2O->AddElement(O, natoms=1);
  
  // world material
  G4Element* N  = new G4Element("Nitrogen", "N", 7, 14.01*g/mole);    
  G4Material* Air20 = new G4Material("Air", 1.205*mg/cm3, ncomponents=2, kStateGas, 293.*kelvin, 1.*atmosphere);
    Air20->AddElement(N, fractionmass=0.7);
    Air20->AddElement(O, fractionmass=0.3);

  // vacuum
  G4double atomicNumber = 1.;
  G4double massOfMole = 1.008*g/mole;
  G4double density = 1.e-25*g/cm3;
  G4double temperature = 2.73*kelvin;
  G4double pressure = 3.e-18*pascal;
  fVacuum = new G4Material("Vacuum", atomicNumber, massOfMole, density, kStateGas, temperature, pressure);
  
  // graphite
  G4Isotope* C12 = new G4Isotope("C12", 6, 12);  
  G4Element* C   = new G4Element("TS_C_of_Graphite","C", ncomponents=1);
  C->AddIsotope(C12, 100.*perCent);
  G4Material* graphite = new G4Material("graphite", 2.27*g/cm3, ncomponents=1,
                         kStateSolid, 293*kelvin, 1*atmosphere);
  graphite->AddElement(C, natoms=1);  
  
  // Define elements for all materials not found in the NIST database  
  G4NistManager* man = G4NistManager::Instance();
  G4Element* Li = man->FindOrBuildElement("Li");
  G4Element* B = man->FindOrBuildElement("B");
  G4Element* F = man->FindOrBuildElement("F");
  G4Element* Na = man->FindOrBuildElement("Na");
  G4Element* Mg = man->FindOrBuildElement("Mg");
  G4Element* Al = man->FindOrBuildElement("Al");
  G4Element* Si = man->FindOrBuildElement("Si");
  G4Element* K = man->FindOrBuildElement("K");
  G4Element* Ca = man->FindOrBuildElement("Ca");
  G4Element* Ti = man->FindOrBuildElement("Ti");
  G4Element* Cr = man->FindOrBuildElement("Cr");
  G4Element* Mn = man->FindOrBuildElement("Mn");
  G4Element* Fe = man->FindOrBuildElement("Fe");
  G4Element* Ni = man->FindOrBuildElement("Ni");
  G4Element* Sb = man->FindOrBuildElement("Sb");
  G4Element* Xe = man->FindOrBuildElement("Xe");
  G4Element* Cs = man->FindOrBuildElement("Cs");
  G4Element* Bi = man->FindOrBuildElement("Bi");
  

  // stainless steel
  G4Material* StainlessSteel = new G4Material("StainlessSteel", density= 8.06*g/cm3, ncomponents=6);
      StainlessSteel->AddElement(C, fractionmass=0.015); 
      StainlessSteel->AddElement(Si, fractionmass=0.008);
      StainlessSteel->AddElement(Cr, fractionmass=0.18);
      StainlessSteel->AddElement(Mn, fractionmass=0.01);
      StainlessSteel->AddElement(Fe, fractionmass=0.697);
      StainlessSteel->AddElement(Ni, fractionmass=0.09);

  //Aluminium
  G4Material* Aluminium = new G4Material("Aluminiun", density= 2.7*g/cm3, ncomponents=1);
      Aluminium->AddElement(Al, fractionmass=1);
	
  // Fe-56 isotope
  G4Isotope* iso_Fe = new G4Isotope("iso_Fe", Z=26, A=56, a=55.9349363*g/mole);
  G4Element* ele_Fe = new G4Element("ele_Fe", "Fe", ncomponents=1);
  ele_Fe->AddIsotope(iso_Fe,abundance=100.*perCent);
  G4Material* mat_Fe=new G4Material("mat_Fe",7.874*g/cm3, ncomponents = 1);
	mat_Fe->AddElement(ele_Fe, fractionmass = 1 );

	
  // Li-6 isotope
  G4Isotope* iso_Li = new G4Isotope("iso_Li", Z=3, A=6, a=6.015122795*g/mole);
  G4Element* ele_Li = new G4Element("ele_Li", "Li", ncomponents=1);
  ele_Li->AddIsotope(iso_Li,abundance=100.*perCent);
  G4Material* mat_Li=new G4Material("mat_Li",0.534*g/cm3, ncomponents = 1);
    mat_Li->AddElement(ele_Li, fractionmass = 1 );

  //Lithium Polyethylene
  G4Material* polyethylene = man->FindOrBuildMaterial("G4_POLYETHYLENE");
  G4Material* LiPoly = new G4Material("LiPoly", 1.06*g/cm3, ncomponents=2);
    LiPoly->AddElement (Li, 7.54*perCent);
    LiPoly->AddMaterial (polyethylene, 92.46*perCent);

  // Nitrogen, STP
  G4Material* Nitrogen = new G4Material("N2", 1.25053*mg/cm3, ncomponents=1);
  Nitrogen->AddElement(N, 2);

  // Polyurethane (foam insulator)
  G4Material* Polyurethane = new G4Material("Polyurethane", 1.005*g/cm3, ncomponents=4);
  Polyurethane->AddElement(C, 3);
  Polyurethane->AddElement(H, 8);
  Polyurethane->AddElement(N, 2);
  Polyurethane->AddElement(O, 1);

  	
  // world mater
  fWorldMater = Air20;//

  //Room Wall Mater
  fRoomMater = man->FindOrBuildMaterial("G4_CONCRETE");
  
  // insulator container
  fInsulatorContainerMater = Aluminium; //StainlessSteel; //
  
  // insulator
  fInsulatorMater = Polyurethane; //fVacuum;
  
  // LAr container
  fLArcontainerMater = Aluminium; //StainlessSteel; //
  	
  // liquid argon target
  fTargetMater = man->FindOrBuildMaterial("G4_lAr");
  //fTargetMater = fVacuum;
  
  // neutron collimator
  fCollimatorMater = LiPoly;
  
  // neutron detector
  fDetectorMater = H2O;

  // kapton
  fkapton  = man->FindOrBuildMaterial("G4_KAPTON"); //Aluminium; //

  //Beam Line
  fBeamLineMater = StainlessSteel;

  //Beam Line Volume
  fBeamLineVolumeMater = fVacuum;
  
  // Buffer volume
  fBufferMater = Nitrogen; //man->FindOrBuildMaterial("G4_N");

  
 ///G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* DetectorConstruction::MaterialWithSingleIsotope( G4String name,
                           G4String symbol, G4double density, G4int Z, G4int A)
{
 // define a material from an isotope
 //
 G4int ncomponents;
 G4double abundance, massfraction;

 G4Isotope* isotope = new G4Isotope(symbol, Z, A);
 
 G4Element* element  = new G4Element(name, symbol, ncomponents=1);
 element->AddIsotope(isotope, abundance= 100.*perCent);
 
 G4Material* material = new G4Material(name, density, ncomponents=1);
 material->AddElement(element, massfraction=100.*perCent);

 return material;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::ConstructVolumes()
{
  // Cleanup old geometry
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  // world
  G4Box*
  sWorld = new G4Box("Container", 
                   fWorldSizeX/2,fWorldSizeY/2,fWorldSizeZ/2); 

  fLWorld = new G4LogicalVolume(sWorld, 
                             fWorldMater,
                             fWorldMater->GetName()); 

  fPWorld = new G4PVPlacement(0, 
                            G4ThreeVector(),            
                            fLWorld,                    
                            fWorldMater->GetName(),     
                            0,                          
                            false,                      
                            0);    

  //Room
  G4Box*
  sRoom = new G4Box("Room_s", fRoomSizeX/2, fRoomSizeY/2, fRoomSizeZ/2);

  fLogicRoom = new G4LogicalVolume(sRoom, fRoomMater, "Room_l");

  fPhysiRoom = new G4PVPlacement(0, 
                                G4ThreeVector(),
                                fLogicRoom,
                                "Room_p",
                                fLWorld,
                                false,
                                0);

  //Room Volume
  G4Box*
  sRoomVolume = new G4Box("RoomVolume_s", (fRoomSizeX - fRoomThickness)/2, (fRoomSizeY - fRoomThickness)/2, (fRoomSizeZ - fRoomThickness)/2);

  fLogicRoomVolume = new G4LogicalVolume(sRoomVolume, fWorldMater, "RoomVolume_l");

  fPhysiRoomVolume = new G4PVPlacement(0,
                                G4ThreeVector(),
                                fLogicRoomVolume,
                                "RoomVolume_p",
                                fLogicRoom,
                                false,
                                0);

  
  // // Insulator container
  // G4Tubs* 
  // sInsulatorContainer = new G4Tubs("InsulatorContainer_s",                                               
  //                       fInsulatorContainerInnerRadius, fInsulatorContainerOuterRadius, 0.5*fInsulatorContainerLength, 0.,CLHEP::twopi);   

  // fLogicInsulatorContainer = new G4LogicalVolume(sInsulatorContainer,     
  //                                     fInsulatorContainerMater,     
  //                                     "InsulatorContainer_l");      

  // fPhysiInsulatorContainer = new G4PVPlacement(0,                   
  //                           				G4ThreeVector(),          
  //                                   fLogicInsulatorContainer ,      
  //                                   "InsulatorContainer_p",         
  //                                   fLWorld,                  
  //                                   false,                    
  //                                   0);                                                      
  // Insulator
  G4Tubs* 
  sInsulator = new G4Tubs("Insulator_s",                                               
                        fInsulatorInnerRadius, fInsulatorOuterRadius, 0.5*fInsulatorLength, 0.,CLHEP::twopi);   

  fLogicInsulator = new G4LogicalVolume(sInsulator,     
                                      fInsulatorMater,     
                                      "Insulator_l");      

  fPhysiInsulator = new G4PVPlacement(0,                   
                            				G4ThreeVector(),          
                                    fLogicInsulator ,      
                                    "Insulator_p",         
                                    fLogicRoomVolume,       
                                    false,                    
                                    0);                       
                                    
  // LAr container
  G4Tubs* 
  sLArcontainer = new G4Tubs("LArcontainer_s",                                               
                        fLArcontainerInnerRadius, fLArcontainerOuterRadius, 0.5*fLArcontainerLength, 0.,CLHEP::twopi);   

  fLogicLArcontainer = new G4LogicalVolume(sLArcontainer,   
                                      fLArcontainerMater,   
                                      "LArcontainer_l");    

  fPhysiLArcontainer = new G4PVPlacement(0,                 
                            				G4ThreeVector(),        
                                    fLogicLArcontainer ,    
                                    "LArcontainer_p",       
                                    fLogicRoomVolume,     
                                    false,                  
                                    0);                                                      
  
  // liquid argon target
  G4Tubs* 
  sTarget = new G4Tubs("Target_s",                                               
                        0, fTargetRadius, 0.5*fTargetLength, 0.,CLHEP::twopi);   

  fLogicTarget = new G4LogicalVolume(sTarget,                                    
                                      fTargetMater,                              
                                      "Target_l");           

  fPhysiTarget = new G4PVPlacement(0,                        
                            				G4ThreeVector(),         
                                    fLogicTarget ,           
                                    "Target_p",              
                                    fLogicRoomVolume,    //fLWorld, //  
                                    false,                   
                                    0);     
                                    
  // Kapton Window logic volume
  G4Tubs*
  sKapWin = new G4Tubs("KapWin_s", 0, fLArcontainerOuterRadius, 0.5*fKaptonThickness, 0., CLHEP::twopi);

  fLogicKapWin = new G4LogicalVolume(sKapWin, fkapton, "KapWin_l");
  
  
  //Kapton Window L
  fPhysiKapWinL1 = new G4PVPlacement(0,
  									G4ThreeVector(0,0, -fTargetLength*0.5-0.5*fKaptonThickness), 
  	                                fLogicKapWin,           
                                    "KapWinL1_p",
                                    fLogicRoomVolume,      
                                    false,                   
                                    0);     
  
  //Kapton Window R
  fPhysiKapWinR1 = new G4PVPlacement(0, 
  									G4ThreeVector(0,0, fTargetLength*0.5+0.5*fKaptonThickness), 
  	                                fLogicKapWin ,          
                                    "KapWinR1_p",              
                                    fLogicRoomVolume,  
                                    false, 
                                    1);

  // G4RotationMatrix* PipeRot = new G4RotationMatrix;
  // PipeRot->rotateY(90.0*deg);
/*
  //lAr fill line In
  G4Tubs*
  slArFillLineIN = new G4Tubs("lArFillLineIN_s", 0, fPipeRadius, 0.25*fPipeLength, 0., CLHEP::twopi);

  fLogiclArFillLineIN = new G4LogicalVolume(slArFillLineIN, fInsulatorContainerMater, "lArFillLineIN_l");

  fPhysilArFillLineIN = new G4PVPlacement(PipeRot, 
  									G4ThreeVector(fLArcontainerRadius+0.25*fPipeLength, 0, 0.25*fTargetLength), 
  	                                fLogiclArFillLineIN,           
                                    "lArFillLineIN_p",              
                                    fLogicInsulator,      
                                    false,                   
                                    0);

  //lAr fill line Out
  G4Tubs*
  slArFillLineOUT = new G4Tubs("lArFillLineOUT_s", 0, fPipeRadius, 0.25*fPipeLength, 0., CLHEP::twopi);

  fLogiclArFillLineOUT = new G4LogicalVolume(slArFillLineOUT, fInsulatorContainerMater, "lArFillLineOUT_l");

  fPhysilArFillLineOUT = new G4PVPlacement(PipeRot, 
  									G4ThreeVector(fInsulatorContainerRadius + 0.25*fPipeLength, 0, 0.25*fTargetLength), 
  	                                fLogiclArFillLineOUT ,           
                                    "lArFillLineOUT_p",              
                                    fLWorld,      
                                    false,                   
                                    0);


  //lAr Vent Line IN
  G4Tubs*
  slArVentLineIN = new G4Tubs("lArVentLineIN_s", 0, fPipeRadius, 0.25*fPipeLength, 0., CLHEP::twopi);

  fLogiclArVentLineIN = new G4LogicalVolume(slArVentLineIN, fInsulatorContainerMater, "lArVentLineIN_l");

  fPhysilArVentLineIN = new G4PVPlacement(PipeRot, 
  									G4ThreeVector(fLArcontainerRadius + 0.25*fPipeLength, 0, -0.25*fTargetLength), 
  	                                fLogiclArVentLineIN ,           
                                    "lArVentLineIN_p",              
                                    fLogicInsulator,      
                                    false,                   
                                    0);

  //lAr Vent Line Out
  G4Tubs*
  slArVentLineOUT = new G4Tubs("lArVentLineOUT_s", 0, fPipeRadius, 0.25*fPipeLength, 0., CLHEP::twopi);

  fLogiclArVentLineOUT = new G4LogicalVolume(slArVentLineOUT, fInsulatorContainerMater, "lArVentLineOUT_l");

  fPhysilArVentLineOUT = new G4PVPlacement(PipeRot, 
  									G4ThreeVector(fInsulatorContainerRadius + 0.25*fPipeLength, 0, -0.25*fTargetLength), 
  	                                fLogiclArVentLineOUT ,           
                                    "lArVentLineOUT_p",              
                                    fLWorld,      
                                    false,                   
                                    0);

  //lAr Vacuum Line
  G4Tubs*
  sVacuumLine = new G4Tubs("VacuumLine_s", 0, fPipeRadius, 0.25*fPipeLength, 0., CLHEP::twopi);

  fLogicVacuumLine = new G4LogicalVolume(sVacuumLine, fInsulatorContainerMater, "VacuumLine_l");

  fPhysiVacuumLine = new G4PVPlacement(PipeRot, 
  									G4ThreeVector(-fInsulatorContainerRadius - 0.25*fPipeLength, 0, 0), 
  	                                fLogicVacuumLine ,           
                                    "VacuumLine_p",              
                                    fLWorld,      
                                    false,                   
                                    0);
*/
  //Beam line
  G4Tubs*
  sBeamLine = new G4Tubs("BeamLine_s", 0, fBeamLineRadiusOUT, 0.5*fBeamLineLength, 0., CLHEP::twopi);

  fLogicBeamLine = new G4LogicalVolume(sBeamLine, fBeamLineMater, "BeamLine_l");

  fPhysiBeamLine = new G4PVPlacement(0,
                    G4ThreeVector(0, 0, 0.5*fInsulatorContainerLength + fBufferLength + 0.5*fBeamLineLength),
                                    fLogicBeamLine ,           
                                    "BeamLine_p",              
                                    fLogicRoomVolume,      
                                    false,                   
                                    0);

  //Beam line Volume
  G4Tubs*
  sBeamLineV = new G4Tubs("BeamLineV_s", 0, fBeamLineRadiusIN, 0.5*fBeamLineLength, 0., CLHEP::twopi);

  fLogicBeamLineV = new G4LogicalVolume(sBeamLineV, fBeamLineVolumeMater, "BeamLineV_l");

  fPhysiBeamLineV = new G4PVPlacement(0, 
                    G4ThreeVector(0, 0, 0), 
                                    fLogicBeamLineV ,           
                                    "BeamLineV_p",              
                                    fLogicBeamLine,      
                                    false,                   
                                    0);
//
//
//  // neutron collimator
//  G4Tubs*
//  sCollimatorSolid = new G4Tubs("CollimatorSolid",                        
//                   0, fCollimatorSolidRadius, fCollimatorSolidLength/2, 0.,CLHEP::twopi ); 
//  G4Tubs*
//  sCollimatorHollow = new G4Tubs("CollimatorHollow",                        
//                   0, fCollimatorHollowRadius, fCollimatorHollowLength/2, 0.,CLHEP::twopi );  
//  G4ThreeVector zTransCollimator(0, 0, -fCollimatorShieldThickness/2); 
//  G4SubtractionSolid* sCollimator =
//  new G4SubtractionSolid("CollimatorSolid-CollimatorHollow", sCollimatorSolid, sCollimatorHollow, 0, zTransCollimator); 
//  
//  fLogicCollimator = new G4LogicalVolume(sCollimator,                                    
//                                      fCollimatorMater,                              
//                                      "sCollimator_l");           
//
//  fPhysiCollimator = new G4PVPlacement(0,                        
//                            				G4ThreeVector(0., 0., fDetectorPositionZ + fDetectorLength/2 + fCollimatorShieldThickness - fCollimatorSolidLength/2),         
//                                    fLogicCollimator ,           
//                                    "sCollimator_p",              
//                                    fLWorld,      
//                                    false,                   
//                                    0);                      
// 

  G4Tubs* 
  sDetector = new G4Tubs("Detector_s",                                               
                        0, fDetectorRadius, 0.5*fDetectorLength, 0.,CLHEP::twopi);   

  fLogicDetector = new G4LogicalVolume(sDetector,                                    
                                      fDetectorMater,                              
                                      "Detector_l");           

  fPhysiDetector = new G4PVPlacement(0,                        
                            				G4ThreeVector(0., 0., fDetectorPositionZ),         
                                    fLogicDetector ,           
                                    "Detector_p",              
                                    fLogicRoomVolume,      
                                    false,                   
                                    0);  
 
 // buffer Container 
  G4Tubs* 
  sBuffer = new G4Tubs("Buffer_s",                                               
                        fBufferInnerRadius, fBufferOuterRadius, 0.5*fBufferLength, 0.,CLHEP::twopi);   

  fLogicBuffer = new G4LogicalVolume(sBuffer,                                    
                                      fLArcontainerMater,                              
                                      "Buffer_l");           

  fPhysiBufferL = new G4PVPlacement(0,                        
                            				G4ThreeVector(0., 0., -fInsulatorContainerLength/2-fBufferLength/2),         
                                    fLogicBuffer ,           
                                    "Buffer_p",              
                                    fLogicRoomVolume,      
                                    false,                   
                                    0);  
  fPhysiBufferR = new G4PVPlacement(0,                        
                            				G4ThreeVector(0., 0., fInsulatorContainerLength/2+fBufferLength/2),         
                                    fLogicBuffer ,           
                                    "Buffer_p",              
                                    fLogicRoomVolume,      
                                    false,                   
                                    1); 

  //buffer Volume

  G4Tubs* 
  sBufferVol = new G4Tubs("BufferVol_s",                                               
                        0, fBufferInnerRadius, 0.5*fBufferLength, 0.,CLHEP::twopi);   

  fLogicBufferVol = new G4LogicalVolume(sBufferVol,                                    
                                      fBufferMater,                              
                                      "BufferVol_l");           

  fPhysiBufferVolL = new G4PVPlacement(0,                        
                                    G4ThreeVector(0., 0., 0),         
                                    fLogicBufferVol ,           
                                    "BufferVol_p",              
                                    fLogicRoomVolume,      
                                    false,                   
                                    0);  
  fPhysiBufferVolR = new G4PVPlacement(0,                        
                                    G4ThreeVector(0., 0., 0),         
                                    fLogicBufferVol ,           
                                    "BufferVol_p",              
                                    fLogicRoomVolume,      
                                    false,                   
                                    1);


                                               
  
  // // set VisAttributes
  // fLWorld->SetVisAttributes(G4VisAttributes::GetInvisible()); 
  // G4VisAttributes* VisAttInsulatorContainer= new G4VisAttributes(G4Colour(0.0, 1.0, 1.0)); // cyan
  // fLogicInsulatorContainer->SetVisAttributes(VisAttInsulatorContainer);  
  // G4VisAttributes* VisAttInsulator= new G4VisAttributes(G4Colour(1.0, 0.0, 1.0)); // magenta 
  // fLogicInsulator->SetVisAttributes(VisAttInsulator);  
  // G4VisAttributes* VisAttLArcontainer= new G4VisAttributes(G4Colour(0.0, 1.0, 0.0)); // green
  // fLogicLArcontainer->SetVisAttributes(VisAttLArcontainer);
  // G4VisAttributes* VisAttTarget= new G4VisAttributes(G4Colour(0.0, 0.0, 1.0)); // blue
  // fLogicTarget->SetVisAttributes(VisAttTarget);   
  // G4VisAttributes* VisAttDetector= new G4VisAttributes(G4Colour(1.0, 1.0, 0.0)); // yellow
  // fLogicDetector->SetVisAttributes(VisAttDetector);
  // G4VisAttributes* VisAttKapWin= new G4VisAttributes(G4Colour(1.0, 1.0, 0.0)); // yellow
  // fLogicKapWin->SetVisAttributes(VisAttKapWin);
  // G4VisAttributes* VisAttBeamLineV= new G4VisAttributes(G4Colour(1.0,0.0,0.0)); // red
  // fLogicBeamLineV->SetVisAttributes(VisAttBeamLineV);          
    
  PrintParameters();
  
  //always return the root volume
  //
  return fPWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::PrintParameters()
{
  G4cout << "\n The World is " << G4BestUnit(fWorldSizeZ,"Length")
         << " of " << fWorldMater->GetName() 
         << "\n \n" << fWorldMater << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetMaterial(G4String materialChoice)
{
  if(materialChoice == "Vacuum") {
    fLogicTarget->SetMaterial(fVacuum);
  }
  else {
    // search the material by its name
    G4Material* pttoMaterial =
       G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);   
    
    if (pttoMaterial) { 
      if(fTargetMater != pttoMaterial) {
        fTargetMater = pttoMaterial;
        if(fLogicTarget) { fLogicTarget->SetMaterial(pttoMaterial); }
        G4RunManager::GetRunManager()->PhysicsHasBeenModified();
      }
    } else {
      G4cout << "\n--> warning from DetectorConstruction::SetMaterial : "
             << materialChoice << " not found" << G4endl;
    }    
  }
  
            
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetSize(G4double value)
{
  fWorldSizeZ = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

