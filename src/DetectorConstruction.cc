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
#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithABool.hh"

#include "DetectorConstruction.hh"

G4Material * DetectorConstruction::GetMaterialByLocalName(G4String local_name){
  if (local_name == "vacuum_rough"){
    return vacuum_rough_;
  } else if (local_name == "argon_liquid"){
    return argon_liquid_;
  } else {
    G4cout << "ERROR: could not find material with local name " << local_name << "\n";
    return NULL;
  }
}

void DetectorConstruction::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if( command == fTargetMaterialCmd ){
    G4cout << "Setting target material to :  " << newValue << "\n"; 
    target_material_ = GetMaterialByLocalName(newValue);
  }
}


DetectorConstruction::DetectorConstruction()
  :G4VUserDetectorConstruction(), G4UImessenger(), 
   fPWorld(0), fLWorld(0), fWorldMater(0), 
   fDetDir(0)
{ 
  DefineMaterials();  

  G4bool broadcast = false;
  fDetDir = new G4UIdirectory("/artie/det/",broadcast);
  fDetDir->SetGuidance("detector construction commands");
        
  fTargetMaterialCmd = new G4UIcmdWithAString("/artie/det/target_material",this);
  fTargetMaterialCmd->SetGuidance("Select material of the target (Vacuum or Argon)");
  fTargetMaterialCmd->SetParameterName("material",false);
  
  fWorldSizeX = 4*m;
  fWorldSizeY = 4*m;
  fWorldSizeZ = 200*m;

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
  
  // Buffer volume
  fBufferLength = 5*cm;    
  fBufferInnerRadius = fLArcontainerInnerRadius;
  fBufferOuterRadius = fLArcontainerOuterRadius;

  //Beam Line
  fBeamLineLength = fDetectorPositionZ - 0.5*fDetectorLength - 0.5*fInsulatorContainerLength- fBufferLength;
  fBeamLineRadiusOUT = 20.*cm;
  fBeamLineRadiusIN = 18.*cm;
}

DetectorConstruction::~DetectorConstruction()
{ 
  delete fDetDir;
  delete fTargetMaterialCmd;
  //delete fTargetLengthCmd;
  //delete fTargetIn;
  //delete fTargetContainerIn;
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  return ConstructVolumes();
}



void DetectorConstruction::DefineMaterials()
{
  const G4double ROOM_TEMP = 273 * kelvin;

  G4NistManager* man = G4NistManager::Instance();
  G4Element* C = man->FindOrBuildElement("C");
  G4Element* H = man->FindOrBuildElement("H");
  G4Element* N = man->FindOrBuildElement("N");
  G4Element* O = man->FindOrBuildElement("O");
  G4Element* Li = man->FindOrBuildElement("Li");
  air_         = man->FindOrBuildMaterial("G4_AIR");
  water_       = man->FindOrBuildMaterial("G4_WATER");
  argon_gas_   = man->FindOrBuildMaterial("G4_Ar");
  argon_liquid_ = man->FindOrBuildMaterial("G4_lAr");
  stainless_    = man->FindOrBuildMaterial("G4_STAINLESS-STEEL");
  aluminum_    = man->FindOrBuildMaterial("G4_Al");
  //graphite_  = man->FindOrBuildMaterial("G4_GRAPHITE");
  graphite_    = man->FindOrBuildMaterial("G4_GRAPHITE_POROUS");
  concrete_    = man->FindOrBuildMaterial("G4_CONCRETE");
  kapton_      = man->FindOrBuildMaterial("G4_KAPTON");

  // Define Custom Materials
  // Vacuum   --  NEEDS UPDATE TO ROUGH VACUUM
  G4double atomicNumber = 1.;
  G4double massOfMole = 1.008*g/mole;
  G4double density = 1.e-25*g/cm3;
  G4double pressure = 3.e-18*pascal;
  // these will likely differ, but setting the same for now:
  vacuum_rough_ = new G4Material("Vacuum", atomicNumber, massOfMole, density, kStateGas, ROOM_TEMP, pressure);
  vacuum_beam_  = vacuum_rough_;

 //Lithium Polyethylene
 G4Material* polyethylene = man->FindOrBuildMaterial("G4_POLYETHYLENE");
 lipoly_ = new G4Material("LiPoly", 1.06*g/cm3, 2);
 lipoly_->AddElement (Li, 7.54*perCent);
 lipoly_->AddMaterial (polyethylene, 92.46*perCent);

  // Polyurethane (foam insulator)
  polyurethane_ = new G4Material("Polyurethane", 1.005*g/cm3, 4);
  polyurethane_->AddElement(C, 3);
  polyurethane_->AddElement(H, 8);
  polyurethane_->AddElement(N, 2);
  polyurethane_->AddElement(O, 1);
	  
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}



G4VPhysicalVolume* DetectorConstruction::ConstructVolumes()
{

  G4cout << "Building detector with target material " << target_material_->GetName() << "\n";




    // Obsolete:  hard-code materials to defined materials or use a defined switchable material...
  // world mater
  fWorldMater = air_;//

  //Room Wall Mater
  fRoomMater = concrete_;
  
  // insulator container
  fInsulatorContainerMater = stainless_;
  
  // insulator
  fInsulatorMater = polyurethane_; 
  
  // LAr container
  fLArcontainerMater = stainless_; 
  	
  // neutron collimator
  fCollimatorMater = lipoly_;
  
  // neutron detector
  fDetectorMater = water_;

  // kapton
  fkapton  = kapton_;

  //Beam Line
  fBeamLineMater = stainless_;

  //Beam Line Volume
  fBeamLineVolumeMater = vacuum_beam_;
  
  // Buffer volume
  fBufferMater = argon_gas_; 

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
                                      target_material_,                              
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



void DetectorConstruction::PrintParameters()
{
  G4cout << "\n The World is " << G4BestUnit(fWorldSizeZ,"Length")
         << " of " << fWorldMater->GetName() 
         << "\n \n" << fWorldMater << G4endl;
}








