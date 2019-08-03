#include <iostream>

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

#include "Detector.hh"

const double Detector::world_xy = 4*m;
const double Detector::world_z  = 200*m;


G4Material * Detector::GetMaterialByLocalName(G4String local_name){
  if (local_name == "vacuum_high"){
    return vacuum_high_;
  } else if (local_name == "vacuum_rough"){
    return vacuum_rough_;
  } else if (local_name == "air"){
    return air_;
  } else if (local_name == "argon_gas"){
    return argon_gas_;
  } else if (local_name == "argon_liquid"){
    return argon_liquid_;
  } else {
    G4cout << "ERROR: could not find material with local name " << local_name << "\n";
    return NULL;
  }
}

void Detector::SetNewValue(G4UIcommand* command, G4String arg){
  std::istringstream is((const char *) arg);

  if( command == fTargetMaterialCmd ){
    G4cout << "Setting target material to :  " << arg << "\n"; 
    target_material_ = GetMaterialByLocalName(arg);
    return;
  }
  if( command == fWorldMaterialCmd ){
    G4cout << "Setting world material to :  " << arg << "\n"; 
    world_material_ = GetMaterialByLocalName(arg);
    return;
  }
  if( command ==  fTzeroLocationCmd){
    is >> tzero_location_;
    tzero_location_ *= m;
    G4cout << "INFO:  Set tzero location " << arg << " parsed as " << tzero_location_ << "\n";
    return;
  }
  if( command ==  fDetectorEntranceCmd){
    is >> detector_entrance_;
    detector_entrance_ *= m;
    G4cout << "INFO:  Set detector entrance " << arg << " parsed as " << detector_entrance_ << "\n";
    return;
  }
  if( command ==  fTargetLengthCmd){
    is >> target_length_;
    target_length_ *= cm;
    G4cout << "INFO:  Set target length " << arg << " parsed as " << target_length_ << "\n";
    return;
  }
  if( command ==  fTargetRadiusCmd){
    is >> target_radius_;
    target_radius_ *= cm;
    G4cout << "INFO:  Set target radius " << arg << " parsed as " << target_radius_ << "\n";
    return;
  }
  if( command ==  fContainerRadiusCmd){
    is >> container_radius_;
    container_radius_ *= cm;
    G4cout << "INFO:  Set container radius " << arg << " parsed as " << container_radius_ << "\n";
    return;
  }
  if( command ==  fInsulationThicknessCmd){
    is >> insulation_thickness_;
    insulation_thickness_ *= cm;
    G4cout << "INFO:  Set insulation thickness " << arg << " parsed as " << insulation_thickness_ << "\n";
    return;
  }

  if( command ==  fWindowThicknessCmd){
    is >> window_thickness_;
    window_thickness_ *= cm;
    G4cout << "INFO:  Set window thickness " << arg << " parsed as " << window_thickness_ << "\n";
    return;
  }

  if( command ==  fTargetInCmd){
    target_in_ = (arg == "true");
    G4cout << "INFO:  Set target volume in flag " << arg << " parsed as " << target_in_ << "\n";
    return;
  }
  if( command ==  fContainerInCmd){
    container_in_ = (arg == "true");
    G4cout << "INFO:  Set target container in flag " << arg << " parsed as " << container_in_ << "\n";
    return;
  }
  if( command ==  fHallInCmd){
    hall_in_ = (arg == "true");
    G4cout << "INFO:  Set target hall in flag " << arg << " parsed as " << hall_in_ << "\n";
    return;
  }
}


Detector::Detector()
  :G4VUserDetectorConstruction(), G4UImessenger()
{ 
  G4bool broadcast = false;
  fDetDir = new G4UIdirectory("/artie/det/",broadcast);
  fDetDir->SetGuidance("detector construction commands");

  fWorldMaterialCmd = new G4UIcmdWithAString("/artie/det/world_material",this);
  fWorldMaterialCmd->SetGuidance("Select material of the world");
  fWorldMaterialCmd->SetParameterName("material",false);
        
  fTargetMaterialCmd = new G4UIcmdWithAString("/artie/det/target_material",this);
  fTargetMaterialCmd->SetGuidance("Select material of the target");
  fTargetMaterialCmd->SetParameterName("material",false);

  fTzeroLocationCmd = new G4UIcmdWithADouble("/artie/det/tzero_location",this);
  fTzeroLocationCmd->SetGuidance("Specify z position of neutron at t=0 (m)");
  fTzeroLocationCmd->SetParameterName("position",false);

  fDetectorEntranceCmd = new G4UIcmdWithADouble("/artie/det/detector_entrance",this);
  fDetectorEntranceCmd->SetGuidance("Specify z position of the detector entrance (m)");
  fDetectorEntranceCmd->SetParameterName("position",false);

  fTargetLengthCmd = new G4UIcmdWithADouble("/artie/det/target_length",this);
  fTargetLengthCmd->SetGuidance("Specify length of the target (cm)");
  fTargetLengthCmd->SetParameterName("length",false);

  fTargetRadiusCmd = new G4UIcmdWithADouble("/artie/det/target_radius",this);
  fTargetRadiusCmd->SetGuidance("Specify radius of the target inner volume (cm)");
  fTargetRadiusCmd->SetParameterName("radius",false);

  fContainerRadiusCmd = new G4UIcmdWithADouble("/artie/det/container_radius",this);
  fContainerRadiusCmd->SetGuidance("Specify outer radius of the target container (cm)");
  fContainerRadiusCmd->SetParameterName("radius",false);

  fInsulationThicknessCmd = new G4UIcmdWithADouble("/artie/det/insulation_thickness",this);
  fInsulationThicknessCmd->SetGuidance("Specify thickness of the insulation wrapped around target container (cm)");
  fInsulationThicknessCmd->SetParameterName("thickness",false);

  fWindowThicknessCmd = new G4UIcmdWithADouble("/artie/det/window_thickness",this);
  fWindowThicknessCmd->SetGuidance("Specify thickness of each window in the containment vessel (cm)");
  fWindowThicknessCmd->SetParameterName("thickness",false);

  fTargetInCmd = new G4UIcmdWithABool("/artie/det/target_in",this);
  fTargetInCmd->SetGuidance("Is the target volume in place?");
  fTargetInCmd->SetParameterName("in",false);

  fContainerInCmd = new G4UIcmdWithABool("/artie/det/container_in",this);
  fContainerInCmd->SetGuidance("Is the target container in place?");
  fContainerInCmd->SetParameterName("in",false);

  fHallInCmd = new G4UIcmdWithABool("/artie/det/hall_in",this);
  fHallInCmd->SetGuidance("Model the experimental hall??");
  fHallInCmd->SetParameterName("in",false);

  world_    = NULL;
  target_   = NULL;
  detector_ = NULL;

  DefineMaterials();  

  // set default values:
  world_material_ = air_;
  target_material_ = argon_liquid_;

  tzero_location_    = -1 * m;
  detector_entrance_ = 69 * m;

  target_length_ = 200 * cm;
  // Target Container Dimensions from Specs:
  // DN25 (ID = 25 mm), Flange size 2-1/8" (OD)
  // OD = 1-3/8"
  // (1 + 3./8.0) * 2.54 = 3.4925 cm
  target_radius_        = 2.50 * cm / 2.0;
  container_radius_     = 3.49 * cm / 2.0;
  insulation_thickness_ = 10.0 * cm;
  window_thickness_     = 0.00762*cm;

  target_in_ = true;
  container_in_ = true;
}

Detector::~Detector()
{ 
  delete fDetDir;
  delete fTargetMaterialCmd;
  //delete fTargetLengthCmd;
  //delete fTargetIn;
  //delete fTargetContainerIn;
}

void Detector::DefineMaterials()
{
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

  // Vacuum: (as low density air)
  // From ideal gas law:  P/density = CONST
  G4double density_rough = air_->GetDensity() * (0.1);
  G4double density_high  = air_->GetDensity() * (1.0E-11);
  vacuum_high_ = new G4Material("vacuum_high", density_high, 1);
  vacuum_high_->AddMaterial(air_, 1.0);
  vacuum_rough_ = new G4Material("vacuum_rough", density_rough, 1);
  vacuum_rough_->AddMaterial(air_, 1.0);    

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




void Detector::PrintParameters()
{
  G4cout << "World material is " << world_->GetMaterial()->GetName() << "\n";
  G4cout << "Target material is " << target_->GetMaterial()->GetName() << "\n";
  G4cout << "Target length is " << target_length_ / m << " m\n";
  G4double flight_path = detector_entrance_ - tzero_location_;
  G4cout << "Neutron Flight Path is " <<  flight_path / m << " m\n";

}


void Detector::ConstructHall(){
  static const double wall_thickness = 1*m;

  if (hall_in_){
    G4Box* outer = new G4Box("Wall_outer", world_xy/2, world_xy/2, world_z/2);
    G4Box* inner = new G4Box("Wall_inner", world_xy/2-wall_thickness, world_xy/2-wall_thickness, world_z/2-wall_thickness);
    G4SubtractionSolid * wall_s = new G4SubtractionSolid("Wall_s",outer,inner);
    G4LogicalVolume * wall_l = new G4LogicalVolume(wall_s, concrete_, "Wall_l");
    new G4PVPlacement(0, G4ThreeVector(), wall_l, "Wall_p", world_, false, 0);
  }

}

void Detector::ConstructBeamPipe(){
  // the beam pipe goes from negative infinity to the detector with a gap to accomodate the target
  static const double beampipe_radius_inner = 18*cm;
  static const double beampipe_radius_outer = 20*cm;

  static const double gap = 2.5 * m;  // gap to accommodate target
  double pos_a = -world_z / 2.0;
  double pos_b = -gap / 2.0;
  double pos_c = gap / 2.0;
  double pos_d = detector_entrance_;

  double left_halflength  = (pos_b - pos_a) / 2.0;
  double left_position    = (pos_a + pos_b) / 2.0; 
  double right_halflength  = (pos_d - pos_c) / 2.0;
  double right_position    = (pos_d + pos_c) / 2.0; 

  G4Tubs* left_beam_s = new G4Tubs("left_beam_s", 0, beampipe_radius_inner, left_halflength, 0.,CLHEP::twopi);   
  G4LogicalVolume* left_beam_l = new G4LogicalVolume(left_beam_s, vacuum_high_, "left_beam_l");                                    
  new G4PVPlacement(0, G4ThreeVector(0,0,left_position), left_beam_l, "left_beam_p", world_, false, 0);

  G4Tubs* right_beam_s = new G4Tubs("right_beam_s", 0, beampipe_radius_inner, right_halflength, 0.,CLHEP::twopi);   
  G4LogicalVolume* right_beam_l = new G4LogicalVolume(right_beam_s, vacuum_high_, "right_beam_l");                                    
  new G4PVPlacement(0, G4ThreeVector(0,0,right_position), right_beam_l, "right_beam_p", world_, false, 0);


  G4Tubs* left_pipe_s = new G4Tubs("left_pipe_s", beampipe_radius_inner, beampipe_radius_outer, left_halflength, 0.,CLHEP::twopi);   
  G4LogicalVolume* left_pipe_l = new G4LogicalVolume(left_pipe_s, vacuum_high_, "left_pipe_l");                                    
  new G4PVPlacement(0, G4ThreeVector(0,0,left_position), left_pipe_l, "left_pipe_p", world_, false, 0);

  G4Tubs* right_pipe_s = new G4Tubs("right_pipe_s", beampipe_radius_inner, beampipe_radius_outer, right_halflength, 0.,CLHEP::twopi);   
  G4LogicalVolume* right_pipe_l = new G4LogicalVolume(right_pipe_s, vacuum_high_, "right_pipe_l");                                    
  new G4PVPlacement(0, G4ThreeVector(0,0,right_position), right_pipe_l, "right_pipe_p", world_, false, 0);

}

void Detector::ConstructTarget(){
  static const double buffer_length = 5*cm; 

  G4cout << "target_in_:  " << target_in_ << "\n";
  if (target_in_){
    G4Tubs* solid = new G4Tubs("Target_s", 0, target_radius_, 0.5*target_length_, 0.,CLHEP::twopi);   
    target_ = new G4LogicalVolume(solid, target_material_, "Target_l");                                    
    new G4PVPlacement(0, G4ThreeVector(), target_, "Target_p", world_, false, 0);
  }                       
  if (container_in_){
    // The stainless steel containment vessel:
    G4Tubs* container_s  = new G4Tubs("Container_s", target_radius_, container_radius_, 0.5*target_length_, 0.,CLHEP::twopi);
    G4LogicalVolume* container_l = new G4LogicalVolume(container_s, stainless_, "Container_l");   
    new G4PVPlacement(0, G4ThreeVector(), container_l, "Container_p", world_, false, 0);

    // The wrapped insulation:
    G4Tubs* insulation_s = new G4Tubs("Insulation_s", container_radius_, container_radius_ + insulation_thickness_, 0.5*target_length_, 0.,CLHEP::twopi);
    G4LogicalVolume* insulation_l = new G4LogicalVolume(insulation_s, polyurethane_, "Insulation_l");   
    new G4PVPlacement(0, G4ThreeVector(), insulation_l, "Insulation_p", world_, false, 0);

    // The windows:
    G4Tubs* window_s = new G4Tubs("Window_s", 0, container_radius_, 0.5*window_thickness_, 0., CLHEP::twopi);
    G4LogicalVolume* window_l = new G4LogicalVolume(window_s, kapton_, "Window_l");
    new G4PVPlacement(0,G4ThreeVector(0,0, -(target_length_+window_thickness_)*0.5),window_l, "left_window_p", world_, false, 0);
    new G4PVPlacement(0,G4ThreeVector(0,0, +(target_length_+window_thickness_)*0.5),window_l, "rght_window_p", world_, false, 0);

    // The gas buffers:
    G4Tubs* buffer_s = new G4Tubs("Buffer_s", 0, target_radius_, 0.5*buffer_length, 0., CLHEP::twopi);
    G4LogicalVolume* buffer_l = new G4LogicalVolume(buffer_s, argon_gas_, "Buffer_l");
    new G4PVPlacement(0,G4ThreeVector(0,0, -(target_length_+2*window_thickness_+buffer_length)*0.5),buffer_l, "left_buffer_p", world_, false, 0);
    new G4PVPlacement(0,G4ThreeVector(0,0, +(target_length_+2*window_thickness_+buffer_length)*0.5),buffer_l, "rght_buffer_p", world_, false, 0);

  }
}

void Detector::ConstructDetector(){
  static const double detector_length  = 20.*cm;  
  static const double detector_radius  =  2.*cm;  // not known  
  double z_center = detector_entrance_ + detector_length * 0.5;

  G4Tubs* solid = new G4Tubs("Detector_s", 0, detector_radius, 0.5*detector_length, 0.,CLHEP::twopi);   
  detector_ = new G4LogicalVolume(solid, water_, "Detector_l");                                      
  new G4PVPlacement(0, G4ThreeVector(0., 0., z_center), detector_, "Detector_p", world_, false, 0);                                

}

G4VPhysicalVolume* Detector::Construct(){
  // Cleanup old geometry
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();


  G4Box* world_solid = new G4Box("Container", world_xy/2, world_xy/2, world_z/2);
  world_= new G4LogicalVolume(world_solid, world_material_, world_material_->GetName());
  G4VPhysicalVolume* physical_world =  
    new G4PVPlacement(0, G4ThreeVector(), world_, world_material_->GetName(), 0, false, 0);

  ConstructHall();
  ConstructBeamPipe();
  ConstructTarget();
  ConstructDetector();

  PrintParameters();

  return physical_world;
}







