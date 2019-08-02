#ifndef Detector_h
#define Detector_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4UImessenger.hh"
#include "globals.hh"

class G4LogicalVolume;
class G4Material;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADouble;
class G4UIcmdWithABool;

class Detector : public G4VUserDetectorConstruction, public G4UImessenger
{
public:
  Detector();
  ~Detector();

  // Detector Constructor Interface:
  virtual G4VPhysicalVolume* Construct();
  // UI Messenger Interface (for setting parameters):
  virtual void SetNewValue(G4UIcommand*, G4String);
  
  const G4LogicalVolume* GetLogicWorld()     const {return world_;};
  const G4LogicalVolume* GetLogicDetector()  const {return detector_;};
  const G4LogicalVolume* GetLogicTarget()    const {return target_;};
  
  void               PrintParameters();

  // location of neutrons at t=0 (hook provided for primary generator)
  double TzeroLocation() const {return tzero_location_; };


  // fixed values:
  static const double world_xy;
  static const double world_z;

private:
  // Logical Volumes:
  G4LogicalVolume * world_;
  G4LogicalVolume * target_;
  G4LogicalVolume * detector_;
    
  // Defined Materials 
  G4Material * air_; 
  G4Material * water_; 
  G4Material * argon_liquid_;
  G4Material * argon_gas_;
  G4Material * stainless_;
  G4Material * aluminum_;
  G4Material * graphite_;
  G4Material * vacuum_rough_;
  G4Material * vacuum_high_;
  G4Material * lipoly_;
  G4Material * polyurethane_;
  G4Material * kapton_;
  G4Material * concrete_;

  // obtain above materials by our own local friendly names:
  G4Material * GetMaterialByLocalName(G4String local_name);

  // Selectable Materials
  G4Material * world_material_;
  G4Material * target_material_;

  // Configurable Geometery:
  G4double   tzero_location_;       // z position of neutron at t=0
  G4double   detector_entrance_;    // z position of detector entrance
  G4double   target_length_;        // target length in meters
  G4double   target_radius_;        // radius of target ( = inner radius of cotainer)
  G4double   container_radius_;     // out radius of container 
  G4double   insulation_thickness_; // thickness of the foam insulation
  G4double   window_thickness_;     // thickness of each window in containment vessel
  G4bool     target_in_;            // is the target volume in?
  G4bool     container_in_;         // is the target container included with target?
  G4bool     hall_in_;              // include model for experimental hall?

  // UI commands:
  // - include command "/control/manual /artie/det" into .mac file for a detailed description of each command:
  G4UIdirectory*             fDetDir;
  G4UIcmdWithAString*        fWorldMaterialCmd;
  G4UIcmdWithAString*        fTargetMaterialCmd;
  G4UIcmdWithADouble*        fTzeroLocationCmd;
  G4UIcmdWithADouble*        fDetectorEntranceCmd;
  G4UIcmdWithADouble*        fTargetLengthCmd;
  G4UIcmdWithADouble*        fTargetRadiusCmd;  
  G4UIcmdWithADouble*        fContainerRadiusCmd;  
  G4UIcmdWithADouble*        fInsulationThicknessCmd;
  G4UIcmdWithADouble*        fWindowThicknessCmd; 
  G4UIcmdWithABool*          fTargetInCmd;
  G4UIcmdWithABool*          fContainerInCmd;
  G4UIcmdWithABool*          fHallInCmd;

  // Divide Construction Into Steps:
  void DefineMaterials();  
  void ConstructHall();
  void ConstructBeamPipe();
  void ConstructTarget();
  void ConstructDetector();
           
};

#endif
