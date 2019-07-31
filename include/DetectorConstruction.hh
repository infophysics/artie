#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4UImessenger.hh"
#include "globals.hh"

class G4LogicalVolume;
class G4Material;
class G4UIdirectory;
class G4UIcommand;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;

class DetectorConstruction : public G4VUserDetectorConstruction, public G4UImessenger
{
public:
  DetectorConstruction();
  ~DetectorConstruction();

  // Detector Constructor Interface:
  virtual G4VPhysicalVolume* Construct();
  // UI Messenger Interface (for setting parameters):
  virtual void SetNewValue(G4UIcommand*, G4String);
  
  void SetSize     (G4double);              
  void SetMaterial (G4String);            

  const G4VPhysicalVolume* GetWorld()      {return fPWorld;};           
                    
  G4double           GetSize()       {return fWorldSizeX;}; //usless function (J.Wang)     
  G4Material*        GetMaterial()   {return fWorldMater;};
     
  const G4LogicalVolume* GetLogicWorld()           {return fLWorld;};   
  const G4LogicalVolume* GetLogicDetector()        {return fLogicDetector;};  
  const G4LogicalVolume* GetLogicTarget()        {return fLogicTarget;};
  const G4LogicalVolume* GetLogicLArcontainer()        {return fLogicLArcontainer;}; 
  const G4LogicalVolume* GetLogicBeamLineV()     {return fLogicBeamLineV;};   
  
  void               PrintParameters();
  
private:
  
  // inch to cm
  G4double           fInch2cm = 2.54;
  
  // Vacuum
  G4Material* fVacuum;
     
  // world 
  G4VPhysicalVolume* fPWorld;
  G4LogicalVolume*   fLWorld;
  G4double           fWorldSizeX;
  G4double           fWorldSizeY;
  G4double           fWorldSizeZ;
  G4Material*        fWorldMater;   
  
  //Room
  G4VPhysicalVolume* fPhysiRoom;
  G4LogicalVolume*   fLogicRoom;
  G4double           fRoomSizeX;
  G4double           fRoomSizeY;
  G4double           fRoomSizeZ;
  G4Material*        fRoomMater;
  
  //Room Volume
  G4VPhysicalVolume* fPhysiRoomVolume;
  G4LogicalVolume*   fLogicRoomVolume;
  G4double           fRoomThickness;
  
     
  // Insulator container
  G4double           fInsulatorContainerLength; 
  G4double           fInsulatorContainerInnerRadius; 
  G4double           fInsulatorContainerOuterRadius; 
  G4LogicalVolume*   fLogicInsulatorContainer;
  G4VPhysicalVolume* fPhysiInsulatorContainer;
  //G4Element*         fInsulatorContainerMater;
  G4Material*        fInsulatorContainerMater; 
  
  //Insulator
  G4double           fInsulatorLength; 
  G4double           fInsulatorInnerRadius; 
  G4double           fInsulatorOuterRadius; 
  G4LogicalVolume*   fLogicInsulator;
  G4VPhysicalVolume* fPhysiInsulator;
  G4Material*        fInsulatorMater; 
  
  // liquid argon container   
  G4double           fLArcontainerLength; 
  G4double           fLArcontainerInnerRadius; 
  G4double           fLArcontainerOuterRadius; 
  G4LogicalVolume*   fLogicLArcontainer;
  G4VPhysicalVolume* fPhysiLArcontainer;
  G4Material*        fLArcontainerMater; 
     
  // liquid argon target
  G4double           fTargetLength; 
  G4double           fTargetRadius; 
  G4LogicalVolume*   fLogicTarget;
  G4VPhysicalVolume* fPhysiTarget;
  G4Material*        fTargetMater; 
  
  // kapton window
  G4double           fKaptonThickness;
  G4Material*        fkapton;  
  G4LogicalVolume*   fLogicKapWin;
  
  // kapton window R1
  G4VPhysicalVolume* fPhysiKapWinR1;
  
  // kapton window R2 
  G4VPhysicalVolume* fPhysiKapWinR2;
  
  // kapton window L1 
  G4VPhysicalVolume* fPhysiKapWinL1;
  
  // kapton window L2
  G4VPhysicalVolume* fPhysiKapWinL2;
  
  //Fill and vent lines
  G4double           fPipeRadius;
  G4double           fPipeLength;
  //G4RotationMatrix*   PipeRot;
  //Fill Line
  G4LogicalVolume*   fLogiclArFillLineIN;
  G4VPhysicalVolume* fPhysilArFillLineIN;
  G4LogicalVolume*   fLogiclArFillLineOUT;
  G4VPhysicalVolume* fPhysilArFillLineOUT;
  //Vent Line
  G4LogicalVolume*   fLogiclArVentLineIN;
  G4VPhysicalVolume* fPhysilArVentLineIN;
  G4LogicalVolume*   fLogiclArVentLineOUT;
  G4VPhysicalVolume* fPhysilArVentLineOUT;
  
  // Vacuum line
  G4LogicalVolume*   fLogicVacuumLine;
  G4VPhysicalVolume* fPhysiVacuumLine;
  
  //Beam Line
  //G4double           fBeamLineRadius;
  G4double           fBeamLineLength;
  G4LogicalVolume*   fLogicBeamLine;
  G4VPhysicalVolume* fPhysiBeamLine;
  G4Material*        fBeamLineMater;
  
  //Beam Line Volume
  G4double           fBeamLineRadiusOUT;
  G4double           fBeamLineRadiusIN;
  G4LogicalVolume*   fLogicBeamLineV;
  G4VPhysicalVolume* fPhysiBeamLineV;
  G4Material*        fBeamLineVolumeMater;
  
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
  
  // buffer container
  G4double						fBufferLength;
  G4double						fBufferInnerRadius;
  G4double						fBufferOuterRadius;
  G4LogicalVolume*   fLogicBuffer;
  G4Material*        fBufferMater; 
  
  G4VPhysicalVolume* fPhysiBufferL;
  G4VPhysicalVolume* fPhysiBufferR;
  
  // buffer volume
  G4LogicalVolume*   fLogicBufferVol;
  G4VPhysicalVolume* fPhysiBufferVolL;
  G4VPhysicalVolume* fPhysiBufferVolR;
  
  
  G4UIdirectory*             fTestemDir;
  G4UIdirectory*             fDetDir;
  G4UIcmdWithAString*        fMaterCmd;
  G4UIcmdWithADoubleAndUnit* fSizeCmd;
  G4UIcommand*               fIsotopeCmd;    
  
  void               DefineMaterials();
  G4VPhysicalVolume* ConstructVolumes();    
           
};

#endif
