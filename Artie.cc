
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "Randomize.hh"
#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "PrimaryGenerator.hh"
#include "Analysis.hh"

#ifdef G4VIS_USE
 #include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

int main(int argc,char** argv) {
 
  //choose the Random engine
  G4Random::setTheEngine(new CLHEP::RanecuEngine);

  // Initialize:  detector, physics, generator, analysis
  G4RunManager* runManager = new G4RunManager;
  runManager->SetUserInitialization(new DetectorConstruction());  
  runManager->SetUserInitialization(new PhysicsList());  
  runManager->SetUserInitialization(new PrimaryGenerator::Init()); 
  runManager->SetUserInitialization(new Analysis::Init()); 
  
  // get the pointer to the User Interface manager 
  G4UImanager* UI = G4UImanager::GetUIpointer();  

#ifdef G4VIS_USE
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
#endif
  if (argc!=1)   // batch mode  
    {
     G4String command = "/control/execute ";
     G4String fileName = argv[1];
     UI->ApplyCommand(command+fileName);
    }    
  else           //define visualization and UI terminal for interactive mode
    { 
#ifdef G4VIS_USE
      G4UIExecutive * ui = new G4UIExecutive(argc,argv);  
      //visManager->Initialize();
#endif    
     UI->ApplyCommand("/control/execute vis.mac");  
#ifdef G4UI_USE
      //G4UIExecutive * ui = new G4UIExecutive(argc,argv);      
      ui->SessionStart();
      delete ui;
#endif
          
#ifdef G4VIS_USE
     delete visManager;
#endif     
    }

  delete runManager;

  return 0;
}


