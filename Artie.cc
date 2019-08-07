
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "Randomize.hh"
#include "Detector.hh"
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
  // Detect interactive mode (if no arguments) and define UI session
  //
  G4UIExecutive* ui = 0;
  if ( argc == 1 ) {
    ui = new G4UIExecutive(argc, argv);
  }
 
  //choose the Random engine
  G4Random::setTheEngine(new CLHEP::RanecuEngine);

  // Initialize:  detector, physics, generator, analysis
  G4RunManager* runManager = new G4RunManager;
  runManager->SetUserInitialization(new Detector());  
  runManager->SetUserInitialization(new PhysicsList());  
  runManager->SetUserInitialization(new PrimaryGenerator::Init()); 
  runManager->SetUserInitialization(new Analysis::Init()); 

  // Initialize visualization
  //
  //G4VisManager* visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();

  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  // Process macro or start UI session
  if ( ! ui ) { 
    // batch mode
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);
  }
  else { 
    ui->SessionStart();
    delete ui;
  }

  delete runManager;
  delete visManager;

  return 0;
}


