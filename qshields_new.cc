// --------------------------------------------------------------
//      GEANT 4 - qshields 
//
//      Simulation of MIBETA (20 TeO2 Array) experiment
//      Oliviero Cremonesi - 2-Dec-99
// --------------------------------------------------------------
// Comments
// vers. 5.0: geant4 physics lists
// 
// --------------------------------------------------------------

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"
#include "G4UIXm.hh"
#include "G4UIExecutive.hh"
#include "Randomize.hh"
#include "G4Gendec.hh"
#include "qshieldsDebug.hh"

#ifdef G4VIS_USE
#include "qshieldsVisManager.hh"
#endif

#include "qshieldsDetectorConstruction.hh"
#ifdef G4GT10
#include "G4PhysListFactory.hh"
#else
#include "PhysicsList.hh"
#include "Shielding.hh"
#endif
#include "qshieldsPrimaryGeneratorAction.hh"
#include "qshieldsRunAction.hh"
#include "qshieldsEventAction.hh"
//#include "qshieldsStackingAction.hh"
#include "Randomize.hh"
#include "qshieldsConfigurator.hh"
#include "qshieldsTrackingAction.hh"
#include "G4VisExecutive.hh"
#include "G4Version.hh"

#include "TStyle.h"

#define MAX_NP	13

int main(int argc,char** argv) {

  
  // Initializes Configurator (Reads command line options)
  qshieldsConfigurator* theConfigurator = qshieldsConfigurator::GetInstance(argc, argv);

  if( theConfigurator->OutType & 1) G4cout << "Text output File : " << theConfigurator->OutFile[0] << G4endl;
  if( theConfigurator->OutType & 2) G4cout << "Root output File : " << theConfigurator->OutFile[1] << G4endl;
  G4cout << "Number of generated particles/chains : "   << " " << theConfigurator->n_particles << G4endl;


  G4String nev(theConfigurator->NumberOfEvents.c_str());
  G4String macro(theConfigurator->MacroFile.c_str());
   
  // Initialise Gendec .. if required
  if( theConfigurator->GenDecAN ) {
    G4Gendec* theGendec = G4Gendec::GetInstance();
    if( theGendec->dec_init_( theConfigurator->GendecFile.c_str(),false ) < 0 ) {
     G4cout << "Error initialising Gendec! " << theConfigurator->GendecFile << G4endl
     << "*** ABORT! *** " << G4endl;
     exit(1);
    }
  }
   
  // Construct the default run manager
  G4RunManager * runManager = new G4RunManager;

  // set mandatory initialization classes
  qshieldsDetectorConstruction* detector = new qshieldsDetectorConstruction;
  theConfigurator->SetCDC( detector );
  runManager->SetUserInitialization(detector);

  if ( theConfigurator->Mode == 3)
  {
	detector->Construct();
//	qshieldsDebug* debug;
//	debug = new qshieldsDebug();
	return 0;
  }

#ifdef G4GT10
  G4PhysListFactory *aPhysicsList = new G4PhysListFactory();
  runManager->SetUserInitialization(aPhysicsList->GetReferencePhysList(theConfigurator->PhysListName));
#else
//  PhysicsList *aPhysicsList = new PhysicsList(theConfigurator->neutron_gen,false);
//  PhysicsList *aPhysicsList = new PhysicsList(true,false);
  PhysicsList *aPhysicsList = new PhysicsList();
  if (theConfigurator->PhysListName == "shielding") {
    runManager->SetUserInitialization (new Shielding);
  }
  else {
    runManager->SetUserInitialization(aPhysicsList);
  }
  if( theConfigurator->StepFactor > 0. ) aPhysicsList->SetCutScale(theConfigurator->StepFactor);
  aPhysicsList->SelectPhysicsList("Livermore_EM");
  if(theConfigurator->PhysListName == "livermore") aPhysicsList->SelectPhysicsList("Livermore_EM");
  else if(theConfigurator->PhysListName == "penelope") aPhysicsList->SelectPhysicsList("Penelope_EM");
  if(theConfigurator->PhysListName == "standard") aPhysicsList->SelectPhysicsList("Standard_EM");
  if( theConfigurator->GenDecAN ) aPhysicsList->killRadioactiveDecays();
#endif
  
  // choose the Random engine
  CLHEP::HepRandom::setTheEngine(new CLHEP::HepJamesRandom( theConfigurator->RSeed ));
  G4cout << "Random Seed: " << CLHEP::HepRandom::getTheSeed() << G4endl;
  
  // set user action classes
  runManager->SetUserAction(new qshieldsPrimaryGeneratorAction());
  runManager->SetUserAction(new qshieldsRunAction());
  runManager->SetUserAction(new qshieldsEventAction);
  runManager->SetUserAction(new qshieldsTrackingAction);
  //  runManager->SetUserAction(new qshieldsStackingAction());
  
  //Initialize G4 kernel
  runManager->Initialize();
  
  // visualization manager
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();

  G4UImanager * UImanager = G4UImanager::GetUIpointer();

  if (theConfigurator->Mode == 0)   // Define UI session for interactive mode.
    {
      G4cout << "Mode == 0" << G4endl;

      G4UIExecutive * ui = new G4UIExecutive(argc,argv);

      // G4UIterminal is a (dumb) terminal.
      /*
#ifdef G4UI_USE_TCSH
      G4cout << "Actually this one" << G4endl;
      G4UIsession * ui = new G4UIterminal(new G4UItcsh);
#else
      G4cout << "Doing this one" << G4endl;
      G4UIsession * ui = new G4UIterminal();
#endif
      */
      ui->SessionStart();
      delete ui;
    }
    

  else           // Batch mode
      { 
	if ( theConfigurator->Mode == 1 ) {

	  UImanager->ApplyCommand("/control/execute qshieldsPrerun.mac");    
	  G4String command = "/control/execute ";
	  UImanager->ApplyCommand(command+macro);
      }
	else if ( theConfigurator->Mode == 2 ) {
	  UImanager->ApplyCommand("/control/execute qshieldsPrerun.mac");    	
	  if ( theConfigurator->OutType & 1 ) {
	    theConfigurator->DataFile.open ( theConfigurator->OutFile[0].c_str() );
	    for( G4int i=0; i<argc; i++ ) theConfigurator->DataFile << argv[i] << " ";
	    theConfigurator->DataFile << G4endl;
	  }
	  if ( theConfigurator->OutType & 2 ) theConfigurator->ROOTFile = new TFile( theConfigurator->OutFile[1].c_str(), "RECREATE");

	  for( G4int i=0; i<theConfigurator->NumberOfLoops; i++ ) {
	    G4cout <<  "Batch " << i << " of " << nev << " events completed ..." << G4endl;
	    UImanager->ApplyCommand("/run/beamOn "+nev);    
	  }
	  if( theConfigurator->OutType & 1 ) 
	    theConfigurator->DataFile.close ();
	  if( theConfigurator->OutType & 2 ) {
	    // Write and close ROOT File
	    theConfigurator->INFOtree->Fill();
	    theConfigurator->INFOtree->Write();
	    theConfigurator->ROOTTree->Write();
	    theConfigurator->ROOTFile->Close();  
	  }
	}
      }



  // job termination
#ifdef G4VIS_USE
  delete visManager;
#endif
  delete runManager;

  return 0;
}

