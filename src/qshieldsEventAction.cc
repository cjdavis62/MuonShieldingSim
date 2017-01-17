// vers. 5.0: geant4 physics lists
#include "qshieldsEventAction.hh"
#include "qshieldsConfigurator.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"
#include "Randomize.hh"

qshieldsEventAction::qshieldsEventAction() {
}


qshieldsEventAction::~qshieldsEventAction() {
}


void qshieldsEventAction::BeginOfEventAction(const G4Event* evt)
{
 qshieldsConfigurator* theC = qshieldsConfigurator::GetInstance();
 theC->ClearED();

 G4int evtNb = evt->GetEventID() + 1;
 
 if (evtNb%theC->LogNumber == 0)
   { 
    G4cout << "Event: " << evtNb << " [" << theC->dChainNumber <<  "]" <<G4endl;
	theC->OutUpdate();
   }
    
}


void qshieldsEventAction::EndOfEventAction(const G4Event* evt)
{
  qshieldsConfigurator* theC = qshieldsConfigurator::GetInstance();

  theC->PrintStats();

  G4TrajectoryContainer * trajectoryContainer = evt->GetTrajectoryContainer();
  G4int n_trajectories = 0;
  if(trajectoryContainer)
    { n_trajectories = trajectoryContainer->entries(); }

  if(G4VVisManager::GetConcreteInstance())
    {
      for(G4int i=0; i<n_trajectories; i++)  { 
	G4Trajectory* trj = (G4Trajectory *)((*(evt->GetTrajectoryContainer()))[i]);
	trj->DrawTrajectory();
      }
    }
}


