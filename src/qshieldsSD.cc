// --------------------------------------------------------------
//      GEANT 4 - qshields 
//
//      Simulation of CUORICINO experiment
//      Oliviero Cremonesi, Silvia Capelli, Luca Gironi  
//	Milano Bicocca 2009 
// --------------------------------------------------------------

#include "qshieldsSD.hh"

#include "qshieldsConfigurator.hh"
#include "qshieldsDetectorConstruction.hh"
#include "G4ProcessManager.hh"

#include "G4VPhysicalVolume.hh"
#include "G4Step.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4SDManager.hh"

#include "G4ios.hh"

qshieldsSD::qshieldsSD(const G4String name, qshieldsDetectorConstruction* )
:G4VSensitiveDetector(name)
{
  G4String HCname="SDCollection";
  collectionName.insert(HCname);
}

qshieldsSD::~qshieldsSD()
{
}

void qshieldsSD::Initialize(G4HCofThisEvent*)
{
}

G4bool qshieldsSD::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
  G4StepPoint* endPoint = aStep->GetPostStepPoint();
  G4String procName = endPoint->GetProcessDefinedStep()->GetProcessName();
  G4String particleName = aStep->GetTrack()->GetParticleDefinition()->GetParticleName();
  G4ThreeVector Position,sPosition;
  G4int num,mer;
  G4double edep;
  qshieldsConfigurator* theC = qshieldsConfigurator::GetInstance();
 
  G4int meta=0;
  if( theC->GetMetastable(aStep->GetTrack()->GetTrackID()) !=0 ) meta=1;

  G4TouchableHistory* theTouch
    = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());

//-- microsphysics section
  if( theC->outMicro ) {
	  Position = aStep->GetPostStepPoint()->GetPosition();
	  if( Position!= sPosition ) {
		  if( (mer=procName.find ( "phot" ))  >= 0 ) num=0;
		  else if( (mer=procName.find ( "compt" ))  >= 0 ) num=1;
		  else if( (mer=procName.find ( "conv" ))  >= 0 ) num=2;
		  else num=-1;

	  	  if( theTouch->GetVolume()->GetName() == "TeO2" ) {
			if( num >=0 ) {	//-- Photon interactions
				edep = fabs( aStep->GetPreStepPoint()->GetTotalEnergy() - aStep->GetPostStepPoint()->GetTotalEnergy() );
	 			theC->SetInt( meta, num, theTouch->GetReplicaNumber(), Position, edep );
		    }
			if( aStep->GetTotalEnergyDeposit() > 0 && procName.find( "eIoni" )  > 0 ) theC->SetCluster( meta, theTouch->GetReplicaNumber(), Position ) ;

		  }
		  sPosition = Position;
	  }
  }
//-- End micro
  
  if( theC->PrintInfo && num >=0 ) 
		G4cout <<  "qshieldsSD: " << theC->dChainNumber << " " << theTouch->GetReplicaNumber() << " " << num << " " 
		<< procName << " " << Position << " " << theC->FirstInteraction[theTouch->GetReplicaNumber()] << " "
		<< aStep->GetPreStepPoint()->GetTotalEnergy() << " " << aStep->GetPostStepPoint()->GetTotalEnergy() << " "
		<< edep << G4endl;
	
  edep = aStep->GetTotalEnergyDeposit();

  if (edep==0.) return false;      

  G4TouchableHistory* theTouchable
    = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());

  G4int p_type=11;

  G4String pname = aStep->GetTrack()->GetDefinition()->GetParticleName();
  if( pname == "gamma" ) p_type = 0;
  else if( pname == "e-" ) p_type = 1;
  else if( pname == "e+" ) p_type = 2;
  else if( pname == "proton" ) p_type = 3;
  else if( pname == "anti_proton" ) p_type = 4;
  else if( pname == "neutron" ) p_type = 5;
  else if( pname == "anti_neutron" ) p_type = 6;
  else if( pname == "alpha" ) p_type = 7;
  else if( pname == "mu-" ) p_type = 9;
  else if( pname == "mu+" ) p_type = 10;
  else if( aStep->GetTrack()->GetDefinition()->GetParticleType() == "nucleus" ) p_type = 8;

  if( theTouchable->GetVolume()->GetName() == "TeO2_P" )  theC->SetEED( meta, theTouchable->GetReplicaNumber(), edep, p_type );
  else if ( theTouchable->GetVolume()->GetName() == "PTFE_P" ) theC->SetPTFEED( meta, edep );
  else if ( theTouchable->GetVolume()->GetName() == "NTD_P" ) theC->SetNTDED( meta, edep );
  else if ( theTouchable->GetVolume()->GetName() == "CuDetUpPlate_P" 
	  ||  theTouchable->GetVolume()->GetName() == "CuDetUpPlate_P" 
	  ||  theTouchable->GetVolume()->GetName() == "CuDetBottomPlate_P" 
	  ||  theTouchable->GetVolume()->GetName() == "CuRS6Tiles_P" 
	  ||  theTouchable->GetVolume()->GetName() == "FrameC_P" 
	  ||  theTouchable->GetVolume()->GetName() == "FrameT_P" 
	  ||  theTouchable->GetVolume()->GetName() == "FrameB_P" 
	  ||  theTouchable->GetVolume()->GetName() == "Column_P" ) theC->SetCuED( meta, edep );
// PENTape WireTray WirePad
  else  theC->SetOtherED( meta, edep );

  return true;
}

void qshieldsSD::EndOfEvent(G4HCofThisEvent*)
{
}

void qshieldsSD::clear()
{} 

void qshieldsSD::DrawAll()
{} 

void qshieldsSD::PrintAll()
{} 

