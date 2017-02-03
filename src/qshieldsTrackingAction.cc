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
//
// $Id: TrackingAction.cc,v 1.2 2010/10/11 14:31:39 maire Exp $
// GEANT4 tag $Name: geant4-09-04 $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

#include "qshieldsTrackingAction.hh"

#include "G4Track.hh"
#include "G4TrackStatus.hh"
#include "G4Ions.hh"
#include "G4ParticleDefinition.hh"
#include "G4RunManager.hh"
#include "qshieldsConfigurator.hh"

using namespace CLHEP;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

qshieldsTrackingAction::qshieldsTrackingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

qshieldsTrackingAction::~qshieldsTrackingAction()
{}

void qshieldsTrackingAction::PreUserTrackingAction(const G4Track* track)
{
  G4ParticleDefinition* particle = track->GetDefinition();
  G4String name   = particle->GetParticleName();
  Charge = particle->GetPDGCharge();
  qshieldsConfigurator* theC = qshieldsConfigurator::GetInstance();
  G4bool isNucleus=false, isExcited=false, isSecondary=false;

  if (name == "neutron") {
    theC->IsNeutron = true;
  }
  else {
    theC->IsNeutron = false;
  }

  if( particle->GetParticleType() == "nucleus" ) {
	isNucleus = true;
	G4Ions* anIon =  reinterpret_cast<G4Ions*>(particle);
	if( anIon->GetExcitationEnergy() > 0.0 ) isExcited = true;
  } 

  if( isNucleus && isExcited ) theC->SetExcNucl(track->GetTrackID());
  if( track->GetParentID()!=0 ) isSecondary = true;  

  G4int isOk=1;
  if( name == "gamma" && !(theC->ParticleMask & 1 << theC->pId[ 'g' ])) isOk=0;
  else if( name == "neutron" && !(theC->ParticleMask & 1 << theC->pId[ 'n' ])) isOk=-1;
  else if( (name == "e+" ||  name == "e-") && !(theC->ParticleMask & 1 << theC->pId[ 'e' ])) isOk=-2;
  else if( name == "alpha" && !(theC->ParticleMask & 1 << theC->pId[ 'a' ])) isOk=-3;
  else if( name == "proton" && !(theC->ParticleMask & 1 << theC->pId[ 'p' ])) isOk=-4;
  else if( Charge>2 && isNucleus && isSecondary && !(theC->ParticleMask & 1 << theC->pId[ 'r' ])) isOk=-5;
// G4cout << "Tracking: " << name << " " << isOk << " " <<  isNucleus << " " << isSecondary << " " << isExcited << G4endl;
  if(theC->PrintInfo) G4cout << "Pre-tracking:  " << name << " " << track->GetVertexKineticEnergy() << " " << track->GetVertexMomentumDirection() << G4endl;
  if( isOk <= 0 ) {
// 	G4cout << "Tracking: " << name << " " << track->GetKineticEnergy() << " " << isOk << G4endl;
	G4Track* tr = (G4Track*) track;
  	tr->SetTrackStatus(fStopAndKill);
  }

  G4double decTime = particle->GetPDGLifeTime()/s;

  if ( !theC->GenDecAN && isSecondary ) {
		G4int aId = track->GetTrackID();
		G4int aParent = track->GetParentID();
		G4int aTest=aParent;
		while( (aTest=theC->GetMetastable(aTest))>0 ) ;
		if( aTest<0 ) theC->SetMetastable(aId,aParent);

 		if ( isNucleus && Charge>2 && isExcited && decTime>theC->TauLim ) { 
		     	G4double dTime = -decTime*log(1.-G4UniformRand());
				theC->SetMetastable(aId,dTime);
  		}
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void qshieldsTrackingAction::PostUserTrackingAction(const G4Track* track)
{
  qshieldsConfigurator* theC = qshieldsConfigurator::GetInstance();
  G4ParticleDefinition* particle = track->GetDefinition();
  G4String name   = particle->GetParticleName();
  Charge = particle->GetPDGCharge();
  G4bool isNucleus=false, isExcited=false, isSecondary=false, isAtRest=false;

//  if (Charge < 3. ) return;

  if(particle->GetParticleType() == "nucleus") {
	isNucleus = true;
	G4Ions* anIon =  reinterpret_cast<G4Ions*>(particle);
	if( anIon->GetExcitationEnergy() > 0.0 ) isExcited = true;
  } 

  if( track->GetParentID()!=0 ) isSecondary = true;  
  if( track->GetKineticEnergy() <=0.) isAtRest = true;
 
  if ( !theC->GenDecAN && isSecondary && isNucleus && Charge>2 && isAtRest && !isExcited ) {
   		 G4Track* tr = (G4Track*) track;
		 tr->SetTrackStatus(fKillTrackAndSecondaries);
		 theC->SetNucleus(Charge,particle->GetBaryonNumber(),particle->GetPDGStable(),track->GetPosition());
		 theC->dDaughterName = name;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
