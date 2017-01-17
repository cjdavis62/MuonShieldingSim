// First version for geant4.7.1 
// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: qshieldsPrimaryGeneratorAction.hh,v 5.1 2005/08/04 09:31:05 cremona Exp $
// GEANT4 tag $Name:  $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef qshieldsPrimaryGeneratorAction_h
#define qshieldsPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include <vector>
#include <fstream>
#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4GenericIon.hh"
#include "G4GSMu.hh"
#include "G4VisExtent.hh"
#include "G4Transform3D.hh"
#include "G4Gendec.hh"

struct UNSISO{
  G4int A, Z;
  G4double Exc;
  G4int id, nsucc;
  G4double Weight;
  std::vector<G4int> succ;
};

class G4ParticleGun;
class G4Event;
class qshieldsDetectorConstruction;
class qshieldsPrimaryGeneratorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class qshieldsPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
    qshieldsPrimaryGeneratorAction();    
   ~qshieldsPrimaryGeneratorAction();
    void GeneratePrimaries(G4Event*);
    G4double Thermals(G4double x, G4double y) { return sqrt(x)*exp(-x/y); }
    G4ThreeVector GetAlphaDir( G4ThreeVector fDirection, G4int fId ){  
	  G4ThreeVector Direction;
	  if( GenAlpha == fId ) { // previous alpha decay
	       Direction = -DirAlpha;
	       GenAlpha = 0;
	  } else {
		   Direction = fDirection;
		   DirAlpha = fDirection;
	       GenAlpha = -fId;
	  }
	  return( Direction);
    }
    void ResetAlphaDir( ){ GenAlpha = 0; }

private:
  G4ParticleGun* particleGun;	  //pointer a to G4 service class

  G4ThreeVector GetPointInVolume ( const G4VisExtent VisExt );
  G4ThreeVector GetRandomDirection();
  G4ThreeVector GetRandomDirection( G4ThreeVector SA, G4ThreeVector P );
  G4ThreeVector GetRandomDirection( G4ThreeVector olDir );
  G4ThreeVector GetRandomShellPosition (G4double Rmin, G4double Rmax, G4double Alt);
  void GenDoubleBeta (G4int type, G4double* EBeta, G4ThreeVector* A, G4double Q);
  void DoubleBetaDistr (G4int type, G4double Q, std::string filen, std::string filan);
  G4ThreeVector GenPosition, DirAlpha, PreviousPosition;
  G4int GenAlpha;
  G4int k_particle,b_particle;
  G4double PreviousChainTime;
  std::vector<G4double> vX;
  std::vector<G4double> vY;
  std::vector<struct UNSISO> Chain;
  G4int nSpec;
  G4GSMu *aMuon;
  G4double qMu;
  G4int nDecay;
  G4int pDecay;
  bool FixPos;
  G4double wtmp;
  G4double dTime;
  G4double pE;
  G4double MuRadius;
  struct DEC_STACKS ds;
// params for Iachello-Kotila 2nu
  G4double F12Max;
  G4double **g0, estep0, emax0, **g1, estep1, emax1;
  G4int nstep0, nstep1;
};

#endif


