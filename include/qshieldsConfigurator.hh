// First version for geant4.7.1 
#ifndef qshieldsConfigurator_h
#define qshieldsConfigurator_h 1
#include "qshieldsDetectorConstruction.hh"
#include "G4ThreeVector.hh"
#include <vector>
#include <map>
#include <fstream>
#include <TTree.h>
#include <TFile.h>
#define MaxDist 1
#define TNumberOfParticles 12
#define TNumberOfInteractions 3
#define tNumberOfDetectors 988
#define MaxNumberOfCalls 200
#define n_shields 5

using namespace CLHEP;

class qshieldsConfigurator {

public:
  static qshieldsConfigurator* GetInstance();
  static qshieldsConfigurator* GetInstance(int argc, char** argv);
  ~qshieldsConfigurator();
  void SetEED( int meta, G4int n, G4double EED, G4int k ) { EventEnergyDeposit[n][k][meta] += EED; } 
  void SetCuED( int meta, G4double CuED ) { CuEnergyDeposit[meta] += CuED; } 
  void SetPTFEED( int meta, G4double PTFEED ) { PTFEEnergyDeposit[meta] += PTFEED; } 
  void SetNTDED( int meta, G4double NTDED ) { NTDEnergyDeposit[meta] += NTDED; }
  void SetMuED (int meta, G4double MuED ) { MuEnergyDeposit[meta] += MuED; }
  void SetOtherED( int meta, G4double OtherED ) { OtherEnergyDeposit[meta] += OtherED; } 
  void OutUpdate(void) { 
	if( OutType & 2) ROOTTree->AutoSave("SaveSelf");
  }
//-- Microphysics 
  void SetInt( G4int m, G4int aInt, G4int aCris, G4ThreeVector aPos, G4double aDep ) { 
  	if( !outMicro ) return;
	IntDep[aCris][numInt[aCris][m]][m] = aDep;
	IntPos[aCris][numInt[aCris][m]][m] = aPos;	
	numInt[aCris][m]++;	//-- total number of interactions
	if( FirstInteraction[aCris][m]<0 ) FirstInteraction[aCris][m] = aInt;
	nInt[aCris][aInt][m]++;	//-- interaction counter
  } 
  void SetCluster( G4int m, G4int aCris, G4ThreeVector aPos ) { 
	G4double dist,distM=1e20;
	G4int nCl;
	if( !outMicro ) return;
	for( G4int i=0;i<nClus[aCris][m];i++ ) {
	 if( (dist=aPos.diff2(cPos[aCris][i][m]))<distM ) { distM = dist; nCl=i; }
//	 G4cout << "       Configuratorhh  SetCluster test: " << i << " " << dist << " " << distM << " " << cPos[aCris][i] << G4endl;
	
	}
	if( nClus[aCris][m] == 0 || distM>MaxDist ) { nCl=nClus[aCris][m]; nClus[aCris][m]++; }
	else {
		cLen[aCris][nCl][0][m] += sqrt(distM);
		cLen[aCris][nCl][1][m] += distM;
		nLen[aCris][nCl][m]++;
	}
//	G4cout << "Configuratorhh  SetCluster: " << aCris << " " << nClus[aCris] << " " << nCl << " " << distM << " " << aPos << G4endl;
	if( aCris>=0 && aCris<tNumberOfDetectors && nCl>=0 && nCl<MaxNumberOfCalls ) cPos[aCris][nCl][m] = aPos; 
  } 
//-- End Microscopic
  
//-- Metastable states
  G4int GetMetastable(G4int aTest) {
	if( aTest>MaxStack ) return(0);
	else if( aTest>=0 ) return(StackId[aTest]); 
	else  return(0); 
  }

  void SetMetastable(G4int aId, G4int aParent) { 
	  if( aId>=MaxStack ) {
	  	for(G4int i=MaxStack; i<aId; i++ ) StackId.push_back(0);
		StackId.push_back(aParent);
	  	MaxStack = aId+1; 
		if( PrintInfo ) G4cout << "MaxStack updated (set): " << aId << " " << MaxStack << " " << StackId.size() << G4endl;
	  }
	  else StackId[aId] = aParent;
  }

  void SetMetastable(G4int aId,double Time) { 
	  if( aId>=MaxStack ) {
	  	for(G4int i=MaxStack; i<aId; i++ ) StackId.push_back(0);
		StackId.push_back(-1);
	  	MaxStack = aId+1; 
		if( PrintInfo ) G4cout << "MaxStack updated: " << aId << " " << MaxStack << " " << StackId.size() << G4endl;
	  }
	  else StackId[aId] = -1;
	  MetaTime=Time; 
  }
//-- End Metastable states

  void ClearED() { 
	for( G4int i=0; i<2; i++ ) {
     CuEnergyDeposit[i] = 0.;
     PTFEEnergyDeposit[i] = 0.;
     NTDEnergyDeposit[i] = 0.;
     MuEnergyDeposit[i] = 0;
     OtherEnergyDeposit[i] = 0.;
    }
    NumberOfCalls = 0; 
	MetaTime = 0.;
	for( G4int i=0; i<MaxStack; i++ ) StackId[i] = 0;
     for( G4int i=0; i<tNumberOfDetectors; i++ ) 
       for( G4int k=0; k<TNumberOfParticles; k++ )
         for( G4int l=0; l<2; l++ )
           EventEnergyDeposit[i][k][l] = 0.0; 
	 ClearRootVariables();
 	 if( outMicro ) {
	     for( G4int i=0; i<tNumberOfDetectors; i++ ) {
		   for( G4int l=0; l<2; l++ ) {
			for( G4int k=0; k<TNumberOfInteractions; k++ ) nInt[i][k][l] = 0; 
			nClus[i][l] = 0; 
			numInt[i][l] = 0;
			FirstInteraction[i][l] = -1;
			for( G4int k=0; k<MaxNumberOfCalls; k++ ) {
			   nLen[i][k][l]=0;
		       for( G4int j=0; j<TNumberOfInteractions; j++ ) cLen[i][j][k][l]=0.;
			}
		   }
		 }	 	
 	 }
  }

  void ClearRootVariables() { 
	// clear ROOT variables
	  	dEnergy=0.;
	  	for( G4int i=0; i<TNumberOfParticles; i++ ) dPEnergy[i] = 0.;

  } 
  void IncreaseNoC() { NumberOfCalls++; }
  G4int GetNoC() { return NumberOfCalls; }
  void SetCDC( qshieldsDetectorConstruction* det ) { CDet = det; }
  void ResetChain() { dChainNumber++; PTime[0] = 0.; LastChainTime = 0;}
  void ResetChain( G4int fIdx) { dChainNumber=fIdx; PTime[0] = 0.; LastChainTime = 0; }
  void SetChainEvent(G4int fId) { if( fId==0 || SingleDecayGranularity ) ResetChain(); ChainEvent = fId;}
  void SetChainTime( double pt ) { 
	PTime[0] += pt; 
  }
  void SetCurrentParticle( int cpt ) {CurrentParticle = cpt;} 
  void SetMuSphereCenter( G4ThreeVector aVector ) { BottomCenter = aVector; }
  
  void SetParticle( double pe, int pq, int pm, int pz ); 
  void SetParticle( double pe, int pq, int pm, int pz, double px ); 
  void SetParticle( double pe, int pq, int pm ) {
    n_particles = 1;
    PEnergy[0] = pe*MeV;
    PType[0] = pq; 
    PMass[0] = static_cast<double> ( pm );
  }
  void SetNucleus( G4int fZ, G4int fA, G4bool fS, G4ThreeVector fP ){AtomicNumber=fZ; MassNumber=fA; StableNucleus=fS; LastPosition=fP;}
  int GetType( int typ );  
  void SetType( char typ, double energ ); 
  void PrintStats();
  void SetExcNucl(G4int fk) { ExcNucl[fk] = true; }
  G4bool GetExcNucl(G4int fk) { return ExcNucl[fk]; }
  void ResetExcNucl() { for(G4int i=0;i<1000;i++) ExcNucl[i]=false; }

  unsigned long long int dChainNumber;
  long	   RSeed;
  G4int	   SAF;
  G4int	   LogNumber;
  G4int	   SType;	// Source element: 0=PL, <> Surf./Vol. nth element
  unsigned int	   SMode;	// Int.(<)/Tot(0)/Ext(>) Surface
  G4int Mode;
  G4int OutType;
  G4ThreeVector SAngle;
  std::vector<G4ThreeVector> SWDirCos;
  std::vector<G4ThreeVector> SPosition;
  G4double SGNorm;
  G4double UniverseR;
  G4double EventEnergyDeposit[tNumberOfDetectors][TNumberOfParticles][2];
  G4int nInt[tNumberOfDetectors][TNumberOfInteractions][2];
  G4int numInt[tNumberOfDetectors][2];
  G4ThreeVector IntPos[tNumberOfDetectors][MaxNumberOfCalls][2];
  G4double IntDep[tNumberOfDetectors][MaxNumberOfCalls][2];
  G4int nClus[tNumberOfDetectors][2]; 
  G4double cLen[tNumberOfDetectors][MaxNumberOfCalls][2][2]; 
  G4int nLen[tNumberOfDetectors][MaxNumberOfCalls][2]; 
  G4ThreeVector cPos[tNumberOfDetectors][MaxNumberOfCalls][2];
  double CuEnergyDeposit[2];
  double PTFEEnergyDeposit[2];
  double NTDEnergyDeposit[2];
  double MuEnergyDeposit[2];
  double OtherEnergyDeposit[2];
  G4double SOffs;
  G4double SDepth;
  G4double LastChainTime;
  std::vector<double> PEnergy;
  std::vector<double> PTime;
  std::vector<double> PMass;
  std::vector<double> PZeta;
  std::vector<double> PEexc;
  std::vector<int> PType;
  G4double ChainE;
  G4double ChainA;
  G4double ChainZ;
  G4double ChainX;
  G4ThreeVector BottomCenter;

  G4String PhysListName;
  G4int CurrentParticle;
  bool GenDec;
  bool GenDecAN;
  bool Debug;
  G4int v1, v2;
  G4int DebugEvents;
  G4int DoubleBeta;
  G4double QBeta;
  G4int n_particles;
  G4int NumCycles;
  G4int ParticleMask;
  std::string MacroFile;
  std::string OutFile[2];
  std::string CfgFile;  
  std::string NumberOfEvents;
  G4String GendecFile;
  G4String SpecInpFile;
  G4String SpecOutFile;
  bool xyz_out;
  G4int iSpec, oSpec;
  G4int NumberOfLoops;
  G4int TNumberOfDetectors;
  G4int NPs; 
  G4int NPb; 
  qshieldsDetectorConstruction* CDet;
  std::map< char,G4int > pId;
  G4int ExtShieldMask;
  G4double CAngle;
  G4double WiR,WiL,WiZ;
  bool PrintInfo;
  bool PrintGendec;
  G4double BaseCut;
  G4int SurfExpGen;
  std::ofstream DataFile;
  std::ofstream SpecFile;
  std::ofstream SourceFile;
  G4int VolumeNumber,ModuleNumber;
  G4int ggCorrelation;
  G4double ThermalMax;
  G4double EThermalMax;
  G4double EDistrUniform;
  G4double EThermalAvr;
  G4double StepFactor;
  bool neutron_gen;
  G4ThreeVector MuPar;
   bool qshieldsTest;
  G4int NumqshieldsTest;
  G4double GenDecThresh;
  G4int GenerationZone;
  G4int partyp;
  G4double parten;
  G4double tThresh;
  G4int ChainEvent;
  bool SingleDecayGranularity;
  G4double TauLim;
  G4int AtomicNumber;
  G4int MassNumber;
  G4int StableNucleus;
  G4int AStop;
  G4int ZStop;
  G4ThreeVector LastPosition;
  G4bool ExcNucl[1000];
  G4bool outMicro;

  std::vector<G4int> StackId;
  G4int MaxStack;
  double MetaTime;

// DBD from file
  std::string fileE;
  std::string fileA;

///////////////////////
// ROOT TTree output stuff
  TFile *ROOTFile;
  TTree *ROOTTree;
  TTree *CreateROOTTree();
  TTree *INFOtree;
  TTree *SOURCETree;

// OUTPUT N-PLE VARIABLES
  int dChannel;
  double dDTime;
  double dEnergy;
  double dPEnergy[TNumberOfParticles];
  double dEParticle;
  double dDirectionX;
  double dDirectionY;
  double dPositionX;
  double dPositionY;
  double dPositionZ;
  double dDepth;
  double dCuEnergyDeposit;
  double dPTFEEnergyDeposit;
  double dNTDEnergyDeposit;
  double dMuEnergyDeposit;
  double dOtherEnergyDeposit;
//-- Microscopic variables
  int dNumInt;
  double dAvDist;
  int dNInt[TNumberOfInteractions];
  int dNClus;
  double dTLen;
  int FirstInteraction[tNumberOfDetectors][2];
  G4String dFirstInteraction;
  std::string dInteractionName;
  std::string dParticleName;
  std::string dDaughterName;

// HEADER VARIABLES ////
  std::string command;
  unsigned long long int iNtot; //number of simulated events
///////////////////////

////// qshields (CUORE) peculiar   section
  std::string TowerCenter;
  int singleTower;
  int NumberOfPlanes;
  G4int	   SElem;	// Source element: 0=PL, <> Surf./Vol. nth element
  G4double SphereRadius;
  G4int GraphicMask[2];
  G4double minEn;  //S.C.

  G4int NumTotCrystals;

private:
  void CuoreMap ( std::string, bool );

// end of qshields section

private:
  void qshieldsHelp(void);
  std::vector<G4double>  qshieldsConvert( std::string, G4int );
  std::vector<std::string> qshieldsSplit( std::string, G4int );
  qshieldsConfigurator();
  void qshieldsReleaseNotes();
  qshieldsConfigurator(int argc, char** argv);
  static  qshieldsConfigurator* theConfiguration;
  void setparameters();
  void RemoveWhiteSpaces ( std::string& );
  G4int LocatePlusSign ( std::string );
  G4int NumberOfCalls;
};


#endif
