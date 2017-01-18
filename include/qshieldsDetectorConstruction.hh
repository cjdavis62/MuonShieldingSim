// vers. 5.0: geant4 physics lists
// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: qshieldsDetectorConstruction.hh,v 5.0 2005/07/21 20:26:17 cremona Exp $
// GEANT4 tag $Name:  $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef qshieldsDetectorConstruction_h
#define qshieldsDetectorConstruction_h 1

#include "qshieldsSD.hh"
#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4Material.hh"
#include <vector>
#include <map>
#include "G4LogicalVolume.hh"
#include "G4Tubs.hh"
#include "G4PVPlacement.hh"
#include "G4UserLimits.hh"
#include "G4Region.hh"
#include "G4DistributedSource.hh"

class qshieldsDetectorConstruction : public G4VUserDetectorConstruction
{
public:
  
  qshieldsDetectorConstruction();
  ~qshieldsDetectorConstruction();
  
  G4VPhysicalVolume* Construct();
  G4double ZMount;
  std::map<std::string, G4Material*> MaterialMap;
  
private:
  static  qshieldsDetectorConstruction* theDetector;
  void DefineMaterials();
  void DumpMaterialsMap();
  void AddCuts(G4Region*, double);
  qshieldsSD* TSD;  //pointer to the sensitive detector
  qshieldsSD* CuSD;  //pointer to the copper sensitive detector
  qshieldsSD* PTFESD;  //pointer to the PTFE sensitive detector
  qshieldsSD* NTDSD;  //pointer to the NTD sensitive detector
  qshieldsSD* MuSD; //pointer to the Muon Shield sensitive detector
  qshieldsSD* OTHERSD;  //pointer to the PEN tape sensitive detector
  std::vector<G4VPhysicalVolume*> SD_phys; //pointer to the physical Absorber
  G4VPhysicalVolume* Hall_phys;    //pointer to the physical World

  G4UserLimits*    theUserLimitsForRoom; 
  G4UserLimits*    theUserLimitsForDetector; 

  G4double         theMaxTimeCuts;
  G4double         theMaxStepSize;
  G4double         theDetectorStepSize;
  G4double         theMinEkine;
  G4double         theRoomMinEkine;  
  G4double         theRoomTimeCut;
  double GammaCut;
  double ElectronCut;
  double PositronCut;

public:
   const G4VPhysicalVolume* GetphysiWorld() {return Hall_phys;};	    
   const std::vector<G4VPhysicalVolume*> GetSensitivePVolume() {return SD_phys;};
   G4int CrCounter;
//   G4int NumPlanes;
//   G4int NumCrystalsPlane;
   G4int NumTowers;
   G4int NumHoles;
   G4int NumHolesExt;
   G4int NumGiuntiCuBe;
   G4int NumScassiRS4;
   std::vector<G4ThreeVector> GiuntiCuBePos;
   std::vector<G4ThreeVector> ScassiRS4Pos;
   std::vector<G4ThreeVector> FrameZ;			
   std::vector<G4ThreeVector> CrystalC;
   std::vector<G4ThreeVector> ColumnC;
   std::vector<G4ThreeVector> NTDC;
   std::vector<G4ThreeVector> TowersC;
   std::vector<G4ThreeVector> WireTrayC;
   std::vector<G4ThreeVector> WireTrayPadC;
   std::vector<G4ThreeVector> WirePadC;
   std::vector<G4ThreeVector> PTFEC;
   std::vector<G4RotationMatrix> UpRot;
   std::vector<G4RotationMatrix> DoRot;
   std::vector<G4double> WireTrayRot;
   std::vector<G4ThreeVector> HoleP;
   std::vector<G4ThreeVector> HolePExt;
   std::vector<G4ThreeVector> PbRS4LatHoleP;
   std::vector<G4DistributedSource *> aSource;
   G4int SourceN;
   G4double totArea;
   G4double totVolume;
   G4ThreeVector CryoPos;
   std::vector< G4VPhysicalVolume * > fPhysicalVolumeVector;
   std::vector< G4Transform3D> fTransformVector;
};


#endif

