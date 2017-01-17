// vers. 5.0: geant4 physics lists
#include "qshieldsDetectorConstruction.hh"
#include "qshieldsConfigurator.hh"
#include "G4SDManager.hh"
#include "G4Material.hh"
#include "G4Isotope.hh"
#include "G4Box.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4Cons.hh" 
#include "G4Torus.hh" 
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4ProductionCuts.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4BooleanSolid.hh"
#include "G4VSolid.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include <vector>

#include "G4ios.hh"
#define n_shields 5

qshieldsDetectorConstruction::qshieldsDetectorConstruction():TSD(0),CuSD(0),PTFESD(),NTDSD(0),OTHERSD(0) {

  theUserLimitsForRoom     = 0; 
  theUserLimitsForDetector = 0; 
  // default time cut = infinite
  //  - note also number of steps cut in stepping action = MaxNoSteps
  theMaxTimeCuts      = DBL_MAX;
  theMaxStepSize      = DBL_MAX;
  theDetectorStepSize = DBL_MAX;
  theRoomTimeCut      = 1000. * nanosecond;
  theMinEkine         = 250.0*eV; // minimum kinetic energy required in volume
  theRoomMinEkine     = 250.0*eV; // minimum kinetic energy required in volume

}

qshieldsDetectorConstruction::~qshieldsDetectorConstruction() { 

 delete theUserLimitsForRoom;
 delete theUserLimitsForDetector;
 MaterialMap.clear();

}

G4VPhysicalVolume* qshieldsDetectorConstruction::Construct()
{
const double PbCutFactor=1.;
const double CuCutFactor=0.1;
const double FeCutFactor=0.1;

qshieldsConfigurator* theC = qshieldsConfigurator::GetInstance();

GammaCut    = 1*cm*theC->StepFactor;
ElectronCut = 1*cm*theC->StepFactor;
PositronCut = 1*cm*theC->StepFactor;
G4Region* PbRegion = new G4Region("Lead_region");
PbRegion->SetProductionCuts(new G4ProductionCuts());
PbRegion->GetProductionCuts()->SetProductionCut(GammaCut*PbCutFactor,"gamma");
PbRegion->GetProductionCuts()->SetProductionCut(ElectronCut*PbCutFactor,"e-");
PbRegion->GetProductionCuts()->SetProductionCut(PositronCut*PbCutFactor,"e+");

G4Region* CuRegion = new G4Region("Copper_region");
CuRegion->SetProductionCuts(new G4ProductionCuts());
CuRegion->GetProductionCuts()->SetProductionCut(GammaCut*CuCutFactor,"gamma");
CuRegion->GetProductionCuts()->SetProductionCut(ElectronCut*CuCutFactor,"e-");
CuRegion->GetProductionCuts()->SetProductionCut(PositronCut*CuCutFactor,"e+");

G4Region* FeRegion = new G4Region("Steel_region");
FeRegion->SetProductionCuts(new G4ProductionCuts());
FeRegion->GetProductionCuts()->SetProductionCut(GammaCut*FeCutFactor,"gamma");
FeRegion->GetProductionCuts()->SetProductionCut(ElectronCut*FeCutFactor,"e-");
FeRegion->GetProductionCuts()->SetProductionCut(PositronCut*FeCutFactor,"e+");

CrCounter=0;
G4int Cel = 0, Elm = 0;
G4VPhysicalVolume** TMPV;
G4VPhysicalVolume** TMPZ;
G4VPhysicalVolume** TMPK;
G4ThreeVector OVCC; 
G4RotationMatrix Rot;
G4Transform3D RTT;

  DefineMaterials();
  DumpMaterialsMap();

G4ThreeVector POFF; // Centro dell'ultimo volume considerato

G4ThreeVector PC; // Centro dell'ultimo volume considerato
TMPV = new G4VPhysicalVolume* [30];  
TMPZ = new G4VPhysicalVolume* [4];
TMPK = new G4VPhysicalVolume* [2];

	if( theC->SType !=0 )  Cel = abs( theC->SType );
	if( theC->SElem >0 )  Elm = abs( theC->SElem );

  // Hall mother volume

G4double HallR;

//				 
// Sensitive Detectors
//
  G4SDManager* SDman = G4SDManager::GetSDMpointer();

  if(!TSD)
  {
    TSD = new qshieldsSD("TSD",this);
    SDman->AddNewDetector( TSD );
  } 


  if(!CuSD)
  {
    CuSD = new qshieldsSD("CuSD",this);
    SDman->AddNewDetector( CuSD );
  } 
  if(!PTFESD)
  {
    PTFESD = new qshieldsSD("PTFESD",this);
    SDman->AddNewDetector( PTFESD );
  } 
  if(!NTDSD)
  {
    NTDSD = new qshieldsSD("NTDSD",this);
    SDman->AddNewDetector( NTDSD );
  } 
  if(!OTHERSD)
  {
    OTHERSD = new qshieldsSD("OTHERSD",this);
    SDman->AddNewDetector( OTHERSD );
  } 

//				 
// Universe
//
  HallR = theC->UniverseR;
  if( HallR <=0. )  HallR = 20*m;
  if( theC->SphereRadius > HallR ) HallR=2.*theC->SphereRadius;
  HallR += 10. * cm;

  G4Sphere * Hall 
   = new G4Sphere("Hall",
		 0., HallR,
		 0.0, 2.*M_PI,
		 0.0, M_PI );
  G4LogicalVolume * Hall_log 
    = new G4LogicalVolume (Hall, MaterialMap["Air"], "Hall_L", 0,0,0);

  G4VPhysicalVolume * Hall_phys
    = new G4PVPlacement(0,G4ThreeVector(),"Hall_P", 
			Hall_log, 0, false, 0);

  Hall_log->SetVisAttributes ( new G4VisAttributes(false) );
  PC = G4ThreeVector();  // Universe Center


#include "geom_params.icc"
// *** Everything is outside the cryostat ***
#include "qshieldsShieldings.icc"
// *** Everything is inside the cryostat ***
#include "qshieldsCryostat.icc"	    
// *** Detector ***
#include "qshieldsDetector.icc"

  if( Cel != 0 ) {
	for( G4int is=0;is < (int) aSource.size();is++) 
		G4cout << "Source element " << aSource[is]->getSolid()->GetName() << 
			" - Volume: " << aSource[is]->getVolume()/1000. << " cm3 - Area: " << aSource[is]->getArea()/100. << " cm2" << G4endl;

	for( G4int is=1;is < (int) aSource.size();is++) {
		aSource[is]->addVolume(  aSource[is-1]->getVolume() );
		aSource[is]->addArea(  aSource[is-1]->getArea() );
	}
	SourceN = aSource.size();
	totArea = aSource[SourceN-1]->getArea();
	totVolume = aSource[SourceN-1]->getVolume();
	G4cout << "Total source elements - Volume: " << totVolume/1000. << " cm3 - Area: " << totArea/100. << " cm2" << G4endl;
  }
 
  // ......................................................................
  // attach user limits ...................................................

  G4cout << G4endl << "User Limits: " << G4endl 
	 << "\t theMaxTimeCuts:     " << G4BestUnit(theMaxTimeCuts,"Time")  
	 << G4endl
	 << "\t theRoomTimeCut:     " << G4BestUnit(theRoomTimeCut,"Time")  
	 << G4endl
	 << "\t theMaxStepSize:     " << G4BestUnit(theMaxStepSize,"Length")
	 << G4endl
	 << "\t theMinEKine:        " << G4BestUnit(theMinEkine,"Energy")   
	 << G4endl
	 << "\t minRoomMinEKine:    " << G4BestUnit(theRoomMinEkine,"Energy")
	 << G4endl << G4endl;

  if (theUserLimitsForRoom != 0) delete theUserLimitsForRoom;
  if (theUserLimitsForDetector != 0) delete theUserLimitsForDetector;
  return Hall_phys;
}

void qshieldsDetectorConstruction::DefineMaterials(){

  G4double a, iz, z, density, fractionmass;
  G4String name = ""; 
  G4String symbol = "";
  G4int natoms;

  G4Material* Poly;
  G4Material* PolyB;
  G4Material* BAcid;
  G4Material* Air;
  G4Material* Tellurite;
  G4Material* CdWO4;
  G4Material* Copper;
  G4Material* Aluminum;
  G4Material* Germanium;
  G4Material* Lead;
  G4Material* LHe;
  G4Material* Iron;
  G4Material* Teflon;
  G4Material* Vacuum;
  G4Material* Silicon;
  G4Material* Gold;
  G4Material* Cadmium;
  G4Material* Steel;
  G4Material* Tungsten;
  G4Material* CuBe;
  G4Material* ZnMoO4;
  G4Material* ZnSe;
  G4Material* Mylar;
  G4Material* GranSassoRock;
  G4Material* Polyvinyltoluene;
 
//------------------- Elements -----------------------------

  // N & O & Cd

  a = 1.00794*g/mole;
  G4Element* elH = new G4Element(name="Hydrogen", symbol="H", iz=1., a);
  a = 10.811*g/mole;
  G4Element* elB = new G4Element(name="Boron", symbol="B", iz=5., a);
  a = 12.011*g/mole;
  G4Element* elBe = new G4Element(name="Berillium",symbol="Be", iz=4., a );
  a = 10.811*g/mole;
  G4Element* elC = new G4Element(name="Carbon", symbol="C", iz=6., a);
  a = 14.01*g/mole;
  G4Element* elN = new G4Element(name="Nitrogen", symbol="N", iz=7., a);
  a = 16.00*g/mole;
  G4Element* elO = new G4Element(name="Oxygen", symbol="O", iz=8., a);
  a = 18.9984032*g/mole;
  G4Element* elF = new G4Element(name="Fluorine", symbol="F", iz=9., a);
  a = 55.847*g/mole;
  G4Element* elFe = new G4Element(name="Iron", symbol="Fe", iz=26., a);
  a=63.55*g/mole;
  G4Element* elCu = new G4Element(name="Copper", symbol="Cu", iz=29., a);
  a = 112.411*g/mole;
  G4Element* elCd = new G4Element(name="Cadmium",symbol="Cd", iz=48., a );
  a = 127.60*g/mole;
  G4Element* elTe = new G4Element(name="Tellurium", symbol="Te", iz=52., a);
  a = 183.84*g/mole;
  G4Element* elW = new G4Element(name="Tungsten",symbol="W", iz=74., a );
  a = 26.981539*g/mole;
  G4Element* elAl = new G4Element(name="Aluminum",symbol="Al", iz=13., a );
  a = 40.0784*g/mole;
  G4Element* elCa = new G4Element(name="Calcium", symbol="Ca", iz=20., a);
  a = 28.085*g/mole;
  G4Element* elSi = new G4Element(name="Silicon", symbol="Si", iz=14., a);
  a = 24.305*g/mole;
  G4Element* elMg = new G4Element(name="Magnesium", symbol="Mg", iz=12., a);
  a = 39.9083*g/mole;
  G4Element* elK = new G4Element(name="Potassium", symbol="K", iz=19., a);
  
//------------------- Materials -----------------------------
  // Copper

  density = 8.960*g/cm3;
  a = 63.55*g/mole;
  z = 29.0;
  Copper = new G4Material(name="Copper"   , z , a, density );


  // Iron
  
  density = 7.87*g/cm3;
  a = 55.845*g/mole;
  z = 26.0;
  Iron = new G4Material(name="Iron"   , z , a, density );
  
  
  // Germanium

  density = 5.323*g/cm3;
  a = 72.61*g/mole;
  z = 32.0;
  Germanium = new G4Material(name="Germanium"   , z , a, density );

  // Liquid Helium

  density = 0.125*g/cm3;
  a = 4.00*g/mole;
  z = 2.0;
  LHe = new G4Material(name="LHe"   , z , a, density );

  // Aluminum

  density = 2.70*g/cm3;
  a = 26.981539*g/mole;
  z = 13.0;
  Aluminum = new G4Material(name="Aluminum"   , z , a, density );

  // Cadmium

  density = 8.9938*g/cm3;
  a = 112.411*g/mole;
  z = 48.0;
  Cadmium = new G4Material(name="Cadmium"   , z , a, density );

  // Lead

  density = 11.35*g/cm3;
  a = 207.19*g/mole;
  z = 82.0;
  Lead = new G4Material(name="Lead"   , z , a, density );

  // Vacuum

  density     = universe_mean_density;    //from PhysicalConstants.h
  G4double pressure    = 3.e-18*pascal;
  G4double temperature = 2.73*kelvin;
  Vacuum   = new G4Material("Vacuum", 1., 1.01*g/mole, density,kStateGas, temperature, pressure);
   
   
  // Gold

  density = 19.32*g/cm3;
  a = 196.9665*g/mole;
  z = 79.0;
  Gold = new G4Material(name="Gold"   , z , a, density );                            

    //Silicon
   
   density = 2.33*g/cm3;
   a = 28.0855*g/mole;
   z = 14.0;
   Silicon = new G4Material(name="Silicon"   , z , a, density );
   
   // Tungsten
  
  density = 19.25*g/cm3;
  a = 183.84*g/mole;
  z = 74.0;
  Tungsten = new G4Material(name="Tungsten"   , z , a, density );

//------------------- Compounds -----------------------------

  // Air
  density = 1.29e-03*g/cm3;
  Air = new G4Material(name="Air", density, 2);
  Air->AddElement(elN, natoms=4);
  Air->AddElement(elO, natoms=1);

  // Steel
  density = 7.874*g/cm3;
  Steel = new G4Material(name="Steel", density, 2);
  Steel->AddElement(elFe, fractionmass=0.96);
  Steel->AddElement(elC, fractionmass=0.04);

  // Multi-layer AlMylar
  density = 1.3*g/cm3;
  Mylar = new G4Material(name="Mylar", density, 4);
  Mylar->AddElement(elC, fractionmass=0.597);
  Mylar->AddElement(elO, fractionmass=0.045);
  Mylar->AddElement(elH, fractionmass=0.354);
  Mylar->AddElement(elAl, fractionmass=0.004);
  
  // Tellurite
  density = 6.0*g/cm3;
  Tellurite = new G4Material(name="Tellurite", density, 2);
  Tellurite->AddElement(elTe, natoms=1);
  Tellurite->AddElement(elO, natoms=2);
  
  // Teflon

  density = 2.2*g/cm3;
  Teflon = new G4Material(name="Teflon", density, 2);
  Teflon->AddElement(elC, natoms=2);
  Teflon->AddElement(elF, natoms=4);
   
  // Poly

  density = 0.935*g/cm3;
  Poly = new G4Material(name="Poly", density, 2);
  Poly->AddElement(elC, natoms=2);
  Poly->AddElement(elH, natoms=4);
  
  // PolyB

  density=0.95*g/cm3;
  PolyB = new G4Material(name="Borated poly", density, 2);
  PolyB->AddElement(elB, fractionmass=5*perCent);
  PolyB->AddMaterial(Poly, fractionmass=95*perCent);
 
  // BAcid //S.C.

  density=1.435*g/cm3;
  BAcid = new G4Material(name="Boric Acid", density, 3);
  BAcid->AddElement(elH, natoms=3);
  BAcid->AddElement(elB, natoms=1);
  BAcid->AddElement(elO, natoms=3);

  // CuBe

  density = 8.2*g/cm3; //
  CuBe = new G4Material(name="CuBe", density, 2);
  CuBe->AddElement(elCu, fractionmass=0.98);
  CuBe->AddElement(elBe, fractionmass=0.02);

  // ZnSe,5.27,1,2,Zinc,Zn,30,65.38,1,Selenium,Se,34,78.96,1

  G4Element* elZn = new G4Element(name="Zinc",symbol="Zn", iz=30, a=65.38*g/mole );
  G4Element* elSe = new G4Element(name="Selenimu",symbol="Se", iz=34, a=78.96*g/mole );
  density = 5.27*g/cm3; //
  ZnSe = new G4Material(name="ZnSe", density, 2);
  ZnSe->AddElement(elZn, natoms=1);
  ZnSe->AddElement(elSe, natoms=1);

  // ZnMoO4,4.3,1,3,Zinc,Zn,30,65.38,1,Molybdenum,Mo,42,95.96,1,Oxygen,O,8,15.999,4
  G4Element* elMo = new G4Element(name="Molybdenum",symbol="Mo", iz=42, a=95.96*g/mole );
  density = 4.3*g/cm3; //
  ZnMoO4 = new G4Material(name="ZnMoO4", density, 3);
  ZnMoO4->AddElement(elZn, natoms=1);
  ZnMoO4->AddElement(elMo, natoms=1);
  ZnMoO4->AddElement(elO, natoms=4);

  // CdWO4,7.9,1,3,Cadmium,Cd,48,112.41,1,Tungsten,W,74,183.84,1,Oxygen,O,8,15.999,4
  density = 7.9*g/cm3; //
  CdWO4 = new G4Material(name="CdWO4", density, 3);
  CdWO4->AddElement(elCd, natoms=1);
  CdWO4->AddElement(elW, natoms=1);
  CdWO4->AddElement(elO, natoms=4);

  // Gran Sasso Rock Hall A
  density = 2.71 * g/cm3;
  GranSassoRock = new G4Material(name="GranSassoRock", density, 7);
  GranSassoRock->AddElement(elC, fractionmass=12.00*perCent);
  GranSassoRock->AddElement(elO, fractionmass=48.40*perCent);
  GranSassoRock->AddElement(elMg, fractionmass=5.64*perCent);
  GranSassoRock->AddElement(elAl, fractionmass=1.04*perCent);
  GranSassoRock->AddElement(elSi, fractionmass=1.28*perCent);
  GranSassoRock->AddElement(elK, fractionmass=1.04*perCent);
  GranSassoRock->AddElement(elCa, fractionmass=30.60*perCent);

  // Polyvinyltoluene
  density = 1.023 * g/cm3;
  Polyvinyltoluene = new G4Material(name="Polyvinyltoluene", density, 2);
  Polyvinyltoluene->AddElement(elC, natoms=9);
  Polyvinyltoluene->AddElement(elH, natoms=10);

  // Build materials hash table, usefull for the dynamic assignament of the material

  MaterialMap["Air"] = Air;
  MaterialMap["Tellurite"] = Tellurite;
  MaterialMap["Copper"] = Copper;
  MaterialMap["Aluminum"] = Aluminum;
  MaterialMap["Germanium"] = Germanium;
  MaterialMap["Lead"] = Lead;
  MaterialMap["LHe"] = LHe;
  MaterialMap["Teflon"] = Teflon;
  MaterialMap["Iron"] = Iron;
  MaterialMap["Vacuum"] = Vacuum;
  MaterialMap["Silicon"] = Silicon;
  MaterialMap["Gold"] = Gold;
  MaterialMap["Cadmium"] = Cadmium;
  MaterialMap["CdWO4"] = CdWO4;
  MaterialMap["Poly"] = Poly;
  MaterialMap["PolyB"] = PolyB;
  MaterialMap["BAcid"] = BAcid; 
  MaterialMap["Steel"] = Steel;
  MaterialMap["Tungsten"] = Tungsten;
  MaterialMap["Mylar"] = Mylar;
  MaterialMap["CuBe"] = CuBe;
  MaterialMap["GranSassoRock"] = GranSassoRock;
  MaterialMap["Polyvinyltoluene"] = Polyvinyltoluene;
}
  
void qshieldsDetectorConstruction::DumpMaterialsMap(){

  G4cout << "Material map : " << G4endl;
  G4cout << G4endl << "     Name                  Density                    Temperature            Pressure  " << G4endl;
  G4cout << "-----------------------------------------------------------------------------------------------" << G4endl << G4endl;

  std::map<std::string, G4Material*>::iterator i;
  for ( i=MaterialMap.begin(); i!=MaterialMap.end(); i++ ) {
    G4cout << "\t" << (i->second)->GetName() << "\t\t" 
	 << "\t" << (i->second)->GetDensity() << "\t\t" 
	 << "\t" << (i->second)->GetTemperature() << "\t\t" 
	 << "\t" << (i->second)->GetPressure() << "\t\t" 
	 << G4endl;
    
  }
}
