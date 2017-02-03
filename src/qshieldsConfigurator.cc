// --------------------------------------------------------------
//      GEANT 4 - qshields 
//
//      Simulation of CUORE0 experiment
//      Oliviero Cremonesi, Silvia Capelli, Luca Gironi  
//	Milano Bicocca 2009 
// --------------------------------------------------------------

#include <unistd.h>
#include <string.h>
#include <getopt.h>
#include <fstream>
#include <iostream>
#include "qshieldsConfigurator.hh"

#define qVERSION	"4.2"
#define ARGS	"A:a:b:BC:c:D:dE:e:F:f:G:g:Hhi:I:j:J:k:K:L:l:M:m:N:n:O:o:P:p:Qq:r:R:S:s:t:TU:u:VvW:w:X:x:y:YZ:z:012:3:4"

//extern int getopt(int argc, char * const argv[], const char *optstd::string);

extern char *optarg;

static qshieldsConfigurator* TheConfiguration = 0;

qshieldsConfigurator* qshieldsConfigurator::GetInstance()
{
  if ( !TheConfiguration ) TheConfiguration = new qshieldsConfigurator();
  return TheConfiguration;
} 

qshieldsConfigurator* qshieldsConfigurator::GetInstance(int argc, char** argv)
{
  if ( !TheConfiguration ) TheConfiguration = new qshieldsConfigurator(argc,argv);
  return TheConfiguration;
} 

qshieldsConfigurator::qshieldsConfigurator(int argc, char** argv)
{
  G4int  ib, c, errflg = 0;
  optarg = NULL;
  std::vector<G4double> TVtmp;

  char cmnd[256];
  cmnd[0] = '\0'; // initialize string to zero length
  assert(argc >= 1); // sanity check
  for( G4int i=0; i<argc-1; i++ ) {
    strcat(cmnd,argv[i]); strcat(cmnd," "); // append command line elements
  }
  strcat(cmnd,argv[argc-1]); // append last command line element without ending space

  char* s;

  v1 = 0;
  v2 = 0;

  setparameters();
  if ( !TheConfiguration ) TheConfiguration = new qshieldsConfigurator();

  while (!errflg && (c = getopt(argc, argv, ARGS))	  != -1)
    switch (c) {
    case 'a'		 :
      n_particles++;
      PTime.push_back( 0. );
      PType.push_back( 2 );
      PEnergy.push_back( atof (optarg)*MeV );
      PMass.push_back( 0. );
      PEexc.push_back( 0. );
      break;
    case 'A'		 :
      CAngle = atof(optarg);
      break;
    case 'b'		 :
      n_particles += 2;
      TVtmp = qshieldsConvert( std::string(optarg),2 );
      DoubleBeta =  G4int( TVtmp[0] );
      if(DoubleBeta==8) {
        std::vector<std::string> ccns;
        ccns = qshieldsSplit( std::string(optarg),3 );
        fileE = ccns[1];
        fileA = ccns[2];
        G4cout << "Double beta configurator: " << fileE << " " << fileA << G4endl;
      } else QBeta =  TVtmp[1]*MeV ;
      G4cout << "Double beta" << DoubleBeta << " Q: " << QBeta << G4endl;
      for( ib=0; ib<2; ib++ ) {
        PTime.push_back( 0. );
        PType.push_back( -1 );
        PEnergy.push_back( 0. );
      }
      if( SType == 0 ) SType = 14;	// Uniform generation in TeO2 crystals
      break;
    case 'B'		 :
      GenDecAN = true;
      break;
    case 'C'		 :
      VolumeNumber = atoi (optarg);
      G4cout << "****** Single Crystal Source number: " << VolumeNumber <<  "****** " << G4endl;
      //	   if( VolumeNumber<0 ) ModuleNumber = -VolumeNumber;		//SC 23.07.13
      break;
    case 'c'		 :
      Mode = 3;
      Debug = true;
      s = optarg;

      sscanf(s, "%d", &DebugEvents);
      if((s = strchr(s, ',')) == NULL)
        break;
      s++;
      sscanf(s, "%d", &v1);
      break;
    case 'd'		 :
      PrintGendec = true ;
      break;
    case 'D'		 :
      if( SMode>0 && SMode!=2 ) {
        G4cout << "Only identical generation options can be repeated - Disc generation neglected!" << G4endl;
      }
      else {
        SMode = 2;
        TVtmp = qshieldsConvert( std::string(optarg),2 );
        SWDirCos.push_back( G4ThreeVector(TVtmp[0],TVtmp[1]*cm,0) );
        G4cout << "Disc source n. " << SWDirCos.size() << ". Axis and Radius: " << SWDirCos.back().x() << " " << SWDirCos.back().y() << G4endl;
      }
      break;
    case 'E'		 :
      ParticleMask = ParticleMask & (0xFFFF ^ 1 << pId[ optarg[0] ]);
      G4cout << "ParticleMask: " << ParticleMask << " " << optarg[0] << " " << pId[ optarg[0] ] << G4endl;
      break;
    case 'e'		 :
      n_particles++;
      PTime.push_back( 0. );
      PType.push_back( -1 );
      PEnergy.push_back( atof (optarg)*MeV );
      PMass.push_back( 0. );
      PEexc.push_back( 0. );
      break;
    case 'f'		 :
      tThresh = atof (optarg);
      break;
    case 'g'		 :
      n_particles++;
      PTime.push_back( 0. );
      PType.push_back( 0 );
      PEnergy.push_back( atof (optarg)*MeV );
      PMass.push_back( 0. );
      PEexc.push_back( 0. );
      break;
    case 'G'		 :
      n_particles++;
      GendecFile = G4String (optarg);
      GenDec = true;
      PTime.push_back( 0. );
      PType.push_back( 0 );
      PEnergy.push_back( 0. );
      PMass.push_back( 0. );
      PZeta.push_back( 0. );
      PEexc.push_back( 0. );
      break;
    case 'h'	 :
      qshieldsHelp();
      exit (1);
    case 'H'		 :
      PrintInfo = true ;
      break;
    case 'i'		 :
      RSeed = atoi (optarg);
      break;
    case 'I'		 :
      //	   ggCorrelation = 1;
      TauLim = atof (optarg);
      G4cout << "Minimum decay time: " << TauLim << " sec" << G4endl;
      break;
    case 'j'		 :
      StepFactor =  atof(optarg);
      break;
    case 'J'		 :
      SpecOutFile = G4String (optarg);
      SpecFile.open ( SpecOutFile.c_str() );
      oSpec = 1;
      break;
    case 'k'		 :
      CuoreMap( std::string(optarg),true );
      break;
    case 'K'		 :
      SpecInpFile = G4String (optarg);
      iSpec = 1;
      break;
    case 'l'		 :
      PhysListName = G4String(optarg);
      break;
    case 'L'		 :
      LogNumber =  atoi(optarg);
      break;
    case 'm'		 :
      MacroFile = std::string(optarg);
      Mode = 1;
      break;
    case 'M'		 :
      NumberOfLoops = atoi (optarg);
      Mode = 2;
      break;
    case 'n'		 :
      n_particles++;
      PTime.push_back( 0. );
      PType.push_back( +3 );
      PEnergy.push_back( atof (optarg)*eV );
      PMass.push_back( 0. );
      PEexc.push_back( 0. );
      neutron_gen = true;
      break;
    case 'N'		 :
      NumberOfEvents = std::string (optarg);
      Mode = 2;
      break;
    case 'o'		 :
      G4int iType;
      if( optarg[0]=='t' ) iType = 1; // Text output
      else if( optarg[0]=='r' ) iType = 2; // RootTree output
      else {
        G4cout << "*** ERROR in output file specification *** " << G4endl;
        exit(1);
      }
      OutType |= iType;
      OutFile[iType-1] = std::string (&optarg[1]);
      break;
    case 'O'		 :
      if( ParticleMask == 0xFFFF ) ParticleMask = 0;
      ParticleMask = ParticleMask | 1 << pId[ optarg[0] ];
      G4cout << "ParticleMask: " << ParticleMask << " " << optarg[0] << " " << pId[ optarg[0] ] << G4endl;
      break;
    case 'p'		 :
      n_particles++;
      PTime.push_back( 0. );
      PType.push_back( 1 );
      PEnergy.push_back( atof (optarg)*MeV );
      PMass.push_back( 0. );
      PEexc.push_back( 0. );
      break;
    case 'P'		 :
      TVtmp = qshieldsConvert( std::string(optarg),3 );
      SPosition.push_back( G4ThreeVector( TVtmp[0], TVtmp[1], TVtmp[2] ) * cm);
      G4cout << "Source position n. " << SPosition.size() << ": " << SPosition[SPosition.size()-1] << G4endl;
      // UniverseR = sqrt( SPosition.x()*SPosition.x() + SPosition.y()*SPosition.y()+ SPosition.z()*SPosition.z() );
      break;
    case 'Q'		 :
      xyz_out=true;
      break;
    case 'q'		 :
      EDistrUniform =  atof(optarg)*MeV;
      break;
    case 'r'		 :
      n_particles++;
      TVtmp = qshieldsConvert( std::string(optarg),4 );
      PTime.push_back( 0. );
      PType.push_back( -2 );
      PEnergy.push_back( TVtmp[0]*MeV );
      PMass.push_back( TVtmp[1] );
      PZeta.push_back( TVtmp[2] );
      PEexc.push_back( TVtmp[3] );
      break;
    case 'R'		 :
      n_particles++;
      PTime.push_back( 0. );
      PType.push_back( 6 ); //negative muon 
      PEnergy.push_back( 0. );
      PMass.push_back( 0. );
      PEexc.push_back( 0. );
      TVtmp = qshieldsConvert( std::string(optarg),3 );
      MuPar = G4ThreeVector( TVtmp[0], TVtmp[1], TVtmp[2] );
      G4cout << "Muon generation parameters: " << MuPar << G4endl;
      break;
    case 's'		 :
      SAF = 1;
      TVtmp = qshieldsConvert( std::string(optarg),3 );
      SAngle = G4ThreeVector( TVtmp[0],TVtmp[1],TVtmp[2] );
      G4cout << "Source Solid Angle: " << SAngle << G4endl;
      break;
    case 'S'		 :
      if( SMode>0 && SMode!=3 ) {
        G4cout << "Only identical generation options can be repeated - Sphere generation neglected!" << G4endl; 
      }
      else {
        SMode = 3;
        SWDirCos.push_back( G4ThreeVector(atof( optarg )*cm,0,0) );
        SphereRadius = SWDirCos.back().x()*cm;
        G4cout << "Sphere source n. " << SWDirCos.size() << ". Radius: " << SWDirCos.back().x() << G4endl;
      }
      break;
    case 't'		 :
      if( SMode>0 && SMode!=4 ) {
        G4cout << "Only identical generation options can be repeated - Disc generation neglected!" << G4endl;
      }
      else {
        SMode = 4;
        TVtmp = qshieldsConvert( std::string(optarg),3 );
        SWDirCos.push_back( G4ThreeVector(TVtmp[0],TVtmp[1]*cm,TVtmp[2]*cm) );
        G4cout << "Tube source n. " << SWDirCos.size() << ". Axis, Radius and heigth: " << SWDirCos.back().x() << " " << SWDirCos.back().y() << " " << SWDirCos.back().z()*2 << " " << G4endl;
      }
      break;
    case 'T'		 :
      EThermalAvr = 0.025*eV;
      break;
    case 'u'		 :
      CuoreMap( std::string(optarg),false );
      break;
    case 'U'		 :
      unsigned int il;
      for( il=0; il<strlen(optarg); il++ ) if( optarg[il] == ':' ) break;
      if( il<strlen(optarg) ) {
        SElem = atoi(optarg+il+1);
        optarg[il] = 0;
      }
      SType = atoi(optarg);
      break;
    case 'v'		 :
      G4cout << "To use vrml interactive graphic analysis you must: " << G4endl
        << " export G4VRMLFILE_VIEWER=vrmlview" << G4endl
        << "... then, to simply view the detector:" << G4endl
        << " /vis/scene/create" << G4endl
        << " /vis/sceneHandler/create VRML2FILE" << G4endl
        << " /vis/viewer/create" << G4endl
        << " /vis/enable true" << G4endl
        << " /vis/viewer/refresh" << G4endl
        << " /vis/viewer/update" << G4endl << G4endl
        << "... or (DAWN ps interface):" << G4endl
        << " /vis/scene/create" << G4endl
        << " /vis/sceneHandler/create DAWNFILE" << G4endl
        << " /vis/viewer/create" << G4endl
        << " /vis/enable true" << G4endl
        << " /vis/viewer/refresh" << G4endl
        << " /vis/viewer/update" << G4endl << G4endl
        << "to analyze instead a track:" << G4endl
        << " /vis/scene/create" << G4endl
        << " /vis/sceneHandler/create VRML2FILE" << G4endl
        << " /vis/viewer/create" << G4endl
        << " /vis/enable true" << G4endl
        << " /vis/viewer/refresh" << G4endl
        << " /tracking/storeTrajectory 1 " << G4endl
        << " /vis/scene/add/trajectories " << G4endl
        << " /run/beamOn 1 " << G4endl
        << "to test geometry:" << G4endl
        << " /gun/particle chargedgeantino" << G4endl 
        << " /tracking/verbose 1" << G4endl
        << " /run/beamOn 1" << G4endl << G4endl
        << "to draw multiple tracks:" << G4endl
        << " /vis/open VRML2FILE" << G4endl
        << " /vis/drawVolume" << G4endl
        << " /vis/scene/add/trajectories" << G4endl
        << " /vis/scene/endOfEventAction accumulate" << G4endl
        << " /tracking/storeTrajectory 1" << G4endl
        << " /run/beamOn 1" << G4endl;
      exit(1);
      break;
    case 'V'		 :
      G4cout << "CUORE elements map: " << G4endl
        << 1 << ":\t PET shield" << G4endl
        << 2 << ":\t Top PET shield" << G4endl
        << 3 << ":\t Boric Acid (H3B04) external shield" << G4endl
        << 4 << ":\t External PB shields" << G4endl
        << 5 << ":\t First Cu radiation shield (300 K)" << G4endl
        << 6 << ":\t First radiation shield flange (Fe)" << G4endl
        << 7 << ":\t Second Cu radiation shield (40 K)" << G4endl
        << 8 << ":\t Second radiation shield flange (Cu)" << G4endl
        << 9 << ":\t Third Cu radiation shield (4 K)" << G4endl
        << 10 << ":\t Third radiation shield flange (Cu)" << G4endl
        << 11 << ":\t Fourth Cu radiation shield (600 mK)" << G4endl
        << 12 << ":\t Fourth radiation shield flange (Cu) " << G4endl
        << 43 << ":\t CuBe Cilindrical Connectors below R4-Flange" << G4endl
        << 13 << ":\t Internal Lead shield (in 3rd shield)" << G4endl
        << 14 << ":\t Internal Lead shield copper parts (Up Ring and Bottom Crandle)" << G4endl
        << 15 << ":\t Steel Rods of the Internal Lead shield" << G4endl
        << 16 << ":\t Fifth Cu radiation shield (50 mK)" << G4endl
        << 17 << ":\t Fifth radiation shield flange (Cu)" << G4endl
        << 18 << ":\t Sixth Cu radiation shield (10 mK) " << G4endl
        << 19 << ":\t Sixth radiation shield (10mK) flange (Cu)" << G4endl
        //	   << 45 << ":\t Sixth radiation shield (10mK) Poly wrapping -- On internal 10mK side surface" << G4endl
        << 20 << ":\t Cu Top Plate Above Std Pb (h=1.2 cm)" << G4endl
        << 21 << ":\t Std Pb plug (h=25 cm)" << G4endl
        << 23 << ":\t DCS Tubes in Std Pb plug" << G4endl
        << 22 << ":\t Cu Tiles on lower 10 mK shield" << G4endl
        << 24 << ":\t Cu Bottom Plate Below Roman Pb (h=4.6 cm)" << G4endl
        << 25 << ":\t Cu TSP plate (h=5 cm)" << G4endl
        << 26 << ":\t Detector Cu Up plate (h=1 mm)" << G4endl
        //	   << 27 << ":\t Detector Cu Up plate Poly wrapping (h=70 um)" << G4endl
        << 28 << ":\t Detector Cu Bottom plate (h=1 mm) " << G4endl
        //	   << 29 << ":\t Detector Cu Bottom plate Poly wrapping (h=70 um)" << G4endl
        << 30 << ":\t Internal Calibration copper tubes " << G4endl
        << 31 << ":\t Internal Calibration Source String (Teflon Layer)" << G4endl
        << 32 << ":\t Internal Calibration Source String (internal Cu Layer)" << G4endl
        //	   << 33 << ":\t Internal Calibration Source String (Tungsten)" << G4endl
        << 34 << ":\t External Calibration Source String (Teflon Layer)" << G4endl
        << 35 << ":\t External Calibration Source String (internal Cu Layer)" << G4endl
        //	   << 36 << ":\t External Calibration Source String (Tungsten)" << G4endl
        << 37 << ":\t CUORE Detector " << G4endl
        << 38 << ":\t Frames" << G4endl
        << 39 << ":\t Wire Trays (Cu) " << G4endl
        << 44 << ":\t PEN Tape (Poly) (placed inside the Wire Trays)" << G4endl
        << 40 << ":\t Wire Pads (Cu)" << G4endl
        << 41 << ":\t NTD Thermistor" << G4endl
        << 42 << ":\t PTFE crystal supports" << G4endl
        //	   << 43 << ":\t Signal Gold Wires (STILL IN PROGRESS, DO NOT USE THEM!!)" << G4endl;
        << 46 << ":\t SuperInsulation Layers (Mylar) on 40K and 4 K shields" << G4endl
        << 54 << ":\t Internal Calibration Weight Capsules (8), Activity = 0.538 Bq" << G4endl
        << 55 << ":\t Internal Calibration Source Capsules (21), Activity = 1.961 Bq" << G4endl
        << 56 << ":\t Internal Calibration Source Capsules (4), Activity = 0.309 Bq" << G4endl
        << 57 << ":\t Internal Calibration Source Capsules (1), Activity = 0.796 Bq" << G4endl
        << 58 << ":\t External Calibration Weight Capsules (8), Activity = 5.115 Bq" << G4endl
        << 59 << ":\t External Calibration Source Capsules (20), Activity = 9.524 Bq" << G4endl
        << 60 << ":\t External Calibration Source Capsules (5), Activity = 4.745 Bq" << G4endl;

      exit(1);
      break;
    case 'w'		 :
      TVtmp = qshieldsConvert( std::string(optarg),3 );
      WiR =  TVtmp[0]*mm*0.001;
      WiL =  TVtmp[1]*mm;
      WiZ =  TVtmp[2]*mm;
      break;
    case 'W'		 :
      if( SMode>0 && SMode!=1 ) {
        G4cout << "Only identical generation options can be repeated - Wire generation neglected!" << G4endl;
      }
      else {
        SMode = 1;
        TVtmp = qshieldsConvert( std::string(optarg),3 );
        SWDirCos.push_back( G4ThreeVector( TVtmp[0], TVtmp[1], TVtmp[2]*cm) );
        G4cout << "Sphere source n. " << SWDirCos.size() << ". Length: " << SWDirCos.back().z()*2 << ". Orientation: " << 
          G4ThreeVector( TVtmp[0], TVtmp[1], sqrt(1.-TVtmp[0]*TVtmp[0]-TVtmp[1]*TVtmp[1]) ) << G4endl;
      }
      break;
    case 'x'		 :
      ExtShieldMask = atoi (optarg);
      break;
    case 'X'		 :
      SurfExpGen = LocatePlusSign( optarg );
      if( SurfExpGen ) G4cout << "Exponential surface profile\n"<<G4endl;
      TVtmp = qshieldsConvert( std::string(optarg),3 );
      SType = static_cast<int> (-fabs(TVtmp[0]));
      SOffs = TVtmp[1]*cm;
      SDepth = TVtmp[2]*cm;
      break;
    case 'y'		 :
      n_particles++;
      PTime.push_back( 0. );
      PType.push_back( 5 );
      PEnergy.push_back( atof (optarg)*MeV );
      PMass.push_back( 0. );
      PEexc.push_back( 0. );
      break;
    case 'Y'		 :
      outMicro = true;
      break;
    case 'z'		 :       //S.C. per mettere taglio in energia nella propagazione
      minEn = atof(optarg)*MeV;;
      break;
    case 'Z'		 :
      n_particles++;
      PTime.push_back( 0. );
      PType.push_back( +6 );
      PEnergy.push_back( atof (optarg)*MeV );
      PMass.push_back( 0. );
      PEexc.push_back( 0. );
      break;
    case '0'		 :
      qshieldsReleaseNotes();
      exit (1);
      break;
    case '1'		 :
      SingleDecayGranularity = true;
      break;
    case '2'		 :
      qshieldsTest = true;
      NumqshieldsTest =  atoi(optarg);
      NumberOfEvents = std::string ("1");
      Mode = 2;
      break;
    case '3'		 :
      BaseCut = atof (optarg)*mm;
      break;
    case '4'		 :
      singleTower = 5;
      NumberOfPlanes = 2;
      //	   TNumberOfDetectors = 988;
      break;
    default :
      errflg++;
      qshieldsHelp();
      exit (1);
    }

  assert( n_particles || GenDec || GenDecAN || Debug);

  if(GraphicMask[0]+GraphicMask[1]==0) {
    GraphicMask[0] = 0xFFFF;
    GraphicMask[1] = 0xFFFF;
  }

  G4cout << "Microphysics output: " << outMicro << G4endl; 

  if ( Mode == 2 ) {
    if( OutType & 1 ) {
      DataFile.open ( OutFile[0].c_str() );
      for( G4int i=0; i<argc; i++ ) DataFile << argv[i] << " ";
      DataFile << G4endl;
    }
    else if( OutType & 2 ) {
      ROOTFile = new TFile( OutFile[1].c_str(),"RECREATE");
      command = std::string(cmnd);
      INFOtree = new TTree("INFOtree", "");
      INFOtree->Branch("CommandLine",&command);
      INFOtree->Branch("Ntot",&iNtot,"Ntot/l");
      INFOtree->Branch("Nchain",&dChainNumber,"Nchain/l");
      //    INFOtree->Branch("SourceType",&SourceType,"SourceType/I");	
      // ...
      INFOtree->Fill();
      INFOtree->AutoSave("SaveSelf");
      // Create ROOT file

      ROOTTree = CreateROOTTree();
    }
  }

  if( SWDirCos.size() >0 ) { 
    if( SWDirCos.size() > SPosition.size() ) {
      G4cout << "*** ERROR ***: Positions must be specified for all geometric items (Wires, Discs, ...)" << G4endl;
      exit(1);
    }
    else if( SWDirCos.size() < SPosition.size() ) {
      G4cout << "*** WARNING ***: More positions than needed have been specified. Last will be neglected" << G4endl;
      G4cout << "*** Press return to continue ***" << G4endl;
      G4String mystr;
      getline (G4cin, mystr);
    }
  }

  SGNorm = 0.;
  G4int csize = SWDirCos.size();
  if( SMode == 1) for(G4int ii=0; ii<csize; ii++ ) SGNorm += SWDirCos[ii].z();
  else if( SMode == 2) for(G4int ii=0; ii<csize; ii++ ) SGNorm += SWDirCos[ii].y()*SWDirCos[ii].y();
  else if( SMode == 3) for(G4int ii=0; ii<csize; ii++ ) SGNorm += SWDirCos[ii].x()*SWDirCos[ii].x()*SWDirCos[ii].x();
  else if( SMode == 4) for(G4int ii=0; ii<csize; ii++ ) SGNorm += SWDirCos[ii].y()*SWDirCos[ii].z();

  if( GenDecAN ) GenDec=false;  

  if( EThermalAvr != 0.0 ) { 
    // EThermalAvr = kT (cio' che e' passato con l'opzione n)
    // ThermalMax = max della distribuzione maxwelliana dei termici
    // EThermalMax = E_max sfino a cui sono distribuiti i neutroni
    ThermalMax = EThermalAvr/1000.;
    G4double wtmp = ThermalMax;
    G4double tmp;
    G4double stmp = 0.0;
    while( (tmp=sqrt(ThermalMax)*exp(-ThermalMax/EThermalAvr )) > stmp ) {
      stmp = tmp;
      ThermalMax += wtmp;
    }
    ThermalMax = stmp;

    wtmp = EThermalAvr;
    EThermalMax = EThermalAvr;
    while( sqrt(EThermalMax)*exp(-EThermalMax/EThermalAvr ) > 1.e-4*ThermalMax ) {
      EThermalMax += wtmp;
    }
    if( PrintInfo ) G4cout << "Thermal neutron parameters: " << EThermalAvr << " "
      << EThermalMax << " " << ThermalMax << G4endl;
  }
} 

qshieldsConfigurator::qshieldsConfigurator() {
  setparameters();
}

qshieldsConfigurator::~qshieldsConfigurator() {
  PEnergy.clear();
  PTime.clear();
  PType.clear();
  pId.clear();
}

void qshieldsConfigurator::SetType( char typ, double energ ) {
  partyp = typ;
  parten = energ;
}

void qshieldsConfigurator::SetParticle( double pe, int pq, int pm, int pz ) { 
  n_particles = 1;
  PEnergy[0] = pe*MeV;
  PType[0] = pq; 
  PMass[0] = static_cast<double> ( pm );
  PZeta[0] = static_cast<double> ( pz );
}
void qshieldsConfigurator::SetParticle( double pe, int pq, int pm, int pz, double px ) { 
  n_particles = 1;
  PEnergy[0] = pe*MeV;
  PType[0] = pq; 
  PMass[0] = static_cast<double> ( pm );
  PZeta[0] = static_cast<double> ( pz );
  PEexc[0] = px;
}

void qshieldsConfigurator::setparameters() 
{
  EThermalAvr = 0.0;
  EDistrUniform = 0.0;
  neutron_gen = false;
  DoubleBeta = 0;
  LogNumber = 1000;
  SAF = 0;
  GenDec = false;
  GenDecAN = false;
  Debug = false;
  n_particles = 0;
  GraphicMask[0] = GraphicMask[1] = 0;
  ParticleMask = 0xFFFFFFFF;
  //  OutFile = "qshields.dat";
  OutType = 0;  
  CfgFile = "qshields.cfg";  
  GendecFile = G4String("");
  SpecInpFile = G4String("");
  SpecOutFile = G4String("");
  iSpec=0;
  oSpec=0;
  Mode = 0;
  NumberOfLoops = 1;
  TNumberOfDetectors = 988;
  NumberOfPlanes = 13;
  NPs = 0; //numero piani dei cristalli small
  NPb = 0; //numero piani dei cristalli big che stanno sotto quelli small
  NumberOfEvents = "0";
#ifdef G4GT10
  PhysListName = G4String("Shielding_LIV");
#else
  PhysListName = G4String("livermore");
#endif
  UniverseR = 0.;
  SType = 0;
  SMode = 0;
  SOffs = 0.;
  SDepth = 0.;
  RSeed = 0L;
  LastChainTime = 0.;
  pId['g'] = 0;
  pId['x'] = 0;
  pId['e'] = 1;
  pId['b'] = 1;
  pId['i'] = 1;
  pId['p'] = 2;
  pId['a'] = 3;
  pId['n'] = 4;
  pId['r'] = 5;
  pId['m'] = 6;
  dChainNumber = 0;
  iNtot = 0;
  NumCycles = 0;
  ExtShieldMask = 0xFFFFFFFF;
  CAngle = 360.;
  PrintInfo = false;
  PrintGendec = false;
  WiR = WiL = WiZ = 0.0;
  BaseCut = 1.0*mm;
  SurfExpGen = 0;
  VolumeNumber = -1;
  ModuleNumber = -1;
  ggCorrelation = 0;
  StepFactor=1;
  MuPar = G4ThreeVector();
  qshieldsTest = false;
  NumqshieldsTest = 0;
  tThresh = 3.16e7; // 1y, maximum correlation time in radioactive chains. when t>tThresh the decay position is changed
  SingleDecayGranularity = false;
  TauLim = 1.e-6;
  AtomicNumber=0;
  MassNumber=0;
  StableNucleus=true;
  AStop=0;
  ZStop=0;
  ChainE = 0.;
  ChainA = 0.;
  ChainZ = 0.;
  ChainX = 0.;
  xyz_out=false;
  outMicro = false;
  SphereRadius = 0.0;
  minEn = 0.0;
  MaxStack=0;
  singleTower=0;

}

void qshieldsConfigurator::RemoveWhiteSpaces ( std::string& s ) {
  for ( std::string::iterator i=s.begin(); i!=s.end(); i++ ) {
    if ( *i == ' ' ) s.erase ( i );
  }
}

G4int qshieldsConfigurator::LocatePlusSign ( std::string s ) {
  G4int iok=0;
  for ( std::string::iterator i=s.begin(); i!=s.end(); i++ ) {
    if ( *i == '+' ) iok=1; 
  }
  return iok;
}

void qshieldsConfigurator::qshieldsReleaseNotes() 
{
  G4cout << G4endl << "\033[1;34mqshields: Simulation of CUORE experiment. Version " << qVERSION << G4endl 
    << qVERSION << " => Code updated to all changes introduced in MCuoreZ in the past months" << G4endl << G4endl
    << "\033[0m4.1 - After a number of problems with g4 generator we have re-introduced gendec decay generator (-B option)." << G4endl
    << "4.0 Assembly geometry improved" << G4endl
    << "- From now on, the geometry is described in *.icc files" << G4endl
    << "3.0 New class for distributed sources" << G4endl
    << "- From now on, no problems for distributed source generation on any set-up element" << G4endl
    << "2.0.2 Metastable states" << G4endl
    << "- Metastable states are decoupled from main radioactive decay (options -I)" << G4endl
    << "2.0.1 Tracking" << G4endl
    << "- Re-introduced -O option for particle propagation selection" << G4endl
    << "2.0.0 Minor revision" << G4endl
    << "- Solved problem with stable ions propagation" << G4endl
    << "- Versioning and release notes" << G4endl
    << "- Re-introduced pHysics Lists selection (-l option)" << G4endl << G4endl;
}

void qshieldsConfigurator::qshieldsHelp() 
{
  G4cout 
    << "\033[1;32m" 
    << "  qshields: Simulation of the CUORE experiment. Version" << qVERSION << G4endl
    << "  Usage: qshields [Options] " << G4endl << G4endl
    //<< "\033[1;30m" // deprecated by michsakai
    << "\033[1;39m" // use terminal default color
    << "Options: " << G4endl
    //<< "\033[0;30m" // deprecated by michsakai
    << "\033[0;39m"  // use terminal default color
    << "|--------|---------------|------------------------------|---------------------------------------------|" << G4endl
    << "| Option | Argument      | Description                  | Additional information                      |" << G4endl
    << "|--------|---------------|------------------------------|---------------------------------------------|" << G4endl
    << "| m      | macrofile     | Execute a macro file         |                                             |" << G4endl
    << "| o      | cFile         | Output file                  | The leading character 'c' should be either: |" << G4endl
    << "|        |               |                              |       - t -> Text file                      |" << G4endl 
    << "|        |               |                              |       - r -> ROOT file                      |" << G4endl 
    << "| G      | chain-file    | Decay chain generation       | Full path to file must be specified         |" << G4endl
    << "|        |               | input file                   | File format                                 |" << G4endl
    << "|        |               |                              | A Z E_exc(MeV) BR(%) id n list              |" << G4endl
    << "|        |               |                              | id = isotope identifier                     |" << G4endl
    << "|        |               |                              | n = number of successors                    |" << G4endl
    << "|        |               |                              | list = list of successor id                 |" << G4endl
    << "| B      |               | gendec radioactive decays    | Activates the gendec decay manager          |" << G4endl
    << "| f      | time          | decay correlation time       | maximum chain correlation time (default 1y) |" << G4endl
    << "|        |               |                              | when event time distance dt > time the decay|" << G4endl
    << "|        |               |                              | position is changed                         |" << G4endl
    << "| 4      |               | mini-tower                   |                                             |" << G4endl
    << "| M      | number        | number of iterations         | generates number iterations                 |" << G4endl
    << "| N      | number        | nomber of events             | generates  number events per iteration      |" << G4endl
    << "| Q      |               | add source points to output  | only for root output files                  |" << G4endl
    << "\033[0;31m" 
    << "|--- Particle Options (* can be repeated for cascades): ----------------------------------------------|" << G4endl
    << "| e      | energy*       | generates an electron        | energy (MeV)                                |" << G4endl
    << "| p      | energy*       | generates a positron         | energy (MeV)                                |" << G4endl
    << "| Z      | energy*       | generates a muon             | energy (MeV)                                |" << G4endl
    << "| g      | energy*       | generates a gamma            | energy (MeV)                                |" << G4endl
    << "| a      | energy*       | generates an alpha           | energy (MeV)                                |" << G4endl
    << "| n      | energy*       | generates a neutron          | energy (MeV)                                |" << G4endl
    << "| b      | (t,Q/names)   | generate a Double Beta Decay | t is an index specifying the decay mode     |" << G4endl
    << "|        |               |                              |   1 = 0 neutrinos, 0+, 2n, M_nu             |" << G4endl
    << "|        |               |                              |   2 = 2 neutrinos, 0+, 2n                   |" << G4endl
    << "|        |               |                              |   3 = 0 neutrinos, 1 Majoron, 0+, 2n        |" << G4endl
    << "|        |               |                              |   4 = 0 neutrinos, 2 Majorons, 0+, 2n       |" << G4endl
    << "|        |               |                              |   5 = 0 neutrinos, 0+, RH                   |" << G4endl
    << "|        |               |                              |   6 = 2 neutrinos, 2+, 2n                   |" << G4endl
    << "|        |               |                              |   7 = 0 neutrinos, 2+, 2n, RH               |" << G4endl
    << "|        |               |                              |   8 = energy and angle distribution files   |" << G4endl
    << "|        |               |                              | Q = transition energy (t<7)                 |" << G4endl
    << "|        |               |                              | names = list of energy and angle            |" << G4endl
    << "|        |               |                              |         distribution files                  |" << G4endl
    << "| r      | (E,A,Z,Exc)*  | generates an ion             | E = kinetic energy (MeV)                    |" << G4endl
    << "|        |               |                              | Exc = nuclear excitation energy (MeV)       |" << G4endl
    << "| R      | (R,E-,E+)     | generates a muon on emisphere| R = radius of (top) emisphere (cm)          |" << G4endl
    << "|        |               |                              | E-,E+ energy interval (GeV)                 |" << G4endl
    << "|        |               |                              | The default center is the bottom of the     |" << G4endl
    << "|        |               |                              | concrete pedestal. It can be modified       |" << G4endl
    << "|        |               |                              | providing an offset R through the P option  |" << G4endl
    << "|-----------------------------------------------------------------------------------------------------|" << G4endl
    << "| O      | id            | excludes id particles        | id is the letter used to specify a particle |" << G4endl
    << "|        |               |                              | generation (eg e for an electron, see above)|" << G4endl
    << "|        |               |                              | Usually used for decay chains               |" << G4endl
    << "| E      | id            | generates only id particles  | same as for O option                        |" << G4endl
    << "|        |               |                              |                                             |" << G4endl
    << "| T      |               | generates a thermal neutron  |                                             |" << G4endl
    << "| q      | E             | uniform energy interval      | generates uniformly in 0-E (MeV)            |" << G4endl
    << "| K      | file          | energy spectrum              | input file format:                          |" << G4endl 
    << "|        |               |                              | E(MeV),Flux (ascending order)               |" << G4endl
    << "| I      | time          | minimum decay time (sec)     | if time difference in larger than time, then|" << G4endl
    << "|        |               |                              | a radioactive chain is \"paused\" to intend   |" << G4endl
    << "|        |               |                              | that all data are saved and a new isotope   |" << G4endl
    << "|        |               |                              | is generated                                |" << G4endl
    << "| z      | E             | propagates only particles    | E = energy threshold                        |" << G4endl
    << "|        |               | with enegy larger than E     |                                             |" << G4endl
    << "\033[0;34m" 
    << "|--- Source Options (* can be repeated for cascades, + use P option for central position): ------_----|" << G4endl
    << "| P      | (x,y,z)       | pointlike source             | in position (x,y,z) [default (0,0,0)]       |" << G4endl
    << "| W+     | (CX,CY,HL)    | linear source                | HL = line half-length                       |" << G4endl
    << "|        |               |                              | CX/CY = X/Y direcotr cosines                |" << G4endl
    << "| w      | (r,l,z)       | wire source (Z axis)         | r = radius (microns)                        |" << G4endl 
    << "|        |               |                              | l = length (mm)                             |" << G4endl
    << "|        |               |                              | z = Z offset (m)                            |" << G4endl
    << "| D+     | (N,R)         | disc source                  | N = normal direction (1/2/3)                |" << G4endl 
    << "|        |               |                              | R = disc radius (cm)                        |" << G4endl
    << "| t+     | (A,R,HH)      | tube source                  | A = tube axis (1/2/3)                       |" << G4endl 
    << "|        |               |                              | R = tube radius (cm)                        |" << G4endl
    << "|        |               |                              | HH = half-height (cm)                       |" << G4endl
    << "| S+     | R             | surface of a sphere          | R = sphere radius (cm)                      |" << G4endl
    << "| U      | n             | uniform in element volume    | n = CUORE0 element id (use option V)        |" << G4endl
    << "| X      | (n,o,d)       | unifor on element surface    | n = CUORE0 element id (use option V)        |" << G4endl
    << "|        |               |                              | o = (additive) depth offset (cm)            |" << G4endl
    << "|        |               |                              | d = depth (cm)                              |" << G4endl
    << "|        |               |                              |   d<0 => uniform in layer of thickness d    |" << G4endl
    << "|        |               |                              |   +d (explicit \'+\' sign) => exp(-x/d)       |" << G4endl 
    << "| s      | (k,min,max)   | solid angle limits           | k =  direction specifier                    |" << G4endl
    << "|        |               |                              |      1/2/3: x/y/z                           |" << G4endl
    << "|        |               |                              |      12/13/...: xy/xz radial direction      |" << G4endl
    << "|        |               |                              |      123: radial direction                  |" << G4endl
    << "|        |               |                              | min/max: coside director limits             |" << G4endl
    << "\033[0;33m" 
    << "|--- Graphic Options: --------------------------------------------------------------------------------|" << G4endl
    << "| k      | (n1,n2,...)   | draw only listed elements    | n1,n2,... = list of CUORE0 elements (V opt) |" << G4endl
    << "| u      | (n1,n2,...)   | does not draw listed elements| n1,n2,... = list of CUORE0 elements (V opt) |" << G4endl
    << "| A      | angle         | draws angular section        | angle = opening angle                       |" << G4endl 
    << "\033[0;35m" 
    << "|--- Generic Options: --------------------------------------------------------------------------------|" << G4endl
    << "| 3      | cut           | changes cut value            | cut = cut value (mm) [default = 1 mm]       |" << G4endl 
    << "| C      | n             | specifies sub-element        | applies to U and X options (multi-elements) |" << G4endl 
    << "|        |               |                              |    n = source sub-element number            |" << G4endl
    << "| l      | name          | E.M. physics list selection  | name = livermore/standard/penelope          |" << G4endl 
    << "| L      | n             | changes notification freq.   | n = number of event for refreshing output   |" << G4endl 
    << "| i      | n             | initializes random number gen| n = random seed                             |" << G4endl 
    << "| j      | n             | scales propagation cuts      | n = multiplicative factor for base cut      |" << G4endl 
    << "| J      | name          | outputs generated particle   |                                             |" << G4endl 
    << "|        |               | energy spectrum              |                                             |" << G4endl 
    << "| v      |               | prints visualization info    |                                             |" << G4endl 
    << "| V      |               | prints summary of            |                                             |" << G4endl 
    << "|        |               | volume/surface informations  |                                             |" << G4endl 
    << "| H      |               | activates verbose mode       |                                             |" << G4endl 
    << "| Y      |               | activates microscopic info   | output is added to root output files        |" << G4endl 
    << "| 1      |               | single isotope number output | only for decay chains                       |" << G4endl
    << "| x      | n             | external shields mask        | n = bit mask (1:Lead,2:PolyB,3:Concrete)    |" << G4endl 
    << "| c      | N[,v1,v2]     | checks for volume            | outputs root file with generated points     |" << G4endl
    << "|        |               | intersections                | N is required, volumes are optional         |" << G4endl 
    << "|        |               |                              | N = number of particles                     |" << G4endl 
    << "|        |               |                              | v1,v2 = id\'s of volumes to compare          |" << G4endl 
    << "| 0      |               | release notes                |                                             |" << G4endl
    << "| h      |               | prints this help page        |                                             |" << G4endl
    << "|--------|---------------|------------------------------|---------------------------------------------|" << G4endl
    //<< "\033[0;30m" << G4endl // deprecated by michsakai
    << "\033[0m" << G4endl // revert terminal color and text style to default
    ;
}

std::vector<std::string> qshieldsConfigurator::qshieldsSplit( std::string pss, G4int n_el ) 
{
  G4int pst,pc=0;
  std::vector<std::string> pv;

  pss = std::string( pss.substr( 1,pss.length()-2 ).c_str() );

  pst = 0;
  for( unsigned int pen=0; pen<pss.length(); pen++ ) {
    if( pss.substr( pen,1 ) == std::string(",") && pc < n_el-1 ) {
      pv.push_back(  pss.substr( pst, pen-pst ) );
      pst = pen+1;
    }
  }
  pv.push_back( pss.substr( pst, pss.length()-pst ) );

  return( pv );
}
std::vector<G4double> qshieldsConfigurator::qshieldsConvert( std::string pss, G4int n_el ) 
{
  G4int pst,pc=0;
  std::vector<G4double> pv;

  pss = std::string( pss.substr( 1,pss.length()-2 ).c_str() );

  pst = 0;
  for( unsigned int pen=0; pen<pss.length(); pen++ ) {
    if( pss.substr( pen,1 ) == std::string(",") && pc < n_el-1 ) {
      pv.push_back( atof( pss.substr( pst, pen-pst ).c_str() ) );
      pst = pen+1;
    }
  }
  pv.push_back( atof( pss.substr( pst, pss.length()-pst ).c_str() ) );

  return( pv );
}


int qshieldsConfigurator::GetType( int type ) {

  if( type == 2 ) return 7;
  else if( type == -1 ) return 1;
  else if( type == 0 ) return 0;
  else if( type == 3 ) return 5;
  else if( type == 1 ) return 3;
  else if( type == 6 ) return 9;
  else if( type == -6 ) return 10;
  else if( type == -2 ) return 8;
  else return 11;
}

void qshieldsConfigurator::PrintStats() {
  G4double dTime;
  G4int numl;
  G4double wgt;
  G4ThreeVector Bari,Square;

  if( PTime[0] == LastChainTime ) dTime = 0.;
  else dTime = PTime[0];
  dDTime = dTime;

  if ( OutType ==0 ) { 
    if ( dCuEnergyDeposit != 0 ) 
      G4cout << dChainNumber << "\t" << dTime<< "\tCu\t" << dCuEnergyDeposit << G4endl;
    if ( dPTFEEnergyDeposit != 0 ) 
      G4cout << dChainNumber << "\t" << dTime<< "\tPTFE\t" << dPTFEEnergyDeposit << G4endl;
    if ( dNTDEnergyDeposit != 0 ) 
      G4cout << dChainNumber << "\t" << dTime<< "\tNTD\t" << dNTDEnergyDeposit << G4endl;
    if ( dOtherEnergyDeposit != 0 ) 
      G4cout << dChainNumber << "\t" << dTime<< "\tOther\t" << dOtherEnergyDeposit << G4endl;
  }

  G4int metamax;
  if( MetaTime>0 ) metamax=2;
  else metamax=1;
  for(G4int meta=0;meta<metamax;meta++) {
    if( meta>0 ) {
      dDTime = dTime = MetaTime;
    }
    for( G4int i=0; i<TNumberOfDetectors; i++ ) {
      if( OutType ==0 || OutType & 1 ) {
        for( G4int k=0; k<TNumberOfParticles; k++ ) {
          if ( EventEnergyDeposit[i][k][meta] != 0 ) {
            if ( OutType == 0 ) 
              G4cout << dChainNumber << "\t" << dTime<< "\t" << i << "\t" << EventEnergyDeposit[i][k][meta] << "\t" << k << G4endl;
            else
              DataFile << dChainNumber << "\t" << dTime<< "\t" << i << "\t" << EventEnergyDeposit[i][k][meta] << "\t" << k << G4endl;
            LastChainTime = PTime[0] = 0.;
            dTime = 0.;
          }
        }
      }

      if( OutType & 2 ) {
        // Microsphysics section: only for root output
        if( outMicro ) {
          wgt = 0.;
          Bari = G4ThreeVector();
          Square = G4ThreeVector();
          for( G4int l=0; l<numInt[i][meta]; l++ ) {
            Bari += IntDep[i][l][meta]*IntPos[i][l][meta];
            wgt += IntDep[i][l][meta];
            Square += IntDep[i][l][meta]*G4ThreeVector(IntPos[i][l][meta].x()*IntPos[i][l][meta].x(),
                                                       IntPos[i][l][meta].y()*IntPos[i][l][meta].y(),IntPos[i][l][meta].z()*IntPos[i][l][meta].z() );
          }			
          Bari /= wgt;	//-- baricenter of energy depositions
          Bari = Square/wgt-G4ThreeVector(Bari.x()*Bari.x(),Bari.y()*Bari.y(),Bari.z()*Bari.z() );
          dAvDist = sqrt(Bari.mag());
          if( FirstInteraction[i][meta] == 0 ) dFirstInteraction = std::string("Phot");
          else if( FirstInteraction[i][meta] == 1 ) dFirstInteraction = std::string("Comp");
          else if( FirstInteraction[i][meta] == 2 ) dFirstInteraction = std::string("Pair");
          else dFirstInteraction = std::string("Unkn");

          if( nClus[i][meta]>0 ) {
            numl = 0;
            dNClus = nClus[i][meta];
            dTLen = 0.;
            for( G4int j=0;j<nClus[i][meta];j++ ) { 
              if( nLen[i][j][meta]>0 ) {
                dTLen +=  cLen[i][j][0][meta];
                numl++;
              }
            }
          } 			
        }
        //--- End microphysics

        for( G4int k=0; k<TNumberOfParticles; k++ ) {
          dPEnergy[k] = EventEnergyDeposit[i][k][meta];
          dEnergy += dPEnergy[k];
        }
        dCuEnergyDeposit = CuEnergyDeposit[meta];
        dPTFEEnergyDeposit = PTFEEnergyDeposit[meta];
        dNTDEnergyDeposit = NTDEnergyDeposit[meta];
	dMuEnergyDeposit = MuEnergyDeposit[meta];
        dOtherEnergyDeposit = OtherEnergyDeposit[meta];

        dChannel = i;
        if( dEnergy>0.) {
          /*		G4cout << dChainNumber << "\t" << dDTime << "\t" << PTime[0]<< "\t" << i << "\t" << dEnergy;
                for( G4int mm=0; mm<TNumberOfParticles; mm++) G4cout << "\t" << dPEnergy[mm];
                G4cout << G4endl; */	          
          if( outMicro ) {
            for( G4int k=0; k<TNumberOfInteractions; k++ )  dNInt[k] = nInt[i][k][meta];
            if( PrintInfo )	G4cout << G4endl << G4endl << ">>>>>>>>> " << i << " E " << dEnergy 
              << " First " << dFirstInteraction 
                << " Foto " << dNInt[0] 
                << " Compton " << dNInt[1] 
                << " Pairs " << dNInt[2] 
                << " AverageDistance " << dAvDist 
                << " NumberOfClusters " << dNClus 
                << " TrackLength " << dTLen << G4endl << G4endl ;
          }				
          ROOTTree->Fill();
          LastChainTime = PTime[0] = 0.;
        }
        ClearRootVariables();
      }
    }
  }

  if( OutType & 2 && xyz_out ) {
    dChannel=-1;
    dDTime=0;
    ROOTTree->Fill();
  }
}	

TTree *qshieldsConfigurator::CreateROOTTree()
{
  ROOTTree = new TTree("qtree", "Tree with MC Data");

  ROOTTree->Branch("Channel",    			&dChannel,    		"Channel/I");
  ROOTTree->Branch("ChainNumber",			&dChainNumber,		"ChainNumber/l");
  ROOTTree->Branch("Time", 	    			&dDTime,				"Time/D");
  ROOTTree->Branch("DepositedEnergy",     			&dEnergy,			"Energy/D");
  ROOTTree->Branch("DepositedEnergyByParticle",				&dPEnergy,			"PEnergy[12]/D");
  ROOTTree->Branch("EnergyDepositedInCopper",				&dCuEnergyDeposit,			"CuEnergyDeposit/D");
  ROOTTree->Branch("EnergyDepositedInPTFE",				&dPTFEEnergyDeposit,			"PTFEEnergyDeposit/D");
  ROOTTree->Branch("EnergyDepositedInNTD",				&dNTDEnergyDeposit,			"NTDEnergyDeposit/D");
  ROOTTree->Branch("EnergyDepositedInMu",                               &dMuEnergyDeposit,                      "MuEnergyDeposit/D");
  ROOTTree->Branch("ParticleName",	"std::string",			&dParticleName);
  ROOTTree->Branch("DaughterName",	"std::string",			&dDaughterName);
  ROOTTree->Branch("ParticleInputEnergy",	&dEParticle,		"ParticleInputEnergy/D");
  ROOTTree->Branch("ParticlePositionX",		&dPositionX,	"PositionX/D");
  ROOTTree->Branch("ParticlePositionY",		&dPositionY,	"PositionY/D");
  ROOTTree->Branch("ParticlePositionZ",		&dPositionZ,	"PositionZ/D");
  ROOTTree->Branch("ParticleDirectionX",	&dDirectionX,	"DirectionX/D");
  ROOTTree->Branch("ParticleDirectionY",	&dDirectionY,	"DirectionY/D");
  ROOTTree->Branch("IsNeutron", &IsNeutron, "IsNeutron/I");
  if( outMicro) { 
    ROOTTree->Branch("PhotoInteractions",		&dNInt[0],			"PhotoInteractions/I");
    ROOTTree->Branch("ComptonInteractions",		&dNInt[1],			"ComptonInteractions/I");
    ROOTTree->Branch("PairInteractions",		&dNInt[2],			"PairInteractions/I");
    ROOTTree->Branch("FirstInteraction",	"std::string",			&dFirstInteraction);
    ROOTTree->Branch("AverageDistance",	&dAvDist,	"AverageDistance/D");
    ROOTTree->Branch("NumberOfClusters",	&dNClus,	"NumberOfClusters/I");
    ROOTTree->Branch("TotalTrackLength",	&dTLen,	"TotalTrackLength/D");
  }
  if( SType<0 ) ROOTTree->Branch("ParticleDepth",	&dDepth,	"Depth/D");
  //  Long64_t SavFreq = sizeof(*ROOTTree)*1000;
  //  ROOTTree->SetAutoSave( SavFreq );
  return( ROOTTree );
}

void qshieldsConfigurator::CuoreMap( std::string pss, bool yn ) 
{
  G4int pst;
  G4int numel; 

  pst = 1;
  for( unsigned int pen=1; pen<pss.length(); pen++ ) {
    if( pss.substr( pen,1 ) == std::string(",") || pss.substr( pen,1 ) == std::string(")") ) {
      numel = atoi( pss.substr( pst, pen-pst ).c_str() );
      if(yn) GraphicMask[numel/32] |= (1 << (numel-1)%32);
      else GraphicMask[numel/32] &= ~(1 << (numel-1)%32);
      pst = pen+1;
    }
  }
}

