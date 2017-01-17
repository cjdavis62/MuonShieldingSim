// First version for geant4.7.1 
#include "qshieldsPrimaryGeneratorAction.hh"
#include "qshieldsDetectorConstruction.hh"
#include "qshieldsConfigurator.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4VisExtent.hh"
#include "Randomize.hh"
#include <fstream>

#include "G4Gendec.hh"

qshieldsPrimaryGeneratorAction::qshieldsPrimaryGeneratorAction()
	:k_particle(0),b_particle(0),PreviousChainTime(0),nSpec(0), pDecay(0) 
{
  qshieldsConfigurator* theC = qshieldsConfigurator::GetInstance();
  if( theC->MuPar != G4ThreeVector() ) aMuon = new G4GSMu(theC->MuPar.x(),theC->MuPar.y(),theC->MuPar.z(),G4ThreeVector()); 
  MuRadius = theC->MuPar.x()*cm;
  std::ifstream iSpecFile;
  if( theC->iSpec!=0 ) {
    G4double xv,yv,yvs=0.,xvs=0;
    iSpecFile.open ( theC->SpecInpFile.c_str() );
    assert( iSpecFile );
    while( !iSpecFile.eof() ) {
     iSpecFile >> xv >> yv;
     if( theC->EDistrUniform <= 0 || 
       (xv>theC->PEnergy[k_particle] && xv<theC->EDistrUniform) ) {
//      if( nSpec<=0 ) {vX.push_back(xvs); vY.push_back(yvs); nSpec++;}
      if( nSpec<=0 || xv>xvs ) {vX.push_back(xv); vY.push_back(yvs); nSpec++;}
      yvs += yv;
     }
     if( xv>theC->EDistrUniform && xvs<theC->EDistrUniform) 
     	{vX.push_back(xv); vY.push_back(yvs+yv); nSpec++;}
     xvs = xv;
    }
     for( G4int jj=0; jj<nSpec; jj++)  G4cout << "Spectrum Entry n. "<< jj+1 << " " << vX[jj] << " " << vY[jj] << G4endl;
  }
  if( theC->DoubleBeta > 0 ) DoubleBetaDistr( theC->DoubleBeta, theC->QBeta, theC->fileE, theC->fileA );

  if( theC->GenDecAN ) ds.index=0;
  else if( theC->GenDec ) {
/*------------- the new logic of decay chains generation ----------------*\
The user has to provide a list of unstable nuclei correctly ordered. 
For each of them the following must be specified: 
A Z f id n_succ succ_list
A: mass number
Z: atomic number
Exc: Level excitation energy
f: branching ratio (%)
id: unique identifier in the chain (numbered starting from 1)
n_succ: number of successors
succ_list: list of successors (according to their id)

The nuclei are generated in sequence following a whole chain including branches.
When the time difference with respect to the previous decay in the chain is
less than a user specified time threshold, the position is maintained and the 
corresponding time written to file.
\*-----------------------------------------------------------------------*/
    std::ifstream iChainFile;
    struct UNSISO sDecay;
    G4int nprz;
    iChainFile.open ( theC->GendecFile.c_str() );
    assert( iChainFile );
    nDecay = 0;
    pDecay = -1;
    wtmp = 0.;
    while( !iChainFile.eof() ) {
     sDecay.A = 0;
     iChainFile >> sDecay.A;
     if( sDecay.A <= 0 ) break;
     iChainFile >> sDecay.Z;
     iChainFile >> sDecay.Exc;
     iChainFile >> sDecay.Weight; sDecay.Weight /= 100.;
     iChainFile >> sDecay.id;
     iChainFile >> sDecay.nsucc;
     nDecay++;
     assert( sDecay.id == nDecay );
     Chain.push_back(sDecay);
     if( sDecay.nsucc ==  0 ) pDecay = sDecay.id-1;
     for( G4int nn=0; nn<sDecay.nsucc; nn++ ) {
       iChainFile >> nprz; Chain[nDecay-1].succ.push_back(nprz-1);
     }
    }
    for(G4int jj=0; jj<nDecay; jj++) {
     G4cout << jj << " "  << Chain[jj].A << " "
     << Chain[jj].Z << " " << Chain[jj].Weight << " " << Chain[jj].nsucc;
     for(G4int kk=0; kk<Chain[jj].nsucc; kk++ ) G4cout << " " << Chain[jj].succ[kk];
     G4cout << G4endl;
    }
  }
  particleGun  = new G4ParticleGun(1);
}


qshieldsPrimaryGeneratorAction::~qshieldsPrimaryGeneratorAction()
{
  delete particleGun;
}

void qshieldsPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) {

static G4ThreeVector Position, Direction, AbsorberCenter;
G4ThreeVector sDirection;
G4double px,py,pz,Mass=0.,kEnergy,kExc,Energy,Momentum,Charge=0.,ExtEnergy;
G4String particleName="";
G4String Type="";
static G4double EBeta[2];
static G4ThreeVector DirBeta[2];
int rcharge=-2, fcharge=0, numgen=0; 
double depth=0.0;
G4double tx,ty,tz;
 
  qshieldsConfigurator* theC = qshieldsConfigurator::GetInstance();
  qshieldsDetectorConstruction* theD = theC->CDet ;
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
#ifdef G4GT10
  G4IonTable* ionTable = G4IonTable::GetIonTable();
#endif
  G4Gendec* theG;
  if( theC->GenDecAN ) theG = G4Gendec::GetInstance();

// Prima di tutto la posizione
 if( theC->MuPar.x() > 0. ) { // muons on a theC->MuPar.x() emisphere
   k_particle %= theC->n_particles;
   tz = MuRadius+1;
   while(tz > MuRadius) {
    aMuon->SetValue();
    Position = aMuon->GetPosition()*cm;
    Direction = aMuon->GetDirection();
	tz = -Position.z()/Direction.z();
	tx = Position.x()+tz*Direction.x();
	ty = Position.y()+tz*Direction.y();
	tz = sqrt(tx*tx+ty*ty);
   }
   Position += theC->BottomCenter;
   kExc = 0.;
   kEnergy = aMuon->GetEnergy();
   qMu = aMuon->GetCharge();
   ExtEnergy=aMuon->GetExt();
   theC->ResetChain();
   if( theC->qshieldsTest ){
    for(G4int i=0; i<theC->NumqshieldsTest;i++ ) {
      aMuon->SetValue();
      Position = aMuon->GetPosition();
      Direction = aMuon->GetDirection();
      kEnergy = aMuon->GetEnergy();
      qMu = aMuon->GetCharge();
      theC->DataFile << i << "\t" 
	   << Position.x() << "\t" << Position.y() << "\t" << Position.z() << "\t" 
	   << Direction.x() << "\t" << Direction.y() << "\t" << Direction.z() << "\t" 
	   << kEnergy << "\t" << qMu << G4endl; 
    }
    exit(1);
   }
 }
 else {
  qMu = 1;
  if( (k_particle %= theC->n_particles) == 0 ) {
   if ( theC->SType == 0 ) { // Dimensionless
	 Position = G4ThreeVector( 0,0,theD->ZMount); // PointLike: position relative to detector center
    if( theC->SMode == 0 )  // Pointlike
	 Position += theC->SPosition[ static_cast<G4int>( G4UniformRand()*theC->SPosition.size() ) ];
//G4cout << "PrimaryG - Point-like source: " << theC->SPosition << " " << theD->ZMount << G4endl;
    if( theC->SMode == 1 ) { // Wire
	 G4double rval = G4UniformRand()*theC->SGNorm;
	 G4double sval = 0.;
	 G4int i;
	 for( i=0; i<(int)theC->SWDirCos.size(); i++ ) if( rval < (sval += theC->SWDirCos[i].z()) ) break;
	 G4ThreeVector DDirection = theC->SWDirCos[i];
	 DDirection.setZ( sqrt(1. - DDirection.x()*DDirection.x() - DDirection.y()*DDirection.y()) );
     Position += theC->SPosition[i] + DDirection*(1. -2.*G4UniformRand())*theC->SWDirCos[i].z();
    } 
    else if( theC->SMode == 2 ) { // Disc
	 G4double rval = G4UniformRand()*theC->SGNorm;
	 G4double sval = 0.;
	 G4int i;
	 for( i=0; i<(int)theC->SWDirCos.size(); i++ ) if( rval < (sval+=theC->SWDirCos[i].y()*theC->SWDirCos[i].y()) ) break;
     sval = sqrt(G4UniformRand())*theC->SWDirCos[i].y(); 
     G4double Phi = G4UniformRand()*twopi;
	 G4ThreeVector DDirection;
	 if( theC->SWDirCos[i].x() == 1. ) DDirection = G4ThreeVector(0,sval*cos(Phi),sval*sin(Phi));
	 else if( theC->SWDirCos[i].x() == 2. ) DDirection = G4ThreeVector(sval*sin(Phi),0,sval*cos(Phi));
	 else if( theC->SWDirCos[i].x() == 3. ) DDirection = G4ThreeVector(sval*cos(Phi),sval*sin(Phi),0);
     Position += theC->SPosition[i] + DDirection;
    }
    else if( theC->SMode == 3 ) { // Sphere
	 G4double rval = G4UniformRand()*theC->SGNorm;
	 G4double sval = 0.;
	 G4int i;
	 for( i=0; i<(int)theC->SWDirCos.size(); i++ ) if( rval < (sval+=theC->SWDirCos[i].x()*theC->SWDirCos[i].x()*theC->SWDirCos[i].x()) ) break;
     G4double cosTheta = 1. - 2.*G4UniformRand();
     G4double sinTheta = sqrt(1. - cosTheta*cosTheta);
     G4double Phi = G4UniformRand()*twopi;
     Position += theC->SPosition[i] + theC->SWDirCos[i].x()*G4ThreeVector( sinTheta*cos(Phi), sinTheta*sin(Phi), cosTheta );
    }
    else if( theC->SMode == 4 ) { // Tube
	 G4double rval = G4UniformRand()*theC->SGNorm;
	 G4double sval = 0.;
	 G4int i;
	 for( i=0; i<(int)theC->SWDirCos.size(); i++ ) if( rval < (sval+=theC->SWDirCos[i].y()*theC->SWDirCos[i].z()) ) break;
     sval = (1. - 2.*G4UniformRand())*theC->SWDirCos[i].z();
     G4double Phi = G4UniformRand()*twopi;
	 G4ThreeVector TDirection;
	 if( theC->SWDirCos[i].x() == 1. ) TDirection = G4ThreeVector(sval,theC->SWDirCos[i].y()*cos(Phi),theC->SWDirCos[i].y()*sin(Phi));
	 else if( theC->SWDirCos[i].x() == 2. ) TDirection = G4ThreeVector(theC->SWDirCos[i].y()*sin(Phi),sval,theC->SWDirCos[i].y()*cos(Phi));
	 else if( theC->SWDirCos[i].x() == 3. ) TDirection = G4ThreeVector(theC->SWDirCos[i].y()*cos(Phi),theC->SWDirCos[i].y()*sin(Phi),sval);
     Position += theC->SPosition[i] + TDirection;
    }
   } else if ( theC->SType > 0 ) { // Cryostat element Volume
	  double vtmp = G4UniformRand()*theD->totVolume;
	  G4int is;
	  for( is=0; is<theD->SourceN; is++ ) if( theD->aSource[is]->getVolume() > vtmp ) break;
	  G4VisExtent VExt = theD->aSource[is]->getSolid()->GetExtent();
	  Position = GetPointInVolume(VExt);
	  while( theD->aSource[is]->getSolid()->Inside(Position)!=kInside ) Position = GetPointInVolume(VExt); 
	  Position *= theD->aSource[is]->getT3D().getRotation();
	  Position += theD->aSource[is]->getT3D().getTranslation(); 
   } else { // Setup element Surface
	  G4int is;
	  G4ThreeVector sDir;
	  double stmp = G4UniformRand()*theD->totArea;
	  for( is=0; is<theD->SourceN; is++ ) if( theD->aSource[is]->getArea() > stmp ) break;
	  if( theC->SurfExpGen > 0 ) depth=-log(1.-G4UniformRand())*theC->SDepth;
	  else ( theC->SDepth < 0. ? depth = -theC->SDepth*G4UniformRand(): depth = theC->SDepth );
	  Position = theD->aSource[is]->getSolid()->GetPointOnSurface(); 
	  sDir = theD->aSource[is]->getSolid()->SurfaceNormal(Position);
	  if( theC->PrintInfo ) G4cout << "PrimaryGenerator surface generation - depth and additive offset: " << depth << " " << fabs(theC->SOffs) << G4endl;
  	  depth += fabs(theC->SOffs);
	  Position -= depth*sDir; 
	  if( depth > 0 ) {
	   while( theD->aSource[is]->getSolid()->DistanceToOut(Position) < depth*0.99 ) 
	   {
	      if( theC->SurfExpGen > 0 ) depth=-log(1.-G4UniformRand())*theC->SDepth;
	      else ( theC->SDepth < 0. ? depth = -theC->SDepth*G4UniformRand(): depth = theC->SDepth );
	  	Position = theD->aSource[is]->getSolid()->GetPointOnSurface(); 
	  	sDir = theD->aSource[is]->getSolid()->SurfaceNormal(Position);
	  	Position -= depth*sDir;  
	   }
	  }
	  Position *= theD->aSource[is]->getT3D().getRotation();
	  Position += theD->aSource[is]->getT3D().getTranslation(); 
   }

// Generate Radioactive decay
   FixPos = false;
   if( theC->GenDecAN ) {
	     while( (numgen = theG->dec_start_( &ds )) >= 0 && 
	          !( theC->ParticleMask & (1 << theC->pId[ ds.typ ]) )) ;
	     if( ds.typ == 'r' )
	       theC->SetParticle( ds.epart, rcharge, ds.qpart, ds.zpart, 0.0 );
	     else
	       theC->SetParticle( ds.epart, ds.qpart, fcharge );
		 if( ds.index != theC->dChainNumber ) {
			theC->ResetChain(ds.index);
			PreviousChainTime = 0.;
		 }
		 else FixPos = true;
	     dTime = ds.tpart;
	 	 if( !FixPos || fabs(dTime) > theC->tThresh ) {
	      GenPosition = Position;
	     }
	     else {
	      Position = GenPosition;
	     }
		if( dTime > 1.e-6 ) ResetAlphaDir();
	    theC->SetChainTime( dTime );    
   } 
   else if( theC->GenDec ) {
	    if( theC->StableNucleus || pDecay>=nDecay || (theC->MassNumber==Chain[nDecay-1].A && theC->AtomicNumber==Chain[nDecay-1].Z) ) {
			theC->ResetChain();
			pDecay = 1;
			theC->SetParticle( 0., -2, Chain[pDecay-1].A, Chain[pDecay-1].Z, 0. );
	 	    //G4cout << "\taGenerate: " << Chain[pDecay-1].A << " " << Chain[pDecay-1].Z << G4endl;
		}
	    else if( theC->MassNumber<Chain[nDecay-1].A ) {
			pDecay++;
			theC->SetParticle( 0., -2, Chain[pDecay-1].A, Chain[pDecay-1].Z, 0. );
	 	    //G4cout << "\tbGenerate: " << Chain[pDecay-1].A << " " << Chain[pDecay-1].Z << G4endl;
			FixPos = true;
		}
		else {
			pDecay++;
			theC->SetParticle( pE, -2, theC->MassNumber, theC->AtomicNumber, 0. );    
	 	    //G4cout << "\tcGenerate: " << theC->MassNumber << " " << theC->AtomicNumber << G4endl;
			FixPos = true;
		} 
   } 
   else  theC->ResetChain();

   if( theC->DoubleBeta > 0 ) {
    GenDoubleBeta( theC->DoubleBeta, EBeta, DirBeta, theC->QBeta );
    // G4cout << "Electron energies: " << EBeta[0] << " " << EBeta[1] << G4endl;
   }

   if( theC->PrintInfo ) G4cout << "PrimaryGnerator - event start position: " << theC->SType <<  " " << theC->SMode << 
  	(Position-G4ThreeVector( 0,0,theD->ZMount))*0.1 << " " 
        << theC->SType <<  " " << theC->SMode << G4endl;
  }

//  G4PrimaryVertex* aVertex = new G4PrimaryVertex( Position,particleGun->GetParticleTime() );

// ... quindi ciclo sulle particelle da generare: direzione, energia, ...

  if( theC->DoubleBeta > 0 && theC->PType[k_particle] == -1 ) {
    if( b_particle > 2 &&  k_particle < theC->n_particles ) {
       G4cout << " No double beta process with more than two electrons is allowed!" << G4endl;
       exit(1);
    }
    b_particle %= 2;
    kEnergy = EBeta[b_particle];
    Direction = DirBeta[b_particle];
    b_particle++;
  }
  else {
    kExc = theC->PEexc[k_particle]*MeV;
    if( theC->SAF > 0 ) {
      if( theC->PrintInfo ) G4cout << theC->SAngle << " " << Position; 
      Direction = GetRandomDirection ( theC->SAngle, Position );
      if( theC->PrintInfo ) G4cout << Direction << G4endl;
    }
    else if( theC->ggCorrelation > 0 ) {
      if( k_particle == 0 ) Direction = GetRandomDirection ();
      else Direction = GetRandomDirection (Direction);
    }
    else Direction = GetRandomDirection ();

    if( theC->EDistrUniform > 0 ) {
      kEnergy = theC->PEnergy[k_particle]+
        G4UniformRand()*(theC->EDistrUniform-theC->PEnergy[k_particle])  ;
//     G4cout << "Uniform " << kEnergy << G4endl;
    }
    else if( theC->EThermalAvr != 0.0 ) {
     G4double FVal;
     do {
      kEnergy = G4UniformRand()*theC->EThermalMax;
      FVal = G4UniformRand()*theC->ThermalMax;
     } while( FVal > Thermals( kEnergy,theC->EThermalAvr ) );
//     G4cout << "Thermal " << kEnergy << G4endl;
    }
    else if( theC->iSpec!=0 ) {
      G4double val = G4UniformRand()*vY[nSpec-1];
      G4int i;
      for( i=0; i<nSpec; i++ ) if( vY[i] > val ) break;
      kEnergy = vX[i-1] + (val-vY[i-1])*(vX[i]-vX[i-1])/(vY[i]-vY[i-1]);
      if( theC->PrintInfo ) G4cout << "Spectrum Energy: " << kEnergy << G4endl;
    }
    else
      kEnergy = theC->PEnergy[k_particle] ;
  }
 }
  if( theC->oSpec!=0 ) theC->SpecFile << kEnergy*1.e3 << G4endl;

// CICLO sul tipo di particella
   if( !theC->GenDec ) theC->SetType(theC->GetType(theC->PType[k_particle]), kEnergy);
   theC->SetCurrentParticle( theC->PType[k_particle] );
   G4ParticleDefinition* aParticleDefinition = 0;
// 1: FOTONE
  if( theC->PType[k_particle] == 0 ) {
   aParticleDefinition 
    = particleTable->FindParticle( particleName="gamma" );

   Mass = 0.;
   Charge = 0.;
  }
//2: CHARGEDGEANTINO
  else if( theC->PType[k_particle] == 5 ) {
   aParticleDefinition 
    = particleTable->FindParticle( particleName="chargedgeantino" );

   Charge = aParticleDefinition->GetPDGCharge();
   Mass = aParticleDefinition->GetPDGMass();
  }
//3: ELECTRON
  else if( theC->PType[k_particle] == -1 ) {
   aParticleDefinition 
    = particleTable->FindParticle( particleName="e-" );

   Charge = aParticleDefinition->GetPDGCharge();
   Mass = aParticleDefinition->GetPDGMass();
  }
//4: NEUTRON
  else if( theC->PType[k_particle] == 3 ) {
   aParticleDefinition 
    = particleTable->FindParticle( particleName="neutron" );
   Charge = aParticleDefinition->GetPDGCharge();
   Mass = aParticleDefinition->GetPDGMass();
  }
//5: POSITRON
  else if( theC->PType[k_particle] == +1 ) {
   aParticleDefinition 
    = particleTable->FindParticle( particleName="e+" );
   Charge = aParticleDefinition->GetPDGCharge();
   Mass = aParticleDefinition->GetPDGMass();
  }
//6: ALPHA
  else if( theC->PType[k_particle] == +2 ) {
 //   Mass = 4.002602*931.49432*MeV;
   if( theC->GenDecAN ) Direction = GetAlphaDir( Direction,-1 );
   aParticleDefinition 
    = particleTable->FindParticle( particleName="alpha" );
   Charge = aParticleDefinition->GetPDGCharge();
   Mass = aParticleDefinition->GetPDGMass();
  }
//6b: MUON
  else if( theC->PType[k_particle] == +6 ) {
   if( qMu > 0.0 ) 
   	aParticleDefinition 
   	 = particleTable->FindParticle( particleName="mu-" );
   else
   	aParticleDefinition 
   	 = particleTable->FindParticle( particleName="mu+" );
   Charge = aParticleDefinition->GetPDGCharge();
   Mass = aParticleDefinition->GetPDGMass();
  }
//10 CHARGED GEANTINO
  else if( theC->PType[k_particle] == 10 ) {
   	aParticleDefinition 
   	 = particleTable->FindParticle( particleName="geantino" );
   Charge = aParticleDefinition->GetPDGCharge();
   Mass = aParticleDefinition->GetPDGMass();
   kEnergy=1.;
  }
//7: NUCLEAR RECOIL
  else if( theC->PType[k_particle] == -2 ) {
    G4int BarionN = static_cast<G4int> (theC->PMass[k_particle]);
    G4int AtomicN = static_cast<G4int> (theC->PZeta[k_particle]);
    aParticleDefinition
#ifdef G4GT10
		= ionTable->GetIon(  AtomicN, BarionN, kExc );
#else
		= particleTable->GetIon(  AtomicN, BarionN, kExc );
#endif
    particleName = aParticleDefinition->GetParticleName();
    Charge   = 0.*eplus;
    particleGun->SetParticleCharge(Charge);
    Mass = aParticleDefinition->GetPDGMass();

	dTime = 0;
	if( theC->GenDecAN ) Direction = GetAlphaDir( Direction,1 );
	else if( theC->GenDec ) {
      G4double DecTime = aParticleDefinition->GetPDGLifeTime()/s;
      dTime = -DecTime*log(1.-G4UniformRand());
 	  if( !FixPos || fabs(dTime) > theC->tThresh ) {
       GenPosition = Position;
      }
      else if( theC->SMode>0 && theC->SMode<5 ) {
                Position = PreviousPosition;
      }
      else {
		// G4ThreeVector DiV = theC->LastPosition - PreviousPosition;
		// G4cout << "primary ... " << PreviousPosition << " " << theC->LastPosition << " " << sqrt(DiV.x()*DiV.x() + DiV.y()*DiV.y() + DiV.z()*DiV.z()) << G4endl;
		Position = theC->LastPosition;
      }
	  PreviousPosition = Position;
      theC->SetChainTime( dTime );    
	}
  }

// Fill ROOTtree
	theC->dDirectionX = Direction.x();
	theC->dDirectionY = Direction.y();
	theC->dPositionX = Position.x();
	theC->dPositionY = Position.y();
	theC->dPositionZ = Position.z();
	theC->dParticleName = std::string(particleName);
	theC->dEParticle = kEnergy;
	theC->dDepth = depth;
	theC->iNtot++;

  if( theC->PrintInfo ) 
     G4cout << "PrimaryGnerator - event summary: " << " " << particleName << " " << Mass << " " << Charge << " " <<
      kEnergy << " " << Position << " " << Direction << G4endl;

  if( theC->PrintInfo && theC->GenDec ) 
     G4cout << "qshields decay chain PARTICLE: " << theC->dChainNumber << " " << G4endl;
  
  Energy = Mass + kEnergy;
  Momentum = sqrt(Energy*Energy-Mass*Mass);
  px = Direction.x()*Momentum;
  py = Direction.y()*Momentum;
  pz = Direction.z()*Momentum;

  particleGun->SetParticleDefinition(aParticleDefinition);
  particleGun->SetParticleTime(0.0*ns);
  particleGun->SetParticlePosition(Position);
  particleGun->SetParticleMomentumDirection(Direction.unit());
  particleGun->SetParticleEnergy(kEnergy);

  particleGun->GeneratePrimaryVertex(anEvent);

  k_particle++;
}

G4ThreeVector qshieldsPrimaryGeneratorAction::GetRandomDirection() {

  G4ThreeVector retval;

  G4double CosTheta;   
  G4double SinTheta;
  
  G4double Phi; 
  G4double SinPhi;
  G4double CosPhi;
  
  G4double rand;

  //
  //	selecting a random direction
  //

  rand = G4UniformRand();

  CosTheta = 2.0*rand -1.0;
  SinTheta = sqrt (1.-CosTheta*CosTheta);
  rand = G4UniformRand();
  Phi = twopi*rand;
  SinPhi = sin (Phi);
  CosPhi = cos (Phi);
  retval.setX(SinTheta*CosPhi);
  retval.setY(SinTheta*SinPhi);
  retval.setZ(CosTheta);

  return retval;
}
G4ThreeVector qshieldsPrimaryGeneratorAction::GetRandomDirection(G4ThreeVector olDir) {

  G4ThreeVector retval;

// Generiamo l'angolo tra le due direzioni
  G4double cosTheta, sinTheta, fVal;   
  fVal = 0.;
  do {
   cosTheta = cos(twopi*G4UniformRand());
   fVal = 1.+1.62*(cosTheta*cosTheta/8.+pow(cosTheta,4)/24.);
  } while( fVal > 1.27*G4UniformRand() );
  sinTheta = sqrt (1.-cosTheta*cosTheta);

  G4double Phi, sinPhi, cosPhi;
  Phi = twopi*G4UniformRand();
  sinPhi = sin (Phi);
  cosPhi = cos (Phi);
  retval.setX(sinTheta*cosPhi);
  retval.setY(sinTheta*sinPhi);
  retval.setZ(cosTheta);
  
/*   cosTheta = olDir.z();
  sinTheta = sqrt (1.-cosTheta*cosTheta);
  if( sinTheta != 0. ) Phi = acos(olDir.x()/sinTheta);
  else Phi = 0.;
  retval.rotateY(acos(cosTheta));
  retval.rotateZ(Phi);
 */  
  retval.rotateUz(olDir);
//  G4cout << "cosTheta: " <<  retval.x()*olDir.x() + retval.y()*olDir.y() + retval.z()*olDir.z() << G4endl;

  return retval;
}

G4ThreeVector qshieldsPrimaryGeneratorAction::GetRandomDirection( G4ThreeVector SA, G4ThreeVector P ) {

G4int flag;
G4double X,Y,Z;

G4double Phi; 
G4double SinPhi;
G4double CosPhi;
G4double CosTheta=0.;
G4double SinTheta;
  
  flag = static_cast<G4int> ( SA.x() );
  if( flag <= 0 ) return( GetRandomDirection() );
  qshieldsConfigurator* theC = qshieldsConfigurator::GetInstance();

  
  //	selecting a random direction on solid angle SA

  Z = SA.y() + G4UniformRand()*( SA.z() - SA.y() );
  SinTheta = sqrt( 1. - Z*Z );

  Phi = twopi*G4UniformRand();
  SinPhi = sin (Phi);
  CosPhi = cos (Phi);
  X = SinTheta*CosPhi;
  Y = SinTheta*SinPhi;
  
  switch( flag ) {
     case 1:
       return( G4ThreeVector( Z,X,Y ) );
       break;
     case 2:
       return( G4ThreeVector( Y,Z,X ) );
       break;
     case 3:
       return( G4ThreeVector( X,Y,Z ) );
       break;
     case 12:
     case 21:
       CosTheta = 0.;
       SinTheta = 1.;
       CosPhi = P.x()/sqrt(P.x()*P.x()+P.y()*P.y());
       SinPhi = P.y()/sqrt(P.x()*P.x()+P.y()*P.y());
       if( theC->PrintInfo ) G4cout << " (" << CosPhi << "," << SinPhi << ",0)" ; 
       break;
     case 13:
     case 31:
       CosTheta = P.z()/sqrt(P.x()*P.x()+P.z()*P.z());
       SinTheta = P.x()/sqrt(P.x()*P.x()+P.z()*P.z());
       CosPhi = 1.;
       SinPhi = 0.;
       if( theC->PrintInfo ) G4cout << " (" << SinTheta << ",0," << CosTheta << ")" ; 
       break;
     case 23:
     case 32:
       CosTheta = P.z()/sqrt(P.y()*P.y()+P.z()*P.z());
       SinTheta = P.y()/sqrt(P.y()*P.y()+P.z()*P.z());
       CosPhi = 0.;
       SinPhi = 1.;
       if( theC->PrintInfo ) G4cout << " (0," << SinTheta << "," << CosTheta << ")" ; 
       break;
     case 123:
     case 132:
     case 213:
     case 231:
     case 312:
     case 321:
       CosTheta = P.z()/P.mag();
       SinTheta = sqrt( 1. - CosTheta*CosTheta );
       CosPhi = P.x()/sqrt(P.x()*P.x()+P.y()*P.y());
       SinPhi = P.y()/sqrt(P.x()*P.x()+P.y()*P.y());
       if( theC->PrintInfo ) G4cout << " " << P/P.mag() ; 
       break;
     default:
       return( GetRandomDirection() );
  }
  return( G4ThreeVector( X*CosPhi*CosTheta-Y*SinPhi+Z*CosPhi*SinTheta, 
       			X*CosTheta*SinPhi+Y*CosPhi+Z*SinPhi*SinTheta,
			-X*SinTheta+Z*CosTheta ) );
}

// Isotropic generator in a cylindrical shell
G4ThreeVector qshieldsPrimaryGeneratorAction::GetRandomShellPosition (G4double Rmin, G4double Rmax, G4double Alt) {
G4double rand;
  
  rand = G4UniformRand();
  G4double Zeta = (rand-0.5)*Alt;
  rand = G4UniformRand();
  G4double Phi = rand*twopi;
  rand = G4UniformRand();
  G4double R = sqrt( rand*(Rmax*Rmax - Rmin*Rmin ) + Rmin*Rmin );
  G4double X = R*cos ( Phi );
  G4double Y = R*sin ( Phi );
  G4ThreeVector position( X, Y ,Zeta );
  return (position);
}

// Double beta decays single electron distributions
void qshieldsPrimaryGeneratorAction::GenDoubleBeta (G4int type, G4double* E, G4ThreeVector* A, G4double Q) {
const G4double em=0.511*MeV;
G4double rand,rval,pval,et=0.,e1,e2,ec,ev,p1,p2;
G4double csign=1.;

  switch( type ) {
   case 1:	// neutrinoless m-neutrino
      do {	
	E[0] = G4UniformRand()*Q;
        E[1] = Q - E[0];
	e1 = E[0]/em + 1.;
	e2 = E[1]/em + 1.;
	ev= e1*e1 * e2*e2;
	rand = G4UniformRand();
      } while( rand > ev/F12Max );
      break;
   case 2:	// 2 neutrinos
      do {	
	E[0] = G4UniformRand()*Q;
	E[1] = G4UniformRand()*Q;
	if( E[0]+E[1] > Q ) {
	  E[0] = Q - E[0];
	  E[1] = Q - E[1];
	}
	e1 = E[0]/em + 1.;
	e2 = E[1]/em + 1.;
	ec = (Q-E[0]-E[1])/em;
	ev = e1*e1 * e2*e2 * ec*ec*ec*ec*ec;
	rand = G4UniformRand();
      } while( rand > ev/F12Max );	// Distribuzione energie cinetiche degli elettroni
      break;
   case 3:	// 1 Majoron
     do {	
	E[0] = G4UniformRand()*Q;
	E[1] = G4UniformRand()*Q;
	if( E[0]+E[1] > Q ) {
	  E[0] = Q - E[0];
	  E[1] = Q - E[1];
	}
	e1 = E[0]/em + 1.;
	e2 = E[1]/em + 1.;
	ec = (Q-E[0]-E[1])/em;
	ev = e1*e1 * e2*e2 * ec;
	rand = G4UniformRand();
      } while( rand > ev/F12Max );	// Distribuzione energie cinetiche degli elettroni
      break;
   case 4:	// 2 Majorons
     do {	
	E[0] = G4UniformRand()*Q;
	E[1] = G4UniformRand()*Q;
	if( E[0]+E[1] > Q ) {
	  E[0] = Q - E[0];
	  E[1] = Q - E[1];
	}
	e1 = E[0]/em + 1.;
	e2 = E[1]/em + 1.;
	ec = (Q-E[0]-E[1])/em;
	ev = e1*e1 * e2*e2 * ec*ec*ec;
	rand = G4UniformRand();
      } while( rand > ev/F12Max );	// Distribuzione energie cinetiche degli elettroni
      break;
   case 5:	// neutrinoless RH
      do {	
	E[0] = G4UniformRand()*Q;
        E[1] = Q - E[0];
	e1 = E[0]/em + 1.;
	e2 = E[1]/em + 1.;
	et = e1 - e2;
	ev= e1*e1 * e2*e2 * et*et;
	rand = G4UniformRand();
      } while( rand > ev/F12Max );
      csign = -1.;
      break;
   case 6:	// 2 neutrinos, 2+
     csign = 1./3.;
     do {	
	E[0] = G4UniformRand()*Q;
	E[1] = G4UniformRand()*Q;
	if( E[0]+E[1] > Q ) {
	  E[0] = Q - E[0];
	  E[1] = Q - E[1];
	}
	e1 = E[0]/em + 1.;
	e2 = E[1]/em + 1.;
	ec = (Q-E[0]-E[1])/em;
	p1 = (E[0] - E[1])/em;
	ev = e1*e1 * e2*e2 * p1*p1 * ec*ec*ec*ec*ec*ec*ec;
	rand = G4UniformRand();
      } while( rand > ev/F12Max );	// Distribuzione energie cinetiche degli elettroni
      break;
   case 7:	// 0 neutrinos, 2+ RH
     do {	
	et = 1. -2.*G4UniformRand();	// cos( theta )
	E[0] = G4UniformRand()*Q;
	E[1] = Q-E[0];
	e1 = E[0]/em + 1.;
	e2 = E[1]/em + 1.;
        p1 = sqrt(E[0]*(E[0]+2.*em));
        p2 = sqrt(E[1]*(E[1]+2.*em));
	ev = e1*e2 * ( 3.*p1*p1*p2*p2*et*et -
	  p1*p2*et*( 10.*(e1*e2+1.)+p1*p1+p2*p2 ) + 5.*(e1*e2+1.)*(p1*p1+p2*p2) - p1*p1*p2*p2 );
	rand = G4UniformRand();
      } while( rand > ev/F12Max );	// Distribuzione energie cinetiche degli elettroni
      break;
   case 8:	// DBD from file
	  G4int nn1,nn2;
	  G4double de1,de2;
	  do {
		rand = F12Max*G4UniformRand();
		E[0] = emax0*G4UniformRand();
		nn1 = E[0]/estep0;
		de1 = E[0] - nn1*estep0;
		E[1] = emax0*G4UniformRand();
		nn2 = E[1]/estep0;
		de2 = E[1] - nn2*estep0;
		rval = g0[nn1][nn2]+(g0[nn1+1][nn2+1]-g0[nn1][nn2])*(de1+de2)/estep0;
	    nn1 = E[0]/estep1;
	    de1 = E[0] - nn1*estep1;
	    nn2 = E[1]/estep1;
	    de2 = E[1] - nn2*estep1;
	    pval = g1[nn1][nn2]+(g1[nn1+1][nn2+1]-g1[nn1][nn2])*(de1+de2)/estep1;
	  } while( rand> rval || pval==0. );
	  p1 = sqrt(E[0]*(E[0]+2.*em));
	  p2 = sqrt(E[1]*(E[1]+2.*em));

	G4double discr, discr2;
	discr2 = rval*rval - 4*pval*(rval-pval-2.*rval*G4UniformRand());
	discr = sqrt(fabs(discr2));
	et = (-rval+discr)/2./pval;
//	et2 = (-rval-discr)/2./pval;
//	G4cout <<  E[0]*1000 << " " << E[1]*1000 << " " << et << " " << (E[0]+E[1])*1000 << " " << G4endl;   
//	G4cout <<  E[0] << " " << E[1] << " " << et << " "<< et2 << " "<< discr2 << " "
//	<< nn1 << " " << nn2 << " " << g0[nn1][nn2] << " " << g1[nn1][nn2] << " " << rval << " " << pval << G4endl;   
//	  et=csign*(p1*p2/((E[0]+em)*(E[1]+em)));	// beta_1*beta_2 => (1-b_1 b_2 cos_theta)
//	  et=(1.-sqrt(1.-2.*et*(1-et/2.)*G4UniformRand()))/et;  // Coseno di  theta
      break;
   default:
      exit(1);
  }
  
  if( type < 7 ) {
   p1 = sqrt(E[0]*(E[0]+2.*em));
   p2 = sqrt(E[1]*(E[1]+2.*em));
   et=csign*(p1*p2/((E[0]+em)*(E[1]+em)));
   et=(1.-sqrt((1.+et)*(1.+et)-4.*et*G4UniformRand()))/et;  // Coseno di  theta
  }
// G4double CosTheta = et;
  A[0] = A[1] = G4ThreeVector( 0.,0.,1. );
  A[1].rotateY( acos(et) );
  A[1].rotateZ( G4UniformRand()*2.*M_PI );
// 17 Jan 2013: Jeffrey Wilson => reversed et and ec
  ec = G4UniformRand()*2.*M_PI;
  et = acos( 1. - 2.*G4UniformRand() );
  A[0].rotateY( et );
  A[0].rotateZ( ec );
  A[1].rotateY( et );
  A[1].rotateZ( ec );
  //G4cout << E[0] << " " << E[1] << " " << CosTheta << " " << ev/F12Max << G4endl;
  
}

void qshieldsPrimaryGeneratorAction::DoubleBetaDistr  ( G4int type, G4double Q, std::string file0, std::string file1 ) {
const G4double em=0.511*MeV;
G4double ek,ec,et;
std::ifstream inpFile;

  G4double T=Q/em;
  
  switch( type ) {
   case 1:	// 0 neutrinos, 0+, 2n, m_neutrino
      ek = T/2. + 1.;
      F12Max = ek*ek*ek*ek;
      break;
   case 2:	// 2 neutrinos, 0+, 2n
      ek = (2.*T-5.)/9.;
      et = ek+1.;
      ec = T-2.*ek;
      F12Max = et*et*et*et * ec*ec*ec*ec*ec;
      break;
   case 3:	// 0 neutrinos, 1 Majoron, 0+, 2n
      ek = (2.*T-2.)/6.;
      et = ek+1.;
      ec = T-2.*ek;
      F12Max = et*et*et*et * ec;
      break;
   case 4:	// 0 neutrinos, 2 Majorons, 0+, 2n
      ek = (2.*T-3.)/7.;
      et = ek+1.;
      ec = T-2.*ek;
      F12Max = et*et*et*et * ec*ec*ec;
      break;
   case 5:	// 0 neutrinos, 0+, 2n, RH
      ek = T + 1.;
      F12Max = ek*ek * T*T;
      break;
   case 6:	// 2 neutrinos, 2+, 2n
      ek = 4.*T-9.;
      ek = (ek+sqrt(ek*ek+88.*T))/22.;
      et = ek+1.;
      ec = T-ek;
      F12Max = et*et * ek*ek * ec*ec*ec*ec*ec*ec*ec;
      break;
   case 7:	// 0 neutrinos, 2+, 2n RH (max per K1=K2=T/2, cos(theta)=-1)
      ek = T/2. + 1.;
      ec = T/2.*(ek+1.);
      F12Max = ek*ek * (3.*ec*ec + ec*( 10.*(ek*ek+1.)+2.*ec ) + 10.*(ek*ek+1.)*ec - ec*ec );
      break;
	   case 8:	// double beta from files
	// electron distribution g0
		  inpFile.open(file0.c_str());
		  if(!inpFile){
		    cout << "[ERROR] Could not load g_0 file "<< file0 << G4endl;
		    exit(1);
		  }
		  G4int idx,idy;
		  G4double xv,yv,prob;
	// reads distribution parameters and counts elements
		  while (!inpFile.eof()) {
		    inpFile >> idx >> idy >> xv >> yv >> prob;
			if( nstep0 == 0 ) estep0 = xv;
			if( idx > nstep0 ) nstep0 = idx;
			if( idy > nstep0 ) nstep0 = idy;
			if(xv>emax0) emax0 = xv;
		  }
		nstep0 +=10;
		  g0 = new G4double*[nstep0];
		  for (G4int i = 0; i < nstep0; ++i)
		    g0[i] = new G4double[nstep0];
		  inpFile.close();
		  inpFile.open(file0.c_str());
		  while (!inpFile.eof()) {
		    inpFile >> idx >> idy >> xv >> yv >> prob;
		    // G4cout << "due: " << idx << " " << idy << " " << xv << " " << yv << " " << prob << G4endl;
			g0[idx-1][idy-1] = prob;
			if(prob>F12Max) F12Max = prob;
		  }
		  inpFile.close();
	// electron correlation g1
		  inpFile.open(file1.c_str());
		  if(!inpFile){
		    cout << "[ERROR] Could not load g_1 file "<< file1 << G4endl;
		    exit(1);
		  }
	// reads distribution parameters and counts elements
		  while (!inpFile.eof()) {
		    inpFile >> idx >> idy >> xv >> yv >> prob;
			if( nstep1 == 0 ) estep1 = xv;
			if( idx > nstep1 ) nstep1 = idx;
			if( idy > nstep1 ) nstep1 = idy;
			if(xv>emax1) emax1 = xv;
		  }
		nstep1+=10;
		  g1 = new G4double*[nstep1];
		  for (G4int i = 0; i < nstep1; ++i)
		    g1[i] = new G4double[nstep1];
		  //inpFile.seekg(0);
		  inpFile.close();
		  inpFile.open(file1.c_str());
		  while (!inpFile.eof()) {
		    inpFile >> idx >> idy >> xv >> yv >> prob;
		    // G4cout << idx << " " << idy << " " << xv << " " << yv << " " << prob << G4endl;
			g1[idx-1][idy-1] = prob;
		  }
		  inpFile.close();
			// for(G4int m=0;m<100000;m++) GenDoubleBeta( thec->DoubleBeta, ebbeta, dirbeta, thec->QBeta );
		  // exit(1);

	      break;
   default:
      ;
  }
}

G4ThreeVector qshieldsPrimaryGeneratorAction::GetPointInVolume ( const G4VisExtent VisExt ) {
  return( G4ThreeVector( VisExt.GetXmin()+(VisExt.GetXmax()-VisExt.GetXmin())*G4UniformRand(),
	VisExt.GetYmin()+(VisExt.GetYmax()-VisExt.GetYmin())*G4UniformRand(),
	VisExt.GetZmin()+(VisExt.GetZmax()-VisExt.GetZmin())*G4UniformRand() ) );
}


