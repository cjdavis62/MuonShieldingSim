#ifndef G4DistributedSource_h
#define G4DistributedSource_h 1

#include "G4Transform3D.hh"
#include "G4VSolid.hh"

class G4DistributedSource 
{
  public:
  
	G4DistributedSource(G4VSolid *fSolid, G4Transform3D fT3D) { 
		aSolid = fSolid; 
		aT3D = fT3D; 
		Area = fSolid->GetSurfaceArea(); 
		Volume = fSolid->GetCubicVolume(); 
	}
    ~G4DistributedSource();

	G4VSolid * getSolid() {return aSolid;}
	G4Transform3D getT3D() {return aT3D;}
	G4double getArea() {return Area;}
	G4double getVolume() {return Volume;}
	void addArea( G4double fArea ) { Area += fArea;}
	void addVolume( G4double fVolume ) { Volume += fVolume;}

  private:
        G4VSolid* aSolid;
		G4Transform3D aT3D;
		G4double Area, Volume;
	
};

#endif

