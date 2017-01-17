#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VisExtent.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4ThreeVector.hh"
#include "qshieldsConfigurator.hh"
#include "qshieldsDetectorConstruction.hh"

#include <vector>
#include <map>

class qshieldsDebug
{
	public:
		 qshieldsDebug();
		~qshieldsDebug();

	private:
		G4ThreeVector GetPointInVolume( const G4VisExtent VisExt );
        int  CheckThisVolume(G4VPhysicalVolume* vol);

        void CheckAllVolumes();
        void CheckSomeVolumes();

        std::vector<G4VPhysicalVolume*> volumes  ;
        std::vector<G4Transform3D     > transform;
        qshieldsConfigurator* theC;
        qshieldsDetectorConstruction* theD;
        G4ThreeVector Position;
        int nEvents;
        int v1, v2;

        std::map<int, const char*> vNames;
};
