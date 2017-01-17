#ifndef qshieldsSD_h
#define qshieldsSD_h 1

#include "G4VSensitiveDetector.hh"
#include "globals.hh"

class qshieldsDetectorConstruction;
class G4HCofThisEvent;
class G4Step;

class qshieldsSD : public G4VSensitiveDetector
{
  public:
  
      qshieldsSD(const G4String, qshieldsDetectorConstruction* );
     ~qshieldsSD();

      void Initialize(G4HCofThisEvent*);
      G4bool ProcessHits(G4Step*,G4TouchableHistory*);
      void EndOfEvent(G4HCofThisEvent*);
      void clear();
      void DrawAll();
      void PrintAll();

  private:
//        qshieldsDetectorConstruction* Detector;
	
};

#endif

