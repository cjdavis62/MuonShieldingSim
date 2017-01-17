// vers. 5.0: geant4 physics lists
#ifndef qshieldsEventAction_h
#define qshieldsEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"


class qshieldsEventAction : public G4UserEventAction
{
  public:
    qshieldsEventAction();
    virtual ~qshieldsEventAction();

  public:
    virtual void   BeginOfEventAction(const G4Event*);
    virtual void   EndOfEventAction(const G4Event*);
    
};

#endif

    
