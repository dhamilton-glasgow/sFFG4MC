#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleDefinition.hh"
#include "G4Event.hh"

//---------------------------------------------------------------------------

class OutputManager;
class PrimaryGeneratorAction;

//---------------------------------------------------------------------------

class EventAction : public G4UserEventAction
{
  public:
  EventAction( OutputManager*, PrimaryGeneratorAction* );
   ~EventAction();

  public:
    void     BeginOfEventAction(const G4Event*);
    void     EndOfEventAction(const G4Event*);

  private:
  OutputManager*        fOutManager;
  PrimaryGeneratorAction* fPGA;

};
#endif

//---------------------------------------------------------------------------

    
