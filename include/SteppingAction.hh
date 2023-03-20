#ifndef SteppingAction_h
#define SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"

class DetectorConstruction;
class AnalysisManager;

//---------------------------------------------------------------------------

class SteppingAction : public G4UserSteppingAction
{
  public:
  SteppingAction(DetectorConstruction*);
   ~SteppingAction();

    void UserSteppingAction(const G4Step*);
    
  private:
    DetectorConstruction*   detector;
    AnalysisManager*        anaManager;
};

//---------------------------------------------------------------------------

#endif
