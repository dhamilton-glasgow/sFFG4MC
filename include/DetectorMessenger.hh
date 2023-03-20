#ifndef DetectorMessenger_h
#define DetectorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

//---------------------------------------------------------------------------

class DetectorConstruction;
class G4UIdirectory;
class G4UIcommand;
class G4UIcmdWithADoubleAndUnit;

//---------------------------------------------------------------------------

class DetectorMessenger: public G4UImessenger
{
public:
  DetectorMessenger(DetectorConstruction*);
  ~DetectorMessenger();
  
  void SetNewValue(G4UIcommand*, G4String);
    
private:
  DetectorConstruction*        fDetector;
  G4UIdirectory*               fDetectorDir;

  G4UIcmdWithADoubleAndUnit*   fNPSAngleCmd;  
  G4UIcmdWithADoubleAndUnit*   fNPSDistanceCmd;  
  G4UIcmdWithADoubleAndUnit*   fHCALAngleCmd;  
  G4UIcmdWithADoubleAndUnit*   fHCALDistanceCmd;  
  G4UIcmdWithADoubleAndUnit*   fShieldThicknessCmd;  
  G4UIcmdWithADoubleAndUnit*   fWindowThicknessCmd;  
  G4UIcmdWithADoubleAndUnit*   fTargetLengthCmd;  

  G4UIcommand*                 fUpdateCmd;
};

#endif

//---------------------------------------------------------------------------

