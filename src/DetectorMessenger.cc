#include "DetectorMessenger.hh"
#include "DetectorConstruction.hh"

#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"

//---------------------------------------------------------------------------

DetectorMessenger::DetectorMessenger(DetectorConstruction* Detect)
  :fDetector(Detect)
{
  fDetectorDir  = new G4UIdirectory("/sFFG4MC/detector/");
  fDetectorDir->SetGuidance("Detector geometry control");

  fNPSAngleCmd = new G4UIcmdWithADoubleAndUnit("/sFFG4MC/detector/setEarmAngle",this);  
  fNPSAngleCmd->SetGuidance("Set Earm angle (value and unit)");
  fNPSAngleCmd->SetUnitCategory("Angle");

  fNPSDistanceCmd = new G4UIcmdWithADoubleAndUnit("/sFFG4MC/detector/setEarmDistance",this);  
  fNPSDistanceCmd->SetGuidance("Set Earm distance  (value and unit)");
  fNPSDistanceCmd->SetUnitCategory("Length");

  fNPSShieldThicknessCmd = new G4UIcmdWithADoubleAndUnit("/sFFG4MC/detector/setEarmShieldThickness",this);  
  fNPSShieldThicknessCmd->SetGuidance("Set Earm lead shield thickness  (value and unit)");
  fNPSShieldThicknessCmd->SetUnitCategory("Length");

  fHCALAngleCmd = new G4UIcmdWithADoubleAndUnit("/sFFG4MC/detector/setHarmAngle",this);  
  fHCALAngleCmd->SetGuidance("Set Harm angle  (value and unit)");
  fHCALAngleCmd->SetUnitCategory("Angle");

  fHCALDistanceCmd = new G4UIcmdWithADoubleAndUnit("/sFFG4MC/detector/setHarmDistance",this);  
  fHCALDistanceCmd->SetGuidance("Set Harm distance  (value and unit)");
  fHCALDistanceCmd->SetUnitCategory("Length");

  fHCALShieldThicknessCmd = new G4UIcmdWithADoubleAndUnit("/sFFG4MC/detector/setHarmShieldThickness",this);  
  fHCALShieldThicknessCmd->SetGuidance("Set Harm lead shield thickness  (value and unit)");
  fHCALShieldThicknessCmd->SetUnitCategory("Length");

  fWindowThicknessCmd = new G4UIcmdWithADoubleAndUnit("/sFFG4MC/detector/setWindowThickness",this);  
  fWindowThicknessCmd->SetGuidance("Set scattering chamber window thickness (value and unit)");
  fWindowThicknessCmd->SetUnitCategory("Length");

  fTargetLengthCmd = new G4UIcmdWithADoubleAndUnit("/sFFG4MC/detector/setTargetLength",this);  
  fTargetLengthCmd->SetGuidance("Set LH2 target length  (value and unit)");
  fTargetLengthCmd->SetUnitCategory("Length");

  fBeamlineCmd = new G4UIcmdWithAnInteger("/sFFG4MC/detector/setBeamline",this);
  fBeamlineCmd->SetGuidance("Set an integer flag for whether to include scattering chamber and beamline (0 or 1)");

  fUpdateCmd = new G4UIcommand("/sFFG4MC/detector/update",this);
  fUpdateCmd->SetGuidance("Update the detector geometry with changed values.");
  fUpdateCmd->SetGuidance("Must be run before beamOn if detector has been changed.");  
}

//---------------------------------------------------------------------------

DetectorMessenger::~DetectorMessenger()
{;}

//---------------------------------------------------------------------------

void DetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{ 
  if(command == fNPSAngleCmd)
    { fDetector->SetNPSAngle(fNPSAngleCmd->GetNewDoubleValue(newValue)); }

  if(command == fNPSDistanceCmd)
    { fDetector->SetNPSDistance(fNPSDistanceCmd->GetNewDoubleValue(newValue)); }

  if(command == fNPSShieldThicknessCmd)
    { fDetector->SetNPSShieldThickness(fNPSShieldThicknessCmd->GetNewDoubleValue(newValue)); }

  if(command == fHCALAngleCmd)
    { fDetector->SetHCALAngle(fHCALAngleCmd->GetNewDoubleValue(newValue)); }

  if(command == fHCALDistanceCmd)
    { fDetector->SetHCALDistance(fHCALDistanceCmd->GetNewDoubleValue(newValue)); }

  if(command == fHCALShieldThicknessCmd)
    { fDetector->SetHCALShieldThickness(fHCALShieldThicknessCmd->GetNewDoubleValue(newValue)); }

  if(command == fWindowThicknessCmd)
    { fDetector->SetWindowThickness(fWindowThicknessCmd->GetNewDoubleValue(newValue)); }

  if(command == fTargetLengthCmd)
    { fDetector->SetTargetLength(fTargetLengthCmd->GetNewDoubleValue(newValue)); }

  if(command == fBeamlineCmd)
    { fDetector->SetBeamlineOn(fBeamlineCmd->GetNewIntValue(newValue)); }

  if(command == fUpdateCmd)
    fDetector->UpdateGeometry();

}

//---------------------------------------------------------------------------
