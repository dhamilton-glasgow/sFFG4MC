#include "PrimaryGeneratorMessenger.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"

#include "Randomize.hh"

//---------------------------------------------------------------------------

PrimaryGeneratorMessenger::PrimaryGeneratorMessenger( PrimaryGeneratorAction* Gun )
:fAction(Gun)
{
  fGunDir = new G4UIdirectory("/sFFG4MC/generator/");
  fGunDir->SetGuidance("PrimaryGenerator control");

  fSetModeCmd = new G4UIcmdWithAnInteger("/sFFG4MC/generator/Mode",this);
  fSetModeCmd->SetGuidance("Set the mode of the generator (0 for BEAM or 1 for ROOT)");

  fSetNEventsCmd = new G4UIcmdWithAnInteger("/sFFG4MC/generator/Nevents",this);
  fSetNEventsCmd->SetGuidance("Set the number of primary events to be generated");

  fSetInputCmd = new G4UIcmdWithAString("/sFFG4MC/generator/InputFile",this);
  fSetInputCmd->SetGuidance("Set the full name and path of the file with the input ROOT ntuple (in ROOT mode)"); 
}

//---------------------------------------------------------------------------

PrimaryGeneratorMessenger::~PrimaryGeneratorMessenger()
{;}

//---------------------------------------------------------------------------


void PrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{ 
  if( command == fSetModeCmd )
     { fAction->SetMode(fSetModeCmd->GetNewIntValue(newValue));}

  if( command == fSetNEventsCmd )
     { fAction->SetNEvents(fSetNEventsCmd->GetNewIntValue(newValue));}

  if( command == fSetInputCmd )
    { fAction->SetUpROOTInput(static_cast<TString>(newValue));}
   
}

//---------------------------------------------------------------------------


