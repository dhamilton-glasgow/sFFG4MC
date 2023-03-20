#ifndef OutputMessenger_h
#define OutputMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

//---------------------------------------------------------------------------

class OutputManager;
class G4UIdirectory;
class G4UIcmdWithAString;

//---------------------------------------------------------------------------

class OutputMessenger: public G4UImessenger
{
  public:
    OutputMessenger(OutputManager*);
   ~OutputMessenger();

    void SetNewValue(G4UIcommand*, G4String);

  private:
    OutputManager*      fOutputManager;
    G4UIdirectory*      fOutputDir;
    G4UIcmdWithAString* fOutFileCmd;
};
#endif

//---------------------------------------------------------------------------

