#include "SteppingAction.hh"
#include "DetectorConstruction.hh"
#include "VirtualDetectorSD.hh"

#include "G4Track.hh"
#include "G4Gamma.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"

using namespace CLHEP;

//---------------------------------------------------------------------------

SteppingAction::SteppingAction(DetectorConstruction* det)
  :detector(det)
{;}

//---------------------------------------------------------------------------

SteppingAction::~SteppingAction()
{;}

//---------------------------------------------------------------------------

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  G4StepPoint* prePoint = aStep->GetPreStepPoint();
  G4StepPoint* endPoint = aStep->GetPostStepPoint();

  VirtualDetectorSD* virtualDetectorSD = (VirtualDetectorSD*)detector->GetVirtualDetectorSD();

  // Trigger sensitive virtual detectors 
  // (note the for loop starts at 1 as 0 is the default copy number for non-SD volumes)
  // copy number 9999 is blacklisted -- do not trigger virtualSD if the particle originates in 9999
  if( prePoint->GetTouchableHandle()->GetVolume()->GetCopyNo() != 9999 ) {
    for(int i=1; i<=detector->GetNoSD(); i++) {
      
      if ( prePoint->GetTouchableHandle()->GetVolume() != detector->GetDetVol(i) &&
	   endPoint->GetTouchableHandle()->GetVolume() == detector->GetDetVol(i) ) {
	
	if(virtualDetectorSD) 
	  virtualDetectorSD->ProcessHits_constStep(aStep,NULL);
      }
    }
  }

}

//---------------------------------------------------------------------------


