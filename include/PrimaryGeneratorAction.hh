#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleTable.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TLorentzVector.h"

//---------------------------------------------------------------------------

class PrimaryGeneratorMessenger;
class G4ParticleGun;
class G4GeneralParticleSource;
class G4Event;

enum { EPGA_BEAM, EPGA_ROOT };

//---------------------------------------------------------------------------

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  PrimaryGeneratorAction( );    
  ~PrimaryGeneratorAction();

  void GeneratePrimaries (G4Event*);
  void SetUpROOTInput    (TString filename);

  void SetMode(G4int mo)          { fMode       = mo;  }
  void SetNEvents(G4int nev)      { fNevents    = nev; }

  G4int         GetMode()                       { return fMode; }
  G4int         GetNEvents()                    { return fNevents; }

  G4double      GetWeight()                     { return fWeight; }
  G4int         GetFlag()                       { return fFlag; }
  G4int         GetNPrimaryParticles()          { return fNPrimParticles; }
  G4ThreeVector GetVertex( G4int i )            { return G4ThreeVector(fVx[i],  fVy[i],  fVz[i]);  }
  G4ThreeVector GetDirection( G4int i )         { return G4ThreeVector(fPx[i], fPy[i], fPz[i]); }
  G4double      GetEnergy( G4int i )            { return fE[i]; }
  G4ParticleDefinition* GetPrimPDef(  G4int i ) { return fPDefinition[i]; }

private:

  PrimaryGeneratorMessenger* fGunMessenger;
  G4ParticleGun*             fParticleGun;    
  G4GeneralParticleSource*   fParticleSource;
  G4ParticleTable*           fParticleTable;

  G4int                      fMode;   
  G4int                      fNPrimParticles;   
  G4double                   fWeight;
  G4int                      fFlag;
  G4int                      fNevents;
  G4int                      fEvent;
  
  static const G4int         fMaxprim = 50;
  
  G4ParticleDefinition*      fPDefinition[fMaxprim];

  TFile*                     fGenFile;    
  TTree*                     fGenTree;  
  Int_t                      fNGenBranches;   

  Float_t                   fWin;
  Int_t                     fPDG[fMaxprim];
  
  Float_t                   fVx[fMaxprim];
  Float_t                   fVy[fMaxprim];
  Float_t                   fVz[fMaxprim];
  Float_t                   fPx[fMaxprim];
  Float_t                   fPy[fMaxprim];
  Float_t                   fPz[fMaxprim];
  Float_t                   fE[fMaxprim];

};
#endif

//---------------------------------------------------------------------------


