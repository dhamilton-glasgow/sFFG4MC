#include "G4GeneralParticleSource.hh"
#include "G4Event.hh"
#include "G4String.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"

#include "TObjArray.h"
#include "TBranch.h"
#include "TString.h"

#include "PrimaryGeneratorAction.hh"
#include "PrimaryGeneratorMessenger.hh"

using namespace CLHEP;

//---------------------------------------------------------------------------

PrimaryGeneratorAction::PrimaryGeneratorAction()
{
  fMode               = EPGA_BEAM;
  fNevents            = -1;
  fEvent              = 0;
  fWeight             = 0;
  fFlag               = 0;

  fGenFile            = NULL;
  fGenTree            = NULL;
  fNGenBranches       = 0;

  fParticleSource    = new G4GeneralParticleSource();
  fParticleGun       = new G4ParticleGun(1);
  fParticleTable     = G4ParticleTable::GetParticleTable();

  fGunMessenger      = new PrimaryGeneratorMessenger(this);

}

//---------------------------------------------------------------------------

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  if( fGenFile ) {
    fGenFile->Close();
    delete fGenFile;
  }
  delete fParticleGun;
  delete fParticleSource;
  delete fGunMessenger;
}

//---------------------------------------------------------------------------

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  
  switch(fMode) {
    
  case EPGA_BEAM:
    
    fParticleSource->GeneratePrimaryVertex(anEvent);

    fFlag           = 0;
    fNPrimParticles = 1;
    fWeight         = (1./fNevents)*1e-6/1.602e-19;  // event weight per uA beam current
    
    fVx[0]          = fParticleSource->GetParticlePosition().getX();
    fVy[0]          = fParticleSource->GetParticlePosition().getY();
    fVz[0]          = fParticleSource->GetParticlePosition().getZ();
    fPx[0]          = fParticleSource->GetParticleMomentumDirection().getX();
    fPy[0]          = fParticleSource->GetParticleMomentumDirection().getY();
    fPz[0]          = fParticleSource->GetParticleMomentumDirection().getZ();
    fE[0]           = fParticleSource->GetParticleEnergy();
    fPDefinition[0] = fParticleSource->GetParticleDefinition();
    
    break;
    
  case EPGA_ROOT:
    
    if(fGenTree) {

      fGenTree->GetEvent(fEvent++);

      for( Int_t j = 0; j < fNPrimParticles; j++ ) {

	fPDefinition[j] = fParticleTable->FindParticle( fPDG[j] );
	
	fParticleGun->SetParticlePosition          ( G4ThreeVector( fVx[j]*cm, fVy[j]*cm , fVz[j]*cm) );
	fParticleGun->SetParticleDefinition        ( fPDefinition[j] );
	fParticleGun->SetParticleMomentumDirection ( G4ThreeVector(fPx[j], fPy[j], fPz[j]).unit() );
	fParticleGun->SetParticleEnergy            ( fE[j] );
	fParticleGun->GeneratePrimaryVertex(anEvent);
	
	fVx[j]          = fParticleGun->GetParticlePosition().getX();
	fVy[j]          = fParticleGun->GetParticlePosition().getY();
	fVz[j]          = fParticleGun->GetParticlePosition().getZ();
	fPx[j]          = fParticleGun->GetParticleMomentumDirection().getX();
	fPy[j]          = fParticleGun->GetParticleMomentumDirection().getY();
	fPz[j]          = fParticleGun->GetParticleMomentumDirection().getZ();
	fE[j]           = fParticleGun->GetParticleEnergy();
	fPDefinition[j] = fParticleGun->GetParticleDefinition();
	
      }
    }
    break;

  default:
    G4cout << "Unknown mode given to PrimiaryGeneratorAction (0 for BEAM, 1 for ROOT)" << G4endl;
  }			       
}


//---------------------------------------------------------------------------

void PrimaryGeneratorAction::SetUpROOTInput(TString filename)
{
  
  fMode = EPGA_ROOT;
  
  fGenFile = new TFile(filename);
  if(!fGenFile)
    G4cout << "PrimaryGeneratorAction::SetUpRootInput(TString filename) - Didn't find filename" << G4endl;
  
  fGenTree = dynamic_cast<TTree*>(fGenFile->Get("TGen"));
  if(!fGenTree)
    G4cout << "PrimaryGeneratorAction::SetUpRootInput(TString filename) - Didn't find tree TGen" << G4endl;
    
  fNGenBranches       = fGenTree->GetNbranches();
  TObjArray* objarray = fGenTree->GetListOfBranches();

  for( Int_t i = 0; i < fNGenBranches; i++ ) {
    
    TBranch *branch = dynamic_cast<TBranch*>     (objarray->At(i));
    TString  bname  = TString( const_cast<char*> (branch->GetName()) );

    if( bname == "Nparticles" ) branch->SetAddress( &fNPrimParticles );
    if( bname == "weight" )     branch->SetAddress( &fWeight  );
    if( bname == "flag" )       branch->SetAddress( &fFlag );

    if( fNPrimParticles < fMaxprim ) {
      if( bname == "vx"     ) branch->SetAddress( fVx );
      if( bname == "vy"     ) branch->SetAddress( fVy );
      if( bname == "vz"     ) branch->SetAddress( fVz );
      if( bname == "px"     ) branch->SetAddress( fPx );
      if( bname == "py"     ) branch->SetAddress( fPy );
      if( bname == "pz"     ) branch->SetAddress( fPz );
      if( bname == "E"      ) branch->SetAddress( fE  );
      if( bname == "pdg"    ) branch->SetAddress( fPDG );
    }
  }

  G4int nev = fGenTree->GetEntries();
  if( fNevents == -1 || fNevents > nev )
    fNevents = nev;
  
}

//---------------------------------------------------------------------------
