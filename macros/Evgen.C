#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TMath.h"
#include "TSystem.h"
#include "TDatabasePDG.h"
#include "TPDGCode.h"
#include "TH2.h"
#include "TBenchmark.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TStyle.h"

#include <iostream>
using namespace std;

// ----------------------------------------------------------------------------

void    InitOutput( TString );
Int_t   SelectReaction();
Float_t GenerateReaction();
Int_t   DecayReaction();

// ----------------------------------------------------------------------------

TRandom3*       fRand;
TDatabasePDG*   fDBpdg;

TFile*          fROOTFile;
TTree*          fROOTTree;

const Int_t     fMaxparticles = 10;
Int_t           fNparticles;  
Float_t         fWeight;
Int_t           fReactFlag;  
Float_t         fVx[fMaxparticles];
Float_t         fVy[fMaxparticles];
Float_t         fVz[fMaxparticles];
Float_t         fPx[fMaxparticles];
Float_t         fPy[fMaxparticles];
Float_t         fPz[fMaxparticles];
Float_t         fE[fMaxparticles];
Int_t           fPDG[fMaxparticles];

TVector3*       fVertex[fMaxparticles];
TLorentzVector* fP4Lab[fMaxparticles];

// ----------------------------------------------------------------------------

void Evgen( TString outrootfile = "testgen.root", ULong64_t nev = 1000 ) 
{

  fRand = new TRandom3(-1);

  InitOutput( outrootfile );

  for( int i=0 ; i < fMaxparticles ; i++ ) {
    fVertex[i] = new TVector3( 0., 0., 0. );
    fP4Lab[i]  = new TLorentzVector( 0., 0., 0., 1. );
    fPDG[i]    = -9999999;
  }

  fDBpdg = new TDatabasePDG();
  TString pdgtable = gSystem->Getenv("ROOTSYS");
  pdgtable.Append("/etc/pdg_table.txt");
  fDBpdg->ReadPDGTable(pdgtable);

  // ----------------------------------------------------------------------------

  for( ULong64_t i = 0; i < nev; i++ ) {

    if( nev%1000 == 0 ) {
     printf("Event %8lld\r", nev);
     fflush(stdout);
    }

    fReactFlag  = SelectReaction();
    fWeight     = GenerateReaction();
    fNparticles = DecayReaction();

    if( fNparticles < 1 ) continue; 

    for( int i=0 ; i < fNparticles ; i++ ) {

      fVx[i]  = (Float_t)fVertex[i]->X();
      fVy[i]  = (Float_t)fVertex[i]->Y();
      fVz[i]  = (Float_t)fVertex[i]->Z();

      fPx[i]  = (Float_t)fP4Lab[i]->Px();
      fPy[i]  = (Float_t)fP4Lab[i]->Py();
      fPz[i]  = (Float_t)fP4Lab[i]->Pz();
      fE[i]   = (Float_t)fP4Lab[i]->E();
    }
    
    fROOTTree->Fill();
  }

  fROOTTree->Write();
  fROOTFile->Close();

}

// ----------------------------------------------------------------------------

void InitOutput( TString outname )
{
  fROOTFile = new TFile(outname,"RECREATE");
  fROOTTree = new TTree("TGen", "Generator tree");
  fROOTTree->SetAutoSave();

  fROOTTree->Branch("weight",     &fWeight,     "weight/F");
  fROOTTree->Branch("flag",       &fReactFlag,  "flag/I");
  fROOTTree->Branch("Nparticles", &fNparticles, "Nparticles/I");

  fROOTTree->Branch("vx",  fVx,  "vx[Nparticles]/F");
  fROOTTree->Branch("vy",  fVy,  "vy[Nparticles]/F");
  fROOTTree->Branch("vz",  fVz,  "vz[Nparticles]/F");
  fROOTTree->Branch("px",  fPx,  "px[Nparticles]/F");
  fROOTTree->Branch("py",  fPy,  "py[Nparticles]/F");
  fROOTTree->Branch("pz",  fPz,  "pz[Nparticles]/F");
  fROOTTree->Branch("E",   fE,   "E[Nparticles]/F");
  fROOTTree->Branch("pdg", fPDG, "pdg[Nparticles]/I");
}

// ----------------------------------------------------------------------------

Int_t SelectReaction()
{
  return 0;
}
  
// ----------------------------------------------------------------------------

Float_t GenerateReaction( )
{
  return 0.0;
}

// ----------------------------------------------------------------------------

Int_t DecayReaction()
{
  return 1;
}

// ----------------------------------------------------------------------------
