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

void InitOutput( TString );
void GenerateReaction();

// ----------------------------------------------------------------------------

const Int_t     fMaxparticles = 10;
TVector3*       fVertex[fMaxparticles];
TLorentzVector* fP4Lab[fMaxparticles];

TRandom3*     fRand;
TDatabasePDG* fDBpdg;
TFile*        fROOTFile;
TTree*        fROOTTree;
Int_t         fNparticles;  
Float_t       fWeight;
Int_t         fReactFlag;  
Float_t       fVx[fMaxparticles];
Float_t       fVy[fMaxparticles];
Float_t       fVz[fMaxparticles];
Float_t       fPx[fMaxparticles];
Float_t       fPy[fMaxparticles];
Float_t       fPz[fMaxparticles];
Float_t       fE[fMaxparticles];
Int_t         fPDG[fMaxparticles];

enum { kEP, kNReact }; 

Float_t fRasterSize = 0.2;   // cm
Float_t fTarLength  = 10.;   // cm
Float_t fBeamE      = 6.6;   // GeV
Float_t fThetaMin   = 13.5 * TMath::DegToRad(); 
Float_t fThetaMax   = 17.5 * TMath::DegToRad(); 
Int_t   fNgen       = 10000;

// ----------------------------------------------------------------------------

void Evgen( TString outrootfile = "testgen.root" ) 
{
  
  fRand = new TRandom3(-1);

  InitOutput( outrootfile );

  for( int i=0 ; i < fMaxparticles ; i++ ) {
    fVertex[i] = new TVector3( 0., 0., -30. );
    fP4Lab[i]  = new TLorentzVector( 0., 0., fBeamE, fBeamE ); 
    fPDG[i]    = 11;
  }
  
  fDBpdg = new TDatabasePDG();
  TString pdgtable = gSystem->Getenv("ROOTSYS");
  pdgtable.Append("/etc/pdg_table.txt");
  fDBpdg->ReadPDGTable(pdgtable);
  
  // ----------------------------------------------------------------------------

  for( Int_t i = 0; i < fNgen; i++ ) {

    if( i%1000 == 0 ) {
     printf("Event %8d\r", i);
     fflush(stdout);
    }

    fReactFlag = fRand->Integer( kNReact );
    GenerateReaction();
    
    if( fNparticles < 1 ) continue; 
    
    for( int i=0 ; i < fNparticles ; i++ ) {
      
      fVx[i]  = (Float_t)fVertex[i]->X();
      fVy[i]  = (Float_t)fVertex[i]->Y();
      fVz[i]  = (Float_t)fVertex[i]->Z();

      fPx[i]  = 1000*(Float_t)fP4Lab[i]->Px();
      fPy[i]  = 1000*(Float_t)fP4Lab[i]->Py();
      fPz[i]  = 1000*(Float_t)fP4Lab[i]->Pz();
      fE[i]   = 1000*(Float_t)fP4Lab[i]->E();
    }
    
    fROOTTree->Fill();
  }

  // ----------------------------------------------------------------------------
  
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

void GenerateReaction()
{

  Float_t xv = fRand->Uniform( -fRasterSize/2, fRasterSize/2 );
  Float_t yv = fRand->Uniform( -fRasterSize/2, fRasterSize/2 );
  Float_t zv = fRand->Uniform( -fTarLength/2, fTarLength/2 );

  Float_t L = (0.0708*6.022e23)*fTarLength*(60.e-6/1.602e-19);

  // Elastic ep scattering from LH2 target
  if( fReactFlag == kEP ) {
    
    fNparticles = 2;
    fWeight = 1.0;
    
    fVertex[0]->SetXYZ( xv, yv, zv );
    fVertex[1]->SetXYZ( xv, yv, zv );
    fPDG[0] = 11;
    fPDG[1] = 2212;

    TLorentzVector beamP4( 0., 0., fBeamE, fBeamE );
    TLorentzVector targP4( 0., 0., 0., fDBpdg->GetParticle(2212)->Mass() );
    
    Float_t mt = fDBpdg->GetParticle(2212)->Mass(); 
    Float_t me = fDBpdg->GetParticle(fPDG[0])->Mass(); 
    Float_t mp = fDBpdg->GetParticle(fPDG[1])->Mass();
    
    Float_t th     = TMath::ACos( fRand->Uniform( TMath::Cos( fThetaMax ), TMath::Cos( fThetaMin ) ));
    Float_t ph     = fRand->Uniform( -TMath::Pi(), TMath::Pi()  );
    Float_t dOmega = TMath::TwoPi()*(cos(fThetaMin)-cos(fThetaMax));
    TVector3 scatP3( TMath::Sin(th)*TMath::Cos(ph), TMath::Sin(th)*TMath::Sin(ph), TMath::Cos(th) );

    Float_t scatE = (fBeamE*targP4.E())/(fBeamE*(1 -TMath::Cos(th)) + targP4.E()-targP4.Vect().Dot(scatP3));
    Float_t scatP = TMath::Sqrt( scatE*scatE - me*me );
    TLorentzVector scatP4( (scatP*scatP3), scatE ); 
    TLorentzVector qP4 = beamP4 - scatP4; 

    TVector3 rest = targP4.BoostVector();
    beamP4.Boost( -rest );
    targP4.Boost( -rest );
    scatP4.Boost( -rest );
    
    TLorentzVector qrestP4  = beamP4 - scatP4;
    TLorentzVector recoilP4 = targP4 + qrestP4;
    scatP4.Boost( rest );
    recoilP4.Boost( rest );
    
    *fP4Lab[0] = scatP4;
    *fP4Lab[1] = recoilP4;

    // ----------------------------------------------------------------------------
    
    Float_t hbarc = 197.327/1000;
    Float_t alpha = 1/137.036;
    Float_t Q2    = -qP4.M2();
    Float_t tau   = Q2/(4*mp*mp);
    Float_t GE    = (1.0-0.24*tau)/(1.0 + 10.98*tau + 12.82*tau*tau + 21.97*tau*tau*tau );
    Float_t GM    = 2.79*(1.0+0.12*tau)/(1.0 + 10.97*tau + 18.86*tau*tau + 6.55*tau*tau*tau );
    
    Float_t dsigMott  = TMath::Power(hbarc*alpha*TMath::Cos(th/2),2) / 
      (4*beamP4.E()*beamP4.E()*TMath::Power(TMath::Sin(th/2),4));        // nb/sr
    Float_t dsigRosen = dsigMott*(scatP4.E()/beamP4.E()) * 
      ( (GE*GE + tau*GM*GM)/(1+tau) + (2*tau*GM*GM)*TMath::Power(TMath::Tan(th/2),2) );

    fWeight = (1./fNgen) * L * dOmega * dsigRosen*1e-9*1e-24;
  }

  // default
  else {
    fWeight = 1.0;
    fNparticles = 1;
  }
}


// ----------------------------------------------------------------------------
