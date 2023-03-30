#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TMath.h"
#include "TSystem.h"
#include "TDatabasePDG.h"
#include "TPDGCode.h"

#include "cteq/cteqpdf.h"

#include <iostream>
using namespace std;

// ----------------------------------------------------------------------------

cteq_pdf_t *__dis_pdf;

void initcteqpdf();
double dissigma( double, double, double );

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

enum { kEP, kDIS, kNReact }; 

Float_t fRasterSize = 0.2;   // cm
Float_t fTarLength  = 10.;   // cm
Float_t fBeamE      = 6.6;   // GeV
Float_t fThetaMin   = 15.4 * TMath::DegToRad(); 
Float_t fThetaMax   = 15.6 * TMath::DegToRad(); 
Float_t fEMin       = 0.2;
Float_t fEMax       = 5.2;
Int_t   fNgen       = 10000;

// ----------------------------------------------------------------------------

void Evgen( TString outrootfile = "testgen.root" ) 
{

  initcteqpdf();

  fRand = new TRandom3(0);

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

  Float_t L      = (0.0708*6.022e23)*fTarLength*(60.e-6/1.602e-19);
  Float_t dOmega = TMath::TwoPi()*(cos(fThetaMin)-cos(fThetaMax));

  TLorentzVector beamP4( 0., 0., fBeamE, fBeamE );
  TLorentzVector targP4( 0., 0., 0., fDBpdg->GetParticle(2212)->Mass() );
  
  Float_t mt = fDBpdg->GetParticle(2212)->Mass(); 
  Float_t me = fDBpdg->GetParticle(fPDG[0])->Mass(); 
  Float_t th, ph;

  TVector3 scatP3, rest;
  Float_t scatE, scatP, Q2;

  TLorentzVector scatP4, qP4, qrestP4;
 
  Float_t hbarc = 197.327/1000;
  Float_t alpha = 1/137.036;

  switch( fReactFlag ) {
    
  case kEP: { // Elastic ep scattering from LH2 target

    fNparticles = 2;
    fVertex[0]->SetXYZ( xv, yv, zv );
    fVertex[1]->SetXYZ( xv, yv, zv );
    fPDG[0] = 11;
    fPDG[1] = 2212;
    
    // ----------------------------------------------------------------------------
    // Kinematics
    
    Float_t mp = fDBpdg->GetParticle(fPDG[1])->Mass();
    
    th = TMath::ACos( fRand->Uniform( TMath::Cos( fThetaMax ), TMath::Cos( fThetaMin ) ));
    ph = fRand->Uniform( -TMath::Pi(), TMath::Pi()  );
    scatP3.SetXYZ( TMath::Sin(th)*TMath::Cos(ph), TMath::Sin(th)*TMath::Sin(ph), TMath::Cos(th) );

    scatE = (fBeamE*targP4.E())/(fBeamE*(1 -TMath::Cos(th)) + targP4.E()-targP4.Vect().Dot(scatP3));
    scatP = TMath::Sqrt( scatE*scatE - me*me );
    scatP4.SetPxPyPzE( scatP*scatP3.X(), scatP*scatP3.Y(), scatP*scatP3.Z(), scatE ); 
    qP4 = beamP4 - scatP4; 

    rest = targP4.BoostVector();
    beamP4.Boost( -rest );
    targP4.Boost( -rest );
    scatP4.Boost( -rest );
    
    qrestP4  = beamP4 - scatP4;
    TLorentzVector recoilP4 = targP4 + qrestP4;
    
    // ----------------------------------------------------------------------------
    // Cross section
    
    hbarc = 197.327/1000;
    alpha = 1/137.036;
    Q2    = -qP4.M2();

    Float_t tau   = Q2/(4*mp*mp);
    Float_t GE    = (1.0-0.24*tau)/(1.0 + 10.98*tau + 12.82*tau*tau + 21.97*tau*tau*tau );
    Float_t GM    = 2.79*(1.0+0.12*tau)/(1.0 + 10.97*tau + 18.86*tau*tau + 6.55*tau*tau*tau );
    
    Float_t dSigMott  = 1e10*hbarc*hbarc*alpha*alpha/
      (4*beamP4.E()*beamP4.E()*TMath::Power(TMath::Sin(th/2),4))
      *(scatP4.E()/beamP4.E())*TMath::Power(TMath::Cos(th/2),2);  
    Float_t dSigRosen = dSigMott * 
      ( (GE*GE + tau*GM*GM)/(1+tau) + (2*tau*GM*GM)*TMath::Power(TMath::Tan(th/2),2) );

    fWeight = (1.*kNReact/fNgen) * L * dOmega * dSigRosen*1e-36;

    scatP4.Boost( rest );
    recoilP4.Boost( rest );
    *fP4Lab[0] = scatP4;
    *fP4Lab[1] = recoilP4;

    break;
  }
  case kDIS: { // DIS from LH2 target
    
    fNparticles = 1;
    fVertex[0]->SetXYZ( xv, yv, zv );
    fPDG[0] = 11;

    // ----------------------------------------------------------------------------
    // Kinematics

    th = TMath::ACos( fRand->Uniform( TMath::Cos( fThetaMax ), TMath::Cos( fThetaMin ) ));
    ph = fRand->Uniform( -TMath::Pi(), TMath::Pi()  );
    scatP3.SetXYZ( TMath::Sin(th)*TMath::Cos(ph), TMath::Sin(th)*TMath::Sin(ph), TMath::Cos(th) );

    scatE = fRand->Uniform( fEMin, fEMax );
    scatP = TMath::Sqrt( scatE*scatE - me*me );
    scatP4.SetPxPyPzE( scatP*scatP3.X(), scatP*scatP3.Y(), scatP*scatP3.Z(), scatE ); 
 
    qP4 = beamP4 - scatP4; 
    TLorentzVector XP4    = targP4 + qP4;
    TLorentzVector ftotP4 = XP4 + scatP4;
    TLorentzVector itotP4 = beamP4 + targP4;

    rest = targP4.BoostVector();
    beamP4.Boost( -rest );
    targP4.Boost( -rest );
    scatP4.Boost( -rest );

    qrestP4  = beamP4 - scatP4;
    Float_t x = -qrestP4.M2()/(2.0*targP4.Dot( qrestP4 ) );
    
    if( ftotP4.M2() > itotP4.M2() || ftotP4.E() > itotP4.E() || x > 1.0 || x < 0.0 ) {
      scatP4.SetPxPyPzE( 0.0, 0.0, 0.0, 0.0 );
    }
    Float_t threst = TMath::ACos( beamP4.Vect().Unit().Dot(scatP4.Vect().Unit()) );

    // ----------------------------------------------------------------------------
    // Cross section

    Float_t  dSigDIS = dissigma( beamP4.E(), threst, scatP4.E() ); 

    fWeight = (1.*kNReact/fNgen) * L * dOmega * (fEMax - fEMin) * dSigDIS * 1e-33;

    scatP4.Boost( rest );
    *fP4Lab[0] = scatP4;

    break;
  }
  default: {
    fWeight = 1.0;
    fNparticles = 1;
  }
  }
}

// ----------------------------------------------------------------------------
// Use CTEQ6 parameterization

void initcteqpdf(){

  __dis_pdf = cteq_pdf_alloc_id(400); // mode 400 = cteq6.6?

  assert(__dis_pdf);
}

// // ----------------------------------------------------------------------------

double dissigma( double ebeam, double th, double eprime ){

    // Return in nb/(GeV*sr)

    double Q2 = 2.0*eprime*ebeam*(1.0-cos(th));
    double nu = ebeam-eprime;
    double Mp = 0.938;

    double x = Q2/(2.0*Mp*nu);
    double y = nu/ebeam;

    if( ! (0.0 < x && x < 1.0 && 0.0 < y && y < 1.0) ){
	printf("WARNING %s line %d  x = %f, y = %f -> eprime = %f GeV   th = %f deg  ebeam = %f GeV\n", __FILE__,
		__LINE__, x, y, eprime, th*180/3.14159, ebeam );
	return 0.0;;
    }

    double qu = cteq_pdf_evolvepdf(__dis_pdf, 1, x, sqrt(Q2) );
    double qd = cteq_pdf_evolvepdf(__dis_pdf, 2, x, sqrt(Q2) );
    double qubar = cteq_pdf_evolvepdf(__dis_pdf, -1, x, sqrt(Q2) );
    double qdbar = cteq_pdf_evolvepdf(__dis_pdf, -2, x, sqrt(Q2) );

    double quv = qu-qubar;
    double qdv = qd-qdbar;

    double qs = cteq_pdf_evolvepdf(__dis_pdf, 3, x, sqrt(Q2) );

    double F2 = 0.0; 
    double e_u =  2.0/3.0;
    double e_d = -1.0/3.0;

    F2 += x*( e_u*e_u*quv + e_d*e_d*qdv ); 
    F2  += x*(2.0*e_u*e_u*qubar + 2.0*e_d*e_d*(qdbar + qs));

    double F1 = F2/(2.0*x);

    double hbarc = 197.327/1000;
    double alpha = 1/137.036;

    // From PDG
    double ds_dxdy = 4.0*TMath::Pi()*pow(alpha,2)*((1.0-y-pow(x*y*Mp,2)/Q2)*F2+y*y*x*F1)/(x*y*Q2);

    // In GeV^-2
    double ds_dOmega_dE = ds_dxdy*eprime/(2.0*TMath::Pi()*Mp*nu);

    return ds_dOmega_dE*pow(hbarc,2)*1e7; // GeV2 -> nb
 }
