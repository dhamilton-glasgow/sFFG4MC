#include "OutputManager.hh"
#include "PrimaryGeneratorAction.hh"
#include "OutputMessenger.hh"

#include "G4UnitsTable.hh"
#include "G4RunManager.hh"
#include "G4Point3D.hh"
#include "G4Transform3D.hh"
#include "G4SystemOfUnits.hh"

#include "TROOT.h"
#include "TApplication.h"
#include "TFile.h"
#include "TSystem.h"
#include "TTree.h"
#include "TString.h"
#include "TMath.h"

//---------------------------------------------------------------------------

OutputManager::OutputManager()
{
  ZeroArray();

  fOutFileName = TString("output/out_default.root");

  fOutMessenger = new OutputMessenger(this);

}

//---------------------------------------------------------------------------

OutputManager::~OutputManager()
{
  if(fROOTfile) {
    fROOTtree->Write();
    fROOTfile->Close();
  }
}

//---------------------------------------------------------------------------

void OutputManager::InitOutput()
{ 

  fROOTfile = new TFile(fOutFileName,"RECREATE","fROOTfile",1);
  fROOTtree = new TTree("TOut","Output Tree");
  fROOTtree->SetAutoSave();

  // Set Event Branches
  fROOTtree->Branch("Event_weight",   &fEvent_weight, "Event_weight/F");  
  fROOTtree->Branch("Event_flag",     &fEvent_flag,   "Event_flag/I");  
  fROOTtree->Branch("Event_num",      &fEvent_num,    "Event_num/I");  

  // Set Primary Branches
  fROOTtree->Branch("Primary_Nhits",  &fPrimary_Nhits, "Primary_Nhits/I");  
  fROOTtree->Branch("Primary_pdg",    fPrimary_pdg,    "Primary_pdg[Primary_Nhits]/I");
  fROOTtree->Branch("Primary_E",      fPrimary_E,      "Primary_E[Primary_Nhits]/F");
  fROOTtree->Branch("Primary_x",      fPrimary_xpos,   "Primary_x[Primary_Nhits]/F");
  fROOTtree->Branch("Primary_y",      fPrimary_ypos,   "Primary_y[Primary_Nhits]/F");
  fROOTtree->Branch("Primary_z",      fPrimary_zpos,   "Primary_z[Primary_Nhits]/F");
  fROOTtree->Branch("Primary_px",     fPrimary_px,     "Primary_px[Primary_Nhits]/F");
  fROOTtree->Branch("Primary_py",     fPrimary_py,     "Primary_py[Primary_Nhits]/F");
  fROOTtree->Branch("Primary_pz",     fPrimary_pz,     "Primary_pz[Primary_Nhits]/F");

  // Set VirtualDetector Hit Branches
  fROOTtree->Branch("Virtual_Nhits", &fVirtual_Nhits, "Virtual_Nhits/I");  
  fROOTtree->Branch("Virtual_pdg",   fVirtual_pdg,    "Virtual_pdg[Virtual_Nhits]/I");
  fROOTtree->Branch("Virtual_det",   fVirtual_det,    "Virtual_det[Virtual_Nhits]/I");
  fROOTtree->Branch("Virtual_mod",   fVirtual_mod,    "Virtual_mod[Virtual_Nhits]/I");
  fROOTtree->Branch("Virtual_row",   fVirtual_row,    "Virtual_row[Virtual_Nhits]/I");
  fROOTtree->Branch("Virtual_col",   fVirtual_col,    "Virtual_col[Virtual_Nhits]/I");
  fROOTtree->Branch("Virtual_tid",   fVirtual_tid,    "Virtual_tid[Virtual_Nhits]/I");
  fROOTtree->Branch("Virtual_pid",   fVirtual_pid,    "Virtual_pid[Virtual_Nhits]/I");
  fROOTtree->Branch("Virtual_E",     fVirtual_E,      "Virtual_E[Virtual_Nhits]/F" );
  fROOTtree->Branch("Virtual_t",     fVirtual_t,      "Virtual_t[Virtual_Nhits]/F" );
  fROOTtree->Branch("Virtual_x",     fVirtual_xpos,   "Virtual_x[Virtual_Nhits]/F"  );
  fROOTtree->Branch("Virtual_y",     fVirtual_ypos,   "Virtual_y[Virtual_Nhits]/F"  );
  fROOTtree->Branch("Virtual_z",     fVirtual_zpos,   "Virtual_z[Virtual_Nhits]/F"  );
  fROOTtree->Branch("Virtual_px",    fVirtual_px,     "Virtual_px[Virtual_Nhits]/F"  );
  fROOTtree->Branch("Virtual_py",    fVirtual_py,     "Virtual_py[Virtual_Nhits]/F"  );
  fROOTtree->Branch("Virtual_pz",    fVirtual_pz,     "Virtual_pz[Virtual_Nhits]/F"  );
  fROOTtree->Branch("Virtual_vx",    fVirtual_vx,     "Virtual_vx[Virtual_Nhits]/F"  );
  fROOTtree->Branch("Virtual_vy",    fVirtual_vy,     "Virtual_vy[Virtual_Nhits]/F"  );
  fROOTtree->Branch("Virtual_vz",    fVirtual_vz,     "Virtual_vz[Virtual_Nhits]/F"  );

  // Set RealDetector Hit Branches
  fROOTtree->Branch("Real_Nhits", &fReal_Nhits, "Real_Nhits/I");  
  fROOTtree->Branch("Real_det",   fReal_det,    "Real_det[Real_Nhits]/I");
  fROOTtree->Branch("Real_mod",   fReal_mod,    "Real_mod[Real_Nhits]/I");
  fROOTtree->Branch("Real_row",   fReal_row,    "Real_row[Real_Nhits]/I");
  fROOTtree->Branch("Real_col",   fReal_col,    "Real_col[Real_Nhits]/I");
  fROOTtree->Branch("Real_edep",  fReal_Edep,   "Real_edep[Real_Nhits]/F" );
  fROOTtree->Branch("Real_t",     fReal_t,      "Real_t[Real_Nhits]/F" );
  fROOTtree->Branch("Real_x",     fReal_xpos,   "Real_x[Real_Nhits]/F"  );
  fROOTtree->Branch("Real_y",     fReal_ypos,   "Real_y[Real_Nhits]/F"  );
  fROOTtree->Branch("Real_z",     fReal_zpos,   "Real_z[Real_Nhits]/F"  );
}

//---------------------------------------------------------------------------

void OutputManager::ZeroArray()
{
  G4ThreeVector zero(0.,0.,0.);

  fPrimary_PDef    = NULL;
  fPrimary_dir     = (zero);
  fPrimary_vtx     = (zero);
  fPrimary_energy  = 0;

  fPrimary_Nhits   = 0;

  for ( Int_t i = 0; i < fMaxprim; i++ ) {
    fPrimary_pdg[i]    = -INT_MAX;
    fPrimary_E[i]      = -INT_MAX;
    fPrimary_xpos[i]   = -INT_MAX;
    fPrimary_ypos[i]   = -INT_MAX;
    fPrimary_zpos[i]   = -INT_MAX;
    fPrimary_px[i]     = -INT_MAX;
    fPrimary_py[i]     = -INT_MAX;
    fPrimary_pz[i]     = -INT_MAX;
  }

  fVirtual_pdef    = NULL;
  fVirtual_p3      = (zero);
  fVirtual_pospre  = (zero);
  fVirtual_pospost = (zero);
  fVirtual_time    = 0;
  fVirtual_detid   = 0;
  fVirtual_Pid     = 0;
  fVirtual_Tid     = 0;
  fVirtual_energy  = 0;

  fReal_pospre     = (zero);
  fReal_pospost    = (zero);
  fReal_time       = 0;
  fReal_detid      = 0;
  fReal_edep       = 0;

  fVirtual_Nhits   = 0;
  fReal_Nhits      = 0;

  for ( Int_t i = 0; i < fMaxhits; i++ ) {
    fVirtual_pdg[i]     = -INT_MAX;
    fVirtual_tid[i]     = -INT_MAX;
    fVirtual_pid[i]     = -INT_MAX; 
    fVirtual_E[i]       = -INT_MAX;
    fVirtual_t[i]       = -INT_MAX;
    fVirtual_xpos[i]    = -INT_MAX;
    fVirtual_ypos[i]    = -INT_MAX;
    fVirtual_zpos[i]    = -INT_MAX;
    fVirtual_px[i]      = -INT_MAX;
    fVirtual_py[i]      = -INT_MAX;
    fVirtual_pz[i]      = -INT_MAX;
    fVirtual_vx[i]      = -INT_MAX;
    fVirtual_vy[i]      = -INT_MAX;
    fVirtual_vz[i]      = -INT_MAX;
    fVirtual_det[i]     = -INT_MAX;
    fVirtual_mod[i]     = -INT_MAX;
    fVirtual_row[i]     = -INT_MAX;
    fVirtual_col[i]     = -INT_MAX;

    fReal_Edep[i]       = -INT_MAX;
    fReal_t[i]          = -INT_MAX;
    fReal_xpos[i]       = -INT_MAX;
    fReal_ypos[i]       = -INT_MAX;
    fReal_zpos[i]       = -INT_MAX;
    fReal_det[i]        = -INT_MAX;
    fReal_mod[i]        = -INT_MAX;
    fReal_row[i]        = -INT_MAX;
    fReal_col[i]        = -INT_MAX;
  }

}

//---------------------------------------------------------------------------

void OutputManager::FillVirtualArray( Int_t hitn ) 
{

  if( hitn < fMaxhits ) {
    
    fVirtual_pdg[hitn]    = (Int_t)fVirtual_pdef->GetPDGEncoding();
    fVirtual_tid[hitn]    = (Int_t)fVirtual_Tid;
    fVirtual_pid[hitn]    = (Int_t)fVirtual_Pid;
    fVirtual_t[hitn]      = (Float_t)fVirtual_time *ns;                                   
    fVirtual_E[hitn]      = (Float_t)fVirtual_energy *MeV;                                

    Float_t xpost = (Float_t)fVirtual_pospost.getX() /cm;                                                          
    Float_t ypost = (Float_t)fVirtual_pospost.getY() /cm;                                                          
    Float_t zpost = (Float_t)fVirtual_pospost.getZ() /cm;                                                          

    fVirtual_xpos[hitn]   = xpost;
    fVirtual_ypos[hitn]   = ypost;
    fVirtual_zpos[hitn]   = zpost;
    fVirtual_px[hitn]     = (Float_t)fVirtual_p3.getX();                             
    fVirtual_py[hitn]     = (Float_t)fVirtual_p3.getY();                             
    fVirtual_pz[hitn]     = (Float_t)fVirtual_p3.getZ();
    fVirtual_vx[hitn]     = (Float_t)fVirtual_vtx.getX()/cm;                             
    fVirtual_vy[hitn]     = (Float_t)fVirtual_vtx.getY()/cm;                                                          
    fVirtual_vz[hitn]     = (Float_t)fVirtual_vtx.getZ()/cm;                             
    
    if( fVirtual_detid >= 1 && fVirtual_detid <= fNearm ) { 
      fVirtual_det[hitn] = 0;                               // NPS
      fVirtual_mod[hitn] = (fVirtual_detid-1)/(fNearm/6);          // Module
      fVirtual_row[hitn] = (fVirtual_detid-1)%(fNearm/6)/fNearmCol;       // Row
      fVirtual_col[hitn] = (fVirtual_detid-1)%(fNearm/6)/fNearmCol;       // Column
    }
    else if( fVirtual_detid >= (fNearm+1) && fVirtual_detid <= (fNearm+fNhodo) ) { 
      fVirtual_det[hitn] = 1;                               // Hodoscope
      fVirtual_mod[hitn] = (fVirtual_detid-(fNearm+1))/(fNhodo/6);       // Module
      fVirtual_row[hitn] = (fVirtual_detid-(fNearm+1))%(fNhodo/6)/fNhodoCol;    // Row
      fVirtual_col[hitn] = (fVirtual_detid-(fNearm+1))%(fNhodo/6)%fNhodoCol;    // Column
    }
    else if( fVirtual_detid >= (fNearm+fNhodo+1) && fVirtual_detid <= (fNearm+fNhodo+fNharm+1) ) { 
      fVirtual_det[hitn] = 2;                               // HCAL
      fVirtual_mod[hitn] = (fVirtual_detid-(fNearm+fNhodo+1))/(fNharm/6);        // Module
      fVirtual_row[hitn] = (fVirtual_detid-(fNearm+fNhodo+1))%(fNharm/6)/fNharmCol;     // Row
      fVirtual_col[hitn] = (fVirtual_detid-(fNearm+fNhodo+1))%(fNharm/6)%fNharmCol;     // Column
    }

    fVirtual_Nhits++;
  }
}

//---------------------------------------------------------------------------

void OutputManager::FillRealArray( G4int hitn ) 
{

  if( hitn < fMaxhits ) {
    
    fReal_Edep[hitn]   = (Float_t)fReal_edep *MeV;                                   
    fReal_t[hitn]      = (Float_t)fReal_time *ns;                                   

    Float_t xpre  = (Float_t)fReal_pospre.getX() /cm;                             
    Float_t ypre  = (Float_t)fReal_pospre.getY() /cm;                                                          
    Float_t zpre  = (Float_t)fReal_pospre.getZ() /cm;                                                          
    Float_t xpost = (Float_t)fReal_pospost.getX() /cm;                                                          
    Float_t ypost = (Float_t)fReal_pospost.getY() /cm;                                                          
    Float_t zpost = (Float_t)fReal_pospost.getZ() /cm;                                                          

    fReal_xpos[hitn]   = (xpre + xpost)/2.;
    fReal_ypos[hitn]   = (ypre + ypost)/2.;
    fReal_zpos[hitn]   = (zpre + zpost)/2.;

    if( fReal_detid >= 1 && fReal_detid <= fNearm ) { 
      fReal_det[hitn] = 0;                               // NPS
      fReal_mod[hitn] = (fReal_detid-1)/(fNearm/6);          // Module
      fReal_row[hitn] = (fReal_detid-1)%(fNearm/6)/fNearmCol;       // Row
      fReal_col[hitn] = (fReal_detid-1)%(fNearm/6)/fNearmCol;       // Column
    }
    else if( fReal_detid >= (fNearm+1) && fReal_detid <= (fNearm+fNhodo) ) { 
      fReal_det[hitn] = 1;                               // Hodoscope
      fReal_mod[hitn] = (fReal_detid-(fNearm+1))/(fNhodo/6);       // Module
      fReal_row[hitn] = (fReal_detid-(fNearm+1))%(fNhodo/6)/fNhodoCol;    // Row
      fReal_col[hitn] = (fReal_detid-(fNearm+1))%(fNhodo/6)%fNhodoCol;    // Column
    }
    else if( fReal_detid >= (fNearm+fNhodo+1) && fReal_detid <= (fNearm+fNhodo+fNharm+1) ) { 
      fReal_det[hitn] = 2;                               // HCAL
      fReal_mod[hitn] = (fReal_detid-(fNearm+fNhodo+1))/(fNharm/6);        // Module
      fReal_row[hitn] = (fReal_detid-(fNearm+fNhodo+1))%(fNharm/6)/fNharmCol;     // Row
      fReal_col[hitn] = (fReal_detid-(fNearm+fNhodo+1))%(fNharm/6)%fNharmCol;     // Column
    }
    
    fReal_Nhits++;
  }

}

//---------------------------------------------------------------------------

void OutputManager::FillPrimaryArray( G4int primn ) 
{

  if( primn < fMaxprim ) {
    
    fPrimary_pdg[primn]  = (Int_t)fPrimary_PDef->GetPDGEncoding();
    fPrimary_E[primn]    = (Float_t)fPrimary_energy *MeV; 
    fPrimary_xpos[primn] = (Float_t)fPrimary_vtx.getX() /cm;                             
    fPrimary_ypos[primn] = (Float_t)fPrimary_vtx.getY() /cm;                             
    fPrimary_zpos[primn] = (Float_t)fPrimary_vtx.getZ() /cm;                             
    fPrimary_px[primn]   = (Float_t)fPrimary_dir.getX() /cm;                             
    fPrimary_py[primn]   = (Float_t)fPrimary_dir.getY() /cm;                             
    fPrimary_pz[primn]   = (Float_t)fPrimary_dir.getZ() /cm;                             
    
    fPrimary_Nhits++;
    
  }

}

//---------------------------------------------------------------------------

void OutputManager::FillTree(G4double weight, G4int flag, G4int num)
{

  fEvent_weight = (Float_t)weight;
  fEvent_flag   = (Int_t)flag;
  fEvent_num    = (Int_t)num;

  fROOTtree->Fill();
}

//---------------------------------------------------------------------------
