#include <TLorentzVector.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMath.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TSystem.h>

#include <iostream>
#include <fstream>

using namespace std;

const Bool_t ApplyThresh     = true;
const Bool_t ApplyWindow     = true;

const Float_t Earm_threshold = 5;
const Float_t Earm_window    = 1e6;
const Float_t Harm_threshold = 0.5;
const Float_t Harm_window    = 1e6;
const Float_t Hodo_threshold = -1;
const Float_t Hodo_window    = 1e6;

void Analyse_sFFG4MC( Int_t run_no = 14 ) { 
  
  //-----------------------------------------------------------------------------------------------------------------------------

  TChain* TOut = new TChain("TOut");
  
  TOut->Add(Form("out/batch_%d*_beam.root", run_no));

  const Int_t     Maxprim = 50;
  const Int_t     Maxhits = 10000;

  Float_t         Event_weight;
  Int_t           Event_flag;
  Int_t           Primary_Nhits;
  Int_t           Primary_pdg[Maxprim];
  Float_t         Primary_E[Maxprim];
  Float_t         Primary_x[Maxprim];
  Float_t         Primary_y[Maxprim];
  Float_t         Primary_z[Maxprim];
  Float_t         Primary_px[Maxprim];
  Float_t         Primary_py[Maxprim];
  Float_t         Primary_pz[Maxprim];
  Int_t           Virtual_Nhits;
  Int_t           Virtual_pdg[Maxhits];
  Int_t           Virtual_det[Maxhits];
  Int_t           Virtual_mod[Maxhits];
  Int_t           Virtual_row[Maxhits];
  Int_t           Virtual_col[Maxhits];
  Int_t           Virtual_tid[Maxhits];
  Int_t           Virtual_pid[Maxhits];
  Float_t         Virtual_E[Maxhits];
  Float_t         Virtual_t[Maxhits];
  Float_t         Virtual_x[Maxhits];
  Float_t         Virtual_y[Maxhits];
  Float_t         Virtual_z[Maxhits];
  Float_t         Virtual_px[Maxhits];
  Float_t         Virtual_py[Maxhits];
  Float_t         Virtual_pz[Maxhits];
  Float_t         Virtual_vx[Maxhits];
  Float_t         Virtual_vy[Maxhits];
  Float_t         Virtual_vz[Maxhits];
  Int_t           Real_Nhits;
  Int_t           Real_det[Maxhits];
  Int_t           Real_mod[Maxhits];
  Int_t           Real_row[Maxhits];
  Int_t           Real_col[Maxhits];
  Float_t         Real_edep[Maxhits];
  Float_t         Real_t[Maxhits];
  Float_t         Real_x[Maxhits];
  Float_t         Real_y[Maxhits];
  Float_t         Real_z[Maxhits];
  
  TOut->SetBranchAddress("Event_weight",&Event_weight);
  TOut->SetBranchAddress("Event_flag",&Event_flag);
  TOut->SetBranchAddress("Primary_Nhits",&Primary_Nhits);
  TOut->SetBranchAddress("Primary_pdg",Primary_pdg);
  TOut->SetBranchAddress("Primary_E",Primary_E);
  TOut->SetBranchAddress("Primary_x",Primary_x);
  TOut->SetBranchAddress("Primary_y",Primary_y);
  TOut->SetBranchAddress("Primary_z",Primary_z);
  TOut->SetBranchAddress("Primary_px",Primary_px);
  TOut->SetBranchAddress("Primary_py",Primary_py);
  TOut->SetBranchAddress("Primary_pz",Primary_pz);
  TOut->SetBranchAddress("Virtual_Nhits",&Virtual_Nhits);
  TOut->SetBranchAddress("Virtual_pdg",Virtual_pdg);
  TOut->SetBranchAddress("Virtual_det",Virtual_det);
  TOut->SetBranchAddress("Virtual_mod",Virtual_mod);
  TOut->SetBranchAddress("Virtual_row",Virtual_row);
  TOut->SetBranchAddress("Virtual_col",Virtual_col);
  TOut->SetBranchAddress("Virtual_tid",Virtual_tid);
  TOut->SetBranchAddress("Virtual_pid",Virtual_pid);
  TOut->SetBranchAddress("Virtual_E",Virtual_E);
  TOut->SetBranchAddress("Virtual_t",Virtual_t);
  TOut->SetBranchAddress("Virtual_x",Virtual_x);
  TOut->SetBranchAddress("Virtual_y",Virtual_y);
  TOut->SetBranchAddress("Virtual_z",Virtual_z);
  TOut->SetBranchAddress("Virtual_px",Virtual_px);
  TOut->SetBranchAddress("Virtual_py",Virtual_py);
  TOut->SetBranchAddress("Virtual_pz",Virtual_pz);
  TOut->SetBranchAddress("Virtual_vx",Virtual_vx);
  TOut->SetBranchAddress("Virtual_vy",Virtual_vy);
  TOut->SetBranchAddress("Virtual_vz",Virtual_vz);
  TOut->SetBranchAddress("Real_Nhits",&Real_Nhits);
  TOut->SetBranchAddress("Real_det",Real_det);
  TOut->SetBranchAddress("Real_mod",Real_mod);
  TOut->SetBranchAddress("Real_row",Real_row);
  TOut->SetBranchAddress("Real_col",Real_col);
  TOut->SetBranchAddress("Real_edep",Real_edep);
  TOut->SetBranchAddress("Real_t",Real_t);
  TOut->SetBranchAddress("Real_x",Real_x);
  TOut->SetBranchAddress("Real_y",Real_y);
  TOut->SetBranchAddress("Real_z",Real_z);

  gStyle->SetOptFit(0);
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);

  gStyle->SetCanvasColor(0);
  gStyle->SetFrameFillColor(0);

  gStyle->SetPadTopMargin(.05);
  gStyle->SetPadLeftMargin(.18);
  gStyle->SetPadRightMargin(.18);
  gStyle->SetPadBottomMargin(.15);

  gStyle->SetTitleOffset(1.1, "X");
  gStyle->SetTitleOffset(1.5, "Y");
  gStyle->SetTitleFont(42,"X");
  gStyle->SetTitleFont(42,"Y");
  gStyle->SetTitleSize(0.055,"X");
  gStyle->SetTitleSize(0.055,"Y");

  gStyle->SetLabelOffset(0.01, "X");
  gStyle->SetLabelOffset(0.01, "Y");
  gStyle->SetLabelFont(42,"X");
  gStyle->SetLabelFont(42,"Y");
  gStyle->SetLabelSize(0.045,"X");
  gStyle->SetLabelSize(0.045,"Y");

  gStyle->SetNdivisions(105,"X");
  gStyle->SetNdivisions(105,"Y");

  gStyle->SetStripDecimals(kFALSE);

  TFile *outfile = new TFile( Form("anasFF_%d.root", run_no),"RECREATE");

  Long64_t nentries = TOut->GetEntries();
  Int_t ntrees      = TOut->GetNtrees();
  
  cout << "Processing TOut chain with " << ntrees << " trees and " << nentries << " total events" << endl;

  //-----------------------------------------------------------------------------------------------------------------------------
  
  Double_t Eb  = 6.6;
  Double_t Mp  = 0.93827;
  Double_t R2D = 180./M_PI;

  TLorentzVector Tp4(0,0,0,Mp);  // target 4vec
  TLorentzVector kp4(0,0,Eb,Eb); // beam 4vec
  TLorentzVector Qp4, kpp4, Rp4; // q, recoil electron, recoil nucleon

  TH1F* hPrim_N   = new TH1F("hPrim_N",  "", 100,0.,10.);
  TH1F* hPrim_E   = new TH1F("hPrim_E",  "", 100,0.,10000.);
  TH1F* hPrim_z   = new TH1F("hPrim_z",  "", 100,-40.,0.);
  TH1F* hPrim_pdg = new TH1F("hPrim_pdg",  "", 100,-30.,30.);
  TH2F* hPrim_xy  = new TH2F("hPrim_xy", "", 100,-1.0,1.0, 100,-1.,1. );

  TH1F* hEarm_N   = new TH1F("hEarm_N",  "", 100,0.,200.);
  TH1F* hEarm_E   = new TH1F("hEarm_E",  "", 100,0.,200.);
  TH1F* hEarm_z   = new TH1F("hEarm_z",  "", 100,0.,500.);
  TH1F* hEarm_pdg = new TH1F("hEarm_pdg",  "", 100,-30.,30.);
  TH2F* hEarm_xy  = new TH2F("hEarm_xy", "", 100,-500.,500., 100,-500.,500. );

  TH1F* hHodo_N   = new TH1F("hHodo_N",  "", 100,0.,200.);
  TH1F* hHodo_E   = new TH1F("hHodo_E",  "", 100,0.,50.);
  TH1F* hHodo_z   = new TH1F("hHodo_z",  "", 100,0.,500.);
  TH1F* hHodo_pdg = new TH1F("hHodo_pdg",  "", 100,-30.,30.);
  TH2F* hHodo_xy  = new TH2F("hHodo_xy", "", 100,-500.,500., 100,-500.,500. );

  TH1F* hHarm_N   = new TH1F("hHarm_N",  "", 100,0.,200.);
  TH1F* hHarm_E   = new TH1F("hHarm_E",  "", 100,0.,50.);
  TH1F* hHarm_z   = new TH1F("hHarm_z",  "", 100,0.,500.);
  TH1F* hHarm_pdg = new TH1F("hHarm_pdg",  "", 100,-30.,30.);
  TH2F* hHarm_xy  = new TH2F("hHarm_xy", "", 100,-500.,500., 100,-500.,500. );

  TH1F* hRealEarm_N   = new TH1F("hRealEarm_N",  "", 100,0.,200.);
  TH1F* hRealEarm_E   = new TH1F("hRealEarm_E",  "", 100,0.,100.);
  TH1F* hRealEarm_z   = new TH1F("hRealEarm_z",  "", 100,0.,500.);
  TH1F* hRealEarm_t   = new TH1F("hRealEarm_t",  "", 100, 0.,200.);
  TH2F* hRealEarm_xy  = new TH2F("hRealEarm_xy", "", 100,-500.,500., 100,-500.,500. );

  TH1F* hRealHodo_N   = new TH1F("hRealHodo_N",  "", 100,0.,200.);
  TH1F* hRealHodo_E   = new TH1F("hRealHodo_E",  "", 100,0.,0.5);
  TH1F* hRealHodo_z   = new TH1F("hRealHodo_z",  "", 100,0.,500.);
  TH1F* hRealHodo_t   = new TH1F("hRealHodo_t",  "", 100, 0.,200.);
  TH2F* hRealHodo_xy  = new TH2F("hRealHodo_xy", "", 100,-500.,500., 100,-500.,500. );

  TH1F* hRealHarm_N   = new TH1F("hRealHarm_N",  "", 100,0.,200.);
  TH1F* hRealHarm_E   = new TH1F("hRealHarm_E",  "", 100,0.,20.);
  TH1F* hRealHarm_z   = new TH1F("hRealHarm_z",  "", 100,0.,500.);
  TH1F* hRealHarm_t   = new TH1F("hRealHarm_t",  "", 100, 0.,200.);
  TH2F* hRealHarm_xy  = new TH2F("hRealHarm_xy", "", 100,-500.,500., 100,-500.,500. );

  //-----------------------------------------------------------------------------------------------------------------------------

  for(Long64_t ev=0; ev<nentries;ev++) {
    
    TOut->GetEntry(ev);
    
    if( ev%10000 == 0 )
      cout << ev << endl;

    Event_weight = Event_weight * 1/ntrees;   // this only works if nevents is the same for all trees in the chain

    if (Event_weight == 0 ) 
      Event_weight = 1;

    // Primary variables
    hPrim_N->Fill( (Float_t)Primary_Nhits );
    for( int i =0; i < Primary_Nhits; i++ ) {

      hPrim_E->Fill( Primary_E[i] , Event_weight/1e6 );
      hPrim_xy->Fill( Primary_x[i], Primary_y[i] , Event_weight/1e6 );
      hPrim_z->Fill( Primary_z[i] , Event_weight/1e6 );
      hPrim_pdg->Fill( (Float_t)Primary_pdg[i] , Event_weight/1e6 );
    }

    // Virutal variables
    for( int i =0; i < Virtual_Nhits; i++ ) {

      if( Virtual_det[i] == 0 ) {
	hEarm_N->Fill( (Float_t)Virtual_Nhits , Event_weight/1e6 );
	hEarm_E->Fill( Virtual_E[i] , Event_weight/1e6 );
	hEarm_xy->Fill( Virtual_x[i], Virtual_y[i] , Event_weight/1e6 );
	hEarm_z->Fill( Virtual_z[i] , Event_weight/1e6 );
	hEarm_pdg->Fill( (Float_t)Virtual_pdg[i] , Event_weight/1e6 );
      }

      if( Virtual_det[i] == 1 ) {
	hHodo_N->Fill( (Float_t)Virtual_Nhits , Event_weight/1e6 );
	hHodo_E->Fill( Virtual_E[i] , Event_weight/1e6 );
	hHodo_xy->Fill( Virtual_x[i], Virtual_y[i] , Event_weight/1e6 );
	hHodo_z->Fill( Virtual_z[i] , Event_weight/1e6 );
	hHodo_pdg->Fill( (Float_t)Virtual_pdg[i] , Event_weight/1e6 );
      }

      if( Virtual_det[i] == 2 ) {
	hHarm_N->Fill( (Float_t)Virtual_Nhits, Event_weight/1e6 );
	hHarm_E->Fill( Virtual_E[i], Event_weight/1e6 );
	hHarm_xy->Fill( Virtual_x[i], Virtual_y[i], Event_weight/1e6 );
	hHarm_z->Fill( Virtual_z[i], Event_weight/1e6 );
	hHarm_pdg->Fill( (Float_t)Virtual_pdg[i], Event_weight/1e6 );
      }
    }

    // Real variables
    for( int i =0; i < Real_Nhits; i++ ) {
      
      if( Real_det[i] == 0 ) {
	if( ApplyThresh && Real_edep[i] > Earm_threshold ) { 
	  if( ApplyWindow && Real_t[i] < Earm_window ) {
	    hRealEarm_N->Fill( (Float_t)Real_Nhits, Event_weight/1e6 );
	    hRealEarm_E->Fill( Real_edep[i], Event_weight/1e6 );
	    hRealEarm_xy->Fill( Real_x[i], Real_y[i], Event_weight/1e6 );
	    hRealEarm_z->Fill( Real_z[i], Event_weight/1e6 );
	    hRealEarm_t->Fill( (Float_t)Real_t[i], Event_weight/1e6 );
	  }
	}
      }
      
      if( Real_det[i] == 1 ) {
	if( ApplyThresh && Real_edep[i] > Hodo_threshold ) {
	  if( ApplyWindow && Real_t[i] < Hodo_window ) {
	    hRealHodo_N->Fill( (Float_t)Real_Nhits, Event_weight/1e6 );
	    hRealHodo_E->Fill( Real_edep[i], Event_weight/1e6 );
 	    hRealHodo_xy->Fill( Real_x[i], Real_y[i], Event_weight/1e6 );
	    hRealHodo_z->Fill( Real_z[i], Event_weight/1e6 );
	    hRealHodo_t->Fill( (Float_t)Real_t[i], Event_weight/1e6 );
	  }
	}
      }

      if( Real_det[i] == 2 ) {
	if( ApplyThresh && Real_edep[i] > Harm_threshold ) {
	  if( ApplyWindow && Real_t[i] < Harm_window ) {
	    hRealHarm_N->Fill( (Float_t)Real_Nhits, Event_weight/1e6 );
	    hRealHarm_E->Fill( Real_edep[i], Event_weight/1e6 );
	    hRealHarm_xy->Fill( Real_x[i], Real_y[i], Event_weight/1e6 );
	    hRealHarm_z->Fill( Real_z[i], Event_weight/1e6 );
	    hRealHarm_t->Fill( (Float_t)Real_t[i], Event_weight/1e6 );
	  }
	}
      }
    
    }
  }

  //-----------------------------------------------------------------------------------------------------------------------------

  TLatex* tex;

  TCanvas* cPrim = new TCanvas("cPrim","",1200,800);
  cPrim->Divide(3,2);
  
  cPrim->cd(1);
  
  tex = new TLatex( 0.24, 0.7, "sFFG4MC");
  tex->SetNDC(1);
  tex->SetTextFont(42);
  tex->SetTextColor(1);
  tex->SetTextSize(0.095);
  tex->Draw();

  tex = new TLatex( 0.24, 0.5, Form("Run %d", run_no ));
  tex->SetNDC(1);
  tex->SetTextFont(42);
  tex->SetTextColor(1);
  tex->SetTextSize(0.095);
  tex->Draw();
  
  tex = new TLatex( 0.24, 0.3, "Primaries");
  tex->SetNDC(1);
  tex->SetTextFont(42);
  tex->SetTextColor(1);
  tex->SetTextSize(0.095);
  tex->Draw();
    
  cPrim->cd(2);
  hPrim_N->SetLineColor(4);
  hPrim_N->Draw("hist c");
  hPrim_N->GetXaxis()->SetTitle("Num Primary Particles");

  cPrim->cd(3);
  hPrim_pdg->SetLineColor(4);
  hPrim_pdg->Draw("hist c");
  hPrim_pdg->GetXaxis()->SetTitle("PDG code");

  cPrim->cd(4);
  hPrim_E->SetLineColor(4);
  hPrim_E->Draw("hist c");
  hPrim_E->GetXaxis()->SetTitle("Kinetic Energy [MeV]");

  cPrim->cd(5);
  hPrim_xy->Draw("colz");
  hPrim_xy->GetXaxis()->SetTitle("x [cm]");
  hPrim_xy->GetYaxis()->SetTitle("y [cm]");

  cPrim->cd(6);
  hPrim_z->SetLineColor(4);
  hPrim_z->Draw("hist c");
  hPrim_z->GetXaxis()->SetTitle("z [cm]");

  cPrim->Print(Form("tempa-%d.pdf", run_no));  
  cPrim->Close();

  //-----------------------------------------------------------------------------------------------------------------------------

  TCanvas* cVirt1 = new TCanvas("cVirt1","",1200,800);
  cVirt1->Divide(3,2);
  
  cVirt1->cd(1);
  
  tex = new TLatex( 0.24, 0.7, "sFFG4MC");
  tex->SetNDC(1);
  tex->SetTextFont(42);
  tex->SetTextColor(1);
  tex->SetTextSize(0.095);
  tex->Draw();

  tex = new TLatex( 0.24, 0.5, Form("Run %d", run_no ));
  tex->SetNDC(1);
  tex->SetTextFont(42);
  tex->SetTextColor(1);
  tex->SetTextSize(0.095);
  tex->Draw();
  
  tex = new TLatex( 0.24, 0.3, "Virtual EArm");
  tex->SetNDC(1);
  tex->SetTextFont(42);
  tex->SetTextColor(1);
  tex->SetTextSize(0.095);
  tex->Draw();
    
  cVirt1->cd(2);
  hEarm_N->SetLineColor(4);
  hEarm_N->Draw("hist c");
  hEarm_N->GetXaxis()->SetTitle("Num Hits");

  cVirt1->cd(3);
  hEarm_pdg->SetLineColor(4);
  hEarm_pdg->Draw("hist c");
  hEarm_pdg->GetXaxis()->SetTitle("PDG code");

  cVirt1->cd(4)->SetLogy(1);
  hEarm_E->SetLineColor(4);
  hEarm_E->Draw("hist c");
  hEarm_E->GetXaxis()->SetTitle("Kinetic Energy [MeV]");

  cVirt1->cd(5);
  hEarm_xy->Draw("colz");
  hEarm_xy->GetXaxis()->SetTitle("x [cm]");
  hEarm_xy->GetYaxis()->SetTitle("y [cm]");

  cVirt1->cd(6);
  hEarm_z->SetLineColor(4);
  hEarm_z->Draw("hist c");
  hEarm_z->GetXaxis()->SetTitle("z [cm]");

  cVirt1->Print(Form("tempb-%d.pdf", run_no));  
  cVirt1->Close();

  //-----------------------------------------------------------------------------------------------------------------------------

  TCanvas* cVirt2 = new TCanvas("cVirt2","",1200,800);
  cVirt2->Divide(3,2);
  
  cVirt2->cd(1);
  
  tex = new TLatex( 0.24, 0.7, "sFFG4MC");
  tex->SetNDC(1);
  tex->SetTextFont(42);
  tex->SetTextColor(1);
  tex->SetTextSize(0.095);
  tex->Draw();

  tex = new TLatex( 0.24, 0.5, Form("Run %d", run_no ));
  tex->SetNDC(1);
  tex->SetTextFont(42);
  tex->SetTextColor(1);
  tex->SetTextSize(0.095);
  tex->Draw();
  
  tex = new TLatex( 0.24, 0.3, "Virtual Hodo");
  tex->SetNDC(1);
  tex->SetTextFont(42);
  tex->SetTextColor(1);
  tex->SetTextSize(0.095);
  tex->Draw();
    
  cVirt2->cd(2);
  hHodo_N->SetLineColor(4);
  hHodo_N->Draw("hist c");
  hHodo_N->GetXaxis()->SetTitle("Num Hits");

  cVirt2->cd(3);
  hHodo_pdg->SetLineColor(4);
  hHodo_pdg->Draw("hist c");
  hHodo_pdg->GetXaxis()->SetTitle("PDG code");

  cVirt2->cd(4)->SetLogy(1);
  hHodo_E->SetLineColor(4);
  hHodo_E->Draw("hist c");
  hHodo_E->GetXaxis()->SetTitle("Kinetic Energy [MeV]");

  cVirt2->cd(5);
  hHodo_xy->Draw("colz");
  hHodo_xy->GetXaxis()->SetTitle("x [cm]");
  hHodo_xy->GetYaxis()->SetTitle("y [cm]");

  cVirt2->cd(6);
  hHodo_z->SetLineColor(4);
  hHodo_z->Draw("hist c");
  hHodo_z->GetXaxis()->SetTitle("z [cm]");

  cVirt2->Print(Form("tempc-%d.pdf", run_no));  
  cVirt2->Close();

  //-----------------------------------------------------------------------------------------------------------------------------

  TCanvas* cVirt3 = new TCanvas("cVirt3","",1200,800);
  cVirt3->Divide(3,2);
  
  cVirt3->cd(1);
  
  tex = new TLatex( 0.24, 0.7, "sFFG4MC");
  tex->SetNDC(1);
  tex->SetTextFont(42);
  tex->SetTextColor(1);
  tex->SetTextSize(0.095);
  tex->Draw();

  tex = new TLatex( 0.24, 0.5, Form("Run %d", run_no ));
  tex->SetNDC(1);
  tex->SetTextFont(42);
  tex->SetTextColor(1);
  tex->SetTextSize(0.095);
  tex->Draw();
  
  tex = new TLatex( 0.24, 0.3, "Virtual HArm");
  tex->SetNDC(1);
  tex->SetTextFont(42);
  tex->SetTextColor(1);
  tex->SetTextSize(0.095);
  tex->Draw();
    
  cVirt3->cd(2);
  hHarm_N->SetLineColor(4);
  hHarm_N->Draw("hist c");
  hHarm_N->GetXaxis()->SetTitle("Num Hits");

  cVirt3->cd(3);
  hHarm_pdg->SetLineColor(4);
  hHarm_pdg->Draw("hist c");
  hHarm_pdg->GetXaxis()->SetTitle("PDG code");

  cVirt3->cd(4)->SetLogy(1);
  hHarm_E->SetLineColor(4);
  hHarm_E->Draw("hist c");
  hHarm_E->GetXaxis()->SetTitle("Kinetic Energy [MeV]");

  cVirt3->cd(5);
  hHarm_xy->Draw("colz");
  hHarm_xy->GetXaxis()->SetTitle("x [cm]");
  hHarm_xy->GetYaxis()->SetTitle("y [cm]");

  cVirt3->cd(6);
  hHarm_z->SetLineColor(4);
  hHarm_z->Draw("hist c");
  hHarm_z->GetXaxis()->SetTitle("z [cm]");

  cVirt3->Print(Form("tempd-%d.pdf", run_no));  
  cVirt3->Close();

  //-----------------------------------------------------------------------------------------------------------------------------

  TCanvas* cReal1 = new TCanvas("cReal1","",1200,800);
  cReal1->Divide(3,2);
  
  cReal1->cd(1);
  
  tex = new TLatex( 0.24, 0.7, "sFFG4MC");
  tex->SetNDC(1);
  tex->SetTextFont(42);
  tex->SetTextColor(1);
  tex->SetTextSize(0.095);
  tex->Draw();

  tex = new TLatex( 0.24, 0.5, Form("Run %d", run_no ));
  tex->SetNDC(1);
  tex->SetTextFont(42);
  tex->SetTextColor(1);
  tex->SetTextSize(0.095);
  tex->Draw();
  
  tex = new TLatex( 0.24, 0.3, "Real EArm");
  tex->SetNDC(1);
  tex->SetTextFont(42);
  tex->SetTextColor(1);
  tex->SetTextSize(0.095);
  tex->Draw();
    
  cReal1->cd(2);
  hRealEarm_N->SetLineColor(4);
  hRealEarm_N->Draw("hist c");
  hRealEarm_N->GetXaxis()->SetTitle("Num Hits");

  cReal1->cd(3);
  hRealEarm_t->SetLineColor(4);
  hRealEarm_t->Draw("hist c");
  hRealEarm_t->GetXaxis()->SetTitle("Time [ns]");

  cReal1->cd(4)->SetLogy(1);
  hRealEarm_E->SetLineColor(4);
  hRealEarm_E->Draw("hist c");
  hRealEarm_E->GetXaxis()->SetTitle("Energy Deposit [MeV]");

  cReal1->cd(5);
  hRealEarm_xy->Draw("colz");
  hRealEarm_xy->GetXaxis()->SetTitle("x [cm]");
  hRealEarm_xy->GetYaxis()->SetTitle("y [cm]");

  cReal1->cd(6);
  hRealEarm_z->SetLineColor(4);
  hRealEarm_z->Draw("hist c");
  hRealEarm_z->GetXaxis()->SetTitle("z [cm]");

  cReal1->Print(Form("tempe-%d.pdf", run_no));  
  cReal1->Close();

  //-----------------------------------------------------------------------------------------------------------------------------

  TCanvas* cReal2 = new TCanvas("cReal2","",1200,800);
  cReal2->Divide(3,2);
  
  cReal2->cd(1);

  tex = new TLatex( 0.24, 0.7, "sFFG4MC");
  tex->SetNDC(1);
  tex->SetTextFont(42);
  tex->SetTextColor(1);
  tex->SetTextSize(0.095);
  tex->Draw();

  tex = new TLatex( 0.24, 0.5, Form("Run %d", run_no ));
  tex->SetNDC(1);
  tex->SetTextFont(42);
  tex->SetTextColor(1);
  tex->SetTextSize(0.095);
  tex->Draw();
  
  tex = new TLatex( 0.24, 0.3, "Real Hodo");
  tex->SetNDC(1);
  tex->SetTextFont(42);
  tex->SetTextColor(1);
  tex->SetTextSize(0.095);
  tex->Draw();
    
  cReal2->cd(2);
  hRealHodo_N->SetLineColor(4);
  hRealHodo_N->Draw("hist c");
  hRealHodo_N->GetXaxis()->SetTitle("Num Hits");

  cReal2->cd(3);
  hRealHodo_t->SetLineColor(4);
  hRealHodo_t->Draw("hist c");
  hRealHodo_t->GetXaxis()->SetTitle("Time [ns]");

  cReal2->cd(4)->SetLogy(1);
  hRealHodo_E->SetLineColor(4);
  hRealHodo_E->Draw("hist c");
  hRealHodo_E->GetXaxis()->SetTitle("Energy Deposit [MeV]");

  cReal2->cd(5);
  hRealHodo_xy->Draw("colz");
  hRealHodo_xy->GetXaxis()->SetTitle("x [cm]");
  hRealHodo_xy->GetYaxis()->SetTitle("y [cm]");

  cReal2->cd(6);
  hRealHodo_z->SetLineColor(4);
  hRealHodo_z->Draw("hist c");
  hRealHodo_z->GetXaxis()->SetTitle("z [cm]");

  cReal2->Print(Form("tempf-%d.pdf", run_no));  
  cReal2->Close();

  //-----------------------------------------------------------------------------------------------------------------------------

  TCanvas* cReal3 = new TCanvas("cReal3","",1200,800);
  cReal3->Divide(3,2);
  
  cReal3->cd(1);
  
  tex = new TLatex( 0.24, 0.7, "sFFG4MC");
  tex->SetNDC(1);
  tex->SetTextFont(42);
  tex->SetTextColor(1);
  tex->SetTextSize(0.095);
  tex->Draw();

  tex = new TLatex( 0.24, 0.5, Form("Run %d", run_no ));
  tex->SetNDC(1);
  tex->SetTextFont(42);
  tex->SetTextColor(1);
  tex->SetTextSize(0.095);
  tex->Draw();
  
  tex = new TLatex( 0.24 , 0.3, "Real HArm");
  tex->SetNDC(1);
  tex->SetTextFont(42);
  tex->SetTextColor(1);
  tex->SetTextSize(0.095);
  tex->Draw();
    
  cReal3->cd(2);
  hRealHarm_N->SetLineColor(4);
  hRealHarm_N->Draw("hist c");
  hRealHarm_N->GetXaxis()->SetTitle("Num Hits");

  cReal3->cd(3);
  hRealHarm_t->SetLineColor(4);
  hRealHarm_t->Draw("hist c");
  hRealHarm_t->GetXaxis()->SetTitle("Time [ns]");

  cReal3->cd(4)->SetLogy(1);
  hRealHarm_E->SetLineColor(4);
  hRealHarm_E->Draw("hist c");
  hRealHarm_E->GetXaxis()->SetTitle("Energy Deposit [MeV]");

  cReal3->cd(5);
  hRealHarm_xy->Draw("colz");
  hRealHarm_xy->GetXaxis()->SetTitle("x [cm]");
  hRealHarm_xy->GetYaxis()->SetTitle("y [cm]");

  cReal3->cd(6);
  hRealHarm_z->SetLineColor(4);
  hRealHarm_z->Draw("hist c");
  hRealHarm_z->GetXaxis()->SetTitle("z [cm]");

  cReal3->Print(Form("tempg-%d.pdf", run_no));  
  cReal3->Close();
  
  //-----------------------------------------------------------------------------------------------------------------------------
  
  gSystem->Exec(Form("pdfunite  temp*.pdf sFFG4MC_%d.pdf", run_no) );  
  gSystem->Exec("rm temp*.pdf");  
  
  outfile->Write();
  outfile->Close();
   
}

//-----------------------------------------------------------------------------------------------------------------------------
