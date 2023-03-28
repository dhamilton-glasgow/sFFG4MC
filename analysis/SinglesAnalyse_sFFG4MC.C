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

const Float_t Earm_threshold = 28.;
const Float_t Earm_window    = 50;
const Float_t Harm_threshold = 3.0;
const Float_t Harm_window    = 50;
const Float_t Hodo_threshold = 0.;
const Float_t Hodo_window    = 50;

const Int_t   NEarm          = 961;
const Int_t   NHodo          = 7201;
const Int_t   NHarm          = 289;

void SinglesAnalyse_sFFG4MC( Int_t run_no = 1 ) { 
  
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
  gStyle->SetTitleOffset(1.7, "Y");
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

  Long64_t nentries = TOut->GetEntries();
  Int_t ntrees      = TOut->GetNtrees();
  
  cout << "Processing TOut chain with " << ntrees << " trees and " << nentries << " total events" << endl;

  //-----------------------------------------------------------------------------------------------------------------------------
  
  TH1F* hRealEarm_N   = new TH1F("hRealEarm_N",  "", NEarm,0.,NEarm);
  TH1F* hRealEarm_E   = new TH1F("hRealEarm_E",  "", 100,0.,1000.);
  TH1F* hRealEarm_z   = new TH1F("hRealEarm_z",  "", 100,0.,800.);
  TH1F* hRealEarm_t   = new TH1F("hRealEarm_t",  "", 100, 0.,50.);
  TH2F* hRealEarm_xy  = new TH2F("hRealEarm_xy", "", 100,-500.,500., 100,-500.,500. );
  TH1F* hRealEarm_Edet[NEarm];
  for(Int_t i=0; i<NEarm; i++) 
    hRealEarm_Edet[i] = new TH1F( Form("hRealEarm_Edet%d", i), "", 100,0.,1000.);

  TH1F* hRealHarm_N   = new TH1F("hRealHarm_N",  "", NHarm,0.,NHarm);
  TH1F* hRealHarm_E   = new TH1F("hRealHarm_E",  "", 100,0.,100.);
  TH1F* hRealHarm_z   = new TH1F("hRealHarm_z",  "", 100,0.,500.);
  TH1F* hRealHarm_t   = new TH1F("hRealHarm_t",  "", 100, 0.,50.);
  TH2F* hRealHarm_xy  = new TH2F("hRealHarm_xy", "", 100,-500.,500., 100,-500.,500. );
  TH1F* hRealHarm_Edet[NHarm];
  
  for(Int_t i=0; i<NHarm; i++) 
    hRealHarm_Edet[i] = new TH1F( Form("hRealHarm_Edet%d", i), "", 100,0.,1000.);

  TH1F* hRealHodo_N   = new TH1F("hRealHodo_N",  "", 7200,0.,NHodo);
  TH1F* hRealHodo_E   = new TH1F("hRealHodo_E",  "", 100,0.,0.3);
  TH1F* hRealHodo_z   = new TH1F("hRealHodo_z",  "", 100,0.,500.);
  TH1F* hRealHodo_t   = new TH1F("hRealHodo_t",  "", 100, 0.,50.);
  TH2F* hRealHodo_xy  = new TH2F("hRealHodo_xy", "", 100,-500.,500., 100,-500.,500. );
  TH1F* hRealHodo_Edet[NHodo];
  for(Int_t i=0; i<NHodo; i++) 
    hRealHodo_Edet[i] = new TH1F( Form("hRealHodo_Edet%d", i), "", 100,0.,0.3);
  
  //-----------------------------------------------------------------------------------------------------------------------------

  //  for(Long64_t ev=0; ev<10;ev++) {
  for(Long64_t ev=0; ev<nentries;ev++) {
    
    TOut->GetEntry(ev);
    
    if( ev%100000 == 0 ) {
     printf("Event %8lld\r", ev);
     fflush(stdout);
    }

    Event_weight = Event_weight * 1/ntrees;   // this only works if nevents is the same for all trees in the chain

    if (Event_weight == 0 ) 
      Event_weight = 1;

    for( int i =0; i < Real_Nhits; i++ ) {
      
      if( Real_det[i] == 0 ) {

	if( ApplyThresh && Real_edep[i] > Earm_threshold ) { 
	  if( ApplyWindow && Real_t[i] < Earm_window ) {

	    Int_t id = (Real_mod[i]*160)+((Real_row[i]+1)*(Real_col[i]+1));
	    hRealEarm_Edet[id]->Fill( Real_edep[i], Event_weight );

	    hRealEarm_N->Fill( (Float_t)id, Event_weight );
	    hRealEarm_E->Fill( Real_edep[i], Event_weight );
	    hRealEarm_xy->Fill( Real_x[i], Real_y[i], Event_weight );
	    hRealEarm_z->Fill( Real_z[i], Event_weight );
	    hRealEarm_t->Fill( (Float_t)Real_t[i], Event_weight );
	  }
	}
      }
      
      if( Real_det[i] == 1 ) {
	if( ApplyThresh && Real_edep[i] > Hodo_threshold ) {
	  if( ApplyWindow && Real_t[i] < Hodo_window ) {
	    Int_t id = (Real_mod[i]*1200)+((Real_row[i]+1)*(Real_col[i]+1));

	    hRealHodo_Edet[id]->Fill( Real_edep[i], Event_weight );

	    hRealHodo_N->Fill( (Float_t)id, Event_weight );
	    hRealHodo_E->Fill( Real_edep[i], Event_weight );
 	    hRealHodo_xy->Fill( Real_x[i], Real_y[i], Event_weight );
	    hRealHodo_z->Fill( Real_z[i], Event_weight );
	    hRealHodo_t->Fill( (Float_t)Real_t[i], Event_weight );
	  }
	}
      }

      if( Real_det[i] == 2 ) {
	if( ApplyThresh && Real_edep[i] > Harm_threshold ) {
	  if( ApplyWindow && Real_t[i] < Harm_window ) {

	    Int_t id = (Real_mod[i]*48)+((Real_row[i]+1)*(Real_col[i]+1));
	    hRealHarm_Edet[id]->Fill( Real_edep[i], Event_weight );

	    hRealHarm_N->Fill( (Float_t)id, Event_weight );
	    hRealHarm_E->Fill( Real_edep[i], Event_weight );
	    hRealHarm_xy->Fill( Real_x[i], Real_y[i], Event_weight );
	    hRealHarm_z->Fill( Real_z[i], Event_weight );
	    hRealHarm_t->Fill( (Float_t)Real_t[i], Event_weight );
	  }
	}
      }
    
    }
  }

  //-----------------------------------------------------------------------------------------------------------------------------

  TLatex* tex;

  TCanvas* cReal1 = new TCanvas("cReal1","",1200,800);
  cReal1->Divide(3,2);
  
  cReal1->cd(1);
  
  tex = new TLatex( 0.09, 0.9, "sFFG4MC");
  tex->SetNDC(1);
  tex->SetTextFont(42);
  tex->SetTextColor(1);
  tex->SetTextSize(0.075);
  tex->Draw();

  tex = new TLatex( 0.09, 0.7, "Beam background");
  tex->SetNDC(1);
  tex->SetTextFont(42);
  tex->SetTextColor(1);
  tex->SetTextSize(0.075);
  tex->Draw();
  
  tex = new TLatex( 0.09, 0.5, "EArm Singles Rates");
  tex->SetNDC(1);
  tex->SetTextFont(42);
  tex->SetTextColor(1);
  tex->SetTextSize(0.075);
  tex->Draw();

  tex = new TLatex( 0.09, 0.3, Form("Threshold = %3.2f MeV", Earm_threshold));
  tex->SetNDC(1);
  tex->SetTextFont(42);
  tex->SetTextColor(1);
  tex->SetTextSize(0.075);
  tex->Draw();

  cReal1->cd(2)->SetLogy(1);
  hRealEarm_E->SetLineColor(4);
  hRealEarm_E->Draw("hist e4");
  hRealEarm_E->GetXaxis()->SetTitle("Energy Deposit [MeV]");

  tex = new TLatex( 0.24, 0.9, Form("Mean = %3.2e MeV", hRealEarm_E->GetMean()));
  tex->SetNDC(1);
  tex->SetTextFont(42);
  tex->SetTextColor(1);
  tex->SetTextSize(0.055);
  tex->Draw();

  cReal1->cd(3);
  hRealEarm_N->SetLineColor(4);
  hRealEarm_N->Draw("hist e4");
  hRealEarm_N->GetYaxis()->SetTitle("Hit Rate per uA [Hz]");
  hRealEarm_N->GetXaxis()->SetTitle("Detector ID");
  hRealEarm_N->GetYaxis()->SetRangeUser(0, 1.3*hRealEarm_N->GetMaximum());

  TF1* meanhitrate = new TF1("meanhitrate","pol0", 0, 1000);  
  meanhitrate->SetLineColor( 1 );
  hRealEarm_N->Fit(meanhitrate,"Q","",0, 1000);;
  
  tex = new TLatex( 0.24, 0.9, Form("Mean = %3.2e Hz", meanhitrate->GetParameter(0)/1.));
  tex->SetNDC(1);
  tex->SetTextFont(42);
  tex->SetTextColor(1);
  tex->SetTextSize(0.055);
  tex->Draw();

  cReal1->cd(4);

  Float_t det[NHodo], sumedep[NHodo], sumdose[NHodo] ;
  Float_t edet[NHodo], esumedep[NHodo], esumdose[NHodo] ;

  for(Int_t i=0; i<NEarm; i++) {
    det[i]     = i;
    sumedep[i] = (hRealEarm_Edet[i]->Integral(0,100))* 16e-13; // convert MeV to J
    sumdose[i] = (sumedep[i] * 60 * 60 * 100 )/(0.5*2*2*20*8/1000.); // divide by detector mass in kg then convert J/s to rad/hr
    // NB the 0.5 in the denominator accounts for the dose being deposited in the front half of the detector
    edet[i]     = 0.;
    esumedep[i] = 0.;
    
  }

  TGraphErrors* gEdep = new TGraphErrors( NEarm, det, sumedep, edet, esumedep );
  gEdep->SetMarkerColor( 4 );
  gEdep->SetLineColor( 4 );
  gEdep->SetMarkerStyle( 20 );
  gEdep->SetMarkerSize( 0.8 );
  gEdep->Draw("AL");
  gEdep->GetXaxis()->SetTitle( "Detector ID");
  gEdep->GetYaxis()->SetTitle( "Edep Rate per #muA [J/s]");
  gEdep->GetXaxis()->SetRangeUser(0,NEarm);
  gEdep->GetYaxis()->SetRangeUser(0, 1.3*TMath::MaxElement(NEarm,gEdep->GetY()) );

  TF1* meanedeprate = new TF1("meanedeprate","pol0", 0, 1000);  
  meanedeprate->SetLineColor( 1 );
  gEdep->Fit(meanedeprate,"NQ","",0, 1000);;
  
  tex = new TLatex( 0.24, 0.9, Form("Mean = %3.2e J/s", meanedeprate->GetParameter(0)));
  tex->SetNDC(1);
  tex->SetTextFont(42);
  tex->SetTextColor(1);
  tex->SetTextSize(0.055);
  tex->Draw();

  cReal1->cd(5);

  TGraphErrors* gDose = new TGraphErrors( NEarm, det, sumdose, edet, esumdose );
  gDose->SetMarkerColor( 4 );
  gDose->SetLineColor( 4 );
  gDose->SetMarkerStyle( 20 );
  gDose->SetMarkerSize( 0.8 );
  gDose->Draw("AL");
  gDose->GetXaxis()->SetTitle( "Detector ID");
  gDose->GetYaxis()->SetTitle( "Dose Rate per #muA [rad/hr]");
  gDose->GetXaxis()->SetRangeUser(0,NEarm);
  gDose->GetYaxis()->SetRangeUser(0, 1.3*TMath::MaxElement(NEarm,gDose->GetY()) );

  TF1* meandoserate = new TF1("meandoserate","pol0", 0, 1000);  
  meandoserate->SetLineColor( 1 );
  gDose->Fit(meandoserate,"NQ","",0, 1000);;
  
  tex = new TLatex( 0.24, 0.9, Form("Mean = %3.2e rad/hr", meandoserate->GetParameter(0)));
  tex->SetNDC(1);
  tex->SetTextFont(42);
  tex->SetTextColor(1);
  tex->SetTextSize(0.055);
  tex->Draw();

  cReal1->cd(6);

  tex = new TLatex( 0.09, 0.9, "720 hours at 60 #muA");
  tex->SetNDC(1);
  tex->SetTextFont(42);
  tex->SetTextColor(1);
  tex->SetTextSize(0.075);
  tex->Draw();

  tex = new TLatex( 0.09, 0.6, Form("Mean hit rate = %3.2e Hz", 60.*meanhitrate->GetParameter(0)/1.));
  tex->SetNDC(1);
  tex->SetTextFont(42);
  tex->SetTextColor(1);
  tex->SetTextSize(0.055);
  tex->Draw();

  tex = new TLatex( 0.09, 0.5, Form("Mean edep rate = %3.2e J/s", 60.*meanedeprate->GetParameter(0)));
  tex->SetNDC(1);
  tex->SetTextFont(42);
  tex->SetTextColor(1);
  tex->SetTextSize(0.055);
  tex->Draw();

  tex = new TLatex( 0.09, 0.4, Form("Mean dose rate = %3.2e rad/hr", 60.*meandoserate->GetParameter(0)));
  tex->SetNDC(1);
  tex->SetTextFont(42);
  tex->SetTextColor(1);
  tex->SetTextSize(0.055);
  tex->Draw();

  tex = new TLatex( 0.09, 0.3, Form("Total dose = %3.2e rad", 60.*720*meandoserate->GetParameter(0)));
  tex->SetNDC(1);
  tex->SetTextFont(42);
  tex->SetTextColor(1);
  tex->SetTextSize(0.055);
  tex->Draw();
  
  cReal1->Print(Form("SinglesEarm-%2.2f.pdf", Earm_threshold));  
  cReal1->Print(Form("SinglesEarm-%2.2f.png", Earm_threshold));  

  //-----------------------------------------------------------------------------------------------------------------------------

  TCanvas* cReal2 = new TCanvas("cReal2","",1200,800);
  cReal2->Divide(3,2);
  
  cReal2->cd(1);
  
  tex = new TLatex( 0.09, 0.9, "sFFG4MC");
  tex->SetNDC(1);
  tex->SetTextFont(42);
  tex->SetTextColor(1);
  tex->SetTextSize(0.075);
  tex->Draw();

  tex = new TLatex( 0.09, 0.7, "Beam background");
  tex->SetNDC(1);
  tex->SetTextFont(42);
  tex->SetTextColor(1);
  tex->SetTextSize(0.075);
  tex->Draw();
  
  tex = new TLatex( 0.09, 0.5, "HArm Singles Rates");
  tex->SetNDC(1);
  tex->SetTextFont(42);
  tex->SetTextColor(1);
  tex->SetTextSize(0.075);
  tex->Draw();

  tex = new TLatex( 0.09, 0.3, Form("Threshold = %3.2f MeV", Harm_threshold));
  tex->SetNDC(1);
  tex->SetTextFont(42);
  tex->SetTextColor(1);
  tex->SetTextSize(0.075);
  tex->Draw();

  cReal2->cd(2)->SetLogy(1);
  hRealHarm_E->SetLineColor(4);
  hRealHarm_E->Draw("hist e4");
  hRealHarm_E->GetXaxis()->SetTitle("Energy Deposit [MeV]");

  tex = new TLatex( 0.24, 0.9, Form("Mean = %3.2e MeV", hRealHarm_E->GetMean()));
  tex->SetNDC(1);
  tex->SetTextFont(42);
  tex->SetTextColor(1);
  tex->SetTextSize(0.055);
  tex->Draw();

  cReal2->cd(3);
  hRealHarm_N->SetLineColor(4);
  hRealHarm_N->Draw("hist e4");
  hRealHarm_N->GetYaxis()->SetTitle("Hit Rate per uA [Hz]");
  hRealHarm_N->GetXaxis()->SetTitle("Detector ID");
  hRealHarm_N->GetYaxis()->SetRangeUser(0, 1.3*hRealHarm_N->GetMaximum());

  TF1* hmeanhitrate = new TF1("hmeanhitrate","pol0", 0, 300);  
  hmeanhitrate->SetLineColor( 1 );
  hRealHarm_N->Fit(hmeanhitrate,"Q","",0, 300);;
  
  tex = new TLatex( 0.24, 0.9, Form("Mean = %3.2e Hz", hmeanhitrate->GetParameter(0)/1.));
  tex->SetNDC(1);
  tex->SetTextFont(42);
  tex->SetTextColor(1);
  tex->SetTextSize(0.055);
  tex->Draw();


  cReal2->cd(4);

  for(Int_t i=0; i<NHarm; i++) {
    det[i]     = i;
    sumedep[i] = (hRealHarm_Edet[i]->Integral(0,100))* 16e-13; // convert MeV to J
    sumdose[i] = (sumedep[i] * 60 * 60 * 100 )/(0.2*15*15*100*4/1000.); // divide by detector mass in kg then convert J/s to rad/hr
    // NB the 0.2 in the denominator accounts for the dose being deposited in the front 20% of the detector
    edet[i]     = 0.;
    esumedep[i] = 0.;
    
  }

  gEdep = new TGraphErrors( NHarm, det, sumedep, edet, esumedep );
  gEdep->SetMarkerColor( 4 );
  gEdep->SetLineColor( 4 );
  gEdep->SetMarkerStyle( 20 );
  gEdep->SetMarkerSize( 0.8 );
  gEdep->Draw("AL");
  gEdep->GetXaxis()->SetTitle( "Detector ID");
  gEdep->GetYaxis()->SetTitle( "Edep Rate per #muA [J/s]");
  gEdep->GetXaxis()->SetRangeUser(0,NHarm);
  gEdep->GetYaxis()->SetRangeUser(0, 1.3*TMath::MaxElement(NHarm,gEdep->GetY()) );

  TF1* hmeanedeprate = new TF1("hmeanedeprate","pol0", 0, 300);  
  hmeanedeprate->SetLineColor( 1 );
  gEdep->Fit(hmeanedeprate,"NQ","",0, 300);;
  
  tex = new TLatex( 0.24, 0.9, Form("Mean = %3.2e J/s", hmeanedeprate->GetParameter(0)));
  tex->SetNDC(1);
  tex->SetTextFont(42);
  tex->SetTextColor(1);
  tex->SetTextSize(0.055);
  tex->Draw();

  cReal2->cd(5);

  gDose = new TGraphErrors( NHarm, det, sumdose, edet, esumdose );
  gDose->SetMarkerColor( 4 );
  gDose->SetLineColor( 4 );
  gDose->SetMarkerStyle( 20 );
  gDose->SetMarkerSize( 0.8 );
  gDose->Draw("AL");
  gDose->GetXaxis()->SetTitle( "Detector ID");
  gDose->GetYaxis()->SetTitle( "Dose Rate per #muA [rad/hr]");
  gDose->GetXaxis()->SetRangeUser(0,NHarm);
  gDose->GetYaxis()->SetRangeUser(0, 1.3*TMath::MaxElement(NHarm,gDose->GetY()) );

  TF1* hmeandoserate = new TF1("hmeandoserate","pol0", 0, 1000);  
  hmeandoserate->SetLineColor( 1 );
  gDose->Fit(hmeandoserate,"NQ","",0, 1000);;
  
  tex = new TLatex( 0.24, 0.9, Form("Mean = %3.2e rad/hr", hmeandoserate->GetParameter(0)));
  tex->SetNDC(1);
  tex->SetTextFont(42);
  tex->SetTextColor(1);
  tex->SetTextSize(0.055);
  tex->Draw();

  cReal2->cd(6);

  tex = new TLatex( 0.09, 0.9, "720 hours at 60 #muA");
  tex->SetNDC(1);
  tex->SetTextFont(42);
  tex->SetTextColor(1);
  tex->SetTextSize(0.075);
  tex->Draw();

  tex = new TLatex( 0.09, 0.6, Form("Mean hit rate = %3.2e Hz", 60.*hmeanhitrate->GetParameter(0)/1.));
  tex->SetNDC(1);
  tex->SetTextFont(42);
  tex->SetTextColor(1);
  tex->SetTextSize(0.055);
  tex->Draw();

  tex = new TLatex( 0.09, 0.5, Form("Mean edep rate = %3.2e J/s", 60.*hmeanedeprate->GetParameter(0)));
  tex->SetNDC(1);
  tex->SetTextFont(42);
  tex->SetTextColor(1);
  tex->SetTextSize(0.055);
  tex->Draw();

  tex = new TLatex( 0.09, 0.4, Form("Mean dose rate = %3.2e rad/hr", 60.*hmeandoserate->GetParameter(0)));
  tex->SetNDC(1);
  tex->SetTextFont(42);
  tex->SetTextColor(1);
  tex->SetTextSize(0.055);
  tex->Draw();

  tex = new TLatex( 0.09, 0.3, Form("Total dose = %3.2e rad", 60.*720*hmeandoserate->GetParameter(0)));
  tex->SetNDC(1);
  tex->SetTextFont(42);
  tex->SetTextColor(1);
  tex->SetTextSize(0.055);
  tex->Draw();
  
  cReal2->Print(Form("SinglesHarm-%2.2f.pdf", Harm_threshold));  
  cReal2->Print(Form("SinglesHarm-%2.2f.png", Harm_threshold));  

  //-----------------------------------------------------------------------------------------------------------------------------

  TCanvas* cReal3 = new TCanvas("cReal3","",1200,800);
  cReal3->Divide(3,2);
  
  cReal3->cd(1);
  
  tex = new TLatex( 0.09, 0.9, "sFFG4MC");
  tex->SetNDC(1);
  tex->SetTextFont(42);
  tex->SetTextColor(1);
  tex->SetTextSize(0.075);
  tex->Draw();

  tex = new TLatex( 0.09, 0.7, "Beam background");
  tex->SetNDC(1);
  tex->SetTextFont(42);
  tex->SetTextColor(1);
  tex->SetTextSize(0.075);
  tex->Draw();
  
  tex = new TLatex( 0.09, 0.5, "Hodo Singles Rates");
  tex->SetNDC(1);
  tex->SetTextFont(42);
  tex->SetTextColor(1);
  tex->SetTextSize(0.075);
  tex->Draw();

  tex = new TLatex( 0.09, 0.3, Form("Threshold = %3.2f MeV", Hodo_threshold));
  tex->SetNDC(1);
  tex->SetTextFont(42);
  tex->SetTextColor(1);
  tex->SetTextSize(0.075);
  tex->Draw();

  cReal3->cd(2)->SetLogy(1);
  hRealHodo_E->SetLineColor(4);
  hRealHodo_E->Draw("hist e4");
  hRealHodo_E->GetXaxis()->SetTitle("Energy Deposit [MeV]");

  tex = new TLatex( 0.24, 0.9, Form("Mean = %3.2e MeV", hRealHodo_E->GetMean()));
  tex->SetNDC(1);
  tex->SetTextFont(42);
  tex->SetTextColor(1);
  tex->SetTextSize(0.055);
  tex->Draw();

  cReal3->cd(3);
  hRealHodo_N->SetLineColor(4);
  hRealHodo_N->Draw("hist e4");
  hRealHodo_N->GetYaxis()->SetTitle("Hit Rate per uA [Hz]");
  hRealHodo_N->GetXaxis()->SetTitle("Detector ID");
  hRealHodo_N->GetYaxis()->SetRangeUser(0, 1.3*hRealHodo_N->GetMaximum());

  TF1* hodmeannhitrate = new TF1("hodmeannhitrate","pol0", 0, NHodo);  
  hodmeannhitrate->SetLineColor( 1 );
  hRealHodo_N->Fit(hodmeannhitrate,"Q","",0, NHodo);;
  
  tex = new TLatex( 0.24, 0.9, Form("Mean = %3.2e Hz", hodmeannhitrate->GetParameter(0)/1.));
  tex->SetNDC(1);
  tex->SetTextFont(42);
  tex->SetTextColor(1);
  tex->SetTextSize(0.055);
  tex->Draw();

  cReal3->cd(4);

  for(Int_t i=0; i<NHodo; i++) {
    det[i]     = i;
    sumedep[i] = (hRealHodo_Edet[i]->Integral(0,100))* 16e-13; // convert MeV to J
    sumdose[i] = (sumedep[i] * 60 * 60 * 100 )/(3*3*10*1/1000.); // divide by detector mass in kg then convert J/s to rad/hr
    edet[i]     = 0.;
    esumedep[i] = 0.;
    
  }

  gEdep = new TGraphErrors( NHodo, det, sumedep, edet, esumedep );
  gEdep->SetMarkerColor( 4 );
  gEdep->SetLineColor( 4 );
  gEdep->SetMarkerStyle( 20 );
  gEdep->SetMarkerSize( 0.8 );
  gEdep->Draw("AL");
  gEdep->GetXaxis()->SetTitle( "Detector ID");
  gEdep->GetYaxis()->SetTitle( "Edep Rate per #muA [J/s]");
  gEdep->GetXaxis()->SetRangeUser(0,NHodo);
  gEdep->GetYaxis()->SetRangeUser(0, 1.3*TMath::MaxElement(NHodo,gEdep->GetY()) );

  TF1* hodmeannedeprate = new TF1("hodmeannedeprate","pol0", 0, NHodo);  
  hodmeannedeprate->SetLineColor( 1 );
  gEdep->Fit(hodmeannedeprate,"NQ","",0, NHodo);;
  
  tex = new TLatex( 0.24, 0.9, Form("Mean = %3.2e J/s", hodmeannedeprate->GetParameter(0)));
  tex->SetNDC(1);
  tex->SetTextFont(42);
  tex->SetTextColor(1);
  tex->SetTextSize(0.055);
  tex->Draw();

  cReal3->cd(5);

  gDose = new TGraphErrors( NHodo, det, sumdose, edet, esumdose );
  gDose->SetMarkerColor( 4 );
  gDose->SetLineColor( 4 );
  gDose->SetMarkerStyle( 20 );
  gDose->SetMarkerSize( 0.8 );
  gDose->Draw("AL");
  gDose->GetXaxis()->SetTitle( "Detector ID");
  gDose->GetYaxis()->SetTitle( "Dose Rate per #muA [rad/hr]");
  gDose->GetXaxis()->SetRangeUser(0,NHodo);
  gDose->GetYaxis()->SetRangeUser(0, 1.3*TMath::MaxElement(NHodo,gDose->GetY()) );

  TF1* hodmeanndoserate = new TF1("hodmeanndoserate","pol0", 0, NHodo);  
  hodmeanndoserate->SetLineColor( 1 );
  gDose->Fit(hodmeanndoserate,"NQ","",0, NHodo);;
  
  tex = new TLatex( 0.24, 0.9, Form("Mean = %3.2e rad/hr", hodmeanndoserate->GetParameter(0)));
  tex->SetNDC(1);
  tex->SetTextFont(42);
  tex->SetTextColor(1);
  tex->SetTextSize(0.055);
  tex->Draw();

  cReal3->cd(6);

  tex = new TLatex( 0.09, 0.9, "720 hours at 60 #muA");
  tex->SetNDC(1);
  tex->SetTextFont(42);
  tex->SetTextColor(1);
  tex->SetTextSize(0.075);
  tex->Draw();

  tex = new TLatex( 0.09, 0.6, Form("Mean hit rate = %3.2e Hz", 60.*hodmeannhitrate->GetParameter(0)/1.));
  tex->SetNDC(1);
  tex->SetTextFont(42);
  tex->SetTextColor(1);
  tex->SetTextSize(0.055);
  tex->Draw();

  tex = new TLatex( 0.09, 0.5, Form("Mean edep rate = %3.2e #muJ/s", 60.*hodmeannedeprate->GetParameter(0)));
  tex->SetNDC(1);
  tex->SetTextFont(42);
  tex->SetTextColor(1);
  tex->SetTextSize(0.055);
  tex->Draw();

  tex = new TLatex( 0.09, 0.4, Form("Mean dose rate = %3.2e rad/hr", 60.*hodmeanndoserate->GetParameter(0)));
  tex->SetNDC(1);
  tex->SetTextFont(42);
  tex->SetTextColor(1);
  tex->SetTextSize(0.055);
  tex->Draw();

  tex = new TLatex( 0.09, 0.3, Form("Total dose = %3.2e rad", 60.*720*hodmeanndoserate->GetParameter(0)));
  tex->SetNDC(1);
  tex->SetTextFont(42);
  tex->SetTextColor(1);
  tex->SetTextSize(0.055);
  tex->Draw();
  
  cReal3->Print(Form("SinglesHodo-%2.2f.pdf", Hodo_threshold));  
  cReal3->Print(Form("SinglesHodo-%2.2f.png", Hodo_threshold));  
   

  }


//-----------------------------------------------------------------------------------------------------------------------------


