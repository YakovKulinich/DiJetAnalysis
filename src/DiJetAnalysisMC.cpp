#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TROOT.h>
#include <TEnv.h>
#include <TCanvas.h>
#include <TRandom.h>

#include <iostream>

#include "MyRoot.h"

#include "DiJetAnalysisMC.h"

DiJetAnalysisMC::DiJetAnalysisMC()
  : DiJetAnalysis( true, true ) , nJZNmax(3)
{}

DiJetAnalysisMC::DiJetAnalysisMC( bool isData, bool is_pPb )
  : DiJetAnalysis( isData, is_pPb) , nJZNmax(3)
{}

DiJetAnalysisMC::~DiJetAnalysisMC(){}

void DiJetAnalysisMC::RunOverTreeFillHistos( int nEvents, 
					     int startEvent ){  
    setupHistograms();
    processEvents( nEvents, startEvent );
}

void DiJetAnalysisMC::PlotExistingHistos(){
  loadHistograms();
    
  plotSpectra();
  plotEtaPhi();
  plotEtaPt();
}

void DiJetAnalysisMC::setupHistograms(){
  /*
  THmulf* hm = new THmulf( "test","test" );
  hm->AddAxis( "pt1"  , "p_{T}^{1}", 50, 10, 60);
  hm->AddAxis( "pt2"  , "p_{T}^{1}", 50, 10, 60);
  hm->AddAxis( "pt3"  , "p_{T}^{1}", 50, 10, 60);
  hm->AddAxis( "pt4"  , "p_{T}^{1}", 50, 10, 60);
  hm->Sumw2();
  
  TRandom r;
  
  for( int i = 0; i < 10000; i++ ){
    hm->Fill( 1,
	      r.Gaus(35,15), r.Gaus(35,15),
	      r.Gaus(35,15), r.Gaus(35,15));
  }

  TH1* h = hm->Projection( "h", "pt1" );
    
  std::cout << "fNameOut: " << m_fNameOut << std::endl;
  TFile* fout = new TFile( m_fNameOut.c_str(),"RECREATE");
  h->Write();
  
  fout->Close();
  */  

  
}

void DiJetAnalysisMC::processEvents( int nEvents, int startEvent ){}

void DiJetAnalysisMC::loadHistograms(){}

void DiJetAnalysisMC::plotSpectra(){}

void DiJetAnalysisMC::plotEtaPhi(){}

void DiJetAnalysisMC::plotEtaPt(){}
