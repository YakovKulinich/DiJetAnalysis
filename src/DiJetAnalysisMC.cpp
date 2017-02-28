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
#include "JetPair.h"

DiJetAnalysisMC::DiJetAnalysisMC()
{}

DiJetAnalysisMC::DiJetAnalysisMC( bool isData, bool is_pPb )
  : DiJetAnalysis( isData, is_pPb)
{}

DiJetAnalysisMC::~DiJetAnalysisMC(){}

//---------------------------------
//            Read Data
//---------------------------------
void DiJetAnalysisMC::RunOverTreeFillHistos( int nEvents, 
					     int startEvent ){  
    SetupHistograms();
    ProcessEvents( nEvents, startEvent );
    SaveOutputs();
}

//---------------------------------
//            Plot Data
//---------------------------------
void DiJetAnalysisMC::ProcessPlotHistos(){
  LoadHistograms();
  PlotSpectra();
  PlotEtaPhi();
  PlotEtaPt();
}

//---------------------------------
//            Read Data
//---------------------------------

void DiJetAnalysisMC::SetupHistograms(){
  m_jznFname[1] = "/home/yakov/Projects/atlas/data/jz1.root";
  m_jznSigma[1] = 6.7890E+07;
  m_jznEff[1]   = 2.8289E-03;

  // Triggers and Spectra
  int    nPtSpectBins = 100; 
  double ptSpectMin   = 0; double ptSpectMax  = nPtSpectBins;
  
  // Eta-Phi Maps
  int    nEtaBins = 100; 
  double etaMin   = -5; double etaMax  = 5;

  int    nPhiBins = 64; 
  double phiMin   = -constants::PI; double phiMax  = constants::PI; 

  int    nPtBins  = 100;
  double ptMin    = 10;
  double ptMax    = ptMin + static_cast<double>(nPtBins)/2;

  // JES JER etc
  int    nEtaForwardBinsFine   = 12;
  int    nEtaForwardBinsCoarse = 3;
  double etaForwardMin   = -constants::FETAMAX;
  double etaForwardMax   = -constants::FETAMIN;

  double ptTruthWidth  = 5;
  double ptTruthMin    = 10;
  double ptTruthMax    = 100;
  int    nPtTruthBins  = (ptTruthMax - ptTruthMin) / ptTruthWidth;

  int    nRPtBins   = 60;
  double rPtMin     = 0; double rPtMax = 3;

  int    nDAngleBins  = 50;
  double dAngleMin    = -0.5; double dAngleMax = -0.5;
  
  for( auto& jznFname : m_jznFname ){ // loop over JZ samples
    int jzn = jznFname.first;
    
    m_jznRPt[ jzn ] =
      new TH3D( Form("h_rPt_%i", jzn),
		";#eta;#it{p}_{T}^{Truth};#it{p}_{T}^{Reco}/#it{p}_{T}^{Truth}",
		nEtaForwardBinsCoarse, etaForwardMin, etaForwardMax,
		nPtTruthBins, ptTruthMin, ptTruthMax,
		nRPtBins, rPtMin, rPtMax);
    m_jznRPt[ jzn ]->Sumw2();
    v_hists.push_back( m_jznRPt[ jzn] );
    
    m_jznDeta[ jzn ] =
      new TH3D( Form("h_dEta_%i", jzn),
		";#eta;#it{p}_{T}^{Truth};#eta^{Reco}-#eta^{Truth}",
		nEtaForwardBinsCoarse, etaForwardMin, etaForwardMax,
		nPtTruthBins, ptTruthMin, ptTruthMax,
		nDAngleBins, dAngleMin, dAngleMax );
    m_jznDeta[ jzn ]->Sumw2();
    v_hists.push_back( m_jznDeta[ jzn] );
      
    m_jznDphi[ jzn ] =
      new TH3D( Form("h_dPhi_%i", jzn),
		";#eta;#it{p}_{T}^{Truth};#phi^{Reco}-#phi^{Truth}",
		nEtaForwardBinsCoarse, etaForwardMin, etaForwardMax,
		nPtTruthBins, ptTruthMin, ptTruthMax,
		nDAngleBins, dAngleMin, dAngleMax );
    m_jznDphi[ jzn ]->Sumw2();
    v_hists.push_back( m_jznDphi[ jzn] );
  } 
}

void DiJetAnalysisMC::ProcessEvents( int nEvents, int startEvent ){
  
}

// pair reco and truth jets with deltaR parameter 
void pairJets(  std::vector< TLorentzVector >& vR_jets,
		std::vector< TLorentzVector >& vT_jets,
		std::vector< JetPair >&  v_paired_jets ){ 

  // exit if either is empty
  if( !vT_jets.size() || !vR_jets.size() ){ return; }

  // clear vectors from previous time
  v_paired_jets.clear();
  
  for( auto& truthJet : vT_jets){
    // for each truth jet, need to find closest reco jet
    // set deltaRmin to 2*PI, find jet with smallest
    // deltaRmin less than Rmax and pair that to the truth jet
    double           deltaRmin    = 2*constants::PI;
    TLorentzVector* pairedRecoJet = NULL;
    for( auto& recoJet : vR_jets ){
      if( DeltaR( recoJet, truthJet ) <= deltaRmin ) {	
	pairedRecoJet = &recoJet;
	deltaRmin     = DeltaR( recoJet, truthJet );
      }
    }  // end loop over reco jets
    // just for safety. if there is at least
    // one reco jet, it should pair to all truth jets
    if( !pairedRecoJet ) continue;
    // we found one!
    // add the truth and reco jets together
    // also store their deltaR
    JetPair pairedJets( pairedRecoJet, &truthJet, deltaRmin);
    v_paired_jets.emplace_back( pairedJets );  
   } // end loop over truth jets
}


//---------------------------------
//            Plot Data
//---------------------------------

void DiJetAnalysisMC::LoadHistograms(){}

void DiJetAnalysisMC::PlotSpectra(){}

void DiJetAnalysisMC::PlotEtaPhi(){}

void DiJetAnalysisMC::PlotEtaPt(){}

/*
  int    nJESbins   = 40;
  double jesMin     = 0; double jesMax = 2;

  int    nJERbins   = 40;
  double jerMin     = 0; double jerMax = 2;

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
