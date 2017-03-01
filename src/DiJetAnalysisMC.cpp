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

void DiJetAnalysisMC::Initialize(){
  DiJetAnalysis::Initialize();

  v_usedJZN.push_back( 1 );
  m_jznFnameIn[ v_usedJZN.back() ] = "/home/yakov/Projects/atlas/data/jz1.root";
  m_jznSigma  [ v_usedJZN.back() ] = 6.7890E+07;
  m_jznEff    [ v_usedJZN.back() ]   = 2.8289E-03;
}

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
  // Triggers and Spectra
  int    nPtSpectBins = 100; 
  double ptSpectMin   = 0; double ptSpectMax  = nPtSpectBins;
  
  // Eta-Phi Maps
  int    nEtaBins = 100; 
  double etaMin   = constants::ETAMIN; double etaMax  = constants::ETAMAX;
  
  int    nPhiBins = 64; 
  double phiMin   = -constants::PI; double phiMax  = constants::PI; 
  
  double ptWidth  = 2;
  double ptMin    = 10;
  double ptMax    = 100;
  int    nPtBins  = (ptMax - ptMin)/ptWidth;
  
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
  
  for( auto jzn : v_usedJZN ){ // loop over JZ samples
    std::cout << "Making - JZ" << jzn << " histograms " << std::endl;

    m_jznEtaSpect[ jzn ] = 
      new TH2D( Form("h_spectEta_%i", jzn ), 
		";#eta_{Reco};#it{p}_{T}^{Reco} [GeV];dN/d#it{p}_{T}",
		nEtaForwardBinsCoarse, etaForwardMin, etaForwardMax,
		nPtSpectBins, ptSpectMin, ptSpectMax ) ;
    AddHistogram( m_jznEtaSpect[ jzn ] );
    
    m_jznEtaPhi[ jzn ] =
      new TH2D( Form("h_etaPhi_%i", jzn),
		";#eta_{Reco};#phi_{Reco}",
		nEtaBins, etaMin, etaMax,
		nPhiBins, phiMin, phiMax );
    AddHistogram( m_jznEtaPhi[ jzn ] );
      
    m_jznEtaPt[ jzn ] =
      new TH2D( Form("h_etaPt_%i", jzn),
		";#eta_{Reco};#it{p}_{T}^{Reco}",
		nEtaBins, etaMin, etaMax,
		nPtBins, ptMin, ptMax );
    AddHistogram( m_jznEtaPt[ jzn ] );
    
    m_jznRPt[ jzn ] =
      new TH3D( Form("h_rPt_%i", jzn),
		";#eta;#it{p}_{T}^{Truth};#it{p}_{T}^{Reco}/#it{p}_{T}^{Truth}",
		nEtaForwardBinsCoarse, etaForwardMin, etaForwardMax,
		nPtTruthBins, ptTruthMin, ptTruthMax,
		nRPtBins, rPtMin, rPtMax);
    AddHistogram( m_jznRPt[ jzn ] );

    m_jznDeta[ jzn ] =
      new TH3D( Form("h_dEta_%i", jzn),
		";#eta;#it{p}_{T}^{Truth};#eta^{Reco}-#eta^{Truth}",
		nEtaForwardBinsCoarse, etaForwardMin, etaForwardMax,
		nPtTruthBins, ptTruthMin, ptTruthMax,
		nDAngleBins, dAngleMin, dAngleMax );
    AddHistogram( m_jznDeta[ jzn ] );    
      
    m_jznDphi[ jzn ] =
      new TH3D( Form("h_dPhi_%i", jzn),
		";#eta;#it{p}_{T}^{Truth};#phi^{Reco}-#phi^{Truth}",
		nEtaForwardBinsCoarse, etaForwardMin, etaForwardMax,
		nPtTruthBins, ptTruthMin, ptTruthMax,
		nDAngleBins, dAngleMin, dAngleMax );
    AddHistogram( m_jznDphi[ jzn ] );
  } 
}

void DiJetAnalysisMC::ProcessEvents( int nEvents, int startEvent ){
  // collections and variables
  std::vector< TLorentzVector >    vT_jets;
  std::vector< TLorentzVector >* p_vT_jets = &vT_jets;

  std::vector< TLorentzVector >    vR_jets;
  std::vector< TLorentzVector >* p_vR_jets = &vR_jets;

  std::vector< bool >    v_isCleanJet;
  std::vector< bool >* p_v_isCleanJet = &v_isCleanJet;

  double rMax = 0.2;
  
  for( auto jzn : v_usedJZN ){
    
    std::cout << "fNameIn: " << m_jznFnameIn[ jzn ] << std::endl;
  
    m_fIn = TFile::Open( m_jznFnameIn[ jzn ].c_str() );
    m_tree = (TTree*) m_fIn->Get( "tree" );

    // Connect to tree
    m_tree->SetBranchAddress( "vT_jets"     , &p_vT_jets );
    m_tree->SetBranchAddress( "vR_C_jets"   , &p_vR_jets );
    m_tree->SetBranchAddress( "v_isCleanJet", &p_v_isCleanJet );

    // n events
    int maxEvents = m_tree->GetEntries();

    nEvents = nEvents > 0 ? nEvents : maxEvents;
    startEvent = startEvent < maxEvents ? startEvent : maxEvents - 1;
  
    int endEvent = startEvent + nEvents < maxEvents ? startEvent + nEvents : maxEvents;

    // event loop
    for( int ev = startEvent; ev < endEvent; ev++ ){
      m_tree->GetEntry( ev );

      ApplyCleaning ( vR_jets, v_isCleanJet );
      ApplyIsolation( 1.0, vR_jets );

      if( DoPrint(ev) ) {
	std::cout << "\nEvent : " << ev 
		  << "    has : " << vR_jets.size() << " reco jets"
		  << "    and : " << vT_jets.size() << " truth jets" 
		  << std::endl; 
      }

      std::vector< JetPair > v_paired_jets;
      PairJets( vR_jets, vT_jets, v_paired_jets );

	
      for( auto& vp : v_paired_jets ){
	/*
	if( DoPrint(ev) ){
	  std::cout << "-----Pair---- deltaR = " << vp.DeltaR() << std::endl;
	  vp.RecoJet()->Print();
	  vp.TruthJet()->Print();
	}
	*/
	  
	double etaTruth = vp.TruthJet()->Eta();
	double phiTruth = vp.TruthJet()->Phi();
	double  ptTruth = vp.TruthJet()->Pt()/1000.;
	  
	double etaReco = vp.RecoJet()->Eta();
	double phiReco = vp.RecoJet()->Phi();
	double  ptReco = vp.RecoJet()->Pt()/1000.;
	  
	  
	m_jznEtaSpect[ jzn ]->Fill( etaReco, ptReco );
	m_jznEtaPhi  [ jzn ]->Fill( etaReco, phiReco );
	m_jznEtaPt   [ jzn ]->Fill( etaReco, ptReco );
	  
	if( vp.DeltaR() <= rMax  ){
	  m_jznRPt [ jzn ]->Fill( etaTruth, ptTruth, ptReco/ptTruth );
	  m_jznDeta[ jzn ]->Fill( etaTruth, ptTruth, etaReco - etaTruth );
	  m_jznDphi[ jzn ]->Fill( etaTruth, ptTruth, phiReco - phiTruth );
	}
      } // end loop over pairs
      
    } // end loop over events
    std::cout << "DONE WITH JZ" << jzn << std::endl;
    m_fIn->Close();
  } // end loop over a JZ sample
}

// pair reco and truth jets with deltaR parameter 
void DiJetAnalysisMC::PairJets(  std::vector< TLorentzVector >& vR_jets,
				 std::vector< TLorentzVector >& vT_jets,
				 std::vector< JetPair >&  v_paired_jets ){ 
  // clear vectors from previous time
  v_paired_jets.clear();
  
  // exit if either is empty
  if( !vT_jets.size() || !vR_jets.size() ){ return; }
  
  for( auto& truthJet : vT_jets){
    // for each truth jet, need to find closest reco jet
    // set deltaRmin to something large, find jet with smallest
    // deltaRmin less than Rmax and pair that to the truth jet
    double              deltaRmin = 4*constants::ETAMAX;
    TLorentzVector* pairedRecoJet = NULL;
    for( auto& recoJet : vR_jets ){
      // "cleaned" jets had their px, py, pz set to 0
      // to avoid changing the vector size
      if( recoJet.Pt() == 0 ) continue;
      // if we have a "real" reco jet, continue
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
