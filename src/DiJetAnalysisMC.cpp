#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TROOT.h>
#include <TEnv.h>
#include <TCanvas.h>
#include <TRandom.h>

#include <iostream>
#include <sstream>
#include <cmath>
#include <boost/format.hpp>
#include <boost/filesystem.hpp>

#include "MyRoot.h"

#include "DiJetAnalysisMC.h"
#include "JetPair.h"

DiJetAnalysisMC::DiJetAnalysisMC() : m_mcType(0)
{}

DiJetAnalysisMC::DiJetAnalysisMC( bool isData, bool is_pPb,
				  int mcType )
  : DiJetAnalysis( isData, is_pPb), m_mcType( mcType )
{}

DiJetAnalysisMC::~DiJetAnalysisMC(){}

void DiJetAnalysisMC::Initialize(){
  m_labelOut = m_isData ? "_data" : "_mc" ;
  m_labelOut = m_is_pPb ? m_labelOut + "_pPb" : m_labelOut + "_pp";

  if( m_mcType == 0 ){
    m_labelOut += "_pythia8powheg";
    m_mcTypeLabel = "Pythia8+Powheg";
    
    m_vUsedJZN.push_back( 1 );
    m_mJznFnameIn[ m_vUsedJZN.back() ]
      = "/home/yakov/Projects/atlas/data/jz1_pythia8powheg.root";
    m_mJznSigma  [ m_vUsedJZN.back() ] = 1.2794E+08;
    m_mJznEff    [ m_vUsedJZN.back() ] = 1.5857E-03;
    m_mJznSumPowhegWeights[ m_vUsedJZN.back() ] = 3.428918e+08;
    
    m_vUsedJZN.push_back( 2 );
    m_mJznFnameIn[ m_vUsedJZN.back() ] 
      = "/home/yakov/Projects/atlas/data/jz2_pythia8powheg.root";
    m_mJznSigma  [ m_vUsedJZN.back() ] = 1.9648E+07;
    m_mJznEff    [ m_vUsedJZN.back() ] = 1.2948E-04;
    m_mJznSumPowhegWeights[ m_vUsedJZN.back() ] = 4.656080e+06;

    // Get Associated weights for powheg.
    // These need to be used at filling time (per event basis).
    TFile* weight_file = TFile::Open("/home/yakov/Projects/atlas/data/Powheg.reweight.root");
    m_hPowhegWeights = (TH3D*)weight_file->Get("h3_pT_y_phi_rw");
    m_hPowhegWeights->SetDirectory(0);
    weight_file->Close();
  } else if( m_mcType == 1 ){
    m_labelOut += "_pythia8";
    m_mcTypeLabel = "Pythia8";

    m_vUsedJZN.push_back( 0 );
    m_mJznFnameIn[ m_vUsedJZN.back() ]
      = "/home/yakov/Projects/atlas/data/jz0_pythia8.root";
    m_mJznSigma  [ m_vUsedJZN.back() ] = 6.7890E+07;
    m_mJznEff    [ m_vUsedJZN.back() ] = 8.6287E-02;
    m_mJznSumPowhegWeights[ m_vUsedJZN.back() ] = 1;
      
    m_vUsedJZN.push_back( 1 );
    m_mJznFnameIn[ m_vUsedJZN.back() ]
      = "/home/yakov/Projects/atlas/data/jz1_pythia8.root";
    m_mJznSigma  [ m_vUsedJZN.back() ] = 6.7890E+07*1.2;
    m_mJznEff    [ m_vUsedJZN.back() ] = 2.8289E-03;
    m_mJznSumPowhegWeights[ m_vUsedJZN.back() ] = 1;
    
    m_vUsedJZN.push_back( 2 );
    m_mJznFnameIn[ m_vUsedJZN.back() ] 
      = "/home/yakov/Projects/atlas/data/jz2_pythia8.root";
    m_mJznSigma  [ m_vUsedJZN.back() ] = 6.3997E+05;
    m_mJznEff    [ m_vUsedJZN.back() ] = 4.2714E-03;
    m_mJznSumPowhegWeights[ m_vUsedJZN.back() ] = 1;
  } else if( m_mcType == 2 ){
    m_labelOut += "_herwig";
    m_mcTypeLabel = "Herwig";

    m_vUsedJZN.push_back( 0 );
    m_mJznFnameIn[ m_vUsedJZN.back() ]
      = "/home/yakov/Projects/atlas/data/jz0_herwig.root";
    m_mJznSigma  [ m_vUsedJZN.back() ] = 69.2475E+06;
    m_mJznEff    [ m_vUsedJZN.back() ] = 5.713013E-04;
    //m_mJznSigma  [ m_vUsedJZN.back() ] = 3.0028E+02;
    //m_mJznEff    [ m_vUsedJZN.back() ] = 2.0899E-02;
    m_mJznSumPowhegWeights[ m_vUsedJZN.back() ] = 1;
    
    m_vUsedJZN.push_back( 1 );
    m_mJznFnameIn[ m_vUsedJZN.back() ]
      = "/home/yakov/Projects/atlas/data/jz1_herwig.root";
    m_mJznSigma  [ m_vUsedJZN.back() ] = 69.2475E+06;
    m_mJznEff    [ m_vUsedJZN.back() ] = 1.5088E-03;
    //m_mJznSigma  [ m_vUsedJZN.back() ] = 3.1313E+02;
    //m_mJznEff    [ m_vUsedJZN.back() ] = 1.5062E-03;
    m_mJznSumPowhegWeights[ m_vUsedJZN.back() ] = 1;
       
    m_vUsedJZN.push_back( 2 );
    m_mJznFnameIn[ m_vUsedJZN.back() ] 
      = "/home/yakov/Projects/atlas/data/jz2_herwig.root";
    m_mJznSigma  [ m_vUsedJZN.back() ] = 4.2389E+05;
    m_mJznEff    [ m_vUsedJZN.back() ] = 3.21E-03;
    //m_mJznSigma  [ m_vUsedJZN.back() ] = 4.2389E+05;
    //m_mJznEff    [ m_vUsedJZN.back() ] = 3.2383E-03;
    m_mJznSumPowhegWeights[ m_vUsedJZN.back() ] = 1;
  }  

  auto checkWriteDir = []( const char* c_dirOut ){
    boost::filesystem::path dir( c_dirOut );  
    if(!(boost::filesystem::exists(dir))){
      std::cout<< c_dirOut << " doesn't Exist."<<std::endl;
      if (boost::filesystem::create_directory(dir))
	std::cout << "....Successfully Created !" << std::endl;
    }
  };

  // Check if the directories exist.
  // If they don't, create them
  m_dirOut   = "output";
  checkWriteDir( m_dirOut.c_str() );
  m_dirOut   += "/output" + m_labelOut;
  checkWriteDir( m_dirOut.c_str() );

  m_rootFname = m_dirOut + "/myOut" + m_labelOut + ".root";
  
  std::cout << "fNameIn/Out: " << m_rootFname << std::endl;

  //========== Cuts and triggers =============
  
  m_dRmax = 0.2;
  
  // calculate sum of sigma and eff
  m_sumSigmaEff = 0;
  for( auto jzn : m_vUsedJZN ){
    m_sumSigmaEff += ( m_mJznSigma[ jzn ] * m_mJznEff[ jzn ] );
  }

  std::cout << "SumSigmaEff = " << m_sumSigmaEff << std::endl;
  
  for( auto jzn : m_vUsedJZN ){ // loop over JZ samples
    m_fIn = TFile::Open( m_mJznFnameIn[ jzn ].c_str() );
    m_tree = (TTree*) m_fIn->Get( "tree" );

    int nEventsTotal = m_tree->GetEntries();
    double     sigma = m_mJznSigma[ jzn ];      
    double       eff = m_mJznEff  [ jzn ];
  
    m_mJznNev    [ jzn ] = nEventsTotal;
    m_mJznWeights[ jzn ] =
      (1./nEventsTotal) * (1./m_sumSigmaEff) * (sigma * eff);
  
    std::cout << "JZ" << jzn << "   weight = "
	      << m_mJznWeights [jzn] << std::endl;
    m_fIn->Close();
  }

  //========== Set Histogram Binning =============
  // Common for all analysis
  
  // Triggers and Spectra
  m_nPtSpectBins = 50;  
  m_ptSpectMin   = 10;
  m_ptSpectMax   = 2 * m_nPtSpectBins + m_ptSpectMin;
  
  // Eta-Phi Maps
  m_nEtaBins = 100; 
  m_etaMin   = constants::ETAMIN;
  m_etaMax   = constants::ETAMAX;
  
  m_nPhiBins = 64; 
  m_phiMin   = -constants::PI;
  m_phiMax   = constants::PI; 
  
  m_ptWidth  = 2;
  m_ptMin    = 10;
  m_ptMax    = 100;
  m_nPtBins  = (m_ptMax - m_ptMin)/m_ptWidth;
  
  m_nEtaForwardBinsFine   = 12;
  m_nEtaForwardBinsCoarse = 3;
  m_etaForwardMin   = -constants::FETAMAX;
  m_etaForwardMax   = -constants::FETAMIN;

  // truth bins
  m_ptTruthWidth  = 5;
  m_ptTruthMin    = 10;
  m_ptTruthMax    = 100;
  m_nPtTruthBins  =
    (m_ptTruthMax - m_ptTruthMin) / m_ptTruthWidth;

  // Jes Jer
  m_nRPtBins   = 60;
  m_rPtMin     = 0;  m_rPtMax = 3;

  // angular bins
  m_nDAngleBins  = 50;
  m_dAngleMin    = -0.5;  m_dAngleMax = 0.5;
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

  std::string cfNameOut = m_dirOut + "/c_myOut" + m_labelOut + ".root";
  m_fOut = new TFile( cfNameOut.c_str(),"RECREATE");
  
  // check if we have one
  PlotSpectra( m_mJznEtaSpect );
  
  PlotEtaPhiPtMap( m_mJznEtaPhi );
  PlotEtaPhiPtMap( m_mJznEtaPt  );

  PlotVsEtaPt( m_mJznRpt , m_mJznNentries, "rPt");
  PlotVsEtaPt( m_mJznDeta, m_mJznNentries, "dEta");
  PlotVsEtaPt( m_mJznDphi, m_mJznNentries, "dPhi");

  std::cout << "DONE! Closing " << cfNameOut << std::endl;
  m_fOut->Close();
  std::cout << "......Closed  " << cfNameOut << std::endl;
}

//---------------------------------
//            Read Data
//---------------------------------

void DiJetAnalysisMC::SetupHistograms(){
  
  for( auto jzn : m_vUsedJZN ){ // loop over JZ samples
    std::cout << "Making - JZ" << jzn << " histograms " << std::endl;

    m_mJznEtaSpect[ jzn ] = 
      new TH2D( Form("h_etaSpect_jz%i", jzn ), 
		";#eta_{Reco};#it{p}_{T}^{Reco} [GeV];dN/d#it{p}_{T}",
		m_nEtaForwardBinsCoarse,
		m_etaForwardMin, m_etaForwardMax,
		m_nPtSpectBins,
		m_ptSpectMin, m_ptSpectMax ) ;
    AddHistogram( m_mJznEtaSpect[ jzn ] );
    
    m_mJznEtaPhi[ jzn ] =
      new TH2D( Form("h_etaPhi_jz%i", jzn),
		";#eta_{Reco};#phi_{Reco}",
		m_nEtaBins, m_etaMin, m_etaMax,
		m_nPhiBins, m_phiMin, m_phiMax );
    AddHistogram( m_mJznEtaPhi[ jzn ] );
      
    m_mJznEtaPt[ jzn ] =
      new TH2D( Form("h_etaPt_jz%i", jzn),
		";#eta_{Reco};#it{p}_{T}^{Reco}",
		m_nEtaBins, m_etaMin, m_etaMax,
		m_nPtBins , m_ptMin  , m_ptMax );
    AddHistogram( m_mJznEtaPt[ jzn ] );
    
    m_mJznRpt[ jzn ] =
      new TH3D( Form("h_rPt_jz%i", jzn),
		";#eta^{Truth};#it{p}_{T}^{Truth};#it{p}_{T}^{Reco}/#it{p}_{T}^{Truth}",
		m_nEtaForwardBinsFine,
		m_etaForwardMin, m_etaForwardMax,
		m_nPtTruthBins,
		m_ptTruthMin, m_ptTruthMax,
		m_nRPtBins,
		m_rPtMin, m_rPtMax);
    AddHistogram( m_mJznRpt[ jzn ] );
  
    m_mJznDeta[ jzn ] =
      new TH3D( Form("h_dEta_jz%i", jzn),
		";#eta^{Truth};#it{p}_{T}^{Truth};#eta^{Reco}-#eta^{Truth}",
		m_nEtaForwardBinsFine,
		m_etaForwardMin, m_etaForwardMax,
		m_nPtTruthBins,
		m_ptTruthMin, m_ptTruthMax,
		m_nDAngleBins,
		m_dAngleMin, m_dAngleMax );
    AddHistogram( m_mJznDeta[ jzn ] );    
  
    m_mJznDphi[ jzn ] =
      new TH3D( Form("h_dPhi_jz%i", jzn),
		";#eta^{Truth};#it{p}_{T}^{Truth};#phi^{Reco}-#phi^{Truth}",
		m_nEtaForwardBinsFine,
		m_etaForwardMin, m_etaForwardMax,
		m_nPtTruthBins,
		m_ptTruthMin, m_ptTruthMax,
		m_nDAngleBins,
		m_dAngleMin, m_dAngleMax );
    AddHistogram( m_mJznDphi[ jzn ] );

    m_mJznNentries[ jzn ] =
      new TH2D( Form("h_nEntries_jz%i", jzn),
		";#eta^{Truth};#it{p}_{T}^{Truth}",
		m_nEtaForwardBinsFine,
		m_etaForwardMin, m_etaForwardMax,
		m_nPtTruthBins,
		m_ptTruthMin, m_ptTruthMax );
    AddHistogram( m_mJznNentries[ jzn ] );
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

  for( auto jzn : m_vUsedJZN ){
    
    std::cout << "fNameIn: " << m_mJznFnameIn[ jzn ] << std::endl;
  
    m_fIn = TFile::Open( m_mJznFnameIn[ jzn ].c_str() );
    m_tree = (TTree*) m_fIn->Get( "tree" );

    // Connect to tree
    m_tree->SetBranchAddress( "vT_jets"     , &p_vT_jets );
    m_tree->SetBranchAddress( "vR_C_jets"   , &p_vR_jets );
    m_tree->SetBranchAddress( "v_isCleanJet", &p_v_isCleanJet );

    // n events
    int nEventsTotal = m_tree->GetEntries();

    nEvents = nEvents > 0 ? nEvents : nEventsTotal;
    startEvent = startEvent < nEventsTotal ? startEvent : nEventsTotal - 1;
  
    int endEvent = startEvent + nEvents < nEventsTotal ?
					  startEvent + nEvents : nEventsTotal;

    // event loop
    for( int ev = startEvent; ev < endEvent; ev++ ){
      m_tree->GetEntry( ev );

      ApplyCleaning ( vR_jets, v_isCleanJet );
      ApplyIsolation( 1.0, vR_jets );

      if( AnalysisTools::DoPrint(ev) ) {
	std::cout << "\nEvent : " << ev 
		  << "    has : " << vR_jets.size() << " reco jets"
		  << "    and : " << vT_jets.size() << " truth jets" 
		  << std::endl; 
      }

      std::vector< JetPair > v_paired_jets;
      PairJets( vR_jets, vT_jets, v_paired_jets );

	
      for( auto& vp : v_paired_jets ){	  
	double jetEtaTruth = vp.TruthJet()->Eta();
	double jetPhiTruth = vp.TruthJet()->Phi();
	double  jetPtTruth = vp.TruthJet()->Pt()/1000.;
	  
	double jetEtaReco = vp.RecoJet()->Eta();
	double jetPhiReco = vp.RecoJet()->Phi();
	double  jetPtReco = vp.RecoJet()->Pt()/1000.;

	double weight = m_mcType == 0 ?
	  GetJetWeight( jetEtaTruth, jetPhiTruth, jetPtTruth ) : 1;
	
	m_mJznEtaSpect[ jzn ]->Fill( jetEtaReco, jetPtReco, weight);
	m_mJznEtaPhi  [ jzn ]->Fill( jetEtaReco, jetPhiReco, weight);
	m_mJznEtaPt   [ jzn ]->Fill( jetEtaReco, jetPtReco , weight);

	// convert positive eta to negative because
	// in pp it doesnt matter.
	// our histos run in negative eta
	// (due to pPb configuration)
	// the labels will be taken care of so it is ok
	double jetEtaRecoAdj  = AdjustEtaForPP( jetEtaReco  );
	double jetEtaTruthAdj = AdjustEtaForPP( jetEtaTruth );
	
	if( vp.DeltaR() <= m_dRmax ){
	  m_mJznRpt   [ jzn ]->
	    Fill( jetEtaTruthAdj, jetPtTruth, jetPtReco/jetPtTruth    , weight);
	  m_mJznDeta  [ jzn ]->
	    Fill( jetEtaTruthAdj, jetPtTruth, jetEtaReco - jetEtaTruth, weight);
	  m_mJznDphi  [ jzn ]->
	    Fill( jetEtaTruthAdj, jetPtTruth, jetPhiReco - jetPhiTruth, weight );
	  m_mJznNentries[ jzn ]->
	    Fill( jetEtaTruthAdj, jetPtTruth, weight );
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
      if( AnalysisTools::DeltaR( recoJet, truthJet ) <= deltaRmin ) {	
	pairedRecoJet = &recoJet;
	deltaRmin     = AnalysisTools::DeltaR( recoJet, truthJet );
      }
    }  // end loop over reco jets
    // just for safety. if there is at least
    // one reco jet, it should pair to all truth jets
    if( !pairedRecoJet ) continue;
    // we found one!
    // add the truth and reco jets together
    // also store their deltaR
    v_paired_jets.emplace_back( JetPair( pairedRecoJet,
					 &truthJet,
					 deltaRmin) );  
  } // end loop over truth jets
}


//---------------------------------
//            Plot Data
//---------------------------------

void DiJetAnalysisMC::LoadHistograms(){
  m_fIn = TFile::Open( m_rootFname.c_str() ); 

  for( auto& jzn : m_vUsedJZN ){
    m_mJznEtaSpect [ jzn ] =
      (TH2D*)m_fIn->Get( Form("h_etaSpect_jz%i" , jzn ) );
    m_mJznEtaSpect [ jzn ]->SetDirectory(0);
    m_mJznEtaPhi   [ jzn ] =
      (TH2D*)m_fIn->Get( Form("h_etaPhi_jz%i"   , jzn ) );
    m_mJznEtaPhi   [ jzn ]->SetDirectory(0);
    m_mJznEtaPt    [ jzn ] =
      (TH2D*)m_fIn->Get( Form("h_etaPt_jz%i"    , jzn ) );
    m_mJznEtaPt    [ jzn ]->SetDirectory(0);

    m_mJznRpt     [ jzn ] =
      (TH3D*)m_fIn->Get( Form("h_rPt_jz%i"    , jzn ) );
    m_mJznRpt     [ jzn ]->SetDirectory(0);
    m_mJznDeta    [ jzn ]   =
      (TH3D*)m_fIn->Get( Form("h_dEta_jz%i"    , jzn ) );
    m_mJznDeta    [ jzn ]->SetDirectory(0);
    m_mJznDphi    [ jzn ] =
      (TH3D*)m_fIn->Get( Form("h_dPhi_jz%i"    , jzn ) );
    m_mJznDphi    [ jzn ]->SetDirectory(0);

    m_mJznNentries [ jzn ] =
      (TH2D*)m_fIn->Get( Form("h_nEntries_jz%i"    , jzn ) );
    m_mJznNentries [ jzn ]->SetDirectory(0);

  }
  m_fIn->Close();
}

void DiJetAnalysisMC::PlotSpectra( std::map< int, TH2* >& mJznHIN ){
  std::map< int, TH1* > mJznH;

  // check if we have one
  if( !mJznHIN.size() ){ return; }

  for( int xBin = 1;
       xBin <= mJznHIN.begin()->second->GetNbinsX();
       xBin++ ){

    TCanvas c_spect("c_spect","c_spect",800,600);
    c_spect.SetLogy();

    TLegend l_spect(0.15, 0.15, 0.46, 0.28);
    StyleTools::SetLegendStyle( &l_spect, StyleTools::lSS );
    l_spect.SetFillStyle(0);

    int style = 1;
    double max = -1;

        // should all be the same
    double etaMin = mJznHIN.begin()->second->
      GetXaxis()->GetBinLowEdge( xBin );
    double etaMax = mJznHIN.begin()->second->
      GetXaxis()->GetBinUpEdge( xBin );
  
    for( auto& jznHIN : mJznHIN ){ // loop over JZN
      int jzn = jznHIN.first;
    
      TH1* hs =
	jznHIN.second->
	ProjectionY( Form("h_spect_jz%i_%2.0f_Eta_%2.0f",
			  jzn,
			  10*std::abs(etaMin),
			  10*std::abs(etaMax) ),
		     xBin, xBin );
      mJznH[ jzn ] = hs ; 
    
      l_spect.AddEntry(  hs, Form("JZ%i", jzn) );
      StyleTools::SetHStyle( hs, style++, StyleTools::hSS);
      hs->Draw("epsame");
      if( max < hs->GetMaximum() ){ max = hs->GetMaximum(); }

      double power = log10(max);
      power = std::ceil(power); 
      max = pow( 10, power );
      for( auto& jznH : mJznH ){
	jznH.second->SetMaximum( max );
	jznH.second->SetMinimum( 0.1 );
      }
    } // end loop over JZN

    int  nPtBins = mJznHIN.begin()->second->GetNbinsY();
    double ptMin = mJznHIN.begin()->second->GetYaxis()->GetXmin();
    double ptMax = mJznHIN.begin()->second->GetYaxis()->GetXmax();
  
    TH1* h_spect = new TH1D( Form("h_spect_%2.0f_Eta_%2.0f_final",
				  10*std::abs(etaMin),
				  10*std::abs(etaMax) ),
			     mJznHIN.begin()->second->GetTitle(),
			     nPtBins, ptMin, ptMax );
    StyleTools::SetHStyle( h_spect, 0, StyleTools::hSS);
    l_spect.AddEntry( h_spect, "Total" );
  
    CombineJZN( h_spect, mJznH );
    
    h_spect->Draw("epsame");

    
    l_spect.Draw();
    DrawTools::DrawAtlasInternalMCRight( 0, 0, StyleTools::lSS, m_mcTypeLabel ); 
    DrawTools::DrawLeftLatex( 0.5, 0.4,
			      GetEtaLabel( etaMin, etaMax ).c_str(),
			      StyleTools::lSS, 1 );

    c_spect.SaveAs( Form("%s/spectra_%2.0f_Eta_%2.0f%s.pdf",
			 m_dirOut.c_str(),
			 std::abs(etaMin)*10,
			 std::abs(etaMax)*10,
			 m_labelOut.c_str() ) );
    c_spect.SaveAs( Form("%s/spectra_%2.0f_Eta_%2.0f%s.png",
			 m_dirOut.c_str(),
			 std::abs(etaMin)*10,
			 std::abs(etaMax)*10,
			 m_labelOut.c_str() ) );

    
    c_spect.Write( Form("c_spectra_%2.0f_Eta_%2.0f%s",
			std::abs(etaMin)*10,
			std::abs(etaMax)*10,
			m_labelOut.c_str() ) );

  } // end loop over eta
}

void DiJetAnalysisMC::PlotEtaPhiPtMap( std::map< int, TH2* >& mJznHIN ){
  TCanvas c_map("c_map","c_map",800,600);

  for( auto& jznHIN : mJznHIN ){
    jznHIN.second->Draw("col");
    StyleTools::SetHStyle( jznHIN.second, 0, StyleTools::hSS);
    DrawTools::DrawAtlasInternalMCLeft( 0, -0.55, StyleTools::lSS, m_mcTypeLabel );  
    c_map.SaveAs( Form("%s/%s%s.pdf", 
		       m_dirOut.c_str(),
		       jznHIN.second->GetName(),
		       m_labelOut.c_str() ) );
    c_map.SaveAs( Form("%s/%s%s.png", 
		       m_dirOut.c_str(),
		       jznHIN.second->GetName(),
		       m_labelOut.c_str() ) );
  }
}

void DiJetAnalysisMC::PlotVsEtaPt( std::map< int, TH3* >& mJznHIN,
				   std::map< int, TH2* >& mJznNIN,
				   const std::string& type ){
  std::vector< TH1* > vMeans;
  std::vector< TH1* > vSigmas;

  std::vector< TH1* > vMeansFinal;
  std::vector< TH1* > vSigmasFinal;

  // check if we have one
  if( !mJznHIN.size() ){ return; }
  
  for( int xBin = 1;
       xBin <= mJznHIN.begin()->second->GetNbinsX();
       xBin++ ){
    
    // should all be the same
    double etaMin = mJznHIN.begin()->second->
      GetXaxis()->GetBinLowEdge( xBin );
    double etaMax = mJznHIN.begin()->second->
      GetXaxis()->GetBinUpEdge( xBin );
  
    int  nPtBins = mJznHIN.begin()->second->GetNbinsY();
    double ptMin = mJznHIN.begin()->second->GetYaxis()->GetXmin();
    double ptMax = mJznHIN.begin()->second->GetYaxis()->GetXmax();

    std::string yTitleMean;
    std::string yTitleSigma;

    if( !type.compare("rPt") ){
      yTitleMean  = "#it{p}_{T}^{Reco}/#it{p}_{T}^{Truth}";
      yTitleSigma = "#sigma" + yTitleMean;
    } else if( !type.compare("dEta") ){
      yTitleMean  = "#Delta#eta";
      yTitleSigma = "#sigma" + yTitleMean;
    } else if( !type.compare("dPhi") ){
      yTitleMean  = "#Delta#phi";
      yTitleSigma = "#sigma" + yTitleMean;
    } 

    std::map< int, TH1* > mJznMean;
    std::map< int, TH1* > mJznSigma;
    std::map< int, TH1* > mJznN;

    int style = 1;
    
    for( auto& jznHIN : mJznHIN ){   // loop over JZN
      int jzn = jznHIN.first;
    
      TH1* hMean = new TH1D( Form("h_%s_mean_jz%i_%2.0f_Eta_%2.0f",
				  type.c_str(),
				  jzn,
				  10*std::abs(etaMin),
				  10*std::abs(etaMax) ),
			     Form("%s;%s;%s",
				  GetEtaLabel( etaMin, etaMax ).c_str(),
				  jznHIN.second->GetYaxis()->GetTitle(),
				  yTitleMean.c_str() ),
			     nPtBins, ptMin, ptMax );    
      StyleTools::SetHStyle( hMean, style, StyleTools::hSS);
      mJznMean[ jzn ] = hMean;
      vMeans.push_back( hMean );
    
      TH1* hSigma = new TH1D( Form("h_%s_sigma_jz%i_%2.0f_Eta_%2.0f",
				   type.c_str(),
				   jzn,
				   10*std::abs(etaMin),
				   10*std::abs(etaMax) ),
			      Form("%s;%s;%s",
				   GetEtaLabel( etaMin, etaMax ).c_str(),
				   jznHIN.second->GetYaxis()->GetTitle(),
				   yTitleSigma.c_str() ),
			      nPtBins, ptMin, ptMax );    
      StyleTools::SetHStyle( hSigma, style, StyleTools::hSS);
      mJznSigma[ jzn ] = hSigma;
      vSigmas.push_back( hSigma );

      TH1* hN = mJznNIN[ jzn ]->
	ProjectionY( Form("h_%s_N_jz%i_%2.0f_Eta_%2.0f",
			  type.c_str(),
			  jzn,
			  10*std::abs(etaMin),
			  10*std::abs(etaMax) ),
		     xBin, xBin);
      mJznN[ jzn ] = hN;
    
      ProjectEtaPtAndFit( jznHIN.second, hMean, hSigma, xBin, xBin );
      
      // increment style
      style++;

    } // end loop over JZN

    // save eta ranges in title to draw all etas on one canvas later
    // this is for the legends.... 
    TH1* hMeanFinal = new TH1D( Form("h_%s_mean_%2.0f_Eta_%2.0f_final",
				     type.c_str(),
				     10*std::abs(etaMin),
				     10*std::abs(etaMax) ),
				Form("%s;%s;%s",
				     GetEtaLabel( etaMin, etaMax).c_str(),
				     mJznHIN.begin()->second->
				     GetYaxis()->GetTitle(),
				     yTitleSigma.c_str() ),
				nPtBins, ptMin, ptMax );
    StyleTools::SetHStyle( hMeanFinal, 0, StyleTools::hSS);
    vMeansFinal.push_back( hMeanFinal );
  
    TH1* hSigmaFinal = new TH1D( Form("h_%s_sigma_%2.0f_Eta_%2.0f_final",
				      type.c_str(),
				      10*std::abs(etaMin),
				      10*std::abs(etaMax) ),
				 Form("%s;%s;%s",
				      GetEtaLabel( etaMin, etaMax).c_str(),
				      mJznHIN.begin()->second->
				      GetYaxis()->GetTitle(),
				      yTitleSigma.c_str() ),
				 nPtBins, ptMin, ptMax );
    StyleTools::SetHStyle( hSigmaFinal, 0, StyleTools::hSS);
    vSigmasFinal.push_back( hSigmaFinal );
    
    CombineJZN( hMeanFinal , mJznMean , mJznN );
    CombineJZN( hSigmaFinal, mJznSigma, mJznN );
    
    // Draw the canvases for mean and sigma
    DrawCanvas( mJznMean,  hMeanFinal,  etaMin, etaMax, type, true ); 
    DrawCanvas( mJznSigma, hSigmaFinal, etaMin, etaMax, type, false ); 
  } // end loop over eta


  DrawCanvas( vMeansFinal , type, true  );
  DrawCanvas( vSigmasFinal, type, false );
}

//---------------------------
//          Tools 
//---------------------------
void DiJetAnalysisMC::ProjectEtaPtAndFit( TH3* h3,
					  TH1* h1Mean, TH1* h1Sigma,
					  int etaBinLow, int etaBinUp ){
  TCanvas c_proj("c_proj","c_proj",800,600);

  double etaMin = h3->GetXaxis()->GetBinLowEdge( etaBinLow );
  double etaMax = h3->GetXaxis()->GetBinUpEdge( etaBinUp );

  int  nPtBins  = h3->GetNbinsY();

  std::vector< TH1* > v_hProj;
  std::vector< TF1* > v_fit;
  
  // loop over ptBins
  for( int ptBin = 1; ptBin <= nPtBins; ptBin++ ){
    double ptMin = h3->GetYaxis()->GetBinLowEdge(ptBin);
    double ptMax = h3->GetYaxis()->GetBinUpEdge(ptBin);
      
    TH1* hProj = h3->
      ProjectionZ( Form("%s_%2.0f_Eta_%2.0f_%2.0f.Pt.%2.0f",
			h3->GetName(),
			10*std::abs(etaMin),
			10*std::abs(etaMax),
			ptMin, ptMax ),
		   etaBinLow, etaBinUp, ptBin, ptBin );
    StyleTools::SetHStyle( hProj, 0, StyleTools::hSS);
    v_hProj.push_back( hProj );
    hProj->SetTitle("");
    
    TF1* fit  = new TF1( Form("f_%s_%2.0f_Eta_%2.0f_%2.0f.Pt.%2.0f",
			      h3->GetName(),
			      10*std::abs(etaMin),
			      10*std::abs(etaMax),
			      ptMin, ptMax ),
			 "gaus(0)" );
    StyleTools::SetHStyle( fit, 0, StyleTools::hSS);
    v_fit.push_back( fit );

    if( hProj->GetEntries() < 5 ){ continue; }
    
    AnalysisTools::FitGaussian( hProj, fit );
        
    h1Mean->SetBinContent( ptBin, fit->GetParameter(1) );
    h1Mean->SetBinError  ( ptBin, fit->GetParError (1) );

    h1Sigma->SetBinContent ( ptBin, fit->GetParameter(2) );
    h1Sigma->SetBinError   ( ptBin, fit->GetParError (2) );
    
    hProj->Draw();
    fit->Draw("same");

    DrawTools::DrawAtlasInternalMCRight( 0, 0, StyleTools::lSS, m_mcTypeLabel ); 
    DrawTools::DrawLeftLatex( 0.5, 0.8,
			      GetEtaLabel( etaMin, etaMax ).c_str(),
			      StyleTools::lSS, 1 );
    DrawTools::DrawLeftLatex( 0.5, 0.73,
			      Form("%3.0f<%s<%3.1f",
				   ptMin,
				   h1Mean->GetXaxis()->GetTitle(),
				   ptMax ),
			      StyleTools::lSS, 1 );


    c_proj.SaveAs( Form("%s/fits/%s_%2.0f_Eta_%2.0f_%2.0f.Pt.%2.0f%s.pdf",
			m_dirOut.c_str(),
			h3->GetName(),
			std::abs(etaMin)*10,
			std::abs(etaMax)*10,
			ptMin, ptMax,
			m_labelOut.c_str() ) );


    c_proj.Write( Form("c_%s_%2.0f_Eta_%2.0f_%2.0f.Pt.%2.0f%s",
		       h3->GetName(),
		       std::abs(etaMin)*10,
		       std::abs(etaMax)*10,
		       ptMin, ptMax,
		       m_labelOut.c_str() ) ); 
  } // end loop over ptBins
}

void DiJetAnalysisMC::CombineJZN( TH1* h_res,
				  std::map< int, TH1* >& mJznHIN){
  for( auto& jznH : mJznHIN ){
    int jzn      = jznH.first;
    double scale = m_mJznEff[ jzn ] * m_mJznSigma[ jzn ] /
      m_mJznSumPowhegWeights[ jzn ];

    std:: cout << "jz" << jzn << "   scale: " << scale << std::endl;
    h_res->Add( jznH.second, scale );
  }
  h_res->Scale( 1./m_sumSigmaEff );
}

void DiJetAnalysisMC::CombineJZN( TH1* h_res,
				  std::map< int, TH1*>& mJznVIN,
				  std::map< int, TH1*>& mJznNIN ){
  for( int xBin = 1; xBin <= h_res->GetNbinsX(); xBin++ ){

    double valTot   = 0; double valErrorTot   = 0;
    double denomTot = 0;
    
    for( auto jznV : mJznVIN ){
      int jzn = jznV.first;

      double valueBin    = mJznVIN[ jzn ]->GetBinContent( xBin );
      double valueBinErr = mJznVIN[ jzn ]->GetBinError( xBin );
      double nEntriesBin = mJznNIN[ jzn ]->GetBinContent( xBin );
      double weight      = m_mJznWeights[ jzn ] / m_mJznSumPowhegWeights[ jzn ];
      
      double val    = weight * nEntriesBin * valueBin;
      double valErr = weight * nEntriesBin * valueBinErr;
      valTot       += val;
      valErrorTot  += valErr * valErr;

      double denom  = weight * nEntriesBin;
      denomTot     += denom;
    }

    // dont want to set bins where val divided by 0
    if( AnalysisTools::EpsilonEqual( denomTot, 0. ) )
      { continue; };
    
    valTot      /= denomTot;
    valErrorTot /= denomTot * denomTot;
    
    h_res->SetBinContent( xBin, valTot );
    h_res->SetBinError  ( xBin, valErrorTot );
  }
}

double DiJetAnalysisMC::GetJetWeight( double eta, double phi, double pt ){
  int xb=m_hPowhegWeights->GetXaxis()->FindBin(pt);
  int yb=m_hPowhegWeights->GetYaxis()->FindBin(eta);
  int zb=m_hPowhegWeights->GetZaxis()->FindBin(phi);
  float jet_weight=m_hPowhegWeights->GetBinContent(xb,yb,zb);

  return jet_weight;
}


//---------------------------
//          Drawing 
//---------------------------
void DiJetAnalysisMC::DrawCanvas( std::map< int, TH1* >& mJznHIN, TH1* hFinal,
				  double etaMin, double etaMax,
				  const std::string& type, bool isMean ){
  TCanvas c("c","c",800,600);
  
  TLegend leg(0.68, 0.64, 0.99, 0.77);
  StyleTools::SetLegendStyle( &leg, StyleTools::lSS );
  leg.SetFillStyle(0);

  for( auto& jznHIN : mJznHIN ){
    jznHIN.second->SetTitle("");
    jznHIN.second->Draw("epsame");
    SetMinMax( jznHIN.second, type, isMean );
    leg.AddEntry( jznHIN.second, Form("JZ%i", jznHIN.first ) );
  }
  hFinal->Draw("same");
  leg.AddEntry( hFinal, "Total" );
  
  leg.Draw();
  
  DrawTools::DrawAtlasInternalMCRight( 0, 0, StyleTools::lSS, m_mcTypeLabel ); 
  DrawTools::DrawLeftLatex( 0.5, 0.8,
			    GetEtaLabel( etaMin, etaMax ).c_str(),
			    StyleTools::lSS, 1 );
  
  std::string meanOrSigma = isMean ? "mean" : "sigma";
  double y0 = GetLineHeight( type );
  
  double xMin = mJznHIN.begin()->second->GetXaxis()->GetXmin();
  double xMax = mJznHIN.begin()->second->GetXaxis()->GetXmax();
  
  TLine line( xMin, y0, xMax, y0);
  line.Draw();
    
  c.SaveAs( Form("%s/%s_%s_%2.0f_Eta_%2.0f%s.pdf",
		 m_dirOut.c_str(),
		 type.c_str(),
		 meanOrSigma.c_str(),
		 std::abs(etaMin)*10,
		 std::abs(etaMax)*10,
		 m_labelOut.c_str() ) );
  c.SaveAs( Form("%s/%s_%s_%2.0f_Eta_%2.0f%s.png",
		 m_dirOut.c_str(),
		 type.c_str(),
		 meanOrSigma.c_str(),
		 std::abs(etaMin)*10,
		 std::abs(etaMax)*10,
		 m_labelOut.c_str() ) );

  c.Write( Form("c_%s_%s_%2.0f_Eta_%2.0f%s",
		type.c_str(),
		meanOrSigma.c_str(),
		std::abs(etaMin)*10,
		std::abs(etaMax)*10,
		m_labelOut.c_str()) );  
}

void DiJetAnalysisMC::DrawCanvas( std::vector< TH1* >& vHIN,
				  const std::string& type, bool isMean ){
  TCanvas c("c","c",800,600);
  
  TLegend leg(0.68, 0.64, 0.99, 0.77);
  StyleTools::SetLegendStyle( &leg, StyleTools::lSS );
  leg.SetFillStyle(0);

  int style = 0;

  // plot every four on canvas
  int dEta = vHIN.size()/3; // plot every 3
  for( unsigned int etaRange = 0;
       etaRange < vHIN.size();
       etaRange += dEta){
    StyleTools::SetHStyle( vHIN[ etaRange], style++, StyleTools::hSS );
    leg.AddEntry( vHIN[ etaRange ], vHIN[ etaRange ]->GetTitle() );
    vHIN[ etaRange ]->SetTitle("");
    vHIN[ etaRange ]->Draw("epsame");
    SetMinMax( vHIN[ etaRange ], type, isMean );
  }

  leg.Draw();
  
  std::string meanOrSigma = isMean ? "mean" : "sigma";
  double y0 = GetLineHeight( type );
  
  double xMin = vHIN.front()->GetXaxis()->GetXmin();
  double xMax = vHIN.front()->GetXaxis()->GetXmax();
  
  TLine line( xMin, y0, xMax, y0);
  line.Draw();

  c.SaveAs( Form("%s/%s_%s%s.pdf",
		 m_dirOut.c_str(),
		 type.c_str(),
		 meanOrSigma.c_str(),
		 m_labelOut.c_str() ) );
  c.SaveAs( Form("%s/%s_%s%s.png",
		 m_dirOut.c_str(),
		 type.c_str(),
		 meanOrSigma.c_str(),
		 m_labelOut.c_str() ) );

  c.Write( Form("c_%s_%s%s",
		type.c_str(),
		meanOrSigma.c_str(),
		m_labelOut.c_str()) );  

}

//===== MinMax and line drawing =====

void DiJetAnalysisMC::
SetMinMax( TH1* h1, const std::string& type, bool isMean ){
  if( !type.compare("rPt") ){ // JES JER
    if( isMean ){
      h1->SetMinimum(0.75);
      h1->SetMaximum(1.25);
    } else { // sigma
      h1->SetMinimum(0.);
      h1->SetMaximum(0.30);
    }
  } else if( !type.compare("dEta") ||
	     !type.compare("dPhi") ) { // ANGLES
    if( isMean ){
      h1->SetMinimum(-0.075);
      h1->SetMaximum(0.075);      
    }
  } else { // sigma
    h1->SetMinimum(0.015);
    h1->SetMaximum(0.15);
  }
}

double DiJetAnalysisMC::GetLineHeight( const std::string& type ){
  double y0 = 0;
  
  if( !type.compare("rPt") ){ // JES/JER
    y0 = 1;
    y0 = 1;
  } else if( !type.compare("rEta") ||
	     !type.compare("dPhi") ) { // ANGLES
    y0 = 0;
    y0 = 0;
  }

  return y0;
}
