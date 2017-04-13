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
#include <cfloat>
#include <boost/format.hpp>
#include <boost/filesystem.hpp>

#include "MyRoot.h"

#include "DiJetAnalysisMC.h"
#include "JetPair.h"

DiJetAnalysisMC::DiJetAnalysisMC() : DiJetAnalysisMC( true, true, 0 )
{}

DiJetAnalysisMC::DiJetAnalysisMC( bool isData, bool is_pPb, int mcType )
  : DiJetAnalysis( isData, is_pPb), m_mcType( mcType )
{
  //========== Set Histogram Binning =============
    // ------ truth binning --------
  m_ptTruthWidth  = 5;
  m_ptTruthMin    = 10;
  m_ptTruthMax    = 100;
  m_nPtTruthBins  =
    (m_ptTruthMax - m_ptTruthMin) / m_ptTruthWidth;

  // ---- JES/PRes/Etc ----- 
  m_nRPtRecoTruthBins   = 50;
  m_rPtRecoTruthMin     = 0;  m_rPtRecoTruthMax = 2;

  m_nDAngleRecoTruthBins = 50;
  m_dAngleRecoTruthMin   = -0.5; m_dAngleRecoTruthMax = 0.5;

  //==================== Cuts ====================    
  m_dRmax    = 0.2;
}

DiJetAnalysisMC::~DiJetAnalysisMC(){}

void DiJetAnalysisMC::Initialize(){
  m_labelOut = m_isData ? "data" : "mc" ;
  m_labelOut = m_is_pPb ? m_labelOut + "_pPb" : m_labelOut + "_pp";

  if( m_mcType == 0 ){
    m_labelOut += "_pythia8powheg";
    m_mcTypeLabel = "Pythia8+Powheg";
    
    m_vUsedJZN.push_back( "JZ1" );
    m_mJznFnameIn[ m_vUsedJZN.back() ]
      = "/home/yakov/Projects/atlas/data/jz1_pythia8powheg.root";
    m_mJznSigma  [ m_vUsedJZN.back() ] = 1.2794E+08;
    m_mJznEff    [ m_vUsedJZN.back() ] = 1.5857E-03;
    m_mJznSumPowhegWeights[ m_vUsedJZN.back() ] = 3.428918e+08;
    m_mJznPtThreshold[ m_vUsedJZN.back() ] = 20;
    
    m_vUsedJZN.push_back( "JZ2" );
    m_mJznFnameIn[ m_vUsedJZN.back() ] 
      = "/home/yakov/Projects/atlas/data/jz2_pythia8powheg.root";
    m_mJznSigma  [ m_vUsedJZN.back() ] = 1.9648E+07;
    m_mJznEff    [ m_vUsedJZN.back() ] = 1.2948E-04;
    m_mJznSumPowhegWeights[ m_vUsedJZN.back() ] = 4.656080e+06;
    m_mJznPtThreshold[ m_vUsedJZN.back() ] = 60;
    
    // Get Associated weights for powheg.
    // These need to be used at filling time (per event basis).
    TFile* weight_file =
      TFile::Open("/home/yakov/Projects/atlas/data/Powheg.reweight.root");
    m_hPowhegWeights = static_cast< TH3D* >( weight_file->Get("h3_pT_y_phi_rw") );
    m_hPowhegWeights->SetDirectory(0);
    weight_file->Close();
  } else if( m_mcType == 1 ){
    m_labelOut += "_pythia8";
    m_mcTypeLabel = "Pythia8";

    m_vUsedJZN.push_back( "JZ0" );
    m_mJznFnameIn[ m_vUsedJZN.back() ]
      = "/home/yakov/Projects/atlas/data/jz0_pythia8.root";
    m_mJznSigma  [ m_vUsedJZN.back() ] = 6.7890E+07;
    m_mJznEff    [ m_vUsedJZN.back() ] = 8.6287E-02;
    m_mJznSumPowhegWeights[ m_vUsedJZN.back() ] = 1;
    m_mJznPtThreshold[ m_vUsedJZN.back() ] = 10;
          
    m_vUsedJZN.push_back( "JZ1" );
    m_mJznFnameIn[ m_vUsedJZN.back() ]
      = "/home/yakov/Projects/atlas/data/jz1_pythia8.root";
    m_mJznSigma  [ m_vUsedJZN.back() ] = 6.7890E+07*1.2;
    m_mJznEff    [ m_vUsedJZN.back() ] = 2.8289E-03;
    m_mJznSumPowhegWeights[ m_vUsedJZN.back() ] = 1;
    m_mJznPtThreshold[ m_vUsedJZN.back() ] = 20;
    
    m_vUsedJZN.push_back( "JZ2" );
    m_mJznFnameIn[ m_vUsedJZN.back() ] 
      = "/home/yakov/Projects/atlas/data/jz2_pythia8.root";
    m_mJznSigma  [ m_vUsedJZN.back() ] = 6.3997E+05;
    m_mJznEff    [ m_vUsedJZN.back() ] = 4.2714E-03;
    m_mJznSumPowhegWeights[ m_vUsedJZN.back() ] = 1;
    m_mJznPtThreshold[ m_vUsedJZN.back() ] = 60;
  } else if( m_mcType == 2 ){
    m_labelOut += "_herwig";
    m_mcTypeLabel = "Herwig";

    m_vUsedJZN.push_back( "JZ0" );
    m_mJznFnameIn[ m_vUsedJZN.back() ]
      = "/home/yakov/Projects/atlas/data/jz0_herwig.root";
    m_mJznSigma  [ m_vUsedJZN.back() ] = 69.2475E+06;
    m_mJznEff    [ m_vUsedJZN.back() ] = 5.713013E-04;
    //m_mJznSigma  [ m_vUsedJZN.back() ] = 3.0028E+02;
    //m_mJznEff    [ m_vUsedJZN.back() ] = 2.0899E-02;
    m_mJznSumPowhegWeights[ m_vUsedJZN.back() ] = 1;
    m_mJznPtThreshold[ m_vUsedJZN.back() ] = 10;
    
    m_vUsedJZN.push_back( "JZ1" );
    m_mJznFnameIn[ m_vUsedJZN.back() ]
      = "/home/yakov/Projects/atlas/data/jz1_herwig.root";
    m_mJznSigma  [ m_vUsedJZN.back() ] = 69.2475E+06;
    m_mJznEff    [ m_vUsedJZN.back() ] = 1.5088E-03;
    //m_mJznSigma  [ m_vUsedJZN.back() ] = 3.1313E+02;
    //m_mJznEff    [ m_vUsedJZN.back() ] = 1.5062E-03;
    m_mJznSumPowhegWeights[ m_vUsedJZN.back() ] = 1;
    m_mJznPtThreshold[ m_vUsedJZN.back() ] = 20;
    
    m_vUsedJZN.push_back( "JZ2" );
    m_mJznFnameIn[ m_vUsedJZN.back() ] 
      = "/home/yakov/Projects/atlas/data/jz2_herwig.root";
    m_mJznSigma  [ m_vUsedJZN.back() ] = 4.2389E+05;
    m_mJznEff    [ m_vUsedJZN.back() ] = 3.21E-03;
    //m_mJznSigma  [ m_vUsedJZN.back() ] = 4.2389E+05;
    //m_mJznEff    [ m_vUsedJZN.back() ] = 3.2383E-03;
    m_mJznSumPowhegWeights[ m_vUsedJZN.back() ] = 1;
    m_mJznPtThreshold[ m_vUsedJZN.back() ] = 60;
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
  m_dirOut   += "/output_" + m_labelOut;
  checkWriteDir( m_dirOut.c_str() );

  m_rootFname = m_dirOut + "/myOut_" + m_labelOut + ".root";
  
  std::cout << "fNameIn/Out: " << m_rootFname << std::endl;

  //========== Cuts and triggers =============
  
  // calculate sum of sigma and eff
  m_sumSigmaEff = 0;
  for( auto jzn : m_vUsedJZN ){
    m_sumSigmaEff += ( m_mJznSigma[ jzn ] * m_mJznEff[ jzn ] );
  }

  std::cout << "SumSigmaEff = " << m_sumSigmaEff << std::endl;
  
  for( auto jzn : m_vUsedJZN ){ // loop over JZ samples
    m_fIn = TFile::Open( m_mJznFnameIn[ jzn ].c_str() );
    m_tree = static_cast< TTree* >( m_fIn->Get( "tree" ) );

    int nEventsTotal = m_tree->GetEntries();
    double     sigma = m_mJznSigma[ jzn ];      
    double       eff = m_mJznEff  [ jzn ];
  
    m_mJznNev    [ jzn ] = nEventsTotal;
    m_mJznWeights[ jzn ] =
      (1./nEventsTotal) * (1./m_sumSigmaEff) * (sigma * eff);
  
    std::cout << jzn << "   weight = "
	      << m_mJznWeights [jzn] << std::endl;
    m_fIn->Close();
  }
}

//---------------------------
// Fill Tree / Plot Controls
//---------------------------

void DiJetAnalysisMC::RunOverTreeFillHistos( int nEvents, 
					     int startEvent ){  
  SetupHistograms();
  ProcessEvents( nEvents, startEvent );
  SaveOutputsFromTree();
}

void DiJetAnalysisMC::ProcessPlotHistos(){
  LoadHistograms();

  std::string cfNameOut = m_dirOut + "/c_myOut_" + m_labelOut + ".root";
  m_fOut = new TFile( cfNameOut.c_str(),"RECREATE");
  
  // check if we have one  
  PlotEtaPhiPtMap( m_mJznEtaPhiMap );
  PlotEtaPhiPtMap( m_mJznEtaPtMap  );

  std::string sReco  = "Reco";
  std::string sTruth = "Truth";
  
  PlotSpectra( m_mJznEtaSpectReco , "spect", sReco );
  PlotSpectra( m_mJznEtaSpectTruth, "spect", sTruth );
  
  PlotEfficiencies( m_mJznEtaSpectTruthPaired, 
		    m_mJznEtaSpectTruth,
		    m_mJznEtaSpectTruthNent,
		    "eff" );
  
  PlotVsEtaPt( m_mJznRecoTruthRpt , m_mJznRecoTruthRptNent, "recoTruthRpt");
  /*
  PlotVsEtaPt( m_mJznRecoTruthDeta, m_mJznRecoTruthDetaNent, "recoTruthDeta");
  PlotVsEtaPt( m_mJznRecoTruthDphi, m_mJznRecoTruthDphiNent, "recoTruthDphi");
  */

  std::cout << "DONE! Closing " << cfNameOut << std::endl;
  m_fOut->Close();
  std::cout << "......Closed  " << cfNameOut << std::endl;
}

//---------------------------------
//            Read Data
//---------------------------------

void DiJetAnalysisMC::SetupHistograms(){
  
  for( auto jzn : m_vUsedJZN ){ // loop over JZ samples
    std::cout << "Making - " << jzn << " histograms " << std::endl;

    // -------- maps ---------
    m_mJznEtaPhiMap[ jzn ] =
      new TH2D( Form("h_etaPhiMap_%s", jzn.c_str() ),
		";#eta_{Reco};#phi_{Reco}",
		m_nEtaMapBins, m_etaMapMin, m_etaMapMax,
		m_nPhiMapBins, m_phiMapMin, m_phiMapMax );
    AddHistogram( m_mJznEtaPhiMap[ jzn ] );
      
    m_mJznEtaPtMap[ jzn ] =
      new TH2D( Form("h_etaPtMap_%s", jzn.c_str() ),
		";#eta_{Reco};#it{p}_{T}^{Reco}",
		m_nEtaMapBins, m_etaMapMin, m_etaMapMax,
		m_nPtMapBins , m_ptMapMin , m_ptMapMax );
    AddHistogram( m_mJznEtaPtMap[ jzn ] );

    // -------- spect --------
    m_mJznEtaSpectReco[ jzn ] = 
      new TH2D( Form("h_etaSpectReco_%s", jzn.c_str() ), 
		";#eta_{Reco};#it{p}_{T}^{Reco} [GeV]",
		m_nVarFwdEtaBins, 0, 1,
		m_nPtSpectBins,
		m_ptSpectMin, m_ptSpectMax ) ;
    m_mJznEtaSpectReco[ jzn ]->GetXaxis()->
      Set( m_nVarFwdEtaBins, &( m_varFwdEtaBinning[0] ) );
    AddHistogram( m_mJznEtaSpectReco[ jzn ] );

    m_mJznEtaSpectTruth[ jzn ] = 
      new TH2D( Form("h_etaSpectTruth_%s", jzn.c_str() ), 
		";#eta_{Truth};#it{p}_{T}^{Truth} [GeV]",
		m_nVarFwdEtaBins, 0, 1,
		m_nPtSpectBins,
		m_ptSpectMin, m_ptSpectMax ) ;
    m_mJznEtaSpectTruth[ jzn ]->GetXaxis()->
      Set( m_nVarFwdEtaBins, &( m_varFwdEtaBinning[0] ) );
    AddHistogram( m_mJznEtaSpectTruth[ jzn ] );

    // The Nent Histogram is only for the overlay samples
    // where we fill with some weight. Otherwise,
    // SpectTruth = SpectTruthNent (b/c fill with w=1)
    m_mJznEtaSpectTruthNent[ jzn ] = 
      new TH2D( Form("h_etaSpectTruthNent_%s", jzn.c_str() ), 
		";#eta_{Truth};#it{p}_{T}^{Truth} [GeV]",
		m_nVarFwdEtaBins, 0, 1,
		m_nPtSpectBins,
		m_ptSpectMin, m_ptSpectMax ) ;
    m_mJznEtaSpectTruthNent[ jzn ]->GetXaxis()->
      Set( m_nVarFwdEtaBins, &( m_varFwdEtaBinning[0] ) );
    AddHistogram( m_mJznEtaSpectTruthNent[ jzn ] );

    m_mJznEtaSpectTruthPaired[ jzn ] = 
      new TH2D( Form("h_etaSpectTruthPaired_%s", jzn.c_str() ), 
		";#eta_{Truth};#it{p}_{T}^{Truth} [GeV]",
		m_nVarFwdEtaBins, 0, 1,
		m_nPtSpectBins,
		m_ptSpectMin, m_ptSpectMax ) ;
    m_mJznEtaSpectTruthPaired[ jzn ]->GetXaxis()->
      Set( m_nVarFwdEtaBins, &( m_varFwdEtaBinning[0] ) );
    AddHistogram( m_mJznEtaSpectTruthPaired[ jzn ] );

    // --------- recoTruthRpt ---------
    m_mJznRecoTruthRpt[ jzn ] =
      new TH3D( Form("h_recoTruthRpt_%s", jzn.c_str()),
		";#eta^{Truth};#it{p}_{T}^{Truth};#it{p}_{T}^{Reco}/#it{p}_{T}^{Truth}",
		m_nEtaForwardBinsFine,
		m_etaForwardMin, m_etaForwardMax,
		m_nPtTruthBins,
		m_ptTruthMin, m_ptTruthMax,
		m_nRPtRecoTruthBins,
		m_rPtRecoTruthMin, m_rPtRecoTruthMax);
    AddHistogram( m_mJznRecoTruthRpt[ jzn ] );

    m_mJznRecoTruthRptNent[ jzn ] =
      new TH2D( Form("h_recoTruthRptNent_%s", jzn.c_str()),
		";#eta^{Truth};#it{p}_{T}^{Truth}",
		m_nEtaForwardBinsFine,
		m_etaForwardMin, m_etaForwardMax,
		m_nPtTruthBins,
		m_ptTruthMin, m_ptTruthMax );
    AddHistogram( m_mJznRecoTruthRptNent[ jzn ] );    

    // --------- recoTruthDeta ---------
    m_mJznRecoTruthDeta[ jzn ] =
      new TH3D( Form("h_recoTruthDeta_%s", jzn.c_str()),
		";#eta^{Truth};#it{p}_{T}^{Truth};#eta^{Reco}-#eta^{Truth}",
		m_nEtaForwardBinsFine,
		m_etaForwardMin, m_etaForwardMax,
		m_nPtTruthBins,
		m_ptTruthMin, m_ptTruthMax,
		m_nDAngleRecoTruthBins,
		m_dAngleRecoTruthMin, m_dAngleRecoTruthMax );
    AddHistogram( m_mJznRecoTruthDeta[ jzn ] );    
    
    m_mJznRecoTruthDetaNent[ jzn ] =
      new TH2D( Form("h_recoTruthDetaNent_%s", jzn.c_str()),
		";#eta^{Truth};#it{p}_{T}^{Truth}",
		m_nEtaForwardBinsFine,
		m_etaForwardMin, m_etaForwardMax,
		m_nPtTruthBins,
		m_ptTruthMin, m_ptTruthMax );
    AddHistogram( m_mJznRecoTruthDetaNent[ jzn ] );    

    // --------- recoTruthDphi ---------
    m_mJznRecoTruthDphi[ jzn ] =
      new TH3D( Form("h_recoTruthDphi_%s", jzn.c_str()),
		";#eta^{Truth};#it{p}_{T}^{Truth};#phi^{Reco}-#phi^{Truth}",
		m_nEtaForwardBinsFine,
		m_etaForwardMin, m_etaForwardMax,
		m_nPtTruthBins,
		m_ptTruthMin, m_ptTruthMax,
		m_nDAngleRecoTruthBins,
		m_dAngleRecoTruthMin, m_dAngleRecoTruthMax );
    AddHistogram( m_mJznRecoTruthDphi[ jzn ] );

    m_mJznRecoTruthDphiNent[ jzn ] =
      new TH2D( Form("h_recoTruthDphiNent_%s", jzn.c_str()),
		";#eta^{Truth};#it{p}_{T}^{Truth}",
		m_nEtaForwardBinsFine,
		m_etaForwardMin, m_etaForwardMax,
		m_nPtTruthBins,
		m_ptTruthMin, m_ptTruthMax );
    AddHistogram( m_mJznRecoTruthDphiNent[ jzn ] );    
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
  
    m_fIn  = TFile::Open( m_mJznFnameIn[ jzn ].c_str() );
    m_tree = static_cast< TTree*>( m_fIn->Get( "tree" ) );

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

    int nJetsTotal = 0;
    int nJetsForward = 0;

    // event loop
    for( m_ev = startEvent; m_ev < endEvent; m_ev++ ){
      m_tree->GetEntry( m_ev );

      ApplyCleaning ( vR_jets, v_isCleanJet );
      ApplyIsolation( 1.0, vR_jets );
      
      std::vector< JetPair > v_paired_jets;
      PairJets( vR_jets, vT_jets, v_paired_jets );
      
      if( anaTool->DoPrint(m_ev) ) {
	std::cout << "\nEvent : " << m_ev 
		  << "    has : " << vR_jets.size() << " reco jets"
		  << "    and : " << vT_jets.size() << " truth jets"
		  << std::endl; 
      }

      
      // fill for truth jets
      // denominator for efficiency because not all
      // reco jets are reconstructed for a truth jet
      // also count total, and fwd truth jets
      for( auto& tJet : vT_jets ){
	double jetEta    = tJet.Eta();
	double jetEtaAdj = AdjustEtaForPP( jetEta );
	double jetPhi    = tJet.Phi();
	double jetPt     = tJet.Pt()/1000.;

	double weight = GetJetWeight( jetEta, jetPhi, jetPt );

	m_mJznEtaSpectTruth    [ jzn ]->
	  Fill( jetEtaAdj, jetPt, weight);

	m_mJznEtaSpectTruthNent[ jzn ]->
	  Fill( jetEtaAdj, jetPt );
		
	// count how many total truth jets
	// and how many forward truth jets
	// cut on pt corresponding to jz sample
	if( jetPt > m_mJznPtThreshold[ jzn ] ){
	  nJetsTotal++;   
	  if( jetEtaAdj < -constants::FETAMIN ){
	    nJetsForward++;
	  }
	} 
      } // end loop over truth jets

      // loop over pairs
      for( auto& vp : v_paired_jets ){	  
	double jetEtaTruth = vp.TruthJet()->Eta();
	double jetPhiTruth = vp.TruthJet()->Phi();
	double  jetPtTruth = vp.TruthJet()->Pt()/1000.;
	  
	double jetEtaReco = vp.RecoJet()->Eta();
	double jetPhiReco = vp.RecoJet()->Phi();
	double  jetPtReco = vp.RecoJet()->Pt()/1000.;

	// convert positive eta to negative because
	// in pp it doesnt matter since detector is
	// symmetric in eta. i.e. eta < 0 <=> eta > 0
	// so we just take positive etas and put them
	// in the negative bins for pp configuration.
	// this saves some overhead with histograms later
	// our histos run in negative eta
	// (due to pPb configuration)
	// the labels will be taken care of so it is ok
	double jetEtaRecoAdj  = AdjustEtaForPP( jetEtaReco  );
	double jetEtaTruthAdj = AdjustEtaForPP( jetEtaTruth );	  
	
	double weight = GetJetWeight( jetEtaTruth, jetPhiTruth, jetPtTruth );
	
	m_mJznEtaPhiMap[ jzn ]->Fill( jetEtaReco, jetPhiReco, weight);
	m_mJznEtaPtMap [ jzn ]->Fill( jetEtaReco, jetPtReco , weight);
        		
	if( vp.DeltaR() <= m_dRmax ){
	  m_mJznEtaSpectReco  [ jzn ]->
	    Fill( jetEtaRecoAdj,  jetPtReco,  weight);

	  m_mJznEtaSpectTruthPaired[ jzn ]->
	    Fill( jetEtaTruthAdj,  jetPtTruth,  weight);

	  m_mJznRecoTruthRpt     [ jzn ]->
	    Fill( jetEtaTruthAdj, jetPtTruth, jetPtReco/jetPtTruth, weight);
	  m_mJznRecoTruthRptNent [ jzn ]->
	    Fill( jetEtaTruthAdj, jetPtTruth );

	  m_mJznRecoTruthDeta    [ jzn ]->
	    Fill( jetEtaTruthAdj, jetPtTruth, jetEtaReco - jetEtaTruth, weight);
	  m_mJznRecoTruthDetaNent[ jzn ]->
	    Fill( jetEtaTruthAdj, jetPtTruth );

	  m_mJznRecoTruthDphi    [ jzn ]->
	    Fill( jetEtaTruthAdj, jetPtTruth, jetPhiReco - jetPhiTruth, weight );
	  m_mJznRecoTruthDphiNent[ jzn ]->
	    Fill( jetEtaTruthAdj, jetPtTruth  );
	}
      } // end loop over pairs
    } // end loop over events
   
    std::cout << "DONE WITH JZ" << jzn
	      << "   nJetsForward: " << nJetsForward
	      << "   nTotalJets: "   << nJetsTotal
	      << "   nTotalEvents: " << nEventsTotal
	      << std::endl;

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
      if( anaTool->DeltaR( recoJet, truthJet ) <= deltaRmin ) {	
	pairedRecoJet = &recoJet;
	deltaRmin     = anaTool->DeltaR( recoJet, truthJet );
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
//           Plot Data
//---------------------------------
void DiJetAnalysisMC::LoadHistograms(){
  m_fIn = TFile::Open( m_rootFname.c_str() ); 

  for( auto& jzn : m_vUsedJZN ){
    // -------- maps ---------
    m_mJznEtaPhiMap[ jzn ] =
      static_cast< TH2D* >( m_fIn->Get( Form("h_etaPhiMap_%s", jzn.c_str() ) ) );
    m_mJznEtaPhiMap[ jzn ]->SetDirectory(0);
    m_mJznEtaPtMap [ jzn ] =
      static_cast< TH2D* >( m_fIn->Get( Form("h_etaPtMap_%s", jzn.c_str() ) ) );
    m_mJznEtaPtMap [ jzn ]->SetDirectory(0);

    // -------- spect --------
    m_mJznEtaSpectReco [ jzn ] =
      static_cast< TH2D* >( m_fIn->Get( Form("h_etaSpectReco_%s", jzn.c_str() ) ) );
    m_mJznEtaSpectReco [ jzn ]->SetDirectory(0);

    // The Nent Histogram is only for the overlay samples
    // where we fill with some weight. Otherwise,
    // SpectTruth = SpectTruthNent (b/c fill with w=1)
    m_mJznEtaSpectTruth    [ jzn ] =
      static_cast< TH2D* >( m_fIn->Get( Form("h_etaSpectTruth_%s", jzn.c_str() ) ) );
    m_mJznEtaSpectTruth    [ jzn ]->SetDirectory(0);
    m_mJznEtaSpectTruthNent[ jzn ] =
      static_cast< TH2D* >( m_fIn->Get( Form("h_etaSpectTruthNent_%s", jzn.c_str() ) ) );
    m_mJznEtaSpectTruthNent[ jzn ]->SetDirectory(0);
    m_mJznEtaSpectTruthPaired[ jzn ] =
      static_cast< TH2D* >( m_fIn->Get( Form("h_etaSpectTruthPaired_%s", jzn.c_str() ) ) );
    m_mJznEtaSpectTruthPaired[ jzn ]->SetDirectory(0);

    // --------- recoTruthRpt ---------
    m_mJznRecoTruthRpt     [ jzn ] =
      static_cast< TH3D* >( m_fIn->Get( Form("h_recoTruthRpt_%s", jzn.c_str() ) ) );
    m_mJznRecoTruthRpt     [ jzn ]->SetDirectory(0);
    m_mJznRecoTruthRptNent [ jzn ] =
      static_cast< TH2D* >( m_fIn->Get( Form("h_recoTruthRptNent_%s", jzn.c_str() ) ) );
    m_mJznRecoTruthRptNent [ jzn ]->SetDirectory(0);

    // --------- recoTruthDeta ---------
    m_mJznRecoTruthDeta    [ jzn ]   =
      static_cast< TH3D* >( m_fIn->Get( Form("h_recoTruthDeta_%s", jzn.c_str() ) ) );
    m_mJznRecoTruthDeta    [ jzn ]->SetDirectory(0);
    m_mJznRecoTruthDetaNent[ jzn ] =
      static_cast< TH2D* >( m_fIn->Get( Form("h_recoTruthDetaNent_%s", jzn.c_str() ) ) );
    m_mJznRecoTruthDetaNent[ jzn ]->SetDirectory(0);

    // --------- recoTruthDphi ---------
    m_mJznRecoTruthDphi    [ jzn ]   =
      static_cast< TH3D* >( m_fIn->Get( Form("h_recoTruthDphi_%s", jzn.c_str() ) ) );
    m_mJznRecoTruthDphi    [ jzn ]->SetDirectory(0);
    m_mJznRecoTruthDphiNent[ jzn ] =
      static_cast< TH2D* >( m_fIn->Get( Form("h_recoTruthDphiNent_%s", jzn.c_str() ) ) );
    m_mJznRecoTruthDphiNent[ jzn ]->SetDirectory(0);
  }
  m_fIn->Close();
}

void DiJetAnalysisMC::PlotSpectra( std::map< std::string, TH2* >& mJznSpect,
				   const std::string& type,
				   const std::string& level ){
  // check if we have one
  if( !mJznSpect.size() ){ return; }
  
  std::string yAxisTitle = "dN/d#it{p}_{T}";
  
  double ptSpectWidth =
    ( m_ptSpectMax - m_ptSpectMin ) / m_nPtSpectBins;

  // vector of projected spect
  std::vector< TH1* > vSpect;
  // vector of final spectra
  std::vector< TH1* > vSpectFinal;
  
  for( int xBin = 1;
       xBin <= mJznSpect.begin()->second->GetNbinsX();
       xBin++ ){

    // keep track of projections
    std::map< std::string, TH1* > mSpectTemp;
    
    TCanvas c_spect("c_spect","c_spect",800,600);
    c_spect.SetLogy();

    TLegend l_spect(0.78, 0.65, 0.99, 0.78);
    styleTool->SetLegendStyle( &l_spect );
    l_spect.SetFillStyle(0);

    int jznStyle =  1;
    double   max = -1;

    // should all be the same
    double etaMin = mJznSpect.begin()->second->
      GetXaxis()->GetBinLowEdge( xBin );
    double etaMax = mJznSpect.begin()->second->
      GetXaxis()->GetBinUpEdge( xBin );
  
    for( auto& jznHIN : mJznSpect ){ // loop over JZN
      std::string jzn = jznHIN.first;
    
      TH1* h_spect =
	jznHIN.second->
	ProjectionY( Form("h_%s%s_%s_%2.0f_Eta_%2.0f",
			  type.c_str(),
			  level.c_str(),
			  jzn.c_str(),
			  10*std::abs(etaMin),
			  10*std::abs(etaMax) ),
		     xBin, xBin );
      styleTool->SetHStyle( h_spect, jznStyle++ );

      h_spect->GetYaxis()->SetTitle( yAxisTitle.c_str() );
      h_spect->Scale( 1./ptSpectWidth );
      mSpectTemp[ jzn ] = h_spect ; 
      vSpect.push_back( h_spect );
      
      h_spect->Draw("epsame");
      max = h_spect->GetMaximum() > max ?
	h_spect->GetMaximum() : max;
  
      l_spect.AddEntry(  h_spect, Form("%s", jzn.c_str() ) );
      
      double power = log10(max);
      power = std::ceil(power); 
      max = pow( 10, power );
      // all spectra plots will have same max, min
      for( auto& jznSpect : mSpectTemp ){ 
	jznSpect.second->SetMaximum( max );
	jznSpect.second->SetMinimum( 0.1 );
      }
    } // end loop over JZN

    int  nPtBins = mJznSpect.begin()->second->GetNbinsY();
    double ptMin = mJznSpect.begin()->second->GetYaxis()->GetXmin();
    double ptMax = mJznSpect.begin()->second->GetYaxis()->GetXmax();

    TH1* h_spectFinal = new TH1D( Form("h_%s%sFinal_%2.0f_Eta_%2.0f_final",
				       type.c_str(),
				       level.c_str(),
				       10*std::abs(etaMin),
				       10*std::abs(etaMax) ),
				  Form("%s;%s;%s",
				       GetEtaLabel( etaMin, etaMax).c_str(),
				       mJznSpect.begin()->second->
				       GetYaxis()->GetTitle(),
				       yAxisTitle.c_str() ),
				  nPtBins, ptMin, ptMax );
    styleTool->SetHStyle( h_spectFinal, 0 );
    vSpectFinal.push_back( h_spectFinal );
  
    CombineJZN( h_spectFinal, mSpectTemp );
    
    h_spectFinal->Draw("epsame");

    l_spect.AddEntry( h_spectFinal, "Total" );
    l_spect.Draw();
    
    drawTool->DrawAtlasInternalMCRight( 0, 0, m_mcTypeLabel ); 
    drawTool->DrawLeftLatex
      ( 0.45, 0.87, GetEtaLabel( etaMin, etaMax ).c_str(), 1 );

    SaveAsAll( c_spect,
	       type, level, "Eta",
	       std::abs(etaMin)*10,
	       std::abs(etaMax)*10 );
  } // end loop over eta

  bool isLog = true;
  DrawCanvas( vSpectFinal, type, level, isLog );

  for( auto& h : vSpect      ){ delete h; }
  for( auto& h : vSpectFinal ){ delete h; }
}

void DiJetAnalysisMC::PlotVsEtaPt( std::map< std::string, TH3* >& mJznHIN,
				   std::map< std::string, TH2* >& mJznNentIN,
				   const std::string& type ){
  std::vector< TH1* > vMeans;
  std::vector< TH1* > vSigmas;
  std::vector< TH1* > vNent;
  
  std::vector< TH1* > vMeansFinal;
  std::vector< TH1* > vSigmasFinal;

  std::string sMean  = "mean";
  std::string sSigma = "sigma";
  
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

    if( !type.compare("recoTruthRpt") ){
      yTitleMean  = "#it{p}_{T}^{Reco}/#it{p}_{T}^{Truth}";
      yTitleSigma = "#sigma" + yTitleMean;
    } else if( !type.compare("recoTruthDeta") ){
      yTitleMean  = "#Delta#eta";
      yTitleSigma = "#sigma" + yTitleMean;
    } else if( !type.compare("recoTruthDphi") ){
      yTitleMean  = "#Delta#phi";
      yTitleSigma = "#sigma" + yTitleMean;
    } 

    std::map< std::string, TH1* > mJznMean;
    std::map< std::string, TH1* > mJznSigma;
    std::map< std::string, TH1* > mJznNent;

    int style = 1;
    
    for( auto& jznHIN : mJznHIN ){   // loop over JZN
      std::string jzn = jznHIN.first;
    
      TH1* h_mean = new TH1D( Form("h_%s_mean_%s_%2.0f_Eta_%2.0f",
				   type.c_str(),
				   jzn.c_str(),
				   10*std::abs(etaMin),
				   10*std::abs(etaMax) ),
			      Form("%s;%s;%s",
				   GetEtaLabel( etaMin, etaMax ).c_str(),
				   jznHIN.second->GetYaxis()->GetTitle(),
				   yTitleMean.c_str() ),
			      nPtBins, ptMin, ptMax );    
      styleTool->SetHStyle( h_mean, style );
      mJznMean[ jzn ] = h_mean;
      vMeans.push_back( h_mean );
    
      TH1* h_sigma = new TH1D( Form("h_%s_sigma_%s_%2.0f_Eta_%2.0f",
				    type.c_str(),
				    jzn.c_str(),
				    10*std::abs(etaMin),
				    10*std::abs(etaMax) ),
			       Form("%s;%s;%s",
				    GetEtaLabel( etaMin, etaMax ).c_str(),
				    jznHIN.second->GetYaxis()->GetTitle(),
				    yTitleSigma.c_str() ),
			       nPtBins, ptMin, ptMax );    
      styleTool->SetHStyle( h_sigma, style );
      mJznSigma[ jzn ] = h_sigma;
      vSigmas.push_back( h_sigma );

      TH1* hNent = mJznNentIN[ jzn ]->
	ProjectionY( Form("h_%s_N_%s_%2.0f_Eta_%2.0f",
			  type.c_str(),
			  jzn.c_str(),
			  10*std::abs(etaMin),
			  10*std::abs(etaMax) ),
		     xBin, xBin);
      mJznNent[ jzn ] = hNent;
      vNent.push_back( hNent );
      
      ProjectAndFitGaus( jznHIN.second, h_mean, h_sigma,
			 xBin, xBin, jzn, m_mcTypeLabel );
      
      // increment style
      style++;

    } // end loop over JZN

    // save eta ranges in title to draw all etas on one canvas later
    // this is for the legends.... 
    TH1* h_meanFinal = new TH1D( Form("h_%s_mean_%2.0f_Eta_%2.0f_final",
				      type.c_str(),
				      10*std::abs(etaMin),
				      10*std::abs(etaMax) ),
				 Form("%s;%s;%s",
				      GetEtaLabel( etaMin, etaMax).c_str(),
				      mJznHIN.begin()->second->
				      GetYaxis()->GetTitle(),
				      yTitleMean.c_str() ),
				 nPtBins, ptMin, ptMax );
    styleTool->SetHStyle( h_meanFinal, 0 );
    vMeansFinal.push_back( h_meanFinal );
  
    TH1* h_sigmaFinal = new TH1D( Form("h_%s_sigma_%2.0f_Eta_%2.0f_final",
				       type.c_str(),
				       10*std::abs(etaMin),
				       10*std::abs(etaMax) ),
				  Form("%s;%s;%s",
				       GetEtaLabel( etaMin, etaMax).c_str(),
				       mJznHIN.begin()->second->
				       GetYaxis()->GetTitle(),
				       yTitleSigma.c_str() ),
				  nPtBins, ptMin, ptMax );
    styleTool->SetHStyle( h_sigmaFinal, 0 );
    vSigmasFinal.push_back( h_sigmaFinal );

    CombineJZN( h_meanFinal , mJznMean , mJznNent );
    CombineJZN( h_sigmaFinal, mJznSigma, mJznNent );
    
    // Draw the canvases for mean and sigma
    DrawCanvas( mJznMean,  h_meanFinal, type, sMean, etaMin, etaMax ); 
    DrawCanvas( mJznSigma, h_sigmaFinal, type, sSigma, etaMin, etaMax ); 
  } // end loop over eta

  // final means in a few eta bins.
  // since these histograms have 12 bins
  // we draw every 3 bins
  int spacing = 3;
  DrawCanvas( vMeansFinal , type, sMean , spacing );
  DrawCanvas( vSigmasFinal, type, sSigma, spacing );

  for( auto& h : vMeans  ){ delete h; }
  for( auto& h : vSigmas ){ delete h; }
  for( auto& h : vNent   ){ delete h; }
  
  for( auto& h : vMeansFinal  ){ delete h; }
  for( auto& h : vSigmasFinal ){ delete h; }
}

void DiJetAnalysisMC::
PlotEfficiencies( std::map< std::string, TH2* >& mJznSpectPaired,
		  std::map< std::string, TH2* >& mJznSpect,
		  std::map< std::string, TH2* >& mJznSpectNent,
		  const std::string& type ){  
  // if there are none, return
  if( !mJznSpect.size() ){ return; }  

  double lX0 = 0.13;
  double lY0 = 0.75;
  double lX1 = 0.39;
  double lY1 = 0.87;

  double xMin = 10.; 
  double xMax = 90.; 

  std::string gTitle = ";#it{p}_{T}^{Truth} [GeV];#it{#varepsilon}_{Reco}";

  std::vector< TH1* > vSpect;
  std::vector< TGraphAsymmErrors* > vEffGrf;
  std::vector< TGraphAsymmErrors* > vEffGrfFinal;
  
  for( int xBin = 1; xBin <= mJznSpect.begin()->second->GetNbinsX(); xBin++ ){
    TCanvas c_eff("c_eff","c_eff",800,600);
    styleTool->SetCStyleEff( c_eff, xMin, m_effMin, xMax, m_effMax, gTitle );
 
    TLegend l_eff( lX0, lY0, lX1, lY1);
    styleTool->SetLegendStyle( &l_eff );
    l_eff.SetFillStyle(0);

    // should all be the same
    double etaMin = mJznSpect.begin()->second->
      GetXaxis()->GetBinLowEdge( xBin );
    double etaMax = mJznSpect.begin()->second->
      GetXaxis()->GetBinUpEdge( xBin );

    // local inside the loop. used only on a per eta-bin
    // basis. the actual projections are saved to the vectors
    // that are global in this function.
    std::map< std::string, TH1* > mSpect;
    std::map< std::string, TH1* > mSpectPaired;
    std::map< std::string, TH1* > mNent;
    std::map< std::string, TGraphAsymmErrors* > mEff;
    
    for( auto& jznSpect : mJznSpect ){
      std::string jzn = jznSpect.first;

      mSpectPaired[ jzn ] =
	mJznSpectPaired[ jzn ]->
	ProjectionY( Form("h_%s_%s_paired_%2.0f_Eta_%2.0f",
			  type.c_str(),
			  jzn.c_str(),
			  10*std::abs(etaMin),
			  10*std::abs(etaMax) ),
		     xBin, xBin);
      vSpect.push_back( mSpectPaired[ jzn ] );

      mSpect[ jzn ] =
	jznSpect.second->
	ProjectionY( Form("h_%s_%s_%2.0f_Eta_%2.0f",
			  type.c_str(),
			  jzn.c_str(),
			  10*std::abs(etaMin),
			  10*std::abs(etaMax) ),
		     xBin, xBin);
      vSpect.push_back( mSpect[ jzn ] );
      
      mNent[ jzn ] =
	mJznSpectNent[ jzn ]->
	ProjectionY( Form("h_%s_%s_nEnt_%2.0f_Eta_%2.0f",
			  type.c_str(),
			  jzn.c_str(),
			  10*std::abs(etaMin),
			  10*std::abs(etaMax) ),
		     xBin, xBin);
      vSpect.push_back( mNent[ jzn ] ); 
    }

    int jznStyle = 1;
    for( auto& jznSpect : mSpect ){
      std::string   jzn = jznSpect.first;
      TH1*       hSpect = jznSpect.second;
      TH1* hSpectPaired = mSpectPaired[ jzn ];

      TGraphAsymmErrors* g_etaEff = new TGraphAsymmErrors();
      styleTool->SetHStyle( g_etaEff, jznStyle++ );
      
      g_etaEff->SetName ( Form("gr_%s_%s", type.c_str(), jzn.c_str() ) );
      vEffGrf.push_back( g_etaEff );
      mEff[ jzn ] = g_etaEff;
      
      l_eff.AddEntry( g_etaEff, Form("%s", jzn.c_str() ) );

      g_etaEff->Divide( hSpectPaired, hSpect,
			"cl=0.683 b(1,1) mode" );
      
      
      g_etaEff->Draw("p");
    } // end loop over JZN

    TGraphAsymmErrors* g_effFinal =  CombineJZN( mEff, mNent );
    if( !g_effFinal ){ continue; } // something wrong. fix later.
    g_effFinal->SetName( Form("gr_%s_%2.0f_Eta_%2.0f",
			      type.c_str(),
			      10*std::abs(etaMin),
			      10*std::abs(etaMax) ) );
    g_effFinal->SetTitle( GetEtaLabel( etaMin, etaMax ).c_str() );
    styleTool->SetHStyle( g_effFinal, 0 );
    vEffGrfFinal.push_back( g_effFinal );

    g_effFinal->Draw("p");

    l_eff.AddEntry( g_effFinal, "Total" );
    l_eff.Draw();

    TLine line( xMin, 1, xMax, 1);
    line.Draw();

    drawTool->DrawAtlasInternalDataRight( 0, 0, m_is_pPb ); 
    drawTool->DrawRightLatex
      ( 0.88, 0.17, GetEtaLabel( etaMin, etaMax).c_str(), 1 );

    
    SaveAsAll( c_eff, type, "", "Eta", std::abs(etaMin)*10, std::abs(etaMax)*10 );
  } // end loop over eta

  DrawCanvas( vEffGrfFinal, type, gTitle, xMin, xMax );

  for( auto& h : vSpect       ){ delete h; }
  for( auto& g : vEffGrf      ){ delete g; }
  for( auto& g : vEffGrfFinal ){ delete g; }
}

void DiJetAnalysisMC::PlotEtaPhiPtMap( std::map< std::string, TH2* >& mJznHIN ){
  TCanvas c_map("c_map","c_map",800,600);

  for( auto& jznHIN : mJznHIN ){
    jznHIN.second->Draw("col");
    styleTool->SetHStyle( jznHIN.second, 0 );
    drawTool->DrawAtlasInternalMCLeft( 0, -0.55, m_mcTypeLabel  );  
    SaveAsPdfPng( c_map, jznHIN.second->GetName() );
  }
}

//---------------------------
//          Tools 
//---------------------------

void DiJetAnalysisMC::CombineJZN( TH1* h_res,
				  std::map< std::string, TH1* >& mJznHIN){
  for( auto& jznH : mJznHIN ){
    std::string jzn      = jznH.first;
    double scale = m_mJznEff[ jzn ] * m_mJznSigma[ jzn ] /
      m_mJznSumPowhegWeights[ jzn ];
    h_res->Add( jznH.second, scale );
  }
  h_res->Scale( 1./m_sumSigmaEff );
}

void DiJetAnalysisMC::CombineJZN( TH1* h_res,
				  std::map< std::string, TH1*>& mJznVIN,
				  std::map< std::string, TH1*>& mJznNentIN ){
  for( int xBin = 1; xBin <= h_res->GetNbinsX(); xBin++ ){

    double valTot      = 0;
    double valErrorTot = 0;
    double denomTot    = 0;

    for( auto& jznV : mJznVIN ){
      std::string jzn = jznV.first;

      double valueBin    = mJznVIN      [ jzn ]->GetBinContent( xBin );
      double valueBinErr = mJznVIN      [ jzn ]->GetBinError( xBin );
      double nEntriesBin = mJznNentIN   [ jzn ]->GetBinContent( xBin );
      double weight      = m_mJznWeights[ jzn ] / m_mJznSumPowhegWeights[ jzn ];

      if( nEntriesBin < m_nMinEntriesGausFit ){ continue; }
      
      double val    = weight * nEntriesBin * valueBin;
      double valErr = weight * nEntriesBin * valueBinErr;
      valTot       += val;
      valErrorTot  += valErr * valErr;

      double denom  = weight * nEntriesBin;
      denomTot     += denom;
    }

    double valFinal      = valTot / denomTot;
    double valErrorFinal = valErrorTot / ( denomTot * denomTot );
    valErrorFinal        = std::sqrt( valErrorFinal );

    // check if we have NaN from
    // dividing valTot by zero (denomTot)
    if( std::isnan( valFinal ) ){ continue; }
    
    h_res->SetBinContent( xBin, valFinal );
    h_res->SetBinError  ( xBin, valErrorFinal );
  }
}

TGraphAsymmErrors* DiJetAnalysisMC::
CombineJZN( std::map< std::string, TGraphAsymmErrors*>& mJznVIN,
	    std::map< std::string, TH1*>& mJznNentIN ){
  
  // check if we have one
  if( !mJznNentIN.size() ){ return NULL; }

  // temp histos for this scope
  // fill values from TGraph into these histos
  std::map< std::string, TH1* > mJznVal;
  std::map< std::string, TH1* > mJznValErrorLow;
  std::map< std::string, TH1* > mJznValErrorHigh;

  int  nXbins = mJznNentIN.begin()->second->GetNbinsX();
  double xMin = mJznNentIN.begin()->second->GetXaxis()->GetXmin();
  double xMax = mJznNentIN.begin()->second->GetXaxis()->GetXmax();
  
  for( auto& jznNent : mJznNentIN ){
    std::string jzn = jznNent.first;
        
    mJznVal         [ jzn ] =
      new TH1D( Form("h_jznVal_%s", jzn.c_str() ),
		Form("h_jznVal_%s", jzn.c_str() ),
		nXbins, xMin, xMax );

    mJznValErrorLow [ jzn ] =
      new TH1D( Form("h_jznValErrorLow_%s", jzn.c_str() ),
		Form("h_jznValErrorLow_%s", jzn.c_str() ),
		nXbins, xMin, xMax );

    mJznValErrorHigh[ jzn ] =
      new TH1D( Form("h_jznValErrorHigh_%s", jzn.c_str() ),
		Form("h_jznValErrorHigh_%s", jzn.c_str() ),
		nXbins, xMin, xMax );
    
    TGraphAsymmErrors* gPts = mJznVIN[ jzn ];
    double x, y;
    
    for( int i = 0; i < gPts->GetN(); i++ ){
      gPts->GetPoint(i, x, y);
      int bin = mJznVal[ jzn ]->FindBin( x );

      double eYlow  = gPts->GetErrorYlow ( i );
      double eYhigh = gPts->GetErrorYhigh( i );
      
      mJznVal         [ jzn ]->SetBinContent( bin, y      );
      mJznValErrorLow [ jzn ]->SetBinContent( bin, eYlow  );
      mJznValErrorHigh[ jzn ]->SetBinContent( bin, eYhigh );
    }    
  }

  std::vector< double > vX;
  std::vector< double > vY;
  std::vector< double > vEx;
  std::vector< double > vEyLow;
  std::vector< double > vEyHigh;
  
  for( int xBin = 1; xBin <= nXbins; xBin++ ){

    double valTot   = 0;
    double denomTot = 0;
    double valErrorLowTot  = 0;
    double valErrorHighTot = 0;

    double xCent  = mJznNentIN.begin()->second->GetBinCenter( xBin );
    double xWidth = mJznNentIN.begin()->second->GetBinWidth ( xBin );
    
    for( auto& jznV : mJznVal ){
      std::string jzn = jznV.first;

      double nEntriesBin = mJznNentIN   [ jzn ]->GetBinContent( xBin );
      double weight      = m_mJznWeights[ jzn ] / m_mJznSumPowhegWeights[ jzn ];

      double valueBin    = jznV.second->GetBinContent( xBin );
     
      double valueBinErrLow  = mJznValErrorLow [ jzn ]->GetBinContent( xBin );
      double valueBinErrHigh = mJznValErrorHigh[ jzn ]->GetBinContent( xBin );

      if( valueBin == 0 && valueBinErrLow == 0 &&  valueBinErrHigh == 0 )
	{ continue; }
      
      if( nEntriesBin < m_nMinEntriesGausFit )
	{ continue; }
      
      double val         = weight * nEntriesBin * valueBin;
      double valErrLow   = weight * nEntriesBin * valueBinErrLow;
      double valErrHigh  = weight * nEntriesBin * valueBinErrHigh;
      valTot            += val;
      valErrorLowTot    += valErrLow  * valErrLow;
      valErrorHighTot   += valErrHigh * valErrHigh;

      double denom  = weight * nEntriesBin;
      denomTot     += denom;
    }

    double valFinal      = valTot / denomTot;
    
    double valErrorLowFinal  = valErrorLowTot  / ( denomTot * denomTot );
    double valErrorHighFinal = valErrorHighTot / ( denomTot * denomTot );

    valErrorLowFinal  = std::sqrt( valErrorLowFinal );
    valErrorHighFinal = std::sqrt( valErrorHighFinal );
    
    // check if we have NaN from
    // dividing valTot by zero (denomTot)
    if( std::isnan( valFinal ) ){ continue; }

    vX.push_back ( xCent  ); vY.push_back( valFinal );
    vEx.push_back( xWidth * 0.5 );
    vEyLow.push_back ( valErrorLowFinal  );
    vEyHigh.push_back( valErrorHighFinal );
  }

  for( auto& jznH : mJznVal          ){ delete jznH.second; }
  for( auto& jznH : mJznValErrorLow  ){ delete jznH.second; }
  for( auto& jznH : mJznValErrorHigh ){ delete jznH.second; }
    
  return new TGraphAsymmErrors( vX.size(),
				&vX[0] , &vY[0],
				&vEx[0], &vEx[0],
				&vEyLow[0], &vEyHigh[0] );
}

double DiJetAnalysisMC::GetJetWeight( double eta, double phi, double pt ){
  if( m_mcType != 0 ) { return 1; } // only for powheg
  
  int xb = m_hPowhegWeights->GetXaxis()->FindBin(pt);
  int yb = m_hPowhegWeights->GetYaxis()->FindBin(eta);
  int zb = m_hPowhegWeights->GetZaxis()->FindBin(phi);
  float jet_weight = m_hPowhegWeights->GetBinContent(xb,yb,zb);

  return jet_weight;
}


//---------------------------
//          Drawing 
//---------------------------
void DiJetAnalysisMC::DrawCanvas( std::map< std::string, TH1* >& mJznHIN,
				  TH1* hFinal,
				  const std::string& type1,
				  const std::string& type2,
				  double etaMin, double etaMax ){
  TCanvas c("c","c",800,600);
  
  TLegend leg(0.13, 0.14, 0.44, 0.27);
  styleTool->SetLegendStyle( &leg );
  leg.SetFillStyle(0);

  for( auto& jznHIN : mJznHIN ){
    jznHIN.second->SetTitle("");
    jznHIN.second->Draw("epsame");
    SetMinMax( jznHIN.second, type1, type2 );
    leg.AddEntry( jznHIN.second, jznHIN.first.c_str() );
  }
  hFinal->Draw("same");
  leg.AddEntry( hFinal, "Total" );
  
  leg.Draw();
  
  drawTool->DrawAtlasInternalMCRight( 0, 0, m_mcTypeLabel ); 
  drawTool->DrawLeftLatex
    ( 0.45, 0.874,GetEtaLabel( etaMin, etaMax ).c_str(), 1 );
  
  double y0 = GetLineHeight( type1 );
  
  double xMin = mJznHIN.begin()->second->GetXaxis()->GetXmin();
  double xMax = mJznHIN.begin()->second->GetXaxis()->GetXmax();

  TLine line( xMin, y0, xMax, y0);
  // dont draw line for sigma plots
  if( type2.compare("sigma") ) { line.Draw(); }
 
  SaveAsAll( c, type1, type2, "Eta",
	     std::abs(etaMin)*10, std::abs(etaMax)*10 );
}

void DiJetAnalysisMC::DrawCanvas( std::vector< TH1* >& vHIN,
				  const std::string& type1,
				  const std::string& type2,
				  int spacing ){
  TCanvas c("c","c",800,600);
  
  TLegend leg(0.68, 0.64, 0.99, 0.77);
  styleTool->SetLegendStyle( &leg );
  leg.SetFillStyle(0);

  int style = 0;

  // for situations where dont want to
  // plot every single bin 
  // plot every n on canvas
  int recoTruthDeta = vHIN.size()/spacing; // plot every n
  for( unsigned int etaRange = 0;
       etaRange < vHIN.size();
       etaRange += recoTruthDeta){
    styleTool->SetHStyle( vHIN[ etaRange], style++ );
    leg.AddEntry( vHIN[ etaRange ], vHIN[ etaRange ]->GetTitle() );
    vHIN[ etaRange ]->SetTitle("");
    vHIN[ etaRange ]->Draw("epsame");
    SetMinMax( vHIN[ etaRange ], type1, type2 );
  }

  leg.Draw();
  
  double y0 = GetLineHeight( type1 );
  
  double xMin = vHIN.front()->GetXaxis()->GetXmin();
  double xMax = vHIN.front()->GetXaxis()->GetXmax();
  
  TLine line( xMin, y0, xMax, y0);
  line.Draw();

  drawTool->DrawAtlasInternalMCRight( 0, 0, m_mcTypeLabel ); 
  SaveAsAll( c, type1, type2 );
}

void DiJetAnalysisMC::DrawCanvas( std::vector< TH1* >& vHIN,
				  const std::string& type1,
				  const std::string& type2,
				  bool logY ){
  TCanvas c("c","c",800,600);
  if( logY ) { c.SetLogy(); }
  
  TLegend leg(0.68, 0.64, 0.99, 0.77);
  styleTool->SetLegendStyle( &leg );
  leg.SetFillStyle(0);

  double max = -1;
  for( auto& h : vHIN ){
    max = h->GetMaximum() > max ? h->GetMaximum() : max;
  }

  if( logY ){
    double power = log10(max);
    power = std::ceil(power); 
    max   = pow( 10, power );
  }
  
  int style = 0;  
  for( auto& h : vHIN ){
    styleTool->SetHStyle( h, style++ );
    leg.AddEntry( h, h->GetTitle() );
    h->SetTitle("");
    h->Draw("epsame");
    h->SetMaximum( max );
    h->SetMinimum( 0.1 );
  }

  drawTool->DrawAtlasInternalMCRight( 0, 0, m_mcTypeLabel  ); 
  leg.Draw();

  SaveAsAll( c, type1, type2 ); 
}


void DiJetAnalysisMC::DrawCanvas( std::vector< TGraphAsymmErrors* >& vGIN,
				  const std::string& type,
				  const std::string& title,
				  double xMin, double xMax ){
  TCanvas c("c","c",800,600);
  styleTool->SetCStyleEff( c, xMin, m_effMin, xMax, m_effMax,
			    title.c_str() );
   
  TLegend leg(0.64, 0.20, 0.95, 0.34);
  styleTool->SetLegendStyle( &leg );
  leg.SetFillStyle(0);
  
  int style = 0;  
  for( auto& gr : vGIN ){
    styleTool->SetHStyle( gr, style++ );
    leg.AddEntry( gr, gr->GetTitle() );
    gr->SetTitle("");
    gr->Draw("epsame");
  }

  drawTool->DrawAtlasInternalMCRight( 0, 0, m_mcTypeLabel  ); 
  leg.Draw();

  TLine line( xMin, 1, xMax, 1);
  line.Draw();
  
  SaveAsAll( c, type ); 
}

//===== MinMax and line drawing =====

void DiJetAnalysisMC::
SetMinMax( TH1* h1, const std::string& type1, const std::string& type2 ){
  // JES JER
  if( !type1.compare("recoTruthRpt") ){ 
    if( !type2.compare("mean") ){ // sigma
      h1->SetMaximum(1.25);
      h1->SetMinimum(0.75);
    } else if( !type2.compare("sigma") ){ // sigma
      h1->SetMaximum(0.34);
      h1->SetMinimum(0.);
    }
  }
  // ANGLES
  else if( !type1.compare("recoTruthDeta") ||
	   !type1.compare("recoTruthDphi") ) { 
    if( !type2.compare("mean") ){ // mean
      h1->SetMaximum(0.075);      
      h1->SetMinimum(-0.075);
    } else if( !type2.compare("sigma") ){ // sigma
      h1->SetMaximum(0.056);
      h1->SetMinimum(0.);
    } 
  } 
}

double DiJetAnalysisMC::GetLineHeight( const std::string& type ){
  double y0 = 0;
  
  if( !type.compare("recoTruthRpt") ){ // JES/JER
    y0 = 1;
    y0 = 1;
  } else if( !type.compare("rEta") ||
	     !type.compare("recoTruthDphi") ) { // ANGLES
    y0 = 0;
    y0 = 0;
  }

  return y0;
}
