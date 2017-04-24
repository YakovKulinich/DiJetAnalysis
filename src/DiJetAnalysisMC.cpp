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
#include <boost/assign.hpp>

#include "MyRoot.h"

#include "DiJetAnalysisMC.h"
#include "JetPair.h"

DiJetAnalysisMC::DiJetAnalysisMC() : DiJetAnalysisMC( true, true , 0 )
{}

DiJetAnalysisMC::DiJetAnalysisMC( bool isData, bool is_pPb, int mcType )
  : DiJetAnalysis( isData, is_pPb, mcType )
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

    boost::assign::push_back( m_vJznUsed )
      ( 1 )( 2 );
    boost::assign::push_back( m_vJznLabel )
      ( "jz1" )( "jz2" );
    boost::assign::push_back( m_vJznFnameIn )
      ( "/home/yakov/Projects/atlas/data/jz1_pythia8powheg.root" )
      ( "/home/yakov/Projects/atlas/data/jz2_pythia8powheg.root" );
    boost::assign::push_back( m_vJznSigma )
      ( 1.2794E+08 )( 1.9648E+07 );
    boost::assign::push_back( m_vJznEff )
      ( 1.5857E-03 )( 1.2948E-04 );
    boost::assign::push_back( m_vJznSumPowhegWeights )
      ( 3.428918e+08 )( 4.656080e+06 );
    boost::assign::push_back( m_vJznPtThreshold )
      ( 20 )( 60 );
    				 
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

    boost::assign::push_back( m_vJznUsed )
      ( 0 )( 1 )( 2 );
    boost::assign::push_back( m_vJznLabel )
      ( "jz0" )( "jz1" )( "jz2" );
    boost::assign::push_back( m_vJznFnameIn )
      ( "/home/yakov/Projects/atlas/data/jz0_pythia8.root" )
      ( "/home/yakov/Projects/atlas/data/jz1_pythia8.root" )
      ( "/home/yakov/Projects/atlas/data/jz2_pythia8.root" );
    boost::assign::push_back( m_vJznSigma )
      ( 6.7890E+07 )( 6.7890E+07*1.2 )( 6.3997E+05 );
    boost::assign::push_back( m_vJznEff )
      ( 8.6287E-02 )( 2.8289E-03 )( 4.2714E-03 );
    boost::assign::push_back( m_vJznSumPowhegWeights )
      ( 1 )( 1 )( 1 );   
    boost::assign::push_back( m_vJznPtThreshold )
      ( 10 )( 20 )( 60 );

  } else if( m_mcType == 2 ){
    m_labelOut += "_herwig";
    m_mcTypeLabel = "Herwig";

    boost::assign::push_back( m_vJznUsed )
      ( 0 )( 1 )( 2 );
    boost::assign::push_back( m_vJznLabel )
      ( "jz0")( "jz1" )( "jz2" );
    boost::assign::push_back( m_vJznFnameIn )
      ( "/home/yakov/Projects/atlas/data/jz0_herwig.root" )
      ( "/home/yakov/Projects/atlas/data/jz1_herwig.root" )
      ( "/home/yakov/Projects/atlas/data/jz2_herwig.root" );
    boost::assign::push_back( m_vJznSigma )
      ( 69.2475E+06 )( 69.2475E+06 )( 4.2389E+05 );
    boost::assign::push_back( m_vJznEff )
      ( 5.713013E-04 )( 1.5088E-03 )( 3.21E-03 );
    boost::assign::push_back( m_vJznSumPowhegWeights )
      ( 1 )( 1 )( 1 );   
    boost::assign::push_back( m_vJznPtThreshold )
      ( 10 )( 20 )( 60 );
  }  

  m_nJzn = m_vJznUsed.size();  
  
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

  for( uint iG = 0; iG < m_nJzn; iG++ )
    { m_sumSigmaEff += ( m_vJznSigma[iG] * m_vJznEff[iG] ); }

  std::cout << "SumSigmaEff = " << m_sumSigmaEff << std::endl;

  for( uint iG = 0; iG < m_nJzn; iG++ ){
    m_fIn = TFile::Open( m_vJznFnameIn[iG].c_str() );
    m_tree = static_cast< TTree* >( m_fIn->Get( "tree" ) );

    int nEventsTotal = m_tree->GetEntries();
    double     sigma = m_vJznSigma[iG];      
    double       eff = m_vJznEff  [iG];
  
    m_vJznNev.push_back( nEventsTotal );
    m_vJznWeights.push_back
      ( (1./nEventsTotal) * (1./m_sumSigmaEff) * (sigma * eff) );
  
    std::cout << iG << "   weight = "
	      << m_vJznWeights.back() << std::endl;
    
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
  
  PlotEtaPhiPtMap( m_vHjznEtaPhiMap );
  PlotEtaPhiPtMap( m_vHjznEtaPtMap  );

  std::string sReco  = "reco";
  std::string sTruth = "truth";
  
  PlotSpectra( m_vHjznEtaSpectReco , "spect", sReco );
  PlotSpectra( m_vHjznEtaSpectTruth, "spect", sTruth );

  PlotDeltaPhi( m_vHjznDphiReco , m_vJznLabel, sReco, m_mcTypeLabel  );
  PlotDeltaPhi( m_vHjznDphiTruth, m_vJznLabel, sTruth, m_mcTypeLabel );
  
  PlotVsEtaPt( m_vHjznRecoTruthRpt , m_vHjznRecoTruthRptNent, "recoTruthRpt");
  /*
  PlotEfficiencies( m_vHjznEtaSpectTruthPaired, 
		    m_vHjznEtaSpectTruth,
		    m_vHjznEtaSpectTruthNent,
		    "eff" );
  
  // PlotVsEtaPt( m_vHjznRecoTruthDeta, m_vHjznRecoTruthDetaNent, "recoTruthDeta");
  // PlotVsEtaPt( m_vHjznRecoTruthDphi, m_vHjznRecoTruthDphiNent, "recoTruthDphi");
 */

  std::cout << "DONE! Closing " << cfNameOut << std::endl;
  m_fOut->Close();
  std::cout << "......Closed  " << cfNameOut << std::endl;
}


//---------------------------------
//            Read Data
//---------------------------------

void DiJetAnalysisMC::SetupHistograms(){
  
  for( uint iG = 0; iG < m_nJzn; iG++){
    int jzn = m_vJznUsed[iG];

    std::cout << "Making - jz" << jzn << " histograms " << std::endl;
    
    // -------- maps ---------
    m_vHjznEtaPhiMap.
      push_back( new TH2D( Form("h_etaPhiMap_jz%d", jzn ),
			  ";#eta_{Reco};#phi_{Reco}",
			  m_nEtaMapBins, m_etaMapMin, m_etaMapMax,
			  m_nPhiMapBins, m_phiMapMin, m_phiMapMax ) );
    AddHistogram( m_vHjznEtaPhiMap.back() );
      
    m_vHjznEtaPtMap.
      push_back( new TH2D( Form("h_etaPtMap_jz%d", jzn ),
			   ";#eta_{Reco};#it{p}_{T}^{Reco}",
			   m_nEtaMapBins, m_etaMapMin, m_etaMapMax,
			   m_nPtMapBins , m_ptMapMin , m_ptMapMax ) );
    AddHistogram( m_vHjznEtaPtMap.back() );
    
    // -------- spect --------
    m_vHjznEtaSpectReco.push_back
      ( new TH2D( Form("h_etaSpectReco_jz%d", jzn ), 
		  ";#eta_{Reco};#it{p}_{T}^{Reco} [GeV]",
		  m_nVarFwdEtaBins, 0, 1,
		  m_nPtSpectBins,
		  m_ptSpectMin, m_ptSpectMax ) );
    m_vHjznEtaSpectReco.back()->GetXaxis()->
      Set( m_nVarFwdEtaBins, &( m_varFwdEtaBinning[0] ) );
    AddHistogram( m_vHjznEtaSpectReco.back() );
    
    m_vHjznEtaSpectTruth.push_back
      ( new TH2D( Form("h_etaSpectTruth_jz%d", jzn ), 
		  ";#eta_{Truth};#it{p}_{T}^{Truth} [GeV]",
		  m_nVarFwdEtaBins, 0, 1,
		  m_nPtSpectBins,
		  m_ptSpectMin, m_ptSpectMax ) );
    m_vHjznEtaSpectTruth.back()->GetXaxis()->
      Set( m_nVarFwdEtaBins, &( m_varFwdEtaBinning[0] ) );
    AddHistogram( m_vHjznEtaSpectTruth.back() );
    
    // The Nent Histogram is only for the overlay samples
    // where we fill with some weight. Otherwise,
    // SpectTruth = SpectTruthNent (b/c fill with w=1)
    m_vHjznEtaSpectTruthNent.push_back
      ( new TH2D( Form("h_etaSpectTruthNent_jz%d", jzn ), 
		  ";#eta_{Truth};#it{p}_{T}^{Truth} [GeV]",
		  m_nVarFwdEtaBins, 0, 1,
		  m_nPtSpectBins,
		  m_ptSpectMin, m_ptSpectMax ) );
    m_vHjznEtaSpectTruthNent.back()->GetXaxis()->
      Set( m_nVarFwdEtaBins, &( m_varFwdEtaBinning[0] ) );
    AddHistogram( m_vHjznEtaSpectTruthNent.back() );

    m_vHjznEtaSpectTruthPaired.push_back
      ( new TH2D( Form("h_etaSpectTruthPaired_jz%d", jzn ), 
		  ";#eta_{Truth};#it{p}_{T}^{Truth} [GeV]",
		  m_nVarFwdEtaBins, 0, 1,
		  m_nPtSpectBins,
		  m_ptSpectMin, m_ptSpectMax ) );
    m_vHjznEtaSpectTruthPaired.back()->GetXaxis()->
      Set( m_nVarFwdEtaBins, &( m_varFwdEtaBinning[0] ) );
    AddHistogram( m_vHjznEtaSpectTruthPaired.back() );
    
    // --------- recoTruthRpt ---------
    m_vHjznRecoTruthRpt.push_back
      ( new TH3D( Form("h_recoTruthRpt_jz%d", jzn),
		  ";#eta^{Truth};#it{p}_{T}^{Truth};#it{p}_{T}^{Reco}/#it{p}_{T}^{Truth}",
		  m_nEtaForwardBinsFine,
		  m_etaForwardMin, m_etaForwardMax,
		  m_nPtTruthBins,
		  m_ptTruthMin, m_ptTruthMax,
		  m_nRPtRecoTruthBins,
		  m_rPtRecoTruthMin, m_rPtRecoTruthMax) );
    AddHistogram( m_vHjznRecoTruthRpt.back() );
    
    m_vHjznRecoTruthRptNent.push_back
      ( new TH2D( Form("h_recoTruthRptNent_jz%d", jzn),
		  ";#eta^{Truth};#it{p}_{T}^{Truth}",
		  m_nEtaForwardBinsFine,
		  m_etaForwardMin, m_etaForwardMax,
		  m_nPtTruthBins,
		  m_ptTruthMin, m_ptTruthMax ) );
    AddHistogram( m_vHjznRecoTruthRptNent.back() );    
    
    // --------- recoTruthDeta ---------
    m_vHjznRecoTruthDeta.push_back
      ( new TH3D( Form("h_recoTruthDeta_jz%d", jzn),
		  ";#eta^{Truth};#it{p}_{T}^{Truth};#eta^{Reco}-#eta^{Truth}",
		  m_nEtaForwardBinsFine,
		  m_etaForwardMin, m_etaForwardMax,
		  m_nPtTruthBins,
		  m_ptTruthMin, m_ptTruthMax,
		  m_nDAngleRecoTruthBins,
		  m_dAngleRecoTruthMin, m_dAngleRecoTruthMax ) );
    AddHistogram( m_vHjznRecoTruthDeta.back() );    
    
    m_vHjznRecoTruthDetaNent.push_back
      ( new TH2D( Form("h_recoTruthDetaNent_jz%d", jzn),
		  ";#eta^{Truth};#it{p}_{T}^{Truth}",
		  m_nEtaForwardBinsFine,
		  m_etaForwardMin, m_etaForwardMax,
		  m_nPtTruthBins,
		  m_ptTruthMin, m_ptTruthMax ) );
    AddHistogram( m_vHjznRecoTruthDetaNent.back() );    

    // --------- recoTruthDphi ---------
    m_vHjznRecoTruthDphi.push_back
      ( new TH3D( Form("h_recoTruthDphi_jz%d", jzn),
		  ";#eta^{Truth};#it{p}_{T}^{Truth};#phi^{Reco}-#phi^{Truth}",
		  m_nEtaForwardBinsFine,
		  m_etaForwardMin, m_etaForwardMax,
		  m_nPtTruthBins,
		  m_ptTruthMin, m_ptTruthMax,
		  m_nDAngleRecoTruthBins,
		  m_dAngleRecoTruthMin, m_dAngleRecoTruthMax ) );
    AddHistogram( m_vHjznRecoTruthDphi.back() );
    
    m_vHjznRecoTruthDphiNent.push_back
      ( new TH2D( Form("h_recoTruthDphiNent_jz%d", jzn),
		  ";#eta^{Truth};#it{p}_{T}^{Truth}",
		  m_nEtaForwardBinsFine,
		  m_etaForwardMin, m_etaForwardMax,
		  m_nPtTruthBins,
		  m_ptTruthMin, m_ptTruthMax ) );
    AddHistogram( m_vHjznRecoTruthDphiNent.back() );


     // -------- dPhi --------
    std::vector< int    > nDphiBins;
    std::vector< double > dPhiMin;
    std::vector< double > dPhiMax;

    boost::assign::push_back( nDphiBins )
      ( m_nVarEtaBins   )( m_nVarEtaBins )
      ( m_nDphiPtBins   )( m_nDphiPtBins )
      ( m_nDphiDphiBins );
    
    boost::assign::push_back( dPhiMin  )
      (        0       )(       0      )
      ( m_nDphiPtMin   )( m_nDphiPtMin )
      ( m_nDphiDphiMin );

    boost::assign::push_back( dPhiMax  )
      (        1       )(       1      )
      ( m_nDphiPtMax   )( m_nDphiPtMax )
      ( m_nDphiDphiMax );

    uint nDim = nDphiBins.size();
    
    THnSparse* hnReco =
      new THnSparseD( Form("hn_dPhiReco_jz%d", jzn ), "",
		      nDim, &nDphiBins[0], &dPhiMin[0], &dPhiMax[0] );
    hnReco->GetAxis(0)->Set( m_nVarEtaBins, &( m_varEtaBinning[0] ) );
    hnReco->GetAxis(1)->Set( m_nVarEtaBins, &( m_varEtaBinning[0] ) );
    hnReco->GetAxis(2)->Set( m_nVarPtBins , &( m_varPtBinning[0]  ) );
    hnReco->GetAxis(3)->Set( m_nVarPtBins , &( m_varPtBinning[0]  ) );

    m_vHjznDphiReco.push_back( hnReco );
    AddHistogram( hnReco );

    THnSparse* hnTruth =
      new THnSparseD( Form("hn_dPhiTruth_jz%d", jzn ), "",
		      nDim, &nDphiBins[0], &dPhiMin[0], &dPhiMax[0] );
    hnTruth->GetAxis(0)->Set( m_nVarEtaBins, &( m_varEtaBinning[0] ) );
    hnTruth->GetAxis(1)->Set( m_nVarEtaBins, &( m_varEtaBinning[0] ) );
    hnTruth->GetAxis(2)->Set( m_nVarPtBins , &( m_varPtBinning[0]  ) );
    hnTruth->GetAxis(3)->Set( m_nVarPtBins , &( m_varPtBinning[0]  ) );

    m_vHjznDphiTruth.push_back( hnTruth );
    AddHistogram( hnTruth );

    
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
  
  for( uint iG = 0; iG < m_nJzn; iG++){
     
    std::cout << "fNameIn: " << m_vJznFnameIn[iG] << std::endl;
  
    m_fIn  = TFile::Open( m_vJznFnameIn[iG].c_str() );
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

      // Do Dphi analysis
      AnalyzeDeltaPhi( m_vHjznDphiReco [iG], vR_jets );
      AnalyzeDeltaPhi( m_vHjznDphiTruth[iG], vT_jets );
      
      // fill for truth jets
      // denominator for efficiency because not all
      // reco jets are reconstructed for a truth jet
      // also count total, and fwd truth jets
      for( auto& tJet : vT_jets ){
	double jetEta    = tJet.Eta();
	double jetEtaAdj = anaTool->AdjustEtaForPP( jetEta, m_is_pPb );
	double jetPhi    = tJet.Phi();
	double jetPt     = tJet.Pt()/1000.;

	double weight = GetJetWeight( jetEta, jetPhi, jetPt );

	m_vHjznEtaSpectTruth    [iG]->
	  Fill( jetEtaAdj, jetPt, weight);

	m_vHjznEtaSpectTruthNent[iG]->
	  Fill( jetEtaAdj, jetPt );
		
	// count how many total truth jets
	// and how many forward truth jets
	// cut on pt corresponding to jz sample
	if( jetPt > m_vJznPtThreshold[iG] ){
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
	double jetEtaRecoAdj  = anaTool->AdjustEtaForPP( jetEtaReco , m_is_pPb );
	double jetEtaTruthAdj = anaTool->AdjustEtaForPP( jetEtaTruth, m_is_pPb  );	  
	
	double weight = GetJetWeight( jetEtaTruth, jetPhiTruth, jetPtTruth );
	
	m_vHjznEtaPhiMap[iG]->Fill( jetEtaReco, jetPhiReco, weight);
	m_vHjznEtaPtMap [iG]->Fill( jetEtaReco, jetPtReco , weight);
        		
	if( vp.DeltaR() <= m_dRmax ){
	  m_vHjznEtaSpectReco  [iG]->
	    Fill( jetEtaRecoAdj,  jetPtReco,  weight);

	  m_vHjznEtaSpectTruthPaired[iG]->
	    Fill( jetEtaTruthAdj,  jetPtTruth,  weight);

	  m_vHjznRecoTruthRpt     [iG]->
	    Fill( jetEtaTruthAdj, jetPtTruth, jetPtReco/jetPtTruth, weight);
	  m_vHjznRecoTruthRptNent [iG]->
	    Fill( jetEtaTruthAdj, jetPtTruth );

	  m_vHjznRecoTruthDeta    [iG]->
	    Fill( jetEtaTruthAdj, jetPtTruth, jetEtaReco - jetEtaTruth, weight);
	  m_vHjznRecoTruthDetaNent[iG]->
	    Fill( jetEtaTruthAdj, jetPtTruth );

	  m_vHjznRecoTruthDphi    [iG]->
	    Fill( jetEtaTruthAdj, jetPtTruth, jetPhiReco - jetPhiTruth, weight );
	  m_vHjznRecoTruthDphiNent[iG]->
	    Fill( jetEtaTruthAdj, jetPtTruth  );
	}
      } // end loop over pairs
    } // end loop over events
   
    std::cout << "DONE WITH jz"      << m_vJznUsed[iG] 
	      << "   nJetsForward: " << nJetsForward
	      << "   nTotalJets: "   << nJetsTotal
	      << "   nTotalEvents: " << nEventsTotal
	      << std::endl;

    m_fIn->Close();
  } // end loop over a JZ sample
}

//---------------------------------
//           Plot Data
//---------------------------------
void DiJetAnalysisMC::LoadHistograms(){
  m_fIn = TFile::Open( m_rootFname.c_str() ); 
  
  for( uint iG = 0; iG < m_nJzn; iG++){
    int jzn = m_vJznUsed[iG];
    
    // -------- maps ---------
    m_vHjznEtaPhiMap.push_back
      ( static_cast< TH2D* >
	( m_fIn->Get( Form("h_etaPhiMap_jz%d", jzn ))));
    m_vHjznEtaPhiMap.back()->SetDirectory(0);

    m_vHjznEtaPtMap.push_back
      ( static_cast< TH2D* >
	( m_fIn->Get( Form("h_etaPtMap_jz%d", jzn ))));
    m_vHjznEtaPtMap.back()->SetDirectory(0);

    // -------- spect --------
    m_vHjznEtaSpectReco.push_back
      ( static_cast< TH2D* >
	( m_fIn->Get( Form("h_etaSpectReco_jz%d", jzn ))));
    m_vHjznEtaSpectReco.back()->SetDirectory(0);

    // The Nent Histogram is only for the overlay samples
    // where we fill with some weight. Otherwise,
    // SpectTruth = SpectTruthNent (b/c fill with w=1)
    m_vHjznEtaSpectTruth.push_back 
      ( static_cast< TH2D* >
	( m_fIn->Get( Form("h_etaSpectTruth_jz%d", jzn ))));
    m_vHjznEtaSpectTruth.back()->SetDirectory(0);
    
    m_vHjznEtaSpectTruthNent.push_back
      ( static_cast< TH2D* >
	( m_fIn->Get( Form("h_etaSpectTruthNent_jz%d", jzn ))));
    m_vHjznEtaSpectTruthNent.back()->SetDirectory(0);

    m_vHjznEtaSpectTruthPaired.push_back
      ( static_cast< TH2D* >
	( m_fIn->Get( Form("h_etaSpectTruthPaired_jz%d", jzn ))));
    m_vHjznEtaSpectTruthPaired.back()->SetDirectory(0);

    // --------- recoTruthRpt ---------
    m_vHjznRecoTruthRpt.push_back
      ( static_cast< TH3D* >
	( m_fIn->Get( Form("h_recoTruthRpt_jz%d", jzn ))));
    m_vHjznRecoTruthRpt.back()->SetDirectory(0);

    m_vHjznRecoTruthRptNent.push_back
      ( static_cast< TH2D* >
	( m_fIn->Get( Form("h_recoTruthRptNent_jz%d", jzn ))));
    m_vHjznRecoTruthRptNent.back()->SetDirectory(0);
    
    // --------- recoTruthDeta ---------
    m_vHjznRecoTruthDeta.push_back
      ( static_cast< TH3D* >
	( m_fIn->Get( Form("h_recoTruthDeta_jz%d", jzn ))));
    m_vHjznRecoTruthDeta.back()->SetDirectory(0);

    m_vHjznRecoTruthDetaNent.push_back
      ( static_cast< TH2D* >
	( m_fIn->Get( Form("h_recoTruthDetaNent_jz%d", jzn ))));
    m_vHjznRecoTruthDetaNent.back()->SetDirectory(0);
    
    // --------- recoTruthDphi ---------
    m_vHjznRecoTruthDphi.push_back
      ( static_cast< TH3D* >
	( m_fIn->Get( Form("h_recoTruthDphi_jz%d", jzn )))); 
    m_vHjznRecoTruthDphi.back()->SetDirectory(0);

    m_vHjznRecoTruthDphiNent.push_back
      ( static_cast< TH2D* >
	( m_fIn->Get( Form("h_recoTruthDphiNent_jz%d", jzn ))));
    m_vHjznRecoTruthDphiNent.back()->SetDirectory(0);

    // -------- dPhi- --------
    m_vHjznDphiReco.push_back
      ( static_cast< THnSparse *>
	( m_fIn->Get( Form("hn_dPhiReco_jz%d", jzn ))));  

    m_vHjznDphiTruth.push_back
      ( static_cast< THnSparse *>
	( m_fIn->Get( Form("hn_dPhiTruth_jz%d", jzn ))));  
}

  m_fIn->Close();
}

void DiJetAnalysisMC::PlotSpectra( std::vector< TH2*>& vJznSpect,
				   const std::string& type,
				   const std::string& level ){
  if( !vJznSpect.size() ){ return; }

  std::string yAxisTitle = "dN/d#it{p}_{T}";
  
  // these are same for all the other histos
  int nBinsX = vJznSpect.front()->GetNbinsX();

  int  nBinsY = vJznSpect.front()->GetNbinsY();
  double yMin = vJznSpect.front()->GetYaxis()->GetXmin();
  double yMax = vJznSpect.front()->GetYaxis()->GetXmax();

  double ptSpectWidth = ( yMax - yMin ) / nBinsY;
  
  // vector of projected spect
  std::vector< std::vector< TH1* > > vSpect;
  vSpect.resize( m_nJzn );
  // vector of final spectra
  std::vector< TH1* > vSpectFinal;
  
  for( uint iG = 0; iG < m_nJzn; iG++){
    int jzn = m_vJznUsed[iG];
    TH2* hJznSpect = vJznSpect[iG];
   
    for( int xBin = 1; xBin <= nBinsX; xBin++ ){
      double xBinMin = hJznSpect->GetXaxis()->GetBinLowEdge( xBin );
      double xBinMax = hJznSpect->GetXaxis()->GetBinUpEdge ( xBin );

      TH1* hSpect =
	hJznSpect->ProjectionY
	( Form("h_%s_%s_%s_jz%d",
	       type.c_str(),
	       level.c_str(),
	       anaTool->GetName
	       ( xBinMin, xBinMax, "Eta" ).c_str(),
	       jzn ),
	  xBin, xBin );

      hSpect->SetTitle
	( anaTool->GetEtaLabel( xBinMin, xBinMax).c_str() );
      hSpect->GetYaxis()->SetTitle( yAxisTitle.c_str() );
      hSpect->Scale( 1./ptSpectWidth );
      vSpect[iG].push_back( hSpect );      
    } // end loop over xAxis (eta)
  } // end loop over jz samples

  for( int iX = 0; iX < nBinsX; iX++ ){
    TH1* hSpect = vSpect[0][iX]; 

    double xBinMin = hSpect->GetXaxis()->GetBinLowEdge( iX + 1 );
    double xBinMax = hSpect->GetXaxis()->GetBinUpEdge ( iX + 1 );
    
    TH1* hSpectFinal =
      new TH1D( Form("h_%s_%s_%s_final",
		     type.c_str(),
		     level.c_str(),
		     anaTool->GetName(xBinMin, xBinMax,"Eta").c_str() ),
		Form("%s;%s;%s",
		     hSpect->GetTitle(),
		     hSpect->GetXaxis()->GetTitle(),
		     hSpect->GetYaxis()->GetTitle() ),
		nBinsY, yMin, yMax );
    vSpectFinal.push_back( hSpectFinal );

    std::vector< TH1* > vEtaTemp;
    for( uint iG = 0; iG < m_nJzn; iG++){
      vEtaTemp.push_back( vSpect[iG][iX] );
    }
    
    CombineJZN( hSpectFinal, vEtaTemp );
  } // end loop over eta

  bool isLog = true;
  DrawCanvas( vSpectFinal, type, level, isLog );

  for( uint iG = 0; iG < m_nJzn; iG++){
    delete vSpectFinal[iG];
    for( int iX = 0; iX < nBinsX; iX++ ){
      delete vSpect[iG][iX];
    }
  }
}

void DiJetAnalysisMC::PlotVsEtaPt( std::vector< TH3* >& vJznHin,
				   std::vector< TH2* >& vJznNentIn,
				   const std::string& type ){
  // check if we have one
  if( !vJznHin.size() ){ return; }

  // is same for all the other histos
  int nBinsX = vJznHin.front()->GetNbinsX();

  int  nBinsY = vJznHin.front()->GetNbinsY();
  double yMin = vJznHin.front()->GetYaxis()->GetXmin();
  double yMax = vJznHin.front()->GetYaxis()->GetXmax();

  // to not make spelling errors later
  std::string sMean  = "mean";
  std::string sSigma = "sigma";

  std::string yTitleMean;
  std::string yTitleSigma;
  GetTypeTitle( type, yTitleMean, yTitleSigma );

  std::vector< std::vector< std::vector< TH1* > > > vProj;
  std::vector< std::vector< std::vector< TF1* > > > vFit;
  std::vector< std::vector< TH1* > > vMeans;
  std::vector< std::vector< TH1* > > vSigmas;
  std::vector< std::vector< TH1* > > vNent;
  vProj  .resize( m_nJzn );
  vFit   .resize( m_nJzn );
  vMeans .resize( m_nJzn );
  vSigmas.resize( m_nJzn );
  vNent  .resize( m_nJzn );
  
  // Get The mean, sigma histograms ready.
  // Take projection onto PT axis for Nent.
  for( uint iG = 0; iG < m_nJzn; iG++ ){
    int jzn      = m_vJznUsed[iG];
    TH3* hJznHin = vJznHin[iG];
   
    for( int xBin = 1; xBin <= nBinsX; xBin++ ){
      double xBinMin = hJznHin->GetXaxis()->GetBinLowEdge( xBin );
      double xBinMax = hJznHin->GetXaxis()->GetBinUpEdge ( xBin );

      std::string xAxisTitle = hJznHin->GetYaxis()->GetTitle();

      // build mean, sigma, project nev
      TH1* hMean = new TH1D
	( Form("h_%s_%s_%s_jz%d",
	       type.c_str(),
	       anaTool->GetName(xBinMin, xBinMax, "Eta").c_str(),
	       sMean.c_str(),
	       jzn ),
	  Form("%s;%s;%s",
	       anaTool->GetEtaLabel( xBinMin, xBinMax ).c_str(),
	       xAxisTitle.c_str(),
	       yTitleMean.c_str() ),
	  nBinsY, yMin, yMax );
      vMeans[iG].push_back( hMean );
      
      TH1* hSigma = new TH1D
	( Form("h_%s_%s_%s_jz%d",
	       type.c_str(),
	       anaTool->GetName(xBinMin, xBinMax, "Eta").c_str(),
	       sSigma.c_str(),
	       jzn ),
	  Form("%s;%s;%s",
	       anaTool->GetEtaLabel( xBinMin, xBinMax ).c_str(),
	       xAxisTitle.c_str(),
	       yTitleSigma.c_str() ),
	  nBinsY, yMin, yMax );
      vSigmas[iG].push_back( hSigma );
      
      TH1* hNent = vJznNentIn[iG]->ProjectionY
	( Form("h_%s_%s_N_jz%d",
	       type.c_str(),
	       anaTool->GetName(xBinMin, xBinMax, "Eta").c_str(),
	       jzn ), xBin, xBin );
      vNent[iG].push_back( hNent );

      // now loop over the second axis and project onto z axis
      // to get the gaussian distributions. fit them.
      vProj[iG].resize( nBinsX );
      vFit [iG].resize( nBinsX );
      
      for( int yBin = 1; yBin <= nBinsY; yBin++ ){
	double yBinMin = hJznHin->GetYaxis()->GetBinLowEdge(yBin);
	double yBinMax = hJznHin->GetYaxis()->GetBinUpEdge (yBin);

	int iX = xBin - 1;

	TH1* hProj = hJznHin->ProjectionZ
	  ( Form("h_%s_%s_%s_jn%d",  type.c_str(),
		 anaTool->GetName(xBinMin, xBinMax, "Eta").c_str(),
		 anaTool->GetName(yBinMin, yBinMax, "Pt").c_str(),
		 jzn ),
	    xBin, xBin, yBin, yBin );

	styleTool->SetHStyle( hProj, 0 );
	vProj[iG][iX].push_back( hProj );
	
 	TF1* fit = anaTool->FitGaussian( hProj );
	styleTool->SetHStyle( fit, 0 );
	vFit[iG][iX].push_back( fit );

	TCanvas c( "c", hProj->GetName(), 800, 600 );
	
	hProj->SetTitle("");
	hProj->Draw();

	fit->Draw("same");
	  
	if( !m_isData ){
	  drawTool->DrawAtlasInternalMCRight( 0, 0, m_mcTypeLabel );
	  drawTool->DrawLeftLatex( 0.18, 0.74, Form("JZ%d", jzn) );
	} else if( m_isData ){
	  drawTool->DrawAtlasInternalDataRight( 0, 0, m_is_pPb ); 
	}
    
	drawTool->DrawLeftLatex
	  ( 0.18, 0.88, anaTool->GetEtaLabel( xBinMin, xBinMax ) );
	drawTool->DrawLeftLatex
	  ( 0.18, 0.81, anaTool->GetLabel( yMin, yMax, "#it{p}_{T}^{Truth}") );

	SaveAsROOT( c, hProj->GetName() );

	// if the fit is "bad" do not set any entries on
	// final histograms
	if( hProj->GetEntries() < m_nMinEntriesGausFit ){ continue; }

	hMean->SetBinContent( yBin, fit->GetParameter(1) );
	hMean->SetBinError  ( yBin, fit->GetParError (1) );

	hSigma->SetBinContent ( yBin, fit->GetParameter(2) );
	hSigma->SetBinError   ( yBin, fit->GetParError (2) );
      }
    } // end loop over xBin (eta)
  }// end loop over jzn

  // Now, loop over the sigmas and means, and plot as function of
  // jz sample and eta
  std::vector< TH1* > vMeansFinal;
  std::vector< TH1* > vSigmasFinal;

  for( int iX = 0; iX < nBinsX; iX++ ){

    int    xBin    = iX + 1;
    double xBinMin = vJznHin[0]->GetXaxis()->GetBinLowEdge( xBin );
    double xBinMax = vJznHin[0]->GetXaxis()->GetBinUpEdge ( xBin );

    std::string xAxisTitle = vJznHin[0]->GetYaxis()->GetTitle();
    
    std::vector< TH1* > vEtaMeanTemp;
    std::vector< TH1* > vEtaSigmaTemp;
    std::vector< TH1* > vEtaNentTemp;
    
    for( uint iG = 0; iG < m_nJzn; iG++ ){
      vEtaMeanTemp .push_back( vMeans [iG][iX] );
      vEtaSigmaTemp.push_back( vSigmas[iG][iX] );
      vEtaNentTemp .push_back( vNent  [iG][iX] );
    }

    // build mean, sigma, project nev
    TH1* hMeanFinal = new TH1D
      ( Form("h_%s_%s_%s_final",
	     type.c_str(),
	     anaTool->GetName(xBinMin, xBinMax, "Eta").c_str(),
	     sMean.c_str() ),
	Form("%s;%s;%s",
	     anaTool->GetEtaLabel( xBinMin, xBinMax ).c_str(),
	     xAxisTitle.c_str(),
	     yTitleMean.c_str() ),
	nBinsY, yMin, yMax );
    vMeansFinal.push_back( hMeanFinal );
      
    TH1* hSigmaFinal = new TH1D
      ( Form("h_%s_%s_%s_final",
	     type.c_str(),
	     anaTool->GetName(xBinMin, xBinMax, "Eta").c_str(),
	     sSigma.c_str() ),
	Form("%s;%s;%s",
	     anaTool->GetEtaLabel( xBinMin, xBinMax ).c_str(),
	     xAxisTitle.c_str(),
	     yTitleSigma.c_str() ),
	nBinsY, yMin, yMax );
    vSigmasFinal.push_back( hSigmaFinal );

    CombineJZN( hMeanFinal , vEtaMeanTemp , vEtaNentTemp  );
    CombineJZN( hSigmaFinal, vEtaSigmaTemp, vEtaNentTemp  );
  }

  DrawCanvas( vMeansFinal , type, sMean , 3 );
  DrawCanvas( vSigmasFinal, type, sSigma, 3 );

  for( uint iG = 0; iG < m_nJzn; iG++){
    delete vMeansFinal [iG];
    delete vSigmasFinal[iG];
    for( int iX = 0; iX < nBinsX; iX++ ){
      delete vMeans[iG][iX];
      delete vSigmas[iG][iX];
      delete vNent[iG][iX];
      for( int iY = 0; iY < nBinsY; iY++ ){
	delete vProj[iG][iX][iY];
	delete vFit [iG][iX][iY];
      }
    }
  }
}

void DiJetAnalysisMC::PlotDphiTogether(){

  std::string particles = m_is_pPb ? "pPb" : "pp";
  
  TFile* fIn  = TFile::Open( Form("output/output_%s/c_myOut_%s.root",
				  m_labelOut.c_str(), m_labelOut.c_str() ) );
  TFile* fOut = new TFile("output/c_myOut_mc.root","recreate");

  for( uint iG = 0; iG < m_nJzn; iG++ ){
    std::string jznLabel = m_vJznLabel[iG];
    
    for( uint eta1Bin = 0; eta1Bin < m_nVarEtaBins; eta1Bin++ ){
      for( uint eta2Bin = 0; eta2Bin < m_nVarEtaBins; eta2Bin++ ){
	for( uint pt1Bin = 0; pt1Bin < m_nVarPtBins; pt1Bin++ ){
	  for( uint pt2Bin = 0; pt2Bin < m_nVarPtBins; pt2Bin++ ){
	  
	    double eta1Low = m_varEtaBinning[ eta1Bin ];
	    double eta1Up  = m_varEtaBinning[ eta1Bin + 1 ];
	    double eta1Center = eta1Low + 0.5 * ( eta1Up - eta1Low );
	  
	    double eta2Low = m_varEtaBinning[ eta2Bin ];
	    double eta2Up  = m_varEtaBinning[ eta2Bin + 1 ];
	    double eta2Center = eta2Low + 0.5 * ( eta2Up - eta2Low );

	    double pt1Low  = m_varPtBinning[ pt1Bin ];
	    double pt1Up   = m_varPtBinning[ pt1Bin + 1 ];

	    double pt2Low  = m_varPtBinning[ pt2Bin ];


	    if( !anaTool->IsForward( eta1Center ) &&
		!anaTool->IsForward( eta2Center ) )
	      { continue; }
	    if( pt1Low < pt2Low  )
	      { continue; }
	    
	    std::string hTag =
	      Form("%s_%s_%s_%s",
		   anaTool->GetName( eta1Low, eta1Up, "Eta1").c_str(),
		   anaTool->GetName( eta2Low, eta2Up, "Eta2").c_str(),
		   anaTool->GetName( pt1Low , pt1Up , "Pt1" ).c_str(),
		   anaTool->GetName( pt2Low , pt2Low, "Pt2" ).c_str() );

	    std::string hName_reco =
	      Form("h_dPhi_reco_%s_%s", hTag.c_str(), jznLabel.c_str() );
	    std::string hName_truth =
	      Form("h_dPhi_truth_%s_%s", hTag.c_str(), jznLabel.c_str() );
	    
	    TCanvas* c_reco =
	      static_cast<TCanvas*>
	      ( fIn->Get( Form("c_%s_%s", hName_reco.c_str(), m_labelOut.c_str())));
	    TCanvas* c_truth  =
	      static_cast<TCanvas*>
	      ( fIn->Get( Form("c_%s_%s", hName_truth.c_str(), m_labelOut.c_str())));

	    TH1* h_reco =
	      static_cast<TH1D*>( c_reco->GetPrimitive( hName_reco.c_str() ) );
	    styleTool->SetHStyle( h_reco, 0 );
	    TH1* h_truth  =
	      static_cast<TH1D*>( c_truth->GetPrimitive( hName_truth.c_str() ) );
	    styleTool->SetHStyle( h_truth, 1 );
	  
	    TF1* f_reco =
	      static_cast<TF1*>( c_reco->GetPrimitive
				 ( Form("f_%s", hName_reco.c_str())));
	    styleTool->SetHStyle( f_reco, 0 );
	    f_reco->SetLineColor( h_reco->GetLineColor() );
	    TF1* f_truth  =
	      static_cast<TF1*>( c_truth->GetPrimitive
				 ( Form("f_%s", hName_truth.c_str())));
	    styleTool->SetHStyle( f_truth, 1 );
	    f_truth->SetLineColor( h_truth->GetLineColor() );
	  
	    TCanvas c("c","c", 800, 600 );
	    
	    TLegend leg( 0.27, 0.41, 0.38, 0.52 );
	    styleTool->SetLegendStyle( &leg );
	    leg.AddEntry( h_reco, "Reco");
	    leg.AddEntry( h_truth , "Truth");

	    h_reco->Draw("epsame");
	    h_truth->Draw("epsame");

	    f_reco->Draw("same");
	    f_truth->Draw("same");

	    leg.Draw("same");

	    if( h_reco->GetMaximum() > h_truth->GetMaximum() ){
	      h_reco->SetMaximum( h_reco->GetMaximum() * 1.1 );
	      h_truth->SetMaximum ( h_reco->GetMaximum() * 1.1 );
	    } else {
	      h_reco->SetMaximum( h_truth->GetMaximum() * 1.1 );
	      h_truth->SetMaximum ( h_truth->GetMaximum() * 1.1 );
	    }
	    
	    drawTool->DrawLeftLatex
	      ( 0.13, 0.87,anaTool->GetLabel( eta1Low, eta1Up, "#eta_{1}" ) );
	    drawTool->DrawLeftLatex
	      ( 0.13, 0.82,anaTool->GetLabel( eta2Low, eta2Up, "#eta_{2}" ) );
	    drawTool->DrawLeftLatex
	      ( 0.13, 0.76,anaTool->GetLabel( pt1Low, pt1Up, "#it{p}_{T}^{1}" ) );
	    drawTool->DrawLeftLatex
	      ( 0.13, 0.69,anaTool->GetLabel( pt2Low, pt2Low, "#it{p}_{T}^{2}" ) );

	    drawTool->DrawAtlasInternal();
	  
	    SaveAsROOT( c, Form("h_dPhi_%s", hTag.c_str() ) );

	    delete h_reco;
	    delete h_truth;
	    delete f_reco;
	    delete f_truth;
	    delete c_reco;
	    delete c_truth;
	  }
	}
      } 
    }
  }
  
  std::cout << "DONE! Closing " << fOut->GetName() << std::endl;
  fOut->Close();
  std::cout << "......Closed  " << fOut->GetName() << std::endl;
}

void DiJetAnalysisMC::PlotEtaPhiPtMap( std::vector< TH2* >& vJznHin ){
  TCanvas c_map("c_map","c_map",800,600);

  for( uint iG = 0; iG < m_nJzn; iG++ ){
    vJznHin[iG]->Draw("col");
    styleTool->SetHStyle( vJznHin[iG], 0 );
    drawTool->DrawAtlasInternalMCLeft( 0, -0.55, m_mcTypeLabel  );  
    SaveAsPdfPng( c_map, vJznHin[iG]->GetName() );
  }
}

//---------------------------
//          Tools 
//---------------------------

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
 
void DiJetAnalysisMC::CombineJZN( TH1* h_res,
				  std::vector< TH1* >& vJznHin){
  for( uint iG = 0; iG < m_nJzn; iG++ ){
    double scale = m_vJznEff[iG] * m_vJznSigma[iG] /
      m_vJznSumPowhegWeights[iG];
    h_res->Add( vJznHin[iG], scale );
  }
  h_res->Scale( 1./m_sumSigmaEff );
}

void DiJetAnalysisMC::CombineJZN( TH1* h_res,
				  std::vector< TH1*>& vJznVIN,
				  std::vector< TH1*>& vJznNentIN ){
  for( int xBin = 1; xBin <= h_res->GetNbinsX(); xBin++ ){

    double valTot      = 0;
    double valErrorTot = 0;
    double denomTot    = 0;

    for( uint iG = 0; iG < m_nJzn; iG++ ){
      double valueBin    = vJznVIN      [iG]->GetBinContent( xBin );
      double valueBinErr = vJznVIN      [iG]->GetBinError( xBin );
      double nEntriesBin = vJznNentIN   [iG]->GetBinContent( xBin );
      double weight      = m_vJznWeights[iG] / m_vJznSumPowhegWeights[iG];

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

double DiJetAnalysisMC::GetJetWeight( double eta, double phi, double pt ){
  if( m_mcType != 0 ) { return 1; } // only for powheg
  
  int xb = m_hPowhegWeights->GetXaxis()->FindBin(pt);
  int yb = m_hPowhegWeights->GetYaxis()->FindBin(eta);
  int zb = m_hPowhegWeights->GetZaxis()->FindBin(phi);
  float jet_weight = m_hPowhegWeights->GetBinContent(xb,yb,zb);

  return jet_weight;
}

void DiJetAnalysisMC::GetTypeTitle( const std::string& type,
				    std::string& yTitleMean,
				    std::string& yTitleSigma ){ 
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
}
