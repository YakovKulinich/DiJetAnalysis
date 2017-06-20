#include <TROOT.h>
#include <TEnv.h>
#include <TGraphAsymmErrors.h>
#include <TLorentzVector.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TF1.h>
#include <THnSparse.h>

#include <iostream>
#include <sstream>
#include <cmath>
#include <cfloat>
#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <boost/assign.hpp>

#include "MyRoot.h"

#include "DiJetAnalysisMC.h"
#include "DeltaPhiProj.h"

TH3* DiJetAnalysisMC::m_hPowhegWeights = NULL;

DiJetAnalysisMC::DiJetAnalysisMC() : DiJetAnalysisMC( false, 0 )
{}

DiJetAnalysisMC::DiJetAnalysisMC( bool is_pPb, int mcType )
  : DiJetAnalysis( is_pPb, false, mcType )
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

  //=============== Histo Names ==================    
  m_etaSpectRecoName  = "etaSpect_reco";
  m_etaSpectTruthName = "etaSpect_truth";

  m_dPhiRecoName      = "dPhi_reco";
  m_dPhiTruthName     = "dPhi_truth";
}

DiJetAnalysisMC::~DiJetAnalysisMC(){}

void DiJetAnalysisMC::Initialize(){
  m_labelOut = m_is_pPb ? "pPb_mc" : "pp_mc" ;
  
  std::string mcMenu;
  if     ( m_mcType == 0 ){ mcMenu = "2015.pp.pythia8powheg"; }
  else if( m_mcType == 1 ){ mcMenu = "2015.pp.pythia8"      ; }
  else if( m_mcType == 2 ){ mcMenu = "2015.pp.herwig"       ; }
  
  std::string labelOut =
    GetConfig()->GetValue( Form("labelOut.%s",mcMenu.c_str()),"");
  m_labelOut += "_" + labelOut;
  m_mcTypeLabel = GetConfig()->GetValue( Form("mcTypeLabel.%s",mcMenu.c_str()),"");

  std::vector< double > vJznUsedD = anaTool->vectoriseD
    ( GetConfig()->GetValue( Form("jznUsed.%s",mcMenu.c_str()),"" )," ");
  // double -> int 
  for( auto d : vJznUsedD ){ m_vJznUsed.push_back(d); }
  m_vJznLabel   =  anaTool->vectorise
    ( GetConfig()->GetValue( Form("jznLabel.%s"  ,mcMenu.c_str()),"" )," ");
  m_vJznFnameIn =  anaTool->vectorise
    ( GetConfig()->GetValue( Form("jznFnameIn.%s",mcMenu.c_str()),"" )," ");
  m_vJznSigma   = anaTool->vectoriseD
    ( GetConfig()->GetValue( Form("jznSigma.%s"  ,mcMenu.c_str()),"" )," ");
  m_vJznEff     = anaTool->vectoriseD
    ( GetConfig()->GetValue( Form("jznEff.%s"    ,mcMenu.c_str()),"" )," ");
  m_vJznSumOverlayWeights = anaTool->vectoriseD
    ( GetConfig()->GetValue( Form("jznSumOverlayWeights.%s",mcMenu.c_str()),"" )," ");
  m_vJznPtThreshold       = anaTool->vectoriseD
    ( GetConfig()->GetValue( Form("jznPtThreshold.%s"      ,mcMenu.c_str()),"" )," ");

  std::string weightFile =
    GetConfig()->GetValue( Form("weightFile.%s",mcMenu.c_str()), "");
  if( !weightFile.empty() ){
    // Get Associated weights for powheg.
    // These need to be used at filling time (per event basis).
    TFile* file =
      TFile::Open( weightFile.c_str() );
    if( !m_hPowhegWeights ) // because its a static member 
      { m_hPowhegWeights = static_cast< TH3D* >( file->Get("h3_pT_y_phi_rw") ); }
    m_hPowhegWeights->SetDirectory(0);
    file->Close();
  }

  m_nJzn = m_vJznUsed.size();  
  
  // Check if the directories exist.
  // If they don't, create them
  m_dirOut   = "output";
  anaTool->CheckWriteDir( m_dirOut.c_str() );
  m_dirOut   += "/output_" + m_labelOut;
  anaTool->CheckWriteDir( m_dirOut.c_str() );

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
  
  // PlotEtaPhiPtMap( m_vHjznEtaPhiMap );
  // PlotEtaPhiPtMap( m_vHjznEtaPtMap  );

  // add a slice "all" to collection
  // rest of plots are for combined slices
  m_vJznLabel.push_back( m_allName );
  
  m_hAllEtaSpectReco  = CombineSamples( m_vHjznEtaSpectReco , m_etaSpectRecoName  );
  m_hAllEtaSpectTruth = CombineSamples( m_vHjznEtaSpectTruth, m_etaSpectTruthName );
  PlotSpectra( m_vHjznEtaSpectReco , m_vJznLabel, m_etaSpectRecoName );
  PlotSpectra( m_vHjznEtaSpectTruth, m_vJznLabel, m_etaSpectTruthName  );


  for( auto & h : m_vHjznDphiReco ){ std::cout << h->GetName() << std::endl; }
  for( auto & h : m_vHjznDphiTruth ){ std::cout << h->GetName() << std::endl; }
  
  
  m_hAllDphiReco  = CombineSamples( m_vHjznDphiReco , m_dPhiRecoName  );
  m_hAllDphiTruth = CombineSamples( m_vHjznDphiTruth, m_dPhiTruthName );
  PlotDeltaPhi( m_vHjznDphiReco , m_vJznLabel, m_dPhiRecoName  );
  PlotDeltaPhi( m_vHjznDphiTruth, m_vJznLabel, m_dPhiTruthName );
  
  PlotVsEtaPt( m_vHjznRecoTruthRpt , m_vHjznRecoTruthRptNent , "recoTruthRpt"  );
  PlotVsEtaPt( m_vHjznRecoTruthDeta, m_vHjznRecoTruthDetaNent, "recoTruthDeta" );
  PlotVsEtaPt( m_vHjznRecoTruthDphi, m_vHjznRecoTruthDphiNent, "recoTruthDphi" );

  std::cout << "DONE! Closing " << cfNameOut << std::endl;
  m_fOut->Close();
  std::cout << "......Closed  " << cfNameOut << std::endl;
}


//---------------------------------
//            Read Data
//---------------------------------

void DiJetAnalysisMC::SetupHistograms(){
  
  for( uint iG = 0; iG < m_nJzn; iG++){
    std::string jzn = m_vJznLabel[iG];
    
    std::cout << "Making -" << jzn << " histograms " << std::endl;
    
    // -------- maps ---------
    m_vHjznEtaPhiMap.
      push_back( new TH2D( Form("h_etaPhiMap_%s", jzn.c_str() ),
			  ";#eta_{Reco};#phi_{Reco}",
			  m_nEtaMapBins, m_etaMapMin, m_etaMapMax,
			  m_nPhiMapBins, m_phiMapMin, m_phiMapMax ) );
    AddHistogram( m_vHjznEtaPhiMap.back() );
      
    m_vHjznEtaPtMap.
      push_back( new TH2D( Form("h_etaPtMap_%s", jzn.c_str() ),
			   ";#eta_{Reco};#it{p}_{T}^{Reco}",
			   m_nEtaMapBins, m_etaMapMin, m_etaMapMax,
			   m_nPtMapBins , m_ptMapMin , m_ptMapMax ) );
    AddHistogram( m_vHjznEtaPtMap.back() );
    
    // -------- spect --------
    m_vHjznEtaSpectReco.push_back
      ( new TH2D( Form("h_%s_%s", m_etaSpectRecoName.c_str(), jzn.c_str() ), 
		  ";#eta_{Reco};#it{p}_{T}^{Reco} [GeV]",
		  m_nVarFwdEtaBins, 0, 1,
		  m_nPtSpectBins,
		  m_ptSpectMin, m_ptSpectMax ) );
    m_vHjznEtaSpectReco.back()->GetXaxis()->
      Set( m_nVarFwdEtaBins, &( m_varFwdEtaBinning[0] ) );
    AddHistogram( m_vHjznEtaSpectReco.back() );
    
    m_vHjznEtaSpectTruth.push_back
      ( new TH2D( Form("h_%s_%s", m_etaSpectTruthName.c_str(), jzn.c_str() ), 
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
      ( new TH2D( Form("h_%sNent_%s", m_etaSpectTruthName.c_str(), jzn.c_str() ), 
		  ";#eta_{Truth};#it{p}_{T}^{Truth} [GeV]",
		  m_nVarFwdEtaBins, 0, 1,
		  m_nPtSpectBins,
		  m_ptSpectMin, m_ptSpectMax ) );
    m_vHjznEtaSpectTruthNent.back()->GetXaxis()->
      Set( m_nVarFwdEtaBins, &( m_varFwdEtaBinning[0] ) );
    AddHistogram( m_vHjznEtaSpectTruthNent.back() );

    m_vHjznEtaSpectTruthPaired.push_back
      ( new TH2D( Form("h_%sPaired_%s", m_etaSpectTruthName.c_str(), jzn.c_str() ), 
		  ";#eta_{Truth};#it{p}_{T}^{Truth} [GeV]",
		  m_nVarFwdEtaBins, 0, 1,
		  m_nPtSpectBins,
		  m_ptSpectMin, m_ptSpectMax ) );
    m_vHjznEtaSpectTruthPaired.back()->GetXaxis()->
      Set( m_nVarFwdEtaBins, &( m_varFwdEtaBinning[0] ) );
    AddHistogram( m_vHjznEtaSpectTruthPaired.back() );
    
    // --------- recoTruthRpt ---------
    m_vHjznRecoTruthRpt.push_back
      ( new TH3D( Form("h_recoTruthRpt_%s", jzn.c_str() ),
		  ";#eta^{Truth};#it{p}_{T}^{Truth};#it{p}_{T}^{Reco}/#it{p}_{T}^{Truth}",
		  m_nEtaForwardBinsFine,
		  m_etaForwardMin, m_etaForwardMax,
		  m_nPtTruthBins,
		  m_ptTruthMin, m_ptTruthMax,
		  m_nRPtRecoTruthBins,
		  m_rPtRecoTruthMin, m_rPtRecoTruthMax) );
    AddHistogram( m_vHjznRecoTruthRpt.back() );
    
    m_vHjznRecoTruthRptNent.push_back
      ( new TH2D( Form("h_recoTruthRptNent_%s", jzn.c_str() ),
		  ";#eta^{Truth};#it{p}_{T}^{Truth}",
		  m_nEtaForwardBinsFine,
		  m_etaForwardMin, m_etaForwardMax,
		  m_nPtTruthBins,
		  m_ptTruthMin, m_ptTruthMax ) );
    AddHistogram( m_vHjznRecoTruthRptNent.back() );    
    
    // --------- recoTruthDeta ---------
    m_vHjznRecoTruthDeta.push_back
      ( new TH3D( Form("h_recoTruthDeta_%s", jzn.c_str() ),
		  ";#eta^{Truth};#it{p}_{T}^{Truth};#eta^{Reco}-#eta^{Truth}",
		  m_nEtaForwardBinsFine,
		  m_etaForwardMin, m_etaForwardMax,
		  m_nPtTruthBins,
		  m_ptTruthMin, m_ptTruthMax,
		  m_nDAngleRecoTruthBins,
		  m_dAngleRecoTruthMin, m_dAngleRecoTruthMax ) );
    AddHistogram( m_vHjznRecoTruthDeta.back() );    
    
    m_vHjznRecoTruthDetaNent.push_back
      ( new TH2D( Form("h_recoTruthDetaNent_%s", jzn.c_str() ),
		  ";#eta^{Truth};#it{p}_{T}^{Truth}",
		  m_nEtaForwardBinsFine,
		  m_etaForwardMin, m_etaForwardMax,
		  m_nPtTruthBins,
		  m_ptTruthMin, m_ptTruthMax ) );
    AddHistogram( m_vHjznRecoTruthDetaNent.back() );    

    // --------- recoTruthDphi ---------
    m_vHjznRecoTruthDphi.push_back
      ( new TH3D( Form("h_recoTruthDphi_%s", jzn.c_str() ),
		  ";#eta^{Truth};#it{p}_{T}^{Truth};#phi^{Reco}-#phi^{Truth}",
		  m_nEtaForwardBinsFine,
		  m_etaForwardMin, m_etaForwardMax,
		  m_nPtTruthBins,
		  m_ptTruthMin, m_ptTruthMax,
		  m_nDAngleRecoTruthBins,
		  m_dAngleRecoTruthMin, m_dAngleRecoTruthMax ) );
    AddHistogram( m_vHjznRecoTruthDphi.back() );
    
    m_vHjznRecoTruthDphiNent.push_back
      ( new TH2D( Form("h_recoTruthDphiNent_%s", jzn.c_str() ),
		  ";#eta^{Truth};#it{p}_{T}^{Truth}",
		  m_nEtaForwardBinsFine,
		  m_etaForwardMin, m_etaForwardMax,
		  m_nPtTruthBins,
		  m_ptTruthMin, m_ptTruthMax ) );
    AddHistogram( m_vHjznRecoTruthDphiNent.back() );


    // -------- dPhi --------
    m_nDphiDim     = m_nDphiBins.size();
     
    THnSparse* hnReco =
      new THnSparseD( Form("h_%s_%s", m_dPhiRecoName.c_str(), jzn.c_str() ), "",
		      m_nDphiDim, &m_nDphiBins[0],
		      &m_dPhiMin[0], &m_dPhiMax[0] );
    hnReco->GetAxis(0)->Set( m_nVarYstarBinsA, &( m_varYstarBinningA[0] ) );
    hnReco->GetAxis(1)->Set( m_nVarYstarBinsB, &( m_varYstarBinningB[0] ) );
    hnReco->GetAxis(2)->Set( m_nVarPtBins , &( m_varPtBinning[0]  ) );
    hnReco->GetAxis(3)->Set( m_nVarPtBins , &( m_varPtBinning[0]  ) );
    m_vHjznDphiReco.push_back( hnReco );
    AddHistogram( hnReco );

    THnSparse* hnRecoNent =
      new THnSparseD( Form("h_%sNent_%s",m_dPhiRecoName.c_str(), jzn.c_str() ), "",
		      m_nDphiDim, &m_nDphiBins[0],
		      &m_dPhiMin[0], &m_dPhiMax[0]);
    hnRecoNent->GetAxis(0)->Set( m_nVarYstarBinsA, &( m_varYstarBinningA[0] ) );
    hnRecoNent->GetAxis(1)->Set( m_nVarYstarBinsB, &( m_varYstarBinningB[0] ) );
    hnRecoNent->GetAxis(2)->Set( m_nVarPtBins , &( m_varPtBinning[0]  ) );
    hnRecoNent->GetAxis(3)->Set( m_nVarPtBins , &( m_varPtBinning[0]  ) );
    m_vHjznDphiRecoNent.push_back( hnRecoNent );
    AddHistogram( hnRecoNent );
    
    THnSparse* hnTruth =
      new THnSparseD( Form("h_%s_%s", m_dPhiTruthName.c_str(), jzn.c_str() ), "",
		      m_nDphiDim, &m_nDphiBins[0],
		      &m_dPhiMin[0], &m_dPhiMax[0] );
    hnTruth->GetAxis(0)->Set( m_nVarYstarBinsA, &( m_varYstarBinningA[0] ) );
    hnTruth->GetAxis(1)->Set( m_nVarYstarBinsB, &( m_varYstarBinningB[0] ) );
    hnTruth->GetAxis(2)->Set( m_nVarPtBins , &( m_varPtBinning[0]  ) );
    hnTruth->GetAxis(3)->Set( m_nVarPtBins , &( m_varPtBinning[0]  ) );
    m_vHjznDphiTruth.push_back( hnTruth );
    AddHistogram( hnTruth );

    THnSparse* hnTruthNent =
      new THnSparseD( Form("h_%sNent_%s", m_dPhiTruthName.c_str(), jzn.c_str() ), "",
		      m_nDphiDim, &m_nDphiBins[0],
		      &m_dPhiMin[0], &m_dPhiMax[0]);
    hnTruthNent->GetAxis(0)->Set( m_nVarYstarBinsA, &( m_varYstarBinningA[0] ) );
    hnTruthNent->GetAxis(1)->Set( m_nVarYstarBinsB, &( m_varYstarBinningB[0] ) );
    hnTruthNent->GetAxis(2)->Set( m_nVarPtBins , &( m_varPtBinning[0]  ) );
    hnTruthNent->GetAxis(3)->Set( m_nVarPtBins , &( m_varPtBinning[0]  ) );
    m_vHjznDphiTruthNent.push_back( hnTruthNent );
    AddHistogram( hnTruthNent );
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

    // event loop

    for( m_ev = startEvent; m_ev < endEvent; m_ev++ ){
      m_tree->GetEntry( m_ev );
            
      if( anaTool->DoPrint(m_ev) ) {
	std::cout << "\nEvent : " << m_ev 
		  << "    has : " << vR_jets.size() << " reco jets"
		  << "    and : " << vT_jets.size() << " truth jets"
		  << std::endl; 
      }

      
      ApplyCleaning ( vR_jets, v_isCleanJet );
      ApplyIsolation( vR_jets, 1.0 );

      std::sort( vR_jets.begin(), vR_jets.end(), anaTool->sortByDecendingPt );
      
      std::vector< TLorentzVector > vRR_paired_jets;
      std::vector< TLorentzVector > vRT_paired_jets;
      PairJets( vR_jets, vT_jets, vRR_paired_jets, vRT_paired_jets ); 
      
      // Do Dphi analysis
      AnalyzeDeltaPhi
	( m_vHjznDphiReco [iG], m_vHjznDphiRecoNent[iG],  vRR_paired_jets, GetJetWeight );
      AnalyzeDeltaPhi
	( m_vHjznDphiTruth[iG], m_vHjznDphiTruthNent[iG], vT_jets, GetJetWeight );

      std::vector< TLorentzVector > vTR_paired_jets;
      std::vector< TLorentzVector > vTT_paired_jets;
      PairJets( vT_jets, vR_jets, vTT_paired_jets, vTR_paired_jets );
      
      // loop over truth jets
      // denominator for efficiency because not all
      // reco jets are reconstructed for a truth jet
      for( auto& tJet : vT_jets ){
	double jetEta    = tJet.Eta();
	double jetEtaAdj = AdjustEtaForPP( jetEta );
	double jetPhi    = tJet.Phi();
	double jetPt     = tJet.Pt()/1000.;

	double weight = GetJetWeight( jetEta, jetPhi, jetPt );

	m_vHjznEtaSpectTruth    [iG]->
	  Fill( jetEtaAdj, jetPt, weight);

	m_vHjznEtaSpectTruthNent[iG]->
	  Fill( jetEtaAdj, jetPt );
      } // end loop over truth jets

      // do JER/JES, angular scales and resolution.
      AnalyzeScaleResolution( vTR_paired_jets, vTT_paired_jets, iG );
    } // end loop over events
   
    std::cout << "DONE WITH jz" << m_vJznUsed[iG] << std::endl;

    m_fIn->Close();
  } // end loop over a JZ sample
}

void DiJetAnalysisMC::AnalyzeScaleResolution( const std::vector< TLorentzVector >& vR_jets,
					      const std::vector< TLorentzVector >& vT_jets,
					      const int iG ){
  // loop over for
  for( uint iJet = 0; iJet < vR_jets.size(); iJet ++ ){ 
    const TLorentzVector& recoJet  = vR_jets[iJet];
    const TLorentzVector& truthJet = vT_jets[iJet];
    
    double jetEtaReco  = recoJet.Eta();
    double jetPhiReco  = recoJet.Phi();
    double  jetPtReco  = recoJet.Pt()/1000.;

    double jetEtaTruth = truthJet.Eta();
    double jetPhiTruth = truthJet.Phi();
    double  jetPtTruth = truthJet.Pt()/1000.;
	  
    // convert positive eta to negative because
    // in pp it doesnt matter since detector is
    // symmetric in eta. i.e. eta < 0 <=> eta > 0
    // so we just take positive etas and put them
    // in the negative bins for pp configuration.
    // this saves some overhead with histograms later
    // our histos run in negative eta
    // (due to pPb configuration)
    // the labels will be taken care of so it is ok
    double jetEtaRecoAdj  = AdjustEtaForPP( jetEtaReco );
    double jetEtaTruthAdj = AdjustEtaForPP( jetEtaTruth);	  
	
    double weight = GetJetWeight( jetEtaTruth, jetPhiTruth, jetPtTruth );
	
    m_vHjznEtaPhiMap[iG]->Fill( jetEtaReco, jetPhiReco, weight);
    m_vHjznEtaPtMap [iG]->Fill( jetEtaReco, jetPtReco , weight);

    double deltaR = anaTool->DeltaR( recoJet, truthJet );
    
    if( deltaR <= m_dRmax ){
      m_vHjznEtaSpectReco       [iG]->
	Fill( jetEtaRecoAdj,  jetPtReco,  weight);

      m_vHjznEtaSpectTruthPaired[iG]->
	Fill( jetEtaTruthAdj,  jetPtTruth,  weight);

      m_vHjznRecoTruthRpt       [iG]->
	Fill( jetEtaTruthAdj, jetPtTruth, jetPtReco/jetPtTruth, weight);
      m_vHjznRecoTruthRptNent   [iG]->
	Fill( jetEtaTruthAdj, jetPtTruth );

      m_vHjznRecoTruthDeta      [iG]->
	Fill( jetEtaTruthAdj, jetPtTruth, jetEtaReco - jetEtaTruth, weight);
      m_vHjznRecoTruthDetaNent  [iG]->
	Fill( jetEtaTruthAdj, jetPtTruth );

      m_vHjznRecoTruthDphi      [iG]->
	Fill( jetEtaTruthAdj, jetPtTruth, jetPhiReco - jetPhiTruth, weight );
      m_vHjznRecoTruthDphiNent  [iG]->
	Fill( jetEtaTruthAdj, jetPtTruth  );
    }
  } // end loop over pairs
}


//---------------------------
//          Tools 
//---------------------------

// pair reco and truth jets with deltaR parameter 
void DiJetAnalysisMC::PairJets( std::vector< TLorentzVector >& vA_jets,
			        std::vector< TLorentzVector >& vB_jets,
				std::vector< TLorentzVector >& vA_paired_jets,
				std::vector< TLorentzVector >& vB_paired_jets ){
  // exit if either is empty
  if( !vA_jets.size() || !vB_jets.size() ){ return; }
  
  for( auto& aJet : vA_jets){
    // for each A jet, need to find closest B jet
    // set deltaRmin to something large, find jet with smallest
    // deltaRmin less than Rmax and pair that to the a jet
    double           deltaRmin = 4 * constants::ETAMAX;
    TLorentzVector* pairedBJet = NULL;
    for( auto& bJet : vB_jets ){
      // "cleaned" jets had their px, py, pz set to 0
      // to avoid changing the vector size, skip if
      // its pt is zero.
      if( bJet.Pt() == 0 ){ continue; }
      double deltaR = anaTool->DeltaR( bJet, aJet );
      if( deltaR <= deltaRmin ) {	
	pairedBJet = &bJet;
	deltaRmin  = deltaR;
      }
    }  // end loop over b jets
    // just for safety. if there is at least
    // one b jet, it should pair to all a jets
    if( !pairedBJet ) continue;

    vA_paired_jets.push_back( aJet );
    vB_paired_jets.push_back( *pairedBJet );
  } // end loop over truth jets
}
 
TH1* DiJetAnalysisMC::CombineSamples( std::vector< TH1* >& vSampleHin,
				      const std::string& name ){
  if( !vSampleHin.size() ){ return NULL; }

  TH1* h_res = static_cast< TH1D* >
    ( vSampleHin[0]->Clone( Form("h_%s_%s", name.c_str() , m_allName.c_str() ) ) );
  h_res->Reset();
  
  for( uint iG = 0; iG < vSampleHin.size(); iG++ ){
    double scale = m_vJznEff[iG] * m_vJznSigma[iG] /
      m_vJznSumOverlayWeights[iG];
    h_res->Add( vSampleHin[iG], scale );
  }
  h_res->Scale( 1./m_sumSigmaEff );

  vSampleHin.push_back( h_res );
  
  return h_res;
}

TH2* DiJetAnalysisMC::CombineSamples( std::vector< TH2* >& vSampleHin,
				      const std::string& name ){
  if( !vSampleHin.size() ){ return NULL; }

  TH2* h_res = static_cast< TH2D* >
    ( vSampleHin[0]->Clone( Form("h_%s_%s", name.c_str() , m_allName.c_str() ) ) );
  h_res->Reset();
  
  for( uint iG = 0; iG < vSampleHin.size(); iG++ ){
    double scale = m_vJznEff[iG] * m_vJznSigma[iG] /
      m_vJznSumOverlayWeights[iG];
    h_res->Add( vSampleHin[iG], scale );
  }
  h_res->Scale( 1./m_sumSigmaEff );

  vSampleHin.push_back( h_res );
  
  return h_res;
}

THnSparse* DiJetAnalysisMC::CombineSamples( std::vector< THnSparse* >& vSampleHin,
				      const std::string& name ){
  if( !vSampleHin.size() ){ return NULL; }

  THnSparse* h_res = static_cast< THnSparseD* >
    ( vSampleHin[0]->Clone( Form("h_%s_%s", name.c_str() , m_allName.c_str() ) ) );
  h_res->Reset();
  
  for( uint iG = 0; iG < vSampleHin.size(); iG++ ){
    double scale = m_vJznEff[iG] * m_vJznSigma[iG] /
      m_vJznSumOverlayWeights[iG];
    h_res->Add( vSampleHin[iG], scale );
  }
  h_res->Scale( 1./m_sumSigmaEff );

  vSampleHin.push_back( h_res );
  
  return h_res;
}

void DiJetAnalysisMC::CombineSamples( TH1* h_res,
				      std::vector< TH1* >& vSampleHin,
				      std::vector< TH1* >& vSampleNentIn,
				      const std::string& name ){
  for( int xBin = 1; xBin <= h_res->GetNbinsX(); xBin++ ){

    double valTot      = 0;
    double valErrorTot = 0;
    double denomTot    = 0;

    for( uint iG = 0; iG < vSampleHin.size(); iG++ ){
      if( vSampleHin[iG]->GetEntries() == 0 ){ continue; }
      
      double valueBin    = vSampleHin      [iG]->GetBinContent( xBin );
      double valueBinErr = vSampleHin      [iG]->GetBinError( xBin );
      double nEntriesBin = vSampleNentIn   [iG]->GetBinContent( xBin );
      double weight      = m_vJznWeights[iG] / m_vJznSumOverlayWeights[iG];

      if( nEntriesBin < m_nMinEntriesFit || valueBin == 0 ){ continue; }
      
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
  if( !m_hPowhegWeights ) { return 1; } 
  
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
    yTitleSigma = "#sigma(" + yTitleMean + ")";
  } else if( !type.compare("recoTruthDeta") ){
    yTitleMean  = "#Delta#eta";
    yTitleSigma = "#sigma(" + yTitleMean + ")";
  } else if( !type.compare("recoTruthDphi") ){
    yTitleMean  = "#Delta#phi";
    yTitleSigma = "#sigma(" + yTitleMean + ")";
  } 
}

void DiJetAnalysisMC::GetInfoBoth( std::string& outSuffix,
			   	   std::string& name_a  , std::string& name_b  ,
				   std::string& label_a , std::string& label_b ,
				   std::string& suffix_a, std::string& suffix_b ){
  outSuffix = m_labelOut;
  name_a    = m_dPhiRecoName  + "_" + m_allName;
  name_b    = m_dPhiTruthName + "_" + m_allName;
  label_a   = "Reco";
  label_b   = "Truth";
  suffix_a  = m_labelOut;
  suffix_b  = m_labelOut;
}


//---------------------------------
//           Plotting
//---------------------------------
void DiJetAnalysisMC::LoadHistograms(){
  m_fIn = TFile::Open( m_rootFname.c_str() ); 
  
  for( uint iG = 0; iG < m_nJzn; iG++){
    std::string jzn = m_vJznLabel[iG];
    
    // -------- maps ---------
    m_vHjznEtaPhiMap.push_back
      ( static_cast< TH2D* >
	( m_fIn->Get( Form("h_etaPhiMap_%s", jzn.c_str() ))));
    m_vHjznEtaPhiMap.back()->SetDirectory(0);

    m_vHjznEtaPtMap.push_back
      ( static_cast< TH2D* >
	( m_fIn->Get( Form("h_etaPtMap_%s", jzn.c_str() ))));
    m_vHjznEtaPtMap.back()->SetDirectory(0);

    // -------- spect --------
    m_vHjznEtaSpectReco.push_back
      ( static_cast< TH2D* >
	( m_fIn->Get
	  ( Form("h_%s_%s", m_etaSpectRecoName.c_str(), jzn.c_str() ))));
    m_vHjznEtaSpectReco.back()->SetDirectory(0);

    // The Nent Histogram is only for the overlay samples
    // where we fill with some weight. Otherwise,
    // SpectTruth = SpectTruthNent (b/c fill with w=1)
    m_vHjznEtaSpectTruth.push_back 
      ( static_cast< TH2D* >
	( m_fIn->
	  Get( Form("h_%s_%s", m_etaSpectTruthName.c_str(), jzn.c_str() ))));
    m_vHjznEtaSpectTruth.back()->SetDirectory(0);
    
    m_vHjznEtaSpectTruthNent.push_back
      ( static_cast< TH2D* >
	( m_fIn->Get
	  ( Form("h_%sNent_%s", m_etaSpectTruthName.c_str(), jzn.c_str() ))));
    m_vHjznEtaSpectTruthNent.back()->SetDirectory(0);

    m_vHjznEtaSpectTruthPaired.push_back
      ( static_cast< TH2D* >
	( m_fIn->Get
	  ( Form("h_%sPaired_%s", m_etaSpectTruthName.c_str(), jzn.c_str() ))));
    m_vHjznEtaSpectTruthPaired.back()->SetDirectory(0);

    // --------- recoTruthRpt ---------
    m_vHjznRecoTruthRpt.push_back
      ( static_cast< TH3D* >
	( m_fIn->Get( Form("h_recoTruthRpt_%s", jzn.c_str() ))));
    m_vHjznRecoTruthRpt.back()->SetDirectory(0);

    m_vHjznRecoTruthRptNent.push_back
      ( static_cast< TH2D* >
	( m_fIn->Get( Form("h_recoTruthRptNent_%s", jzn.c_str() ))));
    m_vHjznRecoTruthRptNent.back()->SetDirectory(0);
    
    // --------- recoTruthDeta ---------
    m_vHjznRecoTruthDeta.push_back
      ( static_cast< TH3D* >
	( m_fIn->Get( Form("h_recoTruthDeta_%s", jzn.c_str() ))));
    m_vHjznRecoTruthDeta.back()->SetDirectory(0);

    m_vHjznRecoTruthDetaNent.push_back
      ( static_cast< TH2D* >
	( m_fIn->Get( Form("h_recoTruthDetaNent_%s", jzn.c_str() ))));
    m_vHjznRecoTruthDetaNent.back()->SetDirectory(0);
    
    // --------- recoTruthDphi ---------
    m_vHjznRecoTruthDphi.push_back
      ( static_cast< TH3D* >
	( m_fIn->Get( Form("h_recoTruthDphi_%s", jzn.c_str() )))); 
    m_vHjznRecoTruthDphi.back()->SetDirectory(0);

    m_vHjznRecoTruthDphiNent.push_back
      ( static_cast< TH2D* >
	( m_fIn->Get( Form("h_recoTruthDphiNent_%s", jzn.c_str() ))));
    m_vHjznRecoTruthDphiNent.back()->SetDirectory(0);

    // -------- dPhi- --------
    m_vHjznDphiReco.push_back
      ( static_cast< THnSparse *>
	( m_fIn->Get( Form("h_%s_%s", m_dPhiRecoName.c_str(), jzn.c_str() ))));  

    m_vHjznDphiRecoNent.push_back
      ( static_cast< THnSparse *>
	( m_fIn->Get( Form("h_%sNent_%s", m_dPhiRecoName.c_str(), jzn.c_str() ))));  
    
    m_vHjznDphiTruth.push_back
      ( static_cast< THnSparse *>
	( m_fIn->Get( Form("h_%s_%s", m_dPhiTruthName.c_str(), jzn.c_str() ))));  

    m_vHjznDphiTruthNent.push_back
      ( static_cast< THnSparse *>
	( m_fIn->Get( Form("h_%sNent_%s", m_dPhiTruthName.c_str(), jzn.c_str() ))));  

  }

  m_fIn->Close();
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
  for( uint iG = 0; iG < m_nJzn; iG++ ){
    std::string jzn  = m_vJznLabel[iG];
    TH3* hJznHin     = vJznHin[iG];
   
    for( int xBin = 1; xBin <= nBinsX; xBin++ ){
      double xBinMin = hJznHin->GetXaxis()->GetBinLowEdge( xBin );
      double xBinMax = hJznHin->GetXaxis()->GetBinUpEdge ( xBin );

      std::string xAxisTitle = hJznHin->GetYaxis()->GetTitle();

      // build mean, sigma, project nev
      TH1* hMean = new TH1D
	( Form("h_%s_%s_%s_%s",
	       type.c_str(),
	       jzn.c_str(),
	       sMean.c_str(),
	       anaTool->GetName(xBinMin, xBinMax, "Eta").c_str() ),
	  Form("%s;%s;%s",
	       anaTool->GetEtaLabel( xBinMin, xBinMax ).c_str(),
	       xAxisTitle.c_str(),
	       yTitleMean.c_str() ),
	  nBinsY, yMin, yMax );
      vMeans[iG].push_back( hMean );
      
      TH1* hSigma = new TH1D
	( Form("h_%s_%s_%s_%s",
	       type.c_str(),
	       jzn.c_str(),
	       sSigma.c_str(),
	       anaTool->GetName(xBinMin, xBinMax, "Eta").c_str()),
	  Form("%s;%s;%s",
	       anaTool->GetEtaLabel( xBinMin, xBinMax ).c_str(),
	       xAxisTitle.c_str(),
	       yTitleSigma.c_str() ),
	  nBinsY, yMin, yMax );
      vSigmas[iG].push_back( hSigma );
      
      TH1* hNent = vJznNentIn[iG]->ProjectionY
	( Form("h_%s_%s_N_%s",
	       type.c_str(),
	       jzn.c_str(),
	       anaTool->GetName(xBinMin, xBinMax, "Eta").c_str() )
	  , xBin, xBin );
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
	  ( Form("h_%s_%s_%s_%s",  type.c_str(), jzn.c_str(), 
		 anaTool->GetName(xBinMin, xBinMax, "Eta").c_str(),
		 anaTool->GetName(yBinMin, yBinMax, "Pt").c_str() ),
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
	  
	DrawAtlasRight();
	drawTool->DrawLeftLatex( 0.18, 0.74, Form("%s", jzn.c_str() ) );
    
	drawTool->DrawLeftLatex
	  ( 0.18, 0.88, anaTool->GetEtaLabel( xBinMin, xBinMax ) );
	drawTool->DrawLeftLatex
	  ( 0.18, 0.81, anaTool->GetLabel( yMin, yMax, "#it{p}_{T}^{Truth}") );

	SaveAsROOT( c, hProj->GetName() );

	// if the fit is "bad" do not set any entries on
	// final histograms
	if( hProj->GetEntries() < m_nMinEntriesFit ){ continue; }

	// FIX THIS. Check quality of fit.
	if( fit->GetParError(1) > 0.25 ){ continue; }
	if( fit->GetParError(2) > 0.03 ){ continue; }
	
	hMean->SetBinContent( yBin, fit->GetParameter(1) );
	hMean->SetBinError  ( yBin, fit->GetParError (1) );

	hSigma->SetBinContent ( yBin, fit->GetParameter(2) );
	hSigma->SetBinError   ( yBin, fit->GetParError (2) );
      }
    } // end loop over xBin (eta)
  }// end loop over jzn

  // Now, loop over the sigmas and means, and plot
  // as function of jz sample and eta
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
      ( Form("h_%s_%s_%s",
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
      ( Form("h_%s_%s_%s",
	     type.c_str(),
	     anaTool->GetName(xBinMin, xBinMax, "Eta").c_str(),
	     sSigma.c_str() ),
	Form("%s;%s;%s",
	     anaTool->GetEtaLabel( xBinMin, xBinMax ).c_str(),
	     xAxisTitle.c_str(),
	     yTitleSigma.c_str() ),
	nBinsY, yMin, yMax );
    vSigmasFinal.push_back( hSigmaFinal );

    CombineSamples( hMeanFinal , vEtaMeanTemp , vEtaNentTemp  );
    CombineSamples( hSigmaFinal, vEtaSigmaTemp, vEtaNentTemp  );
  }

  DrawCanvas( vMeansFinal , Form("h_%s", type.c_str() ), sMean , 4 );
  DrawCanvas( vSigmasFinal, Form("h_%s", type.c_str() ), sSigma, 4 );

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

//---------------------------
//        Drawing
//---------------------------
void DiJetAnalysisMC::DrawCanvas( std::vector< TH1* >& vHIN,
				const std::string& type1,
				const std::string& type2,
				int spacing ){
  TCanvas c("c","c",800,600);
  
  TLegend leg(0.64, 0.61, 0.99, 0.82);
  styleTool->SetLegendStyle( &leg, 0.8 );
  leg.SetFillStyle(0);

  int style = 0;

  // for situations where dont want to
  // plot every single bin 
  // plot every n on canvas
  int dX = vHIN.size()/spacing; // plot every n
  for( uint xRange = 2; xRange < vHIN.size(); xRange += dX){
    styleTool->SetHStyle( vHIN[ xRange], style++ );
    leg.AddEntry( vHIN[ xRange ], vHIN[ xRange ]->GetTitle() );
    vHIN[ xRange ]->SetTitle("");
    vHIN[ xRange ]->Draw("epsame");
    SetMinMax( vHIN[ xRange ], type1, type2 );
  }

  leg.Draw();
  
  double y0 = GetLineHeight( type1 );
  
  double xMin = vHIN.front()->GetXaxis()->GetXmin();
  double xMax = vHIN.front()->GetXaxis()->GetXmax();
  
  TLine line( xMin, y0, xMax, y0);
  line.Draw();

  DrawAtlasRight();
  
  SaveAsAll( c, type1, type2 );
}

//===== MinMax and line drawing =====

void DiJetAnalysisMC::
SetMinMax( TH1* h1, const std::string& type1, const std::string& type2 ){
  // JES JER
  if( type1.find("recoTruthRpt") != std::string::npos ){ 
    if( type2.find("mean") != std::string::npos ){ // sigma
      h1->SetMaximum(1.25);
      h1->SetMinimum(0.75);
    } else if( !type2.compare("sigma") ){ // sigma
      h1->SetMaximum(0.34);
      h1->SetMinimum(0.);
    }
  }
  // ANGLES
  else if( type1.find("recoTruthDeta") != std::string::npos ||
	   type1.find("recoTruthDphi") != std::string::npos ) { 
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
  
  if( type.find("recoTruthRpt") != std::string::npos ){ // JES/JER
    y0 = 1;
    y0 = 1;
  } else if( type.find("rEta")          != std::string::npos ||
	     type.find("recoTruthDphi") != std::string::npos ) { // ANGLES
    y0 = 0;
    y0 = 0;
  }

  return y0;
}

void DiJetAnalysisMC::DrawAtlasRight( double x0, double y0, double scale )
{ drawTool->DrawAtlasInternalMCRight( x0, y0, m_mcTypeLabel, scale ); } 


void DiJetAnalysisMC::DrawAtlasRightBoth( double x0, double y0, double scale )
{ DrawAtlasRight(); }
