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
  m_etaSpectRecoName    = m_etaSpectName + "_" + m_recoName;
  m_etaSpectTruthName   = m_etaSpectName + "_" + m_truthName;

  m_dPhiRecoPtTruthName =
    m_dPhiName + "_" + m_recoName + "_" + m_s_pt + "_" + m_truthName;
  m_dPhiTruthPtRecoName =
    m_dPhiName + "_" + m_truthName + "_" + m_s_pt + "_" + m_recoName;
}

DiJetAnalysisMC::~DiJetAnalysisMC(){}

void DiJetAnalysisMC::Initialize(){
  // This is similar to parent class Initialize()
  // but some differences due to MC labels
  std::string system    = m_is_pPb ? m_s_pPb : m_s_pp;

  m_labelOut = m_is_pPb ? m_s_pPb + "_" + m_sMC : m_s_pp + "_" + m_sMC;    
  
  std::string mcMenu;
  if     ( m_mcType == 0 ){ mcMenu = "2015.pp.pythia8powheg"; }
  else if( m_mcType == 1 ){ mcMenu = "2015.pp.pythia8"      ; }
  else if( m_mcType == 2 ){ mcMenu = "2015.pp.herwig"       ; }
  
  std::string labelOut =
    GetConfig()->GetValue( Form("labelOut.%s",mcMenu.c_str()),"");
  m_labelOut += "_" + labelOut;
  m_mcTypeLabel = GetConfig()->GetValue( Form("mcTypeLabel.%s",mcMenu.c_str()),"");
  
  // Check if the directories exist.
  // If they don't, create them
  m_dirOut   = m_sOutput;
  anaTool->CheckWriteDir( m_dirOut.c_str() );
  m_dirOut   += "/" + m_sOutput + "_" + m_labelOut;
  anaTool->CheckWriteDir( m_dirOut.c_str() );
  
  m_rootFname = m_dirOut + "/" + m_myOutName +"_" + m_labelOut + ".root";
  
  std::cout << "fNameIn/Out: " << m_rootFname << std::endl;

  std::string unfoldingMC = GetConfig()->GetValue( "unfoldingMC", "" );
  
  m_fNameUnfoldingMC = Form( "%s/%s_%s_%s_%s/c_%s_%s_%s_%s_%s.root",
			     m_sOutput.c_str(), m_sOutput.c_str(),
			     system.c_str(), m_sMC.c_str(), unfoldingMC.c_str(),
			     m_myOutName.c_str(), system.c_str(), m_sMC.c_str(),
			     unfoldingMC.c_str(), m_unfoldingFileSuffix.c_str() );
  
  // Get Info On jzn slices
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

  std::string cfNameOut = m_dirOut + "/c_" + m_myOutName + "_" + m_labelOut + ".root";
  
  m_fOut = new TFile( cfNameOut.c_str(),"RECREATE");
  
  // MakeEtaPhiPtMap( m_vHjznEtaPhiMap );
  // MakeEtaPhiPtMap( m_vHjznEtaPtMap  );

  // add a slice "all" to collection
  // rest of plots are for combined slices
  m_vJznLabel.push_back( m_allName );
  
  m_hAllEtaSpectReco  = CombineSamples( m_vHjznEtaSpectReco , m_etaSpectRecoName  );
  m_hAllEtaSpectTruth = CombineSamples( m_vHjznEtaSpectTruth, m_etaSpectTruthName );
  MakeSpectra( m_vHjznEtaSpectReco , m_vJznLabel, m_etaSpectRecoName );
  MakeSpectra( m_vHjznEtaSpectTruth, m_vJznLabel, m_etaSpectTruthName  );

  MakeScaleRes( m_vHjznRecoTruthRpt , m_vHjznRecoTruthRptNent , "recoTruthRpt"  );
  MakeScaleRes( m_vHjznRecoTruthDeta, m_vHjznRecoTruthDetaNent, "recoTruthDeta" );
  MakeScaleRes( m_vHjznRecoTruthDphi, m_vHjznRecoTruthDphiNent, "recoTruthDphi" );

  // Make response matrix
  m_hAllDphiRespMat = CombineSamples( m_vHjznDphiRespMat, m_dPhiRespMatName );
  MakeResponseMatrix( m_vHjznDphiRespMat, m_vJznLabel, m_dPhiRespMatName );
  
  m_hAllDphiReco    = CombineSamples( m_vHjznDphiReco   , m_dPhiRecoName    );
  m_hAllDphiTruth   = CombineSamples( m_vHjznDphiTruth  , m_dPhiTruthName   );
  MakeDeltaPhi( m_vHjznDphiReco , m_vJznLabel, m_dPhiRecoName  );
  MakeDeltaPhi( m_vHjznDphiTruth, m_vJznLabel, m_dPhiTruthName );

  m_hAllDphiRecoPtTruth = CombineSamples( m_vHjznDphiRecoPtTruth, m_dPhiRecoPtTruthName );
  m_hAllDphiTruthPtReco = CombineSamples( m_vHjznDphiTruthPtReco, m_dPhiTruthPtRecoName );
  MakeDeltaPhi( m_vHjznDphiRecoPtTruth, m_vJznLabel, m_dPhiRecoPtTruthName );
  MakeDeltaPhi( m_vHjznDphiTruthPtReco, m_vJznLabel, m_dPhiTruthPtRecoName );

  // For MC need to close the file first,
  // to write dphi response, truth histos.
  // Should add this in two steps but unfolding
  // in MC is really for a test, in data there is no
  // need to do this.
  std::cout << "DONE! Closing " << cfNameOut << std::endl;
  m_fOut->Close(); delete m_fOut;
  std::cout << "......Closed  " << cfNameOut << std::endl;
  // Copy the file. Writing to a file and reading from it once it
  // has been written once open works too,
  // but this is "cleaner". To have separate read and write files.
  TFile::Cp( cfNameOut.c_str(), m_fNameUnfoldingMC.c_str() );
  
  // Open file for updating. Will add the unfolded
  // results on top of the previous ones.
  m_fOut = new TFile( cfNameOut.c_str(),"UPDATE");

  std::cout << "----- Unfolding Data ------" << std::endl;
  // make a vector with just the unfolded result.
  // this is to send it to MakeDeltaPhi(..) to have
  // unfolded results plotted separately
  std::vector< THnSparse*  > m_vHDphiUnfolded;
  std::vector< std::string > m_vLabelUnfolded;

  // open the MC file used for unfolding info.
  // passed to unfolding function.
  TFile* fInMC = TFile::Open( m_fNameUnfoldingMC.c_str() );
  // switch dir back to output file for writing.
  m_fOut->cd(); 

  // make unfolded THnSparse with similar naming convention
  // as the other histograms. At this point, don't care about
  // doing this for all triggers. Altohugh, this can be
  // repeated in a loop with m_allName subsitituted for trigger,
  // and subsequently added to the vectors above.  
  THnSparse* m_hAllDphiRecoUnfolded =
    UnfoldDeltaPhi( m_hAllDphiReco, fInMC, m_dPhiRecoUnfoldedName );
  m_vHDphiUnfolded.push_back( m_hAllDphiRecoUnfolded );
  m_vLabelUnfolded.push_back( m_allName );
  
  MakeDeltaPhi( m_vHDphiUnfolded, m_vLabelUnfolded, m_dPhiRecoUnfoldedName );
  
  std::cout << "DONE! Closing " << cfNameOut << std::endl;
  m_fOut->Close(); delete m_fOut;
  std::cout << "......Closed  " << cfNameOut << std::endl;
}

void DiJetAnalysisMC::PlotHistosTogether(){
  MakeDphiTogether();
  // MakeDphiRecoTruth();
}

//---------------------------------
//            Read Data
//---------------------------------

void DiJetAnalysisMC::SetupHistograms(){
  
  for( uint iG = 0; iG < m_nJzn; iG++){
    std::string jzn = m_vJznLabel[iG];
    
    std::cout << "Making - " << jzn << " histograms " << std::endl;
    
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
    THnSparse* hnDphiReco =
      new THnSparseD( Form("h_%s_%s", m_dPhiRecoName.c_str(), jzn.c_str() ), "",
		      m_nDphiDim, &m_nDphiBins[0],
		      &m_dPhiMin[0], &m_dPhiMax[0] );
    hnDphiReco->GetAxis(0)->Set( m_nVarYstarBinsA, &( m_varYstarBinningA[0] ) );
    hnDphiReco->GetAxis(1)->Set( m_nVarYstarBinsB, &( m_varYstarBinningB[0] ) );
    hnDphiReco->GetAxis(2)->Set( m_nVarPtBins    , &( m_varPtBinning[0]     ) );
    hnDphiReco->GetAxis(3)->Set( m_nVarPtBins    , &( m_varPtBinning[0]     ) );
    hnDphiReco->GetAxis(4)->Set( m_nVarDphiBins  , &( m_varDphiBinning[0]   ) );
    m_vHjznDphiReco.push_back( hnDphiReco );
    AddHistogram( hnDphiReco );

    THnSparse* hnDphiTruth =
      new THnSparseD( Form("h_%s_%s", m_dPhiTruthName.c_str(), jzn.c_str() ), "",
		      m_nDphiDim, &m_nDphiBins[0],
		      &m_dPhiMin[0], &m_dPhiMax[0] );
    hnDphiTruth->GetAxis(0)->Set( m_nVarYstarBinsA, &( m_varYstarBinningA[0] ) );
    hnDphiTruth->GetAxis(1)->Set( m_nVarYstarBinsB, &( m_varYstarBinningB[0] ) );
    hnDphiTruth->GetAxis(2)->Set( m_nVarPtBins    , &( m_varPtBinning[0]     ) );
    hnDphiTruth->GetAxis(3)->Set( m_nVarPtBins    , &( m_varPtBinning[0]     ) );
    hnDphiTruth->GetAxis(4)->Set( m_nVarDphiBins  , &( m_varDphiBinning[0]   ) );
    m_vHjznDphiTruth.push_back( hnDphiTruth );
    AddHistogram( hnDphiTruth );

    // --- dPhi truth reco together ----
    THnSparse* hnDphiRecoPtTruth =
      new THnSparseD( Form("h_%s_%s", m_dPhiRecoPtTruthName.c_str(), jzn.c_str() ), "",
		      m_nDphiDim, &m_nDphiBins[0],
		      &m_dPhiMin[0], &m_dPhiMax[0] );
    hnDphiRecoPtTruth->GetAxis(0)->Set( m_nVarYstarBinsA, &( m_varYstarBinningA[0] ) );
    hnDphiRecoPtTruth->GetAxis(1)->Set( m_nVarYstarBinsB, &( m_varYstarBinningB[0] ) );
    hnDphiRecoPtTruth->GetAxis(2)->Set( m_nVarPtBins    , &( m_varPtBinning[0]     ) );
    hnDphiRecoPtTruth->GetAxis(3)->Set( m_nVarPtBins    , &( m_varPtBinning[0]     ) );
    hnDphiRecoPtTruth->GetAxis(4)->Set( m_nVarDphiBins  , &( m_varDphiBinning[0]   ) );
    m_vHjznDphiRecoPtTruth.push_back( hnDphiRecoPtTruth );
    AddHistogram( hnDphiRecoPtTruth );    

    THnSparse* hnDphiTruthPtReco =
      new THnSparseD( Form("h_%s_%s", m_dPhiTruthPtRecoName.c_str(), jzn.c_str() ), "",
		      m_nDphiDim, &m_nDphiBins[0],
		      &m_dPhiMin[0], &m_dPhiMax[0] );
    hnDphiTruthPtReco->GetAxis(0)->Set( m_nVarYstarBinsA, &( m_varYstarBinningA[0] ) );
    hnDphiTruthPtReco->GetAxis(1)->Set( m_nVarYstarBinsB, &( m_varYstarBinningB[0] ) );
    hnDphiTruthPtReco->GetAxis(2)->Set( m_nVarPtBins    , &( m_varPtBinning[0]     ) );
    hnDphiTruthPtReco->GetAxis(3)->Set( m_nVarPtBins    , &( m_varPtBinning[0]     ) );
    hnDphiTruthPtReco->GetAxis(4)->Set( m_nVarDphiBins  , &( m_varDphiBinning[0]   ) );
    m_vHjznDphiTruthPtReco.push_back( hnDphiTruthPtReco );
    AddHistogram( hnDphiTruthPtReco );    
    
    // -------- Dphi Response Matrix --------    
    // copy vectors from 5d dhpi histos
    m_nDphiRespMatBins = m_nDphiBins;
    m_dPhiRespMatMin   = m_dPhiMin;
    m_dPhiRespMatMax   = m_dPhiMax;

    // add one more dimension (6th)
    m_nDphiRespMatBins.push_back( m_nVarDphiBins );
    m_dPhiRespMatMin  .push_back( 0 );
    m_dPhiRespMatMax  .push_back( 1 );
	
    m_nDphiRespMatDim  = m_nDphiRespMatBins.size(); // just adding truth dPhi axis
    
    THnSparse* hnRespMat =
      new THnSparseD( Form("h_%s_%s", m_dPhiRespMatName.c_str(), jzn.c_str() ), "",
		      m_nDphiRespMatDim, &m_nDphiRespMatBins[0],
		      &m_dPhiRespMatMin[0], &m_dPhiRespMatMax[0] );
    hnRespMat->GetAxis(0)->Set( m_nVarYstarBinsA, &( m_varYstarBinningA[0] ) );
    hnRespMat->GetAxis(1)->Set( m_nVarYstarBinsB, &( m_varYstarBinningB[0] ) );
    hnRespMat->GetAxis(2)->Set( m_nVarPtBins    , &( m_varPtBinning[0]     ) );
    hnRespMat->GetAxis(3)->Set( m_nVarPtBins    , &( m_varPtBinning[0]     ) );
    hnRespMat->GetAxis(4)->Set( m_nVarDphiBins  , &( m_varDphiBinning[0]   ) );
    hnRespMat->GetAxis(5)->Set( m_nVarDphiBins  , &( m_varDphiBinning[0]   ) );
    m_vHjznDphiRespMat.push_back( hnRespMat );
    AddHistogram( hnRespMat );
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
    int nEventsTotal, endEvent;

    nEventsTotal = m_tree->GetEntries();
    nEvents      = nEvents > 0 ? nEvents : nEventsTotal;
    startEvent   = startEvent < nEventsTotal ?
				startEvent : nEventsTotal - 1;
    endEvent     = startEvent + nEvents < nEventsTotal ?
					  startEvent + nEvents : nEventsTotal;

    // -------- EVENT LOOP ---------
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
	( m_vHjznDphiReco [iG], vRR_paired_jets, GetJetWeight );
      AnalyzeDeltaPhi
	( m_vHjznDphiTruth[iG], vT_jets, GetJetWeight );

      // Do Dphi analysis on reco/truth mix
      AnalyzeDeltaPhiTruthReco
	( m_vHjznDphiRecoPtTruth[iG], m_vHjznDphiTruthPtReco[iG],
	  vRR_paired_jets, vRT_paired_jets, GetJetWeight );
      
      std::vector< TLorentzVector > vTR_paired_jets;
      std::vector< TLorentzVector > vTT_paired_jets;
      PairJets( vT_jets, vR_jets, vTT_paired_jets, vTR_paired_jets );

      AnalyzeResponseMatrix
	( m_vHjznDphiRespMat[iG], vTR_paired_jets, vTT_paired_jets, GetJetWeight );
      
      // loop over truth jets
      // denominator for efficiency because not all
      // reco jets are reconstructed for a truth jet
      for( auto& tJet : vT_jets ){
	double jetEta    = tJet.Eta();
	double jetEtaAdj = AdjustEtaForPP( jetEta );
	double jetPhi    = tJet.Phi();
	double jetPt     = tJet.Pt()/1000.;

	double weight    = GetJetWeight( jetEta, jetPhi, jetPt );

	m_vHjznEtaSpectTruth    [iG]->
	  Fill( jetEtaAdj, jetPt, weight);
      } // end loop over truth jets

      // do JER/JES, angular scales and resolution.
      AnalyzeScaleResolution( vTR_paired_jets, vTT_paired_jets, iG );
    } // -------- END EVENT LOOP ---------
   
    std::cout << "DONE WITH " << m_vJznLabel[iG] << std::endl;

    m_fIn->Close(); delete m_fIn;
  } // end loop over a JZ sample
}

//---------------------------
//          Analysis
//---------------------------

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
    double jetEtaRecoAdj  = AdjustEtaForPP( jetEtaReco  );
    double jetEtaTruthAdj = AdjustEtaForPP( jetEtaTruth );	  

    double weightReco  = GetJetWeight( jetEtaReco , jetPhiReco , jetPtReco  );
    double weightTruth = GetJetWeight( jetEtaTruth, jetPhiTruth, jetPtTruth );
	
    m_vHjznEtaPhiMap[iG]->Fill( jetEtaReco, jetPhiReco, weightReco );
    m_vHjznEtaPtMap [iG]->Fill( jetEtaReco, jetPtReco , weightReco );

    m_vHjznEtaSpectReco       [iG]->
      Fill( jetEtaRecoAdj,  jetPtReco,  weightReco );

    m_vHjznEtaSpectTruthPaired[iG]->
      Fill( jetEtaTruthAdj,  jetPtTruth,  weightTruth );

    m_vHjznRecoTruthRpt       [iG]->
      Fill( jetEtaTruthAdj, jetPtTruth, jetPtReco/jetPtTruth, weightTruth );
    m_vHjznRecoTruthRptNent   [iG]->
      Fill( jetEtaTruthAdj, jetPtTruth );

    m_vHjznRecoTruthDeta      [iG]->
      Fill( jetEtaTruthAdj, jetPtTruth, jetEtaReco - jetEtaTruth, weightTruth );
    m_vHjznRecoTruthDetaNent  [iG]->
      Fill( jetEtaTruthAdj, jetPtTruth );

    m_vHjznRecoTruthDphi      [iG]->
      Fill( jetEtaTruthAdj, jetPtTruth, jetPhiReco - jetPhiTruth, weightTruth );
    m_vHjznRecoTruthDphiNent  [iG]->
      Fill( jetEtaTruthAdj, jetPtTruth  );
  } // end loop over pairs
}

void  DiJetAnalysisMC::AnalyzeResponseMatrix( THnSparse* hn,
					      const std::vector< TLorentzVector >& vR_jets,
					      const std::vector< TLorentzVector >& vT_jets,
					      WeightFcn weightFcn ){

  const TLorentzVector* recoJet1  = NULL; const TLorentzVector* recoJet2  = NULL;
  const TLorentzVector* truthJet1 = NULL; const TLorentzVector* truthJet2 = NULL;

  if( !GetDiJets( vR_jets, recoJet1 , recoJet2  ) )
    { return; }
  if( !GetDiJets( vT_jets, truthJet1, truthJet2 ) )
    { return; }

  double truthJet1_pt    = truthJet1->Pt()/1000.;
  double truthJet1_eta   = truthJet1->Eta();
  double truthJet1_phi   = truthJet1->Phi();
  double truthJet1_ystar = GetYstar( *truthJet1 );

  double truthJet2_pt    = truthJet2->Pt()/1000.;
  double truthJet2_ystar = GetYstar( *truthJet2 );
  
  double recoDeltaPhi    = anaTool->DeltaPhi( *recoJet2 , *recoJet1  );
  double truthDeltaPhi   = anaTool->DeltaPhi( *truthJet2, *truthJet1 );
  
  std::vector< double > x;
  x.resize( hn->GetNdimensions() );
    
  // wont change unless we have a weightFcn
  // and then it varys depending on eta, phi, pt.
  double weight = weightFcn ?
    weightFcn( truthJet1_eta, truthJet1_phi, truthJet1_pt ) : 1;   

  x[0] = truthJet1_ystar;  
  x[1] = truthJet2_ystar;
  x[2] = truthJet1_pt ;
  x[3] = truthJet2_pt ;
  x[4] = truthDeltaPhi;
  x[5] = recoDeltaPhi;
  hn    ->Fill( &x[0], weight );

  // for pp, fill twice. once for each side since
  // it is symmetric in pp. For pPb, continue
  if( m_is_pPb ){ return; }
  
  x[0] = -truthJet1_ystar;  
  x[1] = -truthJet2_ystar;
  hn    ->Fill( &x[0], weight );
}

std::pair<double,double> DiJetAnalysisMC::AnalyzeDeltaPhiTruthReco
( THnSparse* hnA, THnSparse* hnB,
  const std::vector< TLorentzVector >& v_jets_A,
  const std::vector< TLorentzVector >& v_jets_B,
  WeightFcn weightFcn ){

  const TLorentzVector* jet1A = NULL; const TLorentzVector* jet2A = NULL;
  const TLorentzVector* jet1B = NULL; const TLorentzVector* jet2B = NULL;

  if( !GetDiJets( v_jets_A, jet1A, jet2A ) )
    { return std::make_pair( -1, -1 ); }

  if( !GetDiJets( v_jets_B, jet1B, jet2B ) )
    { return std::make_pair( -1, -1 ); }

  double jet1A_pt    = jet1A->Pt()/1000.;
  double jet1A_eta   = jet1A->Eta();
  double jet1A_phi   = jet1A->Phi();
  double jet1A_ystar = GetYstar( *jet1A );

  double jet2A_pt    = jet2A->Pt()/1000.;
  double jet2A_ystar = GetYstar( *jet2A );
  
  double deltaPhiA = anaTool->DeltaPhi( *jet2A, *jet1A );

  double jet1B_pt    = jet1B->Pt()/1000.;
  double jet1B_eta   = jet1B->Eta();
  double jet1B_phi   = jet1B->Phi();
  double jet1B_ystar = GetYstar( *jet1B );

  double jet2B_pt    = jet2B->Pt()/1000.;
  double jet2B_ystar = GetYstar( *jet2B );
  
  double deltaPhiB = anaTool->DeltaPhi( *jet2B, *jet1B );

  std::vector< double > xA;
  std::vector< double > xB;
  xA.resize( hnA->GetNdimensions() );
  xB.resize( hnB->GetNdimensions() );
  
  // wont change unless we have a weightFcn
  // and then it varys depending on eta, phi, pt.
  double weightA = weightFcn ? weightFcn( jet1B_eta, jet1B_phi, jet1B_pt ) : 1;   
  double weightB = weightFcn ? weightFcn( jet1A_eta, jet1A_phi, jet1A_pt ) : 1;

  xA[0] = jet1A_ystar;  
  xA[1] = jet2A_ystar;
  xA[2] = jet1B_pt ;
  xA[3] = jet2B_pt ;
  xA[4] = deltaPhiA;
  hnA->Fill( &xA[0], weightA );

  xB[0] = jet1A_ystar;  
  xB[1] = jet2A_ystar;
  xB[2] = jet1A_pt ;
  xB[3] = jet2A_pt ;
  xB[4] = deltaPhiB;
  hnB->Fill( &xB[0], weightB );
  
  // for pp, fill twice. once for each side since
  // it is symmetric in pp. For pPb, continue
  if( m_is_pPb ){ return std::make_pair( deltaPhiA, deltaPhiB ); }
  
  xA[0] = -jet1A_ystar;  
  xA[1] = -jet2A_ystar;
  hnA->Fill( &xA[0], weightA );

  xB[0] = -jet1A_ystar;  
  xB[1] = -jet2A_ystar;
  hnB->Fill( &xB[0], weightB );
 
  return std::make_pair( deltaPhiA, deltaPhiB );
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
    if( !pairedBJet ){ continue; }

    // Require jets to be within some radius
    // If they aren't, dont add a pair.
    if( deltaRmin > m_dRmax ){ continue; }
    
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

void DiJetAnalysisMC::GetInfoBothRecoTruth
( std::string& name_a  , std::string& name_b  ,
  std::string& label_a , std::string& label_b ){
  name_a  = m_dPhiRecoPtTruthName + "_" + m_allName;
  name_b  = m_dPhiTruthPtRecoName + "_" + m_allName;
  label_a = "|#Delta#phi|_{Reco} #it{p}_{T}^{Truth}";
  label_b = "|#Delta#phi|_{Truth} #it{p}_{T}^{Reco}";
}

void DiJetAnalysisMC::GetInfoUnfolding( std::string& measuredName,
					std::string& measuredLabel ){
  measuredName  = m_dPhiRecoName;
  measuredLabel = "Reco";;
}

//---------------------------------
//     Get Quantities / Plot 
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
    
    m_vHjznDphiTruth.push_back
      ( static_cast< THnSparse *>
	( m_fIn->Get( Form("h_%s_%s", m_dPhiTruthName.c_str(), jzn.c_str() ))));  

    // --- dPhi truth reco together ----
    m_vHjznDphiRecoPtTruth.push_back
      ( static_cast< THnSparse *>
	( m_fIn->Get( Form("h_%s_%s", m_dPhiRecoPtTruthName.c_str(), jzn.c_str() ))));  
    
    m_vHjznDphiTruthPtReco.push_back
      ( static_cast< THnSparse *>
	( m_fIn->Get( Form("h_%s_%s", m_dPhiTruthPtRecoName.c_str(), jzn.c_str() ))));  
    
    // -------- Dphi Response Matrix --------    
    m_vHjznDphiRespMat.push_back
      ( static_cast< THnSparse *>
	( m_fIn->Get( Form("h_%s_%s", m_dPhiRespMatName.c_str(), jzn.c_str() ))));  
  }

  m_fIn->Close(); delete m_fIn;
}

void DiJetAnalysisMC::MakeScaleRes( std::vector< TH3* >& vJznHin,
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

void DiJetAnalysisMC::MakeDphiRecoTruth(){
  std::string outSuffix;
  std::string name_a  , name_b  , name_c , name_d ;
  std::string label_a , label_b , label_c, label_d;
  std::string suffix_a, suffix_b;

  GetInfoBoth( outSuffix, name_a, name_b, label_a, label_b, suffix_a, suffix_b );
  GetInfoBothRecoTruth( name_c, name_d, label_c, label_d );
  
  // Check if the directories exist.
  // If they don't, create them
  std::string outDir = m_sOutput;
  anaTool->CheckWriteDir( outDir.c_str() );
  outDir += "/" + m_allName;
  anaTool->CheckWriteDir( outDir.c_str() );
  outDir += "/" + outSuffix;
  anaTool->CheckWriteDir( outDir.c_str() );
  
  TAxis* axis0 = m_dPP->GetTAxis( 0 );
  TAxis* axis1 = m_dPP->GetTAxis( 1 );
  TAxis* axis2 = m_dPP->GetTAxis( 2 );
  TAxis* axis3 = m_dPP->GetTAxis( 3 );
  
  int nAxis0Bins = axis0->GetNbins();
  int nAxis1Bins = axis1->GetNbins();
  int nAxis2Bins = axis2->GetNbins();
  int nAxis3Bins = axis3->GetNbins();

  TFile* fIn = TFile::Open
    ( Form("%s/%s_%s/c_%s_%s.root",
	   m_sOutput.c_str(), m_sOutput.c_str(),suffix_a.c_str(),
	   m_myOutName.c_str(),suffix_a.c_str() ) );
  TFile* fOut  = new TFile
    ( Form("%s/%s/%s/c_%s_%s.root",
	   m_sOutput.c_str(), m_allName.c_str(), outSuffix.c_str(),
	   m_myOutName.c_str(), outSuffix.c_str() ), "update");

  for( int axis0Bin = 1; axis0Bin <= nAxis0Bins; axis0Bin++ ){
    double axis0Low, axis0Up;
    anaTool->GetBinRange
      ( axis0, axis0Bin, axis0Bin, axis0Low, axis0Up );
    
    for( int axis1Bin = 1; axis1Bin <= nAxis1Bins; axis1Bin++ ){
      // check we are in correct ystar and pt bins
      if( !m_dPP->CorrectPhaseSpace
	  ( std::vector<int>{ axis0Bin, axis1Bin, 0, 0 } ) )
	{ continue; }

      double axis1Low, axis1Up;
      anaTool->GetBinRange
	( axis1, axis1Bin, axis1Bin, axis1Low, axis1Up );

      // get widths canvases
      std::string hTagCW =
	Form ("%s_%s",
	      	anaTool->GetName( axis0Low, axis0Up, m_dPP->GetAxisName(0) ).c_str(),
		anaTool->GetName( axis1Low, axis1Up, m_dPP->GetAxisName(1) ).c_str() );
    
      for( int axis2Bin = 1; axis2Bin <= nAxis2Bins; axis2Bin++ ){
	double axis2Low , axis2Up;
	anaTool->GetBinRange
	  ( axis2, axis2Bin, axis2Bin, axis2Low, axis2Up );

	for( int axis3Bin = 1; axis3Bin <= nAxis3Bins; axis3Bin++ ){
	  // check we are in correct ystar and pt bins
	  if( !m_dPP->CorrectPhaseSpace
	      ( std::vector<int>{ axis0Bin, axis1Bin, axis2Bin, axis3Bin } ) )
	    { continue; }

	  double axis3Low , axis3Up;
	  anaTool->GetBinRange
	    ( axis3, axis3Bin, axis3Bin, axis3Low, axis3Up );
	    
	  std::string hTag =
	    Form("%s_%s_%s_%s",
		 anaTool->GetName( axis0Low, axis0Up, m_dPP->GetAxisName(0) ).c_str(),
		 anaTool->GetName( axis1Low, axis1Up, m_dPP->GetAxisName(1) ).c_str(),
		 anaTool->GetName( axis2Low, axis2Up, m_dPP->GetAxisName(2) ).c_str(),
		 anaTool->GetName( axis3Low, axis3Up, m_dPP->GetAxisName(3) ).c_str() );
	  
	  std::string hName_a = Form("h_%s_%s", name_a.c_str(), hTag.c_str() );
	  std::string hName_b = Form("h_%s_%s", name_b.c_str(), hTag.c_str() );
	  std::string hName_c = Form("h_%s_%s", name_c.c_str(), hTag.c_str() );
	  std::string hName_d = Form("h_%s_%s", name_d.c_str(), hTag.c_str() );

	  std::cout << hName_c << std::endl;
	  
	  TH1* h_a = static_cast<TH1D*>( fIn->Get( hName_a.c_str() ) );
	  TH1* h_b = static_cast<TH1D*>( fIn->Get( hName_b.c_str() ) );
	  TH1* h_c = static_cast<TH1D*>( fIn->Get( hName_c.c_str() ) );
	  TH1* h_d = static_cast<TH1D*>( fIn->Get( hName_d.c_str() ) );
	  styleTool->SetHStyle( h_a, 0 );
	  styleTool->SetHStyle( h_b, 5 );
	  styleTool->SetHStyle( h_c, 1 );
	  styleTool->SetHStyle( h_d, 6 );
	  
	  TF1* f_a = static_cast<TF1*>( fIn->Get( Form("f_%s", hName_a.c_str())));
	  TF1* f_b = static_cast<TF1*>( fIn->Get( Form("f_%s", hName_b.c_str())));
	  TF1* f_c = static_cast<TF1*>( fIn->Get( Form("f_%s", hName_c.c_str())));
	  TF1* f_d = static_cast<TF1*>( fIn->Get( Form("f_%s", hName_d.c_str())));
	  styleTool->SetHStyle( f_a, 0 );
	  styleTool->SetHStyle( f_b, 5 );
	  styleTool->SetHStyle( f_a, 1 );
	  styleTool->SetHStyle( f_b, 6 );
	  f_a->SetLineColor( h_a->GetLineColor() );
	  f_b->SetLineColor( h_b->GetLineColor() );
	  f_c->SetLineColor( h_c->GetLineColor() );
	  f_d->SetLineColor( h_d->GetLineColor() );

	  double chi2NDF_a = f_a->GetChisquare()/f_a->GetNDF();
	  double chi2NDF_b = f_b->GetChisquare()/f_b->GetNDF();
	  double chi2NDF_c = f_a->GetChisquare()/f_d->GetNDF();
	  double chi2NDF_d = f_b->GetChisquare()/f_c->GetNDF();
	  
	  TCanvas c("c","c", 800, 600 );
	    
	  TLegend leg( 0.27, 0.38, 0.38, 0.57 );
	  styleTool->SetLegendStyle( &leg , 0.85 );

	  bool save = false;
	  
	  if( h_a->GetEntries() ){
	    h_a->SetMinimum(0);
	    h_a->GetYaxis()->SetNdivisions(505);
	    leg.AddEntry( h_a, Form("%s #Chi^{2}/NDF=%4.2f", label_a.c_str(), chi2NDF_a));
	    h_a->Draw("epsame");
	    f_a->Draw("same");
	    save = true;
	  }

	  if( h_b->GetEntries() ){
	    h_b->SetMinimum(0);
	    h_b->GetYaxis()->SetNdivisions(505);
	    leg.AddEntry( h_b, Form("%s #Chi^{2}/NDF=%4.2f", label_b.c_str(), chi2NDF_b));
	    h_b->Draw("epsame");
	    f_b->Draw("same");
	    save = true;
	  }

	  if( h_c->GetEntries() ){
	    h_c->SetMinimum(0);
	    h_c->GetYaxis()->SetNdivisions(505);
	    leg.AddEntry( h_c, Form("%s #Chi^{2}/NDF=%4.2f", label_c.c_str(), chi2NDF_c));
	    h_c->Draw("epsame");
	    f_c->Draw("same");
	    save = true;
	  }

	  if( h_d->GetEntries() ){
	    h_d->SetMinimum(0);
	    h_d->GetYaxis()->SetNdivisions(505);
	    leg.AddEntry( h_d, Form("%s #Chi^{2}/NDF=%4.2f", label_d.c_str(), chi2NDF_d));
	    h_d->Draw("epsame");
	    f_d->Draw("same");
	    save = true;
	  }

	  leg.Draw("same");

	  if( h_a->GetMaximum() > h_b->GetMaximum() ){
	    h_a->SetMaximum( h_a->GetMaximum() * 1.1 );
	    h_b->SetMaximum( h_a->GetMaximum() * 1.1 );
	  } else {
	    h_a->SetMaximum( h_b->GetMaximum() * 1.1 );
	    h_b->SetMaximum( h_b->GetMaximum() * 1.1 );
	  }

	  DrawTopLeftLabels
	    ( m_dPP, axis0Low, axis0Up, axis1Low, axis1Up,
	      axis2Low, axis2Up, axis3Low, axis3Up, 0.8 );

	  DrawAtlasRightBoth();

	  if( save ) c.SaveAs( Form("%s/%s/%s/h_%s_%s_%s.png",
				    m_sOutput.c_str(), m_allName.c_str(),
				    outSuffix.c_str(), m_dPhiName.c_str(),
				    hTag.c_str(),outSuffix.c_str() ) );

	  SaveAsROOT( c, Form("h_%s_%s", m_dPhiName.c_str(), hTag.c_str() ) );
	  
	  delete  f_a; delete  f_b; delete f_c; delete f_d;
	  delete  h_a; delete  h_b; delete h_c; delete h_d;
	} // end loop over axis3
      } // end loop over axis2
    } // end loop over ystar2
  } // end loop over ystar2
  
  std::cout << "DONE! Closing " << fOut->GetName() << std::endl;
  fOut->Close();
  std::cout << "......Closed  " << fOut->GetName() << std::endl;
}

void DiJetAnalysisMC::MakeResponseMatrix( std::vector< THnSparse* >& vhn,
					  const std::vector< std::string >& vLabel,
					  const std::string& name ){  
  std::vector< TH2* > vDphiRespMat;
  
  // ---- loop over group  ----
  // ---- (jzn or trigger) ----
  for( uint iG = 0; iG < vhn.size(); iG++ ){      
    THnSparse* hn = vhn[iG];
  
    std::string label = vLabel[iG];

    // in data only draw for all
    if( label.compare( m_allName ) ){ continue; } 

    TAxis* axis0 = hn->GetAxis( m_dPP->GetAxisI(0) );
    TAxis* axis1 = hn->GetAxis( m_dPP->GetAxisI(1) );
    TAxis* axis2 = hn->GetAxis( m_dPP->GetAxisI(2) );
    TAxis* axis3 = hn->GetAxis( m_dPP->GetAxisI(3) );
    
    int nAxis0Bins = axis0->GetNbins();
    int nAxis1Bins = axis1->GetNbins();
    int nAxis2Bins = axis2->GetNbins();
    int nAxis3Bins = axis3->GetNbins();
    
    // ---- loop over ystars ----
    for( int axis0Bin = 1; axis0Bin <= nAxis0Bins; axis0Bin++ ){
      // set ranges
      axis0->SetRange( axis0Bin, axis0Bin );
      double axis0Low, axis0Up;
      anaTool->GetBinRange
	( axis0, axis0Bin, axis0Bin, axis0Low, axis0Up );
      for( int axis1Bin = 1; axis1Bin <= nAxis1Bins; axis1Bin++ ){
	// set ranges
	axis1->SetRange( axis1Bin, axis1Bin );
	double axis1Low, axis1Up;
	anaTool->GetBinRange
	  ( axis1, axis1Bin, axis1Bin, axis1Low, axis1Up );
	// ---- loop over axis2 ----
	for( int axis2Bin = 1; axis2Bin <= nAxis2Bins; axis2Bin++ ){
	  // set ranges
	  axis2->SetRange( axis2Bin, axis2Bin );
	  double axis2Low , axis2Up;
	  anaTool->GetBinRange
	    ( axis2, axis2Bin, axis2Bin, axis2Low, axis2Up );
	  // ---- loop over axis3 ----
	  for( int axis3Bin = 1; axis3Bin <= nAxis3Bins; axis3Bin++ ){
	    // check we are in correct ystar and pt bins
	    if( !m_dPP->CorrectPhaseSpace
		( std::vector<int>{ axis0Bin, axis1Bin, axis2Bin, axis3Bin } ) )
	      { continue; }

	    axis3->SetRange( axis3Bin, axis3Bin );

	    double axis3Low , axis3Up;
	    anaTool->GetBinRange
	      ( axis3, axis3Bin, axis3Bin, axis3Low, axis3Up );

	    std::string hTag =
	      Form( "%s_%s_%s_%s",
		    anaTool->GetName( axis0Low, axis0Up, m_dPP->GetAxisName(0) ).c_str(),
		    anaTool->GetName( axis1Low, axis1Up, m_dPP->GetAxisName(1) ).c_str(),
		    anaTool->GetName( axis2Low, axis2Up, m_dPP->GetAxisName(2) ).c_str(),
		    anaTool->GetName( axis3Low, axis3Up, m_dPP->GetAxisName(3) ).c_str() ); 

	    // Take projection onto the dPhi axis
	    TH2* hDphiRespMat = hn->Projection( 4, 5 );
	    hDphiRespMat->SetName
	      ( Form( "h_%s_%s_%s", name.c_str(), label.c_str(), hTag.c_str() ) );
	    styleTool->SetHStyle( hDphiRespMat, 0 );
	    vDphiRespMat.push_back( hDphiRespMat );
	    
	    TCanvas c( "c", hDphiRespMat->GetName(), 800, 600 );
	    c.SetLogz();
	    
	    hDphiRespMat->Draw("col");
	    hDphiRespMat->SetTitle("");
	    
	    DrawTopLeftLabels
	      ( m_dPP, axis0Low, axis0Up, axis1Low, axis1Up,
		axis2Low, axis2Up, axis3Low, axis3Up, 0.8 );
	    
	    DrawAtlasRight();
	    
	    SaveAsROOT( c, hDphiRespMat->GetName() );
	    hDphiRespMat->Write();
	  } // end loop over axis3
	} // end loop over axis2
      } // end loop over axis1     
    } // end loop over axis0
  } // end loop over iG
  for( auto rm : vDphiRespMat ){ delete rm; }
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
