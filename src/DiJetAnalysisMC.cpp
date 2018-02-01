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
#include "UncertaintyProvider.h"

DiJetAnalysisMC::DiJetAnalysisMC()
  : DiJetAnalysisMC( false, 0, 0 ) {}

DiJetAnalysisMC::DiJetAnalysisMC( bool is_pPb )
  : DiJetAnalysisMC( is_pPb, 0 , 0 ) {}

DiJetAnalysisMC::DiJetAnalysisMC( bool is_pPb, int mcType )
  : DiJetAnalysisMC( is_pPb, mcType, 0 ){}

DiJetAnalysisMC::DiJetAnalysisMC( bool is_pPb, int mcType, int uncertComp )
  : DiJetAnalysis( is_pPb, false, mcType, uncertComp ){
  
  //========== Set Histogram Binning =============
  // ------ truth binning --------
  m_ptTruthWidth  = 5;
  m_ptTruthMin    = 20;
  m_ptTruthMax    = 120;
  m_nPtTruthBins  =
    (m_ptTruthMax - m_ptTruthMin) / m_ptTruthWidth;

  // ---- JES/PRes/Etc ----- 
  // --- variable eta/ystar binning ---
  m_nRPtRecoTruthBins   = 100;
  m_rPtRecoTruthMin     = 0;  m_rPtRecoTruthMax = 2;

  m_nDAngleRecoTruthBins = 100;
  m_dAngleRecoTruthMin   = -0.5; m_dAngleRecoTruthMax = 0.5;
  
  // --- Dphi Resp Matrix binning ---  

  boost::assign::push_back( m_vNdPhiRespMatBins )
    ( m_nVarYstarBins     )( m_nVarYstarBins     )
    ( m_nVarPtBinsRespMat )( m_nVarPtBinsRespMat )
    ( m_nVarPtBinsRespMat )( m_nVarPtBinsRespMat )
    ( m_nVarDphiBins      )( m_nVarDphiBins      );

  boost::assign::push_back( m_vNdPhiRespMatRebBins )
    ( m_nVarYstarBins     )( m_nVarYstarBins     )
    ( m_nVarPtBinsRespMat )( m_nVarPtBinsRespMat )
    ( m_nVarPtBinsRespMat )( m_nVarPtBinsRespMat )
    ( m_nVarDphiRebinnedBins  )( m_nVarDphiRebinnedBins  );
  
  boost::assign::push_back( m_vDphiRespMatMin  )
    ( 0 )( 0 )( 0 )( 0 )( 0 )( 0 )( 0 )( 0 );

  boost::assign::push_back( m_vDphiRespMatMax  )
    ( 1 )( 1 )( 1 )( 1 )( 1 )( 1 )( 1 )( 1 );

  m_nDphiRespMatDim = m_vNdPhiRespMatBins.size();
  
  //==================== Cuts ====================    
  m_dRmax = 0.2;

  //=============== Histo Names ==================    
  m_dPhiRecoPairedTruthName =
    m_dPhiName + "_" + m_recoName  + "_" + m_pairedName + "_" + m_truthName;

  m_dPhiRespMatName        = m_dPhiName + "_" + m_respMatName;
  m_dPhiRespMatRebName     = m_dPhiName + "_" + m_respMatName + "_" + m_sReb;

  m_dPhiRecoUnfoldedName   = m_dPhiRecoName + "_" + m_unfoldedName;
  
  m_ptRespMatName          = m_s_pt     + "_" + m_respMatName;

  //============== tool ===============
  m_uncertaintyProvider = new UncertaintyProvider( m_uncertComp, m_is_pPb );

}

DiJetAnalysisMC::~DiJetAnalysisMC(){
  delete m_uncertaintyProvider;
}

void DiJetAnalysisMC::Initialize()
{
  // Initalize things common to everything
  DiJetAnalysis::Initialize();
  
  std::string mcMenu = GetMCMenu();

  // Get Info On jzn slices
  std::vector< double > vJznUsedD = anaTool->vectoriseD
    ( GetConfig()->GetValue( Form("jznUsed.%s", mcMenu.c_str()),"" )," ");
  // double -> int 
  for( auto d : vJznUsedD ){ m_vJznUsed.push_back(d); }
  m_vJznLabels  =  anaTool->vectorise
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
    TFile* fIn  = TFile::Open( m_vJznFnameIn[iG].c_str() );
    TTree* tree = static_cast< TTree* >( fIn->Get( "tree" ) );

    int nEventsTotal = tree->GetEntries();
    double     sigma = m_vJznSigma[iG];      
    double       eff = m_vJznEff  [iG];
  
    m_vJznNev.push_back( nEventsTotal );
    m_vJznWeights.push_back
      ( ( 1./nEventsTotal * m_sumSigmaEff ) * ( sigma * eff ) );
  
    std::cout << iG << "   weight = "
	      << m_vJznWeights.back() << std::endl;
    
    fIn->Close();
  }
}

void DiJetAnalysisMC::AdditionalSuffix( std::string& label ){
  label += "_" + m_mcTypeName;
}

std::string DiJetAnalysisMC::GetMCMenu(){
  return m_is_pPb ?
    m_s_pPb + "." + m_mcTypeName :
    m_s_pp  + "." + m_mcTypeName;  
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

void DiJetAnalysisMC::ProcessPerformance(){

  LoadHistograms();

  TFile* fOut = new TFile( m_fNamePerf.c_str(),"RECREATE");
  
  // add a slice "all" to collection
  // rest of plots include combined slices
  m_vJznLabels.push_back( m_allName );

  m_hAllEtaPtMap = CombineSamples( m_vHjznEtaPtMap, "etaPtMap" );
  MakeEtaPhiPtMap( m_vHjznEtaPtMap , m_vJznLabels, "etaPtMap" );
  
  // only do this for the default sample
  if( !m_uncertComp ){
    // need to cd into the original fout, because scale and res files are created 
    MakeScaleRes( m_vHjznRecoTruthRpt , m_vHjznRecoTruthRptNent , "recoTruthRpt"  );
    fOut->cd();
    MakeScaleRes( m_vHjznRecoTruthDeta, m_vHjznRecoTruthDetaNent, "recoTruthDeta" );
    fOut->cd();
    MakeScaleRes( m_vHjznRecoTruthDphi, m_vHjznRecoTruthDphiNent, "recoTruthDphi" );
    fOut->cd();
  }
  
  m_hAllEtaSpectReco    = CombineSamples( m_vHjznEtaSpectReco   , m_etaSpectRecoName    );
  m_hAllEtaSpectTruth   = CombineSamples( m_vHjznEtaSpectTruth  , m_etaSpectTruthName   );
  m_hAllYstarSpectReco  = CombineSamples( m_vHjznYstarSpectReco , m_ystarSpectRecoName  );
  m_hAllYstarSpectTruth = CombineSamples( m_vHjznYstarSpectTruth, m_ystarSpectTruthName );
  MakeSpectra( m_vHjznEtaSpectReco   , m_vJznLabels, m_etaSpectRecoName    );
  MakeSpectra( m_vHjznEtaSpectTruth  , m_vJznLabels, m_etaSpectTruthName   );
  MakeSpectra( m_vHjznYstarSpectReco , m_vJznLabels, m_ystarSpectRecoName  );
  MakeSpectra( m_vHjznYstarSpectTruth, m_vJznLabels, m_ystarSpectTruthName );

  // make ystar spectra response matrix
  m_hAllYstarSpectRespMat = CombineSamples( m_vHjznYstarSpectRespMat, m_ystarSpectRespMatName );

  MakeSpectCFactorsRespMat( m_vHjznYstarSpectReco,  m_vHjznYstarSpectTruth, m_vHjznYstarSpectRespMat,
			    m_vJznLabels, m_ystarSpectCfactorsName, m_ystarSpectRespMatName );
     
  std::cout << "DONE! Closing " << fOut->GetName() << std::endl;
  fOut->Close();
  delete fOut;
  std::cout << "......Closed  " << std::endl;
}

void DiJetAnalysisMC::UnfoldPerformance(){

  // for now
  // Copy File with original dPhi, spectra, etc,
  // into file where histos with corrections are
  // going to be appended. 
  std::cout << "Copy " << m_fNamePerf << " -> " << m_fNamePerfUF << std::endl;
  TFile::Cp( m_fNamePerf.c_str(), m_fNamePerfUF.c_str() );

  // Open two for reading one for updating.
  // open the MC file used for unfolding info.
  // open the data file used for measured info.
  // passed to unfolding function.
  TFile* fInMC   = TFile::Open( m_fNamePerf.c_str() );
  TFile* fInData = TFile::Open( m_fNamePerf.c_str() );
  TFile* fOut    = new TFile( m_fNamePerfUF.c_str(),"UPDATE");

  std::cout << "----- Unfolding MC ------" << std::endl;
  // make a vector with just the unfolded result.
  // this is to send it to MakeDeltaPhi(..) to have
  // unfolded results plotted separately
  std::vector< TH2* >        m_vHspectUnfolded;
  std::vector< std::string > m_vLabelsUnfolded;

  // make unfolded THnSparse with similar naming convention
  // as the other histograms. At this point, don't care about
  // doing this for all triggers. Altohugh, this can be
  // repeated in a loop with m_allName subsitituted for trigger,
  // and subsequently added to the vectors above.  
  TH2* m_hAllspectRecoUnfolded =
    UnfoldSpectra( fInData, fInMC, m_ystarSpectRecoUnfoldedName );
  m_vHspectUnfolded.push_back( m_hAllspectRecoUnfolded );
  m_vLabelsUnfolded.push_back( m_allName );

  // unfold on MC, just used for testing purpose.
  // so there is no comb subt or normalization or scaling
  MakeSpectra( m_vHspectUnfolded, m_vLabelsUnfolded, m_ystarSpectRecoUnfoldedName );

  std::cout << "DONE! Closing " << fOut->GetName() << std::endl;
  fOut->Close(); delete fOut;
  std::cout << "......Closed  " << std::endl;
}

void DiJetAnalysisMC::ProcessPhysics(){

  LoadHistograms();

  // THIS HAS TO BE CHANGED
  TFile* fInMCPerf  = TFile::Open( m_fNamePerfUF.c_str() );
  
  TFile* fOut = new TFile( m_fNamePhys.c_str(),"RECREATE");

  // add a slice "all" to collection
  // rest of plots include combined slices
  m_vJznLabels.push_back( m_allName );
  
  m_hAllDphiReco  = CombineSamples( m_vHjznDphiReco, m_dPhiRecoName   );
  MakeDeltaPhi( m_vHjznDphiReco , m_vJznLabels, m_dPhiRecoName,
		fInMCPerf, m_ystarSpectRecoName );
  fOut->cd();

  m_hAllDphiTruth = CombineSamples( m_vHjznDphiTruth, m_dPhiTruthName );
  MakeDeltaPhi( m_vHjznDphiTruth, m_vJznLabels, m_dPhiTruthName,
		fInMCPerf, m_ystarSpectTruthName );
  fOut->cd();

  m_hAllDphiRecoPairedTruth = CombineSamples( m_vHjznDphiRecoPairedTruth, m_dPhiRecoPairedTruthName );
  /*
  MakeDeltaPhi( m_vHjznDphiRecoPairedTruth, m_vJznLabels, m_dPhiRecoPairedTruthName,
		fInMCPerf, m_ystarSpectTruthName );
  fOut->cd();
  */
  
  m_hAllDphiRespMat    = CombineSamples( m_vHjznDphiRespMat   , m_dPhiRespMatName    );
  m_hAllDphiRespMatReb = CombineSamples( m_vHjznDphiRespMatReb, m_dPhiRespMatRebName );
  
  MakeDphiCFactorsRespMat( m_vHjznDphiTruth, m_vHjznDphiRecoPairedTruth, m_vHjznDphiRespMatReb,
			   m_vJznLabels, m_dPhiCfactorsName, m_dPhiRespMatName );

  /*
    MakePtRespMat( m_vHjznDphiRespMatReb, m_vJznLabels, m_ptRespMatName );
  */

  std::cout << "DONE! Closing " << fOut->GetName() << std::endl;
  fOut->Close();
  delete fOut;
  std::cout << "......Closed  " << std::endl;
}

void DiJetAnalysisMC::UnfoldPhysics(){
  
  // for now
  // Copy File with original dPhi, spectra, etc,
  // into file where histos with corrections are
  // going to be appended. 
  std::cout << "Copy " << m_fNamePhys << " -> " << m_fNamePhysUF << std::endl;
  TFile::Cp( m_fNamePhys.c_str(), m_fNamePhysUF.c_str() );

  // Open two for reading one for updating.
  // open the MC file used for unfolding info.
  // open the data file used for measured info.
  // passed to unfolding function.
  TFile* fInMC     = TFile::Open( m_fNamePhys.c_str() );
  TFile* fInData   = TFile::Open( m_fNamePhys.c_str() );
  TFile* fInMCPerf = TFile::Open( m_fNamePerfUF.c_str() );
  TFile* fOut      = new TFile( m_fNamePhysUF.c_str(),"UPDATE");
  
  std::cout << "----- Unfolding MC ------" << std::endl;
  // make a vector with just the unfolded result.
  // this is to send it to MakeDeltaPhi(..) to have
  // unfolded results plotted separately
  std::vector< THnSparse*  > m_vHDphiUnfolded;
  std::vector< std::string > m_vLabelsUnfolded;

  // make unfolded THnSparse with similar naming convention
  // as the other histograms. At this point, don't care about
  // doing this for all triggers. Altohugh, this can be
  // repeated in a loop with m_allName subsitituted for trigger,
  // and subsequently added to the vectors above.  
  THnSparse* m_hAllDphiRecoUnfolded =
    UnfoldDeltaPhi( fInData, fInMC, m_dPhiRecoUnfoldedName,
		    fInMCPerf, m_ystarSpectRecoUnfoldedName );
  m_vHDphiUnfolded .push_back( m_hAllDphiRecoUnfolded );
  m_vLabelsUnfolded.push_back( m_allName );

  // unfold on MC, just used for testing purpose.
  // make deltaPhi, give flag (true) that its unfolded response
  // so there is no comb subt or normalization or scaling
  MakeDeltaPhi( m_vHDphiUnfolded, m_vLabelsUnfolded, m_dPhiRecoUnfoldedName,
	        fInMCPerf, m_ystarSpectRecoUnfoldedName );
  
  std::cout << "DONE! Closing " << fOut->GetName() << std::endl;
  fOut->Close(); delete fOut;
  std::cout << "......Closed  " << std::endl;
}

void DiJetAnalysisMC::MakeResultsTogether(){
  
  DiJetAnalysis::MakeResultsTogether();
  
  // buggy with filenames
  // leave this out for now
  // more for performance studies
  // MakeDphiRecoTruth();
}

//---------------------------------
//            Read Data
//---------------------------------

void DiJetAnalysisMC::SetupHistograms(){
  
  for( uint iG = 0; iG < m_nJzn; iG++){
    std::string jzn = m_vJznLabels[iG];
    
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
		  m_nPtSpectBins, m_ptSpectMin, m_ptSpectMax ) );
    m_vHjznEtaSpectReco.back()->GetXaxis()->
      Set( m_nVarFwdEtaBins, &( m_varFwdEtaBinning[0] ) );
    AddHistogram( m_vHjznEtaSpectReco.back() );
    
    m_vHjznEtaSpectTruth.push_back
      ( new TH2D( Form("h_%s_%s", m_etaSpectTruthName.c_str(), jzn.c_str() ), 
		  ";#eta_{Truth};#it{p}_{T}^{Truth} [GeV]",
		  m_nVarFwdEtaBins, 0, 1,
		  m_nPtSpectBins, m_ptSpectMin, m_ptSpectMax ) );
    m_vHjznEtaSpectTruth.back()->GetXaxis()->
      Set( m_nVarFwdEtaBins, &( m_varFwdEtaBinning[0] ) );
    AddHistogram( m_vHjznEtaSpectTruth.back() );

    m_vHjznYstarSpectReco.push_back
      ( new TH2D( Form("h_%s_%s", m_ystarSpectRecoName.c_str(), jzn.c_str() ), 
		  ";#it{y}_{1}*;#it{p}_{T1}^{Reco} [GeV]",
		  m_nVarYstarBins, 0, 1,
		  m_nVarPtBinsRespMat, 0, 1 ) );
    m_vHjznYstarSpectReco.back()->GetXaxis()->
      Set( m_nVarYstarBins, &( m_varYstarBinning[0] ) );
    m_vHjznYstarSpectReco.back()->GetYaxis()->
      Set( m_nVarPtBinsRespMat, &( m_varPtBinningRespMat[0] ) );
    AddHistogram( m_vHjznYstarSpectReco.back() );
    
    m_vHjznYstarSpectTruth.push_back
      ( new TH2D( Form("h_%s_%s", m_ystarSpectTruthName.c_str(), jzn.c_str() ), 
		  ";#it{y}_{1}*;#it{p}_{T1}^{Truth} [GeV]",
		  m_nVarYstarBins, 0, 1,
		  m_nVarPtBinsRespMat, 0, 1 ) );
    m_vHjznYstarSpectTruth.back()->GetXaxis()->
      Set( m_nVarYstarBins, &( m_varYstarBinning[0] ) );
    m_vHjznYstarSpectTruth.back()->GetYaxis()->
      Set( m_nVarPtBinsRespMat, &( m_varPtBinningRespMat[0] ) );
   AddHistogram( m_vHjznYstarSpectTruth.back() );

    // --- spectra response matrix ----
    m_vHjznYstarSpectRespMat.push_back
      ( new TH3D( Form("h_%s_%s", m_ystarSpectRespMatName.c_str(), jzn.c_str() ), 
		  ";#it{y}_{1}*;#it{p}_{T1}^{Reco} [GeV];#it{p}_{T1}^{Truth} [GeV]",
		  m_nVarYstarBins, 0, 1,
		  m_nVarPtBinsRespMat, 0, 1,
		  m_nVarPtBinsRespMat, 0, 1 ) );
    m_vHjznYstarSpectRespMat.back()->GetXaxis()->
      Set( m_nVarYstarBins, &( m_varYstarBinning[0] ) );
    m_vHjznYstarSpectRespMat.back()->GetYaxis()->
      Set( m_nVarPtBinsRespMat, &( m_varPtBinningRespMat[0] ) );
    m_vHjznYstarSpectRespMat.back()->GetZaxis()->
      Set( m_nVarPtBinsRespMat, &( m_varPtBinningRespMat[0] ) );
    AddHistogram( m_vHjznYstarSpectRespMat.back() );
    
    // --------- recoTruthRpt ---------
    m_vHjznRecoTruthRpt.push_back
      ( new TH3D( Form("h_recoTruthRpt_%s", jzn.c_str() ),
		  ";#it{y}_{Truth};#it{p}_{T}^{Truth};#it{p}_{T}^{Reco}/#it{p}_{T}^{Truth}",
		  m_nVarYstarBins, 0, 1,
		  m_nPtTruthBins,m_ptTruthMin, m_ptTruthMax,
		  m_nRPtRecoTruthBins, m_rPtRecoTruthMin, m_rPtRecoTruthMax) );
    m_vHjznRecoTruthRpt.back()->GetXaxis()->
      Set( m_nVarYstarBins, &( m_varYstarBinning[0] ) );
    AddHistogram( m_vHjznRecoTruthRpt.back() );
    
    m_vHjznRecoTruthRptNent.push_back
      ( new TH2D( Form("h_recoTruthRptNent_%s", jzn.c_str() ),
		  ";#it{y}_{Truth};#it{p}_{T}^{Truth}",
		  m_nVarYstarBins, 0, 1,
		  m_nPtTruthBins, m_ptTruthMin, m_ptTruthMax ) );
    m_vHjznRecoTruthRptNent.back()->GetXaxis()->
      Set( m_nVarYstarBins, &( m_varYstarBinning[0] ) );
    AddHistogram( m_vHjznRecoTruthRptNent.back() );    
    
    // --------- recoTruthDeta ---------
    m_vHjznRecoTruthDeta.push_back
      ( new TH3D( Form("h_recoTruthDeta_%s", jzn.c_str() ),
		  ";#it{y}_{Truth};#it{p}_{T}^{Truth};#eta^{Reco}-#eta^{Truth}",
		  m_nVarYstarBins, 0, 1,
		  m_nPtTruthBins,m_ptTruthMin, m_ptTruthMax,
		  m_nDAngleRecoTruthBins, m_dAngleRecoTruthMin, m_dAngleRecoTruthMax ) );
    m_vHjznRecoTruthDeta.back()->GetXaxis()->
      Set( m_nVarYstarBins, &( m_varYstarBinning[0] ) );
    AddHistogram( m_vHjznRecoTruthDeta.back() );    
    
    m_vHjznRecoTruthDetaNent.push_back
      ( new TH2D( Form("h_recoTruthDetaNent_%s", jzn.c_str() ),
		  ";#it{y}_{Truth};#it{p}_{T}^{Truth}",
		  m_nVarYstarBins, 0, 1,
		  m_nPtTruthBins, m_ptTruthMin, m_ptTruthMax ) );
    m_vHjznRecoTruthDetaNent.back()->GetXaxis()->
      Set( m_nVarYstarBins, &( m_varYstarBinning[0] ) );
    AddHistogram( m_vHjznRecoTruthDetaNent.back() );    

    // --------- recoTruthDphi ---------
    m_vHjznRecoTruthDphi.push_back
      ( new TH3D( Form("h_recoTruthDphi_%s", jzn.c_str() ),
		  ";#it{y}_{Truth};#it{p}_{T}^{Truth};#phi^{Reco}-#phi^{Truth}",
		  m_nVarYstarBins, 0, 1,
		  m_nPtTruthBins, m_ptTruthMin, m_ptTruthMax,
		  m_nDAngleRecoTruthBins, m_dAngleRecoTruthMin, m_dAngleRecoTruthMax ) );
    m_vHjznRecoTruthDphi.back()->GetXaxis()->
      Set( m_nVarYstarBins, &( m_varYstarBinning[0] ) );
    AddHistogram( m_vHjznRecoTruthDphi.back() );
    
    m_vHjznRecoTruthDphiNent.push_back
      ( new TH2D( Form("h_recoTruthDphiNent_%s", jzn.c_str() ),
		  ";#it{y}_{Truth};#it{p}_{T}^{Truth}",
		  m_nVarYstarBins, 0, 1,
		  m_nPtTruthBins, m_ptTruthMin, m_ptTruthMax ) );
    m_vHjznRecoTruthDphiNent.back()->GetXaxis()->
      Set( m_nVarYstarBins, &( m_varYstarBinning[0] ) );
    AddHistogram( m_vHjznRecoTruthDphiNent.back() );

    // -------- dPhi --------
    THnSparse* hnDphiReco =
      new THnSparseD( Form("h_%s_%s", m_dPhiRecoName.c_str(), jzn.c_str() ), "",
		      m_nDphiDim, &m_vNdPhiBins[0], &m_vDphiMin[0], &m_vDphiMax[0] );
    hnDphiReco->GetAxis(0)->Set( m_nVarYstarBins, &( m_varYstarBinning[0] ) );
    hnDphiReco->GetAxis(1)->Set( m_nVarYstarBins, &( m_varYstarBinning[0] ) );
    hnDphiReco->GetAxis(2)->Set( m_nVarPtBins   , &( m_varPtBinning[0]    ) );
    hnDphiReco->GetAxis(3)->Set( m_nVarPtBins   , &( m_varPtBinning[0]    ) );
    hnDphiReco->GetAxis(4)->Set( m_nVarDphiBins , &( m_varDphiBinning[0]  ) );
    m_vHjznDphiReco.push_back( hnDphiReco );
    AddHistogram( hnDphiReco );
    
    THnSparse* hnDphiTruth =
      new THnSparseD( Form("h_%s_%s", m_dPhiTruthName.c_str(), jzn.c_str() ), "",
		      m_nDphiDim, &m_vNdPhiBins[0], &m_vDphiMin[0], &m_vDphiMax[0] );
    hnDphiTruth->GetAxis(0)->Set( m_nVarYstarBins, &( m_varYstarBinning[0] ) );
    hnDphiTruth->GetAxis(1)->Set( m_nVarYstarBins, &( m_varYstarBinning[0] ) );
    hnDphiTruth->GetAxis(2)->Set( m_nVarPtBins   , &( m_varPtBinning[0]    ) );
    hnDphiTruth->GetAxis(3)->Set( m_nVarPtBins   , &( m_varPtBinning[0]    ) );
    hnDphiTruth->GetAxis(4)->Set( m_nVarDphiBins , &( m_varDphiBinning[0]  ) );
    m_vHjznDphiTruth.push_back( hnDphiTruth );
    AddHistogram( hnDphiTruth );

    THnSparse* hnDphiRecoPairedTruth =
      new THnSparseD( Form("h_%s_%s", m_dPhiRecoPairedTruthName.c_str(), jzn.c_str() ), "",
		      m_nDphiDim, &m_vNdPhiBins[0], &m_vDphiMin[0], &m_vDphiMax[0] );
    hnDphiRecoPairedTruth->GetAxis(0)->Set( m_nVarYstarBins, &( m_varYstarBinning[0] ) );
    hnDphiRecoPairedTruth->GetAxis(1)->Set( m_nVarYstarBins, &( m_varYstarBinning[0] ) );
    hnDphiRecoPairedTruth->GetAxis(2)->Set( m_nVarPtBins   , &( m_varPtBinning[0]    ) );
    hnDphiRecoPairedTruth->GetAxis(3)->Set( m_nVarPtBins   , &( m_varPtBinning[0]    ) );
    hnDphiRecoPairedTruth->GetAxis(4)->Set( m_nVarDphiBins , &( m_varDphiBinning[0]  ) );
    m_vHjznDphiRecoPairedTruth.push_back( hnDphiRecoPairedTruth );
    AddHistogram( hnDphiRecoPairedTruth );
    
    // -------- Dphi Response Matrix --------    
    THnSparse* hnDphiRespMat =
      new THnSparseD( Form("h_%s_%s", m_dPhiRespMatName.c_str(), jzn.c_str() ), "",
		      m_nDphiRespMatDim, &m_vNdPhiRespMatBins[0],
		      &m_vDphiRespMatMin[0], &m_vDphiRespMatMax[0] );
    hnDphiRespMat->GetAxis(0)->Set( m_nVarYstarBins    , &( m_varYstarBinning[0]     ) );
    hnDphiRespMat->GetAxis(1)->Set( m_nVarYstarBins    , &( m_varYstarBinning[0]     ) );
    hnDphiRespMat->GetAxis(2)->Set( m_nVarPtBinsRespMat, &( m_varPtBinningRespMat[0] ) );
    hnDphiRespMat->GetAxis(3)->Set( m_nVarPtBinsRespMat, &( m_varPtBinningRespMat[0] ) );
    hnDphiRespMat->GetAxis(4)->Set( m_nVarPtBinsRespMat, &( m_varPtBinningRespMat[0] ) );
    hnDphiRespMat->GetAxis(5)->Set( m_nVarPtBinsRespMat, &( m_varPtBinningRespMat[0] ) );
    hnDphiRespMat->GetAxis(6)->Set( m_nVarDphiBins     , &( m_varDphiBinning[0]      ) );
    hnDphiRespMat->GetAxis(7)->Set( m_nVarDphiBins     , &( m_varDphiBinning[0]      ) );

    m_vHjznDphiRespMat.push_back( hnDphiRespMat );
    AddHistogram( hnDphiRespMat );
    
    // change some of the axis names
    hnDphiRespMat->GetAxis(2)->SetTitle( "Reco #it{p}_{T}^{1}"  );
    hnDphiRespMat->GetAxis(3)->SetTitle( "Truth #it{p}_{T}^{1}" );
    hnDphiRespMat->GetAxis(4)->SetTitle( "Reco #it{p}_{T}^{2}"  );
    hnDphiRespMat->GetAxis(5)->SetTitle( "Truth #it{p}_{T}^{2}" );
    hnDphiRespMat->GetAxis(6)->SetTitle( "|#Delta#phi_{Reco}|"  );
    hnDphiRespMat->GetAxis(7)->SetTitle( "|#Delta#phi_{Truth}|" );
    
    // --- Rebinned Dphi Response Matrix -----    
    THnSparse* hnDphiRespMatReb =
      new THnSparseD( Form("h_%s_%s", m_dPhiRespMatRebName.c_str(), jzn.c_str() ), "",
		      m_nDphiRespMatDim, &m_vNdPhiRespMatRebBins[0],
		      &m_vDphiRespMatMin[0], &m_vDphiRespMatMax[0] );
    hnDphiRespMatReb->GetAxis(0)->Set( m_nVarYstarBins    , &( m_varYstarBinning[0]     ) );
    hnDphiRespMatReb->GetAxis(1)->Set( m_nVarYstarBins    , &( m_varYstarBinning[0]     ) );
    hnDphiRespMatReb->GetAxis(2)->Set( m_nVarPtBinsRespMat, &( m_varPtBinningRespMat[0] ) );
    hnDphiRespMatReb->GetAxis(3)->Set( m_nVarPtBinsRespMat, &( m_varPtBinningRespMat[0] ) );
    hnDphiRespMatReb->GetAxis(4)->Set( m_nVarPtBinsRespMat, &( m_varPtBinningRespMat[0] ) );
    hnDphiRespMatReb->GetAxis(5)->Set( m_nVarPtBinsRespMat, &( m_varPtBinningRespMat[0] ) );
    hnDphiRespMatReb->GetAxis(6)->Set( m_nVarDphiRebinnedBins, &( m_varDphiRebinnedBinning[0] ) );
    hnDphiRespMatReb->GetAxis(7)->Set( m_nVarDphiRebinnedBins, &( m_varDphiRebinnedBinning[0] ) );

    m_vHjznDphiRespMatReb.push_back( hnDphiRespMatReb );
    AddHistogram( hnDphiRespMatReb );
    
    // change some of the axis names
    hnDphiRespMatReb->GetAxis(2)->SetTitle( "Reco #it{p}_{T}^{1}"  );
    hnDphiRespMatReb->GetAxis(3)->SetTitle( "Truth #it{p}_{T}^{1}" );
    hnDphiRespMatReb->GetAxis(4)->SetTitle( "Reco #it{p}_{T}^{2}"  );
    hnDphiRespMatReb->GetAxis(5)->SetTitle( "Truth #it{p}_{T}^{2}" );
    hnDphiRespMatReb->GetAxis(6)->SetTitle( "|#Delta#phi_{Reco}|"  );
    hnDphiRespMatReb->GetAxis(7)->SetTitle( "|#Delta#phi_{Truth}|" );    
  } 
}

void DiJetAnalysisMC::ProcessEvents( int nEventsIn, int startEventIn ){ 
  // collections and variables
  std::vector< TLorentzVector >    vT_jets;
  std::vector< TLorentzVector >* p_vT_jets = &vT_jets;
  
  std::vector< TLorentzVector >    vR_jets;
  std::vector< TLorentzVector >* p_vR_jets = &vR_jets;
  
  std::vector< bool >    v_isCleanJet;
  std::vector< bool >* p_v_isCleanJet = &v_isCleanJet;

  std::vector< std::vector< float > >    v_sysUncert;
  std::vector< std::vector< float > >* p_v_sysUncert = & v_sysUncert;
  
  for( uint iG = 0; iG < m_nJzn; iG++){
     
    std::cout << "fNameIn: " << m_vJznFnameIn[iG] << std::endl;
  
    TFile* fIn  = TFile::Open( m_vJznFnameIn[iG].c_str() );
    TTree* tree = static_cast< TTree*>( fIn->Get( "tree" ) );

    // Connect to tree
    tree->SetBranchAddress( "vT_jets"      , &p_vT_jets      );
    tree->SetBranchAddress( "vR_C_jets"    , &p_vR_jets      );
    tree->SetBranchAddress( "v_isCleanJet" , &p_v_isCleanJet );
    tree->SetBranchAddress( "v_sysUncert"  , &p_v_sysUncert  );

    // n events
    int nEvents, startEvent, nEventsTotal, endEvent;
    
    nEventsTotal = tree->GetEntries();
    nEvents      = nEventsIn > 0 ? nEventsIn : nEventsTotal;
    startEvent   = startEventIn < nEventsTotal ?
				  startEventIn : nEventsTotal - 1;
    endEvent     = startEvent + nEvents < nEventsTotal ?
					  startEvent + nEvents : nEventsTotal;

    std::cout << startEvent << " " << endEvent << " " << nEventsTotal << " " << nEvents << std::endl;
    
    // -------- EVENT LOOP ---------
    for( m_ev = startEvent; m_ev < endEvent; m_ev++ ){
      tree->GetEntry( m_ev );
            
      if( anaTool->DoPrint(m_ev) ) {
	std::cout << "\nEvent : " << m_ev 
		  << "    has : " << vR_jets.size() << " reco jets"
		  << "    and : " << vT_jets.size() << " truth jets"
		  << std::endl; 
      }

      ApplyCleaning ( vR_jets, v_isCleanJet );
      ApplyIsolation( vR_jets, 1.0 );
      ApplyIsolation( vT_jets, 1.0 );
      
      std::sort( vR_jets.begin(), vR_jets.end(), anaTool->sortByDecendingPt );

      std::vector< TLorentzVector > vRR_paired_jets;
      std::vector< TLorentzVector > vRT_paired_jets;
      PairJets( vR_jets, vT_jets, vRR_paired_jets, vRT_paired_jets ); 

      std::vector< TLorentzVector > vTR_paired_jets;
      std::vector< TLorentzVector > vTT_paired_jets;
      PairJets( vT_jets, vR_jets, vTT_paired_jets, vTR_paired_jets );

      // this holds the weights.
      // in some cases ( unfolding uncertainty )
      // there will be an additional weight.
      std::vector< double > vRR_paired_jet_weights;

      // If not running on default sample.
      // Apply uncertainties to all reco jets.
      if( m_uncertComp ){
	m_uncertaintyProvider->RegisterUFactors  ( &v_sysUncert );
	m_uncertaintyProvider->ApplyUncertainties( vRR_paired_jets, vRT_paired_jets );
      }
      
      // Do Dphi analysis
      AnalyzeDeltaPhi( m_vHjznDphiReco [iG], vR_jets );
      AnalyzeDeltaPhi( m_vHjznDphiTruth[iG], vT_jets );
      AnalyzeDeltaPhi( m_vHjznDphiRecoPairedTruth [iG], vRR_paired_jets );

      AnalyzeDphiRespMat
	( m_vHjznDphiRespMat[iG], m_vHjznDphiRespMatReb[iG], vTR_paired_jets, vTT_paired_jets );

      // loop over truth jets
      // denominator for efficiency because not all
      // reco jets are reconstructed for a truth jet
      // make spectra for truth jets
      AnalyzeSpectra( m_vHjznEtaSpectTruth[iG], vT_jets );
      // make spectra for the reco jets paired to truth
      AnalyzeSpectra( m_vHjznEtaSpectReco [iG], vTR_paired_jets );
      
      // fill single jet spectra 
      if( vR_jets.size() ){
	m_vHjznYstarSpectReco[iG]->
	  Fill( GetYstar( vR_jets.front() ), vR_jets.front().Pt()/1000. );
	// fill both sides for pp
	if( !m_is_pPb ){
	  m_vHjznYstarSpectReco[iG]->
	    Fill( -GetYstar( vR_jets.front() ), vR_jets.front().Pt()/1000. );
	}
      }
      if( vT_jets.size() ){
	m_vHjznYstarSpectTruth[iG]->
	  Fill( GetYstar( vT_jets.front() ), vT_jets.front().Pt()/1000. );
	// fill both sides for pp
	if( !m_is_pPb ){
	m_vHjznYstarSpectTruth[iG]->
	  Fill( -GetYstar( vT_jets.front() ), vT_jets.front().Pt()/1000. );
	}

      }

      // fill single speectra response matrix
      AnalyzeSpectRespMat( m_vHjznYstarSpectRespMat[iG], vR_jets, vT_jets );
      
      // do JER/JES, angular scales and resolution.
      AnalyzeScaleResolution( vTR_paired_jets, vTT_paired_jets, iG );
    } // -------- END EVENT LOOP ---------
   
    std::cout << "DONE WITH " << m_vJznLabels[iG] << std::endl;

    fIn->Close(); delete fIn;
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
	
    double jetWeightReco  = GetJetWeight( recoJet  );
    double jetWeightTruth = GetJetWeight( truthJet );
	
    m_vHjznEtaPhiMap[iG]->Fill( jetEtaReco, jetPhiReco, jetWeightReco );
    m_vHjznEtaPtMap [iG]->Fill( jetEtaReco, jetPtReco , jetWeightReco );

    m_vHjznRecoTruthRpt       [iG]->
      Fill( jetEtaTruth, jetPtTruth, jetPtReco/jetPtTruth, jetWeightTruth );
    m_vHjznRecoTruthRptNent   [iG]->
      Fill( jetEtaTruth, jetPtTruth );

    m_vHjznRecoTruthDeta      [iG]->
      Fill( jetEtaTruth, jetPtTruth, jetEtaReco - jetEtaTruth, jetWeightTruth );
    m_vHjznRecoTruthDetaNent  [iG]->
      Fill( jetEtaTruth, jetPtTruth );

    m_vHjznRecoTruthDphi      [iG]->
      Fill( jetEtaTruth, jetPtTruth, jetPhiReco - jetPhiTruth, jetWeightTruth );
    m_vHjznRecoTruthDphiNent  [iG]->
      Fill( jetEtaTruth, jetPtTruth  );
  } // end loop over pairs
}

void DiJetAnalysisMC::AnalyzeSpectRespMat( TH3* hRespMat,
					   const std::vector< TLorentzVector >& vR_jets,
					   const std::vector< TLorentzVector >& vT_jets ){
  
  if( vR_jets.size() && vT_jets.size() ){

    const TLorentzVector* recoJet  = &vR_jets.front();
    const TLorentzVector* truthJet = &vT_jets.front();

    
    double recoJet_pt     = recoJet->Pt()/1000.;
  
    double truthJet_pt    = truthJet->Pt()/1000.;
    double truthJet_ystar = GetYstar( *truthJet );
    
    hRespMat->Fill(  truthJet_ystar, recoJet_pt, truthJet_pt );

    // for pp fill plus minus ystar
    if( m_is_pPb ){ return; }

    hRespMat->Fill( -truthJet_ystar, recoJet_pt, truthJet_pt );
  }
}

void DiJetAnalysisMC::AnalyzeDphiRespMat( THnSparse* hnDphi,
					  THnSparse* hnDphiReb,
					  const std::vector< TLorentzVector >& vR_jets,
					  const std::vector< TLorentzVector >& vT_jets ){

  const TLorentzVector* recoJet1  = NULL; const TLorentzVector* recoJet2  = NULL;
  const TLorentzVector* truthJet1 = NULL; const TLorentzVector* truthJet2 = NULL;

  if( !GetDiJets( vR_jets, recoJet1 , recoJet2  ) )
    { return; }
  if( !GetDiJets( vT_jets, truthJet1, truthJet2 ) )
    { return; }

  double recoJet1_pt    = recoJet1->Pt()/1000.;
  
  double recoJet2_pt    = recoJet2->Pt()/1000.;
  
  double truthJet1_pt    = truthJet1->Pt()/1000.;
  double truthJet1_ystar = GetYstar( *truthJet1 );

  double truthJet2_pt    = truthJet2->Pt()/1000.;
  double truthJet2_ystar = GetYstar( *truthJet2 );
  
  double recoDeltaPhi    = anaTool->DeltaPhi( *recoJet2 , *recoJet1  );
  double truthDeltaPhi   = anaTool->DeltaPhi( *truthJet2, *truthJet1 );
  
  std::vector< double > xDphi( hnDphi->GetNdimensions(), 0 );
    
  double jetWeight = GetJetWeight( *truthJet1 );

  // fill for Dphi resp mat
  xDphi[0] = truthJet1_ystar;  
  xDphi[1] = truthJet2_ystar;
  xDphi[2] = recoJet1_pt ;
  xDphi[3] = truthJet1_pt;
  xDphi[4] = recoJet2_pt ;
  xDphi[5] = truthJet2_pt;
  xDphi[6] = recoDeltaPhi;
  xDphi[7] = truthDeltaPhi;
  hnDphi   ->Fill( &xDphi[0], jetWeight );
  hnDphiReb->Fill( &xDphi[0], jetWeight );
  
  // for pp, fill twice. once for each side since
  // it is symmetric in pp. For pPb, continue
  if( m_is_pPb ){ return; }
  
  // fill for Dphi resp mat
  xDphi[0] = -truthJet1_ystar;  
  xDphi[1] = -truthJet2_ystar;
  hnDphi   ->Fill( &xDphi[0], jetWeight );
  hnDphiReb->Fill( &xDphi[0], jetWeight );
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

    /*
    // check if this jet is in some useful pt range
    // here we call that useful range the range of the
    // jer jer plots
    double aJetPt    = aJet.Pt()/1000.;
    if( aJetPt <  m_ptTruthMin ||
	aJetPt >= m_ptTruthMax )
      { continue; }
    */
    
    // also ignore jets outside of a certain ystar range.
    double aJetYstar = GetYstar( aJet );
    if( aJetYstar < m_varYstarBinning.front() ||
	aJetYstar > m_varYstarBinning.back() )
      { continue; }
      
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

TH3* DiJetAnalysisMC::CombineSamples( std::vector< TH3* >& vSampleHin,
				      const std::string& name ){
  if( !vSampleHin.size() ){ return NULL; }

  TH3* h_res = static_cast< TH3D* >
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

  // loop over xBins
  for( int xBin = 1; xBin <= h_res->GetNbinsX(); xBin++ ){

    double valTot      = 0;
    double valErrorTot = 0;
    double denomTot    = 0;

    for( uint iG = 0; iG < vSampleHin.size(); iG++ ){
      if( vSampleHin[iG]->GetEntries() == 0 ){ continue; }
      
      double valueBin    = vSampleHin   [iG]->GetBinContent( xBin );
      double valueBinErr = vSampleHin   [iG]->GetBinError( xBin );
      double nEntriesBin = vSampleNentIn[iG]->GetBinContent( xBin );
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
  } // end loop over xBins
}

void DiJetAnalysisMC::CombineSamples( TH1* h_res,
				      std::vector< TH1* >& vSampleHin,
				      std::vector< TH1* >& vSampleNentInT,
				      std::vector< TH1* >& vSampleNentInR,
				      std::vector< TH2* >& vSampleRespMat,
				      const std::string& name ){

  // loop over xBins
  for( int xBin = 1; xBin <= h_res->GetNbinsX(); xBin++ ){

    double valTot      = 0;
    double valErrorTot = 0;
    double denomTot    = 0;

    for( uint iG = 0; iG < vSampleHin.size(); iG++ ){
      if( vSampleHin[iG]->GetEntries() == 0 ){ continue; }

      double valueBin = vSampleHin    [iG]->GetBinContent( xBin );
      double nEntT    = vSampleNentInT[iG]->GetBinContent( xBin );
      double nEntR    = vSampleNentInR[iG]->GetBinContent( xBin );
      double nEntMat  = vSampleRespMat[iG]->GetBinContent( xBin );
      
      double weight   = m_vJznWeights[iG] / m_vJznSumOverlayWeights[iG];

      double newError =
	std::sqrt( std::pow( nEntT, 2 ) / std::pow( nEntR, 3 ) *
		   ( 1 - std::pow( nEntMat, 2 ) / ( nEntT * nEntR ) ) );

      // if one of them is zero, leave it alone
      if( nEntT == 0 || nEntR == 0 ){
	newError = 0;
      }

      double val    = weight * nEntR * valueBin;
      double valErr = weight * nEntR * newError;

      valTot       += val;
      valErrorTot  += valErr * valErr;

      double denom  = weight * nEntR;
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
  } // end loop over xBins
}

void DiJetAnalysisMC::SetCfactorsErrors( TH1* hR, TH1* hT, TH2* hM, TH1* hC ){

  std::cout << hT->GetName() << std::endl;
  // set factors
  for( int xBin = 1; xBin <= hT->GetNbinsX(); xBin++ ){
    int cBin = hC->FindBin( hT->GetBinCenter(xBin) );

    // if one of them is zero, leave it alone
    if( hT->GetBinContent( xBin ) == 0 ||
	hR->GetBinContent( xBin ) == 0 ){
      hC->SetBinContent( cBin, 1 );
      hC->SetBinError  ( cBin, 0 );
      continue;
    }
	      		
    double vR = hR->GetBinContent( xBin  );
    double vT = hT->GetBinContent( xBin  );

    double vM = hM->GetBinContent( cBin, cBin );
	     
    // error on correction factor    
    double newDphiError =
      std::sqrt( std::pow( vT, 2 ) / std::pow( vR, 3 ) *
		 ( 1 - std::pow( vM, 2 ) / ( vT * vR ) ) );

    std::cout << xBin << " " << vT << " " << vR << " "
	      << vM << "  C = " << hC->GetBinContent( xBin )
	      << "   err = " << newDphiError << std::endl;
    
    hC->SetBinError( cBin, newDphiError );
  }
}

double DiJetAnalysisMC::GetJetWeight( const TLorentzVector& jet ){
  if( !m_hPowhegWeights ) { return 1; } 
  
  int xb = m_hPowhegWeights->GetXaxis()->FindBin( jet.Pt()/1000. );
  int yb = m_hPowhegWeights->GetYaxis()->FindBin( jet.Eta()      );
  int zb = m_hPowhegWeights->GetZaxis()->FindBin( jet.Phi()      );
  float jet_weight = m_hPowhegWeights->GetBinContent( xb, yb, zb);

  return jet_weight;
}

double DiJetAnalysisMC::GetUncertaintyWeight( const TLorentzVector& jet1,
					      const TLorentzVector& jet2 )
{

  return m_uncertaintyProvider->GetUncertaintyWeight( jet1, jet2 );
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

void DiJetAnalysisMC::GetSpectUnfoldingInfo( std::string& measuredName,
					     std::string& recoName,
					     std::string& truthName,
					     std::string& respMatName,
					     std::string& unfoldedLabel,
					     std::string& typeLabel ){

  measuredName  = m_ystarSpectRecoName;
  recoName      = m_ystarSpectRecoName;
  truthName     = m_ystarSpectTruthName;
  respMatName   = m_ystarSpectRespMatName;
  unfoldedLabel = "dN/d#it{p}_{T} [GeV]";;
  typeLabel     = "MC";
}

void DiJetAnalysisMC::GetDphiUnfoldingInfo( std::string& measuredName,
					    std::string& truthName,
					    std::string& unfoldedLabel,
					    std::string& typeLabel ){

  measuredName  = m_dPhiRecoName;
  truthName     = m_dPhiTruthName;
  unfoldedLabel = "|#Delta#phi|";
  typeLabel     = "MC";
}

void DiJetAnalysisMC::GetInfoTogether( std::string& name_a , std::string& name_b ,
				       std::string& label_a, std::string& label_b,
				       std::string& fName_a, std::string& fName_b,
				       int option ){

  std::string* pFname = NULL;
  
  switch( option ){
  case 0:
    pFname = &m_fNamePerfUF;
    name_a = m_ystarSpectName;
    name_b = m_ystarSpectName;
    break;
  case 1:
    pFname = &m_fNamePhysUF;
    name_a = m_dPhiName;
    name_b = m_dPhiName;
    break;
  }
  
  int combinationBoth = GetConfig()->GetValue( "combinationBoth", 0 );

  if( combinationBoth == 0 ){
    name_a   +=  "_" + m_recoName;
    name_b   +=  "_" + m_truthName;
    label_a  = "Reco_{MC}";
    label_b  = "Truth_{MC}";
    fName_a  = *pFname;
    fName_b  = *pFname;
  } else if ( combinationBoth == 1 ){
    name_a   +=  "_" + m_recoName;
    name_b   +=  "_" + m_recoName;
    label_a  = "Reco_{MC}";
    label_b  = "UF_{MC}";
    fName_a  = *pFname;
    fName_b  = *pFname;
  } else if ( combinationBoth == 2 ){
    name_a   +=  "_" + m_recoName + "_" + m_unfoldedName;
    name_b   +=  "_" + m_recoName + "_" + m_unfoldedName;
    label_a  = "#it{p}+Pb";
    label_b  = "#it{pp}";
    m_is_pPb = true;  DiJetAnalysis::Initialize();
    fName_a  = *pFname;
    m_is_pPb = false; DiJetAnalysis::Initialize();
    fName_b  = *pFname;
  } else if ( combinationBoth == 3 ){
    name_a   +=  "_" + m_truthName;
    name_b   +=  "_" + m_truthName;
    label_a  = "#it{p}+Pb";
    label_b  = "#it{pp}";
    m_is_pPb = true;  DiJetAnalysis::Initialize();
    fName_a  = *pFname;
    m_is_pPb = false; DiJetAnalysis::Initialize();
    fName_b  = *pFname;
  }
}

void DiJetAnalysisMC::GetPurityEff( TH2* hRespMat, TH1* hPurity, TH1* hEff ){

  for( int xBin = 1 ; xBin <= hRespMat->GetNbinsX(); xBin++ ){
    double sum = 0;
    for( int yBin = 1 ; yBin <= hRespMat->GetNbinsY(); yBin++ ){
      // if( xBin == yBin ){ continue; }
      sum += hRespMat->GetBinContent( xBin, yBin );
    }
    double eff = sum ? hRespMat->GetBinContent( xBin, xBin )/sum : 0 ;    
    hEff->SetBinContent( xBin, eff );
  }

  for( int yBin = 1 ; yBin <= hRespMat->GetNbinsY(); yBin++ ){
    double sum = 0;
    for( int xBin = 1 ; xBin <= hRespMat->GetNbinsX(); xBin++ ){
      // if( xBin == yBin ){ continue; }
      sum += hRespMat->GetBinContent( xBin, yBin );
    }
    double purity = sum ? hRespMat->GetBinContent( yBin, yBin )/sum : 0 ;    
    hPurity->SetBinContent( yBin, purity );
  }
}

//---------------------------------
//     Get Quantities / Plot 
//---------------------------------
void DiJetAnalysisMC::LoadHistograms(){

  TFile* fIn = TFile::Open( m_fNameRaw.c_str() ); 
  
  for( uint iG = 0; iG < m_nJzn; iG++){
    std::string jzn = m_vJznLabels[iG];
    
    // -------- maps ---------
    m_vHjznEtaPhiMap.push_back
      ( static_cast< TH2D* >
	( fIn->Get( Form("h_etaPhiMap_%s", jzn.c_str() ))));
    m_vHjznEtaPhiMap.back()->SetDirectory(0);

    m_vHjznEtaPtMap.push_back
      ( static_cast< TH2D* >
	( fIn->Get( Form("h_etaPtMap_%s", jzn.c_str() ))));
    m_vHjznEtaPtMap.back()->SetDirectory(0);

    // -------- spect --------
    m_vHjznEtaSpectReco.push_back
      ( static_cast< TH2D* >
	( fIn->Get
	  ( Form("h_%s_%s", m_etaSpectRecoName.c_str(), jzn.c_str() ))));
    m_vHjznEtaSpectReco.back()->SetDirectory(0);

    m_vHjznEtaSpectTruth.push_back 
      ( static_cast< TH2D* >
	( fIn->
	  Get( Form("h_%s_%s", m_etaSpectTruthName.c_str(), jzn.c_str() ))));
    m_vHjznEtaSpectTruth.back()->SetDirectory(0);

    m_vHjznYstarSpectReco.push_back
      ( static_cast< TH2D* >
	( fIn->Get
	  ( Form("h_%s_%s", m_ystarSpectRecoName.c_str(), jzn.c_str() ))));
    m_vHjznYstarSpectReco.back()->SetDirectory(0);

    m_vHjznYstarSpectTruth.push_back 
      ( static_cast< TH2D* >
	( fIn->
	  Get( Form("h_%s_%s", m_ystarSpectTruthName.c_str(), jzn.c_str() ))));
    m_vHjznYstarSpectTruth.back()->SetDirectory(0);

    // --- spectra response matrix ----
    m_vHjznYstarSpectRespMat.push_back 
      ( static_cast< TH3D* >
	( fIn->
	  Get( Form("h_%s_%s", m_ystarSpectRespMatName.c_str(), jzn.c_str() ))));
    m_vHjznYstarSpectRespMat.back()->SetDirectory(0);
    
    // --------- recoTruthRpt ---------
    m_vHjznRecoTruthRpt.push_back
      ( static_cast< TH3D* >
	( fIn->Get( Form("h_recoTruthRpt_%s", jzn.c_str() ))));
    m_vHjznRecoTruthRpt.back()->SetDirectory(0);

    m_vHjznRecoTruthRptNent.push_back
      ( static_cast< TH2D* >
	( fIn->Get( Form("h_recoTruthRptNent_%s", jzn.c_str() ))));
    m_vHjznRecoTruthRptNent.back()->SetDirectory(0);
    
    // --------- recoTruthDeta ---------
    m_vHjznRecoTruthDeta.push_back
      ( static_cast< TH3D* >
	( fIn->Get( Form("h_recoTruthDeta_%s", jzn.c_str() ))));
    m_vHjznRecoTruthDeta.back()->SetDirectory(0);

    m_vHjznRecoTruthDetaNent.push_back
      ( static_cast< TH2D* >
	( fIn->Get( Form("h_recoTruthDetaNent_%s", jzn.c_str() ))));
    m_vHjznRecoTruthDetaNent.back()->SetDirectory(0);
    
    // --------- recoTruthDphi ---------
    m_vHjznRecoTruthDphi.push_back
      ( static_cast< TH3D* >
	( fIn->Get( Form("h_recoTruthDphi_%s", jzn.c_str() )))); 
    m_vHjznRecoTruthDphi.back()->SetDirectory(0);

    m_vHjznRecoTruthDphiNent.push_back
      ( static_cast< TH2D* >
	( fIn->Get( Form("h_recoTruthDphiNent_%s", jzn.c_str() ))));
    m_vHjznRecoTruthDphiNent.back()->SetDirectory(0);

    // -------- dPhi- --------
    m_vHjznDphiReco.push_back
      ( static_cast< THnSparse *>
	( fIn->Get( Form("h_%s_%s", m_dPhiRecoName.c_str(), jzn.c_str() ))));  
    	
    m_vHjznDphiTruth.push_back
      ( static_cast< THnSparse *>
	( fIn->Get( Form("h_%s_%s", m_dPhiTruthName.c_str(), jzn.c_str() ))));  
    
    m_vHjznDphiRecoPairedTruth.push_back
      ( static_cast< THnSparse *>
	( fIn->Get( Form("h_%s_%s", m_dPhiRecoPairedTruthName.c_str(), jzn.c_str() ))));  
          
    // -------- Dphi Response Matrix --------    
    m_vHjznDphiRespMat.push_back
      ( static_cast< THnSparse *>
	( fIn->Get( Form("h_%s_%s", m_dPhiRespMatName.c_str(), jzn.c_str() ))));  
      
    // -------- Dphi Response Matrix --------    
    m_vHjznDphiRespMatReb.push_back
      ( static_cast< THnSparse *>
	( fIn->Get( Form("h_%s_%s", m_dPhiRespMatRebName.c_str(), jzn.c_str() ))));  

  }
  
  fIn->Close(); delete fIn;
}

void DiJetAnalysisMC::MakeScaleRes( std::vector< TH3* >& vJznHin,
				    std::vector< TH2* >& vJznNentIn,
				    const std::string& type ){

  // check if we have one
  if( !vJznHin.size() ){ return; }

  // is same for all the other histos
  int  nBinsX = vJznHin.front()->GetNbinsX();
  double xMin = vJznHin.front()->GetXaxis()->GetXmin();
  double xMax = vJznHin.front()->GetXaxis()->GetXmax();

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
    std::string jzn  = m_vJznLabels[iG];
    TH3* hJznHin     = vJznHin[iG];
   
    for( int xBin = 1; xBin <= nBinsX; xBin++ ){
      double xLow, xUp;
      anaTool->GetBinRange
	( hJznHin->GetXaxis(), xBin, xBin, xLow, xUp );
      
      std::string xAxisTitle = hJznHin->GetYaxis()->GetTitle();
      
      // build mean, sigma, project nev
      TH1* hMean = new TH1D
	( Form("h_%s_%s_%s_%s",
	       type.c_str(), sMean.c_str(), jzn.c_str(),
	       anaTool->GetName(xLow, xUp, "Ystar").c_str() ),
	  Form("%s;%s;%s",
	       anaTool->GetYstarLabel( xLow, xUp, 1 ).c_str(),
	       xAxisTitle.c_str(), yTitleMean.c_str() ),
	  nBinsY, yMin, yMax );
      vMeans[iG].push_back( hMean );
      
      TH1* hSigma = new TH1D
	( Form("h_%s_%s_%s_%s",
	       type.c_str(), sSigma.c_str(), jzn.c_str(),
	       anaTool->GetName(xLow, xUp, "Ystar").c_str()),
	  Form("%s;%s;%s",
	       anaTool->GetYstarLabel( xLow, xUp, 1 ).c_str(),
	       xAxisTitle.c_str(), yTitleSigma.c_str() ),
	  nBinsY, yMin, yMax );
      vSigmas[iG].push_back( hSigma );
      
      TH1* hNent = vJznNentIn[iG]->ProjectionY
	( Form("h_%s_%s_N_%s", type.c_str(), jzn.c_str(),
	       anaTool->GetName(xLow, xUp, "Ystar").c_str() )
	  , xBin, xBin );
      vNent[iG].push_back( hNent );

      // now loop over the second axis and project onto z axis
      // to get the gaussian distributions. fit them.
      vProj[iG].resize( nBinsX );
      vFit [iG].resize( nBinsX );
      
      for( int yBin = 1; yBin <= nBinsY; yBin++ ){
	double yLow, yUp;
	anaTool->GetBinRange
	  ( hJznHin->GetYaxis(), yBin, yBin, yLow, yUp );

	int iX = xBin - 1;

	TH1* hProj = hJznHin->ProjectionZ
	  ( Form("h_%s_%s_%s_%s",  type.c_str(), jzn.c_str(), 
		 anaTool->GetName(xLow, xUp, "Ystar").c_str(),
		 anaTool->GetName(yLow, yUp, "Pt").c_str() ),
	    xBin, xBin, yBin, yBin );

	styleTool->SetHStyle( hProj, 0 );
	vProj[iG][iX].push_back( hProj );

 	TF1* fit = anaTool->FitGaussian( hProj );
	styleTool->SetHStyle( fit, 0 );
	vFit[iG][iX].push_back( fit );

	TCanvas c( hProj->GetName(), hProj->GetName(), 800, 600 );
	
	hProj->SetTitle("");
	hProj->Draw();

	fit->Draw("same");
	  
	DrawAtlasRight();
	drawTool->DrawLeftLatex( 0.18, 0.74, Form("%s", jzn.c_str() ) );
    
	drawTool->DrawLeftLatex
	  ( 0.18, 0.88, anaTool->GetYstarLabel( xLow, xUp, 1 ) );
	drawTool->DrawLeftLatex
	  ( 0.18, 0.81, anaTool->GetLabel( yMin, yMax, "#it{p}_{T}^{Truth} [GeV]") );

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

  std::string xAxisTitle = vJznHin[0]->GetXaxis()->GetTitle();
  std::string yAxisTitle = vJznHin[0]->GetYaxis()->GetTitle();

  // make 2D histos with scale and resolution (mean, sigma)
  TH2D* hMeanAll = new TH2D
    ( Form("h_%s_%s", type.c_str(), sMean.c_str() ),
      Form(";%s;%s;%s",xAxisTitle.c_str(),
	   yAxisTitle.c_str(), yTitleMean.c_str() ),
      nBinsX, xMin, xMax,
      nBinsY, yMin, yMax );

  TH2D* hSigmaAll = new TH2D
    ( Form("h_%s_%s",
	   type.c_str(), sSigma.c_str() ),
      Form(";%s;%s;%s",
	   xAxisTitle.c_str(), yAxisTitle.c_str(),
	   yTitleSigma.c_str() ),
      nBinsX, xMin, xMax,
      nBinsY, yMin, yMax );
  
  for( int iX = 0; iX < nBinsX; iX++ ){
    int xBin = iX + 1;

    double xLow, xUp;
    anaTool->GetBinRange
      ( vJznHin[0]->GetXaxis(), xBin, xBin, xLow, xUp );
        
    std::vector< TH1* > vYstarMeanTemp;
    std::vector< TH1* > vYstarSigmaTemp;
    std::vector< TH1* > vYstarNentTemp;
    
    for( uint iG = 0; iG < m_nJzn; iG++ ){
      vYstarMeanTemp .push_back( vMeans [iG][iX] );
      vYstarSigmaTemp.push_back( vSigmas[iG][iX] );
      vYstarNentTemp .push_back( vNent  [iG][iX] );
    }

    // Make the final, combined mean, sigma
    // Note - 
    // we use yAxisTitle from 3D histo
    // as xAxis on the final plots.
    TH1* hMeanFinal = new TH1D
      ( Form("h_%s_%s_%s",
	     type.c_str(), sMean.c_str(),
	     anaTool->GetName(xLow, xUp, "Ystar").c_str() ),
	Form("%s;%s;%s",
	     anaTool->GetYstarLabel( xLow, xUp, 1 ).c_str(),
	     yAxisTitle.c_str(), yTitleMean.c_str() ),
	nBinsY, yMin, yMax );
    vMeansFinal.push_back( hMeanFinal );
    
    TH1* hSigmaFinal = new TH1D
      ( Form("h_%s_%s_%s",
	     type.c_str(), sSigma.c_str(),
	     anaTool->GetName( xLow, xUp, "Ystar" ).c_str() ),
	Form("%s;%s;%s",
	     anaTool->GetYstarLabel( xLow, xUp, 1 ).c_str(),
	     yAxisTitle.c_str(), yTitleSigma.c_str() ),
	nBinsY, yMin, yMax );
    vSigmasFinal.push_back( hSigmaFinal );

    CombineSamples( hMeanFinal , vYstarMeanTemp , vYstarNentTemp );
    CombineSamples( hSigmaFinal, vYstarSigmaTemp, vYstarNentTemp );

    for( int yBin = 1; yBin <= nBinsY; yBin++ ){
      hMeanAll ->SetBinContent( xBin, yBin, hMeanFinal ->GetBinContent(yBin) );
      hSigmaAll->SetBinContent( xBin, yBin, hSigmaFinal->GetBinContent(yBin) );
    }
  } // end loop over xBin

  DrawCanvas( vMeansFinal ,
	      Form( "h_%s", type.c_str() ), sMean );
  DrawCanvas( vSigmasFinal,
	      Form( "h_%s", type.c_str() ), sSigma );

  std::string fOutName =
    m_dirOut + "/" + type + "_" + m_labelOut + ".root";
  
  TFile* fOut = new TFile( fOutName.c_str(), "RECREATE" );
  
  hMeanAll->Write();
  hSigmaAll->Write();

  fOut->Close();
  
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

void DiJetAnalysisMC::MakeSpectCFactorsRespMat( std::vector< TH2* >& vHspectReco,
						std::vector< TH2* >& vHspectTruth,
					        std::vector< TH3* >& vHspectRespMat,
						const std::vector< std::string >& vLabels,
						const std::string& nameCFactors,
						const std::string& nameRespMat ){

  std::string axisLabel, axisLabelTex;
  GetSpectraLabels( axisLabel, axisLabelTex, nameRespMat );  

  std::vector< TH1* > vSpect;
  std::vector< TH1* > vCfactors;
  std::vector< TH2* > vRespMat;
  std::vector< TF1* > vCfits;

  // these are needed to later combine the samples
  std::vector< std::vector< TH1* > > vCfactorsJzn;
  std::vector< std::vector< TH1* > > vNentTJzn;
  std::vector< std::vector< TH1* > > vNentRJzn;
  std::vector< std::vector< TH2* > > vRespMatJzn;

  std::vector< TH1* > vCfactorsAll;
  
  vCfactorsJzn.resize( m_nJzn );
  vNentTJzn   .resize( m_nJzn );
  vNentRJzn   .resize( m_nJzn );
  vRespMatJzn .resize( m_nJzn );
  
  for( uint iG = 0; iG < vHspectRespMat.size(); iG++ ){

    std::string label = vLabels[iG];

    TH2* hReco    = vHspectReco   [ iG ];
    TH2* hTruth   = vHspectTruth  [ iG ];
    TH3* hRespMatFull = vHspectRespMat[ iG ];
    
    for( int xBin = 1; xBin <= hRespMatFull->GetNbinsX(); xBin++ ){

      double xMin, xMax;
      anaTool->GetBinRange
	( hRespMatFull->GetXaxis(), xBin, xBin, xMin, xMax );

      hRespMatFull->GetXaxis()->SetRange( xBin, xBin );

      std::string hTag = anaTool->GetName( xMin, xMax, axisLabel);
      
      TH2* hRespMat = static_cast< TH2D* >( hRespMatFull->Project3D("zy") );
      hRespMat->SetName( Form("h_%s_%s_%s", nameRespMat.c_str(), label.c_str(), hTag.c_str() ) ) ;
      styleTool->SetHStyle( hRespMat, 0 );
      vRespMat.push_back( hRespMat );

      TCanvas cSpectRespMat( hRespMat->GetName(), hRespMat->GetName(), 800, 600 );
      cSpectRespMat.SetLogz();
      
      hRespMat->Draw("colz");
      hRespMat->SetTitle("");

      drawTool->DrawLeftLatex
	( 0.2, 0.85, anaTool->GetLabel( xMin, xMax, axisLabelTex ).c_str() );
      
      DrawAtlasRight();

      hRespMat->Write();

      if( !label.compare( m_allName ) ){
	SaveAsAll( cSpectRespMat, hRespMat->GetName() );
      }

      TH1* hR = static_cast< TH1D* >
	( hReco->ProjectionY( Form("h_R_%s_%s_%s_%s",
				   m_ystarSpectName.c_str(), m_sCounts.c_str(),
				   label.c_str(), hTag.c_str() ) , xBin, xBin ) );

      TH1* hT = static_cast< TH1D* >
	( hTruth->ProjectionY( Form("h_T_%s_%s_%s_%s",
				    m_ystarSpectName.c_str(), m_sCounts.c_str(),
				    label.c_str(), hTag.c_str() ) , xBin, xBin ) );

      TH1* hC = static_cast< TH1D* >
	  ( hT->Clone( Form("h_%s_%s_%s",
			    nameCFactors.c_str(), label.c_str(), hTag.c_str() ) ) );
      styleTool->SetHStyleRatio( hC );	  
      hC->Divide( hR );
      vCfactors.push_back( hC );

      hC->SetYTitle( "Correction Factor" );

      SetCfactorsErrors( hR, hT, hRespMat, hC );
      
      // Do the fitting and drawing only for
      // the jzn individual samples. Do this later
      // for the combined sample.
      if( label.compare( m_allName ) ){
	      
	TF1* cFit = anaTool->FitPol2
	  ( hC, hC->GetXaxis()->GetXmin(), hC->GetXaxis()->GetXmax() );
	vCfits.push_back( cFit );
	    
	cFit->SetLineStyle( 2 );
	cFit->SetLineColor( kBlue );

	TCanvas cCFactors( hC->GetName(), hC->GetName(), 800, 600 );
	cCFactors.SetLogz();
	    
	hC->Draw("ep");
	cFit->Draw("same");

	drawTool->DrawLeftLatex
	  ( 0.2, 0.85, anaTool->GetLabel( xMin, xMax, axisLabelTex ).c_str() );
      	    
	DrawAtlasRight();

	hC->Write();
	cFit->Write();

	// keep track of only JZ slice counts and response matrixes.
	// these vectors are the ones that will be used to make the
	// the final correction factors.
	// For CFactor histograms, only keep track of ones with
	// "All". These factors will be changed on CombineSamples.
	vCfactorsJzn[ iG ].push_back( hC );
	vNentTJzn   [ iG ].push_back( hT );
	vNentRJzn   [ iG ].push_back( hR );
	vRespMatJzn [ iG ].push_back( hRespMat );

	// SaveAsAll( cCFactors, hC->GetName() );
      } else {
	vCfactorsAll.push_back( hC );
      }
    }
  }

  // now combine to make the final samples
  for( uint spectBin = 0; spectBin < vCfactorsAll.size(); spectBin++ ){

    // prepare all jzn samples for each ystar1,2 pt1,2 bin
    std::vector< TH1* > vCfactorsTmp;
    std::vector< TH1* > vNentTTmp;
    std::vector< TH1* > vNentRTmp;
    std::vector< TH2* > vRespMatTmp;

    for( uint iG = 0; iG < m_nJzn; iG++ ){
      vCfactorsTmp.push_back( vCfactorsJzn[ iG ][ spectBin ] );
      vNentTTmp   .push_back( vNentTJzn   [ iG ][ spectBin ] );
      vNentRTmp   .push_back( vNentRJzn   [ iG ][ spectBin ] );
      vRespMatTmp .push_back( vRespMatJzn [ iG ][ spectBin ] );
    }

    TH1* hCAll = vCfactorsAll[ spectBin ]; 
    CombineSamples( hCAll, vCfactorsTmp, vNentTTmp, vNentRTmp, vRespMatTmp );

    hCAll->Write();
  }
  

  
  for( auto& h : vSpect    ){ delete h; }
  for( auto& h : vCfactors ){ delete h; }
  for( auto& h : vRespMat  ){ delete h; }
}

void DiJetAnalysisMC::MakeDphiCFactorsRespMat( std::vector< THnSparse* >& vHnT,
					       std::vector< THnSparse* >& vHnR,
					       std::vector< THnSparse* >& vHnDphiRespMat,
					       const std::vector< std::string >& vLabels,
					       const std::string& nameCFactors,
					       const std::string& nameRespMat ){

  // these are just for garbage collection
  std::vector< TH1* > vProj;
  std::vector< TH1* > vCfactors;
  std::vector< TH2* > vRespMat;
  std::vector< TF1* > vCfits;
  std::vector< TH1* > vCfitsQuality;
  
  // these are needed to later combine the samples
  std::vector< std::vector< TH1* > > vCfactorsJzn;
  std::vector< std::vector< TH1* > > vNentTJzn;
  std::vector< std::vector< TH1* > > vNentRJzn;
  std::vector< std::vector< TH2* > > vRespMatJzn;

  std::vector< TH1* > vCfactorsAll;
  
  vCfactorsJzn.resize( m_nJzn );
  vNentTJzn   .resize( m_nJzn );
  vNentRJzn   .resize( m_nJzn );
  vRespMatJzn .resize( m_nJzn );
  
  // ---- loop over group  ----
  // ---- (jzn or trigger) ----
  for( uint iG = 0; iG < vHnT.size(); iG++ ){      
    std::string label = vLabels[iG];

    TH1D* hCfactorsProb =
      new TH1D( Form( "hCfactorsProb_%s", label.c_str() ), ";Fit Probability", 25, 0, 1 );
    vCfitsQuality.push_back( hCfactorsProb );
	
    TH1D* hCfactorsChi2 =
      new TH1D( Form( "hCfactorsChi2_%s", label.c_str() ), ";#Chi^{2}/NDF", 25, 0, 10 );
    vCfitsQuality.push_back( hCfactorsChi2 );

    TH2D* hMiiTiRi = new TH2D( Form( "hMiiTiRi_%s", label.c_str() ),
			       ";M_{ii}/R_{i};M_{ii}/T_{i}",
			       50, 0, 1,
			       50, 0, 1 );
    vRespMat.push_back( hMiiTiRi );

    THnSparse* hnT    = vHnT[iG];
    THnSparse* hnR    = vHnR[iG];
    THnSparse* hnDphi = vHnDphiRespMat[iG];
    
    // just in case, reset the histogram's ranges
    // in case they were set in previous functions
    anaTool->ResetAxisRanges( hnT    );
    anaTool->ResetAxisRanges( hnR    );
    anaTool->ResetAxisRanges( hnDphi );
    
    TAxis* axisT0 = hnT->GetAxis( m_dPP->GetAxisI(0) ); int nAxis0Bins = axisT0->GetNbins();
    TAxis* axisT1 = hnT->GetAxis( m_dPP->GetAxisI(1) ); int nAxis1Bins = axisT1->GetNbins();
    TAxis* axisT2 = hnT->GetAxis( m_dPP->GetAxisI(2) ); int nAxis2Bins = axisT2->GetNbins();
    TAxis* axisT3 = hnT->GetAxis( m_dPP->GetAxisI(3) ); int nAxis3Bins = axisT3->GetNbins();
    
    TAxis* axisR0 = hnR->GetAxis( m_dPP->GetAxisI(0) ); 
    TAxis* axisR1 = hnR->GetAxis( m_dPP->GetAxisI(1) ); 
    TAxis* axisR2 = hnR->GetAxis( m_dPP->GetAxisI(2) ); 
    TAxis* axisR3 = hnR->GetAxis( m_dPP->GetAxisI(3) );

    for( int axis0Bin = 1; axis0Bin <= nAxis0Bins; axis0Bin++ ){ 
      // check we are in correct ystar and pt bins
      if( !m_dPP->CorrectPhaseSpace
	  ( std::vector<int>{ axis0Bin, 0, 0, 0 } ) )
	{ continue; }

      axisT0->SetRange( axis0Bin, axis0Bin );
      axisR0->SetRange( axis0Bin, axis0Bin );
      double axis0Low, axis0Up;
      anaTool->GetBinRange( axisT0, axis0Bin, axis0Bin, axis0Low, axis0Up );
	   
      for( int axis1Bin = 1; axis1Bin <= nAxis1Bins; axis1Bin++ ){	
	axisT1->SetRange( axis1Bin, axis1Bin ); 
	axisR1->SetRange( axis1Bin, axis1Bin );
	double axis1Low, axis1Up;
	anaTool->GetBinRange( axisT1, axis1Bin, axis1Bin, axis1Low, axis1Up );

	for( int axis2Bin = 1; axis2Bin <= nAxis2Bins; axis2Bin++ ){
	  axisT2->SetRange( axis2Bin, axis2Bin ); 
	  axisR2->SetRange( axis2Bin, axis2Bin );
	  double axis2Low, axis2Up;
	  anaTool->GetBinRange( axisT2, axis2Bin, axis2Bin, axis2Low, axis2Up );

	  for( int axis3Bin = 1; axis3Bin <= nAxis3Bins; axis3Bin++ ){
	    // check we are in correct ystar and pt bins
	    if( !m_dPP->CorrectPhaseSpace
		( std::vector<int>{ axis0Bin, axis1Bin, axis2Bin, axis3Bin } ) )
	      { continue; }
	    
	    axisT3->SetRange( axis3Bin, axis3Bin ); 
	    axisR3->SetRange( axis3Bin, axis3Bin ); 	    
	    double axis3Low, axis3Up;
	    anaTool->GetBinRange( axisT3, axis3Bin, axis3Bin, axis3Low, axis3Up );
	    
	    std::string hTag =
	      Form( "%s_%s_%s_%s",
		    anaTool->GetName( axis0Low, axis0Up, m_dPP->GetAxisName(0) ).c_str(),
		    anaTool->GetName( axis1Low, axis1Up, m_dPP->GetAxisName(1) ).c_str(),
		    anaTool->GetName( axis2Low, axis2Up, m_dPP->GetAxisName(2) ).c_str(),
		    anaTool->GetName( axis3Low, axis3Up, m_dPP->GetAxisName(3) ).c_str() ); 

	    //---------------------------------------------------
	    //                 Response Matrixes
	    //---------------------------------------------------
	    // -------------- make dphi resp matrix -------------
	    // Here the axis are in order. ystar1 = axis0, etc.
	    // we need to get the mapped axis for this to be true.
	    // because axis0Bin as of now could be a pT axis.
	    // The mapped is back to original. Once we get that
	    // vector, the [0] is ystar1, [1] is ystar2, etc.
	    // Before, axis0, axis1 are arbitrary (chosen from m_dPP)
	    std::vector< int > vMappedaAxisBins = m_dPP->GetMappedBins
	      ( std::vector<int>{ axis0Bin, axis1Bin, axis2Bin, axis3Bin } );

	    int ystar1Bin =  vMappedaAxisBins[0];
	    int ystar2Bin =  vMappedaAxisBins[1];
	    // because of underflow bin we add 1 to whatever 
  	    // bin it actually is first actual bin is 2 not 1.
	    int pt1Bin    =  vMappedaAxisBins[2] + 1; 
	    int pt2Bin    =  vMappedaAxisBins[3] + 1;

	    // set the ranges, project onto dPhi x dPhi axis.
	    // for 
	    hnDphi->GetAxis(0)->SetRange( ystar1Bin, ystar1Bin  );
	    hnDphi->GetAxis(1)->SetRange( ystar2Bin, ystar2Bin  );
	    hnDphi->GetAxis(2)->SetRange( pt1Bin   , pt1Bin     );
	    hnDphi->GetAxis(3)->SetRange( pt1Bin   , pt1Bin     );
	    hnDphi->GetAxis(4)->SetRange( pt2Bin   , pt2Bin     );
	    hnDphi->GetAxis(5)->SetRange( pt2Bin   , pt2Bin     );

	    uint dPhiProjDimX = hnDphi->GetNdimensions() - 2;
	    uint dPhiProjDimY = hnDphi->GetNdimensions() - 1;
	    
	    TH2* hDphiRespMat = static_cast<TH2D*>
	      ( hnDphi->Projection( dPhiProjDimY, dPhiProjDimX ) );
	    hDphiRespMat->SetName
	      ( Form( "h_%s_%s_%s", nameRespMat.c_str(), label.c_str(), hTag.c_str() ) );
	    styleTool->SetHStyle( hDphiRespMat, 0 );
	    vRespMat.push_back( hDphiRespMat );
	    
	    TCanvas cDphiRespMat( hDphiRespMat->GetName(), hDphiRespMat->GetName(), 800, 600 );
	    cDphiRespMat.SetLogz();

	    hDphiRespMat->Draw("colz");
	    hDphiRespMat->SetTitle("");

	    hDphiRespMat->GetXaxis()->SetRange
	      ( m_dPhiRebinnedZoomLowBin, m_dPhiRebinnedZoomHighBin );
	    hDphiRespMat->GetYaxis()->SetRange
	      ( m_dPhiRebinnedZoomLowBin, m_dPhiRebinnedZoomHighBin );
	    
	    DrawTopLeftLabels
	      ( m_dPP, axis0Low, axis0Up, axis1Low, axis1Up,
		axis2Low, axis2Up, axis3Low, axis3Up );
	    
	    DrawAtlasRight();

	    hDphiRespMat->Write();

	    if( !label.compare( m_allName ) )
	      { SaveAsROOT( cDphiRespMat, hDphiRespMat->GetName() ); }
	
	    //---------------------------------------------------
	    //   Counts, Rebinned, Normalized dPhi Histograms
	    //---------------------------------------------------
	    // Take projection onto the dPhi axis
	    TH1* hT = hnT->Projection( 4 );
	    TH1* hR = hnR->Projection( 4 );

	    hT->SetName( Form( "h_T_%s_%s_%s_%s",
			       m_dPhiName.c_str(), m_sCounts.c_str(), label.c_str(), hTag.c_str() ) );
	    hR->SetName( Form( "h_R_%s_%s_%s_%s",
			       m_dPhiName.c_str(), m_sCounts.c_str(), label.c_str(), hTag.c_str() ) );

	    // write the unnormalized histograms.
	    vProj.push_back( hT );
	    vProj.push_back( hR );

	    hT->Write();
	    hR->Write();

	    // rebin before dividing and normalizing.
	    TH1* hTreb = static_cast< TH1D* >
	      ( hT->Rebin
		( m_nVarDphiRebinnedBins,
		  Form( "h_T_%s_%s_%s_%s",
			m_dPhiName.c_str(), m_sReb.c_str(), label.c_str(), hTag.c_str() ),
		  &m_varDphiRebinnedBinning[0] ) );   
	    TH1* hRreb = static_cast< TH1D* >
	      ( hR->Rebin
		( m_nVarDphiRebinnedBins,
		  Form( "h_R_%s_%s_%s_%s",
			m_dPhiName.c_str(), m_sReb.c_str(), label.c_str(), hTag.c_str() ),
		  &m_varDphiRebinnedBinning[0] ) );
   
	    vProj.push_back( hTreb );
	    vProj.push_back( hRreb );
   
	    hTreb->Write();
	    hRreb->Write();
	    
	    //---------------------------------------------------
	    //               Correction Factors
	    //---------------------------------------------------
	    
	    // Make the ratio histogram. This is used to unfold bin-by-bin.
	    // later these are combined on a per-jz sample basis.
	    TH1* hC = static_cast< TH1D* >
	      ( hTreb->Clone
		( Form( "h_%s_%s_%s", nameCFactors.c_str(), label.c_str(), hTag.c_str())));
	    styleTool->SetHStyleRatio( hC );	  
	    hC->Divide( hRreb );
	    vCfactors.push_back( hC );

	    hC->SetYTitle( "Correction Factor" );

	    SetCfactorsErrors( hRreb, hTreb, hDphiRespMat, hC );
	    	    
	    // Smooth the cfactors
	    hC->GetXaxis()->SetRange
	      ( m_dPhiRebinnedZoomLowBin, m_dPhiRebinnedZoomHighBin );
	    // hC->Smooth( 1, "R" );

	    // Do the fitting and drawing only for
	    // the jzn individual samples. Do this later
	    // for the combined sample.
	    if( label.compare( m_allName ) ){
	      
	      TF1* cFit = anaTool->FitPol2( hC, 2.2, m_dPhiZoomHigh );
	      vCfits.push_back( cFit );
	    
	      cFit->SetLineStyle( 2 );
	      cFit->SetLineColor( kBlue );

	      if( cFit->GetProb() ){
		hCfactorsChi2->Fill( cFit->GetChisquare()/cFit->GetNDF() );
		hCfactorsProb->Fill( cFit->GetProb() );
	      }
	      
	      TCanvas cCFactors( hC->GetName(), hC->GetName(), 800, 600 );
	      cCFactors.SetLogz();
	    
	      hC->Draw("ep");
	      cFit->Draw("same");
	    
	      DrawTopLeftLabels
		( m_dPP, axis0Low, axis0Up, axis1Low, axis1Up,
		  axis2Low, axis2Up, axis3Low, axis3Up );
	    
	      DrawAtlasRight();

	      hC->Write();
	      cFit->Write();

	      // keep track of only JZ slice counts and response matrixes.
	      // these vectors are the ones that will be used to make the
	      // the final correction factors.
	      // For CFactor histograms, only keep track of ones with
	      // "All". These factors will be changed on CombineSamples.
	      vCfactorsJzn[ iG ].push_back( hC );
	      vNentTJzn   [ iG ].push_back( hTreb );
	      vNentRJzn   [ iG ].push_back( hRreb );
	      vRespMatJzn [ iG ].push_back( hDphiRespMat );

	      // SaveAsROOT( cCFactors, hC->GetName() );
	    } else {
	      vCfactorsAll.push_back( hC );
	    }
	  } // end loop over axis4
	} // end loop over axis2
      } // end loop over axis1
    } // end loop over axis0
    
    hCfactorsProb->Write();
    hCfactorsChi2->Write();
    hMiiTiRi     ->Write();
  } // end loop over jzn sample

  // now combine to make the final samples
  for( uint ystarPtBin = 0; ystarPtBin < vCfactorsAll.size(); ystarPtBin++ ){

    // prepare all jzn samples for each ystar1,2 pt1,2 bin
    std::vector< TH1* > vCfactorsTmp;
    std::vector< TH1* > vNentTTmp;
    std::vector< TH1* > vNentRTmp;
    std::vector< TH2* > vRespMatTmp;

    for( uint iG = 0; iG < m_nJzn; iG++ ){
      vCfactorsTmp.push_back( vCfactorsJzn[ iG ][ ystarPtBin ] );
      vNentTTmp   .push_back( vNentTJzn   [ iG ][ ystarPtBin ] );
      vNentRTmp   .push_back( vNentRJzn   [ iG ][ ystarPtBin ] );
      vRespMatTmp .push_back( vRespMatJzn [ iG ][ ystarPtBin ] );
    }

    TH1* hCAll = vCfactorsAll[ ystarPtBin ]; 
    CombineSamples( hCAll, vCfactorsTmp, vNentTTmp, vNentRTmp, vRespMatTmp );

    hCAll->Write();
  }
  
  for( auto& f : vCfits        ){ delete f; }
  for( auto& h : vProj         ){ delete h; }
  for( auto& f : vCfactors     ){ delete f; }
  for( auto& r : vRespMat      ){ delete r; }
  for( auto& h : vCfitsQuality ){ delete h; }
}


void DiJetAnalysisMC::MakePtRespMat( std::vector< THnSparse* >& vhnPt,
				     const std::vector< std::string >& vLabels,
				     const std::string& namePt ){  
  
  // ----------- pT response matrix -------------
  std::vector< TH2* > vRespMat;
  std::vector< TH1* > vPurityEff;
 
  // ---- loop over group  ----
  // ---- (jzn or trigger) ----
  for( uint iG = 0; iG < vhnPt.size(); iG++ ){      
    // in data only draw for all
    std::string label = vLabels[iG];

    if( label.compare( m_allName ) ){ continue; } 

    THnSparse* hnPt   = vhnPt  [iG];

    // just in case, reset the axis ranges if this
    // histogram was used in previous functions
    anaTool->ResetAxisRanges( hnPt );
    
    // want the ystar axes (0, 1);
    // dont need to do anything fancy like
    // get mapped axis.
    TAxis* axis0yStar = hnPt->GetAxis(0);
    TAxis* axis1yStar = hnPt->GetAxis(1);

    TAxis* dPhiRaxis  = hnPt->GetAxis( 6 );
    TAxis* dPhiTaxis  = hnPt->GetAxis( 7 );

    int nAxis0BinsYstar = axis0yStar->GetNbins();
    int nAxis1BinsYstar = axis1yStar->GetNbins();

    int nDphiRbins      = dPhiRaxis->GetNbins();
    int nDphiTbins      = dPhiTaxis->GetNbins();

    // -------------- do resp mat for Pt    
    for( int axis0Bin = 1; axis0Bin <= nAxis0BinsYstar; axis0Bin++ ){ 
      // check we are in correct ystar and pt bins
      if( !m_dPP->CorrectPhaseSpace
	  ( std::vector<int>{ axis0Bin, 0, 0, 0 } ) )
	{ continue; }

      axis0yStar->SetRange( axis0Bin, axis0Bin );
      double axis0Low, axis0Up;
      anaTool->GetBinRange
	( axis0yStar, axis0Bin, axis0Bin, axis0Low, axis0Up );

      for( int axis1Bin = 1; axis1Bin <= nAxis1BinsYstar; axis1Bin++ ){

	axis1yStar->SetRange( axis1Bin, axis1Bin );
	double axis1Low, axis1Up;
	anaTool->GetBinRange
	  ( axis1yStar, axis1Bin, axis1Bin, axis1Low, axis1Up );
	
	// Make pT resp matrix
	std::string hTag =
	  Form( "%s_%s",
		anaTool->GetName( axis0Low, axis0Up, m_dPP->GetDefaultAxisName(0) ).c_str(),
		anaTool->GetName( axis1Low, axis1Up, m_dPP->GetDefaultAxisName(1) ).c_str() ); 

	std::vector< std::string > vPtLabel{ m_s_pt1, m_s_pt2 };
	std::vector< int > vPtProjDimY{ 3, 5 };
	std::vector< int > vPtProjDimX{ 2, 4 };

	// Reset the dphi bin ranges so we can still get
	// total pT migration. Do not call ResetAxisRanges
	// because we are in some range in ystar atm.
	dPhiRaxis->SetRange( 1, nDphiRbins );
	dPhiTaxis->SetRange( 1, nDphiTbins );
	
	// do the same for pt1 and pt2
	for( uint ptN = 0; ptN < vPtLabel.size(); ptN++ ){
	  // Take projection onto the two pT axis
	  TH2* hPtRespMat = static_cast< TH2D* >
	    ( hnPt->Projection( vPtProjDimY[ ptN ], vPtProjDimX[ ptN ] ) );
	  hPtRespMat->SetName
	    ( Form( "h_%s_%s_%s_%s", namePt.c_str(), vPtLabel[ ptN ].c_str(),
		    label.c_str(), hTag.c_str() ) );
	  styleTool->SetHStyle( hPtRespMat, 0 );
	  vRespMat.push_back( hPtRespMat );

	  // kind of a pointless bin
	  hPtRespMat->SetBinContent( hPtRespMat->FindBin( 20., 20. ), 0 );
     	    
	  TCanvas c1( hPtRespMat->GetName(), hPtRespMat->GetName(), 800, 600 );
	  c1.SetLogz();

	  gStyle->SetPaintTextFormat(".3f");
	  anaTool->AverageOver( hPtRespMat, "row" );
	  hPtRespMat->Draw("col text");
	  hPtRespMat->SetTitle("");

	  drawTool->DrawLeftLatex
	    ( 0.13, 0.86, CT::AnalysisTools::GetLabel
	      ( axis0Low, axis0Up, m_dPP->GetDefaultAxisLabel(0) ) );
	  drawTool->DrawLeftLatex
	    ( 0.13, 0.79, CT::AnalysisTools::GetLabel
	      ( axis1Low, axis1Up, m_dPP->GetDefaultAxisLabel(1) ) );  
	
	  DrawAtlasRight();
	
	  SaveAsAll( c1, hPtRespMat->GetName() );
	  hPtRespMat->Write();

	  TH1* hPtPurity = static_cast<TH1D*>
	    (hPtRespMat->ProjectionX
	     ( Form( "h_%s_%s_%s_%s_%s", namePt.c_str(), vPtLabel[ptN].c_str(),
		     m_purityName.c_str(), label.c_str(), hTag.c_str())) );
	  TH1* hPtEff    = static_cast<TH1D*>
	    ( hPtRespMat->ProjectionY
	      ( Form( "h_%s_%s_%s_%s_%s", namePt.c_str(), vPtLabel[ptN].c_str(),
		      m_effName.c_str(), label.c_str(), hTag.c_str())) );
	  vPurityEff.push_back( hPtPurity );
	  vPurityEff.push_back( hPtEff    );
	  
	  hPtPurity->Reset();
	  hPtEff->Reset();

	  styleTool->SetHStyle( hPtPurity, 0 );
	  styleTool->SetHStyle( hPtEff   , 1 );
	  
	  GetPurityEff( hPtRespMat, hPtPurity, hPtEff );
	  
	  TCanvas c2( hPtPurity->GetName(), hPtPurity->GetName(), 800, 600 );
	  TLegend leg( 0.60, 0.22, 0.90, 0.33 );
	  styleTool->SetLegendStyle( &leg );
	  leg.SetFillStyle(0);

	  leg.AddEntry( hPtPurity,  "Purity" );
	  leg.AddEntry( hPtEff, "Efficiency" );

	  hPtPurity->SetMaximum(1.0);
	  hPtPurity->SetMinimum(0.0);
	  hPtEff   ->SetMaximum(1.0);
	  hPtEff   ->SetMinimum(0.0);
	  
	  hPtPurity->Draw("ep same");
	  hPtEff   ->Draw("ep same");
	  leg.Draw();

	  drawTool->DrawLeftLatex
	    ( 0.13, 0.86, CT::AnalysisTools::GetLabel
	      ( axis0Low, axis0Up, m_dPP->GetDefaultAxisLabel(0) ) );
	  drawTool->DrawLeftLatex
	    ( 0.13, 0.79, CT::AnalysisTools::GetLabel
	      ( axis1Low, axis1Up, m_dPP->GetDefaultAxisLabel(1) ) );  
	  
	  DrawAtlasRight();
	
	  SaveAsROOT( c2, hPtPurity->GetName() );

	  // Now loop over delta Phi bins, and make these pt projections
	  // over some combination (top 4x4). Dont forget to reset
	  // the delta Phi ranges after this!

	  // keep this same for both of them for now.
	  // assume the have > 4 bins.
	  int dPhiBinRange = dPhiRaxis->GetNbins() > 4 ? 4 : dPhiRaxis->GetNbins();
	  int dPhiBinBegin = nDphiRbins - dPhiBinRange + 1;
	  for( int dPhiRbin = dPhiBinBegin; dPhiRbin <= nDphiRbins; dPhiRbin++ ){
	    dPhiRaxis->SetRange( dPhiRbin, dPhiRbin );
	    	  
	    for( int dPhiTbin = dPhiBinBegin; dPhiTbin <= nDphiTbins; dPhiTbin++ ){
	      dPhiTaxis->SetRange( dPhiTbin, dPhiTbin );
	      
	      for( uint ptN = 0; ptN < vPtLabel.size(); ptN++ ){
		// Make pT dPhi resp matrix
		std::string hTagDphi =
		  Form( "%s_dPhiR_%d_dPhiT_%d", hTag.c_str(), dPhiRbin, dPhiTbin ); 
		
		TH2* hPtRespMatDphi = static_cast< TH2D* >
		  ( hnPt->Projection( vPtProjDimY[ ptN ], vPtProjDimX[ ptN ] ) );	      	      
		hPtRespMatDphi->SetName
		  ( Form( "h_%s_%s_%s_%s_%s", namePt.c_str(), m_dPhiName.c_str(),
			  vPtLabel[ ptN ].c_str(), label.c_str(), hTagDphi.c_str() ) );
		styleTool->SetHStyle( hPtRespMatDphi, 0 );
		vRespMat.push_back( hPtRespMatDphi );

		// kind of a pointless bin
		hPtRespMatDphi->SetBinContent( hPtRespMatDphi->FindBin( 20., 20. ), 0 );
	      
		TCanvas c1( hPtRespMatDphi->GetName(), hPtRespMatDphi->GetName(), 800, 600 );
		c1.SetLogz();

		gStyle->SetPaintTextFormat(".3f");
		anaTool->AverageOver( hPtRespMatDphi, "row" );
		hPtRespMatDphi->Draw("col text");
		hPtRespMatDphi->SetTitle("");

		drawTool->DrawLeftLatex
		  ( 0.13, 0.86, CT::AnalysisTools::GetLabel
		    ( axis0Low, axis0Up, m_dPP->GetDefaultAxisLabel(0) ) );
		drawTool->DrawLeftLatex
		  ( 0.13, 0.79, CT::AnalysisTools::GetLabel
		    ( axis1Low, axis1Up, m_dPP->GetDefaultAxisLabel(1) ) );  
		drawTool->DrawLeftLatex
		  ( 0.4, 0.86, Form( "#delta#phi_{Reco} Bin %d" , dPhiRbin ) );
		drawTool->DrawLeftLatex
		  ( 0.4, 0.79, Form( "#delta#phi_{Truth} Bin %d", dPhiTbin ) );  
		DrawAtlasRight();
	
		SaveAsAll( c1, hPtRespMatDphi->GetName() );
		hPtRespMatDphi->Write();
	      } // end loop over ptN
	    } // end loop over dPhiTbin
	  } // end loop over dPhiRbin
	  
	} // end loop over ptN
      } // end loop over ystar2
    } // end loop over ystar1
  } // end loop over iG
  for( auto rm : vRespMat     ){ delete rm; }
  for( auto pe : vPurityEff   ){ delete pe; }
}

//---------------------------
//        Drawing
//---------------------------
void DiJetAnalysisMC::DrawCanvas( std::vector< TH1* >& vHIN,
				  const std::string& type1,
				  const std::string& type2 ){
  TCanvas c("c","c",800,600);
  
  TLegend leg(0.60, 0.61, 0.99, 0.82);
  styleTool->SetLegendStyle( &leg );
  leg.SetFillStyle(0);

  int style = 0;

  // for situations where dont want to
  // plot every single bin 
  // plot every n on canvas
  for( uint xRange = 0; xRange < vHIN.size(); xRange++ ){
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
  
  SaveAsAll( c, Form("%s_%s", type1.c_str(), type2.c_str() ) );
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
  } else if( type.find("recoTruthDeta") != std::string::npos ||
	     type.find("recoTruthDphi") != std::string::npos ) { // ANGLES
    y0 = 0;
    y0 = 0;
  }

  return y0;
}

void DiJetAnalysisMC::DrawAtlasRight( double x0, double y0, double scale )
{ drawTool->DrawAtlasInternalMCRight( x0, y0, m_mcTypeLabel, m_is_pPb, scale ); } 


void DiJetAnalysisMC::DrawAtlasRightBoth( double x0, double y0, double scale )
{ drawTool->DrawAtlasInternalMCRight( x0, y0, m_mcTypeLabel, 3, scale); }
