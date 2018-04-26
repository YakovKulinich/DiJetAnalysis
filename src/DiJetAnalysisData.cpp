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
#include <THnSparse.h>
#include <TF1.h>
#include <TVirtualFitter.h>
#include <TGaxis.h>

#include <TLine.h>

#include <fstream>
#include <iostream>
#include <cmath>

#include <boost/assign.hpp>
#include <boost/format.hpp>

#include "DiJetAnalysisData.h"
#include "DeltaPhiProj.h"

DiJetAnalysisData::DiJetAnalysisData()
  : DiJetAnalysisData( false, 0, 0 ){}

DiJetAnalysisData::DiJetAnalysisData( bool is_pPb )
  : DiJetAnalysisData( is_pPb, 0, 0 ){}

DiJetAnalysisData::DiJetAnalysisData( bool is_pPb, int mcType )
  : DiJetAnalysisData( is_pPb, mcType, 0 ){}

DiJetAnalysisData::DiJetAnalysisData( bool is_pPb, int mcType, int uncertComp )
  : DiJetAnalysis( is_pPb, true, mcType, uncertComp ){
  
  //========== Set Histogram Binning =============

  //==================== Cuts ====================

  //================== Settings ==================
  m_mbTriggerName      = "";
  
  m_mbTriggerI         = -1;
  m_lowestCentTriggerI = -1;
  
  m_centMbCorrection   = 0;

  m_nJetsRunName = "nJetsRun";
}

DiJetAnalysisData::~DiJetAnalysisData(){}

void DiJetAnalysisData::Initialize(){
  // Initalize things common to everything
  DiJetAnalysis::Initialize();
  
  m_fNameIn = m_is_pPb ?
    "/home/yakov/Projects/atlas/data/pPb/pPb.root" :
    "/home/yakov/Projects/atlas/data/pp/pp.root"  ;

  // directory of where unfolding MC files are

  std::string system = m_is_pPb ? m_s_pPb : m_s_pp;

  std::string unfoldingMCdir =
    Form( "%s/%s_%s_%s_%s",
	  m_sOutput.c_str(), m_sOutput.c_str(), system.c_str(),
	  m_sMC.c_str(), m_mcTypeName.c_str());
  
  // name of file that is used as input for performance unfolding
  m_fNamePerfUnfoldingMC
    = Form( "%s/%s_%s_%s_%s_%s_%s.root",
	    unfoldingMCdir.c_str(), m_myOutName.c_str(),
	    system.c_str(), m_sMC.c_str(), m_mcTypeName.c_str(),
	    m_sPerf.c_str(), m_uncertSuffix.c_str() );

  // name of file that is used as input for physics unfolding
  m_fNamePhysUnfoldingMC
    = Form( "%s/%s_%s_%s_%s_%s_%s.root",
	    unfoldingMCdir.c_str(), m_myOutName.c_str(),
	    system.c_str(), m_sMC.c_str(), m_mcTypeName.c_str(),
	    m_sPhys.c_str(), m_uncertSuffix.c_str() );  
}

//---------------------------
// Fill Tree / Plot Controls
//---------------------------

void DiJetAnalysisData::RunOverTreeFillHistos( int nEvents, 
					       int startEvent ){  
  LoadTriggers();
  SetupHistograms();
  ProcessEvents( nEvents, startEvent );
  SaveOutputsFromTree();
}

void DiJetAnalysisData::ProcessPerformance(){

  LoadTriggers();
  LoadHistograms();

  TFile* fOut = new TFile( m_fNameDefPerf.c_str(),"RECREATE");

  // add a trigger "all" to collection
  // rest of plots include combined triggers
  m_vTriggers.push_back( m_allName );

  m_hAllEtaPtMap  = CombineSamples( m_vHtriggerEtaPtMap, "etaPtMap" );
  m_hAllEtaPhiMap = CombineSamples( m_vHtriggerEtaPhiMap, "etaPhiMap" );
  MakeEtaPhiPtMap( m_vHtriggerEtaPtMap , m_vTriggers, "etaPtMap" );
  MakeEtaPhiPtMap( m_vHtriggerEtaPhiMap , m_vTriggers, "etaPhiMap" );

  m_hAllYstarSpect = CombineSamples( m_vHtriggerYstarSpect, m_ystarSpectName );
  MakeSpectra( m_vHtriggerYstarSpect, m_vTriggers, m_ystarSpectName );

  m_hAllYstarSpectFine = CombineSamples( m_vHtriggerYstarSpectFine, m_ystarSpectFineName );
  MakeSpectra( m_vHtriggerYstarSpectFine, m_vTriggers, m_ystarSpectFineName );

  MakeEfficiencies( m_vHtriggerEtaSpectSim, m_vHtriggerEtaSpectDenom, m_etaEffName );

  // false at end => dont apply PS
  m_hAllNjetsRun = CombineSamples( m_vHtriggerNjetsRun, m_nJetsRunName, false );
  MakeNjetsRun( m_vHtriggerNjetsRun, m_vTriggers, m_nJetsRunName );

  std::cout << "DONE! Closing " << fOut->GetName() << std::endl;
  fOut->Close(); delete fOut;
  std::cout << "......Closed  " << std::endl;
}

void DiJetAnalysisData::UnfoldPerformance(){

  // for now
  // Copy File with original dPhi, spectra, etc,
  // into file where histos with corrections are
  // going to be appended. 
  std::cout << "Copy " << m_fNameDefPerf << " -> " << m_fNamePerfUF << std::endl;
  TFile::Cp( m_fNameDefPerf.c_str(), m_fNamePerfUF.c_str() );

  // Open two for reading one for updating.
  // open the MC file used for unfolding info.
  // open the data file used for measured info.
  // passed to unfolding function.
  TFile* fInMC   = TFile::Open( m_fNamePerfUnfoldingMC.c_str() );
  TFile* fInData = TFile::Open( m_fNameDefPerf.c_str()         );
  TFile* fOut    = new TFile( m_fNamePerfUF.c_str(),"UPDATE"   );

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
  TH2* m_hAllspectUnfolded =
    UnfoldSpectra( fInData, fInMC, m_ystarSpectUnfoldedName );
  m_vHspectUnfolded.push_back( m_hAllspectUnfolded );
  m_vLabelsUnfolded.push_back( m_allName );

  // unfold on MC, just used for testing purpose.
  // so there is no comb subt or normalization or scaling
  MakeSpectra( m_vHspectUnfolded, m_vLabelsUnfolded, m_ystarSpectUnfoldedName );

  std::cout << "DONE! Closing " << fOut->GetName() << std::endl;
  fOut->Close(); delete fOut;
  std::cout << "......Closed  " << std::endl;

}

void DiJetAnalysisData::ProcessPhysics(){

  LoadTriggers();
  LoadHistograms();

  // THIS HAS TO BE CHANGED
  TFile* fInPerf  = TFile::Open( m_fNameDefPerf.c_str() );

  TFile* fOut = new TFile( m_fNameDefPhys.c_str(),"RECREATE");

  // add a trigger "all" to collection
  // rest of plots include combined triggers
  m_vTriggers.push_back( m_allName );

  m_hAllDphi = CombineSamples( m_vHtriggerDphi, m_dPhiName );
  MakeDeltaPhi( m_vHtriggerDphi, m_vTriggers, m_dPhiName, fInPerf, m_ystarSpectName );
  
  std::cout << "DONE! Closing " << fOut->GetName() << std::endl;
  fOut->Close(); delete fOut;
  std::cout << "......Closed  " << std::endl;
}

void DiJetAnalysisData::UnfoldPhysics(){

  // for now
  // Copy File with original dPhi, spectra, etc,
  // into file where histos with corrections are
  // going to be appended. 
  std::cout << "Copy " << m_fNameDefPhys << " -> " << m_fNamePhysUF << std::endl;
  TFile::Cp( m_fNameDefPhys.c_str(), m_fNamePhysUF.c_str() );

  // Open two for reading one for updating.
  // open the MC file used for unfolding info.
  // open teh data file used for measured info.
  // passed to unfolding function.
  TFile* fInMC     = TFile::Open( m_fNamePhysUnfoldingMC.c_str() );
  TFile* fInData   = TFile::Open( m_fNameDefPhys.c_str()         );
  TFile* fInMCPerf = TFile::Open( m_fNamePerfUF.c_str()          );
  TFile* fOut      = new TFile( m_fNamePhysUF.c_str(), "UPDATE"  );

  std::cout << m_fNamePhysUnfoldingMC << std::endl;
  
  std::cout << "----- Unfolding Data ------" << std::endl;
  // make a vector with just the unfolded result.
  // this is to send it to MakeDeltaPhi(..) to have
  // unfolded results plotted separately
  std::vector< THnSparse*  > m_vHDphiUnfolded;
  std::vector< std::string > m_vLabelUnfolded;
  
  // make unfolded THnSparse with similar naming convention
  // as the other histograms. At this point, don't care about
  // doing this for all triggers. Altohugh, this can be
  // repeated in a loop with m_allName subsitituted for trigger,
  // and subsequently added to the vectors above.  
  THnSparse* m_hAllDphiUnfolded =
    UnfoldDeltaPhi( fInData, fInMC, m_dPhiUnfoldedName,
		    fInMCPerf, m_ystarSpectUnfoldedName );
  m_vHDphiUnfolded.push_back( m_hAllDphiUnfolded );
  m_vLabelUnfolded.push_back( m_allName );

  // make deltaPhi, give flag (true) that its unfolded response
  // so there is no comb subt or normalization or scaling
  MakeDeltaPhi( m_vHDphiUnfolded, m_vLabelUnfolded, m_dPhiUnfoldedName,
		fInMCPerf, m_ystarSpectUnfoldedName );
  
  std::cout << "DONE! Closing " << fOut->GetName() << std::endl;
  fOut->Close(); delete fOut;
  std::cout << "......Closed  " << std::endl;

  /*
  fOut = new TFile( m_fNamePhysUF.c_str(),"UPDATE");
  CompareCfactorsRBnRB( fOut );
  fOut->Close(); delete fOut;
  */
}

void DiJetAnalysisData::ProcessSystematics(){

  TFile* fOut  = new TFile( m_fNameSYS.c_str(), "recreate");
 
  MakeSystematicsGraphs( fOut, m_widthName );
  MakeSystematicsGraphs( fOut, m_yieldName );

  std::cout << "DONE! Closing " << fOut->GetName() << std::endl;
  fOut->Close(); delete fOut;
  std::cout << "......Closed  " << std::endl;
}

//---------------------------------
//            Fill Tree
//---------------------------------

void DiJetAnalysisData::LoadTriggers(){
  std::string triggerMenu = GetConfig()->GetValue("triggerMenu", " ");

  m_vTriggers =  anaTool->vectorise
    ( GetConfig()->GetValue( Form("triggers.%s",triggerMenu.c_str()), ""), " ");

  std::vector< std::string> vTriggersRef = anaTool->vectorise
    ( GetConfig()->GetValue( Form("triggersRef.%s",triggerMenu.c_str()), "" ), " ");
  
  // Get Index for reference trigger
  for( uint iR = 0 ; iR < vTriggersRef.size(); iR++ ){
      for( uint iG = 0 ; iG < vTriggersRef.size(); iG++ ){
	if( !m_vTriggers[iG].compare( vTriggersRef[iR] ) )
	  { m_vTriggersRefIndex.push_back(iG); }
      }
  }

  for( uint iG = 0 ; iG < m_vTriggersRefIndex.size(); iG++ ){
    int refTrigIndex = m_vTriggersRefIndex[iG];
    std::cout << " For: "         << m_vTriggers[iG]
	      << " Ref Trig is: " << m_vTriggers[refTrigIndex]
	      << std::endl;      
  }

  m_vTriggersPrescale =  anaTool->vectoriseD
    ( GetConfig()->GetValue
      ( Form("triggersPrescale.%s",triggerMenu.c_str()), "" ), " ");

  m_vTriggersTholdPt =  anaTool->vectoriseD
    ( GetConfig()->GetValue
      ( Form("triggersThold.%s",triggerMenu.c_str()), "" ), " ");

  m_vTriggersEffPtLow = anaTool->vectoriseD
    ( GetConfig()->GetValue
      ( Form("triggersEffPtLow.%s",triggerMenu.c_str()), "" ), " ");

  m_vTriggersEffPtHigh = anaTool->vectoriseD
    ( GetConfig()->GetValue
      ( Form("triggersEffPtHigh.%s",triggerMenu.c_str()), "" ), " ");

  m_vTriggersEtaMin = anaTool->vectoriseD
    ( GetConfig()->GetValue
      ( Form("triggersEtaMin.%s",triggerMenu.c_str()), "" ), " ");

  m_vTriggersEtaMax = anaTool->vectoriseD
    ( GetConfig()->GetValue
      ( Form("triggersEtaMax.%s",triggerMenu.c_str()), "" ), " ");

  // total number of real triggers.
  // later, we add All to vTriggers
  // only for reading and writing purposes
  m_nTriggers = m_vTriggers.size();
  
  // find MB trigger, and lowest pp trigger
  // if it exists. This is to correct
  // MB threshold in central region
  std::string lowestCentTrigger =
    GetConfig()->GetValue( Form("lowestCentTrigger.%s",triggerMenu.c_str()), "");
  
  for( uint iG = 0; iG < m_nTriggers; iG++ ){
    std::string trigger = m_vTriggers[iG];
    if( trigger.find("_mb_") != std::string::npos ){ 
      m_mbTriggerI    = iG;
      m_mbTriggerName = trigger;
      std::cout << "Found " << trigger << " at " << m_mbTriggerI << std::endl;
    } else if( !trigger.compare( lowestCentTrigger ) ){ 
      m_lowestCentTriggerI  = iG;
      std::cout << "Found " << trigger << " at "
		<< m_lowestCentTriggerI << std::endl;
    }
  }

  if( m_lowestCentTriggerI >= 0 && m_mbTriggerI >= 0){
    m_centMbCorrection =
      m_vTriggersEffPtLow[ m_lowestCentTriggerI ] -
      m_vTriggersEffPtHigh[ m_mbTriggerI ];
    std::cout << "m_centMbCorrection: " << m_centMbCorrection << std::endl;
  }

  std::vector< int >    vRuns = anaTool->vectoriseI
    ( GetConfig()->GetValue
      ( Form("runs.%s",triggerMenu.c_str()), "" ), " ");
  std::vector< double > vLumi = anaTool->vectoriseD
    ( GetConfig()->GetValue
      ( Form("lumi.%s",triggerMenu.c_str()), "" ), " ");
  m_nRuns = vRuns.size();

  for( uint r = 0; r < m_nRuns; r++ ){
    int     runNumber = vRuns[r];
    double luminosity = vLumi[r];
    m_mRunBin [ runNumber ] = r + 1;
    m_mRunLumi[ runNumber ] = luminosity;
  }
  
  for( const auto& m : m_mRunBin ){
    std::cout << m.first << " : " << m.second << std::endl;
  }
  std::cout << "---" << std::endl;
  for( const auto& m : m_mRunLumi ){
    std::cout << m.first << " : " << m.second << std::endl;
  }
}

void DiJetAnalysisData::SetupHistograms(){

  for( auto& trigger : m_vTriggers ){
    std::cout << "Making - " << trigger << " histograms " << std::endl;          
    
    // -------- maps ---------
    m_vHtriggerEtaPhiMap.push_back
      ( new TH2D( Form("h_etaPhiMap_%s", trigger.c_str() ), 
		   ";#eta;#phi",
		   m_nEtaMapBins, m_etaMapMin, m_etaMapMax,
		   m_nPhiMapBins, m_phiMapMin, m_phiMapMax ) ) ;
    AddHistogram( m_vHtriggerEtaPhiMap.back() );
     
    m_vHtriggerEtaPtMap.push_back
      ( new TH2D( Form("h_etaPtMap_%s", trigger.c_str() ), 
		   ";#eta;#it{p}_{T} [GeV]",
		   m_nEtaMapBins, m_etaMapMin, m_etaMapMax,
		   m_nPtMapBins , m_ptMapMin , m_ptMapMax) ) ;
    AddHistogram( m_vHtriggerEtaPtMap.back() );

    m_vHtriggerEtaPhiPtMap.push_back
      ( new TH3D( Form("h_etaPhiPtMap_%s", trigger.c_str() ), 
		  ";#eta;#it{p}_{T} [GeV]",
		  m_nEtaMapBins, m_etaMapMin, m_etaMapMax,
		  m_nPhiMapBins, m_phiMapMin, m_phiMapMax,
		  22, 28, 50 ) ) ;
    AddHistogram( m_vHtriggerEtaPhiPtMap.back() );
    
    // -------- spect --------
    m_vHtriggerYstarSpect.push_back
      ( new TH2D( Form("h_%s_%s", m_ystarSpectName.c_str(), trigger.c_str() ), 
		  ";#it{y}*;#it{p}_{T1} [GeV];",
		  m_nVarYstarBins, 0, 1,
		  m_nVarPtBinsUfOf, 0, 1 ) ) ;
    m_vHtriggerYstarSpect.back()->GetXaxis()->
      Set( m_nVarYstarBins, &( m_varYstarBinning[0] ) );
    m_vHtriggerYstarSpect.back()->GetYaxis()->
      Set( m_nVarPtBinsUfOf, &( m_varPtBinningUfOf[0] ) );
    AddHistogram( m_vHtriggerYstarSpect.back() );

    m_vHtriggerYstarSpectFine.push_back
      ( new TH2D( Form("h_%s_%s", m_ystarSpectFineName.c_str(), trigger.c_str() ), 
		  ";#it{y}*;#it{p}_{T} [GeV];",
		  m_nVarYstarBins, 0, 1,
		  m_nPtSpectBins,m_ptSpectMin, m_ptSpectMax ) ) ;
    m_vHtriggerYstarSpectFine.back()->GetXaxis()->
      Set( m_nVarYstarBins, &( m_varYstarBinning[0] ) );
    AddHistogram( m_vHtriggerYstarSpectFine.back() );

    m_vHtriggerNjetsRun.push_back
      ( new TH3D( Form("h_%s_%s", m_nJetsRunName.c_str(), trigger.c_str() ), 
		  ";#it{y}*;#it{p}_{T} [GeV];Run Number",
		  m_nVarYstarBins, 0, 1,
		  m_nVarPtBins, 0, 1,
		  m_nRuns, 0, m_nRuns ) ) ;
    m_vHtriggerNjetsRun.back()->GetXaxis()->
      Set( m_nVarYstarBins, &( m_varYstarBinning[0] ) );
    m_vHtriggerNjetsRun.back()->GetYaxis()->
      Set( m_nVarPtBins, &( m_varPtBinning[0] ) );
    AddHistogram( m_vHtriggerNjetsRun.back() );

    m_vHtriggerNjetsRunFine.push_back
      ( new TH3D( Form("h_%s_fine_%s", m_nJetsRunName.c_str(), trigger.c_str() ), 
		  ";#it{y}*;#it{p}_{T} [GeV];Run Number",
		  m_nVarYstarBins, 0, 1,
		  m_nPtSpectBins,m_ptSpectMin, m_ptSpectMax,
		  m_nRuns, 0, m_nRuns ) ) ;
    m_vHtriggerNjetsRunFine.back()->GetXaxis()->
      Set( m_nVarYstarBins, &( m_varYstarBinning[0] ) );
    AddHistogram( m_vHtriggerNjetsRunFine.back() );

    // ----- efficiencies ----
    m_vHtriggerEtaSpectSim.push_back
      ( new TH2D ( Form("h_%sSim_%s", m_etaSpectName.c_str(), trigger.c_str() ), 
		   ";#eta;#it{p}_{T} [GeV]",
		   m_nVarFwdEtaBins, 0, 1, 50, 0, 100 ) ) ;
    m_vHtriggerEtaSpectSim.back()->GetXaxis()->
      Set( m_nVarFwdEtaBins, &( m_varFwdEtaBinning[0] ) );
    AddHistogram( m_vHtriggerEtaSpectSim.back() );

    m_vHtriggerEtaSpectDenom.push_back
      ( new TH2D( Form("h_%sDenom_%s", m_etaSpectName.c_str(), trigger.c_str() ), 
		  ";#eta;#it{p}_{T} [GeV];",
		  m_nVarFwdEtaBins, 0, 1, 50, 0, 100 ) ) ;
    m_vHtriggerEtaSpectDenom.back()->GetXaxis()->
      Set( m_nVarFwdEtaBins, &( m_varFwdEtaBinning[0] ) );
    AddHistogram( m_vHtriggerEtaSpectDenom.back() );
    
    // -------- dPhi --------
    m_nDphiDim    = m_vNdPhiBins.size();
   
    THnSparse* hn =
      new THnSparseD( Form("h_%s_%s", m_dPhiName.c_str(), trigger.c_str() ), "",
		      m_nDphiDim, &m_vNdPhiBins[0], &m_vDphiMin[0], &m_vDphiMax[0] );
    hn->GetAxis(0)->Set( m_nVarYstarBins, &( m_varYstarBinning[0] ) );
    hn->GetAxis(1)->Set( m_nVarYstarBins, &( m_varYstarBinning[0] ) );
    hn->GetAxis(2)->Set( m_nVarPtBins   , &( m_varPtBinning[0]    ) );
    hn->GetAxis(3)->Set( m_nVarPtBins   , &( m_varPtBinning[0]    ) );
    hn->GetAxis(4)->Set( m_nVarDphiBins , &( m_varDphiBinning[0]  ) );
    
    m_vHtriggerDphi.push_back( hn );
    AddHistogram( hn );
  }
}

void DiJetAnalysisData::ProcessEvents( int nEvents, int startEvent ){  
  // collections and variables
  std::vector< TLorentzVector >    vTrig_jets;
  std::vector< TLorentzVector >* p_vTrig_jets = &vTrig_jets;

  std::vector< TLorentzVector >    vR_jets;
  std::vector< TLorentzVector >* p_vR_jets = &vR_jets;

  std::vector< bool >    v_isCleanJet;
  std::vector< bool >* p_v_isCleanJet = &v_isCleanJet;

  std::vector< bool  > vTriggerFired;

  // because vector bool doesnt return lvalue
  std::map< int, bool > mTriggerFired; 
  
  int runNumber = 0;
  int LBN       = 0;
  //----------------------------------------
  //  Open file and tree, Fill histograms
  //----------------------------------------

  int nthreads = 4;
  ROOT::EnableImplicitMT(nthreads);

  std::cout << "fNameIn: " << m_fNameIn << std::endl;
  
  TFile* fIn  = TFile::Open( m_fNameIn.c_str() );
  TTree* tree = static_cast<TTree*>( fIn->Get( "tree" ) );

  // Connect to tree
  tree->SetBranchAddress( "vTrig_jets"  , &p_vTrig_jets );
  tree->SetBranchAddress( "vR_C_jets"   , &p_vR_jets );
  tree->SetBranchAddress( "v_isCleanJet", &p_v_isCleanJet );

  vTriggerFired.resize   ( m_nTriggers );

  for( uint iG = 0; iG < m_nTriggers; iG++ ){
    tree->SetBranchAddress
      ( Form("passed_%s", m_vTriggers[iG].c_str() ), &mTriggerFired[iG] );
  }

  tree->SetBranchAddress( "runNumber", &runNumber );
  tree->SetBranchAddress( "LBN"      , &LBN );  

  std::vector< int > mTrigPassed; 
  std::vector< int > mTrigFailed;  
  mTrigPassed.resize( m_nTriggers );
  mTrigFailed.resize( m_nTriggers );

  std::vector< int > mTrigJetsForward;
  std::vector< int > mTrigJetsTotal;
  mTrigJetsForward.resize( m_nTriggers );
  mTrigJetsTotal.resize  ( m_nTriggers );
  
  // n events
  int nEventsTotal = tree->GetEntries();

  nEvents = nEvents > 0 ? nEvents : nEventsTotal;
   startEvent = startEvent < nEventsTotal ? startEvent : nEventsTotal - 1;
  
  int endEvent = startEvent + nEvents < nEventsTotal ?
					startEvent + nEvents : nEventsTotal;

  // -------- EVENT LOOP ---------
  for( m_ev = startEvent; m_ev < endEvent; m_ev++ ){
    tree->GetEntry( m_ev );

    // THESE ARE LB THAT WERENT EXCLUDED INITIALLY. OLD GRL FOR 2015 pp
    if( !m_is_pPb && runNumber == 286411 && LBN > 1129 && LBN < 1131 ){ continue; }
    if( !m_is_pPb && runNumber == 286367 && LBN > 75   && LBN < 78   ){ continue; }
    if( !m_is_pPb && runNumber == 286367 && LBN > 121                ){ continue; }
    if( !m_is_pPb && runNumber == 286364 && LBN > 618                ){ continue; }
        
    if( anaTool->DoPrint(m_ev) ) {
      std::cout << "\nEvent : " << m_ev << "    runN : " << runNumber
		<< "    has : " << vR_jets.size() << " jets" 
		<< "    and : " << vTrig_jets.size() << " trig jets" 
		<< std::endl; 
    }

    ApplyCleaning ( vR_jets, v_isCleanJet );
    ApplyIsolation( vR_jets, 1.0 );
    
    std::sort( vR_jets.begin(), vR_jets.end(), anaTool->sortByDecendingPt );
    
    // some runs and lbn are bad
    // for efficiency plots
    bool goodRunLBN = true;
    // for pPb
    // skip that run for efficiencies and spectra
    // had mostly MB events
    if( runNumber == 312649 ){ goodRunLBN = false; }
    if( runNumber == 312796 && LBN <= 194 ){ goodRunLBN = false; }
    // in pp
    // skip these LB in the following run
    if( runNumber == 286282 && LBN <= 726 ){ goodRunLBN = false; }

    bool triggerPassed = false;
    // loop over passed triggers
    for( uint iG = 0; iG < m_nTriggers; iG++ ){

      // check if we have that trigger
      // if we dont - continue
      if( !mTriggerFired[iG] ) {
	mTrigFailed[iG]++;
	continue;	
      }
      
      const TLorentzVector* jet1 = NULL; const TLorentzVector* jet2 = NULL;
      // dPhi for all triggers, matched
      if( GetDiJets( vR_jets, jet1, jet2 ) ){
	if( JetInTrigRange( jet1, iG ) ) {
	  AnalyzeDeltaPhi( m_vHtriggerDphi[iG], vR_jets );
	}
      }

      // fill leading jet spectra for later.
      // used in yields normalization
      if( JetInTrigRange( jet1, iG ) ){
	double jetPt1    = jet1->Pt()/1000.;
	double jetYstar1 = GetYstar( *jet1 ); 

	triggerPassed = true;
	
	m_vHtriggerYstarSpect[iG]->
	  Fill( jetYstar1, jetPt1 );
	m_vHtriggerYstarSpectFine[iG]->
	  Fill( jetYstar1, jetPt1 );
	if( !m_is_pPb ){
	  m_vHtriggerYstarSpect[iG]->
	    Fill( -jetYstar1, jetPt1 );
	  m_vHtriggerYstarSpectFine[iG]->
	    Fill( -jetYstar1, jetPt1 );
	}
      }

      // loop over jets 
      for( auto& jet : vR_jets ){

	double jetPt    = jet.Pt()/1000.;
	double jetEta   = jet.Eta();
	double jetPhi   = jet.Phi();
	double jetYstar = GetYstar( jet );

	if( triggerPassed ){
	  if( jetPt > 28 && jetPt < 35 ){
	    m_vHtriggerEtaPhiMap[iG]->Fill( jetEta, jetPhi );
	    m_vHtriggerEtaPtMap [iG]->Fill( jetEta, jetPt  ); 
	  }
	  m_vHtriggerEtaPhiPtMap[iG]->Fill( jetEta, jetPhi, jetPt );
	}

	// check if jet is in trigger range
	if( !JetInTrigRange( &jet, iG ) ) continue;
	
	double rn = m_mRunBin[ runNumber ] - 0.5; 
	m_vHtriggerNjetsRun    [iG]->Fill( jetYstar, jetPt, rn );
	m_vHtriggerNjetsRunFine[iG]->Fill( jetYstar, jetPt, rn );
      } // end loop over jets
    } // end loop over iG

    // EFFICIENCIES - only for good run and LBN
    if( !goodRunLBN ){ continue; }
 
    AnalyzeEff( vR_jets, vTrig_jets, mTriggerFired ); 
    
  } // -------- END EVENT LOOP ---------
  std::cout << "DONE! Has: " << nEventsTotal << " events." << std::endl;

  fIn->Close(); delete fIn;
}

//---------------------------
//       Analysis
//---------------------------

bool DiJetAnalysisData::JetInTrigPtRange( const TLorentzVector* jet, int iG,
					  double ptHighExtra ){  

  if( !jet ){ return false; }
  
  double jetPt = jet->Pt() / 1000.;
  if( jetPt > m_vTriggersEffPtLow[iG] &&
      ( jetPt < m_vTriggersEffPtHigh[iG] + ptHighExtra ||
	m_vTriggersEffPtHigh[iG] < 0 ) )
    { return true; }
  return false;
}

bool DiJetAnalysisData::JetInTrigEtaRange( const TLorentzVector* jet, int iG ){

  if( !jet ){ return false; }
  
  if( iG == m_mbTriggerI ){ return true; }
  double jetEta = jet->Eta();
  double etaMin = m_vTriggersEtaMin[iG];
  double etaMax = m_vTriggersEtaMax[iG];
  // for pPb check only one side
  if( m_is_pPb && jetEta > -etaMax && jetEta < -etaMin )
    { return true; }
  // for pp check both sides
  else if ( std::abs(jetEta) < etaMax && std::abs(jetEta) > etaMin )
    { return true; }
  return false;
}

bool DiJetAnalysisData::JetInTrigRange( const TLorentzVector* jet, int iG ){

  if( !jet ){ return false; }
  
  double jetEta = jet->Eta();

  // in central triggers, we begin with j20
  // this corresponds to one efficiency
  // in forward, we begin with j25
  // this corresponds to another efficiency.
  // so we adjust here where the mb range goes to.
  bool applyCentMbCorrection =
    ( iG == m_mbTriggerI && IsCentralDetector( jetEta ) ) ? true : false;
  
  if( JetInTrigPtRange
      ( jet, iG, applyCentMbCorrection ? m_centMbCorrection : 0 ) &&
      JetInTrigEtaRange( jet, iG ) )
    { return true; }
  
  return false;
}

bool DiJetAnalysisData::TrigJetAboveThold( const TLorentzVector* jet, int iG ){

  if( !jet ){ return false; }
  
  if( jet->Pt()/1000. >= m_vTriggersTholdPt[iG] ){ return true; }
  return false;
}

bool DiJetAnalysisData::TrigJetInTrigRange( const TLorentzVector* jet, int iG ){

  if( !jet ){ return false; }
  
  if( TrigJetAboveThold( jet, iG ) && JetInTrigEtaRange( jet, iG ) )
    { return true; }
  return false;
}

void DiJetAnalysisData::AnalyzeEff( std::vector< TLorentzVector >& vR_jets,
				    std::vector< TLorentzVector >& vTrig_jets,
				    std::map< int, bool >& mTriggerFired ){

  if( vTrig_jets.empty() ){ return; }

  // sort trigger jets
  std::sort( vTrig_jets.begin(), vTrig_jets.end(), anaTool->sortByDecendingPt );

  // take highest pt trigger jet in
  // fwd and central eta ranges 
  // make sure to find at least one 
  const TLorentzVector* tJetFwd  = NULL; const TLorentzVector* tJetCent = NULL;    
  if( !GetFwdCentJets( vTrig_jets, tJetFwd, tJetCent ) ){ return; }

  // take highest pt reco jet in
  // fwd and central eta ranges.
  // make sure to find at least one 
  const TLorentzVector* rJetFwd  = NULL; const TLorentzVector* rJetCent = NULL;
  if( !GetFwdCentJets( vR_jets, rJetFwd, rJetCent ) ){ return; }
  
  // first fill histogram for denominator. I.e. did a trigger fire
  for( uint iG = 0; iG < m_nTriggers; iG++ ){      
    if( !mTriggerFired[iG] ){ continue; }
    FillHistoWithJets( rJetFwd, rJetCent, m_vHtriggerEtaSpectDenom[iG] );
  }

  // now fill histos for triggers that could have passed (emulate)
  // with the condition that the reference trigger fired too
  for( uint iG = 0; iG < m_nTriggers; iG++ ){      
    // check if this triggers reference trigger fired
    int refTrigIndex = m_vTriggersRefIndex[iG];
    if( !mTriggerFired[ refTrigIndex ] ){ continue; }
    
    // emulate that trigger. Check if it should fire
    // only if trigger jet is above threshold
    // and in eta range for that trigger
    // if we had both a forward and central jet for example, and
    // a forward and central trigger, this will pass on two occasions
    if( !( tJetFwd  ? TrigJetInTrigRange( tJetFwd , iG ) : false ) &&
	!( tJetCent ? TrigJetInTrigRange( tJetCent, iG ) : false ) )
      { continue; }

    // if we have one or the other, fill the
    // simulated spectra histograms
    FillHistoWithJets( rJetFwd, rJetCent, m_vHtriggerEtaSpectSim[iG] );
  } // end loop over iG
}

//---------------------------
//       Tools
//---------------------------  

void DiJetAnalysisData::CleanEfficiency( TGraphAsymmErrors* g, int iG ){
  double x, y;
  for( int i = 0; i < g->GetN(); i++ ){
    if( g->GetPoint( i, x, y ) == -1 ){ continue; }
    if( x > ( m_vTriggersTholdPt[iG] - 15 ) ){ continue; }
    g->SetPoint     ( i, x, 0 );
    g->SetPointError( i, 0, 0, 0, 0 );
  }
}

TH1* DiJetAnalysisData::CombineSamples( std::vector< TH1* >& vSampleHin,
					const std::string& name ){

  if( !vSampleHin.size() ){ return NULL; }

  TH1* h_res = static_cast< TH1D* >
    ( vSampleHin[0]->Clone( Form("h_%s_%s", name.c_str() , m_allName.c_str() ) ) );
  h_res->Reset();
  
  for( uint iG = 0; iG < vSampleHin.size(); iG++ ){
    double scale = m_vTriggersPrescale[iG];
    h_res->Add( vSampleHin[iG], scale );
  }
  vSampleHin.push_back( h_res );
  
  return h_res;
}

TH2* DiJetAnalysisData::CombineSamples( std::vector< TH2* >& vSampleHin,
					const std::string& name ){

  if( !vSampleHin.size() ){ return NULL; }

  TH2* h_res = static_cast< TH2D* >
    ( vSampleHin[0]->Clone( Form("h_%s_%s", name.c_str() , m_allName.c_str() ) ) );
  h_res->Reset();
    
  for( uint iG = 0; iG < vSampleHin.size(); iG++ ){
    double scale = m_vTriggersPrescale[iG];
    h_res->Add( vSampleHin[iG], scale );
  }
  vSampleHin.push_back( h_res );
    
  return h_res;
}

TH3* DiJetAnalysisData::CombineSamples( std::vector< TH3* >& vSampleHin,
					const std::string& name,
					bool applyPS ){

  if( !vSampleHin.size() ){ return NULL; }

  TH3* h_res = static_cast< TH3D* >
    ( vSampleHin[0]->Clone( Form("h_%s_%s", name.c_str() , m_allName.c_str() ) ) );
  h_res->Reset();
    
  for( uint iG = 0; iG < vSampleHin.size(); iG++ ){
    double scale = applyPS ? m_vTriggersPrescale[iG] : 1;
    h_res->Add( vSampleHin[iG], scale );
  }
  vSampleHin.push_back( h_res );
    
  return h_res;
}

THnSparse* DiJetAnalysisData::CombineSamples( std::vector< THnSparse* >& vSampleHin,
					      const std::string& name  ){

  if( !vSampleHin.size() ){ return NULL; }
  
  THnSparse* h_res = static_cast< THnSparseD* >
    ( vSampleHin[0]->Clone( Form("h_%s_%s", name.c_str() , m_allName.c_str() ) ) );
  h_res->Reset();
  
  std::vector< int > x { 1, 2, 1, 1, 6, 1 };
  std::vector< int > xR{ 1, 2, 1, 1, 7, 1 }; 
  std::vector< int > xL{ 1, 2, 1, 1, 5, 1 }; 
  
  for( uint iG = 0; iG < vSampleHin.size(); iG++ ){
    THnSparse* hTrig = vSampleHin[iG];
    double scale = m_vTriggersPrescale[iG];
    anaTool->SetZeroEntryError( hTrig );
    h_res->Add( hTrig, scale );
    std::cout << m_vTriggers[ iG ] << " - " << hTrig->GetBinContent( &x[0] ) << " "
	      << hTrig->GetBinError( &x[0] ) << " " << scale << std::endl;
    std::cout << "   L - " << hTrig->GetBinContent( &xL[0] ) << " "
	      << hTrig->GetBinError( &xL[0] ) << std::endl;
    std::cout << "   R - " << hTrig->GetBinContent( &xR[0] ) << " "
	      << hTrig->GetBinError( &xR[0] ) << std::endl;
  }
  vSampleHin.push_back( h_res );

  std::cout << h_res->GetBinContent( &x[0] ) << " "
	    << h_res->GetBinError  ( &x[0] ) << std::endl;
  std::cout << "   L - " << h_res->GetBinContent( &xL[0] ) << " "
	    << h_res->GetBinError  ( &xL[0] ) << std::endl;
  std::cout << "   R - " <<  h_res->GetBinContent( &xR[0] ) << " "
	    << h_res->GetBinError  ( &xR[0] ) << std::endl;

  return h_res;
}

void DiJetAnalysisData::GetSpectUnfoldingInfo( std::string& measuredName,
					       std::string& recoName,
					       std::string& truthName,
					       std::string& respMatName,
					       std::string& unfoldedLabel,
					       std::string& typeLabel ){

  measuredName  = m_ystarSpectName;
  recoName      = m_ystarSpectRecoName;
  truthName     = m_ystarSpectTruthName;
  respMatName   = m_ystarSpectRespMatName;
  unfoldedLabel = "dN/d#it{p}_{T} [GeV]";;
  typeLabel     = "Data";
}

void DiJetAnalysisData::GetDphiUnfoldingInfo( std::string& measuredName,
					      std::string& truthName,
					      std::string& unfoldedLabel,
					      std::string& typeLabel ){

  measuredName  = m_dPhiName;
  truthName     = m_dPhiTruthName;
  unfoldedLabel = "|#Delta#phi|";
  typeLabel     = "Data";
}

void DiJetAnalysisData::GetInfoTogether( std::string& name_a , std::string& name_b ,
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
  case 2:
    pFname = &m_fNameSYS;
    name_a = m_dPhiName;
    name_b = m_dPhiName;
    break;
  }
  
  int combinationBoth = GetConfig()->GetValue( "combinationBoth", 0 );

  // for the final plots, name_a and name_b
  // are same in both pp and pPb files
  // so instead, they are the name of the nominal
  // and systematic TGraphs.
  if( combinationBoth == 0 ){
    name_a    += "_" + m_unfoldedName;
    name_b    += "_" + m_systematicsName;
    label_a   = "#it{p}+Pb";
    label_b   = "#it{pp}";
    m_is_pPb  = true;  Initialize();
    fName_a   = *pFname;
    m_is_pPb  = false; Initialize();
    fName_b   = *pFname;
  } else if( combinationBoth == 1 ){
    name_a    += "_" + m_unfoldedName;
    name_b    += "_" + m_unfoldedName;
    label_a   = "#it{p}+Pb";
    label_b   = "#it{pp}";
    m_is_pPb  = true;  Initialize();
    fName_a   = *pFname;
    m_is_pPb  = false; Initialize();
    fName_b   = *pFname;
  } else if( combinationBoth == 2 ){
    label_a   = "Rec #it{p}+Pb";
    label_b   = "Rec #it{pp}";
    m_is_pPb  = true;  Initialize();
    fName_a   = *pFname;
    m_is_pPb  = false; Initialize();
    fName_b   = *pFname;
  } else if ( combinationBoth == 3 ){
    name_a    += "_" + m_unfoldedName;
    label_a   = "UF #it{pp}";
    label_b   = "Rec #it{pp}";
    m_is_pPb  = false; Initialize();
    fName_a   = *pFname;
    fName_b   = *pFname;
  } else if ( combinationBoth == 4 ){
    name_a    += "_" + m_unfoldedName;
    label_a   = "UF #it{p}+Pb";
    label_b   = "Rec #it{p}+Pb";
    m_is_pPb  = true; Initialize();
    fName_a   = *pFname;
    fName_b   = *pFname;
  }
}

//---------------------------------
//     Get Quantities / Plot 
//---------------------------------

void DiJetAnalysisData::LoadHistograms( int opt ){
  
  TFile* fIn = TFile::Open( m_fNameDefRaw.c_str() ); 

  for( auto& trigger : m_vTriggers ){
    // -------- maps ---------
    m_vHtriggerEtaPhiMap.push_back
      ( static_cast< TH2D* >
	( fIn-> Get( Form("h_etaPhiMap_%s",trigger.c_str() ))));
    m_vHtriggerEtaPhiMap.back()->SetDirectory(0);
    
    m_vHtriggerEtaPtMap.push_back
      ( static_cast< TH2D* >( fIn->Get( Form("h_etaPtMap_%s", trigger.c_str() ))));
    m_vHtriggerEtaPtMap.back()->SetDirectory(0);
    
    // -------- spect --------
    m_vHtriggerYstarSpect.push_back
      ( static_cast< TH2D* >
	( fIn->Get( Form("h_%s_%s", m_ystarSpectName.c_str(), trigger.c_str() ))));
    m_vHtriggerYstarSpect.back()->SetDirectory(0);

    m_vHtriggerYstarSpectFine.push_back
      ( static_cast< TH2D* >
	( fIn->Get( Form("h_%s_%s", m_ystarSpectFineName.c_str(), trigger.c_str() ))));
    m_vHtriggerYstarSpectFine.back()->SetDirectory(0);

    m_vHtriggerNjetsRun.push_back
      ( static_cast< TH3D* >
	( fIn->Get( Form("h_%s_%s", m_nJetsRunName.c_str(), trigger.c_str() ))));
    m_vHtriggerNjetsRun.back()->SetDirectory(0);

    m_vHtriggerNjetsRunFine.push_back
      ( static_cast< TH3D* >
	( fIn->Get( Form("h_%s_fine_%s", m_nJetsRunName.c_str(), trigger.c_str() ))));
    m_vHtriggerNjetsRunFine.back()->SetDirectory(0);
// ----- efficiencies ----
    m_vHtriggerEtaSpectSim.push_back
      ( static_cast< TH2D* >
	( fIn->Get( Form("h_%sSim_%s", m_etaSpectName.c_str(), trigger.c_str() ))));
    m_vHtriggerEtaSpectSim.back()->SetDirectory(0);

    m_vHtriggerEtaSpectDenom.push_back
      ( static_cast< TH2D* >
	( fIn->Get( Form("h_%sDenom_%s", m_etaSpectName.c_str(), trigger.c_str() ))));
    m_vHtriggerEtaSpectDenom.back()->SetDirectory(0);
    // -------- dPhi- --------
    m_vHtriggerDphi.push_back
      ( static_cast< THnSparse* >
	( fIn->Get( Form("h_%s_%s", m_dPhiName.c_str(), trigger.c_str() ))));
  }
  
  fIn->Close(); delete fIn;
}
					  
void DiJetAnalysisData::MakeEfficiencies( std::vector< TH2* >& vTrigSpect,
					  std::vector< TH2* >& vTrigSpectRef,
					  const std::string& type ){

  // to use for reference on axis etc.
  TH2* hExample = NULL;
  if( vTrigSpect.size() ){
    hExample = vTrigSpect.front();
  } else{
    return;
  }

  double lX0 = 0.18;
  double lY0 = 0.70;
  double lX1 = 0.40;
  double lY1 = 0.90;

  double legScale = 0.75;
  
  if( m_is_pPb ){
    legScale = 0.7;
    lY0 = 0.80;
    lX1 = 0.30;
    lY1 = 0.90;
  }
  
  // us m_hAllEtaSpect because its always there
  double xMin = 0;
  double xMax = 100;

  // use this as reference because
  // it should be in every file
  int nXbins = hExample->GetNbinsX();
  
  std::vector< std::vector< TH1* > > vSpect;
  std::vector< std::vector< TH1* > > vSpectRef;
  std::vector< std::vector< TGraphAsymmErrors* > > vEffGrf;
  
  vSpect   .resize( m_nTriggers );
  vSpectRef.resize( m_nTriggers );
  vEffGrf  .resize( m_nTriggers );
  
  for( uint iG = 0; iG < m_nTriggers; iG++ ){
    std::string trigger = m_vTriggers[iG];
    
    for( int xBin = 1; xBin <= nXbins; xBin++ ){
      // should all be the same
      double etaMin, etaMax;
      anaTool->GetBinRange
	( hExample->GetXaxis(), xBin, xBin, etaMin, etaMax );
            
      std::string etaLabel = anaTool->GetEtaLabel
	( etaMin, etaMax, m_is_pPb );
      std::string etaName  = anaTool->GetName
	( etaMin, etaMax, "Eta" );
      
      TH1* h_etaSpect = vTrigSpect[iG]->
	ProjectionY( Form("h_%s_%s_%s",
			  type.c_str(),
			  etaName.c_str(),
			  trigger.c_str() ),
		     xBin, xBin );
      h_etaSpect->SetTitle( etaLabel.c_str() );
      vSpect[iG].push_back( h_etaSpect );

      TH1* h_etaSpectRef = vTrigSpectRef[iG]->
	ProjectionY( Form("h_%s_%s_ref_%s",
			  type.c_str(),
			  etaName.c_str(),
			  trigger.c_str() ),
		     xBin, xBin );
      h_etaSpectRef->SetTitle( etaLabel.c_str() );
      vSpectRef[iG].push_back( h_etaSpectRef );

      TGraphAsymmErrors* g_etaEff = new TGraphAsymmErrors();
      g_etaEff->SetName( Form("gr_%s_%s_%s",
			      type.c_str(),
			      etaName.c_str(),
			      trigger.c_str() ) );
      g_etaEff->SetTitle( etaLabel.c_str() );
      vEffGrf[iG].push_back( g_etaEff );      
    } // end loop over eta
  } // end loop over triggers

  //------------------------------------------------
  //--------- Now Divide by Appropriate Ref  -------
  //------------------------------------------------
  for( uint iG = 0; iG < m_nTriggers; iG++ ){
    for( int iX = 0; iX < nXbins; iX++ ){
      int  refTrigIndex = m_vTriggersRefIndex[iG];
      TH1* hSpect    = vSpect[iG][iX];
      TH1* hSpectRef = vSpectRef[refTrigIndex][iX];
      vEffGrf[iG][iX]->Divide( hSpect, hSpectRef );
      CleanEfficiency( vEffGrf[iG][iX], iG );
    }
  }
  
  //------------------------------------------------
  // ----------------- Now Plot --------------------
  //------------------------------------------------
  std::string gLabel = ";#it{p}_{T} [GeV];#it{#varepsilon}_{Trigger}";
  
  //------------------------------------------------
  //---------- Draw Triggers in Eta Bins -----------
  //------------------------------------------------
  for( int iX = 0; iX < nXbins; iX++ ){
    // no longer have title.
    // Need to make output files
    // with eta names
    int    xBin = iX + 1;

    double etaMin, etaMax;
    anaTool->GetBinRange
      ( hExample->GetXaxis(), xBin, xBin, etaMin, etaMax );
    double etaCenter = hExample->GetXaxis()->GetBinCenter ( xBin );
    
    // temporary, dont draw the 3.1->3.2 bin
    if( std::abs(etaCenter) < 3.2 && std::abs(etaCenter) > 3.1 ){ continue; }
    
    // for pPb, dont draw at anything above -3.2
    if( m_is_pPb && etaCenter > -constants::FETAMIN ){ continue; }
    
    std::string cName  = anaTool->GetName( etaMin, etaMax, "Eta" );
    std::string cLabel = anaTool->GetEtaLabel( etaMin, etaMax, m_is_pPb );
    
    TCanvas c( cLabel.c_str(), cLabel.c_str(), 800, 600 );
    styleTool->SetCStyleGraph( c, xMin, m_effMin, xMax, m_effMax, gLabel );
    
    TLegend leg( lX0, lY0, lX1, lY1 );
    styleTool->SetLegendStyle( &leg, legScale );
    
    int style = 0;
    for( uint iG = 0; iG < m_nTriggers; iG++ ){
      std::string trigger = m_vTriggers[iG];

      // dont draw MB trigger
      if( !trigger.compare( m_mbTriggerName ) )
	{ continue; }
      if( !trigger.compare( m_allName ) )
	{ continue; }

      // for pp, dont draw central triggers below -3.2
      // or forward triggers above -3.2
      if( !m_is_pPb && etaCenter < -constants::FETAMIN &&
	  trigger.find("320eta490") == std::string::npos ){ continue; }
      else if( !m_is_pPb && etaCenter > -constants::FETAMIN &&
	       trigger.find("320eta490") != std::string::npos){ continue; }
      
      TGraphAsymmErrors* g = vEffGrf[iG][iX];
      styleTool->SetHStyle( g, style++ );      
      g->Draw("p");
      leg.AddEntry( g, trigger.c_str(), "lp" );
    } // end loop over iG
    
    leg.Draw("same");

    TLine line( xMin, 1, xMax, 1);
    line.Draw();    
    
    DrawAtlasRight( 0, 0, 0.85 );
    if( !m_is_pPb && iX == 1 ){
      drawTool->DrawRightLatex( 0.65, 0.725, "-3.2<|#eta|<3.2" );
    } else {
      drawTool->DrawRightLatex( 0.65, 0.725, cLabel, 0.85 );
    }
    
    SaveAsAll( c, Form("%s_%s", type.c_str(), cName.c_str() ) );
  } // end loop over iX

  // delete
  for( uint iG = 0; iG < m_nTriggers; iG++ ){
    for( int iX = 0; iX < nXbins; iX++ ) {
      delete vSpect[iG][iX];
      delete vSpectRef[iG][iX];
      delete vEffGrf[iG][iX];
    }
  } 
}

void DiJetAnalysisData::MakeNjetsRun( std::vector< TH3* >& vhSample,
				      const std::vector< std::string >& vLabels,
				      const std::string& name ){
  
  for( uint iG = 0; iG < m_nTriggers; iG++ ){
    
    std::string label = vLabels[ iG ];

    if( !label.compare( m_allName ) ){ continue; }
    if( label.find("_mb_") != std::string::npos ){ continue; }
    // if( !m_is_pPb && label.compare("HLT_j25_320eta490_L1TE5") ){ continue; }
    // if( !m_is_pPb && label.compare("HLT_j35_320eta490_L1TE10") ){ continue; }
    // if( !m_is_pPb && label.compare("HLT_j45_320eta490") ){ continue; }
    // if( !m_is_pPb && label.compare("HLT_j30_L1TE5") ){ continue; }
    // if( !m_is_pPb && label.compare("HLT_j30_L1TE5") ){ continue; }
    
    TH3* hSample = vhSample[ iG ];
    
    int xBin = 2, yBin = 1; 
    double xMin, xMax;
    anaTool->GetBinRange
      ( hSample->GetXaxis(), xBin, xBin, xMin, xMax );
    double yMin, yMax;
    anaTool->GetBinRange
      ( hSample->GetYaxis(), yBin, yBin, yMin, yMax );
      
    std::string hTag = Form( "%s_%s", anaTool->GetName( xMin, xMax, "Ystar1").c_str(),
			     anaTool->GetName( yMin, yMax, "Pt1").c_str() );

    std::string hName = "h_" + name + "_" + label + "_" + hTag;
    hSample->GetXaxis()->SetRange( xBin, xBin );
    hSample->GetYaxis()->SetRange( yBin, yBin );

    TH1* hNjetsRun = static_cast< TH1D* >
      ( hSample->ProjectionZ( hName.c_str(), xBin, xBin, yBin, yBin ) );

    double nJetsTot = 0;
    
    std::cout << "+++++++++++ " << label << " - " << hName << std::endl;
    
    for( auto& rB : m_mRunBin ){
      int run = rB.first;
      int bin = rB.second;
      double lumi = m_mRunLumi[ run ];

      double nJets      = hNjetsRun->GetBinContent( bin );
      double nJetsError = hNjetsRun->GetBinError  ( bin );

      nJetsTot += nJets;
      
      double nJetsScaled      = nJets     /lumi;
      double nJetsErrorScaled = nJetsError/lumi;
      
      std::cout << bin << " " << run << " " << lumi << " " << nJets << " " << nJetsScaled << std::endl;

      if( !nJets ){ continue; }

      hNjetsRun->GetXaxis()->SetBinLabel( bin, Form("%d", run ) );
      hNjetsRun->SetBinContent( bin, nJetsScaled      );
      hNjetsRun->SetBinError  ( bin, nJetsErrorScaled );
    }
    std::cout << "Had: " << nJetsTot << " jets total." << std::endl;

    TCanvas c( "c", "c", 800, 600 );
    styleTool->SetHStyle( hNjetsRun, 0 );
    hNjetsRun->SetMinimum( 0 );
    hNjetsRun->SetMaximum( hNjetsRun->GetMaximum() * 1.6 );

    hNjetsRun->SetYTitle( "N_{jets}/L_{run}");
    hNjetsRun->Draw();

    drawTool->DrawLeftLatex( 0.2, 0.3, label );


    drawTool->DrawLeftLatex
      ( 0.18, 0.86, anaTool->GetLabel( xMin, xMax, "#it{y}_{1}*" ), 1.0 );
    drawTool->DrawLeftLatex
      ( 0.18, 0.78, anaTool->GetLabel( yMin, yMax, "#it{p}_{T}^{1}"), 1.0 );
    
    DrawAtlasRight();

    SaveAsAll( c, hName );
    
    hNjetsRun->Write();
  }
}

void DiJetAnalysisData::MakeSystematicsGraphs( TFile* fOut, const std::string& name ){

  std::vector< int > v_uc;
  std::map< int, TFile* > mFinUC;

  // see if we are dealing with widths of yields.
  bool isYield = name.find( m_yieldName ) != std::string::npos ? true : false;
  
  // get map of sys factor to TFile*
  TFile* fInNominal = GetListOfSystUncert( v_uc, mFinUC );
  
  std::string allUnfoldedName    =
    m_dPhiName + "_" + m_unfoldedName    + "_" + name + "_" + m_allName;
  std::string allSystematicsName =
    m_dPhiName + "_" + m_systematicsName + "_" + name + "_" + m_allName;

  std::vector< TH1* > vHdef;
  std::vector< TH1* > vHsyst;
  std::vector< TGraphAsymmErrors* > vG;
  
  TAxis* axis0 = m_dPP->GetTAxis(0); int nAxis0Bins = axis0->GetNbins();
  TAxis* axis1 = m_dPP->GetTAxis(1); int nAxis1Bins = axis1->GetNbins();
  TAxis* axis2 = m_dPP->GetTAxis(2); int nAxis2Bins = axis2->GetNbins();
  TAxis* axis3 = m_dPP->GetTAxis(3); int nAxis3Bins = axis3->GetNbins();

  // for canvas, since its tgraphasymmerrors.
  double x0 = axis3->GetXmin();
  double x1 = axis3->GetXmax(); 

  double y0 = isYield ? m_dPhiYieldMin : m_dPhiWidthMin;
  double y1 = isYield ? m_dPhiYieldMax : m_dPhiWidthMax;

  double pDx1 = 0.2;
  
  std::string xTitle = m_dPP->GetAxisLabel(3);
  std::string yTitle = isYield ? m_sYieldTitle : m_sWidthTitle;
  std::string gTitle = ";" + xTitle + ";" + yTitle;
  
  std::string yTitleUncert = "#delta " + yTitle + " %";

  std::string sRatio       = isYield ? m_sYieldRatioTitle : m_sWidthRatioTitle;
  std::string gTitleRatio  = ";" + xTitle + ";Data/MC";

  TH1* hDef  = new TH1D( "hDef", gTitle.c_str(), nAxis3Bins, x0, x1 );
  styleTool->SetHStyle( hDef, 0 );
  hDef->SetMaximum( y1 );
  hDef->SetMinimum( y0 );
  hDef->GetXaxis()->SetTitle("");
  hDef->GetYaxis()->SetNdivisions( 502 );

  TH1* hDefR = new TH1D( "hDefR", gTitleRatio.c_str(), nAxis3Bins, x0, x1 );
  styleTool->SetHStyleRatio( hDefR, 0 );
  hDefR->SetMaximum( 1.55 );
  hDefR->SetMinimum( 0.45 );
  
  // lines to be drawn along axis3. this is
  // x-axis that widths are plotted as function of 
  TLine line( x0, 1, x1, 1 );
  line.SetLineWidth( 2 );

  TLine lineP25( x0, 1.25, x1, 1.25 );
  lineP25.SetLineStyle( 2  );
  lineP25.SetLineColor( 12 );
  lineP25.SetLineWidth( 2  );
	  
  TLine lineN25( x0, 0.75, x1, 0.75 );
  lineN25.SetLineStyle( 2  );
  lineN25.SetLineColor( 12 );
  lineN25.SetLineWidth( 2  );

  
  // ------------------------------------------------
  //    MC STUFF - TO COMPARE DATA RESULTS TO MC
  // ------------------------------------------------
  TFile* fMCp8 = m_is_pPb ?
    TFile::Open( "output/output_pPb_mc_pythia8/myOut_pPb_mc_pythia8_phys_0.root" ) :
    TFile::Open( "output/output_pp_mc_pythia8/myOut_pp_mc_pythia8_phys_0.root" );

  std::string hMCname = isYield ? "dPhi_truth_yield_All" : "dPhi_truth_width_All";
  
  // change back to the fOut, for writing.
  fOut->cd();
  
  for( int axis0Bin = 1; axis0Bin <= nAxis0Bins; axis0Bin++ ){
    double axis0Low, axis0Up;
    anaTool->GetBinRange
      ( axis0, axis0Bin, axis0Bin, axis0Low, axis0Up );
    
    for( int axis1Bin = 1; axis1Bin <= nAxis1Bins; axis1Bin++ ){
      if( !m_dPP->CorrectPhaseSpace
	  ( std::vector<int>{ axis0Bin, axis1Bin, 0, 0 } ) )
	{ continue; }
      
      double axis1Low, axis1Up;
      anaTool->GetBinRange
	( axis1, axis1Bin, axis1Bin, axis1Low, axis1Up );

      TCanvas cAll( "cAll", "cAll", 800, 700 );
      styleTool->SetCStyleGraph( cAll, x0, y0, x1, y1, gTitle.c_str() );
      TPad pad1All( "pad1All", "", 0.0, 0.35, 1.0, 1.0 );
      pad1All.SetBottomMargin(0.023);
      pad1All.Draw();
      TPad pad2All( "pad2All", "", 0.0, 0.0, 1.0, 0.34 );
      pad2All.SetTopMargin(0.05);
      pad2All.SetBottomMargin(0.27);
      pad2All.Draw();

      pad1All.cd();
      hDef->Draw();
      styleTool->HideAxis( hDef, "x" );
      pad2All.cd();
      hDefR->Draw();

      if( axis1Bin == 3 && isYield ){
	hDef->SetMaximum( 2E-3 );
	hDef->SetMinimum( 5E-7 );
      }
      
      // for yields, use log scale on yaxis
      if( isYield ){ pad1All.SetLogy(); }

      // this is bs. works though. vary the last factor 
      double deltaYleg = ( axis1Bin - 1 ) * 0.04;

      double legSX0MC = isYield ? 0.34 : 0.19;
      double legSX0Data  = isYield ? 0.45 : 0.30;
      double legSX1MC = isYield ? 0.50 : 0.35;
      double legSX1Data  = isYield ? 0.73 : 0.47;

      double legSY0 = isYield ? 0.025 : 0.025;
      double legSY1 = isYield ? 0.13 + deltaYleg : 0.13 + deltaYleg;
      
      TLegend legMC( legSX0MC, legSY0, legSX1MC, legSY1 );
      styleTool->SetLegendStyle( &legMC, 0.85 );
      TLegend legData( legSX0Data,  legSY0,  legSX1Data, legSY1 );
      styleTool->SetLegendStyle( &legData, 0.85 );
      

      // tag for widths canvas
      std::string hTagC =
	Form ("%s_%s",
	      anaTool->GetName( axis0Low, axis0Up, m_dPP->GetAxisName(0) ).c_str(),
	      anaTool->GetName( axis1Low, axis1Up, m_dPP->GetAxisName(1) ).c_str() );

      int style = 0;
      // ---- loop over axis2 ----
      for( int axis2Bin = 1; axis2Bin <= nAxis2Bins; axis2Bin++ ){
	if( !m_dPP->CorrectPhaseSpace
	    ( std::vector<int>{ axis0Bin, axis1Bin, axis2Bin, 0 } ) )
	  { continue; }

	double axis2Low , axis2Up;
	anaTool->GetBinRange
	  ( axis2, axis2Bin, axis2Bin, axis2Low, axis2Up );

	std::string hTag =
	  Form( "%s_%s_%s",
	        anaTool->GetName( axis0Low, axis0Up, m_dPP->GetAxisName(0) ).c_str(),
		anaTool->GetName( axis1Low, axis1Up, m_dPP->GetAxisName(1) ).c_str(),
		anaTool->GetName( axis2Low, axis2Up, m_dPP->GetAxisName(2) ).c_str() ); 
	
	std::string hNominalName = "h_" + allUnfoldedName + "_" + hTag;
	std::string gNominalName = "g_" + allUnfoldedName + "_" + hTag;
	
	TH1* hNominal = static_cast<TH1D*>( fInNominal->Get( hNominalName.c_str() ) );
       	vHdef.push_back( hNominal );

	// Make all the systematics histograms for each y1, pt1, pt2 bin.
	// since we look one bin at a time and process all systematics,
	// need to have these histograms ready. Store them to a map
	// with the key as the systematic component number.
	std::map< int, TH1* > mHsystTmp; 
	for( auto uc : v_uc ){

	  // skip uc = 0 (default)
	  if( !uc ){ continue; }

	  std::string tmpUncertSuffix = uc > 0 ? Form("P%d", uc) : Form("N%d", -1 * uc) ;
	  std::string hSystematicName = "h_" + allUnfoldedName + "_" + hTag + "_" + tmpUncertSuffix;

	  TH1* hSystematic = static_cast<TH1D*>( hNominal->Clone( hSystematicName.c_str() ) );
	  vHsyst.push_back( hSystematic );
	  hSystematic->Reset();
	  mHsystTmp[ uc ] = hSystematic;
	}

	double systMax = 60;
	double systMin = -50;
	
	// Make all the systematics histograms for each y1, pt1, pt2 bin.
	// These are more specific (PJER, PJER, HIJES, TOT, etc )
	std::string hPTOTName = "h_" + allSystematicsName + "_PTOT_" + hTag;
	TH1* hPTOT = static_cast< TH1D* >( hNominal->Clone( hPTOTName.c_str() ) );
	vHsyst.push_back( hPTOT );
	hPTOT->Reset();
	styleTool->SetHStyle( hPTOT, 0 );
	
	std::string hNTOTName = "h_" + allSystematicsName + "_NTOT_" + hTag;
	TH1* hNTOT = static_cast< TH1D* >( hNominal->Clone( hNTOTName.c_str() ) );
	vHsyst.push_back( hNTOT );
	hNTOT->Reset();
	styleTool->SetHStyle( hNTOT, 0 );

	std::string hPJESName = "h_" + allSystematicsName + "_PJES_" + hTag;
	TH1* hPJES = static_cast< TH1D* >( hNominal->Clone( hPJESName.c_str() ) );
	vHsyst.push_back( hPJES );
	hPJES->Reset();
	styleTool->SetHStyle( hPJES, 1 );

	std::string hNJESName = "h_" + allSystematicsName + "_NJES_" + hTag;
	TH1* hNJES = static_cast< TH1D* >( hNominal->Clone( hNJESName.c_str() ) );
	vHsyst.push_back( hNJES );
	hNJES->Reset();
	styleTool->SetHStyle( hNJES, 1 );
	
	std::string hPHIJESName = "h_" + allSystematicsName + "_PHIJES_" + hTag;
	TH1* hPHIJES = static_cast< TH1D* >( hNominal->Clone( hPHIJESName.c_str() ) );
	vHsyst.push_back( hPHIJES );
	hPHIJES->Reset();
	styleTool->SetHStyle( hPHIJES, 2 );

	std::string hNHIJESName = "h_" + allSystematicsName + "_NHIJES_" + hTag;
	TH1* hNHIJES = static_cast< TH1D* >( hNominal->Clone( hNHIJESName.c_str() ) );
	vHsyst.push_back( hNHIJES );
	hNHIJES->Reset();
	styleTool->SetHStyle( hNHIJES, 2 );

	std::string hPJERName = "h_" + allSystematicsName + "_PJER_" + hTag;
	std::string hNJERName = "h_" + allSystematicsName + "_NJER_" + hTag;
	TH1* hPJER = static_cast< TH1D* >( hNominal->Clone( hPJERName.c_str() ) );
	vHsyst.push_back( hPJER );
	hPJER->Reset();
	styleTool->SetHStyle( hPJER, 4 );

	std::string hPANGName = "h_" + allSystematicsName + "_PANG_" + hTag;
	std::string hNANGName = "h_" + allSystematicsName + "_NANG_" + hTag;
	TH1* hPANG = static_cast< TH1D* >( hNominal->Clone( hPANGName.c_str() ) );
	vHsyst.push_back( hPANG );
	hPANG->Reset();
	styleTool->SetHStyle( hPANG, 3 );
		
	std::string hPUNFName = "h_" + allSystematicsName + "_PUNF_" + hTag;
	std::string hNUNFName = "h_" + allSystematicsName + "_NUNF_" + hTag;
	TH1* hPUNF = static_cast< TH1D* >( hNominal->Clone( hPUNFName.c_str() ) );
	vHsyst.push_back( hPUNF );
	hPUNF->Reset();
	styleTool->SetHStyle( hPUNF, 0 );
	hPUNF->SetLineStyle( 7 );

	std::string hPFITName = "h_" + allSystematicsName + "_PFIT_" + hTag;
	std::string hNFITName = "h_" + allSystematicsName + "_NFIT_" + hTag;
	TH1* hPFIT = static_cast< TH1D* >( hNominal->Clone( hPFITName.c_str() ) ); 
	vHsyst.push_back( hPFIT );
	hPFIT->Reset();
	styleTool->SetHStyle( hPFIT, 1 );
	hPFIT->SetLineStyle( 7 );
	
	std::vector< double > pX;
	std::vector< double > eX( nAxis3Bins, pDx1 * 0.5 );
 
	std::vector< double > pY;
	std::vector< double > eYP;
	std::vector< double > eYN;

	std::vector< double > eYPJES;
	std::vector< double > eYNJES;

	//---------------------------------------------------
	//------------------ DO WORK HERE -------------------
	//---------------------------------------------------
	// for( int axis3Bin = nAxis3Bins; axis3Bin >= 1; axis3Bin-- ){
	for( int axis3Bin = 1; axis3Bin <= nAxis3Bins; axis3Bin++ ){
	  
	  std::vector< TH1* > vHunc;

	  // add uncertainties in quadrature;
	  double uncertFinYP = 0;
	  double uncertFinYN = 0;

	  // add uncertainties in quadrature;
	  double uncertFinYPJES = 0;
	  double uncertFinYNJES = 0;

	  double yNominal = hNominal->GetBinContent( axis3Bin );
	  pX.push_back( hNominal->GetBinCenter( axis3Bin ) );
	  
	  // loop over uncertainties
	  for( auto uc : v_uc ){

	    // skip uc = 0 (default)
	    if( !uc ){ continue; }

	    // this is same as hNominalName but leave it separate
	    // to avoid confusion / make more flexible later.
	    std::string hUncertaintyName = "h_" + allUnfoldedName + "_" + hTag;

	    TH1* hUncertainty = static_cast<TH1D*>( mFinUC[ uc ]->Get( hUncertaintyName.c_str() ) );
	    vHunc.push_back( hUncertainty );

	    int sign = uc > 0 ? 1 : -1;

	    double yShifted = hUncertainty->GetBinContent( axis3Bin );

	    double uncertaintySq = std::pow( ( yNominal - yShifted ) / yNominal, 2 );
	    double uncertainty   = std::sqrt( uncertaintySq ); 
	    
	    std::cout << "++++" << uc << " ++++" << axis3Bin << " " << axis1Bin << " "
		      << axis2Bin << " " << sign << " " << yShifted << " "
		      << yNominal << " " << uncertainty << std::endl;

	    mHsystTmp[ uc ]->SetBinContent( axis3Bin, yShifted );
	    
	    if( !m_is_pPb && !isYield && axis1Bin == 3 && axis3Bin == 1 &&
		( axis2Bin == 1 || axis2Bin == 2 ) ){
	      continue;
	    }
	    
	    // for JER negative is same as positive (20)
	    // for Angular Res negative is same as positive (21)
	    // for Fitting negative is same as positive (22)
	    // for ReWeight/Unf negative is same as positive (23)
	    // and we do not have a NN, just use PN
	    if( uc  >= 20 && uc <= 23 ){
	      uncertFinYP += uncertaintySq;
	      uncertFinYN += uncertaintySq;

	      // set these uncertainties in % ( x100 ) directly, since they
	      // do not depend on quadrature sum with anything.
	      if( uc == 20 ){
		hPJER->SetBinContent( axis3Bin, uncertainty * 100 );
	      } else if( uc == 21 ){
		hPANG->SetBinContent( axis3Bin, uncertainty * 100 );
	      } else if( uc == 22 ){
		hPFIT->SetBinContent( axis3Bin, uncertainty * 100 );
	      } else if( uc == 23 ){
		hPUNF->SetBinContent( axis3Bin, uncertainty * 100 );
	      } 
	      continue;
	    } else if( sign > 0 ){
	      // add onto total positive uncertainty
	      uncertFinYP += uncertaintySq ;
	      // here go positive uncertainties
	      // first, check if it is JES, then if HIJES
	      if( uc >= 1 && uc <= 18 ){
		uncertFinYPJES += uncertaintySq;
	      } else if( uc == 19 ){
		hPHIJES->SetBinContent( axis3Bin, uncertainty * 100 );
	      } 
	    } else if( sign < 0 ){
	      // add onto total negative uncertainty
	      uncertFinYN += uncertaintySq;
	      // here go negative uncertainties
	      // first, check if it is JES, then if HIJES
	      if( uc >= -18 && uc <= -1 ){
		uncertFinYNJES += uncertaintySq;
	      } else if( uc == -19 ){
		hNHIJES->SetBinContent( axis3Bin, uncertainty * 100 );
	      } 
	    }
	  } // end loop over uncertainties

	  // clean up, or there are memory problems
	  // because each file writes same histo into
	  // same memory address. so eventually you start
	  // deleting deleted stuff if you dont do it now
	  for( auto& h : vHunc ){ delete h; }

	  uncertFinYP = uncertFinYP >= 0 ? std::sqrt( uncertFinYP ) : 0.0;
	  uncertFinYN = uncertFinYN >= 0 ? std::sqrt( uncertFinYN ) : 0.0;

	  hPTOT->SetBinContent( axis3Bin, uncertFinYP * 100 );
	  hNTOT->SetBinContent( axis3Bin, uncertFinYN * 100 );
	  
	  uncertFinYPJES = uncertFinYPJES >= 0 ? std::sqrt( uncertFinYPJES ) : 0.0;
	  uncertFinYNJES = uncertFinYNJES >= 0 ? std::sqrt( uncertFinYNJES ) : 0.0;

	  hPJES->SetBinContent( axis3Bin, uncertFinYPJES * 100 );
	  hNJES->SetBinContent( axis3Bin, uncertFinYNJES * 100 );
	  
	  std::cout << hNominalName << " " << uncertFinYP << " " << uncertFinYN << std::endl;
	  
	  pY .push_back( yNominal );
	  eYP.push_back( yNominal * uncertFinYP );
	  eYN.push_back( yNominal * uncertFinYN );
	} // end loop over axis3
	
	std::string gSystematicsName = "g_" + allSystematicsName + "_" + hTag;

	TGraphAsymmErrors* gNominal     = new TGraphAsymmErrors( hNominal );
	TGraphAsymmErrors* gSystematics = new TGraphAsymmErrors
	  ( nAxis3Bins, &(pX[0]), &(pY[0]), &(eX[0]), &(eX[0]), &(eYN[0]), &(eYP[0]) );
	
	gNominal    ->SetName( gNominalName.c_str() );
	gSystematics->SetName( gSystematicsName.c_str() );

	vG.push_back( gNominal );
	vG.push_back( gSystematics );
	
	styleTool->SetHStyle( gNominal    , style );
	styleTool->SetHStyle( gSystematics, style );

	// Get MC
	TH1* hMCp8 =
	  static_cast< TH1D* >( fMCp8->Get( Form("h_%s_%s", hMCname.c_str(), hTag.c_str() ) ) );
	TGraphAsymmErrors* gMCp8 = new TGraphAsymmErrors( hMCp8 );
	styleTool->SetHStyle( gMCp8, style + 4 );
	
	// add some displacement along x
	double* xDef      = gNominal->GetX();
	double* eXDefLow  = gNominal->GetEXlow();
	double* eXDefHigh = gNominal->GetEXhigh();
	double* xSys      = gSystematics->GetX();
	double* xDefMC      = gMCp8->GetX();
	double* eXDefLowMC  = gMCp8->GetEXlow();
	double* eXDefHighMC = gMCp8->GetEXhigh();
	for( int iX = 0; iX < nAxis3Bins; iX++ ){
	  *(      xDef + iX ) -= style * pDx1;
	  *(  eXDefLow + iX ) = 0;
	  *( eXDefHigh + iX ) = 0;
	  *(      xDefMC + iX ) -= style * pDx1;
	  *(  eXDefLowMC + iX ) = 0;
	  *( eXDefHighMC + iX ) = 0;
	  *(      xSys + iX ) -= style * pDx1;
	}
	
	gSystematics->SetTitle("");
       	gSystematics->GetXaxis()->SetRangeUser( x0, x1 );
       	gSystematics->GetYaxis()->SetRangeUser( y0, y1 );
	
	gNominal->SetTitle("");
       	gNominal->GetXaxis()->SetRangeUser( x0, x1 );
       	gNominal->GetYaxis()->SetRangeUser( y0, y1 );

	legData.AddEntry
	  ( gNominal, Form( "Data, %s", anaTool->GetLabel
			    ( axis2Low, axis2Up, m_dPP->GetAxisLabel(2) ).c_str() ), "lp" );	
	legMC.AddEntry( gMCp8, "MC", "lp" );
	
	// switch back to pad1
	pad1All.cd();
	// draw systematics first
	gSystematics->Draw("2");
	gMCp8->Draw("p");
	gNominal->Draw("p");

	// get ratio and draw
	pad2All.cd();
	TGraphAsymmErrors* gR = anaTool->GetRatioWithSys( gNominal, gSystematics, gMCp8 );
	gR->SetName( Form( "h_data_MC_ratio_%s", hTag.c_str()));
	styleTool->SetHStyle( gR, style + 4 );
	gR->Draw("p");
	line.Draw();
	lineN25.Draw();
	lineP25.Draw();
	
	//+++++++++++++++++++++++++++++++++++++++++++
	// Draw Uncertainties
	//+++++++++++++++++++++++++++++++++++++++++++
	TCanvas cSyst( "cSyst", "cSyst", 800, 600 );
	TLegend legSyst( 0.2, 0.2, 0.9, 0.3 );
	styleTool->SetLegendStyle( &legSyst, 0.90 );
	legSyst.SetNColumns( 4 );
	
	legSyst.AddEntry( hPTOT, "Total", "l" );
	legSyst.AddEntry( hPJES, "JES", "l" );
	legSyst.AddEntry( hPHIJES, "HI JES", "l" );
	legSyst.AddEntry( hPJER, "JER", "l" );
	legSyst.AddEntry( hPANG, "JAR", "l" );
	legSyst.AddEntry( hPUNF, "Unfolding", "l" );
	if( !isYield ){
	  legSyst.AddEntry( hPFIT, "Fitting", "l" );
	}


	TH1* hNJER = static_cast< TH1D* >
	  ( hPJER->Clone( hNJERName.c_str() ) );
	TH1* hNANG = static_cast< TH1D* >
	  ( hPANG->Clone( hNANGName.c_str() ) );
	TH1* hNUNF = static_cast< TH1D* >
	  ( hPUNF->Clone( hNUNFName.c_str() ) );
	TH1* hNFIT = static_cast< TH1D* >
	  ( hPFIT->Clone( hNFITName.c_str() ) );
	
	int maxBin = hPTOT->GetMaximumBin();
	int minBin = hNTOT->GetMaximumBin();

	double newSystMax = systMax;
	double newSystMin = systMin;
	
	if( hPTOT->GetBinContent( maxBin ) > 20 &&
	    hPTOT->GetBinContent( maxBin ) < 25 ){
	  newSystMax += 10;
	} else if( hPTOT->GetBinContent( maxBin ) > 25 &&
		   hPTOT->GetBinContent( maxBin ) < 30 ){
	  newSystMax += 15;
	} else if( hPTOT->GetBinContent( maxBin ) > 30 &&
		   hPTOT->GetBinContent( maxBin ) < 50 ){
	  newSystMax += 25;
	}

	if( hNTOT->GetBinContent( minBin ) > 25 &&
	    hNTOT->GetBinContent( minBin ) < 40 ){
	  newSystMin -= 10;
	} else if( hNTOT->GetBinContent( minBin ) > 40 &&
		   hNTOT->GetBinContent( minBin ) < 50 ){
	  newSystMin -= 15;
	}
	
	hPJES->SetYTitle( yTitleUncert.c_str() );
	hPJES->SetNdivisions( 505 );
	hPJES->SetMaximum( newSystMax );
	hPJES->SetMinimum( newSystMin );
	hPJES->Draw( "histo same" );
	hNJES->Scale( -1. );
	hNJES->Draw( "histo same" );
        hPHIJES->SetMaximum( newSystMax );
	hPHIJES->SetMinimum( newSystMin );
	hPHIJES->Draw( "histo same" );
	hNHIJES->Scale( -1. );
	hNHIJES->Draw( "histo same" );
	hPANG->SetMaximum( newSystMax );
	hPANG->SetMinimum( newSystMin );
	hPANG->Draw( "histo same" );
	hNANG->Scale( -1. );
	hNANG->Draw( "histo same" );
	hPJER->SetMaximum( newSystMax );
	hPJER->SetMinimum( newSystMin );
	hPJER->Draw( "histo same" );
	hNJER->Scale( -1. );
	hNJER->Draw( "histo same" );
	hPUNF->SetMaximum( newSystMax );
	hPUNF->SetMinimum( newSystMin );
	hPUNF->Draw( "histo same" );
	hNUNF->Scale( -1. );
	hNUNF->Draw( "histo same" );
	if( !isYield ){
	  hPFIT->SetMaximum( newSystMax );
	  hPFIT->SetMinimum( newSystMin );
	  hPFIT->Draw( "histo same" );
	  hNFIT->Scale( -1. );
	  hNFIT->Draw( "histo same" );
	}
	// Draw Last
	hPTOT->SetMaximum( newSystMax );
	hPTOT->SetMinimum( newSystMin );
	hPTOT->Draw( "histo same" );
	hNTOT->Scale( -1. );
	hNTOT->Draw( "histo same" );
		
	// Draw the final canvas with all of the graphs.
	legSyst.Draw();

	drawTool->DrawLeftLatex
	  ( 0.46, 0.866, anaTool->GetLabel( axis0Low, axis0Up, m_dPP->GetAxisLabel(0) ), 0.85 );
	drawTool->DrawLeftLatex
	  ( 0.18 , 0.87 , anaTool->GetLabel( axis1Low, axis1Up, m_dPP->GetAxisLabel(1) ), 0.85 );
	drawTool->DrawLeftLatex
	  ( 0.18 , 0.795, anaTool->GetLabel( axis2Low, axis2Up, m_dPP->GetAxisLabel(2) ), 0.85 );
	DrawAtlasRight( 0, 0, 0.9 );
	
	std::string hNameSystFinal =
	"h_" + name + "_" + m_systematicsName + "_" + m_sFinal + "_" + hTag;
 
	SaveAsAll( cSyst, hNameSystFinal );

	// INCREMENT STYLE
	style++;
	// !!! Important !!!
	// Switch back to cAll first, becase a canvas for
	// each uncertainty was created for each axis2Bin.
	cAll.cd();
      } // end loop over axis2

      pad1All.cd();
      // Draw the final canvas with all of the graphs.
      legData.Draw();
      legMC.Draw();
      
      DrawTopLeftLabels
	( m_dPP, axis0Low, axis0Up, axis1Low, axis1Up, 0, 0, 0, 0);

      DrawAtlasRight();

      std::string hNameFinal =
	"h_" + name + "_" + m_sFinal + "_" + hTagC;

      SaveAsAll( cAll, hNameFinal );
    } // end loop over axis1     
  } // end loop over axis0

  for( auto& h : vHdef  ){ h->Write(); delete h; }
  for( auto& h : vHsyst ){ h->Write(); delete h; }
  for( auto& g : vG     ){ g->Write(); delete g; }
  for( auto& p : mFinUC ){ p.second->Close();    }
}

void DiJetAnalysisData::MakeFinalPlotsTogether( TFile* fOut, const std::string& name ){

  std::vector< int > v_uc;
  std::map< int, TFile* > mFinUC;

  // dont need those files open
  for( auto& p : mFinUC ){ p.second->Close(); }

  // Get Infos, same as in DiJetAnalysis::MakeDphiTogether(..)
  std::string name_def, name_syst;
  std::string label_a, label_b ;
  std::string fName_a, fName_b;

  GetInfoTogether( name_def, name_syst, label_a, label_b, fName_a, fName_b, 2 );

  TFile* fInA = TFile::Open( fName_a.c_str() );
  TFile* fInB = TFile::Open( fName_b.c_str() );
  
  // see if we are dealing with widths of yields.
  bool isYield = name.find( m_yieldName ) != std::string::npos ? true : false;
  
  // get map of sys factor to TFile*
  GetListOfSystUncert( v_uc, mFinUC );
  
  // change back to the fOut, for writing.
  fOut->cd();

  std::vector< TH1* > vHdef;
  std::vector< TH1* > vHsyst;
  std::vector< TH1* > vR;
  std::vector< TGraphAsymmErrors* > vG;
  std::vector< TGraphAsymmErrors* > vGfinal;
  
  std::string allUnfoldedName    =
    m_dPhiName + "_" + m_unfoldedName    + "_" + name + "_" + m_allName;
  std::string allSystematicsName =
    m_dPhiName + "_" + m_systematicsName + "_" + name + "_" + m_allName;

  TAxis* axis0 = m_dPP->GetTAxis(0); int nAxis0Bins = axis0->GetNbins();
  TAxis* axis1 = m_dPP->GetTAxis(1); int nAxis1Bins = axis1->GetNbins();
  TAxis* axis2 = m_dPP->GetTAxis(2); int nAxis2Bins = axis2->GetNbins();
  TAxis* axis3 = m_dPP->GetTAxis(3); int nAxis3Bins = axis3->GetNbins();

  // for canvas, since its tgraphasymmerrors.
  double x0 = axis3->GetXmin();
  double x1 = axis3->GetXmax(); 

  double y0 = isYield ? m_dPhiYieldMin : 0.0;
  double y1 = isYield ? m_dPhiYieldMax : 0.9;

  double pDx1 = 0.20;
  double pDx2 = 0.23;

  // lines to be drawn along axis3. this is
  // x-axis that widths are plotted as function of 
  TLine line( x0, 1, x1, 1 );
  line.SetLineWidth( 2 );

  TLine lineP25( x0, 1.25, x1, 1.25 );
  lineP25.SetLineStyle( 2  );
  lineP25.SetLineColor( 12 );
  lineP25.SetLineWidth( 2  );
	  
  TLine lineN25( x0, 0.75, x1, 0.75 );
  lineN25.SetLineStyle( 2  );
  lineN25.SetLineColor( 12 );
  lineN25.SetLineWidth( 2  );

  TLine lineP50( x0, 1.50, x1, 1.50 );
  lineP50.SetLineStyle( 2  );
  lineP50.SetLineColor( 12 );
  lineP50.SetLineWidth( 1  );
	  
  TLine lineN50( x0, 0.50, x1, 0.50 );
  lineN50.SetLineStyle( 2  );
  lineN50.SetLineColor( 12 );
  lineN50.SetLineWidth( 1  );

  TLine lineP1S( x0, 1.34, x1, 1.34 );
  lineP1S.SetLineStyle( 2  );
  lineP1S.SetLineColor( 12 );
  lineP1S.SetLineWidth( 1  );
	  
  TLine lineN1S( x0, 0.66, x1, 0.66 );
  lineN1S.SetLineStyle( 2  );
  lineN1S.SetLineColor( 12 );
  lineN1S.SetLineWidth( 1  );

  std::string yTitle = isYield ? m_sYieldTitle : m_sWidthTitle;
  std::string xTitle = m_dPP->GetAxisLabel(3);
  std::string gTitle = ";" + xTitle + ";" + yTitle;

  std::string sRatio       = isYield ? m_sYieldRatioTitle : m_sWidthRatioTitle;
  std::string gTitleRatio  = ";" + xTitle + "; " + sRatio;

  double y0R = 0.5;
  double y1R = 1.6;
  
  TH1* hDef  = new TH1D( "hDef", gTitle.c_str(), nAxis3Bins, x0, x1 );
  styleTool->SetHStyle( hDef, 0 );
  hDef->SetMaximum( y1 );
  hDef->SetMinimum( y0 );
  hDef->GetXaxis()->SetTitle("");
  hDef->GetYaxis()->SetNdivisions( 502 );

  TH1* hDefR = new TH1D( "hDefR", gTitleRatio.c_str(), nAxis3Bins, x0, x1 );
  styleTool->SetHStyleRatio( hDefR, 0 );
  hDefR->SetMaximum( 1.55 );
  hDefR->SetMinimum( 0.45 );
  
  for( int axis0Bin = 1; axis0Bin <= nAxis0Bins; axis0Bin++ ){
    double axis0Low, axis0Up;
    anaTool->GetBinRange
      ( axis0, axis0Bin, axis0Bin, axis0Low, axis0Up );
    
    for( int axis1Bin = 1; axis1Bin <= nAxis1Bins; axis1Bin++ ){
      if( !m_dPP->CorrectPhaseSpace
	  ( std::vector<int>{ axis0Bin, axis1Bin, 0, 0 } ) )
	{ continue; }

      if( axis1Bin == 3 && isYield ){
	hDef->SetMaximum( 2E-3 );
	hDef->SetMinimum( 5E-7 );
      }
      
      double axis1Low, axis1Up;
      anaTool->GetBinRange
	( axis1, axis1Bin, axis1Bin, axis1Low, axis1Up );

      // tag for widths canvas
      std::string hTagC =
	Form ("%s_%s",
	      anaTool->GetName( axis0Low, axis0Up, m_dPP->GetAxisName(0) ).c_str(),
	      anaTool->GetName( axis1Low, axis1Up, m_dPP->GetAxisName(1) ).c_str() );
      
      // Now Draw everything.
      TCanvas cF( "cF", "cF", 800, 700 );
      TPad pad1F( "pad1", "", 0.0, 0.35, 1.0, 1.0 );
      pad1F.SetBottomMargin(0.023);
      pad1F.Draw();
      TPad pad2F( "pad2F", "", 0.0, 0.0, 1.0, 0.34 );
      pad2F.SetTopMargin(0.05);
      pad2F.SetBottomMargin(0.27);
      pad2F.Draw();
     
      pad2F.cd();
      hDefR->Draw();
      pad1F.cd();
      if( isYield ){
	hDef ->GetYaxis()->SetTitleOffset( 1.3 );
      }
      hDef->Draw();
      pad2F.cd();

      // this is bs. works though. vary the last factor 
      double deltaYleg = ( axis1Bin - 1 ) * 0.04;

      double legSX0pPb = isYield ? 0.33 : 0.18;
      double legSX0pp  = isYield ? 0.45 : 0.30;
      double legSX1pPb = isYield ? 0.52 : 0.37;
      double legSX1pp  = isYield ? 0.73 : 0.47;

      double legSY0 = isYield ? 0.025 : 0.025;
      double legSY1 = isYield ? 0.13 + deltaYleg : 0.13 + deltaYleg;
      
      TLegend legSpPb( legSX0pPb, legSY0, legSX1pPb, legSY1 );
      styleTool->SetLegendStyle( &legSpPb, 0.85 );
      TLegend legSpp ( legSX0pp,  legSY0,  legSX1pp, legSY1 );
      styleTool->SetLegendStyle( &legSpp, 0.85 );
      
      int style = 0;
      
      // ---- loop over axis2 ----
      for( int axis2Bin = 1; axis2Bin <= nAxis2Bins; axis2Bin++ ){
	if( !m_dPP->CorrectPhaseSpace
	    ( std::vector<int>{ axis0Bin, axis1Bin, axis2Bin, 0 } ) )
	  { continue; }

	double axis2Low , axis2Up;
	anaTool->GetBinRange
	  ( axis2, axis2Bin, axis2Bin, axis2Low, axis2Up );

	std::string hTag =
	  Form( "%s_%s_%s",
	        anaTool->GetName( axis0Low, axis0Up, m_dPP->GetAxisName(0) ).c_str(),
		anaTool->GetName( axis1Low, axis1Up, m_dPP->GetAxisName(1) ).c_str(),
		anaTool->GetName( axis2Low, axis2Up, m_dPP->GetAxisName(2) ).c_str() ); 

	TCanvas cAll( "cAll", "cAll", 800, 600 );
	styleTool->SetCStyleGraph( cAll, x0, y0, x1, y1, gTitle.c_str() );

	TLegend legAll( 0.55, 0.25, 0.65, 0.38 );
	styleTool->SetLegendStyle( &legAll, 1.0 );

	// ---------------------------------------
	//      Nominal Graphs and Histograms
	// ---------------------------------------
	std::string hNominalName = "h_" + allUnfoldedName + "_" + hTag;
	std::string gNominalName = "g_" + allUnfoldedName + "_" + hTag;

	TH1* hNominalA = static_cast<TH1D*>( fInA->Get( hNominalName.c_str() ) );
	TH1* hNominalB = static_cast<TH1D*>( fInB->Get( hNominalName.c_str() ) );
	hNominalA->SetName( Form( "%s_A", hNominalName.c_str() ) );
	hNominalB->SetName( Form( "%s_B", hNominalName.c_str() ) );
	vHdef.push_back( hNominalA );
	vHdef.push_back( hNominalB );
	
	TGraphAsymmErrors* gNominalA =
	  static_cast< TGraphAsymmErrors* >( fInA->Get( gNominalName.c_str() ) );
	TGraphAsymmErrors* gNominalB =
	  static_cast< TGraphAsymmErrors* >( fInB->Get( gNominalName.c_str() ) );
	gNominalA->SetName( Form( "%s_A", gNominalName.c_str() ) );
	gNominalB->SetName( Form( "%s_B", gNominalName.c_str() ) );
	styleTool->SetHStyle( gNominalA, 0 );
	styleTool->SetHStyle( gNominalB, 1 );
	vG.push_back( gNominalA );
	vG.push_back( gNominalB );

	legAll.AddEntry( gNominalA, label_a.c_str(), "lp" );
	legAll.AddEntry( gNominalB, label_b.c_str(), "lp" );
	
	// ---------------------------------------
	//    Systematics Graphs and Histograms
	// ---------------------------------------
	std::string hSystematicsName = "h_" + allSystematicsName + "_" + hTag;
	std::string gSystematicsName = "g_" + allSystematicsName + "_" + hTag;

	TGraphAsymmErrors* gSystematicsA =
	  static_cast< TGraphAsymmErrors* >( fInA->Get( gSystematicsName.c_str() ) );
	TGraphAsymmErrors* gSystematicsB =
	  static_cast< TGraphAsymmErrors* >( fInB->Get( gSystematicsName.c_str() ) );
	gSystematicsA->SetName( Form( "%s_A", gSystematicsName.c_str() ) );
	gSystematicsB->SetName( Form( "%s_B", gSystematicsName.c_str() ) );
	styleTool->SetHStyle( gSystematicsA, 0 );
	styleTool->SetHStyle( gSystematicsB, 1 );
	vG.push_back( gSystematicsA );
	vG.push_back( gSystematicsB );

	// offset the graph points
	// and clean some "bad" points
	for( int i = 0; i < nAxis3Bins; i++ ){
	  double x0, y0;
	  
	  x0 = hNominalB->GetBinCenter ( i  + 1 );
	  y0 = hNominalB->GetBinContent( i  + 1 );

	  double x, y;
	  gNominalB->GetPoint( i, x, y );
	  
	  std::cout << " dddddddddddd " << i << " " << x0 << " " << x << std::endl; 
	  gNominalB    ->SetPoint( i, x0 + pDx2, y0 );
	  gSystematicsB->SetPoint( i, x0 + pDx2, y0 );
	  x0 = hNominalA->GetBinCenter ( i  + 1 );
	  y0 = hNominalA->GetBinContent( i  + 1 );
	  gNominalA    ->SetPoint( i, x0 + pDx2, y0 );
	  gSystematicsA->SetPoint( i, x0 + pDx2, y0 );
	  
	}
	
	// now draw everything
	// pad1All.cd();
	// for yields, use log scale on yaxis
	if( isYield ){ cAll.SetLogy(); }
	
	gSystematicsA->Draw("2");
	gNominalA->Draw("p");
	gSystematicsB->Draw("2");
	gNominalB->Draw("p");
	
	legAll.Draw();

	DrawTopLeftLabels
	  ( m_dPP, axis0Low, axis0Up, axis1Low, axis1Up,
	    axis2Low, axis2Up );

	DrawAtlasRightBoth( 0, 0, 1.0, true );
	
	std::string hNameFinal =
	  "h_" + name + "_" + m_sFinal + "_" + hTag;

	SaveAsPdfPng( cAll, hNameFinal, true );
	SaveAsROOT  ( cAll, hNameFinal );

	//-------------------------------------------------------
	//  Make Ratio Plots With Statistical Uncertainties
	//-------------------------------------------------------

	std::string hNominalRName = "h_" + allUnfoldedName + "_" + m_sRatio + "_" + hTag;
	std::string gNominalRName = "g_" + allUnfoldedName + "_" + m_sRatio + "_" + hTag;

	TH1* hNominalR = static_cast< TH1D* >( hNominalA->Clone( hNominalRName.c_str() ) );
	hNominalR->Divide( hNominalB );
	vR.push_back( hNominalR );
	
	std::map< int, TH1* > mHsystTmpA;
	std::map< int, TH1* > mHsystTmpB; 
	for( auto uc : v_uc ){

	  // skip uc = 0 (default)
	  if( !uc ){ continue; }

	  std::string tmpUncertSuffix = uc > 0 ? Form("P%d", uc) : Form("N%d", -1 * uc) ;
	  std::string hSystematicName = "h_" + allUnfoldedName + "_" + hTag + "_" + tmpUncertSuffix;

	  TH1* hSystematicA = static_cast<TH1D*>( fInA->Get( hSystematicName.c_str() ) );
	  TH1* hSystematicB = static_cast<TH1D*>( fInB->Get( hSystematicName.c_str() ) );
	  mHsystTmpA[ uc ] = hSystematicA;
	  mHsystTmpB[ uc ] = hSystematicB;
	}

	std::vector< double > pX;
        std::vector< double > eX( nAxis3Bins, pDx1 * 0.5 );
 
	std::vector< double > pY;
	std::vector< double > eYP;
	std::vector< double > eYN;

	for( int axis3Bin = 1; axis3Bin <= nAxis3Bins; axis3Bin++ ){

	  // add uncertainties in quadrature;
	  double uncertFinYP = 0;
	  double uncertFinYN = 0;

	  // for FIT(22) and UNC(23)
	  double uncertUnCorrYP = 0;
	  double uncertUnCorrYN = 0;

	  // nominal ratio
	  double ratioABnominal = hNominalR->GetBinContent( axis3Bin );
	  pX.push_back( hNominalR->GetBinCenter( axis3Bin ) );
	  
	  for( auto uc : v_uc ){
	    
	    // skip uc = 0 (default)
	    if( !uc ){ continue; }
	     
	    // for FIT, UNC, errors uncorrelated, add in quadrature
	    // for pp and pPb, and combine to a total error which
	    // is then added to the other correlated systematics.
	    if( uc == 22 || uc == 23 ){
	      double eYA = mHsystTmpA[ uc ]->GetBinContent( axis3Bin );
	      double yA  = hNominalA->GetBinContent( axis3Bin );
	      double deltaEyA =  eYA - yA;
	      
	      double eYB = mHsystTmpB[ uc ]->GetBinContent( axis3Bin );
	      double yB  = hNominalB->GetBinContent( axis3Bin );
	      double deltaEyB =  eYB - yB;

	      double quadSumAB = std::sqrt( std::pow( deltaEyA, 2 ) +
					    std::pow( deltaEyB, 2 ) );
	      
	      uncertUnCorrYP += std::pow( quadSumAB, 2 );
	      uncertUnCorrYN += std::pow( quadSumAB, 2 );
	    } else {
	      double shiftA  =  mHsystTmpA[ uc ]->GetBinContent( axis3Bin );
	      double shiftB  =  mHsystTmpB[ uc ]->GetBinContent( axis3Bin );
	      
	      double ratioABshifted = shiftA / shiftB;
	      
	      double uncertaintySq =
		std::pow( ( ratioABnominal - ratioABshifted ) / ratioABnominal, 2 );

	      int sign = uc > 0 ? 1 : -1;

	      std::cout << "++++" << uc << " ++++" << axis3Bin << " " << axis1Bin << " "
			<< axis2Bin << " " << sign << " " << ratioABshifted << " "
			<< ratioABnominal << " " << std::sqrt( uncertaintySq ) << std::endl;
	      
	      if( sign > 0 ){
		uncertFinYP += uncertaintySq;
		// for ANG, JER, FIT, UNF, need to symmetereize.
		if( uc >= 20 && uc <= 23 ){
		  uncertFinYN += uncertaintySq;
		}
		continue;
	      } else{
		uncertFinYN += uncertaintySq;
		continue;
	      }
	    }
	  }

	  uncertFinYP = uncertUnCorrYP >= 0 ? uncertFinYP + uncertUnCorrYP : 0;
	  uncertFinYN = uncertUnCorrYN >= 0 ? uncertFinYN + uncertUnCorrYN : 0;
	  
	  // now get the final uncertainties.
	  uncertFinYP = uncertFinYP >= 0 ? std::sqrt( uncertFinYP ) : 0.0;
	  uncertFinYN = uncertFinYN >= 0 ? std::sqrt( uncertFinYN ) : 0.0;

	  std::cout << hNominalRName << " " << uncertFinYP << " " << uncertFinYN << std::endl;
	  
	  pY .push_back( ratioABnominal );
	  eYP.push_back( ratioABnominal * uncertFinYP );
	  eYN.push_back( ratioABnominal * uncertFinYN );
	} // end loop over axis3

	std::string gSystematicsRName = "g_" + allSystematicsName + "_" + m_sRatio + "_" + hTag;

	TGraphAsymmErrors* gNominalR     = new TGraphAsymmErrors( hNominalR );
	TGraphAsymmErrors* gSystematicsR = new TGraphAsymmErrors
	  ( nAxis3Bins, &(pX[0]), &(pY[0]), &(eX[0]), &(eX[0]), &(eYN[0]), &(eYP[0]) );

	gNominalR    ->SetName( gNominalRName.c_str() );
	gSystematicsR->SetName( gSystematicsRName.c_str() );

	// get rid of errors on nominal graph
	double* eXNominalLow  = gNominalR->GetEXlow();
	double* eXNominalHigh = gNominalR->GetEXhigh();
	for( int iX = 0; iX < nAxis3Bins; iX++ ){	  
	  *(  eXNominalLow + iX ) = 0;
	  *( eXNominalHigh + iX ) = 0;
	}

	vGfinal.push_back( gNominalR );
	vGfinal.push_back( gSystematicsR );
	
	styleTool->SetHStyle( gNominalR    , 0 );
	styleTool->SetHStyle( gSystematicsR, 0 );
	
	gSystematicsR->SetTitle("");
	gSystematicsR->GetXaxis()->SetRangeUser( x0, x1 );
	gSystematicsR->GetYaxis()->SetRangeUser( y0R, y1R );
	
	gNominalR->SetTitle("");
	gNominalR->GetXaxis()->SetRangeUser( x0, x1 );
	gNominalR->GetYaxis()->SetRangeUser( y0R, y1R );

	TCanvas cRatio( "cRatio", "cRatio", 800, 600 );
	styleTool->SetCStyleGraph( cRatio, x0, y0R, x1, y1R, gTitleRatio.c_str() );

	// draw systematics first
	gSystematicsR->SetFillStyle(0);
	gSystematicsR->Draw("2");
	gNominalR->Draw("p");

	line.Draw();
	
	lineP25.Draw();
	lineN25.Draw();

	drawTool->DrawLeftLatex
	  ( 0.18, 0.86, anaTool->GetLabel( axis1Low, axis1Up, m_dPP->GetAxisLabel(1) ), 1.0 );
	drawTool->DrawLeftLatex
	  ( 0.18, 0.78, anaTool->GetLabel( axis2Low, axis2Up, m_dPP->GetAxisLabel(2) ), 1.0 );
	drawTool->DrawLeftLatex
	  ( 0.18, 0.275, anaTool->GetLabel( axis0Low, axis0Up, m_dPP->GetAxisLabel(0) ), 1.0 );
	
	DrawAtlasRightBoth( 0, 0, 1.0, true );
	  
	std::string hNameFinalR =
	  "h_" + name + "_" + m_sRatio + "_" +m_sFinal + "_" + hTag;

	SaveAsPdfPng( cRatio, hNameFinalR, true );
	SaveAsROOT  ( cRatio, hNameFinalR );

	//-------------------------------------------
	// now draw everything together
	// for yields, use log scale on yaxis
	// Also, change the xshift on plots first
	//-------------------------------------------
	// offset the graph points
	for( int i = 0; i < gSystematicsB->GetN(); i++ ){
	  double x0, y0;
	  x0 = hNominalB->GetBinCenter ( i  + 1 );
	  y0 = hNominalB->GetBinContent( i  + 1 );
	  gNominalB    ->SetPoint( i, x0 + ( axis2Bin - 2 ) * pDx2, y0 );
	  gSystematicsB->SetPoint( i, x0 + ( axis2Bin - 2 ) * pDx2, y0 );
	  x0 = hNominalA->GetBinCenter ( i  + 1 );
	  y0 = hNominalA->GetBinContent( i  + 1 );
	  gNominalA    ->SetPoint( i, x0 + ( axis2Bin - 2 ) * pDx2, y0 );
	  gSystematicsA->SetPoint( i, x0 + ( axis2Bin - 2 ) * pDx2, y0 );
	  x0 = hNominalR->GetBinCenter ( i  + 1 );
	  y0 = hNominalR->GetBinContent( i  + 1 );
	  gNominalR    ->SetPoint( i, x0 + ( axis2Bin - 2 ) * pDx2, y0 );
	  gSystematicsR->SetPoint( i, x0 + ( axis2Bin - 2 ) * pDx2, y0 );


	  // clean some bad points
	  if( !isYield && axis1Bin == 3 && i == 0 &&
	      ( axis2Bin == 1 || axis2Bin == 2 ) ){

	    double x, y;
	    gNominalA->GetPoint( i, x, y );
	    std::cout << " lllllllllll " << gNominalA->GetName()
		      << " " << i << " " << x << std::endl;

	    gNominalA    ->SetPoint( i, -10, y0 );
	    gSystematicsA->SetPoint( i, -10, y0 );
	    gNominalB    ->SetPoint( i, -10, y0 );
	    gSystematicsB->SetPoint( i, -10, y0 );
	    gNominalR    ->SetPoint( i, -10, y0 );
	    gSystematicsR->SetPoint( i, -10, y0 );
	    continue;
	  }
	}

	pad1F.cd();
	
	if( isYield ){ pad1F.SetLogy(); }

	styleTool->HideAxis( hDef, "x" );

	legSpPb.AddEntry( gNominalA, label_a.c_str() , "lp" );
	legSpp .AddEntry( gNominalB, Form("%s, %s", label_b.c_str(), anaTool->GetLabel
					  ( axis2Low, axis2Up, m_dPP->GetAxisLabel(2) ).c_str() ), "lp" );

	styleTool->SetHStyle( gSystematicsA, style );
	styleTool->SetHStyle( gSystematicsB, style + 4 );
	styleTool->SetHStyle( gNominalA, style );
	styleTool->SetHStyle( gNominalB, style + 4 );
	
	gSystematicsA->SetFillStyle(0);
	gSystematicsB->Draw("2");
	gSystematicsA->Draw("2");
	gNominalB->Draw("p");
	gNominalA->Draw("p");

	pad2F.cd();
	
	styleTool->SetHStyle( gSystematicsR, style );
	styleTool->SetHStyle( gNominalR, style );
	gSystematicsR->SetFillStyle(0);
	gSystematicsR->Draw("2");
	gNominalR->Draw("p");

	style++;
      } // end loop over axis2
      pad1F.cd();

      DrawTopLeftLabels
	( m_dPP, axis0Low, axis0Up, axis1Low, axis1Up );
	
      
      DrawAtlasRightBoth( 0, 0, 1.0, true );

      legSpPb.Draw();
      legSpp .Draw();
	
      pad2F.cd();

      line.Draw();
      lineP25.Draw();
      lineN25.Draw();
    	
      std::string hNameFinalC =
	"h_" + name + "_" + m_sFinal + "_" + hTagC;
    
      SaveAsPdfPng( cF, hNameFinalC, true );
      SaveAsROOT  ( cF, hNameFinalC );
    } // end loop over axis1
  } // end loop over axis0

  for( auto& g : vGfinal ){ g->Write(); delete g; }

  for( auto& h : vHdef  ){ delete h; }
  for( auto& h : vHsyst ){ delete h; }
  for( auto& r : vR     ){ delete r; }
  for( auto& g : vG     ){ delete g; }
}

// CLEAN THIS UP! ITS NOT WORTH THE HEADACHE NOW
// BUT STILL, NOT GOOD (03.30.18)
void DiJetAnalysisData::CompareCfactorsRBnRB( TFile* fOut ){

  std::vector< TH1* > vC;
  std::vector< TH1* > vR;
  
  TFile* fRB  =
    TFile::Open( m_fNamePhysUF.c_str() );
  TFile* fDF = m_is_pPb ?
    TFile::Open( "data/output_pPb_data_nr/myOut_pPb_data_phys_UF_0.root" ): 
    TFile::Open( "data/output_pp_data_nr/myOut_pp_data_phys_UF_0.root"   );
  
  fOut->cd();
  
  TAxis* axis0 = m_dPP->GetTAxis(0); int nAxis0Bins = axis0->GetNbins();
  TAxis* axis1 = m_dPP->GetTAxis(1); int nAxis1Bins = axis1->GetNbins();
  TAxis* axis2 = m_dPP->GetTAxis(2); int nAxis2Bins = axis2->GetNbins();
  TAxis* axis3 = m_dPP->GetTAxis(3); int nAxis3Bins = axis3->GetNbins();

  // lines to be drawn along axis3. this is
  // x-axis that widths are plotted as function of 
  double dPhiXmin = m_dPhiZoomLow; double dPhiXmax = m_dPhiZoomHigh;
  
  TLine lineDphi( dPhiXmin, 1, dPhiXmax, 1 );
  lineDphi.SetLineWidth( 2 );

  TLine lineDphiP05( dPhiXmin, 1.05, dPhiXmax, 1.05 );
  lineDphiP05.SetLineStyle( 2  );
  lineDphiP05.SetLineColor( 12 );
  lineDphiP05.SetLineWidth( 1  );
	  
  TLine lineDphiN05( dPhiXmin, 0.95, dPhiXmax, 0.95 );
  lineDphiN05.SetLineStyle( 2  );
  lineDphiN05.SetLineColor( 12 );
  lineDphiN05.SetLineWidth( 1  );

  TLine lineDphiP25( dPhiXmin, 1.25, dPhiXmax, 1.25 );
  lineDphiP25.SetLineStyle( 2  );
  lineDphiP25.SetLineColor( 12 );
  lineDphiP25.SetLineWidth( 2  );
	  
  TLine lineDphiN25( dPhiXmin, 0.75, dPhiXmax, 0.75 );
  lineDphiN25.SetLineStyle( 2  );
  lineDphiN25.SetLineColor( 12 );
  lineDphiN25.SetLineWidth( 2  );

  double xMin = m_varYstarBinning.front();
  double xMax = m_varYstarBinning.back();
  
  TLine line( xMin, 1, xMax, 1 );
  line.SetLineWidth( 2 );

  TLine lineP05( xMin, 1.05, xMax, 1.05 );
  lineP05.SetLineStyle( 2  );
  lineP05.SetLineColor( 12 );
  lineP05.SetLineWidth( 1  );
	  
  TLine lineN05( xMin, 0.95, xMax, 0.95 );
  lineN05.SetLineStyle( 2  );
  lineN05.SetLineColor( 12 );
  lineN05.SetLineWidth( 1  );

  TLine lineP25( xMin, 1.25, xMax, 1.25 );
  lineP25.SetLineStyle( 2  );
  lineP25.SetLineColor( 12 );
  lineP25.SetLineWidth( 2  );
	  
  TLine lineN25( xMin, 0.75, xMax, 0.75 );
  lineN25.SetLineStyle( 2  );
  lineN25.SetLineColor( 12 );
  lineN25.SetLineWidth( 2  );

  
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

      for( int axis2Bin = 1; axis2Bin <= nAxis2Bins; axis2Bin++ ){

	// check we are in correct ystar and pt bins
	if( !m_dPP->CorrectPhaseSpace
	    ( std::vector<int>{ axis0Bin, axis1Bin, axis2Bin, 0 } ) )
	  { continue; }

	double axis2Low , axis2Up;
	anaTool->GetBinRange
	  ( axis2, axis2Bin, axis2Bin, axis2Low, axis2Up );

	// ----------- widths -----------
	std::string hTagW =
	  Form( "%s_%s_%s",
	     	anaTool->GetName( axis0Low, axis0Up, m_dPP->GetAxisName(0) ).c_str(),
		anaTool->GetName( axis1Low, axis1Up, m_dPP->GetAxisName(1) ).c_str(),
		anaTool->GetName( axis2Low, axis2Up, m_dPP->GetAxisName(2) ).c_str() ); 
	  
	std::string hNameW =
	  "h_" + m_dPhiUnfoldedName + "_" + m_widthName + "_" + m_allName + "_" + hTagW;

	TCanvas cWidthsCmp( "cWidthsCmp", "cWidthsCmp", 800, 800 );

	TPad padW1("padW1", "", 0.0, 0.35, 1.0, 1.0 );
	padW1.SetBottomMargin(0.0);
	padW1.Draw();
	  
	TPad padW2("padW2", "", 0.0, 0.0, 1.0, 0.34 );
	padW2.SetTopMargin(0.05);
	padW2.SetBottomMargin(0.25);
	padW2.Draw();

	TLegend legWidths( 0.50, 0.1, 0.89, 0.28 );
	styleTool->SetLegendStyle( &legWidths );

	padW1.cd();

	TH1* hDphiWidthsRB = static_cast< TH1D* >( fRB->Get( hNameW.c_str() ) );
	TH1* hDphiWidthsDF = static_cast< TH1D* >( fDF->Get( hNameW.c_str() ) );
	styleTool->SetHStyle( hDphiWidthsRB, 0 );
	styleTool->SetHStyle( hDphiWidthsDF, 5 );
	// hDphiWidthsRB->SetMarkerSize( hDphiWidthsRB->GetMarkerSize() * 1.5 );
	vC.push_back( hDphiWidthsRB );
	vC.push_back( hDphiWidthsDF );

	hDphiWidthsRB->SetMinimum( m_dPhiWidthMin );
	hDphiWidthsDF->SetMaximum( m_dPhiWidthMax );

	hDphiWidthsRB->Draw("epsame X0");
	hDphiWidthsDF->Draw("epsame X0");

	legWidths.AddEntry( hDphiWidthsRB, "Re-Binned" );
	legWidths.AddEntry( hDphiWidthsDF, "Default"   );
		
	legWidths.Draw();

	DrawTopLeftLabels
	  ( m_dPP, axis0Low, axis0Up, axis1Low, axis1Up,
	    axis2Low, axis2Up, 0, 0 );

	DrawAtlasRight();
	  
	padW2.cd();

	std::string hNameWR = hNameW + "_" + m_sRatio;
	TH1* hWR = static_cast< TH1D* >( hDphiWidthsRB->Clone( hNameWR.c_str() ) );
	styleTool->SetHStyle( hWR, 0 );
	vR.push_back( hWR );
	  
	hWR->SetMaximum( 1.5 );
	hWR->SetMinimum( 0.5 );
	hWR->SetYTitle( "Ratio" );

	hWR->Divide( hDphiWidthsDF );

	hWR->Draw("ep X0" );
	  
	line.Draw();
	lineP25.Draw();
	lineN25.Draw();
	  
	SaveAsAll( cWidthsCmp, hNameW );

	
	// ----------- yields -----------
	std::string hTagY =
	  Form( "%s_%s_%s",
	     	anaTool->GetName( axis0Low, axis0Up, m_dPP->GetAxisName(0) ).c_str(),
		anaTool->GetName( axis1Low, axis1Up, m_dPP->GetAxisName(1) ).c_str(),
		anaTool->GetName( axis2Low, axis2Up, m_dPP->GetAxisName(2) ).c_str() ); 
	  
	std::string hNameY =
	  "h_" + m_dPhiUnfoldedName + "_" + m_yieldName + "_" + m_allName + "_" + hTagY;

	TCanvas cYieldsCmp( "cYieldsCmp", "cYieldsCmp", 800, 800 );

	TPad padY1("padY1", "", 0.0, 0.35, 1.0, 1.0 );
	padY1.SetBottomMargin(0.0);
	padY1.Draw();
	  
	TPad padY2("padY2", "", 0.0, 0.0, 1.0, 0.34 );
	padY2.SetTopMargin(0.05);
	padY2.SetBottomMargin(0.25);
	padY2.Draw();

	TLegend legYields( 0.50, 0.1, 0.89, 0.28 );
	styleTool->SetLegendStyle( &legYields );

	padY1.cd();
	padY1.SetLogy();
	
	TH1* hDphiYieldsRB = static_cast< TH1D* >( fRB->Get( hNameY.c_str() ) );
	TH1* hDphiYieldsDF = static_cast< TH1D* >( fDF->Get( hNameY.c_str() ) );
	styleTool->SetHStyle( hDphiYieldsRB, 0 );
	styleTool->SetHStyle( hDphiYieldsDF, 5 );
	// hDphiYieldsRB->SetMarkerSize( hDphiYieldsRB->GetMarkerSize() * 1.5 );
	vC.push_back( hDphiYieldsRB );
	vC.push_back( hDphiYieldsDF );

	hDphiYieldsRB->SetMinimum( m_dPhiYieldMin );
	hDphiYieldsDF->SetMaximum( m_dPhiYieldMax );

	hDphiYieldsRB->Draw("epsame X0");
	hDphiYieldsDF->Draw("epsame X0");

	legYields.AddEntry( hDphiYieldsRB, "Re-Binned" );
	legYields.AddEntry( hDphiYieldsDF, "Default"   );
			    
	  
	legYields.Draw();

	DrawTopLeftLabels
	  ( m_dPP, axis0Low, axis0Up, axis1Low, axis1Up,
	    axis2Low, axis2Up, 0, 0 );

	DrawAtlasRight();
	  
	padY2.cd();

	std::string hNameYR = hNameY + "_" + m_sRatio;
	TH1* hYR = static_cast< TH1D* >( hDphiYieldsRB->Clone( hNameYR.c_str() ) );
	styleTool->SetHStyle( hYR, 0 );
	vR.push_back( hYR );
	  
	hYR->SetMaximum( 1.5 );
	hYR->SetMinimum( 0.5 );
	hYR->SetYTitle( "Ratio" );

	hYR->Divide( hDphiYieldsDF );

	hYR->Draw("ep X0" );
	  
	line.Draw();
	lineP25.Draw();
	lineN25.Draw();
	  
	SaveAsAll( cYieldsCmp, hNameY );

	
	for( int axis3Bin = 1; axis3Bin <= nAxis3Bins; axis3Bin++ ){
	  // check we are in correct ystar and pt bins
	  if( !m_dPP->CorrectPhaseSpace
	      ( std::vector<int>{ axis0Bin, axis1Bin, axis2Bin, axis3Bin } ) )
	    { continue; }

	  double axis3Low , axis3Up;
	  anaTool->GetBinRange
	    ( axis3, axis3Bin, axis3Bin, axis3Low, axis3Up );
	  
	  //-----------------------------------------------------
	  //         CORRECTION FACTORS AND THEIR RATIO
	  //-----------------------------------------------------
	  
	  std::string hTag =
	    anaTool->GetName( axis0Low, axis0Up, m_dPP->GetAxisName(0) ) + "_" + 
	    anaTool->GetName( axis1Low, axis1Up, m_dPP->GetAxisName(1) ) + "_" + 
	    anaTool->GetName( axis2Low, axis2Up, m_dPP->GetAxisName(2) ) + "_" +  
	    anaTool->GetName( axis3Low, axis3Up, m_dPP->GetAxisName(3) ); 
	  
	  std::string hNameCF = "h_" + m_dPhiCfactorsName + "_" + m_allName + "_" + hTag;

	  TH1* hCFRB = static_cast< TH1D* >( fRB->Get( hNameCF.c_str() ) );
	  TH1* hCFDF = static_cast< TH1D* >( fDF->Get( hNameCF.c_str() ) );
	  styleTool->SetHStyle( hCFRB, 0 );
	  styleTool->SetHStyle( hCFDF, 5 );
	  // hCFDF->SetMarkerSize( hDphiWidthsRB->GetMarkerSize() * 1.5 );
	  vC.push_back( hCFRB  );
	  vC.push_back( hCFDF );
	  
	  TLegend legCF( 0.7, 0.2, 0.8, 0.4 );
	  styleTool->SetLegendStyle( &legCF );
	  legCF.AddEntry( hCFRB , "Rebinned" );
	  legCF.AddEntry( hCFDF, "Default"  );
	  
	  TCanvas cCF( "cCF", "cCF", 800, 800 );

	  hCFRB->Draw("ep X0 same" );
	  hCFDF->Draw("ep X0 same" );

	  legCF.Draw();

	  lineDphi.Draw();
	  lineDphiP25.Draw();
	  lineDphiN25.Draw();
	  
	  DrawTopLeftLabels
	    ( m_dPP, axis0Low, axis0Up, axis1Low, axis1Up,
	      axis2Low, axis2Up, axis3Low, axis3Up );
	    
	  DrawAtlasRight();
	  
	  SaveAsAll( cCF, hNameCF );

	  //-----------------------------------------------------
	  //         DELTA PHI DISTRIBUTIONS AND RATIO
	  //-----------------------------------------------------
	  
	  std::string hNameDphi = "h_" + m_dPhiUnfoldedName + "_" + m_allName + "_" + hTag;
	  
	  TH1* hDphiRB = static_cast< TH1D* >( fRB->Get( hNameDphi.c_str() ) );
	  TH1* hDphiDF = static_cast< TH1D* >( fDF->Get( hNameDphi.c_str() ) );
	  styleTool->SetHStyle( hDphiRB, 0 );
	  styleTool->SetHStyle( hDphiDF, 5 );
	  // hDphiDF->SetMarkerSize( hDphiWidthsRB->GetMarkerSize() * 1.5 );
	  vC.push_back( hDphiRB  );
	  vC.push_back( hDphiDF );
	  
	  TLegend legDphi( 0.7, 0.20, 0.8, 0.4 );
	  styleTool->SetLegendStyle( &legDphi );
	  legDphi.AddEntry( hDphiRB, "Weighted"   );
	  legDphi.AddEntry( hDphiDF, "UnWeighted" );
	  
	  TCanvas cDphi( "cDphi", "cDphi", 800, 800 );
	  TPad padDphi1("padDphi1", "", 0.0, 0.35, 1.0, 1.0 );
	  padDphi1.SetBottomMargin(0.0);
	  padDphi1.Draw();
	  
	  TPad padDphi2("padDphi2", "", 0.0, 0.0, 1.0, 0.34 );
	  padDphi2.SetTopMargin(0.05);
	  padDphi2.SetBottomMargin(0.25);
	  padDphi2.Draw();

	  padDphi1.cd();
	  padDphi1.SetLogy();
	  hDphiRB->Draw("ep X0 same" );
	  hDphiDF->Draw("ep X0 same" );

	  legDphi.Draw();

	  DrawTopLeftLabels
	    ( m_dPP, axis0Low, axis0Up, axis1Low, axis1Up,
	      axis2Low, axis2Up, axis3Low, axis3Up );
	    
	  DrawAtlasRight();
	  
	  padDphi2.cd();

	  std::string hNameDphiR = hNameDphi + "_" + m_sRatio;

	  TH1* hRdPhi = static_cast< TH1D* >( hDphiRB->Clone( hNameDphiR.c_str() ) );
	  styleTool->SetHStyleRatio( hRdPhi );
	  vR.push_back( hRdPhi );
	  
	  hRdPhi->SetMaximum( 1.2 );
	  hRdPhi->SetMinimum( 0.8 );

	  hRdPhi->Divide( hDphiDF );
	  hRdPhi->SetYTitle( "Weighted/UnWeighted" );

	  hRdPhi->Draw( "hist p X0" );

	  lineDphi.Draw();
	  lineDphiP05.Draw();
	  lineDphiN05.Draw();
	  
	  SaveAsAll( cDphi, hNameDphi );
	}
      }
    }
  }

  // for( auto& c : vC ){ delete c; }
  // for( auto& r : vR ){ delete r; }

  std::string nameRB = m_dPhiUnfoldedName;
  std::string nameDF = m_dPhiUnfoldedName;

  TH1* h_chi2RB = static_cast< TH1D* >
    ( fRB->Get( Form( "h_dPhiChiS_%s", nameRB.c_str() ) ) );
  styleTool->SetHStyle( h_chi2RB, 0 );
  
  TH1* h_chi2DF = static_cast< TH1D* >
    ( fDF->Get( Form( "h_dPhiChiS_%s", nameDF.c_str() ) ) );
  styleTool->SetHStyle( h_chi2DF, 1 );

  TH1* h_probRB = static_cast< TH1D* >
    ( fRB->Get( Form( "h_dPhiProb_%s", nameRB.c_str() ) ) );
  styleTool->SetHStyle( h_probRB, 0 );
  
  TH1* h_probDF = static_cast< TH1D* >
    ( fDF->Get( Form( "h_dPhiProb_%s", nameDF.c_str() ) ) );
  styleTool->SetHStyle( h_probDF, 1 );

  TCanvas c( "c", "c", 1200, 600 );
  c.Divide( 2, 1 );

  TLegend leg( 0.6, 0.6, 0.7, 0.7 );
  styleTool->SetLegendStyle( &leg );
  
  c.cd(1);
  h_chi2RB->SetMaximum( h_chi2RB->GetMaximum() > h_chi2DF->GetMaximum() ?
			h_chi2RB->GetMaximum() * 1.1 : h_chi2DF->GetMaximum() * 1.1 );
  h_chi2RB->Draw("hist C same");
  h_chi2DF->Draw("hist C same");
  DrawAtlasRightBoth();

  c.cd(2);
  h_probRB->SetMaximum( h_probRB->GetMaximum() > h_probDF->GetMaximum() ?
			h_probRB->GetMaximum() * 1.1 : h_probDF->GetMaximum() * 1.1 );
  h_probRB->Draw("hist C same");
  h_probDF->Draw("hist C same");

  leg.AddEntry( h_probRB, "Re-Binned" );
  leg.AddEntry( h_probDF, "Default"   );

  leg.Draw();
  
  DrawAtlasRightBoth();

  SaveAsAll( c, "h_chi2_prob" );
}

//---------------------------
//        Drawing
//---------------------------

void DiJetAnalysisData::DrawAtlasRight( double x0, double y0, double scale )
{ drawTool->DrawAtlasInternalDataRight( x0, y0, m_is_pPb, scale ); } 

void DiJetAnalysisData::DrawAtlasRightBoth( double x0, double y0, double scale, bool fin ){
  drawTool->DrawAtlasInternal( scale );
  if( fin ){
    drawTool->DrawRightLatex
      ( 0.875, 0.805, drawTool->GetLumipPb(), scale, 1 );
    drawTool->DrawRightLatex
      ( 0.875, 0.735, drawTool->GetLumipp(), scale, 1 );
  } else { 
    drawTool->DrawRightLatex( 0.87, 0.78, "Data" );	  	  
  }
} 
