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

  //=============== Histo Names ==================    
  m_ystarSpectFineRunsName = m_ystarSpectFineName + "_" + m_sRuns;
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

  m_hAllEtaPtMap = CombineSamples( m_vHtriggerEtaPtMap, "etaPtMap" );
  MakeEtaPhiPtMap( m_vHtriggerEtaPtMap , m_vTriggers, "etaPtMap" );

  m_hAllYstarSpect = CombineSamples( m_vHtriggerYstarSpect, m_ystarSpectName );
  MakeSpectra( m_vHtriggerYstarSpect, m_vTriggers, m_ystarSpectName );

  m_hAllYstarSpectFine = CombineSamples( m_vHtriggerYstarSpectFine, m_ystarSpectFineName );
  MakeSpectra( m_vHtriggerYstarSpectFine, m_vTriggers, m_ystarSpectFineName );

  MakeEfficiencies( m_vHtriggerEtaSpectSim, m_vHtriggerEtaSpectDenom, m_effName );
  
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
		  ";#eta;#it{p}_{T} [GeV];",
		  m_nVarYstarBins, 0, 1,
		  m_nPtSpectBins,m_ptSpectMin, m_ptSpectMax ) ) ;
    m_vHtriggerYstarSpectFine.back()->GetXaxis()->
      Set( m_nVarYstarBins, &( m_varYstarBinning[0] ) );
    AddHistogram( m_vHtriggerYstarSpectFine.back() );

    m_vHtriggerYstarSpectFineRuns.push_back
      ( new TH3D( Form("h_%s_%s", m_ystarSpectFineRunsName.c_str(), trigger.c_str() ), 
		  ";#eta;#it{p}_{T} [GeV];",
		  m_nVarYstarBins, 0, 1,
		  m_nPtSpectBins,m_ptSpectMin, m_ptSpectMax,
		  m_nRuns, 0, m_nRuns ) ) ;
    m_vHtriggerYstarSpectFineRuns.back()->GetXaxis()->
      Set( m_nVarYstarBins, &( m_varYstarBinning[0] ) );
    AddHistogram( m_vHtriggerYstarSpectFineRuns.back() );

    // ----- efficiencies ----
    m_vHtriggerEtaSpectSim.push_back
      ( new TH2D ( Form("h_%sSim_%s", m_etaSpectName.c_str(), trigger.c_str() ), 
		   ";#eta;#it{p}_{T} [GeV]",
		   m_nVarFwdEtaBins, 0, 1,
		   m_nPtSpectBins, m_ptSpectMin, m_ptSpectMax ) ) ;
    m_vHtriggerEtaSpectSim.back()->GetXaxis()->
      Set( m_nVarFwdEtaBins, &( m_varFwdEtaBinning[0] ) );
    AddHistogram( m_vHtriggerEtaSpectSim.back() );

    m_vHtriggerEtaSpectDenom.push_back
      ( new TH2D( Form("h_%sDenom_%s", m_etaSpectName.c_str(), trigger.c_str() ), 
		  ";#eta;#it{p}_{T} [GeV];",
		  m_nVarFwdEtaBins, 0, 1,
		  m_nPtSpectBins, m_ptSpectMin, m_ptSpectMax ) ) ;
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
    
    // SPECTRA AND ETA PHI PT MAPS
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

	m_vHtriggerYstarSpect[iG]->
	  Fill( jetYstar1, jetPt1 );
	m_vHtriggerYstarSpectFine[iG]->
	  Fill( jetYstar1, jetPt1 );
	m_vHtriggerYstarSpectFineRuns[iG]->
	  Fill( jetYstar1, jetPt1, m_mRunBin[runNumber] - 1 );	
	if( !m_is_pPb ){
	  m_vHtriggerYstarSpect[iG]->
	    Fill( -jetYstar1, jetPt1 );
	  m_vHtriggerYstarSpectFine[iG]->
	    Fill( -jetYstar1, jetPt1 );
	  m_vHtriggerYstarSpectFineRuns[iG]->
	    Fill( -jetYstar1, jetPt1, m_mRunBin[runNumber] - 1 );
	}
      }

      // loop over jets 
      for( auto& jet : vR_jets ){

	double jetPt = jet.Pt()/1000.;

	// check if jet is in trigger range
	if( !JetInTrigRange( &jet, iG ) ) continue;
	
	// ETA-PHI
	double jetEta   = jet.Eta();
	double jetPhi   = jet.Phi();
	// double jetYstar = GetYstar( jet );
	
	m_vHtriggerEtaPhiMap[iG]->Fill( jetEta, jetPhi );
	m_vHtriggerEtaPtMap [iG]->Fill( jetEta, jetPt  ); 
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
					const std::string& name ){

  if( !vSampleHin.size() ){ return NULL; }

  TH3* h_res = static_cast< TH3D* >
    ( vSampleHin[0]->Clone( Form("h_%s_%s", name.c_str() , m_allName.c_str() ) ) );
  h_res->Reset();
    
  for( uint iG = 0; iG < vSampleHin.size(); iG++ ){
    double scale = m_vTriggersPrescale[iG];
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
	( fIn-> Get( Form("h_etaPhiMap_%s",
			    trigger.c_str() ))));
    m_vHtriggerEtaPhiMap.back()->SetDirectory(0);
    
    m_vHtriggerEtaPtMap.push_back
      ( static_cast< TH2D* >
	( fIn->Get( Form("h_etaPtMap_%s", trigger.c_str() ))));
    m_vHtriggerEtaPtMap.back()->SetDirectory(0);
    
    // -------- spect --------
    m_vHtriggerYstarSpect.push_back
      ( static_cast< TH2D* >
	( fIn->Get
	  ( Form("h_%s_%s", m_ystarSpectName.c_str(), trigger.c_str() ))));
    m_vHtriggerYstarSpect.back()->SetDirectory(0);

    m_vHtriggerYstarSpectFine.push_back
      ( static_cast< TH2D* >
	( fIn->Get
	  ( Form("h_%s_%s", m_ystarSpectFineName.c_str(), trigger.c_str() ))));
    m_vHtriggerYstarSpectFine.back()->SetDirectory(0);

    m_vHtriggerYstarSpectFineRuns.push_back
      ( static_cast< TH3D* >
	( fIn->Get
	  ( Form("h_%s_%s", m_ystarSpectFineRunsName.c_str(), trigger.c_str() ))));
    m_vHtriggerYstarSpectFineRuns.back()->SetDirectory(0);

    // ----- efficiencies ----
    m_vHtriggerEtaSpectSim.push_back
      ( static_cast< TH2D* >
	( fIn->Get
	  ( Form("h_%sSim_%s", m_etaSpectName.c_str(), trigger.c_str() ))));
    m_vHtriggerEtaSpectSim.back()->SetDirectory(0);

    m_vHtriggerEtaSpectDenom.push_back
      ( static_cast< TH2D* >
	( fIn->Get
	  ( Form("h_%sDenom_%s", m_etaSpectName.c_str(), trigger.c_str() ))));
    m_vHtriggerEtaSpectDenom.back()->SetDirectory(0);
    // -------- dPhi- --------
    m_vHtriggerDphi.push_back
      ( static_cast< THnSparse* >
	( fIn->Get( Form("h_%s_%s", m_dPhiName.c_str(), trigger.c_str() ))));
  }
  
  fIn->Close(); delete fIn;
}

void DiJetAnalysisData::MakeRunSpectra( std::vector< TH3* >& vSampleSpect,
					const std::vector< std::string>& vLabels,
					const std::string& name ){

  /*
  if( !vSampleSpect.size() ){ return; }

  std::string axisLabel, axisLabelTex;
  GetSpectraLabels( axisLabel, axisLabelTex, name );

  bool isEta = name.find( m_sEta ) != std::string::npos ? true : false;

  // use this as reference because
  // it should be in every file
  TH3*  hRef = vSampleSpect[0];
  int nXbins = hRef->GetNbinsX();

  uint nSamples = vSampleSpect.size();
  
  std::vector< std::vector< std::vector< TH1* > > > vSpect;
  std::vector< std::vector< std::vector< TH1* > > > vSpectCounts;
  vSpect      .resize( nSamples );
  vSpectCounts.resize( nSamples );

  for( auto& v : vSpect       ){ v.resize( m_nRuns ); }
  for( auto& v : vSpectCounts ){ v.resize( m_nRuns ); }

  std::string yAxisTitle = "dN/d#it{p}_{T} [GeV]";

  double max = -1;
  
  for( uint iG = 0; iG < nSamples; iG++){

    std::string label = vLabels[iG];

    if( label.compare( m_allName ) ){ continue; }

    for( uint iRun = 1; iRun <= m_nRuns; iRun++ ){
    
      for( int xBin = 1; xBin <= nXbins; xBin++ ){
      
	double xMin, xMax;
	anaTool->GetBinRange
	  ( hRef->GetXaxis(), xBin, xBin, xMin, xMax );
      
	std::string hTag = anaTool->GetName( xMin, xMax, axisLabel);
       
	TH1* hSpectCounts =
	  vSampleSpect[iG]->
	  ProjectionY( Form("h_%s_%s_%s_%s",
			    name.c_str(), m_sCounts.c_str(), label.c_str(), hTag.c_str() ),
		       xBin, xBin );
      
	TH1* hSpect = static_cast< TH1D* >
	  ( hSpectCounts->Clone
	    ( Form("h_%s_%s_%s",
		   name.c_str(), label.c_str(), hTag.c_str() ) ) );
      
	hSpect->SetTitle( anaTool->GetLabel( xMin, xMax, axisLabelTex ).c_str() );
	hSpect->SetYTitle( yAxisTitle.c_str() );
      
	vSpect      [iG].push_back( hSpect );
	vSpectCounts[iG].push_back( hSpectCounts );
      
	// scale by width of bins to get dN/dpT
	hSpect->Scale( 1./m_mRunLumi(), "width" );
      
	hSpect->Write();
	hSpectCounts->Write();
      
	// get min max from the final histograms
	if( label.compare( m_allName ) ){ continue; }
	if( max < hSpect->GetMaximum() ){ max = hSpect->GetMaximum(); }
      
      } // end loop over xBin
    } // end loop over runs
  } // end loop over iG

  // set maxima globally for all spectra hists.
  // easier to compare. Set on log scale.

  max = anaTool->GetLogMaximum( max );
  
  double lX0, lY0, lX1, lY1;
  
  //------------------------------------------------
  //------- For Each xAxis Bin, draw an IGs --------
  //------------------------------------------------

  if( m_is_pPb ){ lX0 = 0.45; lY0 = 0.54; lX1 = 0.76; lY1 = 0.67; }
  else          { lX0 = 0.25; lY0 = 0.20; lX1 = 0.45; lY1 = 0.40; }
  
  for( int iX = 0; iX < nXbins; iX++ ){
    int    xBin = iX + 1;
    double xMin, xMax;
    anaTool->GetBinRange( hRef->GetXaxis(), xBin, xBin, xMin, xMax );
    double xCenter = hRef->GetXaxis()->GetBinCenter ( xBin );

    // for pPb, dont draw at anything above -3.2
    if( isEta && m_is_pPb && xCenter > -constants::FETAMIN ){ continue; }
     
    std::string cName  = anaTool->GetName ( xMin, xMax, axisLabel    );
    std::string cLabel = anaTool->GetLabel( xMin, xMax, axisLabelTex );
    
    TCanvas c( cLabel.c_str(), cLabel.c_str(), 800, 600 );
    c.SetLogy();
    
    TLegend leg( lX0, lY0, lX1, lY1 );
    styleTool->SetLegendStyle( &leg, 0.7 );
    
    int style = 1;
    for( uint iG = 0; iG < nSamples; iG++ ){
      std::string label = vLabels[iG];

      // for pp, dont draw central triggers below -3.2
      // or forward triggers above -3.2
      // for mb, and total draw everything
      bool isMb  = label.find("_mb_") != std::string::npos
	? true : false;
      bool isAll = !label.compare( m_allName )
	? true : false;

      TH1* h = vSpect[iG][iX];
           
      if( !h->GetEntries() && !isMb && !isAll ){ continue; }
      else if( !h->GetEntries() && !isMb && !isAll ){ continue; }
      
      if( iG == nSamples ){ style = 0; }
      styleTool->SetHStyle( h, style++ );
      h->Draw("epsame");
      h->SetMinimum( m_ptSpectYaxisMin );
      h->SetMaximum( max );
      leg.AddEntry( h, label.c_str() );
    } // end loop over iG
    
    leg.Draw("same");

    DrawAtlasRight();    
    drawTool->DrawRightLatex( 0.87, 0.73, cLabel );

    SaveAsAll( c, Form("%s_%s", name.c_str(), cName.c_str() ) );
  } // end loop over iX
  
  //------------------------------------------------
  //--------- For each IG, draw xAxisBins ----------
  //------------------------------------------------

  if( m_is_pPb ){ lX0 = 0.70; lY0 = 0.54; lX1 = 0.85; lY1 = 0.71; }
  else          { lX0 = 0.20; lY0 = 0.23; lX1 = 0.47; lY1 = 0.40; }
  
  for( uint iG = 0; iG < nSamples; iG++ ){
    std::string label = vLabels[iG];

    if( label.compare( m_allName ) ){ continue; };
    
    std::string cName  = label;
    std::string cLabel = label;

    TCanvas c( "c", cLabel.c_str(), 800, 600 );
    c.SetLogy();
    
    TLegend leg( lX0, lY0, lX1, lY1 );
    styleTool->SetLegendStyle( &leg, 0.7 );

    int style = 0;
    for( int iX = 0; iX < nXbins; iX++ ){
      int       xBin = iX + 1;
      double xCenter = hRef->GetXaxis()->GetBinCenter ( xBin );

      // for pPb, dont draw at anything above -3.2
      if( m_is_pPb && xCenter > -constants::FETAMIN ){ continue; }

      // for pp, dont draw central triggers below -3.2
      // or forward triggers above -3.2
      // for mb, and total draw everything
      bool isMb  = label.find("_mb_") != std::string::npos
	? true : false;
      bool isAll = !label.compare( m_allName )
	? true : false;

      TH1* h = vSpect[iG][iX];
      
      if( h->GetEntries() && !isMb && !isAll ){ continue; }
      else if( h->GetEntries() && !isMb && !isAll ){ continue; }

      styleTool->SetHStyle( h, style++ );
      h->Draw("epsame");
      h->SetMinimum( m_ptSpectYaxisMin );
      h->SetMaximum( max );
      leg.AddEntry( h, h->GetTitle() );
    } // end loop over iX
    leg.Draw("same");

    DrawAtlasRight();    

    SaveAsAll( c, Form("%s_%s", name.c_str(), cName.c_str() ) );
  } // end loop over iG  

  // delete
  for( uint iG = 0; iG < nSamples; iG++ ){
    for( int iX = 0; iX < nXbins; iX++ ){
      delete vSpect      [iG][iX];
      delete vSpectCounts[iG][iX];
    }
  }
  */
}
					  
void DiJetAnalysisData::MakeEfficiencies( std::vector< TH2* >& vTrigSpect,
					  std::vector< TH2* >& vTrigSpectRef,
					  const std::string& type ){
  double lX0, lY0, lX1, lY1;

  if( m_is_pPb ){ lX0 = 0.13; lY0 = 0.75; lX1 = 0.39; lY1 = 0.87; }
  else          { lX0 = 0.13; lY0 = 0.68; lX1 = 0.39; lY1 = 0.89; }
  
  // us m_hAllEtaSpect because its always there
  double xMin = 0;
  double xMax = 100;

  // to use for reference on axis etc.
  TH2* hExample = NULL;
  if( vTrigSpect.size() ){
    hExample = vTrigSpect.front();
  } else{
    return;
  }
  
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
  //-------- Draw Eta Bins for each Trigger --------
  //------------------------------------------------
  for( uint iG = 0; iG < m_nTriggers; iG++ ){
    std::string trigger = m_vTriggers[iG];
    
    // dont draw MB trigger
    if( !trigger.compare( m_mbTriggerName ) ){ continue; }
    if( !trigger.compare( m_allName       ) ){ continue; }

    std::string cName  = trigger;
    std::string cLabel = trigger;

    TCanvas c( cLabel.c_str(), cLabel.c_str(), 800, 600 );
    styleTool->SetCStyleGraph( c, xMin, m_effMin, xMax, m_effMax, gLabel );
    
    TLegend leg( lX0, lY0, lX1, lY1 );
    styleTool->SetLegendStyle( &leg  );

    int style = 0;
    for( int iX = 0; iX < nXbins; iX++ ){
      int         xBin = iX + 1;
      double etaCenter = hExample->GetXaxis()->GetBinCenter ( xBin );

      // temporary, dont draw the 3.1->3.2 bin
      if( std::abs(etaCenter) < 3.2 && std::abs(etaCenter) > 3.1 ){ continue; }
      // for pPb, dont draw at anything above -3.2
      if( m_is_pPb && etaCenter > -constants::FETAMIN ){ continue; }

      // for pp, dont draw central triggers below -3.2
      // or forward triggers above -3.2
      if( !m_is_pPb && etaCenter < -constants::FETAMIN &&
	  trigger.find("320eta490") == std::string::npos ){ continue; }
      else if( !m_is_pPb && etaCenter > -constants::FETAMIN &&
	       trigger.find("320eta490") != std::string::npos){ continue; }
      
      TGraphAsymmErrors* g = vEffGrf[iG][iX];
      styleTool->SetHStyle( g, style++ );
      g->Draw("p");
      leg.AddEntry( g, g->GetTitle() );
      g->SetTitle("");
    } // end loop over iX
    
    leg.Draw("same");

    TLine line( xMin, 1, xMax, 1);
    line.Draw();
    
    DrawAtlasRight();
    drawTool->DrawRightLatex( 0.88, 0.76, cLabel );
    
    SaveAsAll( c, Form("%s_%s", type.c_str(), cName.c_str() ) );
  } // end loop over iG

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
    styleTool->SetLegendStyle( &leg  );
    
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
      leg.AddEntry( g, trigger.c_str() );
    } // end loop over iG
    
    leg.Draw("same");

    TLine line( xMin, 1, xMax, 1);
    line.Draw();    
    
    DrawAtlasRight();
    drawTool->DrawRightLatex( 0.88, 0.76, cLabel );

    SaveAsAll( c, Form("%s_%s", type.c_str(), cName.c_str() ) );
  } // end loop over iX

  // delete
  for( uint iG = 0; iG < m_nTriggers; iG++ ){
    for( int iX = 0; iX < nXbins; iX++ )
      {
	delete vSpect[iG][iX];
	delete vSpectRef[iG][iX];
	delete vEffGrf[iG][iX];
      }
  } 
}

void DiJetAnalysisData::MakeSystematicsGraphs( TFile* fOut, const std::string& name ){

  std::vector< int > v_uc;
  std::map< int, TFile* > mFinUC;

  // see if we are dealing with widths of yields.
  bool isYield = name.find( m_yieldName ) != std::string::npos ? true : false;
  
  // get map of sys factor to TFile*
  TFile* fInNominal = GetListOfSystUncert( v_uc, mFinUC );
  // change back to the fOut, for writing.
  fOut->cd();
  
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

  double pDx = 0.05;
  
  std::string yTitle = isYield ?
    "Pair Jet Per Jet_{1} Yield" : "|#Delta#phi| width";
  std::string xTitle = m_dPP->GetAxisLabel(3);
  std::string gTitle = ";" + xTitle + ";" + yTitle;

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

      TCanvas c( "c", "c", 800, 600 );
      styleTool->SetCStyleGraph
	( c, x0, y0, x1, y1, gTitle.c_str() );

      // for yields, use log scale on yaxis
      if( isYield ){ c.SetLogy(); }

      double legX0, legX1, legY0, legY1;

      // this is bs. works though. vary the last factor 
      double deltaYleg = ( axis1Bin - 1 ) * 0.075;
 
      if( isYield ){
	legX0 = 0.30;
	legY0 = 0.22;
	legX1 = 0.71;
	legY1 = 0.29 + deltaYleg;
      } else {
	legX0 = 0.54;
	legY0 = 0.22;
	legX1 = 0.89;
	legY1 = 0.29 + deltaYleg;
      }
      
      TLegend leg( legX0, legY0, legX1, legY1 );
      styleTool->SetLegendStyle( &leg );

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
	
	TH1D* hNominal = static_cast<TH1D*>( fInNominal->Get( hNominalName.c_str() ) );
	vHdef.push_back( hNominal );

	// Make all the systematics histograms for this y1, pt1, pt2 bin.
	// since we look one bin at a time and process all systematics,
	// need to have these histograms ready
	std::map< int, TH1D* > mHsystTmp; 

	for( auto uc : v_uc ){

	  // skip uc = 0 (default)
	  if( !uc ){ continue; }

	  std::string    uncertSuffix = uc > 0 ? Form("P%d", uc) : Form("N%d", -1 * uc) ;
	  std::string hSystematicName = "h_" + allUnfoldedName + "_" + hTag + "_" + m_uncertSuffix;

	  TH1D* hSystematic = static_cast<TH1D*>( hNominal->Clone( hSystematicName.c_str() ) );
	  vHsyst.push_back( hSystematic );
	  hSystematic->Reset();
	  mHsystTmp[ uc ] = hSystematic;
	}
	
	std::vector< double > pX;
	std::vector< double > eX( nAxis3Bins, 0.1 * hNominal->GetBinWidth(1) );
 
	std::vector< double > pY;
	std::vector< double > eYP;
	std::vector< double > eYN;

	std::vector< double > eYPJES;
	std::vector< double > eYNJES;

	std::vector< double > eYPJER;
	std::vector< double > eYNJER;
	
	//---------------------------------------------------
	//------------------ DO WORK HERE -------------------
	//---------------------------------------------------
	for( int axis3Bin = 1; axis3Bin <= nAxis3Bins; axis3Bin++ ){

	  std::vector< TH1* > vHunc;
	  
	  std::vector< double > eYPtmp;
	  std::vector< double > eYNtmp;

	  std::vector< double > eYPtmpJES;
	  std::vector< double > eYNtmpJES;

	  std::vector< double > eYPtmpJER;
	  std::vector< double > eYNtmpJER;
	  
	  double yNominal = hNominal->GetBinContent( axis3Bin );
	  pX.push_back( hNominal->GetBinCenter( axis3Bin ) );
	  
	  // loop over uncertainties
	  for( auto uc : v_uc ){

	    // skip uc = 0 (default)
	    if( !uc ){ continue; }

	    // this is same as hNominalName but leave it separate
	    // to avoid confusion / make more flexible later.
	    std::string hUncertaintyName = "h_" + allUnfoldedName + "_" + hTag;

	    TH1D* hUncertainty = static_cast<TH1D*>
	      ( mFinUC[ uc ]->Get( hUncertaintyName.c_str() ) );
	    vHunc.push_back( hUncertainty );

	    int sign = uc > 0 ? 1 : -1;

	    double yShifted = hUncertainty->GetBinContent( axis3Bin );

	    double uncertainty  = ( yNominal - yShifted )/yNominal;

	    std::cout << "++++" << uc << " ++++" << axis3Bin << " " << axis1Bin << " " << axis2Bin << " " 
		      << sign << " " << yShifted << " " << yNominal << " " << uncertainty << std::endl;

	    mHsystTmp[ uc ]->SetBinContent( axis3Bin, uncertainty );
	    
	    // for JER negative is same as positive (20)
	    // for Fitting negative is same as positive (22)
	    // for ReWeight negative is same as positive (23)
	    // and we do not have a NN, just use PN
	    if( uc  == 20 || uc == 22 || uc == 23 ){
	      eYPtmp.push_back( uncertainty );
	      eYNtmp.push_back( uncertainty );
	      continue;
	    }
	    if( sign > 0 ){
	      eYPtmp.push_back( uncertainty );
	    } else {
	      eYNtmp.push_back( uncertainty );
	    }
	  } // end loop over uncertainties

	  // add uncertainties in quadrature;
	  double uncertaintyFinalYP = 0;
	  double uncertaintyFinalYN = 0;

	  // clean up, or there are memory problems
	  // because each file writes same histo into
	  // same memory address. so eventually you start
	  // deleting deleted stuff if you dont do it now
	  for( auto& h : vHunc ){ delete h; }

	  for( auto u : eYPtmp ){
	    uncertaintyFinalYP += std::pow( u , 2 );
	  }

	  for( auto u : eYNtmp ){
	    uncertaintyFinalYN += std::pow( u , 2 );
	  }

	  uncertaintyFinalYP =
	    uncertaintyFinalYP >= 0 ? std::sqrt( uncertaintyFinalYP ) : 0.0;

	  uncertaintyFinalYN =
	    uncertaintyFinalYN >= 0 ? std::sqrt( uncertaintyFinalYN ) : 0.0;

	  std::cout << hNominalName << " "
		    << uncertaintyFinalYP << " "
		    << uncertaintyFinalYN << std::endl;
	  
	  pY .push_back( yNominal );
	  eYP.push_back( yNominal * uncertaintyFinalYP );
	  eYN.push_back( yNominal * uncertaintyFinalYN );
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
	style++;

	// add some displacement along x
	double* xDef      = gNominal->GetX();
	double* eXDefLow  = gNominal->GetEXlow();
	double* eXDefHigh = gNominal->GetEXhigh();
	double* xSys      = gSystematics->GetX();
	for( int iX = 0; iX < nAxis3Bins; iX++ ){
	  *(      xDef + iX ) += style * pDx;
	  *(  eXDefLow + iX ) = 0;
	  *( eXDefHigh + iX ) = 0;
	  *(      xSys + iX ) += style * pDx;
	}
	
	gSystematics->SetTitle("");
       	gSystematics->GetXaxis()->SetRangeUser( x0, x1 );
       	gSystematics->GetYaxis()->SetRangeUser( y0, y1 );
	
	gNominal->SetTitle("");
       	gNominal->GetXaxis()->SetRangeUser( x0, x1 );
       	gNominal->GetYaxis()->SetRangeUser( y0, y1 );

	leg.AddEntry
	  ( gNominal, anaTool->GetLabel( axis2Low, axis2Up, m_dPP->GetAxisLabel(2) ).c_str(), "lp" );	
	
	// draw systematics first
	gSystematics->Draw("2");
	gNominal->Draw("p");
      } // end loop over axis2

      // Draw the final canvas with all of the graphs.
      leg.Draw();
      
      DrawTopLeftLabels
	( m_dPP, axis0Low, axis0Up, axis1Low, axis1Up, 0, 0, 0, 0);

      DrawAtlasRight();

      std::string hNameFinal =
	"h_" + name + "_" + m_sFinal + "_" + hTagC;

      SaveAsAll( c, hNameFinal );
    } // end loop over axis1     
  } // end loop over axis0

  for( auto& h : vHdef  ){ h->Write(); delete h; }
  for( auto& h : vHsyst ){ h->Write(); delete h; }
  for( auto& g : vG     ){ g->Write(); delete g; }
}

void DiJetAnalysisData::MakeFinalPlotsTogether( TFile* fOut, const std::string& name ){

  /*
  std::vector< int > v_uc;
  std::map< int, TFile* > mFinUC;

  // see if we are dealing with widths of yields.
  bool isYield = name.find( m_yieldName ) != std::string::npos ? true : false;
  
  // get map of sys factor to TFile*
  TFile* fInNominal = GetListOfSystUncert( v_uc, mFinUC );
  // change back to the fOut, for writing.
  fOut->cd();
  
  std::string allUnfoldedName    =
    m_dPhiName + "_" + m_unfoldedName    + "_" + name + "_" + m_allName;
  std::string allSystematicsName =
    m_dPhiName + "_" + m_systematicsName + "_" + name + "_" + m_allName;

  std::vector< TH1* > vHdef;
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

  double pDx = 0.05;
  
  std::string yTitle = isYield ?
    "Pair Jet Per Jet_{1} Yield" : "|#Delta#phi| width";
  std::string xTitle = m_dPP->GetAxisLabel(3);
  std::string gTitle = ";" + xTitle + ";" + yTitle;

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

      TCanvas c( "c", "c", 800, 600 );
      styleTool->SetCStyleGraph
	( c, x0, y0, x1, y1, gTitle.c_str() );

      // for yields, use log scale on yaxis
      if( isYield ){ c.SetLogy(); }

      double legX0, legX1, legY0, legY1;

      // this is bs. works though. vary the last factor 
      double deltaYleg = ( axis1Bin - 1 ) * 0.075;
 
      if( isYield ){
	legX0 = 0.30;
	legY0 = 0.22;
	legX1 = 0.71;
	legY1 = 0.29 + deltaYleg;
      } else {
	legX0 = 0.54;
	legY0 = 0.22;
	legX1 = 0.89;
	legY1 = 0.29 + deltaYleg;
      }
      
      TLegend leg( legX0, legY0, legX1, legY1 );
      styleTool->SetLegendStyle( &leg );

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
	
	TH1D* hNominal = static_cast<TH1D*>( fInNominal->Get( hNominalName.c_str() ) );
	vHdef.push_back( hNominal );

	std::vector< double > pX;
	std::vector< double > eX( nAxis3Bins, 0.1 * hNominal->GetBinWidth(1) );
 
	std::vector< double > pY;
	std::vector< double > eYP;
	std::vector< double > eYN;

	std::vector< double > eYPJES;
	std::vector< double > eYNJES;

	std::vector< double > eYPJER;
	std::vector< double > eYNJER;
		
	//---------------------------------------------------
	//------------------ DO WORK HERE -------------------
	//---------------------------------------------------
	for( int axis3Bin = 1; axis3Bin <= nAxis3Bins; axis3Bin++ ){

	  std::vector< TH1* > vHunc;
	  
	  std::vector< double > eYPtmp;
	  std::vector< double > eYNtmp;

	  std::vector< double > eYPtmpJES;
	  std::vector< double > eYNtmpJES;

	  std::vector< double > eYPtmpJER;
	  std::vector< double > eYNtmpJER;
	  
	  double yNominal = hNominal->GetBinContent( axis3Bin );
	  pX.push_back( hNominal->GetBinCenter( axis3Bin ) );
	  
	  // loop over uncertainties
	  for( auto uc : v_uc ){

	    // skip uc = 0 (default)
	    if( !uc ){ continue; }

	    // this is same as hNominalName but leave it separate
	    // to avoid confusion / make more flexible later.
	    std::string hUncertaintyName = "h_" + allUnfoldedName + "_" + hTag;

	    TH1D* hUncertainty = static_cast<TH1D*>
	      ( mFinUC[uc]->Get( hUncertaintyName.c_str() ) );
	    vHunc.push_back( hUncertainty );

	    int sign = uc > 0 ? 1 : -1;

	    double yShifted = hUncertainty->GetBinContent( axis3Bin );

	    double uncertainty  = ( yNominal - yShifted )/yNominal;

	    std::cout << "++++" << uc << " ++++" << axis3Bin << " " << axis1Bin << " " << axis2Bin << " " 
		      << sign << " " << yShifted << " " << yNominal << " " << uncertainty << std::endl;

	    // for JER negative is same as positive
	    // and we do not have a N20, just use P20
	    if( uc  == 20 ){
	      eYPtmp.push_back( uncertainty );
	      eYNtmp.push_back( uncertainty );
	      continue;
	    }
	    // for Fitting negative is same as positive
	    // and we do not have a N20, just use P20
	    if( uc  == 22 ){
	      eYPtmp.push_back( uncertainty );
	      eYNtmp.push_back( uncertainty );
	      continue;
	    }
	    // for ReWeight negative is same as positive
	    // and we do not have a N20, just use P20
	    if( uc  == 23 ){
	      eYPtmp.push_back( uncertainty );
	      eYNtmp.push_back( uncertainty );
	      continue;
	    }
	    
	    if( sign > 0 ){
	      eYPtmp.push_back( uncertainty );
	    } else {
	      eYNtmp.push_back( uncertainty );
	    }
	  } // end loop over uncertainties

	  // add uncertainties in quadrature;
	  double uncertaintyFinalYP = 0;
	  double uncertaintyFinalYN = 0;

	  // clean up, or there are memory problems
	  // because each file writes same histo into
	  // same memory address. so eventually you start
	  // deleting deleted stuff if you dont do it now
	  for( auto& h : vHunc ){ delete h; }

	  for( auto u : eYPtmp ){
	    uncertaintyFinalYP += std::pow( u , 2 );
	  }

	  for( auto u : eYNtmp ){
	    uncertaintyFinalYN += std::pow( u , 2 );
	  }

	  uncertaintyFinalYP =
	    uncertaintyFinalYP >= 0 ? std::sqrt( uncertaintyFinalYP ) : 0.0;

	  uncertaintyFinalYN =
	    uncertaintyFinalYN >= 0 ? std::sqrt( uncertaintyFinalYN ) : 0.0;

	  std::cout << hNominalName << " "
		    << uncertaintyFinalYP << " "
		    << uncertaintyFinalYN << std::endl;
	  
	  pY .push_back( yNominal );
	  eYP.push_back( yNominal * uncertaintyFinalYP );
	  eYN.push_back( yNominal * uncertaintyFinalYN );
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
	style++;

	// add some displacement along x
	double* xDef      = gNominal->GetX();
	double* eXDefLow  = gNominal->GetEXlow();
	double* eXDefHigh = gNominal->GetEXhigh();
	double* xSys      = gSystematics->GetX();
	for( int iX = 0; iX < nAxis3Bins; iX++ ){
	  *(      xDef + iX ) += style * pDx;
	  *(  eXDefLow + iX ) = 0;
	  *( eXDefHigh + iX ) = 0;
	  *(      xSys + iX ) += style * pDx;
	}
	
	gSystematics->SetTitle("");
       	gSystematics->GetXaxis()->SetRangeUser( x0, x1 );
       	gSystematics->GetYaxis()->SetRangeUser( y0, y1 );
	
	gNominal->SetTitle("");
       	gNominal->GetXaxis()->SetRangeUser( x0, x1 );
       	gNominal->GetYaxis()->SetRangeUser( y0, y1 );

	leg.AddEntry
	  ( gNominal, anaTool->GetLabel( axis2Low, axis2Up, m_dPP->GetAxisLabel(2) ).c_str(), "lp" );	
	
	// draw systematics first
	gSystematics->Draw("2");
	gNominal->Draw("p");
	
	gSystematics->Write();
	gNominal->Write();
      } // end loop over axis2

      // Draw the final canvas with all of the graphs.
      leg.Draw();
      
      DrawTopLeftLabels
	( m_dPP, axis0Low, axis0Up, axis1Low, axis1Up, 0, 0, 0, 0);

      DrawAtlasRight();

      std::string hNameFinal =
	"h_" + name + "_" + m_sFinal + "_" + hTagC;

      SaveAsAll( c, hNameFinal );
      
    } // end loop over axis1     
  } // end loop over axis0

  for( auto& h : vHdef ){ delete h; }
  for( auto& g : vG    ){ delete g; }
  */
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

void DiJetAnalysisData::DrawAtlasRightBoth( double x0, double y0, double scale ){
  drawTool->DrawAtlasInternal( scale );
  drawTool->DrawRightLatex( 0.87, 0.87, "Data" );	  	  
} 
