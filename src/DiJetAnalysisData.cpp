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
  : DiJetAnalysis( is_pPb, true, mcType, uncertComp )
{
  //========== Set Histogram Binning =============

  //==================== Cuts ====================

  //================== Settings ==================
  m_mbTriggerName = "";
  
  m_mbTriggerI         = -1;
  m_lowestCentTriggerI = -1;
  
  m_centMbCorrection   = 0;

  //=============== Histo Names ==================    
}

DiJetAnalysisData::~DiJetAnalysisData(){}

void DiJetAnalysisData::Initialize(){
  // Initalize things common to everything
  DiJetAnalysis::Initialize();
  
  m_fNameIn = m_is_pPb ?
    "/home/yakov/Projects/atlas/data/pPb/pPb.root" :
    "/home/yakov/Projects/atlas/data/pp/pp.root"  ;
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

  // MakeEtaPhiPtMap( m_vHtriggerEtaPhiMap );
  // MakeEtaPhiPtMap( m_vHtriggerEtaPtMap  );

  // add a trigger "all" to collection
  // rest of plots include combined triggers
  m_vTriggers.push_back( m_allName );
  
  m_hAllEtaSpect = CombineSamples( m_vHtriggerEtaSpect, m_etaSpectName );
  MakeSpectra( m_vHtriggerEtaSpect, m_vTriggers, m_etaSpectName );

  m_hAllYstarSpect = CombineSamples( m_vHtriggerYstarSpect, m_ystarSpectName );
  MakeSpectra( m_vHtriggerYstarSpect, m_vTriggers, m_ystarSpectName );

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
  TFile* fInData = TFile::Open( m_fNamePerf.c_str()            );
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
  TFile* fInMCPerf  = TFile::Open( m_fNamePerfUF.c_str() );

  TFile* fOut = new TFile( m_fNamePhys.c_str(),"RECREATE");

  // add a trigger "all" to collection
  // rest of plots include combined triggers
  m_vTriggers.push_back( m_allName );

  m_hAllDphi = CombineSamples( m_vHtriggerDphi, m_dPhiName );
  MakeDeltaPhi( m_vHtriggerDphi, m_vTriggers, m_dPhiName,
		fInMCPerf, m_ystarSpectUnfoldedName );
  
  std::cout << "DONE! Closing " << fOut->GetName() << std::endl;
  fOut->Close(); delete fOut;
  std::cout << "......Closed  " << std::endl;
}

void DiJetAnalysisData::UnfoldPhysics(){

  // for now
  // Copy File with original dPhi, spectra, etc,
  // into file where histos with corrections are
  // going to be appended. 
  std::cout << "Copy " << m_fNamePhys << " -> " << m_fNamePhysUF << std::endl;
  TFile::Cp( m_fNamePhys.c_str(), m_fNamePhysUF.c_str() );

  // Open two for reading one for updating.
  // open the MC file used for unfolding info.
  // open teh data file used for measured info.
  // passed to unfolding function.
  TFile* fInMC     = TFile::Open( m_fNamePhysUnfoldingMC.c_str()  );
  TFile* fInData   = TFile::Open( m_fNameDefPhys.c_str()          );
  TFile* fInMCPerf = TFile::Open( m_fNamePerfUF.c_str() );
  TFile* fOut      = new TFile( m_fNamePhysUF.c_str(),"UPDATE"    );

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
}

void DiJetAnalysisData::ProcessSystematics(){
  MakeSystematicsGraphs();
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
		   ";#eta;#it{p}_{T}",
		   m_nEtaMapBins, m_etaMapMin, m_etaMapMax,
		   m_nPtMapBins , m_ptMapMin , m_ptMapMax) ) ;
    AddHistogram( m_vHtriggerEtaPtMap.back() );

    // -------- spect --------
    m_vHtriggerEtaSpect.push_back
      ( new TH2D( Form("h_%s_%s", m_etaSpectName.c_str(), trigger.c_str() ), 
		  ";#eta;#it{p}_{T} [GeV];",
		  m_nVarFwdEtaBins, 0, 1,
		  m_nPtSpectBins,m_ptSpectMin, m_ptSpectMax ) ) ;
    m_vHtriggerEtaSpect.back()->GetXaxis()->
      Set( m_nVarFwdEtaBins, &( m_varFwdEtaBinning[0] ) );
    AddHistogram( m_vHtriggerEtaSpect.back() );

    m_vHtriggerYstarSpect.push_back
      ( new TH2D( Form("h_%s_%s", m_ystarSpectName.c_str(), trigger.c_str() ), 
		  ";#it{y}*;#it{p}_{T} [GeV];",
		  m_nVarYstarBins, 0, 1,
		  m_nVarPtBinsRespMat, 0, 1 ) ) ;
    m_vHtriggerYstarSpect.back()->GetXaxis()->
      Set( m_nVarYstarBins, &( m_varYstarBinning[0] ) );
    m_vHtriggerYstarSpect.back()->GetYaxis()->
      Set( m_nVarPtBinsRespMat, &( m_varPtBinningRespMat[0] ) );
   
    AddHistogram( m_vHtriggerYstarSpect.back() );

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
	if( JetInTrigRange( *jet1, iG ) ) {
	  AnalyzeDeltaPhi( m_vHtriggerDphi[iG], vR_jets );
	}
      }
  
      // loop over jets 
      for( auto& jet : vR_jets ){
	// cleaned jets have px, py, pz set to 0
	if( jet.Pt() <= 0 ) continue;
	
	// ETA-PHI
	double jetEta = jet.Eta();
	double jetPhi = jet.Phi();
	double jetPt  = jet.Pt()/1000.;
	
	m_vHtriggerEtaPhiMap[iG]->Fill( jetEta, jetPhi );
	m_vHtriggerEtaPtMap [iG]->Fill( jetEta, jetPt ); 

	// check if the jet is in appropriate range
	// for the trigger fired
	//	if( !JetInTrigRange( jet, iG ) ){ continue; };

	double jetYstar = GetYstar( jet );
	
	// fill spectra 
	m_vHtriggerEtaSpect[iG]->
	  Fill( jetEta  , jetPt );
	m_vHtriggerYstarSpect[iG]->
	  Fill( jetYstar, jetPt );

	// for pp fill both, pPb continue
	if( m_is_pPb ){ continue; }

	m_vHtriggerEtaSpect[iG]->
	  Fill( -jetEta  , jetPt );
	m_vHtriggerYstarSpect[iG]->
	  Fill( -jetYstar, jetPt );

      } // end loop over jets
    } // end loop over iG
    // EFFICIENCIES - only for good run and LBN
    if( goodRunLBN ){ AnalyzeEff( vR_jets, vTrig_jets, mTriggerFired ); }    
  } // -------- END EVENT LOOP ---------
  std::cout << "DONE! Has: " << nEventsTotal << " events." << std::endl;

  fIn->Close(); delete fIn;
}

//---------------------------
//       Analysis
//---------------------------

bool DiJetAnalysisData::JetInTrigPtRange( const TLorentzVector& jet, int iG,
					  double ptHighExtra ){  
  double jetPt = jet.Pt() / 1000.;
  if( jetPt > m_vTriggersEffPtLow[iG] &&
      ( jetPt < m_vTriggersEffPtHigh[iG] + ptHighExtra ||
	m_vTriggersEffPtHigh[iG] < 0 ) )
    { return true; }
  return false;
}

bool DiJetAnalysisData::JetInTrigEtaRange( const TLorentzVector& jet, int iG ){
  if( iG == m_mbTriggerI ){ return true; }
  double jetEta = jet.Eta();
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

bool DiJetAnalysisData::JetInTrigRange( const TLorentzVector& jet, int iG ){
  double jetEta = jet.Eta();

  // in central triggers, we begin with j20
  // this corresponds to ~30 efficiency
  // in forward, we begin with j10
  // this corresponds to ~20 efficiency.
  // so we adjust here where the mb range goes to.
  bool applyCentMbCorrection =
    ( iG == m_mbTriggerI && IsCentralDetector( jetEta ) ) ? true : false;
  
  if( JetInTrigPtRange
      ( jet, iG, applyCentMbCorrection ? m_centMbCorrection : 0 ) &&
      JetInTrigEtaRange( jet, iG ) )
    { return true; }
  
  return false;
}

bool DiJetAnalysisData::TrigJetAboveThold( const TLorentzVector& jet, int iG ){
  if( jet.Pt()/1000. >= m_vTriggersTholdPt[iG] ){ return true; }
  return false;
}

bool DiJetAnalysisData::TrigJetInTrigRange( const TLorentzVector& jet, int iG ){
  if( TrigJetAboveThold( jet, iG ) && JetInTrigEtaRange( jet, iG ) )
    { return true; }
  return false;
}

void DiJetAnalysisData::AnalyzeEff( std::vector< TLorentzVector >& vR_jets,
				    std::vector< TLorentzVector >& vTrig_jets,
				    std::map< int, bool >& mTriggerFired ){
  if( vTrig_jets.empty() ){ return; }

  // sort trigger jets
  std::sort( vTrig_jets.begin(), vTrig_jets.end(),
	     anaTool->sortByDecendingPt );

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
    if( !( tJetFwd  ? TrigJetInTrigRange( *tJetFwd , iG ) : false ) &&
	!( tJetCent ? TrigJetInTrigRange( *tJetCent, iG ) : false ) )
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

THnSparse* DiJetAnalysisData::CombineSamples( std::vector< THnSparse* >& vSampleHin,
					      const std::string& name  ){
  if( !vSampleHin.size() ){ return NULL; }
  
  THnSparse* h_res = static_cast< THnSparseD* >
    ( vSampleHin[0]->Clone( Form("h_%s_%s", name.c_str() , m_allName.c_str() ) ) );
  h_res->Reset();
  
  for( uint iG = 0; iG < vSampleHin.size(); iG++ ){
    double scale = m_vTriggersPrescale[iG];
    h_res->Add( vSampleHin[iG], scale );
  }
  vSampleHin.push_back( h_res );
  
  return h_res;
}

void DiJetAnalysisData::GetInfoTogether( std::string& name_a  , std::string& name_b  ,
					 std::string& label_a , std::string& label_b ,
					 std::string& fName_a , std::string& fName_b ){

  int combinationBoth = GetConfig()->GetValue( "combinationBoth", 0 );

  if( combinationBoth == 0 ){
    // here a and b are for same file. one is for
    // default points, one is for systematics.
    name_a    = m_dPhiUnfoldedName    + "_" + m_allName;
    name_b    = m_dPhiSystematicsName + "_" + m_allName;
    label_a   = "#it{p}+Pb";
    label_b   = "#it{pp}";
    m_is_pPb  = true;  Initialize();
    fName_a   = m_fNameSYS;
    m_is_pPb  = false; Initialize();
    fName_b   = m_fNameSYS;
  } else if( combinationBoth == 1 ){
    name_a    = m_dPhiUnfoldedName + "_" + m_allName;
    name_b    = m_dPhiUnfoldedName + "_" + m_allName;
    label_a   = "#it{p}+Pb";
    label_b   = "#it{pp}";
    m_is_pPb  = true;  Initialize();
    fName_a   = m_fNamePhysUF;
    m_is_pPb  = false; Initialize();
    fName_b   = m_fNamePhysUF;
  } else if( combinationBoth == 2 ){
    name_a    = m_dPhiName + "_"  + m_allName;
    name_b    = m_dPhiName + "_"  + m_allName;
    label_a   = "#it{p}+Pb";
    label_b   = "#it{pp}";
    m_is_pPb  = true;  Initialize();
    fName_a   = m_fNamePhysUF;
    m_is_pPb  = false; Initialize();
    fName_b   = m_fNamePhysUF;
  } else if ( combinationBoth == 3 ){
    name_a    = m_dPhiUnfoldedName + "_" + m_allName;
    name_b    = m_dPhiName + "_" + m_allName;
    label_a   = "Unfolded #it{pp}";
    label_b   = "#it{pp}";
    m_is_pPb  = false; Initialize();
    fName_a   = m_fNamePhysUF;
    fName_b   = m_fNamePhysUF;
  } else if ( combinationBoth == 4 ){
    name_a    = m_dPhiUnfoldedName + "_" + m_allName;
    name_b    = m_dPhiName + "_" + m_allName;
    label_a   = "Unfolded #it{p}+Pb";
    label_b   = "#it{p}+Pb";
    m_is_pPb  = true; Initialize();
    fName_a   = m_fNamePhysUF;
    fName_b   = m_fNamePhysUF;
  }
}

void DiJetAnalysisData::GetSpectUnfoldingInfo( std::string& measuredName,
					       std::string& truthName,
					       std::string& respMatName,
					       std::string& unfoldedLabel,
					       std::string& typeLabel ){

  measuredName  = m_ystarSpectName;
  truthName     = m_ystarSpectTruthName;
  respMatName   = m_ystarSpectRespMatName;
  unfoldedLabel = "dN/d#it{p}_{T}";;
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

//---------------------------------
//     Get Quantities / Plot 
//---------------------------------

void DiJetAnalysisData::LoadHistograms(){
  
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
    m_vHtriggerEtaSpect.push_back
      ( static_cast< TH2D *>
	( fIn->Get
	  ( Form("h_%s_%s", m_etaSpectName.c_str(), trigger.c_str() ))));
    m_vHtriggerEtaSpect.back()->SetDirectory(0);

    m_vHtriggerYstarSpect.push_back
      ( static_cast< TH2D *>
	( fIn->Get
	  ( Form("h_%s_%s", m_ystarSpectName.c_str(), trigger.c_str() ))));
    m_vHtriggerYstarSpect.back()->SetDirectory(0);
 
    // ----- efficiencies ----
    m_vHtriggerEtaSpectSim.push_back
      ( static_cast< TH2D *>
	( fIn->Get
	  ( Form("h_%sSim_%s", m_etaSpectName.c_str(), trigger.c_str() ))));
    m_vHtriggerEtaSpectSim.back()->SetDirectory(0);

    m_vHtriggerEtaSpectDenom.push_back
      ( static_cast< TH2D *>
	( fIn->Get
	  ( Form("h_%sDenom_%s", m_etaSpectName.c_str(), trigger.c_str() ))));
    m_vHtriggerEtaSpectDenom.back()->SetDirectory(0);
    // -------- dPhi- --------
    m_vHtriggerDphi.push_back
      ( static_cast< THnSparse *>
	( fIn->Get( Form("h_%s_%s", m_dPhiName.c_str(), trigger.c_str() ))));
  }
  
  fIn->Close(); delete fIn;
}

void DiJetAnalysisData::MakeEfficiencies( std::vector< TH2* >& vTrigSpect,
					  std::vector< TH2* >& vTrigSpectRef,
					  const std::string& type ){
  double lX0, lY0, lX1, lY1;

  if( m_is_pPb ){ lX0 = 0.13; lY0 = 0.75; lX1 = 0.39; lY1 = 0.87; }
  else          { lX0 = 0.13; lY0 = 0.68; lX1 = 0.39; lY1 = 0.89; }
  
  // us m_hAllEtaSpect because its always there
  double xMin = m_hAllEtaSpect->GetYaxis()->GetXmin() - 10;
  double xMax = 100;

  // use this as reference because
  // it should be in every file
  int nXbins = m_hAllEtaSpect->GetNbinsX();
  
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
	( m_hAllEtaSpect->GetXaxis(), xBin, xBin, etaMin, etaMax );
            
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
      double etaCenter = m_hAllEtaSpect->GetXaxis()->GetBinCenter ( xBin );

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
      ( m_hAllEtaSpect->GetXaxis(), xBin, xBin, etaMin, etaMax );
    double etaCenter = m_hAllEtaSpect->GetXaxis()->GetBinCenter ( xBin );
    
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

void DiJetAnalysisData::MakeSystematicsGraphs(){

  std::vector< int > v_uc;
  std::map< int, TFile* > mFinUC;

  TFile* fInDefault = GetListOfSystUncert( v_uc, mFinUC );
  
  TFile* fOut  = new TFile( m_fNameSYS.c_str(), "recreate");
  
  std::string allUnfoldedName    = m_dPhiUnfoldedName    + "_" + m_allName;
  std::string allSystematicsName = m_dPhiSystematicsName + "_" + m_allName;

  std::vector< TH1* > vHdef;
  std::vector< TGraphAsymmErrors* > vGw;
  
  TAxis* axis0 = m_dPP->GetTAxis(0); int nAxis0Bins = axis0->GetNbins();
  TAxis* axis1 = m_dPP->GetTAxis(1); int nAxis1Bins = axis1->GetNbins();
  TAxis* axis2 = m_dPP->GetTAxis(2); int nAxis2Bins = axis2->GetNbins();
  TAxis* axis3 = m_dPP->GetTAxis(3); int nAxis3Bins = axis3->GetNbins();

  // for canvas, since its tgraphasymmerrors.
  double x0 = axis3->GetXmin();
  double x1 = axis3->GetXmax(); 
  double y0 = m_dPhiWidthMin;
  double y1 = m_dPhiWidthMax;

  double pDx = 0.05;
  
  std::string yTitle = "|#Delta#phi| width";
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

      TCanvas cW( "c", "c", 800, 600 );
      styleTool->SetCStyleGraph
	( cW, x0, y0, x1, y1, gTitle.c_str() );

      TLegend legW( 0.25, 0.13, 0.9, 0.26 );
      styleTool->SetLegendStyle( &legW );
      legW.SetNColumns(2);

      // tag for widths canvas
      std::string hTagCW =
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

	std::string hDefaultName = "h_" + allUnfoldedName + "_" + hTag;
	std::string gDefaultName = "g_" + allUnfoldedName + "_" + hTag;
	
	TH1D* hDefault = static_cast<TH1D*>( fInDefault->Get( hDefaultName.c_str() ) );
	vHdef.push_back( hDefault );
	
	std::vector< double > pX;
	std::vector< double > eX( nAxis3Bins, 0.1 * hDefault->GetBinWidth(1) );
 
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
	  
	  double yDefault =   hDefault->GetBinContent( axis3Bin );
	  pX.push_back(       hDefault->GetBinCenter( axis3Bin ) );
	  
	  // loop over uncertainties
	  for( auto uc : v_uc ){
	    // skip uc = 0 (default)
	    if( !uc ){ continue; }
	    std::string hUncertaintyName = "h_" + allUnfoldedName + "_" + hTag;
	    TH1D* hUncertainty = static_cast<TH1D*>
	      ( mFinUC[uc]->Get( hUncertaintyName.c_str() ) );
	    vHunc.push_back( hUncertainty );

	    int sign = uc > 0 ? 1 : -1;

	    double yShifted = hUncertainty->GetBinContent( axis3Bin );

	    double uncertainty  = ( yDefault - yShifted )/yDefault;

	    std::cout << "++++" << uc << " ++++" << axis3Bin << " " << axis1Bin << " " << axis2Bin << " " 
		      << sign << " " << yShifted << " " << yDefault << " " << uncertainty << std::endl;

	    // for JER negative is same as positive
	    if( uc  == 20 ){
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

	  
	  std::cout << hDefaultName << " "
		    << uncertaintyFinalYP << " "
		    << uncertaintyFinalYN << std::endl;
	  
	  pY .push_back( yDefault );
	  eYP.push_back( yDefault * uncertaintyFinalYP );
	  eYN.push_back( yDefault * uncertaintyFinalYN );
	} // end loop over axis3

	std::string gSystematicsName = "g_" + allSystematicsName + "_" + hTag;

	TGraphAsymmErrors* gDefault     = new TGraphAsymmErrors( hDefault );
	TGraphAsymmErrors* gSystematics = new TGraphAsymmErrors
	  ( nAxis3Bins, &(pX[0]), &(pY[0]), &(eX[0]), &(eX[0]), &(eYN[0]), &(eYP[0]) );

	// add some displacement along x
	double* xDef      = gDefault->GetX();
	double* eXDefLow  = gDefault->GetEXlow();
	double* eXDefHigh = gDefault->GetEXhigh();
	double* xSys      = gSystematics->GetX();
	for( int iX = 0; iX < nAxis3Bins; iX++ ){
	  *(      xDef + iX ) += style * pDx;
	  *(  eXDefLow + iX ) = 0;
	  *( eXDefHigh + iX ) = 0;
	  *(      xSys + iX ) += style * pDx;
	}
	
	styleTool->SetHStyle( gDefault    , style );
	styleTool->SetHStyle( gSystematics, style );
	style++;
	
	gDefault    ->SetName( gDefaultName.c_str() );
	gSystematics->SetName( gSystematicsName.c_str() );

	vGw.push_back( gDefault );
	vGw.push_back( gSystematics );
	
	gSystematics->SetTitle("");
       	gSystematics->GetXaxis()->SetRangeUser( x0, x1 );
       	gSystematics->GetYaxis()->SetRangeUser( y0, y1 );
	
	gDefault->SetTitle("");
       	gDefault->GetXaxis()->SetRangeUser( x0, x1 );
       	gDefault->GetYaxis()->SetRangeUser( y0, y1 );

	legW.AddEntry
	    ( gDefault,
	      anaTool->GetLabel( axis2Low, axis2Up, m_dPP->GetAxisLabel(2) ).c_str() );	
	
	// draw systematics first
	gSystematics->Draw("2");
	gDefault->Draw("p");
	
	gSystematics->Write();
	gDefault->Write();
      } // end loop over axis2

      // Draw the final canvas with all of the graphs.
      legW.Draw();
      
      DrawTopLeftLabels
	( m_dPP, axis0Low, axis0Up, axis1Low, axis1Up, 0, 0, 0, 0);

      DrawAtlasRight();

      SaveAsAll( cW, Form("h_%s_%s_%s", m_dPhiName.c_str(),
			  m_sFinal.c_str(), hTagCW.c_str()) ); 
    } // end loop over axis1     
  } // end loop over axis0

  for( auto& h : vHdef ){ delete h; }
  for( auto& g : vGw   ){ delete g; }
  
  std::cout << "DONE! Closing " << fOut->GetName() << std::endl;
  fOut->Close(); delete fOut;
  std::cout << "......Closed  " << std::endl;
}

//---------------------------
//        Drawing
//---------------------------

void DiJetAnalysisData::DrawAtlasRight( double x0, double y0, double scale )
{ drawTool->DrawAtlasInternalDataRight( x0, y0, m_is_pPb, scale ); } 

void DiJetAnalysisData::DrawAtlasRightBoth( double x0, double y0, double scale ){
  drawTool->DrawAtlasInternal( scale );
  drawTool->DrawRightLatex( 0.88, 0.87, "Data" );	  	  
} 
