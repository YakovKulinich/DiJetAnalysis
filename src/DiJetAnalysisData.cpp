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

#include <iostream>
#include <cmath>

#include <boost/assign.hpp>
#include <boost/format.hpp>

#include "DiJetAnalysisData.h"
#include "DeltaPhiProj.h"

DiJetAnalysisData::DiJetAnalysisData() : DiJetAnalysisData( false )
{}

DiJetAnalysisData::DiJetAnalysisData( bool is_pPb )
  : DiJetAnalysis( is_pPb, true )
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
  DiJetAnalysis::Initialize();
  
  m_fNameIn = m_is_pPb ?
    "/home/yakov/Projects/atlas/data/pPb.root" :
    "/home/yakov/Projects/atlas/data/pp.root"  ;
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

void DiJetAnalysisData::ProcessPlotHistos(){
  LoadTriggers();
  LoadHistograms();

  std::string cfNameOut = m_dirOut + "/c_myOut_" + m_labelOut + ".root";
  m_fOut = new TFile( cfNameOut.c_str(),"RECREATE");

  // MakeEtaPhiPtMap( m_vHtriggerEtaPhiMap );
  // MakeEtaPhiPtMap( m_vHtriggerEtaPtMap  );

  // add a trigger "all" to collection
  // rest of plots are for combined triggers
  m_vTriggers.push_back( m_allName );
  
  m_hAllEtaSpect = CombineSamples( m_vHtriggerEtaSpect, m_etaSpectName );
  MakeSpectra( m_vHtriggerEtaSpect, m_vTriggers, m_etaSpectName );

  MakeEfficiencies( m_vHtriggerEtaSpectSim, m_vHtriggerEtaSpectDenom, m_effName );
  
  m_hAllDphi = CombineSamples( m_vHtriggerDphi, m_dPhiName );
  MakeDeltaPhi( m_vHtriggerDphi, m_vTriggers, m_dPhiName );

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
  THnSparse* m_hAllDphiUnfolded =
    UnfoldDeltaPhi( m_hAllDphi, fInMC, m_dPhiUnfoldedName );
  m_vHDphiUnfolded.push_back( m_hAllDphiUnfolded );
  m_vLabelUnfolded.push_back( m_allName );
  
  std::cout << "DONE! Closing " << cfNameOut << std::endl;
  m_fOut->Close(); delete m_fOut;
  std::cout << "......Closed  " << cfNameOut << std::endl;
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
		  m_nPtSpectBins,
		  m_ptSpectMin, m_ptSpectMax ) ) ;
    m_vHtriggerEtaSpect.back()->GetXaxis()->
      Set( m_nVarFwdEtaBins, &( m_varFwdEtaBinning[0] ) );
    AddHistogram( m_vHtriggerEtaSpect.back() );

    // ----- efficiencies ----
    m_vHtriggerEtaSpectSim.push_back
      ( new TH2D ( Form("h_%sSim_%s", m_etaSpectName.c_str(), trigger.c_str() ), 
		   ";#eta;#it{p}_{T} [GeV]",
		   m_nVarFwdEtaBins, 0, 1,
		   m_nPtSpectBins,
		   m_ptSpectMin, m_ptSpectMax ) ) ;
    m_vHtriggerEtaSpectSim.back()->GetXaxis()->
      Set( m_nVarFwdEtaBins, &( m_varFwdEtaBinning[0] ) );
    AddHistogram( m_vHtriggerEtaSpectSim.back() );

    m_vHtriggerEtaSpectDenom.push_back
      ( new TH2D( Form("h_%sDenom_%s", m_etaSpectName.c_str(), trigger.c_str() ), 
		  ";#eta;#it{p}_{T} [GeV];",
		  m_nVarFwdEtaBins, 0, 1,
		  m_nPtSpectBins,
		  m_ptSpectMin, m_ptSpectMax ) ) ;
    m_vHtriggerEtaSpectDenom.back()->GetXaxis()->
      Set( m_nVarFwdEtaBins, &( m_varFwdEtaBinning[0] ) );
    AddHistogram( m_vHtriggerEtaSpectDenom.back() );
    
    // -------- dPhi --------
    m_nDphiDim    = m_nDphiBins.size();
   
    THnSparse* hn =
      new THnSparseD( Form("h_%s_%s", m_dPhiName.c_str(), trigger.c_str() ), "",
		      m_nDphiDim, &m_nDphiBins[0],
		      &m_dPhiMin[0], &m_dPhiMax[0] );
    hn->GetAxis(0)->Set( m_nVarYstarBinsA, &( m_varYstarBinningA[0] ) );
    hn->GetAxis(1)->Set( m_nVarYstarBinsB, &( m_varYstarBinningB[0] ) );
    hn->GetAxis(2)->Set( m_nVarPtBins    , &( m_varPtBinning[0]     ) );
    hn->GetAxis(3)->Set( m_nVarPtBins    , &( m_varPtBinning[0]     ) );
    hn->GetAxis(4)->Set( m_nVarDphiBins  , &( m_varDphiBinning[0]   ) );
    
    m_vHtriggerDphi.push_back( hn );
    AddHistogram( hn );
    
    THnSparse* hnNent =
      new THnSparseD( Form("h_%sNent_%s", m_dPhiName.c_str(), trigger.c_str() ), "",
		      m_nDphiDim, &m_nDphiBins[0],
		      &m_dPhiMin[0], &m_dPhiMax[0] );
    hnNent->GetAxis(0)->Set( m_nVarYstarBinsA, &( m_varYstarBinningA[0] ) );
    hnNent->GetAxis(1)->Set( m_nVarYstarBinsB, &( m_varYstarBinningB[0] ) );
    hnNent->GetAxis(2)->Set( m_nVarPtBins    , &( m_varPtBinning[0]     ) );
    hnNent->GetAxis(3)->Set( m_nVarPtBins    , &( m_varPtBinning[0]     ) );
    hnNent->GetAxis(4)->Set( m_nVarDphiBins  , &( m_varDphiBinning[0]   ) );
	
    m_vHtriggerDphiNent.push_back( hnNent );
    AddHistogram( hnNent );
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
  
  m_fIn = TFile::Open( m_fNameIn.c_str() );
  m_tree = static_cast<TTree*>( m_fIn->Get( "tree" ) );

  // Connect to tree
  m_tree->SetBranchAddress( "vTrig_jets"  , &p_vTrig_jets );
  m_tree->SetBranchAddress( "vR_C_jets"   , &p_vR_jets );
  m_tree->SetBranchAddress( "v_isCleanJet", &p_v_isCleanJet );

  vTriggerFired.resize   ( m_nTriggers );

  for( uint iG = 0; iG < m_nTriggers; iG++ ){
    m_tree->SetBranchAddress
      ( Form("passed_%s", m_vTriggers[iG].c_str() ), &mTriggerFired[iG] );
  }

  m_tree->SetBranchAddress( "runNumber", &runNumber );
  m_tree->SetBranchAddress( "LBN"      , &LBN );  

  std::vector< int > mTrigPassed; 
  std::vector< int > mTrigFailed;  
  mTrigPassed.resize( m_nTriggers );
  mTrigFailed.resize( m_nTriggers );

  std::vector< int > mTrigJetsForward;
  std::vector< int > mTrigJetsTotal;
  mTrigJetsForward.resize( m_nTriggers );
  mTrigJetsTotal.resize  ( m_nTriggers );
  
  // n events
  int nEventsTotal = m_tree->GetEntries();

  nEvents = nEvents > 0 ? nEvents : nEventsTotal;
   startEvent = startEvent < nEventsTotal ? startEvent : nEventsTotal - 1;
  
  int endEvent = startEvent + nEvents < nEventsTotal ?
					startEvent + nEvents : nEventsTotal;

  // -------- EVENT LOOP ---------
  for( m_ev = startEvent; m_ev < endEvent; m_ev++ ){
    m_tree->GetEntry( m_ev );

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
	  AnalyzeDeltaPhi( m_vHtriggerDphi[iG], m_vHtriggerDphiNent[iG], vR_jets );
	}
      }
  
      // loop over jets 
      for( auto& jet : vR_jets ){
	// cleaned jets have px, py, pz set to 0
	if( jet.Pt() == 0 ) continue;
	
	// ETA-PHI
	double jetEta = jet.Eta();
	double jetPhi = jet.Phi();
	double jetPt = jet.Pt()/1000.;
	
	m_vHtriggerEtaPhiMap[iG]->Fill( jetEta, jetPhi );
	m_vHtriggerEtaPtMap [iG]->Fill( jetEta, jetPt ); 

	// convert positive eta to negative because
	// in pp it doesnt matter.
	// our histos run in negative eta
	// (due to pPb configuration)
	// the labels will be taken care of so it is ok
	double jetEtaAdj = AdjustEtaForPP( jetEta );

	// fill spectra 
	m_vHtriggerEtaSpect[iG]->Fill( jetEtaAdj, jetPt );
	
	// check if the jet is in appropriate range
	// for the trigger fired
	if( !JetInTrigRange( jet, iG ) ){ continue; };
      } // end loop over jets
    } // end loop over iG
    // EFFICIENCIES - only for good run and LBN
    if( goodRunLBN ){ AnalyzeEff( vR_jets, vTrig_jets, mTriggerFired ); }    
  } // -------- END EVENT LOOP ---------
  std::cout << "DONE! Has: " << nEventsTotal << " events." << std::endl;

  m_fIn->Close(); delete m_fIn;
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

void DiJetAnalysisData::GetInfoBoth( std::string& outSuffix,
				     std::string& name_a   , std::string& name_b  ,
				     std::string& label_a  , std::string& label_b ,
				     std::string& suffix_a , std::string& suffix_b ){
  outSuffix = "data";
  name_a    = m_dPhiName + "_" + m_allName;
  name_b    = m_dPhiName + "_" + m_allName;
  label_a   = "#it{p}+Pb";
  label_b   = "#it{pp}";
  suffix_a  = "pPb_data";
  suffix_b  = "pp_data";
}

void DiJetAnalysisData::GetInfoUnfolding( std::string& measuredName,
					  std::string& measuredLabel ){
  measuredName  = m_dPhiName;
  measuredLabel = "Data";
}

//---------------------------------
//     Get Quantities / Plot 
//---------------------------------

void DiJetAnalysisData::LoadHistograms(){
  m_fIn = TFile::Open( m_rootFname.c_str() ); 

  for( auto& trigger : m_vTriggers ){
    // -------- maps ---------
    m_vHtriggerEtaPhiMap.push_back
      ( static_cast< TH2D* >
	( m_fIn-> Get( Form("h_etaPhiMap_%s",
			    trigger.c_str() ))));
    m_vHtriggerEtaPhiMap.back()->SetDirectory(0);
    
    m_vHtriggerEtaPtMap.push_back
      ( static_cast< TH2D* >
	( m_fIn->Get( Form("h_etaPtMap_%s", trigger.c_str() ))));
    m_vHtriggerEtaPtMap.back()->SetDirectory(0);
    
    // -------- spect --------
    m_vHtriggerEtaSpect.push_back
      ( static_cast< TH2D *>
	( m_fIn->Get
	  ( Form("h_%s_%s", m_etaSpectName.c_str(), trigger.c_str() ))));
    m_vHtriggerEtaSpect.back()->SetDirectory(0);
    
    // ----- efficiencies ----
    m_vHtriggerEtaSpectSim.push_back
      ( static_cast< TH2D *>
	( m_fIn->Get
	  ( Form("h_%sSim_%s", m_etaSpectName.c_str(), trigger.c_str() ))));
    m_vHtriggerEtaSpectSim.back()->SetDirectory(0);

    m_vHtriggerEtaSpectDenom.push_back
      ( static_cast< TH2D *>
	( m_fIn->Get
	  ( Form("h_%sDenom_%s", m_etaSpectName.c_str(), trigger.c_str() ))));
    m_vHtriggerEtaSpectDenom.back()->SetDirectory(0);
    // -------- dPhi- --------
    m_vHtriggerDphi.push_back
      ( static_cast< THnSparse *>
	( m_fIn->Get( Form("h_%s_%s", m_dPhiName.c_str(), trigger.c_str() ))));

    m_vHtriggerDphiNent.push_back
      ( static_cast< THnSparse *>
	( m_fIn->Get( Form("h_%sNent_%s", m_dPhiName.c_str(), trigger.c_str() ))));
  }
  
  m_fIn->Close(); delete m_fIn;
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
      double etaMin = m_hAllEtaSpect->GetXaxis()->GetBinLowEdge( xBin );
      double etaMax = m_hAllEtaSpect->GetXaxis()->GetBinUpEdge ( xBin );
            
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
    if( !trigger.compare( m_mbTriggerName ) )
      { continue; }
    if( !trigger.compare( m_allName ) )
      { continue; }

    std::string cName  = trigger;
    std::string cLabel = trigger;

    TCanvas c( "c", cLabel.c_str(), 800, 600 );
    styleTool->SetCStyleEff( c, xMin, m_effMin, xMax, m_effMax, gLabel );
    
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
    
    SaveAsAll( c, type, cName );
  } // end loop over iG

  //------------------------------------------------
  //---------- Draw Triggers in Eta Bins -----------
  //------------------------------------------------
  for( int iX = 0; iX < nXbins; iX++ ){
    // no longer have title.
    // Need to make output files
    // with eta names
    int         xBin = iX + 1;
    double    etaMin = m_hAllEtaSpect->GetXaxis()->GetBinLowEdge( xBin );
    double    etaMax = m_hAllEtaSpect->GetXaxis()->GetBinUpEdge ( xBin );
    double etaCenter = m_hAllEtaSpect->GetXaxis()->GetBinCenter ( xBin );
    
    // temporary, dont draw the 3.1->3.2 bin
    if( std::abs(etaCenter) < 3.2 && std::abs(etaCenter) > 3.1 ){ continue; }
    
    // for pPb, dont draw at anything above -3.2
    if( m_is_pPb && etaCenter > -constants::FETAMIN ){ continue; }
    
    std::string cName  = anaTool->GetName( etaMin, etaMax, "Eta" );
    std::string cLabel = anaTool->GetEtaLabel( etaMin, etaMax, m_is_pPb );
    
    TCanvas c( "c", cLabel.c_str(), 800, 600 );
    styleTool->SetCStyleEff( c, xMin, m_effMin, xMax, m_effMax, gLabel );
    
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

    SaveAsAll( c, type, cName );
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

//---------------------------
//        Drawing
//---------------------------

void DiJetAnalysisData::DrawAtlasRight( double x0, double y0, double scale )
{ drawTool->DrawAtlasInternalDataRight( x0, y0, m_is_pPb, scale ); } 

void DiJetAnalysisData::DrawAtlasRightBoth( double x0, double y0, double scale ){
  drawTool->DrawAtlasInternal( scale );
  drawTool->DrawRightLatex( 0.88, 0.87, "Data", 0.8 );	  	  
} 
