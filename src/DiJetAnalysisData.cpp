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

#include <TLine.h>

#include <iostream>
#include <cmath>

#include <boost/assign.hpp>
#include <boost/format.hpp>

#include "MyRoot.h"

#include "DiJetAnalysisData.h"
#include "DeltaPhiProj.h"

DiJetAnalysisData::DiJetAnalysisData() : DiJetAnalysisData( false )
{}

DiJetAnalysisData::DiJetAnalysisData( bool is_pPb )
  : DiJetAnalysis( is_pPb, true )
{
  //========== Set Histogram Binning =============

  //==================== Cuts ====================    
}

DiJetAnalysisData::~DiJetAnalysisData(){}

void DiJetAnalysisData::Initialize(){
  DiJetAnalysis::Initialize();
  
  m_fNameIn = m_is_pPb ?
    "/home/yakov/Projects/atlas/data/pPb.root" :
    "/home/yakov/Projects/atlas/data/pp.root"  ;

  m_mbTriggerName = "";
  m_allName       = "All";
  
  m_mbTriggerI         = -1;
  m_lowestCentTriggerI = -1;
  
  m_centMbCorrection   = 0;
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

  PlotEtaPhiPtMap( m_vHtriggerEtaPhiMap );
  PlotEtaPhiPtMap( m_vHtriggerEtaPtMap  );

  PlotSpectra( m_vHtriggerEtaSpect, "spect" );

  PlotEfficiencies( m_vHtriggerEtaSpectSim, m_vHtriggerEtaSpectDenom, "eff" );

  PlotDeltaPhi(  m_vHtriggerDphi, m_vHtriggerDphiNent, m_vTriggers );
  
  std::cout << "DONE! Closing " << cfNameOut << std::endl;
  m_fOut->Close();
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
  
  // This is an artificial one. We never loop over m_vTriggers
  // unless reading / writing histos.
  // So its ok to put this at the end, just as long as it gets
  // included in the writing/reading of histos.
  // this does not alter m_nTriggers
  m_vTriggers.push_back("All");
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
      ( new TH2D( Form("h_etaSpect_%s", trigger.c_str() ), 
		  ";#eta;#it{p}_{T} [GeV];",
		  m_nVarFwdEtaBins, 0, 1,
		  m_nPtSpectBins,
		  m_ptSpectMin, m_ptSpectMax ) ) ;
    m_vHtriggerEtaSpect.back()->GetXaxis()->
      Set( m_nVarFwdEtaBins, &( m_varFwdEtaBinning[0] ) );
    AddHistogram( m_vHtriggerEtaSpect.back() );

    if( !trigger.compare( m_allName ) )
      { m_hAllEtaSpect = m_vHtriggerEtaSpect.back(); }

    // ----- efficiencies ----
    m_vHtriggerEtaSpectSim.push_back
      ( new TH2D ( Form("h_etaSpectSim_%s", trigger.c_str() ), 
		   ";#eta;#it{p}_{T} [GeV]",
		   m_nVarFwdEtaBins, 0, 1,
		   m_nPtSpectBins,
		   m_ptSpectMin, m_ptSpectMax ) ) ;
    m_vHtriggerEtaSpectSim.back()->GetXaxis()->
      Set( m_nVarFwdEtaBins, &( m_varFwdEtaBinning[0] ) );
    AddHistogram( m_vHtriggerEtaSpectSim.back() );

    m_vHtriggerEtaSpectDenom.push_back
      ( new TH2D( Form("h_etaSpectDenom_%s", trigger.c_str() ), 
		  ";#eta;#it{p}_{T} [GeV];",
		  m_nVarFwdEtaBins, 0, 1,
		  m_nPtSpectBins,
		  m_ptSpectMin, m_ptSpectMax ) ) ;
    m_vHtriggerEtaSpectDenom.back()->GetXaxis()->
      Set( m_nVarFwdEtaBins, &( m_varFwdEtaBinning[0] ) );
    AddHistogram( m_vHtriggerEtaSpectDenom.back() );
    
    // -------- dPhi --------
    m_nDphiDim     = m_nDphiBins.size();
   
    THnSparse* hn =
      new THnSparseD( Form("hn_dPhi_%s", trigger.c_str() ), "",
		      m_nDphiDim, &m_nDphiBins[0],
		      &m_dPhiMin[0], &m_dPhiMax[0] );
    hn->GetAxis(0)->Set( m_nVarYstarBinsA, &( m_varYstarBinningA[0] ) );
    hn->GetAxis(1)->Set( m_nVarYstarBinsB, &( m_varYstarBinningB[0] ) );
    hn->GetAxis(2)->Set( m_nVarPtBins , &( m_varPtBinning[0]  ) );
    hn->GetAxis(3)->Set( m_nVarPtBins , &( m_varPtBinning[0]  ) );

    m_vHtriggerDphi.push_back( hn );
    AddHistogram( hn );
    
    THnSparse* hnNent =
      new THnSparseD( Form("hn_dPhiNent_%s", trigger.c_str() ), "",
		      m_nDphiDim, &m_nDphiBins[0],
		      &m_dPhiMin[0], &m_dPhiMax[0] );
    hnNent->GetAxis(0)->Set( m_nVarYstarBinsA, &( m_varYstarBinningA[0] ) );
    hnNent->GetAxis(1)->Set( m_nVarYstarBinsB, &( m_varYstarBinningB[0] ) );
    hnNent->GetAxis(2)->Set( m_nVarPtBins , &( m_varPtBinning[0]  ) );
    hnNent->GetAxis(3)->Set( m_nVarPtBins , &( m_varPtBinning[0]  ) );
    
    m_vHtriggerDphiNent.push_back( hnNent );
    AddHistogram( hnNent );
    
    if( !trigger.compare( m_allName ) ){
      m_hAllDphi     = hn;
      m_hAllDphiNent = hnNent;
    }
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

  // event loop
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
      
      /*
      // dPhi for trigger without matching
      AnalyzeDeltaPhi
	( m_vHtriggerDphi[iG], m_vHtriggerDphiNent[iG],
	  vR_jets, 1 );
      */
      
      const TLorentzVector* jet1 = NULL; const TLorentzVector* jet2 = NULL;
      // dPhi for all triggers, matched
      if( GetDiJets( vR_jets, jet1, jet2 ) ){
	if( JetInTrigRange( *jet1, iG ) ) {
	    AnalyzeDeltaPhi( m_hAllDphi, m_hAllDphiNent,
			     vR_jets, m_vTriggersPrescale[iG]);
	    AnalyzeDeltaPhi( m_vHtriggerDphi[iG], m_vHtriggerDphiNent[iG],
			     vR_jets, 1 );
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
	
	m_hAllEtaSpect->Fill( jetEtaAdj, jetPt, m_vTriggersPrescale[iG] );
      } // end loop over jets
    } // end loop over iG
    // EFFICIENCIES - only for good run and LBN
    if( goodRunLBN ){ AnalyzeEff( vR_jets, vTrig_jets, mTriggerFired ); }    
  } // end event loop
  std::cout << "DONE! Has: " << nEventsTotal << " events." << std::endl;

  m_fIn->Close();
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

void DiJetAnalysisData::CleanEfficiency( TGraphAsymmErrors* g, int iG ){
  double x, y;
  for( int i = 0; i < g->GetN(); i++ ){
    if( g->GetPoint( i, x, y ) == -1 ){ continue; }
    if( x > ( m_vTriggersTholdPt[iG] - 15 ) ){ continue; }
    g->SetPoint     ( i, x, 0 );
    g->SetPointError( i, 0, 0, 0, 0 );
  }
}

//---------------------------------
//            Plotting
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
	( m_fIn->Get( Form("h_etaSpect_%s", trigger.c_str() ))));
    m_vHtriggerEtaSpect.back()->SetDirectory(0);
    
    if( !trigger.compare( m_allName ) )
      { m_hAllEtaSpect = m_vHtriggerEtaSpect.back(); }

    // ----- efficiencies ----
    m_vHtriggerEtaSpectSim.push_back
      ( static_cast< TH2D *>
	( m_fIn->Get( Form("h_etaSpectSim_%s", trigger.c_str() ))));
    m_vHtriggerEtaSpectSim.back()->SetDirectory(0);

    m_vHtriggerEtaSpectDenom.push_back
      ( static_cast< TH2D *>
	( m_fIn->Get( Form("h_etaSpectDenom_%s", trigger.c_str() ))));
    m_vHtriggerEtaSpectDenom.back()->SetDirectory(0);
    // -------- dPhi- --------
    m_vHtriggerDphi.push_back
      ( static_cast< THnSparse *>
	( m_fIn->Get( Form("hn_dPhi_%s", trigger.c_str() ))));

    m_vHtriggerDphiNent.push_back
      ( static_cast< THnSparse *>
	( m_fIn->Get( Form("hn_dPhiNent_%s", trigger.c_str() ))));
    
    if( !trigger.compare(m_allName) ){
      m_hAllDphi     = m_vHtriggerDphi.back();
      m_hAllDphiNent = m_vHtriggerDphiNent.back();
    }
  }
  
  m_fIn->Close();
}

void DiJetAnalysisData::PlotSpectra( std::vector< TH2* >& vTrigSpect,
				     const std::string& type){
  std::string yAxisTitle = "dN/d#it{p}_{T}";

  double ptSpectWidth =
    ( m_ptSpectMax - m_ptSpectMin ) / m_nPtSpectBins;

  // use this as reference because
  // it should be in every file
  int nXbins = m_hAllEtaSpect->GetNbinsX();
  
  std::vector< std::vector< TH1* > > vSpect;
  vSpect.resize( m_nTriggers + 1 );
  
  double max = -1;
  for( uint iG = 0; iG < m_nTriggers + 1; iG++){
    std::string trigger = m_vTriggers[iG];

    for( int xBin = 1; xBin <= nXbins; xBin++ ){
      double etaMin = m_hAllEtaSpect->GetXaxis()->GetBinLowEdge( xBin );
      double etaMax = m_hAllEtaSpect->GetXaxis()->GetBinUpEdge ( xBin );

      std::string etaName = anaTool->GetName( etaMin, etaMax, "Eta");
      
      TH1* h_etaSpect =
	vTrigSpect[iG]->
	ProjectionY( Form("h_%s_%s_%s",
			  type.c_str(),
			  etaName.c_str(),
			  trigger.c_str() ),
		     xBin, xBin );
      h_etaSpect->SetTitle
	( anaTool->GetEtaLabel( etaMin, etaMax, m_is_pPb ).c_str() );
      h_etaSpect->GetYaxis()->SetTitle( yAxisTitle.c_str() );
      h_etaSpect->Scale( 1./ptSpectWidth );
      vSpect[iG].push_back( h_etaSpect );
      
      if( max < h_etaSpect->GetMaximum() )
	{ max = h_etaSpect->GetMaximum(); }
    } // end loop over xBin
  } // end loop over iG

  // set maxima globally for all spectra hists.
  // easier to compare. Set on log scale.
  double power = log10(max);
  power = std::ceil(power);
  max = pow( 10, power );
  
  //------------------------------------------------
  //------- Draw Eta as Fucntion of Triggers -------
  //------------------------------------------------
  double lX0, lY0, lX1, lY1;
  
  if( m_is_pPb ){ lX0 = 0.70; lY0 = 0.54; lX1 = 0.85; lY1 = 0.71; }
  else          { lX0 = 0.20; lY0 = 0.23; lX1 = 0.47; lY1 = 0.40; }
  
  for( uint iG = 0; iG < m_nTriggers + 1; iG++ ){
    std::string trigger = m_vTriggers[iG];

    if( trigger.compare( m_allName ) ){ continue; };
    
    std::string cName  = trigger;
    std::string cLabel = trigger;

    TCanvas c( "c", cLabel.c_str(), 800, 600 );
    c.SetLogy();
    
    TLegend leg( lX0, lY0, lX1, lY1 );
    styleTool->SetLegendStyle( &leg, 0.65 );

    int style = 0;
    for( int iX = 0; iX < nXbins; iX++ ){
      int         xBin = iX + 1;
      double etaCenter = m_hAllEtaSpect->GetXaxis()->GetBinCenter ( xBin );

      // for pPb, dont draw at anything above -3.2
      if( m_is_pPb && etaCenter > -constants::FETAMIN ){ continue; }

      // for pp, dont draw central triggers below -3.2
      // or forward triggers above -3.2
      // for mb, and total draw everything
      bool isMb  = !trigger.compare(m_mbTriggerName) ? true : false;
      bool isAll = !trigger.compare(m_allName)       ? true : false;
      
      if( !m_is_pPb && etaCenter < -constants::FETAMIN &&
	  trigger.find("320eta490") == std::string::npos &&
	  !isMb && !isAll ){ continue; }
      else if( !m_is_pPb && etaCenter > -constants::FETAMIN &&
	       trigger.find("320eta490") != std::string::npos &&
	       !isMb && !isAll ){ continue; }

      TH1* h = vSpect[iG][iX];
      styleTool->SetHStyle( h, style++ );
      h->Draw("epsame");
      h->SetMinimum( 1 );
      h->SetMaximum( max );
      leg.AddEntry( h, h->GetTitle() );
      h->SetTitle("");
    } // end loop over iX
    
    leg.Draw("same");
    drawTool->DrawAtlasInternalDataRight( 0, 0, m_is_pPb ); 
    drawTool->DrawRightLatex( 0.88, 0.76, cLabel );
    SaveAsAll( c, type, cName );
  } // end loop over iG
  
  //------------------------------------------------
  //------- Draw Triggers as Fucntion of Eta -------
  //------------------------------------------------
  if( m_is_pPb ){ lX0 = 0.45; lY0 = 0.54; lX1 = 0.76; lY1 = 0.67; }
  else          { lX0 = 0.35; lY0 = 0.67; lX1 = 0.65; lY1 = 0.88; }
  
  for( int iX = 0; iX < nXbins; iX++ ){
    // no longer have title.
    // Need to make output files
    // with eta names
    int      xBin = iX + 1;
    double etaMin = m_hAllEtaSpect->GetXaxis()->GetBinLowEdge( xBin );
    double etaMax = m_hAllEtaSpect->GetXaxis()->GetBinUpEdge ( xBin );
    double etaCenter = m_hAllEtaSpect->GetXaxis()->GetBinCenter ( xBin );

    // for pPb, dont draw at anything above -3.2
    if( m_is_pPb && etaCenter > -constants::FETAMIN ){ continue; }
     
    std::string cName  = anaTool->GetName( etaMin, etaMax, "Eta" );
    std::string cLabel = anaTool->GetEtaLabel( etaMin, etaMax, m_is_pPb );
    
    TCanvas c( "c", cLabel.c_str(), 800, 600 );
    c.SetLogy();
    
    TLegend leg( lX0, lY0, lX1, lY1 );
    styleTool->SetLegendStyle( &leg  );
    
    int style = 1;
    for( uint iG = 0; iG < m_nTriggers + 1; iG++ ){
      std::string trigger = m_vTriggers[iG];
      
      // for pp, dont draw central triggers below -3.2
      // or forward triggers above -3.2
      // for mb, and total draw everything
      if( !m_is_pPb && etaCenter < -constants::FETAMIN &&
	  trigger.find("320eta490") == std::string::npos &&
	  trigger.compare(m_mbTriggerName) &&
	  trigger.compare(m_allName) ){ continue; }
      else if( !m_is_pPb && etaCenter > -constants::FETAMIN &&
	       trigger.find("320eta490") != std::string::npos &&
	       trigger.compare(m_mbTriggerName) &&
	       trigger.compare(m_allName) ){ continue; }
      
      TH1* h = vSpect[iG][iX];
      if( iG == m_nTriggers ){ style = 0; }
      styleTool->SetHStyle( h, style++ );
      h->Draw("epsame");
      h->SetMinimum( 1 );
      h->SetMaximum( max );
      leg.AddEntry( h, m_vTriggers[iG].c_str() );
    } // end loop over iG
    
    leg.Draw("same");
    drawTool->DrawAtlasInternalDataRight( 0, 0, m_is_pPb ); 
    drawTool->DrawRightLatex( 0.88, 0.76, cLabel );

    SaveAsAll( c, type, cName );
  } // end loop over iX

  // delete
  for( uint iG = 0; iG < m_nTriggers; iG++ ){
    for( int iX = 0; iX < nXbins; iX++ )
      { delete vSpect[iG][iX]; }
  }
}

void DiJetAnalysisData::PlotEfficiencies( std::vector< TH2* >& vTrigSpect,
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
    
    drawTool->DrawAtlasInternalDataRight( 0, 0, m_is_pPb ); 
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
    
    drawTool->DrawAtlasInternalDataRight( 0, 0, m_is_pPb ); 
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

void DiJetAnalysisData::PlotEtaPhiPtMap( std::vector< TH2* >& vTrigHin ){
  TCanvas c_map("c_map","c_map",800,600);
  
  for( uint iG = 0; iG < m_nTriggers; iG++){
    TH2* h = vTrigHin[iG];
    h->Draw("col");
    styleTool->SetHStyle( h, 0 );
    drawTool->DrawAtlasInternalDataRight
      ( 0, -0.55, m_is_pPb );  
    
    SaveAsAll( c_map, h->GetName() );
  }
}

void DiJetAnalysisData::PlotDphiTogether(){
  // Check if the directories exist.
  // If they don't, create them
  std::string outDir = "output";
  anaTool->CheckWriteDir( outDir.c_str() );
  outDir += "/all";
  anaTool->CheckWriteDir( outDir.c_str() );
  outDir += "/data";
  anaTool->CheckWriteDir( outDir.c_str() );

  TAxis* axis0 = m_dPP->GetTAxis( 0 );
  TAxis* axis1 = m_dPP->GetTAxis( 1 );
  TAxis* axis2 = m_dPP->GetTAxis( 2 );
  TAxis* axis3 = m_dPP->GetTAxis( 3 );
  
  int nAxis0Bins = axis0->GetNbins();
  int nAxis1Bins = axis1->GetNbins();
  int nAxis2Bins = axis2->GetNbins();
  int nAxis3Bins = axis3->GetNbins();
 
  std::string trigger_pPb = m_allName;
  std::string trigger_pp  = m_allName; 
  
  TFile* fIn_pPb = TFile::Open("output/output_pPb_data/c_myOut_pPb_data.root");
  TFile* fIn_pp  = TFile::Open("output/output_pp_data/c_myOut_pp_data.root");
  TFile* fOut    = new TFile  ("output/all/data/c_myOut_data.root","recreate");
  
  for( int axis0Bin = 1; axis0Bin <= nAxis0Bins; axis0Bin++ ){
     // set ranges
    double axis0Low, axis0Up;
    anaTool->GetBinRange
      ( axis0, axis0Bin, axis0Bin, axis0Low, axis0Up );
    
    for( int axis1Bin = 1; axis1Bin <= nAxis1Bins; axis1Bin++ ){
      double axis1Low, axis1Up;
      anaTool->GetBinRange
	( axis1, axis1Bin, axis1Bin, axis1Low, axis1Up );

      // get widths canvases
      std::string hTagCW =
	Form ("%s_%s",
	      	anaTool->GetName( axis0Low, axis0Up, m_dPP->GetAxisName(0) ).c_str(),
		anaTool->GetName( axis1Low, axis1Up, m_dPP->GetAxisName(1) ).c_str() );
    
      // Make canvas+leg for widths
      TCanvas cW("cW","cW", 800, 600 );

      TLegend legW( 0.33, 0.13, 0.87, 0.26 );
      styleTool->SetLegendStyle( &legW );
      legW.SetNColumns(2);
      
      int style = 0;

      std::vector< TH1* > vHw;
      
      for( int axis2Bin = 1; axis2Bin <= nAxis2Bins; axis2Bin++ ){
	double axis2Low , axis2Up;
	anaTool->GetBinRange
	  ( axis2, axis2Bin, axis2Bin, axis2Low, axis2Up );

	// switch back to cW canvas
	// because a new one was creted in axis3 loop
	cW.cd();
	
	// get widths histos
	std::string hTagW =
	  Form("%s_%s_%s",
	       anaTool->GetName( axis0Low, axis0Up, m_dPP->GetAxisName(0) ).c_str(),
	       anaTool->GetName( axis1Low, axis1Up, m_dPP->GetAxisName(1) ).c_str(),
	       anaTool->GetName( axis2Low, axis2Up, m_dPP->GetAxisName(2) ).c_str() );
	
	std::string hNameW_pPb = Form("h_dPhi_%s_%s", hTagW.c_str(), trigger_pPb.c_str() );
	std::string hNameW_pp  = Form("h_dPhi_%s_%s", hTagW.c_str(), trigger_pp.c_str() );

	TH1* hW_pPb = static_cast<TH1D*>( fIn_pPb->Get( hNameW_pPb.c_str() ) );
	TH1* hW_pp  = static_cast<TH1D*>( fIn_pp ->Get( hNameW_pp. c_str() ) );
	styleTool->SetHStyle( hW_pPb, style );
	styleTool->SetHStyle( hW_pp , style + 5 );
	hW_pPb->SetMarkerSize( hW_pPb->GetMarkerSize() * 1.5 );
	hW_pp-> SetMarkerSize( hW_pp-> GetMarkerSize() * 1.5 );
	vHw.push_back( hW_pPb ); vHw.push_back( hW_pp  );
	
	style++;

	if( hW_pPb->GetMean() ){
	  legW.AddEntry
	    ( hW_pPb, Form( "p+Pb %s", anaTool->GetLabel
			    ( axis2Low, axis2Up, m_dPP->GetAxisLabel(2) ).c_str() ) );	
	  hW_pPb->Draw("ep same X0");
	}

	if( hW_pp->GetMean() ){
	  legW.AddEntry
	    ( hW_pp , Form( "pp %s" , anaTool->GetLabel
			    ( axis2Low, axis2Up , m_dPP->GetAxisLabel(2) ).c_str() ) );	
	  hW_pp->Draw("ep same X0");
	}

	for( int axis3Bin = 1; axis3Bin <= nAxis3Bins; axis3Bin++ ){
          double axis3Low , axis3Up;
	  anaTool->GetBinRange
	    ( axis3, axis3Bin, axis3Bin, axis3Low, axis3Up );
	    
	  std::string hTag =
	    Form("%s_%s_%s_%s",
		 anaTool->GetName( axis0Low, axis0Up, m_dPP->GetAxisName(0) ).c_str(),
		 anaTool->GetName( axis1Low, axis1Up, m_dPP->GetAxisName(1) ).c_str(),
		 anaTool->GetName( axis2Low, axis2Up, m_dPP->GetAxisName(2) ).c_str(),
		 anaTool->GetName( axis3Low, axis3Up, m_dPP->GetAxisName(3) ).c_str() );
	  
	  std::string hName_pPb = Form("h_dPhi_%s_%s", hTag.c_str(), trigger_pPb.c_str() );
	  std::string hName_pp  = Form("h_dPhi_%s_%s", hTag.c_str(), trigger_pp.c_str() );
	    
	  TH1* h_pPb = static_cast<TH1D*>( fIn_pPb->Get( hName_pPb.c_str() ) );
	  TH1* h_pp  = static_cast<TH1D*>( fIn_pp ->Get( hName_pp.c_str() ) );
	  styleTool->SetHStyle( h_pPb, 0 );
	  styleTool->SetHStyle( h_pp , 1 );
	  
	  TF1* f_pPb = static_cast<TF1*>( fIn_pPb->Get( Form("f_%s", hName_pp.c_str())));
	  TF1* f_pp  = static_cast<TF1*>( fIn_pp ->Get( Form("f_%s", hName_pp.c_str())));
	  styleTool->SetHStyle( f_pPb, 0 );
	  styleTool->SetHStyle( f_pp , 1 );
	  f_pPb->SetLineColor( h_pPb->GetLineColor() );
	  f_pp->SetLineColor( h_pp->GetLineColor() );

	  double chi2NDF_pPb = f_pPb->GetChisquare()/f_pPb->GetNDF();
	  double chi2NDF_pp  = f_pp ->GetChisquare()/f_pp ->GetNDF();
	  
	  TCanvas c("c","c", 800, 600 );
	    
	  TLegend leg( 0.27, 0.41, 0.38, 0.52 );
	  styleTool->SetLegendStyle( &leg , 0.85 );

	  bool save = false;
	  
	  if( h_pPb->GetEntries() ){
	    h_pPb->SetMinimum(0);
	    h_pPb->GetXaxis()->SetRangeUser( 2, constants::PI );
	    leg.AddEntry( h_pPb, Form("p+Pb #Chi^{2}/NDF=%4.2f", chi2NDF_pPb ) );
	    h_pPb->Draw("epsame");
	    f_pPb->Draw("same");
	    save = true;
	  }

	  if( h_pp->GetEntries() ){
	    h_pp->SetMinimum(0);
	    h_pp->GetXaxis()->SetRangeUser( 2, constants::PI );
	    leg.AddEntry( h_pp ,  Form("pp #Chi^{2}/NDF=%4.2f", chi2NDF_pp ) );
	    h_pp->Draw("epsame");
	    f_pp->Draw("same");
	    save = true;
	  }

	  leg.Draw("same");

	  if( h_pPb->GetMaximum() > h_pp->GetMaximum() ){
	    h_pPb->SetMaximum( h_pPb->GetMaximum() * 1.1 );
	    h_pp->SetMaximum ( h_pPb->GetMaximum() * 1.1 );
	  } else {
	    h_pPb->SetMaximum( h_pp->GetMaximum() * 1.1 );
	    h_pp->SetMaximum ( h_pp->GetMaximum() * 1.1 );
	  }

	  drawTool->DrawTopLeftLabels
	    ( m_dPP, axis0Low, axis0Up, axis1Low, axis1Up,
	      axis2Low, axis2Up, axis3Low, axis3Up, 0.8 );

	  drawTool->DrawRightLatex( 0.88, 0.87, "Data", 0.8 );	  
	  drawTool->DrawAtlasInternal();

	  if( save ){ c.SaveAs( Form("output/all/data/h_dPhi_%s_data.png", hTag.c_str() )); }
	  // if( save ){ c.SaveAs( Form("output/all/data/h_dPhi_%s_data.pdf", hTag.c_str() )); }
	  SaveAsROOT( c, Form("h_dPhi_%s", hTag.c_str() ) );
	  
	  delete  f_pPb; delete  f_pp;
	  delete  h_pPb; delete  h_pp;
	} // end loop over axis3
      } // end loop over axis2

      // back to cW canvas
      cW.cd();

      legW.Draw("same");

      drawTool->DrawTopLeftLabels
	( m_dPP, axis0Low, axis0Up, axis1Low, axis1Up,
	  0, 0, 0, 0, 0.8 );
      
      drawTool->DrawRightLatex( 0.88, 0.87, "Data", 0.8 );     
      drawTool->DrawAtlasInternal();

      cW.SaveAs( Form("output/all/data/h_dPhi_%s_data.png", hTagCW.c_str() ));
      // cW.SaveAs( Form("output/all/data/h_dPhi_%s_data.pdf", hTagCW.c_str() ));
      SaveAsROOT( cW, Form("h_dPhi_%s", hTagCW.c_str() ) );

      for( auto & hW : vHw ){ delete hW; }
     } // end loop over ystar2
  } // end loop over ystar2
  
  std::cout << "DONE! Closing " << fOut->GetName() << std::endl;
  fOut->Close();
  std::cout << "......Closed  " << fOut->GetName() << std::endl;
}
