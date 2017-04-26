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

DiJetAnalysisData::DiJetAnalysisData() : DiJetAnalysisData( true, true )
{}

DiJetAnalysisData::DiJetAnalysisData
( bool isData, bool is_pPb )
  : DiJetAnalysis( isData, is_pPb )
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

  m_allName = "All";
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

  PlotEfficiencies( m_vHtriggerEtaSpectSim, m_vHtriggerEtaSpect, "eff" );
  
  PlotDeltaPhi( m_vHtriggerDphi, m_vTriggers );
  
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

  m_vTriggersTholdPt =  anaTool->vectoriseD
    ( GetConfig()->GetValue( Form("triggersThold.%s",triggerMenu.c_str()), "" ), " ");

  m_vTriggersEffPtLow = anaTool->vectoriseD
    ( GetConfig()->GetValue( Form("triggersEffLow.%s",triggerMenu.c_str()), "" ), " ");

  m_vTriggersEffPtHigh = anaTool->vectoriseD
    ( GetConfig()->GetValue( Form("triggersEffHigh.%s",triggerMenu.c_str()), "" ), " ");

  m_nTriggers = m_vTriggers.size();
  
  // See if we can find the min bias trigger
  for( uint iG = 0; iG < m_nTriggers; iG++ ){
    std::string trigger = m_vTriggers[iG];
    if( trigger.find("_mb_") != std::string::npos ){ 
      m_mbTriggerI    = iG;
      m_mbTriggerName = trigger;
      std::cout << "Found " << trigger << " at " << m_mbTriggerI << std::endl;
      break;
    }
  }

  // This is an artificial one. We never loop over m_vTriggers
  // unless reading / writing histos.
  // So its ok to put this at the end, just as long as it gets
  // included in the writing/reading of histos.
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
    
    m_vHtriggerEtaSpectSim.push_back
      ( new TH2D ( Form("h_etaSpectSim_%s", trigger.c_str() ), 
		   ";#eta;#it{p}_{T} [GeV]",
		   m_nVarFwdEtaBins, 0, 1,
		   m_nPtSpectBins,
		   m_ptSpectMin, m_ptSpectMax ) ) ;
    m_vHtriggerEtaSpectSim.back()->GetXaxis()->
      Set( m_nVarFwdEtaBins, &( m_varFwdEtaBinning[0] ) );
    AddHistogram( m_vHtriggerEtaSpectSim.back() );

    // -------- dPhi --------
    std::vector< int    > nDphiBins;
    std::vector< double > dPhiMin;
    std::vector< double > dPhiMax;

    boost::assign::push_back( nDphiBins )
      ( m_nVarYstarBinsA )( m_nVarYstarBinsB )
      ( m_nDphiPtBins    )( m_nDphiPtBins    )
      ( m_nDphiDphiBins  );
    
    boost::assign::push_back( dPhiMin  )
      (        0       )(       0      )
      ( m_nDphiPtMin   )( m_nDphiPtMin )
      ( m_nDphiDphiMin );

    boost::assign::push_back( dPhiMax  )
      (        1       )(       1      )
      ( m_nDphiPtMax   )( m_nDphiPtMax )
      ( m_nDphiDphiMax );

    uint nDim = nDphiBins.size();
    
    THnSparse* hn =
      new THnSparseD( Form("hn_dPhi_%s", trigger.c_str() ), "",
		      nDim, &nDphiBins[0], &dPhiMin[0], &dPhiMax[0] );
    hn->GetAxis(0)->Set( m_nVarYstarBinsA, &( m_varYstarBinningA[0] ) );
    hn->GetAxis(1)->Set( m_nVarYstarBinsB, &( m_varYstarBinningB[0] ) );
    hn->GetAxis(2)->Set( m_nVarPtBins , &( m_varPtBinning[0]  ) );
    hn->GetAxis(3)->Set( m_nVarPtBins , &( m_varPtBinning[0]  ) );

    m_vHtriggerDphi.push_back( hn );
    AddHistogram( hn );

    if( !trigger.compare( m_allName ) )
      { m_hAllDphi = hn; }
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
  std::vector< float > vTriggerPrescale;

  // because vector bool doesnt return lvalue
  std::map< int, bool > mTriggerFired; 

  int runNumber = 0;
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
  vTriggerPrescale.resize( m_nTriggers );

  for( uint iG = 0; iG < m_nTriggers; iG++ ){
    m_tree->SetBranchAddress
      ( Form("passed_%s", m_vTriggers[iG].c_str() ), &mTriggerFired[iG] );
    m_tree->SetBranchAddress
      ( Form("prescale_%s", m_vTriggers[iG].c_str() ),&vTriggerPrescale[iG] );
  }

  m_tree->SetBranchAddress( "runNumber", &runNumber );  

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
    if( runNumber == 312649 ){ continue; }

    ApplyCleaning ( vR_jets, v_isCleanJet );
    ApplyIsolation( 1.0, vR_jets );

    if( anaTool->DoPrint(m_ev) ) {
      std::cout << "\nEvent : " << m_ev << "    runN : " << runNumber
		<< "    has : " << vR_jets.size() << " jets" 
		<< "    and : " << vTrig_jets.size() << " trig jets" 
		<< std::endl; 
    }
    
    std::sort( vR_jets.begin(), vR_jets.end(),
	       anaTool->sortByDecendingPt );

    // SPECTRA AND ETA PHI PT MAPS
    // loop over passed triggers
    for( uint iG = 0; iG < m_nTriggers; iG++ ){

      // check if we have that trigger
      // if we dont - continue
      if( !mTriggerFired[iG] ) {
	mTrigFailed[iG]++;
	continue;	
      }

      // dPhi for trigger without matching
      AnalyzeDeltaPhi( m_vHtriggerDphi[iG], vR_jets );

      TLorentzVector jet1, jet2;
      // dPhi for all triggers, matched
      if( GetDiJets( vR_jets, jet1, jet2 ) ){
	double jet1_eta = jet1.Eta();
	double jet2_eta = jet2.Eta();

	if( ( IsForwardDetector( jet1_eta ) &&
	      IsInTriggerRange( jet1, iG ) ) ||
	    ( IsForwardDetector( jet2_eta ) &&
	      IsInTriggerRange( jet2, iG ) ) )
	  { AnalyzeDeltaPhi( m_hAllDphi, vR_jets,
			     vTriggerPrescale[iG]); }
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
	double jetEtaAdj = anaTool->AdjustEtaForPP( jetEta, m_is_pPb );
	
	if( jetPt >= m_vTriggersTholdPt[iG] )
	  { mTrigJetsTotal[iG] ++; }
	
	// saves computation time
	// can theoreticall leave out
	// will go to overflow bins in histos
	if( !IsForwardDetector( jetEtaAdj) ){ continue; }
	m_vHtriggerEtaSpect[iG]->Fill( jetEtaAdj, jetPt );

	// check if the jet is in appropriate range
	// for the trigger fired
	if( !IsInTriggerRange( jet, iG ) ){ continue; };
	
	m_hAllEtaSpect->Fill( jetEtaAdj, jetPt, vTriggerPrescale[iG] );
	mTrigJetsForward[iG]++;
      } // end loop over jets
      mTrigPassed[iG]++;
    } // end loop over iG
    // EFFICIENCIES
    AnalyzeEff( vTrig_jets, vR_jets, mTriggerFired );    
  } // end event loop
  std::cout << "DONE! Has: " << nEventsTotal << " events." << std::endl;

  for( unsigned iG = 0; iG < m_vTriggers.size(); iG++ ){
    std::cout <<  m_vTriggers[iG]
	      << "  has: " << mTrigJetsForward[iG]
	      << "  FWD jets, " << mTrigJetsTotal[iG]
	      << "  TOTAL jets, with " << mTrigPassed[iG]
	      << "  evPassed and "<< mTrigFailed[iG]
	      << "  failed. Total = "
	      << mTrigPassed[iG] + mTrigFailed[iG]
	      << std::endl;
  }

  m_fIn->Close();
}

//---------------------------
//       Analysis
//---------------------------
bool DiJetAnalysisData::IsInTriggerRange( TLorentzVector& jet, uint iG ){
  double jetPt  = jet.Pt()/1000.;
  double jetEta = jet.Eta();
  if( jetPt > m_vTriggersEffPtLow[iG] &&
      ( jetPt < m_vTriggersEffPtHigh[iG] || m_vTriggersEffPtHigh[iG] < 0 ) )
    { return true; }
  return false;
}

void DiJetAnalysisData::AnalyzeEff( std::vector< TLorentzVector >& vTrig_jets,
				    std::vector< TLorentzVector >& vR_jets,
				    std::map< int, bool >& mTriggerFired ){
  if( vTrig_jets.empty() ){ return; }

  // sort trigger jets
  std::sort( vTrig_jets.begin(), vTrig_jets.end(),
	     anaTool->sortByDecendingPt );

  // take highest pt trigger jet in forward eta range
  // for trigger jets, we make this at 3.2.
  // For reco jets we use 3.3.
  TLorentzVector* tJet = NULL;
  for( auto& jet : vTrig_jets ){
    if( anaTool->AdjustEtaForPP
	( jet.Eta(), m_is_pPb) > -constants::TETAMIN )
      { continue; }
    tJet = &jet;
    break; // found highest forward trigger jet
  }    

  // didnt find a trigger jet in forward eta range
  if( !tJet ) return;

  // now fill histos for triggers that could have passed
  for( uint iG = 0; iG < m_nTriggers; iG++ ){
    // check if this triggers reference trigger fired
    int refTrigIndex = m_vTriggersRefIndex[iG];
    if( !mTriggerFired[refTrigIndex] ){ continue; }
    
    double tJetPt = tJet->Pt()/1000;
   
    // now check if we should fill for that trigger
    // only if trigger jet is above threshold
    // for that trigger
    if( tJetPt < m_vTriggersTholdPt[iG] ){ continue; } 

    // fill jets for that trigger
    for( auto& jet : vR_jets ){
      // avoid cleaned jets
      if( jet.Pt() == 0 ) continue;

      double jetEtaAdj = anaTool->AdjustEtaForPP( jet.Eta(), m_is_pPb );

      // eta cut. want forward trigger jet
      if( !IsForwardDetector( jetEtaAdj ) ){ continue; }
      
      m_vHtriggerEtaSpectSim[iG]->
	Fill( jetEtaAdj, jet.Pt()/1000. );
    } // end loop over jets
  } // end loop over iG
}

//---------------------------------
//            Plot Data
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
    
    m_vHtriggerEtaSpectSim.push_back
      ( static_cast< TH2D *>
	( m_fIn->Get( Form("h_etaSpectSim_%s", trigger.c_str() ))));
    m_vHtriggerEtaSpectSim.back()->SetDirectory(0);

    // -------- dPhi- --------
    m_vHtriggerDphi.push_back
      ( static_cast< THnSparse *>
	( m_fIn->Get( Form("hn_dPhi_%s", trigger.c_str() ))));

    if( !trigger.compare(m_allName ) )
      { m_hAllDphi = m_vHtriggerDphi.back(); }
  }
  
  m_fIn->Close();
}

void DiJetAnalysisData::PlotSpectra( std::vector< TH2* >& vTrigSpect,
				     const std::string& type){
  std::string yAxisTitle = "dN/d#it{p}_{T}";

  double ptSpectWidth =
    ( m_ptSpectMax - m_ptSpectMin ) / m_nPtSpectBins;
	          
  double lX0, lY0, lX1, lY1;
  if( m_is_pPb ){ lX0 = 0.45; lY0 = 0.54; lX1 = 0.76; lY1 = 0.67; }
  else          { lX0 = 0.56; lY0 = 0.54; lX1 = 0.86; lY1 = 0.67; }

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

      std::string etaName = anaTool->GetName(etaMin,etaMax,"Eta");
      
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
  for( uint iG = 0; iG < m_nTriggers + 1; iG++ ){
    std::string cName  = m_vTriggers[iG];
    std::string cLabel = m_vTriggers[iG];

    TCanvas c( "c", cLabel.c_str(), 800, 600 );
    c.SetLogy();
    
    TLegend leg( 0.7, lY0, lX1, lY1 );
    styleTool->SetLegendStyle( &leg  );

    int style = 0;
    for( int iX = 0; iX < nXbins; iX++ ){
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
    drawTool->DrawLeftLatex( 0.15, 0.17, cLabel );
    SaveAsAll( c, type, cName );
  } // end loop over iG
  
  //------------------------------------------------
  //------- Draw Triggers as Fucntion of Eta -------
  //------------------------------------------------
  for( int iX = 0; iX < nXbins; iX++ ){
    // no longer have title.
    // Need to make output files
    // with eta names
    int      xBin = iX + 1;
    double etaMin = m_hAllEtaSpect->GetXaxis()->GetBinLowEdge( xBin );
    double etaMax = m_hAllEtaSpect->GetXaxis()->GetBinUpEdge ( xBin );
    
    std::string cName  = anaTool->GetName( etaMin, etaMax, "Eta" );
    std::string cLabel = anaTool->GetEtaLabel( etaMin, etaMax, m_is_pPb );
    
    TCanvas c( "c", cLabel.c_str(), 800, 600 );
    c.SetLogy();
    
    TLegend leg( lX0, lY0, lX1, lY1 );
    styleTool->SetLegendStyle( &leg  );
    
    int style = 1;
    for( uint iG = 0; iG < m_nTriggers + 1; iG++ ){
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
    drawTool->DrawLeftLatex( 0.15, 0.17, cLabel );

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
  double lX0 = 0.13;
  double lY0 = 0.75;
  double lX1 = 0.39;
  double lY1 = 0.87;

  // us m_hAllEtaSpect because its always there
  double xMin = m_hAllEtaSpect->GetYaxis()->GetXmin() - 10;
  double xMax = 60.; // m_hAllEtaSpect->GetYaxis()->GetXmax();

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
    // dont draw MB trigger
    if( !m_vTriggers[iG].compare( m_mbTriggerName ) )
      { continue; }
    if( !m_vTriggers[iG].compare( m_allName ) )
      { continue; }
    std::string cName  = m_vTriggers[iG];
    std::string cLabel = m_vTriggers[iG];

    TCanvas c( "c", cLabel.c_str(), 800, 600 );
    styleTool->SetCStyleEff( c, xMin, m_effMin, xMax, m_effMax, gLabel );
    
    TLegend leg( lX0, lY0, lX1, lY1 );
    styleTool->SetLegendStyle( &leg  );

    int style = 0;
    for( int iX = 0; iX < nXbins; iX++ ){
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
    drawTool->DrawLeftLatex( 0.4, 0.17, cLabel );
    
    SaveAsAll( c, type, cName );
  } // end loop over iG

  //------------------------------------------------
  //---------- Draw Triggers in Eta Bins -----------
  //------------------------------------------------
  for( int iX = 0; iX < nXbins; iX++ ){
    // no longer have title.
    // Need to make output files
    // with eta names
    int      xBin = iX + 1;
    double etaMin = m_hAllEtaSpect->GetXaxis()->GetBinLowEdge( xBin );
    double etaMax = m_hAllEtaSpect->GetXaxis()->GetBinUpEdge ( xBin );
    
    std::string cName  = anaTool->GetName( etaMin, etaMax, "Eta" );
    std::string cLabel = anaTool->GetEtaLabel( etaMin, etaMax, m_is_pPb );
    
    TCanvas c( "c", cLabel.c_str(), 800, 600 );
    styleTool->SetCStyleEff( c, xMin, m_effMin, xMax, m_effMax, gLabel );
    
    TLegend leg( lX0, lY0, lX1, lY1 );
    styleTool->SetLegendStyle( &leg  );
    
    int style = 0;
    for( uint iG = 0; iG < m_nTriggers; iG++ ){
      // dont draw MB trigger
      if( !m_vTriggers[iG].compare( m_mbTriggerName ) )
	{ continue; }
      if( !m_vTriggers[iG].compare( m_allName ) )
	{ continue; }

      TGraphAsymmErrors* g = vEffGrf[iG][iX];
      styleTool->SetHStyle( g, style++ );
      g->Draw("p");
      leg.AddEntry( g, m_vTriggers[iG].c_str() );
    } // end loop over iG
    
    leg.Draw("same");

    TLine line( xMin, 1, xMax, 1);
    line.Draw();    
    
    drawTool->DrawAtlasInternalDataRight( 0, 0, m_is_pPb ); 
    drawTool->DrawLeftLatex( 0.4, 0.17, cLabel );

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
  std::string configName = "config/configJetsFwd";
  
  TEnv config_pPb;
  TEnv config_pp;

  std::string triggerMenu;
  
  config_pPb.ReadFile( Form("%s_pPb_data.cfg", configName.c_str() ), EEnvLevel(0));
  config_pp.ReadFile ( Form("%s_pp_data.cfg" , configName.c_str() ), EEnvLevel(0));

  triggerMenu = config_pPb.GetValue("triggerMenu", "");
  std::string trigger_pPb =
    config_pPb.GetValue( Form("finalTrigger.%s",triggerMenu.c_str()), "");

  triggerMenu = config_pp.GetValue("triggerMenu", "" );
  std::string trigger_pp =
    config_pp.GetValue( Form("finalTrigger.%s",triggerMenu.c_str()), "");

  std::cout << " -- " << trigger_pPb << std::endl;
  std::cout << " -- " << trigger_pp << std::endl;
  
  TFile* fIn_pPb = TFile::Open("output/output_pPb_data/c_myOut_pPb_data.root");
  TFile* fIn_pp  = TFile::Open("output/output_pp_data/c_myOut_pp_data.root");
  TFile* fOut    = new TFile  ("output/all/c_myOut_data.root","recreate");
  
  for( uint ystar1Bin = 0; ystar1Bin < m_nVarYstarBinsA; ystar1Bin++ ){
    for( uint ystar2Bin = 0; ystar2Bin < m_nVarYstarBinsB; ystar2Bin++ ){
      for( uint pt1Bin = 0; pt1Bin < m_nVarPtBins; pt1Bin++ ){
	for( uint pt2Bin = 0; pt2Bin < m_nVarPtBins; pt2Bin++ ){
	  
	  double ystar1Low = m_varYstarBinningA[ ystar1Bin ];
	  double ystar1Up  = m_varYstarBinningA[ ystar1Bin + 1 ];
	  double ystar1Center = ystar1Low + 0.5 * ( ystar1Up - ystar1Low );
	  
	  double ystar2Low = m_varYstarBinningB[ ystar2Bin ];
	  double ystar2Up  = m_varYstarBinningB[ ystar2Bin + 1 ];
	  double ystar2Center = ystar2Low + 0.5 * ( ystar2Up - ystar2Low );

	  double pt1Low  = m_varPtBinning[ pt1Bin ];
	  double pt1Up   = m_varPtBinning[ pt1Bin + 1 ];

	  double pt2Low  = m_varPtBinning[ pt2Bin ];


	  if( !IsForwardYstar( ystar1Center ) &&
	      !IsForwardYstar( ystar2Center ) )
	    { continue; }
	  if( pt1Low < pt2Low  )
	    { continue; }
	    
	  std::string hTag =
	    Form("%s_%s_%s_%s",
		 anaTool->GetName( ystar1Low, ystar1Up, "Ystar1").c_str(),
		 anaTool->GetName( ystar2Low, ystar2Up, "Ystar2").c_str(),
		 anaTool->GetName( pt1Low , pt1Up , "Pt1" ).c_str(),
		 anaTool->GetName( pt2Low , pt2Low, "Pt2" ).c_str() );

	  std::string hName_pPb =
	    Form("h_dPhi_%s_%s", hTag.c_str(), trigger_pPb.c_str() );
	  std::string hName_pp =
	    Form("h_dPhi_%s_%s", hTag.c_str(), trigger_pp.c_str() );
	    
	  TCanvas* c_pPb =
	    static_cast<TCanvas*>
	    ( fIn_pPb->Get( Form("c_%s_pPb_data", hName_pPb.c_str() ) ) );
	  TCanvas* c_pp  =
	    static_cast<TCanvas*>
	    ( fIn_pp->Get( Form("c_%s_pp_data", hName_pp.c_str() ) ) );
	  
	  TH1* h_pPb =
	    static_cast<TH1D*>( c_pPb->GetPrimitive( hName_pPb.c_str() ) );
	  styleTool->SetHStyle( h_pPb, 0 );
	  TH1* h_pp  =
	    static_cast<TH1D*>( c_pp->GetPrimitive( hName_pp.c_str() ) );
	  styleTool->SetHStyle( h_pp, 1 );
	  
	  TF1* f_pPb =
	    static_cast<TF1*>( c_pPb->GetPrimitive
			       ( Form("f_%s", hName_pp.c_str())));
	  styleTool->SetHStyle( f_pPb, 0 );
	  f_pPb->SetLineColor( h_pPb->GetLineColor() );
	  TF1* f_pp  =
	    static_cast<TF1*>( c_pp->GetPrimitive
			       ( Form("f_%s", hName_pp.c_str())));
	  styleTool->SetHStyle( f_pp, 1 );
	  f_pp->SetLineColor( h_pp->GetLineColor() );
	  
	  TCanvas c("c","c", 800, 600 );
	    
	  TLegend leg( 0.27, 0.41, 0.38, 0.52 );
	  styleTool->SetLegendStyle( &leg );
	  leg.AddEntry( h_pPb, "p+Pb");
	  leg.AddEntry( h_pp , "pp");

	  h_pPb->Draw("epsame");
	  h_pp->Draw("epsame");

	  f_pPb->Draw("same");
	  f_pp->Draw("same");

	  leg.Draw("same");

	  if( h_pPb->GetMaximum() > h_pp->GetMaximum() ){
	    h_pPb->SetMaximum( h_pPb->GetMaximum() * 1.1 );
	    h_pp->SetMaximum ( h_pPb->GetMaximum() * 1.1 );
	  } else {
	    h_pPb->SetMaximum( h_pp->GetMaximum() * 1.1 );
	    h_pp->SetMaximum ( h_pp->GetMaximum() * 1.1 );
	  }
	    
	  drawTool->DrawLeftLatex
	    ( 0.13, 0.87,anaTool->GetLabel( ystar1Low, ystar1Up, "#it{y}*_{1}" ) );
	  drawTool->DrawLeftLatex
	    ( 0.13, 0.82,anaTool->GetLabel( ystar2Low, ystar2Up, "#it{y}*_{1}" ) );
	  drawTool->DrawLeftLatex
	    ( 0.13, 0.76,anaTool->GetLabel( pt1Low, pt1Up , "#it{p}_{T}^{1}" ) );
	  drawTool->DrawLeftLatex
	    ( 0.13, 0.69,anaTool->GetLabel( pt2Low, pt2Low, "#it{p}_{T}^{2}" ) );

	  drawTool->DrawAtlasInternal();

	  c.SaveAs( Form("output/all/data/h_dPhi_%s.png", hTag.c_str() ));
	  SaveAsROOT( c, Form("h_dPhi_%s", hTag.c_str() ) );

	  delete h_pPb;
	  delete h_pp;
	  delete f_pPb;
	  delete f_pp;
	  delete c_pPb;
	  delete c_pp;
	}
      }
    } 
  }
  
  std::cout << "DONE! Closing " << fOut->GetName() << std::endl;
  fOut->Close();
  std::cout << "......Closed  " << fOut->GetName() << std::endl;
}
