#include <TLine.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TEnv.h>
#include <TCanvas.h>
#include <TEfficiency.h>
#include <TGraphAsymmErrors.h>

#include <iostream>
#include <cmath>

#include <boost/assign.hpp>
#include <boost/format.hpp>

#include "MyRoot.h"

#include "DiJetAnalysisData.h"

DiJetAnalysisData::DiJetAnalysisData() : DiJetAnalysisData( true, true )
{}

DiJetAnalysisData::DiJetAnalysisData( bool isData,
				      bool is_pPb )
  : DiJetAnalysis( isData, is_pPb)
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

  PlotEtaPhiPtMap( m_vTriggerEtaPhiMap );
  PlotEtaPhiPtMap( m_vTriggerEtaPtMap  );

  PlotSpectra( m_vTriggerEtaSpect, "spect" );

  PlotEfficiencies( m_vTriggerEtaSpectSim, m_vTriggerEtaSpect, "eff" );
  
  PlotDeltaPhi( m_vDphi, "dPhi" );
  
  std::cout << "DONE! Closing " << cfNameOut << std::endl;
  m_fOut->Close();
  std::cout << "......Closed  " << cfNameOut << std::endl;
}

//---------------------------------
//            Fill Tree
//---------------------------------

void DiJetAnalysisData::LoadTriggers(){
  std::string configName   = "config/configJetsFwd";
  std::string configSuffix = m_is_pPb ? "_pPb.cfg" : "_pp.cfg";
  configName += configSuffix;
  
  TEnv config;
  config.ReadFile( configName.c_str(), EEnvLevel(0));

  std::string triggerMenu = config.GetValue("triggerMenu", " ");
  m_vTriggers =
    anaTool->vectorise
    ( config.GetValue( Form("triggers.%s",triggerMenu.c_str()), ""), " ");

  std::vector< std::string> vRefTriggers =
    anaTool->vectorise
    ( config.GetValue( Form("refTriggers.%s",triggerMenu.c_str()), "" ), " ");
  m_vRefTriggerIndex.reserve( vRefTriggers.size() );
    
  // Get Index for reference trigger
  for( uint iR = 0 ; iR < vRefTriggers.size(); iR++ ){
      for( uint iT = 0 ; iT < vRefTriggers.size(); iT++ ){
	if( !m_vTriggers[iT].compare( vRefTriggers[iR] ) )
	  { m_vRefTriggerIndex.push_back(iT); }
      }
  }
  
  for( uint iT = 0 ; iT < m_vRefTriggerIndex.size(); iT++ ){
    int refTrigIndex = m_vRefTriggerIndex[iT];
    std::cout << " For: "         << m_vTriggers[iT]
	      << " Ref Trig is: " << m_vTriggers[refTrigIndex]
	      << std::endl;      
  }
  
  m_vTholdPtTriggers =
    anaTool->vectoriseD
    ( config.GetValue( Form("tholdTriggers.%s",triggerMenu.c_str()), "" ), " ");

  m_vEffPtTriggers =
    anaTool->vectoriseD
    ( config.GetValue( Form("effTriggers.%s",triggerMenu.c_str()), "" ), " ");
  
  m_nTriggers = m_vTriggers.size();

  for( uint iT = 0; iT < m_nTriggers; iT++ ){
    std::string trigger = m_vTriggers[iT];
    if( trigger.find("_mb_") != std::string::npos ){ 
      m_mbTriggerI = iT;
      std::cout << "Found " << trigger << " at " << m_mbTriggerI << std::endl;
      break;
    }
  }
}

void DiJetAnalysisData::SetupHistograms(){

  for( auto& trigger : m_vTriggers ){
    std::cout << "Making - " << trigger << " histograms " << std::endl;          
    
    // -------- maps ---------
    m_vTriggerEtaPhiMap.push_back
      ( new TH2D( Form("h_etaPhiMap_%s", trigger.c_str() ), 
		   ";#eta;#phi",
		   m_nEtaMapBins, m_etaMapMin, m_etaMapMax,
		   m_nPhiMapBins, m_phiMapMin, m_phiMapMax ) ) ;
    AddHistogram( m_vTriggerEtaPhiMap.back() );
     
    m_vTriggerEtaPtMap.push_back
      ( new TH2D( Form("h_etaPtMap_%s", trigger.c_str() ), 
		   ";#eta;#it{p}_{T}",
		   m_nEtaMapBins, m_etaMapMin, m_etaMapMax,
		   m_nPtMapBins , m_ptMapMin , m_ptMapMax) ) ;
    AddHistogram( m_vTriggerEtaPtMap.back() );

    // -------- spect --------
    m_vTriggerEtaSpect.push_back
      ( new TH2D( Form("h_etaSpect_%s", trigger.c_str() ), 
		  ";#eta;#it{p}_{T} [GeV];",
		  m_nVarFwdEtaBins, 0, 1,
		  m_nPtSpectBins,
		  m_ptSpectMin, m_ptSpectMax ) ) ;
    m_vTriggerEtaSpect.back()->GetXaxis()->
      Set( m_nVarFwdEtaBins, &( m_varFwdEtaBinning[0] ) );
    AddHistogram( m_vTriggerEtaSpect.back() );
      
    m_vTriggerEtaSpectSim.push_back
      ( new TH2D ( Form("h_etaSpectSim_%s", trigger.c_str() ), 
		   ";#eta;#it{p}_{T} [GeV]",
		   m_nVarFwdEtaBins, 0, 1,
		   m_nPtSpectBins,
		   m_ptSpectMin, m_ptSpectMax ) ) ;
    m_vTriggerEtaSpectSim.back()->GetXaxis()->
      Set( m_nVarFwdEtaBins, &( m_varFwdEtaBinning[0] ) );
    AddHistogram( m_vTriggerEtaSpectSim.back() );

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
    
    THnSparse* hn =
      new THnSparseD( Form("hn_dPhi_%s", trigger.c_str() ), "",
		      nDim, &nDphiBins[0], &dPhiMin[0], &dPhiMax[0] );
    hn->GetAxis(0)->Set( m_nVarEtaBins, &( m_varEtaBinning[0] ) );
    hn->GetAxis(1)->Set( m_nVarEtaBins, &( m_varEtaBinning[0] ) );
    hn->GetAxis(2)->Set( m_nVarPtBins , &( m_varPtBinning[0]  ) );
    hn->GetAxis(3)->Set( m_nVarPtBins , &( m_varPtBinning[0]  ) );

    m_vDphi.push_back( hn );
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

  for( uint iT = 0; iT < m_nTriggers; iT++ ){
    m_tree->SetBranchAddress
      ( Form("passed_%s", m_vTriggers[iT].c_str() ), &mTriggerFired[iT] );
    m_tree->SetBranchAddress
      ( Form("prescale_%s", m_vTriggers[iT].c_str() ),&vTriggerPrescale[iT] );
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

    bool haveHighPtJet = false; 
    bool mbTrigFired   = false;
    int  nHltTrigNotFired  = 0;
    // SPECTRA AND ETA PHI PT MAPS
    // loop over passed triggers
    for( uint iT = 0; iT < m_nTriggers; iT++ ){

      // loop over jets 
      for( auto& jet : vR_jets ){
	// ETA-PHI
	double jetPt = jet.Pt()/1000.;
	if( jetPt > 35 ) haveHighPtJet = true;
      }

      // check if we have that trigger
      // if we dont - continue
      if( !mTriggerFired[iT] ) {
	mTrigFailed[iT]++;
	if( iT > 0 ){ nHltTrigNotFired++; }
	continue;	
      }
      if( iT == 0 ){ mbTrigFired = true; }

      if( vR_jets.size() >= 2 )
	{ AnalyzeDeltaPhi( m_vDphi[iT], vR_jets ); }
      
      // loop over jets 
      for( auto& jet : vR_jets ){
	// cleaned jets have px, py, pz set to 0
	if( jet.Pt() == 0 ) continue;

	// ETA-PHI
	double jetEta = jet.Eta();
	double jetPhi = jet.Phi();
	double jetPt = jet.Pt()/1000.;

	if( jetPt > 35 ) haveHighPtJet = true;
	
	m_vTriggerEtaPhiMap[iT]->Fill( jetEta, jetPhi );
	m_vTriggerEtaPtMap [iT]->Fill( jetEta, jetPt ); 
	// eta cut

	// convert positive eta to negative because
	// in pp it doesnt matter.
	// our histos run in negative eta
	// (due to pPb configuration)
	// the labels will be taken care of so it is ok
	double jetEtaAdj = AdjustEtaForPP( jetEta );
	
	if( jetPt >= m_vTholdPtTriggers[iT] )
	  { mTrigJetsTotal[iT] ++; }
	
	// saves computation time
	// can theoreticall leave out
	// will go to overflow bins in histos
	if( jetEtaAdj > -constants::FETAMIN ) { continue; }
	m_vTriggerEtaSpect[iT]->Fill( jetEtaAdj, jetPt );
	
	// MB eff is set at zero
	// will always pass taht trigger
	if( jetPt >= m_vEffPtTriggers[iT] )
	  { mTrigJetsForward[iT]++; }
      } // end loop over jets

      mTrigPassed[iT]++;
    } // end loop over iT

    // EFFICIENCIES
    AnalyzeEff( vTrig_jets, vR_jets, mTriggerFired );
    
    if( haveHighPtJet && mbTrigFired && nHltTrigNotFired == 3 ){
      std::cout << "\nReco Jets:" << std::endl;
      for(auto & jet : vR_jets ){
	std::cout << jet.Pt() << "  " << jet.Eta() << std::endl;
      }
      std::cout << "\nTrig Jets:" << std::endl;
      for(auto & jet : vTrig_jets ){
	std::cout << jet.Pt()<< "  " << jet.Eta() << std::endl;
      }
      std::cout << "\nTrig:" << std::endl;
      for(auto & trig : mTriggerFired ){
	std::cout << trig.second << " - " << m_vTriggers[trig.first]
		  << "  " << vTriggerPrescale[trig.first] << std::endl;
      }
      std::cout << "----------" << std::endl;
    }
  } // end event loop
  std::cout << "DONE! Has: " << nEventsTotal << " events." << std::endl;

  for( unsigned iT = 0; iT < m_vTriggers.size(); iT++ ){
    std::cout <<  m_vTriggers[iT]
	      << "  has: " << mTrigJetsForward[iT]
	      << "  FWD jets, " << mTrigJetsTotal[iT]
	      << "  TOTAL jets, with " << mTrigPassed[iT]
	      << "  evPassed and "<< mTrigFailed[iT]
	      << "  failed. Total = "
	      << mTrigPassed[iT] + mTrigFailed[iT]
	      << std::endl;
  }

  m_fIn->Close();
}

//---------------------------
//       Analysis
//---------------------------
void DiJetAnalysisData::AnalyzeEff( std::vector< TLorentzVector >& vTrig_jets,
				    std::vector< TLorentzVector >& vR_jets,
				    std::map< int, bool >& mTriggerFired ){
  if( vTrig_jets.empty() ){ return; }

  // sort trigger jets
  std::sort( vTrig_jets.begin(), vTrig_jets.end(),
	     anaTool->sortByDecendingPt );

  // take highest pt trigger jet in forward eta range
  TLorentzVector* tJet = NULL;
  for( auto& jet : vTrig_jets ){
    if( AdjustEtaForPP( jet.Eta() ) > -constants::FETAMIN )
      { continue; }
    tJet = &jet;
    break; // found highest forward trigger jet
  }    

  // didnt find a trigger jet in forward eta range
  if( !tJet ) return;

  // now fill histos for triggers that could have passed
  for( uint iT = 0; iT < m_nTriggers; iT++ ){
    // check if this triggers reference trigger fired
    int refTrigIndex = m_vRefTriggerIndex[iT];
    if( !mTriggerFired[refTrigIndex] ){ continue; }
    
    double tJetPt = tJet->Pt()/1000;
   
    // now check if we should fill for that trigger
    // only if trigger jet is above threshold
    // for that trigger
    if( tJetPt < m_vTholdPtTriggers[iT] ){ continue; } 

    // fill jets for that trigger
    for( auto& jet : vR_jets ){
      // avoid cleaned jets
      if( jet.Pt() == 0 ) continue;

      double jetEtaAdj = AdjustEtaForPP( jet.Eta() );

      // eta cut. want forward trigger jet
      if( jetEtaAdj > -constants::FETAMIN ){ continue; }

      m_vTriggerEtaSpectSim[iT]->
	Fill( jetEtaAdj, jet.Pt()/1000. );
    } // end loop over jets
  } // end loop over iT
}

//---------------------------------
//            Plot Data
//---------------------------------

void DiJetAnalysisData::LoadHistograms(){
  m_fIn = TFile::Open( m_rootFname.c_str() ); 

  for( auto& trigger : m_vTriggers ){
    // -------- maps ---------
    m_vTriggerEtaPhiMap.
      push_back( static_cast< TH2D* >
		 ( m_fIn-> Get( Form("h_etaPhiMap_%s",
				     trigger.c_str() ))));
    m_vTriggerEtaPhiMap.back()->SetDirectory(0);

    m_vTriggerEtaPtMap.
      push_back( static_cast< TH2D* >
		 ( m_fIn->Get( Form("h_etaPtMap_%s", trigger.c_str() ))));
    m_vTriggerEtaPtMap.back()->SetDirectory(0);

    // -------- spect --------
    m_vTriggerEtaSpect.
      push_back( static_cast< TH2D *>
		 ( m_fIn->Get( Form("h_etaSpect_%s", trigger.c_str() ))));
    m_vTriggerEtaSpect.back()->SetDirectory(0);

    m_vTriggerEtaSpectSim.
      push_back( static_cast< TH2D *>
		 ( m_fIn->Get( Form("h_etaSpectSim_%s", trigger.c_str() ))));
    m_vTriggerEtaSpectSim.back()->SetDirectory(0);

    // -------- dPhi- --------
    m_vDphi.
      push_back( static_cast< THnSparse *>
		 ( m_fIn->Get( Form("hn_dPhi_%s", trigger.c_str() ))));
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
  TH2*   hMB = vTrigSpect[ m_mbTriggerI ];
  int nXbins = hMB->GetNbinsX();

  std::vector< std::vector< TH1* > > vSpect;
  vSpect.resize( m_nTriggers );
  
  double max = -1;
  for( uint iT = 0; iT < m_nTriggers; iT++){
    std::string trigger = m_vTriggers[iT];

    for( int xBin = 1; xBin <= nXbins; xBin++ ){
      double etaMin = hMB->GetXaxis()->GetBinLowEdge( xBin );
      double etaMax = hMB->GetXaxis()->GetBinUpEdge ( xBin );

      std::string etaName = anaTool->GetName(etaMin,etaMax,"Eta");
      
      TH1* h_etaSpect =
	vTrigSpect[iT]->
	ProjectionY( Form("h_%s_%s_%s",
			  type.c_str(),
			  trigger.c_str(),
			  etaName.c_str() ),
		     xBin, xBin );
      h_etaSpect->SetTitle( GetEtaLabel( etaMin, etaMax ).c_str() );
      h_etaSpect->GetYaxis()->SetTitle( yAxisTitle.c_str() );
      h_etaSpect->Scale( 1./ptSpectWidth );
      vSpect[iT].push_back( h_etaSpect );
      
      if( max < h_etaSpect->GetMaximum() )
	{ max = h_etaSpect->GetMaximum(); }
    } // end loop over xBin
  } // end loop over iT

  // set maxima globally for all spectra hists.
  // easier to compare. Set on log scale.
  double power = log10(max);
  power = std::ceil(power);
  max = pow( 10, power );

  //------------------------------------------------
  //------- Draw Eta as Fucntion of Triggers -------
  //------------------------------------------------
  for( uint iT = 0; iT < m_nTriggers; iT++ ){
    std::string cName  = m_vTriggers[iT];
    std::string cLabel = m_vTriggers[iT];

    TCanvas c( "c", cLabel.c_str(), 800, 600 );
    c.SetLogy();
    
    TLegend leg( 0.7, lY0, lX1, lY1 );
    styleTool->SetLegendStyle( &leg  );

    int style = 0;
    for( int iE = 0; iE < nXbins; iE++ ){
      TH1* h = vSpect[iT][iE];
      styleTool->SetHStyle( h, style++ );
      h->Draw("epsame");
      h->SetMinimum( 1 );
      h->SetMaximum( max );
      leg.AddEntry( h, h->GetTitle() );
      h->SetTitle("");
    } // end loop over iE
    
    leg.Draw("same");
    drawTool->DrawAtlasInternalDataRight
      ( 0, 0, m_is_pPb ); 
    drawTool->DrawLeftLatex
      ( 0.15, 0.17, cLabel.c_str(), 1 );
    SaveAsAll( c, type, cName );
  } // end loop over iT
  
  //------------------------------------------------
  //------- Draw Triggers as Fucntion of Eta -------
  //------------------------------------------------
  for( int iE = 0; iE < nXbins; iE++ ){
    // no longer have title.
    // Need to make output files
    // with eta names
    int      xBin = iE + 1;
    double etaMin = hMB->GetXaxis()->GetBinLowEdge( xBin );
    double etaMax = hMB->GetXaxis()->GetBinUpEdge ( xBin );
    
    std::string cName  = anaTool->GetName( etaMin, etaMax, "Eta" );
    std::string cLabel = GetEtaLabel( etaMin, etaMax );
    
    TCanvas c( "c", cLabel.c_str(), 800, 600 );
    c.SetLogy();
    
    TLegend leg( lX0, lY0, lX1, lY1 );
    styleTool->SetLegendStyle( &leg  );
    
    int style = 0;
    for( uint iT = 0; iT < m_nTriggers; iT++ ){
      TH1* h = vSpect[iT][iE];
      styleTool->SetHStyle( h, style++ );
      h->Draw("epsame");
      h->SetMinimum( 1 );
      h->SetMaximum( max );
      leg.AddEntry( h, m_vTriggers[iT].c_str() );
    } // end loop over iT
    
    leg.Draw("same");
    drawTool->DrawAtlasInternalDataRight
      ( 0, 0, m_is_pPb ); 
    drawTool->DrawLeftLatex
      ( 0.15, 0.17, cLabel.c_str(), 1 );

    SaveAsAll( c, type, cName );
  } // end loop over iE

  // delete
  for( uint iT = 0; iT < m_nTriggers; iT++ ){
    for( int iE = 0; iE < nXbins; iE++ )
      { delete vSpect[iT][iE]; }
  }
}

void DiJetAnalysisData::PlotEfficiencies( std::vector< TH2* >& vTrigSpect,
					  std::vector< TH2* >& vTrigSpectRef,
					  const std::string& type ){

  // if there is no mb trigger, return
  if( !vTrigSpect[ m_mbTriggerI ] ){ return; }  

  double lX0 = 0.13;
  double lY0 = 0.75;
  double lX1 = 0.39;
  double lY1 = 0.87;

  double xMin = vTrigSpect[ m_mbTriggerI ]->GetYaxis()->GetXmin() - 10;
  double xMax = 60.; //vTrigSpect[ m_mbTriggerI ]->GetYaxis()->GetXmax();

  // use this as reference because
  // it should be in every file
  TH2*   hMB = vTrigSpect[ m_mbTriggerI ];
  int nXbins = hMB->GetNbinsX();
  
  std::vector< std::vector< TH1* > > vSpect;
  std::vector< std::vector< TH1* > > vSpectRef;
  std::vector< std::vector< TGraphAsymmErrors* > > vEffGrf;
  
  vSpect.resize   ( m_nTriggers );
  vSpectRef.resize( m_nTriggers );
  vEffGrf.resize  ( m_nTriggers );
  
  for( uint iT = 0; iT < m_nTriggers; iT++ ){
    std::string trigger = m_vTriggers[iT];
 
    for( int xBin = 1; xBin <= nXbins; xBin++ ){
      // should all be the same
      double etaMin = hMB->GetXaxis()->GetBinLowEdge( xBin );
      double etaMax = hMB->GetXaxis()->GetBinUpEdge ( xBin );

      std::string etaLabel = GetEtaLabel( etaMin, etaMax );
      std::string etaName  = anaTool->GetName( etaMin, etaMax, "Eta" );
      
      TH1* h_etaSpect = vTrigSpect[iT]->
	ProjectionY( Form("h_%s_%s_%s",
			  type.c_str(),
			  trigger.c_str(),
			  etaName.c_str() ),
		     xBin, xBin );
      h_etaSpect->SetTitle( etaLabel.c_str() );
      vSpect[iT].push_back( h_etaSpect );

      TH1* h_etaSpectRef = vTrigSpectRef[iT]->
	ProjectionY( Form("h_%s_%s_ref_%s",
			  type.c_str(),
			  trigger.c_str(),
			  etaName.c_str() ),
		     xBin, xBin );
      h_etaSpectRef->SetTitle( etaLabel.c_str() );
      vSpectRef[iT].push_back( h_etaSpectRef );

      TGraphAsymmErrors* g_etaEff = new TGraphAsymmErrors();
      g_etaEff->SetName( Form("gr_%s_%s_%s",
			      type.c_str(),
			      trigger.c_str(),
			      etaName.c_str() ) );
      g_etaEff->SetTitle( etaLabel.c_str() );
      vEffGrf[iT].push_back( g_etaEff );      
    } // end loop over eta
  } // end loop over triggers

  //------------------------------------------------
  //--------- Now Divide by Appropriate Ref  -------
  //------------------------------------------------
  for( uint iT = 0; iT < m_nTriggers; iT++ ){
    for( int iE = 0; iE < nXbins; iE++ ){
      int  refTrigIndex = m_vRefTriggerIndex[iT];
      TH1* hSpect    = vSpect[iT][iE];
      TH1* hSpectRef = vSpectRef[refTrigIndex][iE];
      vEffGrf[iT][iE]->Divide( hSpect, hSpectRef );
    }
  }
  
  //------------------------------------------------
  // ----------------- Now Plot --------------------
  //------------------------------------------------
  std::string gLabel = ";#it{p}_{T} [GeV];#it{#varepsilon}_{Trigger}";

  //------------------------------------------------
  //------- Draw Eta as Fucntion of Triggers -------
  //------------------------------------------------
  for( uint iT = 0; iT < m_nTriggers; iT++ ){
    // dont draw MB trigger
    if( !m_vTriggers[iT].compare( m_vTriggers[m_mbTriggerI] ) )
      { continue; }
    std::string cName  = m_vTriggers[iT];
    std::string cLabel = m_vTriggers[iT];

    TCanvas c( "c", cLabel.c_str(), 800, 600 );
    styleTool->SetCStyleEff( c, xMin, m_effMin, xMax, m_effMax, gLabel );
    
    TLegend leg( lX0, lY0, lX1, lY1 );
    styleTool->SetLegendStyle( &leg  );

    int style = 0;
    for( int iE = 0; iE < nXbins; iE++ ){
      TGraphAsymmErrors* g = vEffGrf[iT][iE];
      styleTool->SetHStyle( g, style++ );
      g->Draw("p");
      leg.AddEntry( g, g->GetTitle() );
      g->SetTitle("");
    } // end loop over iE
    
    leg.Draw("same");

    TLine line( xMin, 1, xMax, 1);
    line.Draw();
    
    drawTool->DrawAtlasInternalDataRight
      ( 0, 0, m_is_pPb ); 
    drawTool->DrawLeftLatex
      ( 0.15, 0.17, cLabel.c_str(), 1 );
    
    SaveAsAll( c, type, cName );
  } // end loop over iT
  
  //------------------------------------------------
  //------- Draw Triggers as Fucntion of Eta -------
  //------------------------------------------------
  for( int iE = 0; iE < nXbins; iE++ ){
    // no longer have title.
    // Need to make output files
    // with eta names
    int      xBin = iE + 1;
    double etaMin = hMB->GetXaxis()->GetBinLowEdge( xBin );
    double etaMax = hMB->GetXaxis()->GetBinUpEdge ( xBin );
    
    std::string cName  = anaTool->GetName( etaMin, etaMax, "Eta" );
    std::string cLabel = GetEtaLabel( etaMin, etaMax );
    
    TCanvas c( "c", cLabel.c_str(), 800, 600 );
    styleTool->SetCStyleEff( c, xMin, m_effMin, xMax, m_effMax, gLabel );
    
    TLegend leg( lX0, lY0, lX1, lY1 );
    styleTool->SetLegendStyle( &leg  );
    
    int style = 0;
    for( uint iT = 0; iT < m_nTriggers; iT++ ){
      // dont draw MB trigger
      if( !m_vTriggers[iT].compare( m_vTriggers[m_mbTriggerI] ) )
	{ continue; }
  
      TGraphAsymmErrors* g = vEffGrf[iT][iE];
      styleTool->SetHStyle( g, style++ );
      g->Draw("p");
      leg.AddEntry( g, m_vTriggers[iT].c_str() );
    } // end loop over iT
    
    leg.Draw("same");

    TLine line( xMin, 1, xMax, 1);
    line.Draw();    
    
    drawTool->DrawAtlasInternalDataRight
      ( 0, 0, m_is_pPb ); 
    drawTool->DrawLeftLatex
      ( 0.15, 0.17, cLabel.c_str(), 1 );

    SaveAsAll( c, type, cName );
  } // end loop over iE

  /*
  // delete
  for( uint iT = 0; iT < m_nTriggers; iT++ ){
    for( int iE = 0; iE < nXbins; iE++ )
      { delete vSpect[iT][iE]; }
  }
  */
  
}

void DiJetAnalysisData::PlotDeltaPhi( std::vector< THnSparse* >& vhn,
				      const std::string& type ){
  std::vector< TH1* > vDphi;
    
  for( uint iT = 0; iT < m_nTriggers; iT++ ){
    THnSparse* hn = vhn[iT];
    std::string trigger = m_vTriggers[iT];
    
    TAxis* eta1Axis = hn->GetAxis(0); TAxis* pt1Axis = hn->GetAxis(2);
    TAxis* eta2Axis = hn->GetAxis(1); TAxis* pt2Axis = hn->GetAxis(3); 

    int nEta1Bins = eta1Axis->GetNbins();
    int nEta2Bins = eta2Axis->GetNbins();
    int nPt1Bins  =  pt1Axis->GetNbins();
    int nPt2Bins  =  pt2Axis->GetNbins();
    
    for( int eta1Bin = 1; eta1Bin <= nEta1Bins; eta1Bin++ ){
      eta1Axis->SetRange( eta1Bin, eta1Bin );
      for( int eta2Bin = 1; eta2Bin <= nEta2Bins; eta2Bin++ ){
	eta2Axis->SetRange( eta2Bin, eta2Bin );
	for( int pt1Bin = 1; pt1Bin <= nPt1Bins; pt1Bin++ ){
	  pt1Axis->SetRange( pt1Bin, pt1Bin );
	  for( int pt2Bin = 1; pt2Bin <= nPt2Bins; pt2Bin++ ){
	    pt2Axis->SetRange( pt2Bin, nPt2Bins );
	    
	    // Take projection onto the dPhi axis
	    TH1* hDphi = hn->Projection( 4 );
	    styleTool->SetHStyle( hDphi, 0 );
	    vDphi.push_back( hDphi );

	    double eta1Low, eta1Up, eta2Low, eta2Up;
	    double pt1Low , pt1Up , pt2Low , pt2Up;

	    anaTool->GetBinRange
	      ( eta1Axis, eta1Bin, eta1Bin, eta1Low, eta1Up );
	    anaTool->GetBinRange
	      ( eta2Axis, eta2Bin, eta2Bin, eta2Low, eta2Up );
	    anaTool->GetBinRange
	      ( pt1Axis, pt1Bin, pt1Bin, pt1Low, pt1Up );
	    anaTool->GetBinRange
	      ( pt2Axis, pt2Bin, pt2Bin, pt2Low, pt2Up );

	    hDphi->SetName
	      ( Form( "h_dPhi_%s_%s_%s_%s_%s",
		      trigger.c_str(),
		      anaTool->GetName( eta1Low, eta1Up, "Eta1").c_str(),
		      anaTool->GetName( eta2Low, eta2Up, "Eta2").c_str(),
		      anaTool->GetName( pt1Low , pt1Up , "Pt1" ).c_str(),
		      anaTool->GetName( pt2Low , pt2Low, "Pt2" ).c_str() ) );

	    TCanvas c( "c", hDphi->GetName(), 800, 600 );
	    hDphi->Draw();
	    if( hDphi->GetEntries() )
	      { hDphi->Scale( 1./hDphi->Integral() ); }
	    hDphi->GetYaxis()->SetTitle("Normalized Count");
	    hDphi->SetTitle("");
	    
	    drawTool->DrawLeftLatex
	      ( 0.13, 0.87,
	        GetLabel( eta1Low, eta1Up, "#eta_{1}" ).c_str(), 1 );
	    drawTool->DrawLeftLatex
	      ( 0.13, 0.82,
		GetLabel( eta2Low, eta2Up, "#eta_{2}" ).c_str(), 1 );
	    drawTool->DrawLeftLatex
	      ( 0.13, 0.76,
		GetLabel( pt1Low, pt1Up, "#it{p}_{T}^{1}" ).c_str(), 1 );
	    drawTool->DrawLeftLatex
	      ( 0.13, 0.69,
		GetLabel( pt2Low, pt2Low, "#it{p}_{T}^{2}" ).c_str(), 1 );
	    
	    drawTool->DrawAtlasInternalDataRight
	      ( 0, 0, m_is_pPb ); 

	    SaveAsROOT( c, hDphi->GetName() );
	  }
      	}
      }      
    }
  } // end loop over iT
}

void DiJetAnalysisData::PlotEtaPhiPtMap( std::vector< TH2* >& vTrigHin ){
  TCanvas c_map("c_map","c_map",800,600);
  
  for( uint iT = 0; iT < m_nTriggers; iT++){
    TH2* h = vTrigHin[iT];
    h->Draw("col");
    styleTool->SetHStyle( h, 0 );
    drawTool->DrawAtlasInternalDataRight
      ( 0, -0.55, m_is_pPb );  
    
    SaveAsAll( c_map, h->GetName() );
  }
}

void DiJetAnalysisData::PlotDataTogether(){
  std::string configName = "config/configJetsFwd";
  
  TEnv config_pPb;
  TEnv config_pp;

  std::string triggerMenu;
  
  config_pPb.ReadFile( Form("%s_pPb.cfg", configName.c_str() ), EEnvLevel(0));
  config_pp.ReadFile ( Form("%s_pp.cfg" , configName.c_str() ), EEnvLevel(0));

  triggerMenu = config_pPb.GetValue("triggerMenu", " ");
  std::string trigger_pPb =
    config_pPb.GetValue( Form("finalTrigger.%s",triggerMenu.c_str()), "");

  triggerMenu = config_pp.GetValue("triggerMenu", " " );
  std::string trigger_pp =
    config_pp.GetValue( Form("finalTrigger.%s",triggerMenu.c_str()), "");

  std::cout << " -- " << trigger_pPb << std::endl;
  std::cout << " -- " << trigger_pp << std::endl;
  
  TFile* fIn_pPb = TFile::Open("output/output_data_pPb/c_myOut_data_pPb.root");
  TFile* fIn_pp  = TFile::Open("output/output_data_pp/c_myOut_data_pp.root");
  TFile* fOut    = new TFile  ("output/c_myOut_data.root","recreate");
  
  for( uint eta1Bin = 0; eta1Bin < m_nVarEtaBins; eta1Bin++ ){
    for( uint eta2Bin = 0; eta2Bin < m_nVarEtaBins; eta2Bin++ ){
      for( uint pt1Bin = 0; pt1Bin < m_nVarPtBins; pt1Bin++ ){
	for( uint pt2Bin = 0; pt2Bin < m_nVarPtBins; pt2Bin++ ){

	  double eta1Low = m_varEtaBinning[ eta1Bin ];
	  double eta1Up  = m_varEtaBinning[ eta1Bin + 1 ];
	    
	  double eta2Low = m_varEtaBinning[ eta2Bin ];
	  double eta2Up  = m_varEtaBinning[ eta2Bin + 1 ];

	  double pt1Low  = m_varPtBinning[ pt1Bin ];
	  double pt1Up   = m_varPtBinning[ pt1Bin + 1 ];

	  double pt2Low  = m_varPtBinning[ pt2Bin ];
	    
	  std::string hName =
	    Form("%s_%s_%s_%s",
		 anaTool->GetName( eta1Low, eta1Up, "Eta1").c_str(),
		 anaTool->GetName( eta2Low, eta2Up, "Eta2").c_str(),
		 anaTool->GetName( pt1Low , pt1Up , "Pt1" ).c_str(),
		 anaTool->GetName( pt2Low , pt2Low, "Pt2" ).c_str() );

	  std::string hName_pPb =
	    Form("h_dPhi_%s_%s", trigger_pPb.c_str(), hName.c_str() );
	  std::string hName_pp =
	    Form("h_dPhi_%s_%s", trigger_pp.c_str() , hName.c_str() );
	    
	  TCanvas* c_pPb =
	    static_cast<TCanvas*>
	    ( fIn_pPb->Get( Form("c_%s_data_pPb", hName_pPb.c_str() ) ) );
	  TCanvas* c_pp  =
	    static_cast<TCanvas*>
	    ( fIn_pp->Get( Form("c_%s_data_pp", hName_pp.c_str() ) ) );

	  TH1* h_pPb =
	    static_cast<TH1D*>( c_pPb->GetPrimitive( hName_pPb.c_str() ) );
	  styleTool->SetHStyle( h_pPb, 0 );
	  TH1* h_pp  =
	    static_cast<TH1D*>( c_pp->GetPrimitive( hName_pp.c_str() ) );
	  styleTool->SetHStyle( h_pp, 1 );

	  TCanvas c("c","c", 800, 600 );
	    
	  TLegend leg( 0.27, 0.41, 0.38, 0.52 );
	  styleTool->SetLegendStyle( &leg );
	  leg.AddEntry( h_pPb, "p+Pb");
	  leg.AddEntry( h_pp , "pp");

	  h_pPb->Draw("epsame");
	  h_pp->Draw("epsame");
	  leg.Draw("same");

	  if( h_pPb->GetMaximum() > h_pp->GetMaximum() ){
	    h_pPb->SetMaximum( h_pPb->GetMaximum() * 1.1 );
	    h_pp->SetMaximum ( h_pPb->GetMaximum() * 1.1 );
	  } else {
	    h_pPb->SetMaximum( h_pp->GetMaximum() * 1.1 );
	    h_pp->SetMaximum ( h_pp->GetMaximum() * 1.1 );
	  }
	    
	  drawTool->DrawLeftLatex
	    ( 0.13, 0.87,
	      GetLabel( eta1Low, eta1Up, "#eta_{1}" ).c_str(), 1 );
	  drawTool->DrawLeftLatex
	    ( 0.13, 0.82,
	      GetLabel( eta2Low, eta2Up, "#eta_{2}" ).c_str(), 1 );
	  drawTool->DrawLeftLatex
	    ( 0.13, 0.76,
	      GetLabel( pt1Low, pt1Up, "#it{p}_{T}^{1}" ).c_str(), 1 );
	  drawTool->DrawLeftLatex
	    ( 0.13, 0.69,
	      GetLabel( pt2Low, pt2Low, "#it{p}_{T}^{2}" ).c_str(), 1 );

	  drawTool->DrawAtlasInternal();
	  
	  SaveAsROOT( c, Form("h_dPhi_%s", hName.c_str() ) );

	  delete h_pPb;
	  delete h_pp;
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
