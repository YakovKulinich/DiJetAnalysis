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

  m_vRefTriggers =
    anaTool->vectorise
    ( config.GetValue( Form("refTriggers.%s",triggerMenu.c_str()), "" ), " ");

  m_vTholdPtTriggers =
    anaTool->vectoriseD
    ( config.GetValue( Form("tholdTriggers.%s",triggerMenu.c_str()), "" ), " ");

  m_vEffPtTriggers =
    anaTool->vectoriseD
    ( config.GetValue( Form("effTriggers.%s",triggerMenu.c_str()), "" ), " ");
  
  m_nTriggers = m_vTriggers.size();

  std::cout << "Using " << m_nTriggers << " triggers:" << std::endl;
  for( auto& s : m_vTriggers ){ std::cout << s << std::endl;}

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
      
    m_vTriggerEtaSpectEff.push_back
      ( new TH2D ( Form("h_etaSpectEff_%s", trigger.c_str() ), 
		   ";#eta;#it{p}_{T} [GeV]",
		   m_nVarFwdEtaBins, 0, 1,
		   m_nPtSpectBins,
		   m_ptSpectMin, m_ptSpectMax ) ) ;
    m_vTriggerEtaSpectEff.back()->GetXaxis()->
      Set( m_nVarFwdEtaBins, &( m_varFwdEtaBinning[0] ) );
    AddHistogram( m_vTriggerEtaSpectEff.back() );

    // -------- dPhi --------
    std::vector< int    > nDphiBins;
    std::vector< double > dPhiMin;
    std::vector< double > dPhiMax;


    boost::assign::push_back( nDphiBins )
      ( m_nVarFwdEtaBins )( m_nVarFwdEtaBins )
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
      new THnSparseD( Form("hs_dPhi_%s", trigger.c_str() ), "",
		      nDim, &nDphiBins[0], &dPhiMin[0], &dPhiMax[0] );
    hn->GetAxis(0)->Set( m_nVarFwdEtaBins, &( m_varFwdEtaBinning[0] ) );
    hn->GetAxis(1)->Set( m_nVarFwdEtaBins, &( m_varFwdEtaBinning[0] ) );
  
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

    ApplyCleaning ( vR_jets, v_isCleanJet );
    ApplyIsolation( 1.0, vR_jets );

    if( anaTool->DoPrint(m_ev) ) {
      std::cout << "\nEvent : " << m_ev << "    runN : " << runNumber
		<< "    has : " << vR_jets.size() << " jets" 
		<< std::endl; 
    }

    std::sort( vR_jets.begin(), vR_jets.end(),
	       anaTool->sortByDecendingPt );
    
    // SPECTRA AND ETA PHI PT MAPS
    // loop over passed triggers
    for( uint iT = 0; iT < m_nTriggers; iT++ ){

      // check if we have that trigger
      if( !mTriggerFired[iT] ) {
	mTrigFailed[iT]++;
	continue;	
      }

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

	m_vTriggerEtaPhiMap[iT]->Fill( jetEta, jetPhi );
	m_vTriggerEtaPtMap [iT]->Fill( jetEta, jetPt ); 
	// eta cut

	// convert positive eta to negative because
	// in pp it doesnt matter.
	// our histos run in negative eta
	// (due to pPb configuration)
	// the labels will be taken care of so it is ok
	double jetEtaAdj = AdjustEtaForPP( jetEta );
	
	if( jetPt >= m_vTholdPtTriggers[iT] ){
	  mTrigJetsTotal[iT] ++;
	}
	
	// saves computation time
	// can theoreticall leave out
	// will go to overflow bins in histos
	if( jetEtaAdj > -constants::FETAMIN ) { continue; }

	// above 5 is to exclude MB trigger
	// because its threshold is set at zero
	if( jetPt >= m_vEffPtTriggers[iT] && jetPt > 5 ){
	  mTrigJetsForward[iT]++;
	}
	
	m_vTriggerEtaSpect[iT]->Fill( jetEtaAdj, jetPt );
      } // end loop over jets

      mTrigPassed[iT]++;
    } // end loop over iT

    // EFFICIENCIES
    AnalyzeEff( vTrig_jets, vR_jets, mTriggerFired );
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
  // if we had no MB trigger fired
  // or the trigger jet collection is
  // empty, return.
  if( !mTriggerFired[ m_mbTriggerI ] ){ return; }
  if( vTrig_jets.empty() )            { return; }

  // sort trigger jets
  std::sort( vTrig_jets.begin(), vTrig_jets.end(),
	     anaTool->sortByDecendingPt );

  // take highest pt trigger jet in forward eta range
  TLorentzVector* tJet = NULL;
  for( auto& jet : vTrig_jets ){
    if( AdjustEtaForPP( jet.Eta() ) > -constants::FETAMIN )
      continue;
    tJet = &jet;
    break; // found highest forward trigger jet
  }    

  // didnt find a trigger jet in forward eta range
  if( !tJet ) return;

  // now fill histos for triggers that could have passed
  for( uint iT = 0; iT < m_nTriggers; iT++ ){
    double tJetPt = tJet->Pt()/1000;
    // now check if we should fill for that trigger
    if( tJetPt > m_vTholdPtTriggers[iT] ){ 
      // fill jets for that trigger
      for( auto& jet : vR_jets ){
	double jetEtaAdj = AdjustEtaForPP( jet.Eta() );
	// eta cut. want forward trigger jet
	if( jetEtaAdj > -constants::FETAMIN ) continue;
	m_vTriggerEtaSpectEff[iT]->
	  Fill( jetEtaAdj, jet.Pt()/1000. );
      } // end loop over jets
    }
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
  }
  
  m_fIn->Close();
}

void DiJetAnalysisData::PlotSpectra( std::vector< TH2* >& mTrigSpect,
				     const std::string& type){
  std::string yAxisTitle = "dN/d#it{p}_{T}";

  double ptSpectWidth =
    ( m_ptSpectMax - m_ptSpectMin ) / m_nPtSpectBins;
	          
  double lX0, lY0, lX1, lY1;
  if( m_is_pPb ){ lX0 = 0.45; lY0 = 0.54; lX1 = 0.76; lY1 = 0.67; }
  else          { lX0 = 0.56; lY0 = 0.54; lX1 = 0.86; lY1 = 0.67; }

  // use this as reference because
  // it should be in every file
  TH2* hMB   = mTrigSpect[ m_mbTriggerI ];
  int nXbins = hMB->GetNbinsX();
  
  std::vector< std::vector< TH1* > > vSpect;
  vSpect.resize( m_nTriggers );
  
  double max = -1;
  for( uint iT = 0; iT < m_nTriggers; iT++){
    std::string trigger = m_vTriggers[iT];

    for( int xBin = 1; xBin <= nXbins; xBin++ ){
      double etaMin = hMB->GetXaxis()->GetBinLowEdge( xBin );
      double etaMax = hMB->GetXaxis()->GetBinUpEdge ( xBin );

      TH1* h_etaSpect =
	mTrigSpect[iT]->
	ProjectionY( Form("h_%s_%s_%2.0f_Eta_%2.0f",
			  type.c_str(),
			  trigger.c_str(),
			  10*std::abs(etaMin),
			  10*std::abs(etaMax) ),
		     xBin, xBin );
      h_etaSpect->SetTitle( GetEtaLabel( etaMin, etaMax ).c_str() );
      h_etaSpect->GetYaxis()->SetTitle( yAxisTitle.c_str() );
      h_etaSpect->Scale( 1./ptSpectWidth );
      vSpect [iT].push_back( h_etaSpect );
      
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
    styleTool->SetLegendStyle( &leg, CT::StyleTools::lSS );

    int style = 0;
    for( int iE = 0; iE < nXbins; iE++ ){
      TH1* h = vSpect[iT][iE];
      styleTool->SetHStyle( h, style++, CT::StyleTools::hSS);
      h->Draw("epsame");
      h->SetMinimum( 1 );
      h->SetMaximum( max );
      leg.AddEntry( h, h->GetTitle() );
      h->SetTitle("");
    } // end loop over iE
    
    leg.Draw("same");
    drawTool->DrawAtlasInternalDataRight( 0, 0, CT::StyleTools::lSS, m_is_pPb ); 
    drawTool->DrawLeftLatex( 0.15, 0.17,
			      cLabel.c_str(),
			      CT::StyleTools::lSS, 1 );
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
    styleTool->SetLegendStyle( &leg, CT::StyleTools::lSS );
    
    int style = 0;
    for( uint iT = 0; iT < m_nTriggers; iT++ ){
      TH1* h = vSpect[iT][iE];
      styleTool->SetHStyle( h, style++, CT::StyleTools::hSS);
      h->Draw("epsame");
      h->SetMinimum( 1 );
      h->SetMaximum( max );
      leg.AddEntry( h, m_vTriggers[iT].c_str() );
    } // end loop over iT
    
    leg.Draw("same");
    drawTool->DrawAtlasInternalDataRight( 0, 0, CT::StyleTools::lSS, m_is_pPb ); 
    drawTool->DrawLeftLatex( 0.15, 0.17,
			      cLabel.c_str(),
			      CT::StyleTools::lSS, 1 );

    SaveAsAll( c, type, cName );
  } // end loop over iE

  // delete
  for( uint iT = 0; iT < m_nTriggers; iT++ ){
    for( int iE = 0; iE < nXbins; iE++ )
      { delete vSpect[iT][iE]; }
  }
}

void DiJetAnalysisData::PlotEfficiencies( std::vector< TH2* >& mTrigSpect,
					  const std::string& type ){}

void DiJetAnalysisData::PlotDeltaPhi( std::vector< THnSparse* >& vhn ){
  /*
  std::vector< TH1* > vDphi;
    
  for( uint iT = 0; iT < m_nTriggers; iT++ ){
    THnSparse* hn = vhn[iT];

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
	    pt2Axis->SetRange( pt2Bin, pt2Bin );

	    
	    
	    // Take projection onto the dPhi axis
	    TH1* hDphi = hn->Projection( 4 );
	    styleTool->SetHStyle( hDphi, 0, CT::StyleTools::hSS );
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

	    
	    anaTool->GetBinRange()
	    hDphi->SetName
	      ("h_dPhi_%s_%s_%s_s",
	       anaTool->GetName( eta1Low, eta1Up ),
	       anaTool->GetName( eta2Low, eta2Up ),
	       anaTool->GetName( pt1Low , pt1Up  ),
	       anaTool->GetName( pt2Low , pt2Up  ) );

	    TCanvas c( "c", hDphi->GetName(), 800, 600 );
	    hDphi->Draw();

	    drawTool->DrawLeftLatex
	      ( 0.18, 0.88,
		anaTool->GetUnit( eta1Low, eta1Up, "#eta" ),
		,CT::StyleTools::lSS, 1 );
	    drawTool->DrawLeftLatex
	      ( 0.88, 0.71,
		anaTool->GetUnit( eta2Low, eta2Up, "#eta" ),
		,CT::StyleTools::lSS, 1 );
	    drawTool->DrawLeftLatex
	      ( 0.88, 0.74,
		anaTool->GetUnit( pt1Low, pt1Up, "#it{p}_{T}^{1}" ),
		,CT::StyleTools::lSS, 1 );
	    drawTool->DrawLeftLatex
	      ( 0.88, 0.67,
		anaTool->GetUnit( pt2Low, pt2Up, "#it{p}_{T}^{2}" ),
		,CT::StyleTools::lSS, 1 );
	    
	    drawTool->DrawAtlasInternalDataRight
	      ( 0, 0, CT::StyleTools::lSS, m_is_pPb ); 
	  }
      	}
      }      
    }
  } // end loop over iT
  */
}

void DiJetAnalysisData::PlotEtaPhiPtMap( std::vector< TH2* >& vTrigHin ){
  TCanvas c_map("c_map","c_map",800,600);
  
  for( uint iT = 0; iT < m_nTriggers; iT++){
    TH2* h = vTrigHin[iT];
    h->Draw("col");
    styleTool->SetHStyle( h, 0, CT::StyleTools::hSS);
    drawTool->DrawAtlasInternalDataRight
      ( 0, -0.55, CT::StyleTools::lSS, m_is_pPb );  
    
    SaveAsAll( c_map, h->GetName() );
  }
}

