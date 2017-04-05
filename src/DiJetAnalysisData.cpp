#include <TLine.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TEnv.h>
#include <TCanvas.h>
#include <TEfficiency.h>
#include <TGraphAsymmErrors.h>

#include <iostream>
#include <cmath>

#include "MyRoot.h"

#include "DiJetAnalysisData.h"

DiJetAnalysisData::DiJetAnalysisData() : DiJetAnalysisData( true, true )
{}

DiJetAnalysisData::DiJetAnalysisData( bool isData,
				      bool is_pPb )
  : DiJetAnalysis( isData, is_pPb)
{
  //========== Set Histogram Binning =============
  /*
  // Efficiencies
  m_nPtEffBins   = 50; 
  m_ptEffMin     = 0;
  m_ptEffMax     = m_nPtEffBins;
  */
  //==================== Cuts ====================    
}

DiJetAnalysisData::~DiJetAnalysisData(){}

void DiJetAnalysisData::Initialize(){
  DiJetAnalysis::Initialize();
  
  m_fNameIn = m_is_pPb ?
    "/home/yakov/Projects/atlas/data/pPb.root" :
    "/home/yakov/Projects/atlas/data/pp.root"  ;
}

//---------------------------------
//            Read Data
//---------------------------------
void DiJetAnalysisData::RunOverTreeFillHistos( int nEvents, 
					       int startEvent ){  
  LoadTriggers();
  SetupHistograms();
  ProcessEvents( nEvents, startEvent );
  SaveOutputsFromTree();
}

//---------------------------------
//            Plot Data
//---------------------------------
void DiJetAnalysisData::ProcessPlotHistos(){
  LoadTriggers();
  LoadHistograms();

  std::string cfNameOut = m_dirOut + "/c_myOut_" + m_labelOut + ".root";
  m_fOut = new TFile( cfNameOut.c_str(),"RECREATE");
    
  PlotEtaPhiPtMap( m_mTriggerEtaPhiMap );
  PlotEtaPhiPtMap( m_mTriggerEtaPtMap  );

  PlotSpectra( m_mTriggerEtaSpect, "spect" );

  PlotEfficiencies( m_mTriggerEtaSpectEff, "eff" );
  
  std::cout << "DONE! Closing " << cfNameOut << std::endl;
  m_fOut->Close();
  std::cout << "......Closed  " << cfNameOut << std::endl;
}

//---------------------------------
//            Read Data
//---------------------------------

void DiJetAnalysisData::LoadTriggers(){
  std::string configName   = "config/configJetsFwd";
  std::string configSuffix = m_is_pPb ? "_pPb.cfg" : "_pp.cfg";
  configName += configSuffix;
  
  TEnv config;
  config.ReadFile( configName.c_str(), EEnvLevel(0));

  std::string triggerMenu = config.GetValue("triggerMenu","");
  m_vTriggers =
    AnalysisTools::vectorise( config.
			      GetValue( Form("triggers.%s",
					     triggerMenu.c_str() )
					,"") , " " );
  
  for( auto& trigger : m_vTriggers ){
    if( trigger.find("_mb_") != std::string::npos ){ 
      m_mbTrigger = trigger; 
      std::cout << "Found Min Bias trigger " << m_mbTrigger << std::endl;
      break;
    }
  }

  v_tJetPt.push_back(10);
  v_tJetPt.push_back(15);
  v_tJetPt.push_back(25);
  v_tJetPt.push_back(35);
  v_tJetPt.push_back(45);
  v_tJetPt.push_back(55);
  v_tJetPt.push_back(65);

  std::string triggerLabel;

  // associate pt cut to trigger
  for( auto pt : v_tJetPt ){
    triggerLabel = "j" + std::to_string(pt);
    // find trigger with say, j10_ion in it
    // associate it to a jet pt cut
    for( auto& trigger : m_vTriggers ){
      if( trigger.find(triggerLabel ) != std::string::npos ){
	m_tJetPtTrigger[pt] = trigger;
      }
    }
  }

  for( auto ptTrig = m_tJetPtTrigger.begin();
       ptTrig!= m_tJetPtTrigger.end();
       ptTrig++ ){
    std::cout << "For " << ptTrig->first << " GeV Trigger Jet : " 
	      << ptTrig->second << std::endl;
  }
}

void DiJetAnalysisData::SetupHistograms(){

  for( auto& trigger : m_vTriggers ){
    std::cout << "Making - " << trigger << " histograms " << std::endl;          
    m_mTriggerEtaPhiMap[ trigger ] = 
      new TH2D( Form("h_etaPhiMap_%s", trigger.c_str() ), 
		";#eta;#phi",
		m_nEtaBins, m_etaMin, m_etaMax,
		m_nPhiBins, m_phiMin, m_phiMax ) ;
    AddHistogram( m_mTriggerEtaPhiMap[ trigger ] );
     
    m_mTriggerEtaPtMap[ trigger ] = 
      new TH2D( Form("h_etaPtMap_%s", trigger.c_str() ), 
		";#eta;#it{p}_{T}",
		m_nEtaBins, m_etaMin, m_etaMax,
		m_nPtBins , m_ptMin , m_ptMax) ;
    AddHistogram( m_mTriggerEtaPtMap[ trigger ] );

    m_mTriggerEtaSpect[ trigger ] = 
      new TH2D( Form("h_etaSpect_%s", trigger.c_str() ), 
		";#eta;#it{p}_{T} [GeV];",
		m_nVarEtaBins, 0, 1,
		m_nPtSpectBins,
		m_ptSpectMin, m_ptSpectMax ) ;
    m_mTriggerEtaSpect[ trigger ]->GetXaxis()->
      Set( m_nVarEtaBins, &( m_varEtaBinning[0] ) );
    AddHistogram( m_mTriggerEtaSpect[ trigger ] );
      
    m_mTriggerEtaSpectEff[ trigger ] = 
      new TH2D( Form("h_etaSpectEff_%s", trigger.c_str() ), 
		";#eta;#it{p}_{T} [GeV]",
		m_nVarEtaBins, 0, 1,
		m_nPtSpectBins,
		m_ptSpectMin, m_ptSpectMax ) ;
    m_mTriggerEtaSpectEff[ trigger ]->GetXaxis()->
      Set( m_nVarEtaBins, &( m_varEtaBinning[0] ) );
    AddHistogram( m_mTriggerEtaSpectEff[trigger ] );

    m_mTriggerRunPrescale[ trigger ] = 
      new TH2D( Form("h_runPrescale_%s", trigger.c_str() ), 
		";Run No;Prescale",
		1700, 312600, 314300, 1999, 1, 2000) ;
    AddHistogram( m_mTriggerRunPrescale[ trigger ] );
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

  std::map< std::string, bool  > mTriggerFired;
  std::map< std::string, float > mTriggerPrescale;

  int    runNumber     = 0;
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

  for( auto& trigger : m_vTriggers ){
    m_tree->SetBranchAddress( Form("passed_%s",trigger.c_str() ),
			      &mTriggerFired[trigger] );
    m_tree->SetBranchAddress( Form("prescale_%s", trigger.c_str() ),
			      &mTriggerPrescale[trigger] );
  }

  m_tree->SetBranchAddress( "runNumber"   , &runNumber );  

  std::map< std::string, int > mTrigPassed;
  std::map< std::string, int > mTrigFailed;
  
  std::map< std::string, int > mTrigJetsForward;
  std::map< std::string, int > mTrigJetsTotal;

  // map of minimum pt for counting number of jets
  // in the near 1 efficiency region of that trigger
  std::map< std::string, double > mTrigPtThreshold;
  // pt min for mb is 0
  mTrigPtThreshold[ m_mbTrigger ] = 0; 
  
  for( auto & ptTrig : m_tJetPtTrigger ){
    mTrigPassed     [ ptTrig.second ] = 0;
    mTrigFailed     [ ptTrig.second ] = 0;
    mTrigJetsTotal  [ ptTrig.second ] = 0;
    // set the threshold for near 1 eff at N+5 
    // where N is from jN
    mTrigPtThreshold[ ptTrig.second ] = ptTrig.first + 5;
  }

  
  // n events
  int nEventsTotal = m_tree->GetEntries();

  nEvents = nEvents > 0 ? nEvents : nEventsTotal;
  startEvent = startEvent < nEventsTotal ? startEvent : nEventsTotal - 1;
  
  int endEvent = startEvent + nEvents < nEventsTotal ?
					startEvent + nEvents : nEventsTotal;

  // event loop
  for( int ev = startEvent; ev < endEvent; ev++ ){
    m_tree->GetEntry( ev );

    ApplyCleaning ( vR_jets, v_isCleanJet );
    ApplyIsolation( 1.0, vR_jets );

    if( AnalysisTools::DoPrint(ev) ) {
      std::cout << "\nEvent : " << ev << "    runN : " << runNumber
		<< "    has : " << vR_jets.size() << " jets" 
		<< std::endl; 
    }
    
    // SPECTRA AND ETA PHI PT MAPS
    // loop over passed triggers
    for( auto& trigger : m_vTriggers ){
      // check if we have that trigger
      if( !mTriggerFired[trigger] ) {
	mTrigFailed[ trigger ]++;
	continue;	
      }

      // loop over jets 
      for( auto& jet : vR_jets ){
	// cleaned jets have px, py, pz set to 0
	if( jet.Pt() == 0 ) continue;

	// ETA-PHI
	double jetEta = jet.Eta();
	double jetPhi = jet.Phi();
	double jetPt = jet.Pt()/1000.;

	m_mTriggerEtaPhiMap[ trigger ]->Fill( jetEta, jetPhi );
	m_mTriggerEtaPtMap [ trigger ]->Fill( jetEta, jetPt ); 
	// eta cut

	// convert positive eta to negative because
	// in pp it doesnt matter.
	// our histos run in negative eta
	// (due to pPb configuration)
	// the labels will be taken care of so it is ok
	double jetEtaAdj = AdjustEtaForPP( jetEta );
	
	if( jetPt >= mTrigPtThreshold[ trigger ] ){
	  mTrigJetsTotal[ trigger ] ++;
	}
	
	// saves computation time
	// can theoreticall leave out
	// will go to overflow bins in histos
	if( jetEtaAdj > -constants::FETAMIN ) { continue; }

	if( jetPt >= mTrigPtThreshold[ trigger ] ){
	  mTrigJetsForward[ trigger ] ++;
	}
	
	m_mTriggerEtaSpect[ trigger ]->Fill( jetEtaAdj, jetPt );
	// fill mb efficiency histo also
	if( !trigger.compare(m_mbTrigger) ){
	  m_mTriggerEtaSpectEff[m_mbTrigger]->Fill( jetEtaAdj, jetPt );
	}
      } // end loop over jets
      mTrigPassed[ trigger ]++;
    } // end loop over triggers
    // EFFICIENCIES
    ProcessEfficiencies( vTrig_jets, vR_jets, mTriggerFired );
  } // end event loop
  std::cout << "DONE! Has: " << nEventsTotal << " events." << std::endl;

  for( auto& th : m_mTriggerEtaSpect ){
    std::cout << th.first << "   has: " << mTrigJetsForward[ th.first ]
	      << " FWD jets, " << mTrigJetsTotal[ th.first ]
	      << " TOTAL jets, with " << mTrigPassed[ th.first ]
	      << " evPassed and "<< mTrigFailed[ th.first ]
	      << " failed. Total = "
	      << mTrigPassed[ th.first ] + mTrigFailed[ th.first ]
	      << std::endl;
  }

  m_fIn->Close();
}

void DiJetAnalysisData::
ProcessEfficiencies( std::vector< TLorentzVector >& vTrig_jets,
		     std::vector< TLorentzVector >& vR_jets,
		     std::map< std::string, bool >& mTriggerFired ){
  // only do if had min bias trigger
  bool haveMinBiasTrigger = mTriggerFired[ m_mbTrigger ] ? true : false;
  
  
  if( !haveMinBiasTrigger ) return;
  if( vTrig_jets.empty()  ) return;

  // sort trigger jets
  std::sort( vTrig_jets.begin(), vTrig_jets.end(),
	     AnalysisTools::sortByDecendingPt );

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
  for( auto ptTrig : m_tJetPtTrigger){
    double tJetPt = tJet->Pt()/1000;
    // now check if we should fill for that trigger
    if( tJetPt > ptTrig.first ){ 
      // fill jets for that trigger
      for( auto& jet : vR_jets ){
	// eta cut
	double jetEtaAdj = AdjustEtaForPP( jet.Eta() );
	if( jetEtaAdj > -constants::FETAMIN ) continue;
	m_mTriggerEtaSpectEff[ ptTrig.second ]->
	  Fill( jetEtaAdj, jet.Pt()/1000. );
      } // end loop over jets
    }
  } // end loop over triggers 
}

//---------------------------------
//            Plot Data
//---------------------------------

void DiJetAnalysisData::LoadHistograms(){
  m_fIn = TFile::Open( m_rootFname.c_str() ); 

  for( auto& trigger : m_vTriggers ){
    m_mTriggerEtaPhiMap   [ trigger ] =
      static_cast< TH2D* >
      ( m_fIn-> Get( Form("h_etaPhiMap_%s", trigger.c_str() ) ) );
    m_mTriggerEtaPhiMap   [ trigger ]->SetDirectory(0);
    m_mTriggerEtaPtMap    [ trigger ] =
      static_cast< TH2D* >
      ( m_fIn->Get( Form("h_etaPtMap_%s", trigger.c_str() ) ) );
    m_mTriggerEtaPtMap    [ trigger ]->SetDirectory(0);

    m_mTriggerEtaSpect [ trigger ] =
      static_cast< TH2D *>
      ( m_fIn->Get( Form("h_etaSpect_%s", trigger.c_str() ) ) );
    m_mTriggerEtaSpect [ trigger ]->SetDirectory(0);
    m_mTriggerEtaSpectEff[ trigger ] =
      static_cast< TH2D* >
      ( m_fIn->Get( Form("h_etaSpectEff_%s", trigger.c_str() ) ) );
    m_mTriggerEtaSpectEff[ trigger ]->SetDirectory(0);
  }
  
  m_fIn->Close();
}

void DiJetAnalysisData::
PlotSpectra( std::map< std::string, TH2* >& mTrigSpect,
	     const std::string& type){
  // if there are none, return
  if( !mTrigSpect.size() ){ return; }  

  std::string yAxisTitle = "dN/d#it{p}_{T}";

  double ptSpectWidth =
    ( m_ptSpectMax - m_ptSpectMin ) / m_nPtSpectBins;
	          
  double lX0, lY0, lX1, lY1;

  if( m_is_pPb ){ lX0 = 0.45; lY0 = 0.54; lX1 = 0.76; lY1 = 0.67; }
  else          { lX0 = 0.56; lY0 = 0.54; lX1 = 0.86; lY1 = 0.67; }
  
  std::vector< TH1* > vSpect;
  std::multimap< std::string, TH1* > mmTrigSpect;

  double max = -1;

  for( int xBin = 1; xBin <= mTrigSpect[ m_mbTrigger ]->GetNbinsX(); xBin++ ){
  
    TCanvas c_spect("c_spect","c_spect",800,600);
    c_spect.SetLogy();
    
    TLegend l_etaSpect( lX0, lY0, lX1, lY1 );
    StyleTools::SetLegendStyle( &l_etaSpect, StyleTools::lSS );
    l_etaSpect.SetFillStyle(0);
  
    int trigStyle = 0;
    // should all be the same
    double etaMin = mTrigSpect[ m_mbTrigger ]->
      GetXaxis()->GetBinLowEdge( xBin );
    double etaMax = mTrigSpect[ m_mbTrigger ]->
      GetXaxis()->GetBinUpEdge( xBin );

    for( auto& tH : mTrigSpect ){
      TH1* h_etaSpect =
	tH.second->
	ProjectionY( Form("h_%s_%s_%2.0f_Eta_%2.0f",
			  type.c_str(),
			  tH.first.c_str(),
			  10*std::abs(etaMin),
			  10*std::abs(etaMax) ),
		     xBin, xBin );
      StyleTools::SetHStyle( h_etaSpect, trigStyle++, StyleTools::hSS);

      h_etaSpect->GetYaxis()->SetTitle( yAxisTitle.c_str() );
      h_etaSpect->Scale( 1./ptSpectWidth );
      vSpect.push_back( h_etaSpect );

      int etaStyle = xBin - 1; // increments with thee x/etaBin
      TH1* h_trigSpect = static_cast<TH1*>( h_etaSpect->Clone() );
      h_trigSpect->SetTitle( GetEtaLabel(etaMin, etaMax).c_str() );
      StyleTools::SetHStyle( h_trigSpect, etaStyle, StyleTools::hSS);
      vSpect.push_back( h_trigSpect );
      
      mmTrigSpect.insert( std::pair< std::string, TH1* >
			    ( tH.first, h_trigSpect ) );
      
      h_etaSpect->Draw("epsame");
      if( max < h_etaSpect->GetMaximum() )
	{ max = h_etaSpect->GetMaximum(); }

      l_etaSpect.AddEntry( h_etaSpect, tH.first.c_str() );
    }
    
    double power = log10(max);
    power = std::ceil(power);
    max = pow( 10, power );
    for( auto& h : vSpect ){
      h->SetMinimum( 1 );
      h->SetMaximum( max );
    }
  
    l_etaSpect.Draw();
    DrawTools::DrawAtlasInternalDataRight( 0, 0, StyleTools::lSS, m_is_pPb ); 
    DrawTools::DrawLeftLatex( 0.45, 0.87,
			      GetEtaLabel( etaMin, etaMax).c_str(),
			      StyleTools::lSS, 1 );

    SaveAsAll( c_spect, type, "", "Eta", std::abs(etaMin)*10, std::abs(etaMax)*10 );    
  } // end loop over eta bins
  
  // now draw for each trigger, spectra in various eta bins
  for( auto& trig : m_vTriggers ){ // loop over triggers
    TCanvas c_spect("c_spect","c_spect",800,600);
    c_spect.SetLogy();

    TLegend l_trigSpect( 0.7, lY0, lX1, lY1 );
    StyleTools::SetLegendStyle( &l_trigSpect, StyleTools::lSS );
    l_trigSpect.SetFillStyle(0);

    auto rangeTrigSpect = mmTrigSpect.equal_range( trig );
    for( auto it = rangeTrigSpect.first;
	 it != rangeTrigSpect.second;
	 it++ ){
      std::cout << "For: "<< it->first << " -> "
		<< it->second->GetName() << std::endl; 
      it->second->Draw("epsame");
      it->second->SetMinimum( 1 );
      it->second->SetMaximum( max );
      l_trigSpect.AddEntry( it->second, it->second->GetTitle() );
      it->second->SetTitle("");
    }
    l_trigSpect.Draw("same");
    DrawTools::DrawAtlasInternalDataRight( 0, 0, StyleTools::lSS, m_is_pPb ); 
    DrawTools::DrawLeftLatex( 0.15, 0.17,
			      trig.c_str(),
			      StyleTools::lSS, 1 );
 
    SaveAsAll( c_spect, type, trig );
  }// end loop over triggers

  for( auto& h : vSpect ){ delete h; }
}

void DiJetAnalysisData::
PlotEfficiencies( std::map< std::string, TH2* >& mTrigSpect,
		  const std::string& type ){
  // if there is no mb trigger, return
  if( !mTrigSpect[ m_mbTrigger ] ){ return; }  

  double lX0 = 0.13;
  double lY0 = 0.75;
  double lX1 = 0.39;
  double lY1 = 0.87;

  double xMin = 0.;   //mTrigSpect[ m_mbTrigger ]->GetYaxis()->GetXmin();
  double xMax = 90.;  //mTrigSpect[ m_mbTrigger ]->GetYaxis()->GetXmax();

  std::string gTitle = ";#it{p}_{T} [GeV];#it{#varepsilon}_{Trigger}";
  
  // For later
  // in pp we have heavily prescaled triggers
  // dont divide by MB. Use something "hard coded"
  std::map< int, int > mTrigRefTrig;
  mTrigRefTrig[ 10 ] = -1; // mb
  mTrigRefTrig[ 15 ] = -1; // mb
  mTrigRefTrig[ 25 ] = 10; 
  
  std::vector< TH1* > vMbTrigger;
  std::vector< TH1* > vSpect;
  std::vector< TGraphAsymmErrors* > vEffGrf;
  
  std::multimap< std::string, TGraphAsymmErrors* > mmTrigEffGrf;
  
  for( int xBin = 1; xBin <= mTrigSpect[ m_mbTrigger ]->GetNbinsX(); xBin++ ){
  
    TCanvas c_eff("c_eff","c_eff",800,600);
    StyleTools::SetCStyleEff( c_eff, xMin, m_effMin, xMax, m_effMax, gTitle );
    
    TLegend l_eff( lX0, lY0, lX1, lY1);
    StyleTools::SetLegendStyle( &l_eff, StyleTools::lSS );
    l_eff.SetFillStyle(0);

    int trigStyle = 0;
    // should all be the same
    double etaMin = mTrigSpect[ m_mbTrigger ]->
      GetXaxis()->GetBinLowEdge( xBin );
    double etaMax = mTrigSpect[ m_mbTrigger ]->
      GetXaxis()->GetBinUpEdge( xBin );

    // local inside the loop. used only on a per eta-bin
    // basis. the actual projections are saved to the vectors
    // that are global in this function.
    std::map< std::string, TH1* > mSpectTemp;
    
    TH1* h_mbTrig =
      mTrigSpect[ m_mbTrigger ]->
      ProjectionY( Form("h_%s_%s_%2.0f_Eta_%2.0f",
			type.c_str(),
			m_mbTrigger.c_str(),
			10*std::abs(etaMin),
			10*std::abs(etaMax) ),
		   xBin, xBin );
    vMbTrigger.push_back( h_mbTrig );
  
    for( auto& trigSpect : mTrigSpect ){
      std::string trigger = trigSpect.first;
      mSpectTemp[ trigger ] =
	trigSpect.second->
	ProjectionY( Form("h_%s_%s_%2.0f_Eta_%2.0f",
			  type.c_str(),
			  trigger.c_str(),
			  10*std::abs(etaMin),
			  10*std::abs(etaMax) ),
		     xBin, xBin);
      vSpect.push_back( mSpectTemp[ trigger ] );
    }
  
    for( auto& trigSpect : mSpectTemp ){
      std::string trigger = trigSpect.first;
      
      TGraphAsymmErrors* g_etaEff = new TGraphAsymmErrors();
      g_etaEff->SetName ( Form("gr_%s_%s", type.c_str(), trigger.c_str() ) );
      g_etaEff->SetTitle( gTitle.c_str() );
      vEffGrf.push_back( g_etaEff );
      // dont draw the MB efficiencies, its just 1
      // dont add mb efficiencies to legend
      if( !trigger.compare(m_mbTrigger) ) continue;
      l_eff.AddEntry( g_etaEff, trigger.c_str() );

      TH1*            denom = h_mbTrig;
      std::string denomName = m_mbTrigger;
      
      //-------------------- pp ---------------------
      // in pp we have heavily prescaled triggers
      // dont divide by MB. Use something "hard coded"
      for( auto it = v_tJetPt.begin(); it!= v_tJetPt.end(); it++  ){
	if( !m_tJetPtTrigger[*it].compare( trigger ) ){
	  int refTrigPt = mTrigRefTrig[ *it ];
	  if( refTrigPt < 0 ) {
	    break;
	  }
	  denom     = mSpectTemp[ m_tJetPtTrigger[ refTrigPt ] ];
	  denomName = m_tJetPtTrigger[ refTrigPt ];
	  std::cout << "refTrigPt " << refTrigPt << "  " << denomName << std::endl; 
	}
      }

      std::cout << "-----For: "     << trigger
		<< " dividing by: " << denomName << std::endl;
      
      g_etaEff->Divide( trigSpect.second, denom,
		       "cl=0.683 b(1,1) mode" );
      StyleTools::SetHStyle( g_etaEff, trigStyle++, StyleTools::hSS);

      int etaStyle = xBin - 1; // increments with thee x/etaBin
      TGraphAsymmErrors* g_trigEff =
	static_cast<TGraphAsymmErrors*>( g_etaEff->Clone() );
      g_trigEff->SetTitle( GetEtaLabel(etaMin, etaMax).c_str() );
      StyleTools::SetHStyle( g_trigEff, etaStyle, StyleTools::hSS);
      vEffGrf.push_back( g_trigEff );
      
      mmTrigEffGrf.insert( std::pair< std::string, TGraphAsymmErrors* >
			  ( trigger, g_trigEff ) );
      
      g_etaEff->Draw("p"); 
    } // end loop over triggers
  
    l_eff.Draw();

    TLine line( xMin, 1, xMax, 1);
    line.Draw();

    DrawTools::DrawAtlasInternalDataRight( 0, 0,  StyleTools::lSS, m_is_pPb ); 
    DrawTools::DrawRightLatex( 0.88, 0.17,
			       GetEtaLabel( etaMin, etaMax).c_str(),
			       StyleTools::lSS, 1 );

    SaveAsAll( c_eff, type, "", "Eta", std::abs(etaMin)*10, std::abs(etaMax)*10 );
  } // end loop over eta

  // now draw for each trigger, spectra in various eta bins
  for( auto& trig : m_vTriggers ){ // loop over triggers
    TCanvas c_eff("c_eff","c_eff",800,600);
    StyleTools::SetCStyleEff( c_eff, xMin, m_effMin, xMax, m_effMax, gTitle );
      
    if( !trig.compare( m_mbTrigger ) ){ continue; }
    
    TLegend l_trigEff( lX0, lY0, lX1, lY1);
    StyleTools::SetLegendStyle( &l_trigEff, StyleTools::lSS );
    l_trigEff.SetFillStyle(0);
    
    auto rangeTrigSpect = mmTrigEffGrf.equal_range( trig );
    for( auto it = rangeTrigSpect.first;
	 it != rangeTrigSpect.second;
	 it++ ){
      std::cout << "For: "<< it->first << " -> "
		<< it->second->GetName() << std::endl;

      it->second->Draw("p"); 
      l_trigEff.AddEntry( it->second, it->second->GetTitle() );
    }

    l_trigEff.Draw();

    TLine line( xMin, 1, xMax, 1);
    line.Draw();
    
    DrawTools::DrawAtlasInternalDataRight( 0, 0, StyleTools::lSS, m_is_pPb ); 
    DrawTools::DrawRightLatex( 0.88, 0.17,
			       trig.c_str(),
			       StyleTools::lSS, 1 );
 
    SaveAsAll( c_eff, type, trig );
  }// end loop over triggers

  for( auto& h : vMbTrigger ){ delete h; }
  for( auto& h : vSpect     ){ delete h; }
  for( auto& g : vEffGrf    ){ delete g; }
}

void DiJetAnalysisData::
PlotEtaPhiPtMap( std::map< std::string, TH2* >& mTrigHIN ){
  TCanvas c_map("c_map","c_map",800,600);

  for( auto& trigH : mTrigHIN ){
    trigH.second->Draw("col");
    StyleTools::SetHStyle( trigH.second, 0, StyleTools::hSS);
    DrawTools::DrawAtlasInternalDataRight( 0, -0.55,
					   StyleTools::lSS, m_is_pPb );  

    SaveAsAll( c_map, trigH.second->GetName() );
  }
}
