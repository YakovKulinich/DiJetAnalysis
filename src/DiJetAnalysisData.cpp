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

DiJetAnalysisData::DiJetAnalysisData()
{}

DiJetAnalysisData::DiJetAnalysisData( bool isData,
				      bool is_pPb )
  : DiJetAnalysis( isData, is_pPb)
{}

DiJetAnalysisData::~DiJetAnalysisData(){}

void DiJetAnalysisData::Initialize(){
  DiJetAnalysis::Initialize();
  
  m_fNameIn = m_is_pPb ?
    "/home/yakov/Projects/atlas/data/pPb.root" :
    "/home/yakov/Projects/atlas/data/pp.root"  ;


  // Efficiencies
  m_nPtEffBins   = 50; 
  m_ptEffMin     = 0;
  m_ptEffMax     = m_nPtEffBins;
}

//---------------------------------
//            Read Data
//---------------------------------
void DiJetAnalysisData::RunOverTreeFillHistos( int nEvents, 
					       int startEvent ){  
  LoadTriggers();
  SetupHistograms();
  ProcessEvents( nEvents, startEvent );
  SaveOutputs();
}

//---------------------------------
//            Plot Data
//---------------------------------
void DiJetAnalysisData::ProcessPlotHistos(){
  LoadTriggers();
  LoadHistograms();

  std::string cfNameOut = m_dirOut + "/c_myOut" + m_labelOut + ".root";
  m_fOut = new TFile( cfNameOut.c_str(),"RECREATE");
  
  PlotSpectra( m_mTriggerEtaSpect );

  PlotEfficiencies( m_mTriggerEtaEff );
  
  PlotEtaPhiPtMap( m_mTriggerEtaPhi );
  PlotEtaPhiPtMap( m_mTriggerEtaPt  );
  
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
      mbTrigger = trigger; 
      std::cout << "Found Min Bias trigger " << mbTrigger << std::endl;
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
    m_mTriggerEtaSpect[ trigger ] = 
      new TH2D( Form("h_etaSpect_%s", trigger.c_str() ), 
		";#eta;#it{p}_{T} [GeV];dN/d#it{p}_{T}",
		m_nEtaForwardBinsCoarse,
		m_etaForwardMin, m_etaForwardMax,
		m_nPtSpectBins,
		m_ptSpectMin, m_ptSpectMax ) ;
    AddHistogram( m_mTriggerEtaSpect[ trigger ] );
      
    m_mTriggerEtaEff[ trigger ] = 
      new TH2D( Form("h_etaEff_%s", trigger.c_str() ), 
		";#eta;#it{p}_{T} [GeV]",
		m_nEtaForwardBinsCoarse,
		m_etaForwardMin, m_etaForwardMax,
		m_nPtEffBins,
		m_ptEffMin, m_ptEffMax ) ;
    AddHistogram( m_mTriggerEtaEff[trigger ] );
      
    m_mTriggerRunPrescale[ trigger ] = 
      new TH2D( Form("h_runPrescale_%s", trigger.c_str() ), 
		";Run No;Prescale",
		1700, 312600, 314300, 1999, 1, 2000) ;
    AddHistogram( m_mTriggerRunPrescale[ trigger ] );
    
    m_mTriggerEtaPhi[ trigger ] = 
      new TH2D( Form("h_etaPhi_%s", trigger.c_str() ), 
		";#eta;#phi",
		m_nEtaBins, m_etaMin, m_etaMax,
		m_nPhiBins, m_phiMin, m_phiMax ) ;
    AddHistogram( m_mTriggerEtaPhi[ trigger ] );
     
    m_mTriggerEtaPt[ trigger ] = 
      new TH2D( Form("h_etaPt_%s", trigger.c_str() ), 
		";#eta;#it{p}_{T}",
		m_nEtaBins, m_etaMin, m_etaMax,
		m_nPtBins , m_ptMin , m_ptMax) ;
    AddHistogram( m_mTriggerEtaPt[ trigger ] );
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

  std::map< std::string, bool  > m_mTriggerFired;
  std::map< std::string, float > m_mTriggerPrescale;

  int    runNumber     = 0;

  //----------------------------------------
  //  Open file and tree, Fill histograms
  //----------------------------------------
  std::cout << "fNameIn: " << m_fNameIn << std::endl;
  
  m_fIn = TFile::Open( m_fNameIn.c_str() );
  m_tree = (TTree*) m_fIn->Get( "tree" );

  // Connect to tree
  m_tree->SetBranchAddress( "vTrig_jets"  , &p_vTrig_jets );
  m_tree->SetBranchAddress( "vR_C_jets"   , &p_vR_jets );
  m_tree->SetBranchAddress( "v_isCleanJet", &p_v_isCleanJet );

  for( auto& trigger : m_vTriggers ){
    m_tree->SetBranchAddress( Form("passed_%s",trigger.c_str() ),
			      &m_mTriggerFired[trigger] );
    m_tree->SetBranchAddress( Form("prescale_%s", trigger.c_str() ),
			      &m_mTriggerPrescale[trigger] );
  }

  m_tree->SetBranchAddress( "runNumber"   , &runNumber );  

  // n events
  int maxEvents = m_tree->GetEntries();

  nEvents = nEvents > 0 ? nEvents : maxEvents;
  startEvent = startEvent < maxEvents ? startEvent : maxEvents - 1;
  
  int endEvent = startEvent + nEvents < maxEvents ? startEvent + nEvents : maxEvents;

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
    
    // loop over jets 
    for( auto& jet : vR_jets ){
      // cleaned jets have px, py, pz set to 0
      if( jet.Pt() == 0 ) continue;
      // SPECTRA AND ETA PHI PT MAPS
      // loop over passed triggers
      for( auto& trigger : m_vTriggers ){
	// check if we have that trigger
	if( !m_mTriggerFired[trigger] ) continue;
	// ETA-PHI
	double jetEta = jet.Eta();
	double jetPhi = jet.Phi();
	double jetPt = jet.Pt()/1000.;

	m_mTriggerEtaPhi[ trigger ]->Fill( jetEta, jetPhi );
	m_mTriggerEtaPt [ trigger ]->Fill( jetEta, jetPt ); 
	// eta cut

	// convert positive eta to negative because
	// in pp it doesnt matter.
	// our histos run in negative eta
	// (due to pPb configuration)
	// the labels will be taken care of so it is ok
	double jetEtaAdj = AdjustEtaForPP( jetEta );
	
	// saves computation time
	// can theoreticall leave out
	// will go to overflow bins in histos
	if( jetEtaAdj > -constants::FETAMIN ) continue;

	m_mTriggerEtaSpect[ trigger ]->Fill( jetEtaAdj, jetPt );
	// fill mb efficiency histo also
	if( !trigger.compare(mbTrigger) ){
	  m_mTriggerEtaEff[mbTrigger]->Fill( jetEtaAdj, jetPt );
	}
      } // end loop over triggers
    } // end loop over jets
    // EFFICIENCIES
    ProcessEfficiencies( vTrig_jets, vR_jets, m_mTriggerFired );
  } // end event loop
  std::cout << "DONE!" << std::endl;
  m_fIn->Close();
}

void DiJetAnalysisData::
ProcessEfficiencies( std::vector< TLorentzVector >& vTrig_jets,
		     std::vector< TLorentzVector >& vR_jets,
		     std::map< std::string, bool >& mTriggerFired){
  // only do if had min bias trigger
  bool haveMinBiasTrigger = mTriggerFired[ mbTrigger ] ? true : false;
  
  
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
	m_mTriggerEtaEff[ ptTrig.second ]->
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
    m_mTriggerEtaSpect [ trigger ] =
      (TH2D*)m_fIn->Get( Form("h_etaSpect_%s" , trigger.c_str() ) );
    m_mTriggerEtaSpect [ trigger ]->SetDirectory(0);
    m_mTriggerEtaEff   [ trigger ] =
      (TH2D*)m_fIn->Get( Form("h_etaEff_%s"   , trigger.c_str() ) );
    m_mTriggerEtaEff   [ trigger ]->SetDirectory(0);
    m_mTriggerEtaPhi   [ trigger ] =
      (TH2D*)m_fIn->Get( Form("h_etaPhi_%s"   , trigger.c_str() ) );
    m_mTriggerEtaPhi   [ trigger ]->SetDirectory(0);
    m_mTriggerEtaPt    [ trigger ] =
      (TH2D*)m_fIn->Get( Form("h_etaPt_%s"    , trigger.c_str() ) );
    m_mTriggerEtaPt    [ trigger ]->SetDirectory(0);
  }
  
  m_fIn->Close();
}

void DiJetAnalysisData::
PlotSpectra( std::map< std::string, TH2* >& mTriggerH ){
  // if there are none, return
  if( !mTriggerH[ mbTrigger ] ){ return; }  
  
  std::vector< TH1* > vTriggerSpect;

  for( int xBin = 1; xBin <= mTriggerH[ mbTrigger ]->GetNbinsX(); xBin++ ){
  
    TCanvas c_spect("c_spect","c_spect",800,600);
    c_spect.SetLogy();

    double x0, y0, x1, y1;

    if( m_is_pPb ){ x0 = 0.45; y0 = 0.54; x1 = 0.76; y1 = 0.67; }
    else          { x0 = 0.56; y0 = 0.54; x1 = 0.86; y1 = 0.67; }
    
    TLegend l_spect( x0, y0, x1, y1 );
    StyleTools::SetLegendStyle( &l_spect, StyleTools::lSS );
    l_spect.SetFillStyle(0);

    int style = 0;
    double max = -1;

    // should all be the same
    double etaMin = mTriggerH[ mbTrigger ]->
      GetXaxis()->GetBinLowEdge( xBin );
    double etaMax = mTriggerH[ mbTrigger ]->
      GetXaxis()->GetBinUpEdge( xBin );

  
    for( auto& tH : mTriggerH ){
      TH1* hs =
	tH.second->
	ProjectionY( Form("h_spect_%2.0f_Eta_%2.0f_%s",
			  10*std::abs(etaMin),
			  10*std::abs(etaMax),
			  tH.first.c_str() ),
		     xBin, xBin );
      vTriggerSpect.push_back( hs ); 
    
      l_spect.AddEntry(  hs, tH.first.c_str() );
      StyleTools::SetHStyle( hs, style++, StyleTools::hSS);
      hs->Draw("epsame");
      if( max < hs->GetMaximum() ){ max = hs->GetMaximum(); }
    }

    double power = log10(max);
    power = std::ceil(power);
    max = pow( 10, power );
    for( auto& h : vTriggerSpect ){
      h->SetMaximum( max );
      h->SetMinimum( 1 );
    }
  
    l_spect.Draw();
    DrawTools::DrawAtlasInternalDataRight( 0, 0, StyleTools::lSS, m_is_pPb ); 
    DrawTools::DrawLeftLatex( 0.75, 0.725,
			      GetEtaLabel( etaMin, etaMax).c_str(),
			      StyleTools::lSS, 1 );
  
    c_spect.SaveAs( Form("%s/spectra_%2.0f_Eta_%2.0f%s.pdf",
			 m_dirOut.c_str(),
			 std::abs(etaMin)*10,
			 std::abs(etaMax)*10,
			 m_labelOut.c_str() ) );
    c_spect.SaveAs( Form("%s/spectra_%2.0f_Eta_%2.0f%s.png",
			 m_dirOut.c_str(),
			 std::abs(etaMin)*10,
			 std::abs(etaMax)*10,
			 m_labelOut.c_str() ) );
    c_spect.Write( Form("c_spectra_%2.0f_Eta_%2.0f%s",
			std::abs(etaMin)*10,
			std::abs(etaMax)*10,
			m_labelOut.c_str()) );
  } // end loop over eta bins
}

void DiJetAnalysisData::
PlotEfficiencies( std::map< std::string, TH2* >& mTriggerH ){
  // if there are none, return
  if( !mTriggerH[ mbTrigger ] ){ return; }  

  std::vector< TH1* > vMbTrigger;
  std::vector< TH1* > vTriggerEff;
  std::vector< TGraphAsymmErrors* > vTriggerEffGrf;

  for( int xBin = 1; xBin <= mTriggerH[ mbTrigger ]->GetNbinsX(); xBin++ ){
  
    TCanvas c_eff("c_eff","c_eff",800,600);

    TLegend l_eff(0.38, 0.28, 0.61, 0.41);
    StyleTools::SetLegendStyle( &l_eff, StyleTools::lSS );
    l_eff.SetFillStyle(0);

    int style = 0; bool haveDrawn = false;
    double xMin = 0; double xMax = 0;

    // should all be the same
    double etaMin = mTriggerH[ mbTrigger ]->
      GetXaxis()->GetBinLowEdge( xBin );
    double etaMax = mTriggerH[ mbTrigger ]->
      GetXaxis()->GetBinUpEdge( xBin );

    std::map< std::string, TH1* > mTriggerEff;
    std::map< std::string, TGraphAsymmErrors* > mTriggerEffGrf;

    TH1* h_mbTrig =
      mTriggerH[ mbTrigger ]->
      ProjectionY( Form("h_eff_%2.0f_Eta_%2.0f_%s",
			10*std::abs(etaMin),
			10*std::abs(etaMax),
			mbTrigger.c_str() ),
		   xBin, xBin );
    vMbTrigger.push_back( h_mbTrig );
  
    for( auto& trigger : m_vTriggers ){
      mTriggerEff[ trigger ] =
	mTriggerH[ trigger ]->
	ProjectionY( Form("h_eff_%2.0f_Eta_%2.0f_%s",
			  10*std::abs(etaMin),
			  10*std::abs(etaMax),
			  trigger.c_str() ),
		     xBin, xBin);
      vTriggerEff.push_back( mTriggerEff[ trigger ] );
    }
  
    for( auto& trigger : m_vTriggers ){
      mTriggerEffGrf[ trigger ] = new TGraphAsymmErrors();
      TGraphAsymmErrors* tG = mTriggerEffGrf[ trigger ];
      tG->SetName ( Form("gr_eff_%s", trigger.c_str() ) );
      tG->SetTitle( ";#it{p}_{T} [GeV];#it{#varepsilon}_{Trigger}" );
      vTriggerEffGrf.push_back( mTriggerEffGrf[ trigger ] );
      // dont draw the MB efficiencies, its just 1
      // dont add mb efficiencies to legend
      if( !trigger.compare(mbTrigger) ) continue;
      l_eff.AddEntry( tG, trigger.c_str() );

      TH1*            denom = h_mbTrig;
      std::string denomName = mbTrigger;
    
      //-------------------- pp ---------------------
      // in pp we have heavily prescaled triggers
      // dont divide by MB. Use lower trigger.
      // i.e. for j25, use j15.
      // only for j10 use mn
      if( !m_is_pPb ){
	for( auto it = v_tJetPt.begin(); it!= v_tJetPt.end(); it++  ){
	  if( !m_tJetPtTrigger[*it].compare( trigger ) &&
	      it != v_tJetPt.begin() ){
	    denom     = mTriggerEff[ m_tJetPtTrigger[ *(it-1) ] ];
	    denomName = m_tJetPtTrigger[ *(it-1) ];
	    break;
	  }
	}
      }
      //---------------------------------------------

      std::cout << "-----For: "     << trigger
		<< " dividing by: " << denomName << std::endl;

    
      tG->Divide( mTriggerEff[ trigger ], denom,
		  "cl=0.683 b(1,1) mode" );
      tG->SetMaximum( 1.3 );
      tG->SetMinimum( 0. );
      StyleTools::SetHStyle( tG, style++, StyleTools::hSS);

      if( !haveDrawn ) {
	tG->Draw("ap"); // first one
	haveDrawn = true;
	xMin = tG->GetXaxis()->GetXmin();
	xMax = tG->GetXaxis()->GetXmax();
      } else { tG->Draw("p same"); } // others
    }
  
    l_eff.Draw();

    TLine line( xMin, 1, xMax, 1);
    line.Draw();

    DrawTools::DrawAtlasInternalDataRight( 0, 0,  StyleTools::lSS, m_is_pPb ); 
    DrawTools::DrawLeftLatex( 0.5, 0.8,
			      GetEtaLabel( etaMin, etaMax).c_str(),
			      StyleTools::lSS, 1 );

    c_eff.SaveAs( Form("%s/efficiencies_%2.0f_Eta_%2.0f%s.pdf",
		       m_dirOut.c_str(),
		       std::abs(etaMin)*10,
		       std::abs(etaMax)*10,
		       m_labelOut.c_str() ) );
    c_eff.SaveAs( Form("%s/efficiencies_%2.0f_Eta_%2.0f%s.png",
		       m_dirOut.c_str(),
		       std::abs(etaMin)*10,
		       std::abs(etaMax)*10,
		       m_labelOut.c_str() ) );
    c_eff.Write( Form("c_efficiencies_%2.0f_Eta_%2.0f%s",
		      std::abs(etaMin)*10,
		      std::abs(etaMax)*10,
		      m_labelOut.c_str() ) );
  } // end loop over eta
}

void DiJetAnalysisData::
PlotEtaPhiPtMap( std::map< std::string, TH2* >& mTriggerH ){
  TCanvas c_map("c_map","c_map",800,600);

  for( auto& tH : mTriggerH ){
    tH.second->Draw("col");
    StyleTools::SetHStyle( tH.second, 0, StyleTools::hSS);
    DrawTools::DrawAtlasInternalDataRight( 0, -0.55,
					   StyleTools::lSS, m_is_pPb );  
    c_map.SaveAs( Form("%s/%s%s.pdf", 
		       m_dirOut.c_str(),
		       tH.second->GetName(),
		       m_labelOut.c_str() ) );
    c_map.SaveAs( Form("%s/%s%s.png", 
		       m_dirOut.c_str(),
		       tH.second->GetName(),
		       m_labelOut.c_str() ) );
    c_map.Write( Form("c_%s%s",
		      tH.second->GetName(),
		      m_labelOut.c_str() ) );
  }
}
