#include <TFile.h>
#include <TTree.h>
#include <TLine.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TEnv.h>
#include <TCanvas.h>
#include <TEfficiency.h>
#include <TGraphAsymmErrors.h>

#include <iostream>

#include "MyRoot.h"

#include "DiJetAnalysisData.h"

DiJetAnalysisData::DiJetAnalysisData()
  : DiJetAnalysis( true, true )
{}

DiJetAnalysisData::DiJetAnalysisData( bool isData, bool is_pPb )
  : DiJetAnalysis( isData, is_pPb)
{}

DiJetAnalysisData::~DiJetAnalysisData(){}

void DiJetAnalysisData::RunOverTreeFillHistos( int nEvents, 
					   int startEvent ){  
  loadTriggers();
  setupHistograms();
  processEvents( nEvents, startEvent );
 }

void DiJetAnalysisData::PlotExistingHistos(){
  loadTriggers();
  loadHistograms();
    
  plotSpectra();
  plotEfficiencies();
  plotEtaPhi();
  plotEtaPt();
}

//---------------------------------
//            Read Data
//---------------------------------

void DiJetAnalysisData::loadTriggers(){
  std::string configName   = "config/configJetsFwd";
  std::string configSuffix = m_is_pPb ? "_pPb.cfg" : "_pp.cfg";
  configName += configSuffix;
  
  TEnv config;
  config.ReadFile( configName.c_str(), EEnvLevel(0));

  std::string triggerMenu = config.GetValue("triggerMenu","");
  v_triggers =
    vectorise( config.GetValue( Form("triggers.%s", triggerMenu.c_str() ),"") , " " );

  for( auto& trigger : v_triggers ){
    if( trigger.find("_mb_") != std::string::npos ){ 
      mbTrigger = trigger; 
      std::cout << "Found Min Bias trigger " << mbTrigger << std::endl;
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
    for( auto& trigger : v_triggers ){
      if( trigger.find(triggerLabel ) != std::string::npos ){
	m_tJetPtTrigger[pt] = trigger;
      }
    }
  }

  for( auto ptTrig = m_tJetPtTrigger.begin(); ptTrig!= m_tJetPtTrigger.end(); ptTrig++ ){
    std::cout << "For " << ptTrig->first << " GeV Trigger Jet : " 
	      << ptTrig->second << std::endl;
  }
}

void DiJetAnalysisData::setupHistograms(){
  // Triggers and Spectra
  int    nPtSpectBins = 100; 
  double ptSpectMin   = 0; double ptSpectMax  = nPtSpectBins;

  int    nPtEffBins   = 100; 
  double ptEffMin     = 0; double ptEffMax  = nPtEffBins;
  
  // Eta-Phi Maps
  int    nEtaBins = 100; 
  double etaMin   = -5; double etaMax  = 5;

  int    nPhiBins = 64; 
  double phiMin   = -constants::PI; double phiMax  = constants::PI; 

  int    nPtBins  = 100;
  double ptMin    = 10;
  double ptMax    = ptMin + static_cast<double>(nPtBins)/2;

  int    nEtaForwardBinsFine   = 12;
  int    nEtaForwardBinsCoarse = 3;
  double etaForwardMin   = -constants::FETAMAX;
  double etaForwardMax   = -constants::FETAMIN;

  int style = 0; double scale = 0.6;
  
  for( auto& trigger : v_triggers ){
    std::cout << "Making - " << trigger << " histograms " << std::endl;
    m_triggerSpect[ trigger ] = 
      new TH1D( Form("h_spect_%s", trigger.c_str() ), 
		Form("h_%s;#it{p}_{T} [GeV];d#it{N}/d#it{p}_{T}",
		     trigger.c_str() ),
		nPtSpectBins, ptSpectMin, ptSpectMax ) ;
    m_triggerSpect[ trigger ]->Sumw2();
    v_hists.push_back( m_triggerSpect[ trigger ] );
    SetHStyle( m_triggerSpect[ trigger ], style, scale );

   m_triggerSpectEta[ trigger ] = 
      new TH2D( Form("h_spectEta_%s", trigger.c_str() ), 
		Form("h_%s;#eta;#it{p}_{T} [GeV]",
		     trigger.c_str() ),
		nEtaForwardBinsCoarse, etaForwardMin, etaForwardMax,
		nPtSpectBins, ptSpectMin, ptSpectMax ) ;
    m_triggerSpectEta[ trigger ]->Sumw2();
    v_hists.push_back( m_triggerSpectEta[ trigger ] );
    SetHStyle( m_triggerSpectEta[ trigger ], style, scale );
 
    
    m_triggerEff[ trigger ] = 
      new TH1D( Form("h_eff_%s", trigger.c_str() ), 
		Form("h_%s;#it{p}_{T} [GeV];#it{#varepsilon}_{Trigger}",
		     trigger.c_str() ),
		nPtEffBins, ptEffMin, ptEffMax ) ;
    m_triggerEff[ trigger ]->Sumw2();
    v_hists.push_back( m_triggerEff[ trigger ] );
    SetHStyle( m_triggerEff[ trigger ], style, scale );

    m_triggerEffEta[ trigger ] = 
      new TH2D( Form("h_effEta_%s", trigger.c_str() ), 
		Form("h_%s;#eta;#it{p}_{T} [GeV]",
		     trigger.c_str() ),
		nEtaForwardBinsCoarse, etaForwardMin, etaForwardMax,
		nPtEffBins, ptEffMin, ptEffMax ) ;
    m_triggerEffEta[ trigger ]->Sumw2();
    v_hists.push_back( m_triggerEffEta[ trigger ] );
    SetHStyle( m_triggerEffEta[ trigger ], 0, scale );
    
    m_triggerRunPrescale[ trigger ] = 
      new TH2D( Form("h_runPrescale_%s", trigger.c_str() ), 
		Form("h_%s;Run No;Prescale", trigger.c_str() ),
		1700, 312600, 314300, 1999, 1, 2000) ;
    m_triggerRunPrescale[ trigger ]->Sumw2();
    m_triggerRunPrescale[ trigger ]->GetXaxis()->SetNdivisions(505);  
    m_triggerRunPrescale[ trigger ]->GetYaxis()->SetNdivisions(505);  
    v_hists.push_back( m_triggerRunPrescale[ trigger ] );
    
    m_triggerEtaPhi[ trigger ] = 
      new TH2D( Form("h_etaPhi_%s", trigger.c_str() ), 
		Form("h_%s;#eta;#phi", trigger.c_str() ),
		nEtaBins, etaMin, etaMax,
		nPhiBins, phiMin, phiMax ) ;
    m_triggerEtaPhi[ trigger ]->Sumw2();
    m_triggerEtaPhi[ trigger ]->GetXaxis()->SetNdivisions(505);  
    m_triggerEtaPhi[ trigger ]->GetYaxis()->SetNdivisions(505);  
    v_hists.push_back( m_triggerEtaPhi[ trigger ] );
    SetHStyle( m_triggerEtaPhi[ trigger ], 0, scale );
  
    m_triggerEtaPt[ trigger ] = 
      new TH2D( Form("h_etaPt_%s", trigger.c_str() ), 
		Form("h_%s;#eta;#it{p}_{T}", trigger.c_str() ),
		nEtaForwardBinsFine, etaForwardMin, etaForwardMax,
		nPtBins, ptMin, ptMax) ;
    m_triggerEtaPt[ trigger ]->Sumw2();
    m_triggerEtaPt[ trigger ]->GetXaxis()->SetNdivisions(505);  
    m_triggerEtaPt[ trigger ]->GetYaxis()->SetNdivisions(505);  
    v_hists.push_back( m_triggerEtaPt[ trigger ] );
    SetHStyle( m_triggerEtaPt[ trigger ], 0, scale );

    style++;
  }
}

void DiJetAnalysisData::processEvents( int nEvents, int startEvent ){
  std::string fNameIn = m_is_pPb ?
    "/home/yakov/Projects/atlas/data/pPb.root" :
    "/home/yakov/Projects/atlas/data/pp.root"  ;

  std::cout << "fNameIn: " << fNameIn << std::endl;
  
  // collections and variables
  std::vector< TLorentzVector >    vTrig_jets;
  std::vector< TLorentzVector >* p_vTrig_jets = &vTrig_jets;

  std::vector< TLorentzVector >    vR_jets;
  std::vector< TLorentzVector >* p_vR_jets = &vR_jets;

  std::vector< bool >    v_isCleanJet;
  std::vector< bool >* p_v_isCleanJet = &v_isCleanJet;

  std::map< std::string, bool  > m_triggerFired;
  std::map< std::string, float > m_triggerPrescale;

  int    runNumber     = 0;

  //----------------------------------------
  //  Open file and tree, Fill histograms
  //----------------------------------------
  TFile* fin = TFile::Open( fNameIn.c_str() );
  TTree* tree = (TTree*) fin->Get( "tree" );

  // Connect to tree
  tree->SetBranchAddress( "vTrig_jets"  , &p_vTrig_jets );
  tree->SetBranchAddress( "vR_C_jets"   , &p_vR_jets );
  tree->SetBranchAddress( "v_isCleanJet", &p_v_isCleanJet );

  for( auto& trigger : v_triggers ){
    tree->SetBranchAddress( Form("passed_%s",trigger.c_str() ), &m_triggerFired[trigger] );
    tree->SetBranchAddress( Form("prescale_%s", trigger.c_str() ), &m_triggerPrescale[trigger] );
  }

  tree->SetBranchAddress( "runNumber"   , &runNumber );  

  // n events
  int maxEvents = tree->GetEntries();

  nEvents = nEvents > 0 ? nEvents : maxEvents;
  startEvent = startEvent < maxEvents ? startEvent : maxEvents - 1;
  
  int endEvent = startEvent + nEvents < maxEvents ? startEvent + nEvents : maxEvents;

  // event loop
  for( int ev = startEvent; ev < endEvent; ev++ ){
    tree->GetEntry( ev );

    applyCleaning ( vR_jets, v_isCleanJet );
    applyIsolation( 1.0, vR_jets );

    if( DoPrint(ev) ) {
      std::cout << "\nEvent : " << ev << "    runN : " << runNumber
		<< "    has : " << vR_jets.size() << " jets" 
		<< std::endl; 
    }
    
    // loop over jets 
    for( auto& jet : vR_jets ){
      double pt  = jet.Pt();
      if( pt == 0 ) continue;
      // SPECTRA AND ETA PHI PT MAPS
      // loop over passed triggers
      for( auto& trigger : v_triggers ){
	// check if we have that trigger
	if( !m_triggerFired[trigger] ) continue;
	// ETA-PHI
	double jetEta = jet.Eta();
	double jetPhi = jet.Phi();
	m_triggerEtaPhi[ trigger ]->Fill( jetEta, jetPhi );
	// eta cut
	if( jet.Eta() > -constants::FETAMIN ) continue;
	double jetPt = jet.Pt()/1000.;
	m_triggerSpect   [ trigger ]->Fill( jetPt );
	m_triggerSpectEta[ trigger ]->Fill( jetEta, jetPt );
	m_triggerEtaPt   [ trigger ]->Fill( jetEta, jetPt ); 
	// fill mb efficiency histo also
	if( !trigger.compare(mbTrigger) ){
	  m_triggerEff   [mbTrigger]->Fill( jetPt );
	  m_triggerEffEta[mbTrigger]->Fill( jetEta, jetPt );
	}
      } // end loop over triggers
    } // end loop over jets
    // EFFICIENCIES
    processEfficiencies( vTrig_jets, vR_jets, m_triggerFired );
  } // end event loop
  std::cout << "DONE!" << std::endl;

  //----------------------------------------
  //  Close the input file, 
  //  write histos to output
  //----------------------------------------
  fin->Close();
  std::cout << "fNameOut: " << m_fNameOut << std::endl;
  TFile* fout = new TFile( m_fNameOut.c_str(),"RECREATE");
  for( auto& h  : v_hists  ) { h-> Write(); }
  for( auto& f  : v_functs ) { f-> Write(); }
  for( auto& gr : v_graphs ) { gr->Write(); }

  
  fout->Close();
}

void DiJetAnalysisData::processEfficiencies( std::vector< TLorentzVector >& vTrig_jets,
			  std::vector< TLorentzVector >& vR_jets,
			  std::map< std::string, bool >& m_triggerFired){
  // only do if had min bias trigger
  bool haveMinBiasTrigger = m_triggerFired[ mbTrigger ] ? true : false;
  
  
  if( !haveMinBiasTrigger ) return;
  if( vTrig_jets.empty()  ) return;

  // sort trigger jets
  std::sort( vTrig_jets.begin(), vTrig_jets.end(), sortByDecendingPt );

  // take highest pt trigger jet in forward eta range
  TLorentzVector* tJet = NULL;
  for( auto& jet : vTrig_jets ){
    if( jet.Eta() > -constants::FETAMIN ) continue;
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
	if( jet.Eta() > -constants::FETAMIN ) continue;
	m_triggerEff   [ ptTrig.second ]->
	  Fill( jet.Pt()/1000. );
	m_triggerEffEta[ ptTrig.second ]->
	  Fill( jet.Eta(), jet.Pt()/1000. );
      } // end loop over jets
    }
  } // end loop over triggers 
}

//---------------------------------
//            Plot Data
//---------------------------------

void DiJetAnalysisData::loadHistograms(){
  TFile* fin = TFile::Open( m_fNameOut.c_str() ); 

  for( auto& trigger : v_triggers ){
    m_triggerSpect [ trigger ] = (TH1D*)fin->Get( Form("h_spect_%s" , trigger.c_str() ) );
    m_triggerEff   [ trigger ] = (TH1D*)fin->Get( Form("h_eff_%s"   , trigger.c_str() ) );
    m_triggerEtaPhi[ trigger ] = (TH1D*)fin->Get( Form("h_etaPhi_%s", trigger.c_str() ) );
    m_triggerEtaPt [ trigger ] = (TH1D*)fin->Get( Form("h_etaPt_%s" , trigger.c_str() ) );
  }
}

void DiJetAnalysisData::plotSpectra(){
  TCanvas* c_spect = new TCanvas("c_spect","c_spect",800,600);
  c_spect->SetLogy();

  TLegend* l_spect = new TLegend(0.23, 0.23, 0.54, 0.36);
  SetLegendStyle( l_spect, 0.55 );
  l_spect->SetFillStyle(0);

  int style = 0;
  
  double max = -1;

  for( auto& sh : m_triggerSpect ){
    l_spect->AddEntry( sh.second, sh.first.c_str() );
    sh.second->Draw("epsame");
    SetHStyle( sh.second, style++, 0.6);
    if( max < sh.second->GetMaximum() ){ max = sh.second->GetMaximum(); }
  }

  double power = log10(max);
  power = std::ceil(power);
  max = pow( 10, power );
  for( auto& sh : m_triggerSpect ){
    sh.second->SetMaximum( max );
    sh.second->SetMinimum( 1 );
  }

  l_spect->Draw();

  DrawAtlasInternalDataRight( 0, 0, 0.6, m_is_pPb ); 

  c_spect->SaveAs( Form("output/spectra%s.pdf", m_labelOut.c_str() ) );

}

void DiJetAnalysisData::plotEfficiencies(){
  TCanvas* c_eff = new TCanvas("c_eff","c_eff",800,600);

  TLegend* l_eff = new TLegend(0.38, 0.28, 0.61, 0.41);
  SetLegendStyle( l_eff, 0.55 );
  l_eff->SetFillStyle(0);

  int style = 0; bool haveDrawn = false;
  double xMin = 0; double xMax = 0;
  
  TH1* h_mbTrig = m_triggerEff[ mbTrigger ];

  std::map< std::string, TGraphAsymmErrors* > m_triggerEffGrf;
  std::map< std::string, TGraphAsymmErrors* > m_triggerEffGrfEta;
  
  for( auto& trigger : v_triggers ){
    m_triggerEffGrf[ trigger ] = new TGraphAsymmErrors();
    TGraphAsymmErrors* tG = m_triggerEffGrf[ trigger ];
    tG->SetName ( Form("gr_eff_%s", trigger.c_str() ) );
    tG->SetTitle( Form("gr_%s;#it{p}_{T} [GeV];#it{#varepsilon}_{Trigger}",
			     trigger.c_str() ) );
    // dont draw the MB efficiencies, its just 1
    // dont add mb efficiencies to legend
    if( !trigger.compare(mbTrigger) ) continue;
    l_eff->AddEntry( tG, trigger.c_str() );

    TH1*            denom = h_mbTrig;
    std::string denomName = "";
    
    //-------------------- pp ---------------------
    // in pp we have heavily prescaled triggers
    // dont divide by MB. Use lower trigger.
    // i.e. for j25, use j15.
    // only for j10 use mn
    if( !m_is_pPb ){
      for( auto it = v_tJetPt.begin(); it!= v_tJetPt.end(); it++  ){
	if( !m_tJetPtTrigger[*it].compare( trigger ) &&
	    it != v_tJetPt.begin() ){
	  denom     = m_triggerEff[ m_tJetPtTrigger[ *(it-1) ] ];
	  denomName = m_tJetPtTrigger[ *(it-1) ];
	  break;
	}
      }
    }
    //---------------------------------------------

    std::cout << "-----For: "     << trigger
	      << " dividing by: " << denomName << std::endl;

    
    tG->Divide( m_triggerEff[ trigger ], denom,
		"cl=0.683 b(1,1) mode" );
    tG->SetMaximum( 1.3 );
    tG->SetMinimum( 0. );
    SetHStyle( tG, style++, 0.6);

    if( !haveDrawn ) {
      tG->Draw("ap"); // first one
      haveDrawn = true;
      xMin = tG->GetXaxis()->GetXmin();
      xMax = tG->GetXaxis()->GetXmax();
    } else { tG->Draw("p same"); } // others
  }
  
  l_eff->Draw();

  TLine* line = new TLine( xMin, 1, xMax, 1);
  line->Draw();

  DrawAtlasInternalDataRight( 0, 0,  0.6, m_is_pPb ); 

  c_eff->SaveAs( Form("output/efficiencies%s.pdf", m_labelOut.c_str() ) );
}

void DiJetAnalysisData::plotEtaPhi(){
  TCanvas* c_etaPhi = new TCanvas("c_etaPhi","c_etaPhi",800,600);

  for( auto& sh : m_triggerEtaPhi ){
    sh.second->Draw("col");
    SetHStyle( sh.second, 0, 0.6);
    DrawAtlasInternalDataRight( 0, -0.55, 0.6, m_is_pPb );  
    c_etaPhi->SaveAs( Form("output/etaPhi%s_%s.pdf", 
			   m_labelOut.c_str(), sh.first.c_str() ) );
    std::cout << sh.first << std::endl;
  }
}

void DiJetAnalysisData::plotEtaPt(){
  TCanvas* c_ptEta = new TCanvas("c_ptEta","c_ptEta",800,600);

  for( auto& sh : m_triggerEtaPt ){
    sh.second->Draw("col");
    SetHStyle( sh.second, 0, 0.6);
    DrawAtlasInternalDataLeft( 0, 0, 0.6, m_is_pPb );  
    c_ptEta->SaveAs( Form("output/ptEta%s_%s.pdf", 
			  m_labelOut.c_str(), sh.first.c_str() ) );
  }
}
