#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TROOT.h>
#include <TEnv.h>
#include <TCanvas.h>
#include <TRandom.h>

#include <iostream>

#include "MyRoot.h"

#include "DiJetAnalysisMC.h"
#include "JetPair.h"

DiJetAnalysisMC::DiJetAnalysisMC()
{}

DiJetAnalysisMC::DiJetAnalysisMC( bool isData, bool is_pPb )
  : DiJetAnalysis( isData, is_pPb)
{}

DiJetAnalysisMC::~DiJetAnalysisMC(){}

const double lStyleScale = 0.60;
const double hStyleScale = 0.75;

void DiJetAnalysisMC::Initialize(){
  DiJetAnalysis::Initialize();

  m_dRmax = 0.2;
  
  m_v_usedJZN.push_back( 1 );
  m_m_jznFnameIn[ m_v_usedJZN.back() ] = "/home/yakov/Projects/atlas/data/jz1.root";
  m_m_jznSigma  [ m_v_usedJZN.back() ] = 6.7890E+07;
  m_m_jznEff    [ m_v_usedJZN.back() ] = 2.8289E-03;

  m_v_usedJZN.push_back( 2 );
  m_m_jznFnameIn[ m_v_usedJZN.back() ] = "/home/yakov/Projects/atlas/data/jz2.root";
  m_m_jznSigma  [ m_v_usedJZN.back() ] = 6.3997E+05;
  m_m_jznEff    [ m_v_usedJZN.back() ] = 4.2714E-03;

  // calculate sum of sigma and eff
  m_sumSigmaEff = 0;
  for( auto jzn : m_v_usedJZN ){
    m_sumSigmaEff += ( m_m_jznSigma[ jzn ] * m_m_jznEff[ jzn ] );
  }

  std::cout << "SumSigmaEff = " << m_sumSigmaEff << std::endl;
  
  for( auto jzn : m_v_usedJZN ){ // loop over JZ samples
    m_fIn = TFile::Open( m_m_jznFnameIn[ jzn ].c_str() );
    m_tree = (TTree*) m_fIn->Get( "tree" );

    int nEventsTotal = m_tree->GetEntries();
    double     sigma = m_m_jznSigma[ jzn ];      
    double       eff = m_m_jznEff  [ jzn ];
  
    m_m_jznNev    [ jzn ] = nEventsTotal;
    m_m_jznWeights[ jzn ] =
      (1./nEventsTotal)*(sigma * eff)*(1./m_sumSigmaEff);
  
    std::cout << "JZ" << jzn << "   weight = " << m_m_jznWeights [jzn] << std::endl;
    m_fIn->Close();
  }
}

//---------------------------------
//            Read Data
//---------------------------------
void DiJetAnalysisMC::RunOverTreeFillHistos( int nEvents, 
					     int startEvent ){  
  SetupHistograms();
  ProcessEvents( nEvents, startEvent );
  SaveOutputs();
}

//---------------------------------
//            Plot Data
//---------------------------------
void DiJetAnalysisMC::ProcessPlotHistos(){
  LoadHistograms();

  std::string cfNameOut = m_dirOut + "/c_myOut" + m_labelOut + ".root";
  m_fOut = new TFile( cfNameOut.c_str(),"RECREATE");
  
  int etaBinMax;
  
  // check if we have one
  if( !m_m_jznEtaSpect.size() ){ return; }
  // make some spectra in various eta bins
  etaBinMax = m_m_jznEtaSpect.begin()->second->GetNbinsX();
  for( int etaBin = 1; etaBin <= etaBinMax; etaBin++){
    PlotSpectra( m_m_jznEtaSpect, etaBin, etaBin );
  }
  PlotSpectra( m_m_jznEtaSpect, 1, etaBinMax );
  
  PlotEtaPhiPtMap( m_m_jznEtaPhi );
  PlotEtaPhiPtMap( m_m_jznEtaPt  );

  // check if we have one
  if( !m_m_jznRpt.size() ){ return; }
  // make some spectra in various eta bins
  etaBinMax = m_m_jznRpt.begin()->second->GetNbinsX();
  for( int etaBin = 1; etaBin <= etaBinMax; etaBin++){
    PlotVsEtaPt( etaBin, etaBin,
		 m_m_jznRpt, m_m_jznNentries, 0);
  }
  PlotVsEtaPt( 1, etaBinMax,
	       m_m_jznRpt, m_m_jznNentries, 0);

  
  // check if we have one
  if( !m_m_jznDeta.size() ){ return; }
  // make some spectra in various eta bins
  etaBinMax = m_m_jznDeta.begin()->second->GetNbinsX();
  for( int etaBin = 1; etaBin <= etaBinMax; etaBin++){
    PlotVsEtaPt( etaBin, etaBin,
		 m_m_jznDeta, m_m_jznNentries, 1);
  }
  PlotVsEtaPt( 1, etaBinMax,
	       m_m_jznDeta, m_m_jznNentries, 1);
  
  // check if we have one
  if( !m_m_jznDphi.size() ){ return; }
  // make some spectra in various eta bins
  etaBinMax = m_m_jznDphi.begin()->second->GetNbinsX();
  for( int etaBin = 1; etaBin <= etaBinMax; etaBin++){
    PlotVsEtaPt( etaBin, etaBin,
		 m_m_jznDphi, m_m_jznNentries, 2);
  }
  PlotVsEtaPt( 1, etaBinMax,
	       m_m_jznDphi, m_m_jznNentries, 2);

  std::cout << "DONE! Closing " << cfNameOut << std::endl;
  m_fOut->Close();
  std::cout << "......Closed  " << cfNameOut << std::endl;
}

//---------------------------------
//            Read Data
//---------------------------------

void DiJetAnalysisMC::SetupHistograms(){
  // Triggers and Spectra
  int    nPtSpectBins = 100; 
  double ptSpectMin   = 0; double ptSpectMax  = nPtSpectBins;
  
  // Eta-Phi Maps
  int    nEtaBins = 100; 
  double etaMin   = constants::ETAMIN; double etaMax  = constants::ETAMAX;
  
  int    nPhiBins = 64; 
  double phiMin   = -constants::PI; double phiMax  = constants::PI; 
  
  double ptWidth  = 2;
  double ptMin    = 10;
  double ptMax    = 100;
  int    nPtBins  = (ptMax - ptMin)/ptWidth;
  
  // JES JER etc
  int    nEtaForwardBinsFine   = 12;
  int    nEtaForwardBinsCoarse = 3;
  double etaForwardMin   = -constants::FETAMAX;
  double etaForwardMax   = -constants::FETAMIN;

  double ptTruthWidth  = 5;
  double ptTruthMin    = 10;
  double ptTruthMax    = 100;
  int    nPtTruthBins  = (ptTruthMax - ptTruthMin) / ptTruthWidth;

  int    nRPtBins   = 60;
  double rPtMin     = 0; double rPtMax = 3;

  int    nDAngleBins  = 50;
  double dAngleMin    = -0.5; double dAngleMax = 0.5;
  
  for( auto jzn : m_v_usedJZN ){ // loop over JZ samples
    std::cout << "Making - JZ" << jzn << " histograms " << std::endl;

    m_m_jznEtaSpect[ jzn ] = 
      new TH2D( Form("h_etaSpect_jz%i", jzn ), 
		";#eta_{Reco};#it{p}_{T}^{Reco} [GeV];dN/d#it{p}_{T}",
		nEtaForwardBinsCoarse, etaForwardMin, etaForwardMax,
		nPtSpectBins, ptSpectMin, ptSpectMax ) ;
    AddHistogram( m_m_jznEtaSpect[ jzn ] );
    
    m_m_jznEtaPhi[ jzn ] =
      new TH2D( Form("h_etaPhi_jz%i", jzn),
		";#eta_{Reco};#phi_{Reco}",
		nEtaBins, etaMin, etaMax,
		nPhiBins, phiMin, phiMax );
    AddHistogram( m_m_jznEtaPhi[ jzn ] );
      
    m_m_jznEtaPt[ jzn ] =
      new TH2D( Form("h_etaPt_jz%i", jzn),
		";#eta_{Reco};#it{p}_{T}^{Reco}",
		nEtaBins, etaMin, etaMax,
		nPtBins , ptMin  , ptMax );
    AddHistogram( m_m_jznEtaPt[ jzn ] );
    
    m_m_jznRpt[ jzn ] =
      new TH3D( Form("h_rPt_jz%i", jzn),
		";#eta^{Truth};#it{p}_{T}^{Truth};#it{p}_{T}^{Reco}/#it{p}_{T}^{Truth}",
		nEtaForwardBinsFine, etaForwardMin, etaForwardMax,
		nPtTruthBins, ptTruthMin, ptTruthMax,
		nRPtBins, rPtMin, rPtMax);
    AddHistogram( m_m_jznRpt[ jzn ] );
  
    m_m_jznDeta[ jzn ] =
      new TH3D( Form("h_dEta_jz%i", jzn),
		";#eta^{Truth};#it{p}_{T}^{Truth};#eta^{Reco}-#eta^{Truth}",
		nEtaForwardBinsFine, etaForwardMin, etaForwardMax,
		nPtTruthBins, ptTruthMin, ptTruthMax,
		nDAngleBins, dAngleMin, dAngleMax );
    AddHistogram( m_m_jznDeta[ jzn ] );    
  
    m_m_jznDphi[ jzn ] =
      new TH3D( Form("h_dPhi_jz%i", jzn),
		";#eta^{Truth};#it{p}_{T}^{Truth};#phi^{Reco}-#phi^{Truth}",
		nEtaForwardBinsFine, etaForwardMin, etaForwardMax,
		nPtTruthBins, ptTruthMin, ptTruthMax,
		nDAngleBins, dAngleMin, dAngleMax );
    AddHistogram( m_m_jznDphi[ jzn ] );

    m_m_jznNentries[ jzn ] =
      new TH2D( Form("h_nEntries_jz%i", jzn),
		";#eta^{Truth};#it{p}_{T}^{Truth}",
		nEtaForwardBinsFine, etaForwardMin, etaForwardMax,
		nPtTruthBins, ptTruthMin, ptTruthMax );
    AddHistogram( m_m_jznNentries[ jzn ] );
  } 
}

void DiJetAnalysisMC::ProcessEvents( int nEvents, int startEvent ){
  // collections and variables
  std::vector< TLorentzVector >    vT_jets;
  std::vector< TLorentzVector >* p_vT_jets = &vT_jets;

  std::vector< TLorentzVector >    vR_jets;
  std::vector< TLorentzVector >* p_vR_jets = &vR_jets;

  std::vector< bool >    v_isCleanJet;
  std::vector< bool >* p_v_isCleanJet = &v_isCleanJet;

  for( auto jzn : m_v_usedJZN ){
    
    std::cout << "fNameIn: " << m_m_jznFnameIn[ jzn ] << std::endl;
  
    m_fIn = TFile::Open( m_m_jznFnameIn[ jzn ].c_str() );
    m_tree = (TTree*) m_fIn->Get( "tree" );

    // Connect to tree
    m_tree->SetBranchAddress( "vT_jets"     , &p_vT_jets );
    m_tree->SetBranchAddress( "vR_C_jets"   , &p_vR_jets );
    m_tree->SetBranchAddress( "v_isCleanJet", &p_v_isCleanJet );

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
	std::cout << "\nEvent : " << ev 
		  << "    has : " << vR_jets.size() << " reco jets"
		  << "    and : " << vT_jets.size() << " truth jets" 
		  << std::endl; 
      }

      std::vector< JetPair > v_paired_jets;
      PairJets( vR_jets, vT_jets, v_paired_jets );

	
      for( auto& vp : v_paired_jets ){	  
	double etaTruth = vp.TruthJet()->Eta();
	double phiTruth = vp.TruthJet()->Phi();
	double  ptTruth = vp.TruthJet()->Pt()/1000.;
	  
	double etaReco = vp.RecoJet()->Eta();
	double phiReco = vp.RecoJet()->Phi();
	double  ptReco = vp.RecoJet()->Pt()/1000.;
	  
	  
	m_m_jznEtaSpect[ jzn ]->Fill( etaReco, ptReco );
	m_m_jznEtaPhi  [ jzn ]->Fill( etaReco, phiReco );
	m_m_jznEtaPt   [ jzn ]->Fill( etaReco, ptReco );
	  
	if( vp.DeltaR() <= m_dRmax ){
	  m_m_jznRpt   [ jzn ]->Fill( etaTruth, ptTruth, ptReco/ptTruth );
	  m_m_jznDeta  [ jzn ]->Fill( etaTruth, ptTruth, etaReco - etaTruth );
	  m_m_jznDphi  [ jzn ]->Fill( etaTruth, ptTruth, phiReco - phiTruth );
	  m_m_jznNentries[ jzn ]->Fill( etaTruth, ptTruth );
	}
      } // end loop over pairs
    } // end loop over events
   
    std::cout << "DONE WITH JZ" << std::endl;

    m_fIn->Close();
  } // end loop over a JZ sample
}

// pair reco and truth jets with deltaR parameter 
void DiJetAnalysisMC::PairJets(  std::vector< TLorentzVector >& vR_jets,
				 std::vector< TLorentzVector >& vT_jets,
				 std::vector< JetPair >&  v_paired_jets ){ 
  // clear vectors from previous time
  v_paired_jets.clear();
  
  // exit if either is empty
  if( !vT_jets.size() || !vR_jets.size() ){ return; }
  
  for( auto& truthJet : vT_jets){
    // for each truth jet, need to find closest reco jet
    // set deltaRmin to something large, find jet with smallest
    // deltaRmin less than Rmax and pair that to the truth jet
    double              deltaRmin = 4*constants::ETAMAX;
    TLorentzVector* pairedRecoJet = NULL;
    for( auto& recoJet : vR_jets ){
      // "cleaned" jets had their px, py, pz set to 0
      // to avoid changing the vector size
      if( recoJet.Pt() == 0 ) continue;
      // if we have a "real" reco jet, continue
      if( AnalysisTools::DeltaR( recoJet, truthJet ) <= deltaRmin ) {	
	pairedRecoJet = &recoJet;
	deltaRmin     = AnalysisTools::DeltaR( recoJet, truthJet );
      }
    }  // end loop over reco jets
    // just for safety. if there is at least
    // one reco jet, it should pair to all truth jets
    if( !pairedRecoJet ) continue;
    // we found one!
    // add the truth and reco jets together
    // also store their deltaR
    v_paired_jets.emplace_back( JetPair( pairedRecoJet,
					 &truthJet,
					 deltaRmin) );  
  } // end loop over truth jets
}


//---------------------------------
//            Plot Data
//---------------------------------

void DiJetAnalysisMC::LoadHistograms(){
  m_fIn = TFile::Open( m_rootFname.c_str() ); 

  for( auto& jzn : m_v_usedJZN ){
    m_m_jznEtaSpect [ jzn ] =
      (TH2D*)m_fIn->Get( Form("h_etaSpect_jz%i" , jzn ) );
    m_m_jznEtaSpect [ jzn ]->SetDirectory(0);
    m_m_jznEtaPhi   [ jzn ] =
      (TH2D*)m_fIn->Get( Form("h_etaPhi_jz%i"   , jzn ) );
    m_m_jznEtaPhi   [ jzn ]->SetDirectory(0);
    m_m_jznEtaPt    [ jzn ] =
      (TH2D*)m_fIn->Get( Form("h_etaPt_jz%i"    , jzn ) );
    m_m_jznEtaPt    [ jzn ]->SetDirectory(0);

    m_m_jznRpt     [ jzn ] =
      (TH3D*)m_fIn->Get( Form("h_rPt_jz%i"    , jzn ) );
    m_m_jznRpt     [ jzn ]->SetDirectory(0);
    m_m_jznDeta    [ jzn ]   =
      (TH3D*)m_fIn->Get( Form("h_dEta_jz%i"    , jzn ) );
    m_m_jznDeta    [ jzn ]->SetDirectory(0);
    m_m_jznDphi    [ jzn ] =
      (TH3D*)m_fIn->Get( Form("h_dPhi_jz%i"    , jzn ) );
    m_m_jznDphi    [ jzn ]->SetDirectory(0);

    m_m_jznNentries [ jzn ] =
      (TH2D*)m_fIn->Get( Form("h_nEntries_jz%i"    , jzn ) );
    m_m_jznNentries [ jzn ]->SetDirectory(0);

  }
  m_fIn->Close();
}

void DiJetAnalysisMC::PlotSpectra( std::map< int, TH2* >& m_jznHIN,
				   int etaBinLow, int etaBinUp ){
  TCanvas c_spect("c_spect","c_spect",800,600);
  c_spect.SetLogy();

  TLegend l_spect(0.23, 0.23, 0.54, 0.36);
  StyleTools::SetLegendStyle( &l_spect, lStyleScale );
  l_spect.SetFillStyle(0);

  int style = 1;
  double max = -1;

  std::map< int, TH1* > m_jznH;

  // check if we have one
  if( !m_jznHIN.size() ){ return; }
  // should all be the same
  double etaMin = m_jznHIN.begin()->second->
    GetXaxis()->GetBinLowEdge( etaBinLow );
  double etaMax = m_jznHIN.begin()->second->
    GetXaxis()->GetBinUpEdge( etaBinUp );
  
  for( auto& jznHIN : m_jznHIN ){ // loop over JZN
    int jzn = jznHIN.first;
    
    TH1* hs =
      jznHIN.second->
      ProjectionY( Form("h_%s_%2.0f.Eta.%2.0f",
			jznHIN.second->GetName(),
			10*std::abs(etaMin),
			10*std::abs(etaMax) ),
		   etaBinLow, etaBinUp );
    m_jznH[ jzn] = hs ; 
    
    l_spect.AddEntry(  hs, Form("JZ%i", jzn) );
    StyleTools::SetHStyle( hs, style++, hStyleScale);
    hs->Draw("epsame");
    if( max < hs->GetMaximum() ){ max = hs->GetMaximum(); }

    double power = log10(max);
    power = std::ceil(power); 
    max = pow( 10, power );
    for( auto& jznH : m_jznH ){
      jznH.second->SetMaximum( max );
      jznH.second->SetMinimum( 1 );
    }
  } // end loop over JZN

  int  nPtBins = m_jznHIN.begin()->second->GetNbinsY();
  double ptMin = m_jznHIN.begin()->second->GetYaxis()->GetXmin();
  double ptMax = m_jznHIN.begin()->second->GetYaxis()->GetXmax();
  
  TH1* h_spect = new TH1D( Form("h_%s_%2.0f.Eta.%2.0f_final",
				m_jznHIN.begin()->second->GetName(),
				10*std::abs(etaMin),
				10*std::abs(etaMax) ),
			   m_jznHIN.begin()->second->GetTitle(),
			   nPtBins, ptMin, ptMax );
  StyleTools::SetHStyle( h_spect, 0, hStyleScale);
  l_spect.AddEntry( h_spect, "Total" );
  
  for( auto& jznH : m_jznH ){
    int jzn      = jznH.first;
    double scale = m_m_jznEff[ jzn ] * m_m_jznSigma[ jzn ];
    h_spect->Add( jznH.second, scale );
  }
  h_spect->Scale( 1./m_sumSigmaEff );
  h_spect->Draw("epsame");
  
  l_spect.Draw();
  DrawTools::DrawAtlasInternalMCRight( 0, 0, lStyleScale, m_is_pPb ); 
  DrawTools::DrawLeftLatex( 0.5, 0.81,
			    Form("%3.1f<#eta<%3.1f", etaMin, etaMax),
			    hStyleScale, 1 );
  
  c_spect.SaveAs( Form("%s/spectra_%2.0f.Eta.%2.0f%s.pdf",
		       m_dirOut.c_str(),
		       std::abs(etaMin)*10,
		       std::abs(etaMax)*10,
		       m_labelOut.c_str() ) );
  c_spect.SaveAs( Form("%s/spectra_%2.0f.Eta.%2.0f%s.png",
		       m_dirOut.c_str(),
		       std::abs(etaMin)*10,
		       std::abs(etaMax)*10,
		       m_labelOut.c_str() ) );

    
  c_spect.Write( Form("c_spectra_%2.0f.Eta.%2.0f%s",
		      std::abs(etaMin)*10,
		      std::abs(etaMax)*10,
		      m_labelOut.c_str() ) );
}

void DiJetAnalysisMC::PlotEtaPhiPtMap( std::map< int, TH2* >& m_jznHIN ){
  TCanvas c_map("c_map","c_map",800,600);

  for( auto& jznHIN : m_jznHIN ){
    jznHIN.second->Draw("col");
    StyleTools::SetHStyle( jznHIN.second, 0, hStyleScale);
    DrawTools::DrawAtlasInternalMCLeft( 0, -0.55, lStyleScale, true );  
    c_map.SaveAs( Form("%s/%s%s.pdf", 
		       m_dirOut.c_str(),
		       jznHIN.second->GetName(),
		       m_labelOut.c_str() ) );
    c_map.SaveAs( Form("%s/%s%s.png", 
		       m_dirOut.c_str(),
		       jznHIN.second->GetName(),
		       m_labelOut.c_str() ) );
  }
}

void DiJetAnalysisMC::PlotVsEtaPt( int etaBinLow, int etaBinUp,
				   std::map< int, TH3* >& m_jznHIN,
				   std::map< int, TH2* >& m_jznNIN,
				   int type ){
  TCanvas c_mean("c_mean","c_mean",800,600);
  TCanvas c_sigma("c_sigma","c_sigma",800,600);
  
  TLegend leg(0.68, 0.64, 0.99, 0.77);
  StyleTools::SetLegendStyle( &leg, lStyleScale );
  leg.SetFillStyle(0);

  int style = 1;

  std::map< int, TH1* > m_jznMean;
  std::map< int, TH1* > m_jznSigma;
  std::map< int, TH1* > m_jznN;

  // check if we have one
  if( !m_jznHIN.size() ){ return; }
  // should all be the same
  double etaMin = m_jznHIN.begin()->second->
    GetXaxis()->GetBinLowEdge( etaBinLow );
  double etaMax = m_jznHIN.begin()->second->
    GetXaxis()->GetBinUpEdge( etaBinUp );
  
  int  nPtBins = m_jznHIN.begin()->second->GetNbinsY();
  double ptMin = m_jznHIN.begin()->second->GetYaxis()->GetXmin();
  double ptMax = m_jznHIN.begin()->second->GetYaxis()->GetXmax();

  std::string yTitleMean;
  std::string yTitleSigma;

  if( type == 0 ){
    yTitleMean  = "#it{p}_{T}^{Reco}/#it{p}_{T}^{Truth}";
    yTitleSigma = "#sigma" + yTitleMean;
  } else if( type == 1 ){
    yTitleMean  = "#Delta#eta";
    yTitleSigma = "#sigma" + yTitleMean;
  } else if( type == 2 ){
    yTitleMean  = "#Delta#phi";
    yTitleSigma = "#sigma" + yTitleMean;
  } 
  
  for( auto& jznHIN : m_jznHIN ){   // loop over JZN
    int jzn = jznHIN.first;
    
    TH1* hMean = new TH1D( Form("%s_mean_%2.0f.Eta.%2.0f",
				jznHIN.second->GetName(),
				10*std::abs(etaMin),
				10*std::abs(etaMax) ),
			   Form(";#it{p}_{T}^{Truth};%s",
				yTitleMean.c_str()),
			   nPtBins, ptMin, ptMax );    
    StyleTools::SetHStyle( hMean, style, hStyleScale);
    m_jznMean[ jzn ] = hMean;
    
    TH1* hSigma = new TH1D( Form("%s_sigma_%2.0f.Eta.%2.0f",
				 jznHIN.second->GetName(),
				 10*std::abs(etaMin),
				 10*std::abs(etaMax) ),
			    Form(";#it{p}_{T}^{Truth};%s",
				 yTitleSigma.c_str() ),
			    nPtBins, ptMin, ptMax );    
    StyleTools::SetHStyle( hSigma, style, hStyleScale);
    m_jznSigma[ jzn ] = hSigma; 

    TH1* hN = m_jznNIN[ jzn ]->
      ProjectionY( Form("%s_N_%2.0f.Eta.%2.0f",
			jznHIN.second->GetName(),
			10*std::abs(etaMin),
			10*std::abs(etaMax) ),
		   etaBinLow, etaBinUp);
    m_jznN[ jzn ] = hN;
    
    ProjectEtaPtAndFit( jznHIN.second, hMean, hSigma, etaBinLow, etaBinUp );

    // Draw on Mean canvas
    c_mean.cd();
    leg.AddEntry(  hMean, Form("JZ%i", jzn) );
    hMean->Draw("epsame");

    if( type == 0 ){ // JES JER
      hMean->SetMinimum(0.75);
      hMean->SetMaximum(1.25);
    } else if( type == 1 || type == 2 ) { // ANGLES
      hMean->SetMinimum(-0.2);
      hMean->SetMaximum(0.2);      
    }

    // Draw on Sigma canvas
    c_sigma.cd();
    hSigma->Draw("epsame");

    if( type == 0 ){
      hSigma->SetMinimum(0.);
      hSigma->SetMaximum(0.30);
    } else if ( type == 1 || type == 2 ){
      hSigma->SetMinimum(0.);
      hSigma->SetMaximum(0.125);
    }

    // increment style
    style++;

  } // end loop over JZN

  // The name will have jz1 in it but oh well...
  // final is is important 
  TH1* h_mean = new TH1D( Form("%s_mean_%2.0f.Eta.%2.0f_final",
			       m_jznHIN.begin()->second->GetName(),
			       10*std::abs(etaMin),
			       10*std::abs(etaMax) ),
			  Form(";#it{p}_{T}^{Truth};%s",
			       yTitleMean.c_str() ),
			  nPtBins, ptMin, ptMax );
  StyleTools::SetHStyle( h_mean, 0, hStyleScale);
  
  // The name will have jz1 in it but oh well...
  // final is is important
  TH1* h_sigma = new TH1D( Form("%s_sigma_%2.0f.Eta.%2.0f_final",
				m_jznHIN.begin()->second->GetName(),
				10*std::abs(etaMin),
				10*std::abs(etaMax) ),
			   Form(";#it{p}_{T}^{Truth};%s",
				yTitleSigma.c_str() ),
			   nPtBins, ptMin, ptMax );
  StyleTools::SetHStyle( h_sigma, 0, hStyleScale);

  CombineJZN( h_mean , m_jznMean , m_jznN );
  CombineJZN( h_sigma, m_jznSigma, m_jznN );

  c_mean.cd();
  h_mean->Draw("epsame");
  leg.AddEntry( h_mean, "Total" );
  
  c_sigma.cd();
  h_sigma->Draw("epsame");
  
  // Draw the canvases for mean and sigma
  bool isMean;
  DrawCanvas(  m_jznHIN ,
	       c_mean , leg,
	       etaMin , etaMax,
	       type   , isMean = true ); 

  DrawCanvas(  m_jznHIN  ,
	       c_sigma , leg,
	       etaMin  , etaMax,
	       type    , isMean = false ); 
}

void DiJetAnalysisMC::ProjectEtaPtAndFit( TH3* h3,
					  TH1* h1Mean, TH1* h1Sigma,
					  int etaBinLow, int etaBinUp ){
  TCanvas c_proj("c_proj","c_proj",800,600);

  double etaMin = h3->GetXaxis()->GetBinLowEdge( etaBinLow );
  double etaMax = h3->GetXaxis()->GetBinUpEdge( etaBinUp );

  int  nPtBins  = h3->GetNbinsY();

  std::vector< TH1* > v_hProj;
  std::vector< TF1* > v_fit;
  
  // loop over ptBins
  for( int ptBin = 1; ptBin <= nPtBins; ptBin++ ){
    double ptMin = h3->GetYaxis()->GetBinLowEdge(ptBin);
    double ptMax = h3->GetYaxis()->GetBinUpEdge(ptBin);
      
    TH1* hProj = h3->
      ProjectionZ( Form("%s_%2.0f.Eta.%2.0f_%2.0f.Pt.%2.0f",
			h3->GetName(),
			10*std::abs(etaMin),
			10*std::abs(etaMax),
			ptMin, ptMax ),
		   etaBinLow, etaBinUp, ptBin, ptBin );
    StyleTools::SetHStyle( hProj, 0, hStyleScale);
    v_hProj.push_back( hProj );
    hProj->SetTitle("");
    
    TF1* fit  = new TF1( Form("f_%s_%2.0f.Eta.%2.0f_%2.0f.Pt.%2.0f",
			      h3->GetName(),
			      10*std::abs(etaMin),
			      10*std::abs(etaMax),
			      ptMin, ptMax ),
			 "gaus(0)" );
    StyleTools::SetHStyle( fit, 0, hStyleScale);
    v_fit.push_back( fit );

    if( hProj->GetEntries() < 5 ){ continue; }
    
    AnalysisTools::FitGaussian( hProj, fit );
        
    h1Mean->SetBinContent( ptBin, fit->GetParameter(1) );
    h1Mean->SetBinError  ( ptBin, fit->GetParError (1) );

    h1Sigma->SetBinContent ( ptBin, fit->GetParameter(2) );
    h1Sigma->SetBinError   ( ptBin, fit->GetParError (2) );
    
    hProj->Draw();
    fit->Draw("same");

    DrawTools::DrawAtlasInternalMCRight( 0, 0, lStyleScale, true ); 
    DrawTools::DrawLeftLatex( 0.5, 0.81,
			      Form("%3.1f<#eta<%3.1f", etaMin, etaMax),
			      lStyleScale, 1 );
    DrawTools::DrawLeftLatex( 0.5, 0.74,
			      Form("%3.0f<#it{p}_{T}^{Truth}<%3.1f", ptMin, ptMax ),
			      lStyleScale, 1 );

    c_proj.Write( Form("c_%s_%2.0f.Eta.%2.0f_%2.0f.Pt.%2.0f%s",
		       h3->GetName(),
		       std::abs(etaMin)*10,
		       std::abs(etaMax)*10,
		       ptMin, ptMax,
		       m_labelOut.c_str() ) ); 
  } // end loop over ptBins
}

//---------------------------
//          Tools 
//---------------------------

void DiJetAnalysisMC::CombineJZN( TH1* h_res,
				  std::map< int, TH1*>& m_jznV1,
				  std::map< int, TH1*>& m_jznN1 ){
  for( int xBin = 1; xBin <= h_res->GetNbinsX(); xBin++ ){

    double valTot   = 0; double valErrorTot   = 0;
    double denomTot = 0;
    
    for( auto jznV : m_jznV1 ){
      int jzn = jznV.first;

      double valueBin    = m_jznV1[ jzn ]->GetBinContent( xBin );
      double valueBinErr = m_jznV1[ jzn ]->GetBinError( xBin );
      double nEntriesBin = m_jznN1[ jzn ]->GetBinContent( xBin );
      double weight      = m_m_jznWeights[ jzn ];

      /*if( xBin == 10 ){
	std::cout << "  jz" << jzn << "   weight: " << weight
		<< "  valueBin: " << valueBin
		<< "  nentries: " << nEntriesBin << std::endl;
		}*/
      
      double val    = weight * nEntriesBin * valueBin;
      double valErr = weight * nEntriesBin * valueBinErr;
      valTot       += val;
      valErrorTot  += valErr * valErr;

      double denom  = weight * nEntriesBin;
      denomTot     += denom;
    }
    
    valTot      /= denomTot;
    valErrorTot /= denomTot * denomTot;

    /*   
    std::cout << h_res->GetName() << "  xbin: " << xBin
	      << "   valT: " << valTot << std::endl; 
    */
    
    h_res->SetBinContent( xBin, valTot );
    h_res->SetBinError  ( xBin, valErrorTot );
  }
}
  
void DiJetAnalysisMC::DrawCanvas( std::map< int, TH3* >& m_jznHIN,
				  TCanvas& c, TLegend& leg,
				  double etaMin, double etaMax,
				  int type, bool isMean ){
  // Draw on JES canvas
  c.cd();
  leg.Draw();
  DrawTools::DrawAtlasInternalMCRight( 0, 0, lStyleScale, true ); 
  DrawTools::DrawLeftLatex( 0.5, 0.81,
			    Form("%3.1f<#eta<%3.1f", etaMin, etaMax),
			    lStyleScale, 1 );
  
  double y0;

  std::string meanOrSigma = isMean ? "mean" : "sigma";
  
  if( type == 0 ){ // JES JER
    y0 = 1;
    y0 = 1;
  } else if( type == 1 || type == 2 ) { // ANGLES
    y0 = 0;
    y0 = 0;
  }

  double ptMin = m_jznHIN.begin()->second->GetYaxis()->GetXmin();
  double ptMax = m_jznHIN.begin()->second->GetYaxis()->GetXmax();
  
  TLine line( ptMin, y0, ptMax, y0);
  line.Draw();
    
  c.SaveAs( Form("%s/%s_%s_%2.0f.Eta.%2.0f%s.pdf",
		 m_dirOut.c_str(),
		 m_jznHIN.begin()->second->GetName(),
		 meanOrSigma.c_str(),
		 std::abs(etaMin)*10,
		 std::abs(etaMax)*10,
		 m_labelOut.c_str() ) );
  c.SaveAs( Form("%s/%s_%s_%2.0f.Eta.%2.0f%s.png",
		 m_dirOut.c_str(),
		 m_jznHIN.begin()->second->GetName(),
		 meanOrSigma.c_str(),
		 std::abs(etaMin)*10,
		 std::abs(etaMax)*10,
		 m_labelOut.c_str() ) );

  c.Write( Form("c_%s_%s_%2.0f.Eta.%2.0f%s",
		m_jznHIN.begin()->second->GetName(),
		meanOrSigma.c_str(),
		std::abs(etaMin)*10,
		std::abs(etaMax)*10,
		m_labelOut.c_str()) );  
}
