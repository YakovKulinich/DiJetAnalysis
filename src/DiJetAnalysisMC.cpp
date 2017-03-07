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

void DiJetAnalysisMC::Initialize(){
  DiJetAnalysis::Initialize();

  m_dRmax = 0.2;
  
  m_vUsedJZN.push_back( 1 );
  m_mJznFnameIn[ m_vUsedJZN.back() ] = "/home/yakov/Projects/atlas/data/jz1.root";
  m_mJznSigma  [ m_vUsedJZN.back() ] = 6.7890E+07;
  m_mJznEff    [ m_vUsedJZN.back() ] = 2.8289E-03;

  /*
  m_vUsedJZN.push_back( 2 );
  m_mJznFnameIn[ m_vUsedJZN.back() ] = "/home/yakov/Projects/atlas/data/jz2.root";
  m_mJznSigma  [ m_vUsedJZN.back() ] = 6.3997E+05;
  m_mJznEff    [ m_vUsedJZN.back() ] = 4.2714E-03;
  */
  
  // calculate sum of sigma and eff
  m_sumSigmaEff = 0;
  for( auto jzn : m_vUsedJZN ){
    m_sumSigmaEff += ( m_mJznSigma[ jzn ] * m_mJznEff[ jzn ] );
  }

  std::cout << "SumSigmaEff = " << m_sumSigmaEff << std::endl;
  
  for( auto jzn : m_vUsedJZN ){ // loop over JZ samples
    m_fIn = TFile::Open( m_mJznFnameIn[ jzn ].c_str() );
    m_tree = (TTree*) m_fIn->Get( "tree" );

    int nEventsTotal = m_tree->GetEntries();
    double     sigma = m_mJznSigma[ jzn ];      
    double       eff = m_mJznEff  [ jzn ];
  
    m_mJznNev    [ jzn ] = nEventsTotal;
    m_mJznWeights[ jzn ] =
      (1./nEventsTotal)*(sigma * eff)*(1./m_sumSigmaEff);
  
    std::cout << "JZ" << jzn << "   weight = " << m_mJznWeights [jzn] << std::endl;
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
  
  // check if we have one
  PlotSpectra( m_mJznEtaSpect );
  
  PlotEtaPhiPtMap( m_mJznEtaPhi );
  PlotEtaPhiPtMap( m_mJznEtaPt  );

  PlotVsEtaPt( m_mJznRpt , m_mJznNentries, 0);
  PlotVsEtaPt( m_mJznDeta, m_mJznNentries, 1);
  PlotVsEtaPt( m_mJznDphi, m_mJznNentries, 2);

  std::cout << "DONE! Closing " << cfNameOut << std::endl;
  m_fOut->Close();
  std::cout << "......Closed  " << cfNameOut << std::endl;
}

//---------------------------------
//            Read Data
//---------------------------------

void DiJetAnalysisMC::SetupHistograms(){
  // Triggers and Spectra
  int    nPtSpectBins = 50; 
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
  
  for( auto jzn : m_vUsedJZN ){ // loop over JZ samples
    std::cout << "Making - JZ" << jzn << " histograms " << std::endl;

    m_mJznEtaSpect[ jzn ] = 
      new TH2D( Form("h_etaSpect_jz%i", jzn ), 
		";#eta_{Reco};#it{p}_{T}^{Reco} [GeV];dN/d#it{p}_{T}",
		nEtaForwardBinsCoarse, etaForwardMin, etaForwardMax,
		nPtSpectBins, ptSpectMin, ptSpectMax ) ;
    AddHistogram( m_mJznEtaSpect[ jzn ] );
    
    m_mJznEtaPhi[ jzn ] =
      new TH2D( Form("h_etaPhi_jz%i", jzn),
		";#eta_{Reco};#phi_{Reco}",
		nEtaBins, etaMin, etaMax,
		nPhiBins, phiMin, phiMax );
    AddHistogram( m_mJznEtaPhi[ jzn ] );
      
    m_mJznEtaPt[ jzn ] =
      new TH2D( Form("h_etaPt_jz%i", jzn),
		";#eta_{Reco};#it{p}_{T}^{Reco}",
		nEtaBins, etaMin, etaMax,
		nPtBins , ptMin  , ptMax );
    AddHistogram( m_mJznEtaPt[ jzn ] );
    
    m_mJznRpt[ jzn ] =
      new TH3D( Form("h_rPt_jz%i", jzn),
		";#eta^{Truth};#it{p}_{T}^{Truth};#it{p}_{T}^{Reco}/#it{p}_{T}^{Truth}",
		nEtaForwardBinsFine, etaForwardMin, etaForwardMax,
		nPtTruthBins, ptTruthMin, ptTruthMax,
		nRPtBins, rPtMin, rPtMax);
    AddHistogram( m_mJznRpt[ jzn ] );
  
    m_mJznDeta[ jzn ] =
      new TH3D( Form("h_dEta_jz%i", jzn),
		";#eta^{Truth};#it{p}_{T}^{Truth};#eta^{Reco}-#eta^{Truth}",
		nEtaForwardBinsFine, etaForwardMin, etaForwardMax,
		nPtTruthBins, ptTruthMin, ptTruthMax,
		nDAngleBins, dAngleMin, dAngleMax );
    AddHistogram( m_mJznDeta[ jzn ] );    
  
    m_mJznDphi[ jzn ] =
      new TH3D( Form("h_dPhi_jz%i", jzn),
		";#eta^{Truth};#it{p}_{T}^{Truth};#phi^{Reco}-#phi^{Truth}",
		nEtaForwardBinsFine, etaForwardMin, etaForwardMax,
		nPtTruthBins, ptTruthMin, ptTruthMax,
		nDAngleBins, dAngleMin, dAngleMax );
    AddHistogram( m_mJznDphi[ jzn ] );

    m_mJznNentries[ jzn ] =
      new TH2D( Form("h_nEntries_jz%i", jzn),
		";#eta^{Truth};#it{p}_{T}^{Truth}",
		nEtaForwardBinsFine, etaForwardMin, etaForwardMax,
		nPtTruthBins, ptTruthMin, ptTruthMax );
    AddHistogram( m_mJznNentries[ jzn ] );
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

  for( auto jzn : m_vUsedJZN ){
    
    std::cout << "fNameIn: " << m_mJznFnameIn[ jzn ] << std::endl;
  
    m_fIn = TFile::Open( m_mJznFnameIn[ jzn ].c_str() );
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
	  
	  
	m_mJznEtaSpect[ jzn ]->Fill( etaReco, ptTruth );
	m_mJznEtaPhi  [ jzn ]->Fill( etaReco, phiReco );
	m_mJznEtaPt   [ jzn ]->Fill( etaReco, ptReco );
	  
	if( vp.DeltaR() <= m_dRmax ){
	  m_mJznRpt   [ jzn ]->Fill( etaTruth, ptTruth, ptReco/ptTruth );
	  m_mJznDeta  [ jzn ]->Fill( etaTruth, ptTruth, etaReco - etaTruth );
	  m_mJznDphi  [ jzn ]->Fill( etaTruth, ptTruth, phiReco - phiTruth );
	  m_mJznNentries[ jzn ]->Fill( etaTruth, ptTruth );
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

  for( auto& jzn : m_vUsedJZN ){
    m_mJznEtaSpect [ jzn ] =
      (TH2D*)m_fIn->Get( Form("h_etaSpect_jz%i" , jzn ) );
    m_mJznEtaSpect [ jzn ]->SetDirectory(0);
    m_mJznEtaPhi   [ jzn ] =
      (TH2D*)m_fIn->Get( Form("h_etaPhi_jz%i"   , jzn ) );
    m_mJznEtaPhi   [ jzn ]->SetDirectory(0);
    m_mJznEtaPt    [ jzn ] =
      (TH2D*)m_fIn->Get( Form("h_etaPt_jz%i"    , jzn ) );
    m_mJznEtaPt    [ jzn ]->SetDirectory(0);

    m_mJznRpt     [ jzn ] =
      (TH3D*)m_fIn->Get( Form("h_rPt_jz%i"    , jzn ) );
    m_mJznRpt     [ jzn ]->SetDirectory(0);
    m_mJznDeta    [ jzn ]   =
      (TH3D*)m_fIn->Get( Form("h_dEta_jz%i"    , jzn ) );
    m_mJznDeta    [ jzn ]->SetDirectory(0);
    m_mJznDphi    [ jzn ] =
      (TH3D*)m_fIn->Get( Form("h_dPhi_jz%i"    , jzn ) );
    m_mJznDphi    [ jzn ]->SetDirectory(0);

    m_mJznNentries [ jzn ] =
      (TH2D*)m_fIn->Get( Form("h_nEntries_jz%i"    , jzn ) );
    m_mJznNentries [ jzn ]->SetDirectory(0);

  }
  m_fIn->Close();
}

void DiJetAnalysisMC::PlotSpectra( std::map< int, TH2* >& mJznHIN ){
  std::map< int, TH1* > mJznH;

  // check if we have one
  if( !mJznHIN.size() ){ return; }

  for( int xBin = 1;
       xBin < mJznHIN.begin()->second->GetNbinsX();
       xBin++ ){

    TCanvas c_spect("c_spect","c_spect",800,600);
    c_spect.SetLogy();

    TLegend l_spect(0.15, 0.15, 0.46, 0.28);
    StyleTools::SetLegendStyle( &l_spect, StyleTools::lSS );
    l_spect.SetFillStyle(0);

    int style = 1;
    double max = -1;

    
    // should all be the same
    double etaMin = mJznHIN.begin()->second->
      GetXaxis()->GetBinLowEdge( xBin );
    double etaMax = mJznHIN.begin()->second->
      GetXaxis()->GetBinUpEdge( xBin );
  
    for( auto& jznHIN : mJznHIN ){ // loop over JZN
      int jzn = jznHIN.first;
    
      TH1* hs =
	jznHIN.second->
	ProjectionY( Form("h_%s_%2.0f.Eta.%2.0f",
			  jznHIN.second->GetName(),
			  10*std::abs(etaMin),
			  10*std::abs(etaMax) ),
		     xBin, xBin );
      mJznH[ jzn] = hs ; 
    
      l_spect.AddEntry(  hs, Form("JZ%i", jzn) );
      StyleTools::SetHStyle( hs, style++, StyleTools::lSS);
      hs->Draw("epsame");
      if( max < hs->GetMaximum() ){ max = hs->GetMaximum(); }

      double power = log10(max);
      power = std::ceil(power); 
      max = pow( 10, power );
      for( auto& jznH : mJznH ){
	jznH.second->SetMaximum( max );
	jznH.second->SetMinimum( 0.1 );
      }
    } // end loop over JZN

    l_spect.Draw();
    DrawTools::DrawAtlasInternalMCRight( 0, 0, StyleTools::lSS, m_is_pPb ); 
    DrawTools::DrawLeftLatex( 0.5, 0.4,
			      Form("%3.1f<#eta<%3.1f", etaMin, etaMax),
			      StyleTools::lSS, 1 );
  
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
  } // end loop over eta
  
  /*
  int  nPtBins = mJznHIN.begin()->second->GetNbinsY();
  double ptMin = mJznHIN.begin()->second->GetYaxis()->GetXmin();
  double ptMax = mJznHIN.begin()->second->GetYaxis()->GetXmax();
  
  TH1* h_spect = new TH1D( Form("h_%s_%2.0f.Eta.%2.0f_final",
				mJznHIN.begin()->second->GetName(),
				10*std::abs(etaMin),
				10*std::abs(etaMax) ),
			   mJznHIN.begin()->second->GetTitle(),
			   nPtBins, ptMin, ptMax );
  StyleTools::SetHStyle( h_spect, 0, StyleTools::lSS);
  l_spect.AddEntry( h_spect, "Total" );
  
  for( auto& jznH : mJznH ){
    int jzn      = jznH.first;
    double scale = m_mJznEff[ jzn ] * m_mJznSigma[ jzn ];
    std:: cout << "jz" << jzn << "   eff*sig: " << scale << std::endl;
    h_spect->Add( jznH.second, scale );
  }
  h_spect->Scale( 1./m_sumSigmaEff );
  h_spect->Draw("epsame");
  */
}

void DiJetAnalysisMC::PlotEtaPhiPtMap( std::map< int, TH2* >& mJznHIN ){
  TCanvas c_map("c_map","c_map",800,600);

  for( auto& jznHIN : mJznHIN ){
    jznHIN.second->Draw("col");
    StyleTools::SetHStyle( jznHIN.second, 0, StyleTools::lSS);
    DrawTools::DrawAtlasInternalMCLeft( 0, -0.55, StyleTools::lSS, true );  
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

void DiJetAnalysisMC::PlotVsEtaPt( std::map< int, TH3* >& mJznHIN,
				   std::map< int, TH2* >& mJznNIN,
				   int type ){
  std::map< int, TH1* > mJznMean;
  std::map< int, TH1* > mJznSigma;
  std::map< int, TH1* > mJznN;

  std::vector< TH1* > vMean;
  std::vector< TH1* > vSigma;

  // check if we have one
  if( !mJznHIN.size() ){ return; }

  for( int xBin = 1;
       xBin < mJznHIN.begin()->second->GetNbinsX();
       xBin++ ){
    
    TCanvas c_mean("c_mean","c_mean",800,600);
    TCanvas c_sigma("c_sigma","c_sigma",800,600);
  
    TLegend leg(0.68, 0.64, 0.99, 0.77);
    StyleTools::SetLegendStyle( &leg, StyleTools::lSS );
    leg.SetFillStyle(0);

    int style = 1;
    
    // should all be the same
    double etaMin = mJznHIN.begin()->second->
      GetXaxis()->GetBinLowEdge( xBin );
    double etaMax = mJznHIN.begin()->second->
      GetXaxis()->GetBinUpEdge( xBin );
  
    int  nPtBins = mJznHIN.begin()->second->GetNbinsY();
    double ptMin = mJznHIN.begin()->second->GetYaxis()->GetXmin();
    double ptMax = mJznHIN.begin()->second->GetYaxis()->GetXmax();

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
  
    for( auto& jznHIN : mJznHIN ){   // loop over JZN
      int jzn = jznHIN.first;
    
      TH1* hMean = new TH1D( Form("%s_mean_%2.0f.Eta.%2.0f",
				  jznHIN.second->GetName(),
				  10*std::abs(etaMin),
				  10*std::abs(etaMax) ),
			     Form(";#it{p}_{T}^{Truth};%s",
				  yTitleMean.c_str()),
			     nPtBins, ptMin, ptMax );    
      StyleTools::SetHStyle( hMean, style, StyleTools::lSS);
      mJznMean[ jzn ] = hMean;
    
      TH1* hSigma = new TH1D( Form("%s_sigma_%2.0f.Eta.%2.0f",
				   jznHIN.second->GetName(),
				   10*std::abs(etaMin),
				   10*std::abs(etaMax) ),
			      Form(";#it{p}_{T}^{Truth};%s",
				   yTitleSigma.c_str() ),
			      nPtBins, ptMin, ptMax );    
      StyleTools::SetHStyle( hSigma, style, StyleTools::lSS);
      mJznSigma[ jzn ] = hSigma; 

      TH1* hN = mJznNIN[ jzn ]->
	ProjectionY( Form("%s_N_%2.0f.Eta.%2.0f",
			  jznHIN.second->GetName(),
			  10*std::abs(etaMin),
			  10*std::abs(etaMax) ),
		     xBin, xBin);
      mJznN[ jzn ] = hN;
    
      ProjectEtaPtAndFit( jznHIN.second, hMean, hSigma, xBin, xBin );

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
    TH1* hMeanFinal = new TH1D( Form("%s_mean_%2.0f.Eta.%2.0f_final",
				     mJznHIN.begin()->second->GetName(),
				     10*std::abs(etaMin),
				     10*std::abs(etaMax) ),
				Form(";#it{p}_{T}^{Truth};%s",
				     yTitleMean.c_str() ),
				nPtBins, ptMin, ptMax );
    StyleTools::SetHStyle( hMeanFinal, 0, StyleTools::lSS);
    vMean.push_back( hMeanFinal );
  
    // The name will have jz1 in it but oh well...
    // final is is important
    TH1* hSigmaFinal = new TH1D( Form("%s_sigma_%2.0f.Eta.%2.0f_final",
				      mJznHIN.begin()->second->GetName(),
				      10*std::abs(etaMin),
				      10*std::abs(etaMax) ),
				 Form(";#it{p}_{T}^{Truth};%s",
				      yTitleSigma.c_str() ),
				 nPtBins, ptMin, ptMax );
    StyleTools::SetHStyle( hSigmaFinal, 0, StyleTools::lSS);
    vSigma.push_back( hSigmaFinal );
  
    /*
      CombineJZN( hMeanFinal , mJznMean , mJznN );
      CombineJZN( hSigmaFinal, mJznSigma, mJznN );
  
      c_mean.cd();
      hMeanFinal->Draw("epsame");
      leg.AddEntry( hMeanFinal, "Total" );
  
      c_sigma.cd();
      hSigmaFinal->Draw("epsame");
    */
  
    // Draw the canvases for mean and sigma
    bool isMean;
    DrawCanvas(  mJznHIN ,
		 c_mean , leg,
		 etaMin , etaMax,
		 type   , isMean = true ); 

    DrawCanvas(  mJznHIN  ,
		 c_sigma , leg,
		 etaMin  , etaMax,
		 type    , isMean = false ); 

  } // end loop over eta
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
    StyleTools::SetHStyle( hProj, 0, StyleTools::lSS);
    v_hProj.push_back( hProj );
    hProj->SetTitle("");
    
    TF1* fit  = new TF1( Form("f_%s_%2.0f.Eta.%2.0f_%2.0f.Pt.%2.0f",
			      h3->GetName(),
			      10*std::abs(etaMin),
			      10*std::abs(etaMax),
			      ptMin, ptMax ),
			 "gaus(0)" );
    StyleTools::SetHStyle( fit, 0, StyleTools::lSS);
    v_fit.push_back( fit );

    if( hProj->GetEntries() < 5 ){ continue; }
    
    AnalysisTools::FitGaussian( hProj, fit );
        
    h1Mean->SetBinContent( ptBin, fit->GetParameter(1) );
    h1Mean->SetBinError  ( ptBin, fit->GetParError (1) );

    h1Sigma->SetBinContent ( ptBin, fit->GetParameter(2) );
    h1Sigma->SetBinError   ( ptBin, fit->GetParError (2) );
    
    hProj->Draw();
    fit->Draw("same");

    DrawTools::DrawAtlasInternalMCRight( 0, 0, StyleTools::lSS, true ); 
    DrawTools::DrawLeftLatex( 0.5, 0.8,
			      Form("%3.1f<#eta<%3.1f", etaMin, etaMax),
			      StyleTools::lSS, 1 );
    DrawTools::DrawLeftLatex( 0.5, 0.73,
			      Form("%3.0f<#it{p}_{T}^{Truth}<%3.1f", ptMin, ptMax ),
			      StyleTools::lSS, 1 );


    c_proj.SaveAs( Form("%s/fits/%s_%2.0f.Eta.%2.0f_%2.0f.Pt.%2.0f%s.pdf",
			m_dirOut.c_str(),
			h3->GetName(),
			std::abs(etaMin)*10,
			std::abs(etaMax)*10,
			ptMin, ptMax,
			m_labelOut.c_str() ) );


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
				  std::map< int, TH1*>& mJznV1,
				  std::map< int, TH1*>& mJznN1 ){
  for( int xBin = 1; xBin <= h_res->GetNbinsX(); xBin++ ){

    double valTot   = 0; double valErrorTot   = 0;
    double denomTot = 0;
    
    for( auto jznV : mJznV1 ){
      int jzn = jznV.first;

      double valueBin    = mJznV1[ jzn ]->GetBinContent( xBin );
      double valueBinErr = mJznV1[ jzn ]->GetBinError( xBin );
      double nEntriesBin = mJznN1[ jzn ]->GetBinContent( xBin );
      double weight      = m_mJznWeights[ jzn ];

      double val    = weight * nEntriesBin * valueBin;
      double valErr = weight * nEntriesBin * valueBinErr;
      valTot       += val;
      valErrorTot  += valErr * valErr;

      double denom  = weight * nEntriesBin;
      denomTot     += denom;
    }
    
    valTot      /= denomTot;
    valErrorTot /= denomTot * denomTot;
    
    h_res->SetBinContent( xBin, valTot );
    h_res->SetBinError  ( xBin, valErrorTot );
  }
}

void DiJetAnalysisMC::DrawCanvas( std::map< int, TH3* >& mJznHIN,
				  TCanvas& c, TLegend& leg,
				  double etaMin, double etaMax,
				  int type, bool isMean ){
  // Draw on JES canvas
  c.cd();
  leg.Draw();
  DrawTools::DrawAtlasInternalMCRight( 0, 0, StyleTools::lSS, true ); 
  DrawTools::DrawLeftLatex( 0.5, 0.8,
			    Form("%3.1f<#eta<%3.1f", etaMin, etaMax),
			    StyleTools::lSS, 1 );
  
  double y0;

  std::string meanOrSigma = isMean ? "mean" : "sigma";
  
  if( type == 0 ){ // JES JER
    y0 = 1;
    y0 = 1;
  } else if( type == 1 || type == 2 ) { // ANGLES
    y0 = 0;
    y0 = 0;
  }

  double ptMin = mJznHIN.begin()->second->GetYaxis()->GetXmin();
  double ptMax = mJznHIN.begin()->second->GetYaxis()->GetXmax();
  
  TLine line( ptMin, y0, ptMax, y0);
  line.Draw();
    
  c.SaveAs( Form("%s/%s_%s_%2.0f.Eta.%2.0f%s.pdf",
		 m_dirOut.c_str(),
		 mJznHIN.begin()->second->GetName(),
		 meanOrSigma.c_str(),
		 std::abs(etaMin)*10,
		 std::abs(etaMax)*10,
		 m_labelOut.c_str() ) );
  c.SaveAs( Form("%s/%s_%s_%2.0f.Eta.%2.0f%s.png",
		 m_dirOut.c_str(),
		 mJznHIN.begin()->second->GetName(),
		 meanOrSigma.c_str(),
		 std::abs(etaMin)*10,
		 std::abs(etaMax)*10,
		 m_labelOut.c_str() ) );

  c.Write( Form("c_%s_%s_%2.0f.Eta.%2.0f%s",
		mJznHIN.begin()->second->GetName(),
		meanOrSigma.c_str(),
		std::abs(etaMin)*10,
		std::abs(etaMax)*10,
		m_labelOut.c_str()) );  
}

void DiJetAnalysisMC::DrawCanvas( std::vector< TH1* >& mHIN,
				  TCanvas& c, TLegend& leg,
				  int type, bool isMean ){

}
