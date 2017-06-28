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
#include <TVirtualFitter.h>
#include <TRandom.h>

#include <cmath>
#include <iostream>
#include <sstream>
#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <boost/assign.hpp>

#include "DiJetAnalysis.h"
#include "DeltaPhiProj.h"

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldBinByBin.h"

DiJetAnalysis::DiJetAnalysis() : DiJetAnalysis( true, true, 0 ){}

DiJetAnalysis::DiJetAnalysis( bool is_pPb )
  : DiJetAnalysis( is_pPb, false, 0 ){}

DiJetAnalysis::DiJetAnalysis( bool is_pPb, bool isData )
  : DiJetAnalysis( is_pPb, isData, 0 ){}

DiJetAnalysis::DiJetAnalysis( bool is_pPb, bool isData, int mcType )
  : m_is_pPb( is_pPb ), m_isData( isData ), m_mcType( mcType ),
    m_fIn(NULL), m_fOut(NULL), m_tree(NULL)
{
  //========== Set Histogram Binning =============
  // Common for all analysis
  
  // -------- maps ---------
  m_nEtaMapBins = 100; 
  m_etaMapMin   = constants::ETAMIN;
  m_etaMapMax   = constants::ETAMAX;
  
  m_nPhiMapBins = 64; 
  m_phiMapMin   = -constants::PI;
  m_phiMapMax   = constants::PI; 

  m_nPtMapBins  = 45;
  m_ptMapMin    = 10;
  m_ptMapMax    = 100;
  
  // -------- spect --------
  m_ptSpectMin   = 10;
  m_ptSpectMax   = 110;
  m_nPtSpectBins = ( m_ptSpectMax - m_ptSpectMin ) / 2 ; 
  
  // ---- JES/PRes/Etc ----- 
  m_etaForwardMin = -constants::FETAMAX;
  m_etaForwardMax = -constants::FYSTARMIN;
  m_nEtaForwardBinsFine = ( m_etaForwardMax - m_etaForwardMin ) / 0.2 ;
  
  // ---- forward eta binning ---
  boost::assign::push_back( m_varFwdEtaBinning )
    ( -4.0 )( -3.3)( -3.1 )( -2.7 );
  m_nVarFwdEtaBins = m_varFwdEtaBinning.size() - 1;
  
  // -------- eff ---------
  m_effMin = 0.;
  m_effMax = 1.4;

  // --- variable eta/ystar binning ---
  // ystarB
  boost::assign::push_back( m_varYstarBinningA )
    ( -4.0 )( -2.7 )( -1.8 )( 0 )( 1.8 )( 4.0 );
  m_nVarYstarBinsA = m_varYstarBinningA.size() - 1;
  
  // ystarA
  boost::assign::push_back( m_varYstarBinningB )
    ( -4.0 )( -2.7 )( -1.8 )( 0 )( 1.8 )( 4.0 );
  m_nVarYstarBinsB = m_varYstarBinningB.size() - 1;
  
  // --- variable pt binning ---
  boost::assign::push_back( m_varPtBinning )
    ( 30 )( 38 )( 45 )( 90 );
  m_nVarPtBins = m_varPtBinning.size() - 1;

  // --- variable dphi binning ---
  // coarse bin width = 4 x fine bin width
  m_dPhiFineCoarseFactor = 4;
  // number of bins closer to 0, where we have
  // small statistics.
  m_nDphiCoarseBins = 8;
  // number of bins closer to pi, where we have
  // bigger statistics
  m_nDphiFineBins   = 32;
  // 8coarse bins + 32fine bins = 8*4+32 = 64 bins;
  // bin width for fine bins should be ~0.05
  m_dPhiFineBinWidth = constants::PI /
    ( m_nDphiFineBins + m_nDphiCoarseBins * m_dPhiFineCoarseFactor );

  for( int coarseBin = 0; coarseBin < m_nDphiCoarseBins; coarseBin++ )
    { m_varDphiBinning.push_back
	( coarseBin * m_dPhiFineBinWidth * m_dPhiFineCoarseFactor ); }
  double fineBinStart =
    m_nDphiCoarseBins * m_dPhiFineBinWidth * m_dPhiFineCoarseFactor;
  for( int fineBin = 0; fineBin <= m_nDphiFineBins; fineBin++ )
    { m_varDphiBinning.push_back
	( fineBin * m_dPhiFineBinWidth + fineBinStart ); }

  m_nVarDphiBins = m_varDphiBinning.size() - 1;
  std::cout << " ----------------- " << m_nVarDphiBins << std::endl;
  
  // --- dPhiBins ---  
  boost::assign::push_back( m_nDphiBins )
    ( m_nVarYstarBinsA )( m_nVarYstarBinsB )
    ( m_nVarPtBins     )( m_nVarPtBins     )
    ( m_nVarDphiBins   );
    
  boost::assign::push_back( m_dPhiMin  )
    ( 0 )( 0 )( 0 )( 0 )( 0 );

  boost::assign::push_back( m_dPhiMax  )
    ( 1 )( 1 )( 1 )( 1 )( 1 );

  m_nDphiDim     = m_nDphiBins.size();
  
  m_dPhiDphiMin  = 0;
  m_dPhiDphiMax  = constants::PI;
  
  m_dPhiWidthMin = 0.00;
  m_dPhiWidthMax = 0.5;

  m_dPhiZoomLow  = 2.0;
    
  //========= Set DeltaPhi Axes Order ============
  // The DeltaPhiProj object will have the order
  // onto which to take projections. It also knows
  // what axis is what. The object knows which axis
  // has what name and range. I.e. Pt1 is actually
  // The third axis in the THnSparse. 
  m_dPP = new DeltaPhiProj( YS1, PT1, PT2, YS2 );

  m_dPP->AddTAxis( new TAxis( m_nVarYstarBinsA, &m_varYstarBinningA[0] ) );
  m_dPP->AddTAxis( new TAxis( m_nVarYstarBinsB, &m_varYstarBinningB[0] ) );
  m_dPP->AddTAxis( new TAxis( m_nVarPtBins    , &m_varPtBinning[0]     ) );
  m_dPP->AddTAxis( new TAxis( m_nVarPtBins    , &m_varPtBinning[0]     ) );
  
  //==================== Cuts ====================
  m_nMinEntriesFit = 20;

  m_dPhiThirdJetFraction = 0.4;

  //========= common tools ========
  anaTool   = new CT::AnalysisTools();
  drawTool  = new CT::DrawTools();
  styleTool = new CT::StyleTools();

  //========= config file =========
  std::string configName   = "config/configJetsFwd";
  
  std::string configDataType = m_is_pPb ? "_pPb" : "_pp";
  configName += configDataType;

  std::string configSuffix = m_isData ? "_data.cfg" : "_mc.cfg" ;
  configName += configSuffix;

  std::cout << "Reading " << configName << std::endl;

  m_config = new TEnv();
  m_config->ReadFile( configName.c_str(), EEnvLevel(0));

  //========== settings ===========  
  // name for "All" histos
  // this is either merged JZN or Data from Triggers
  m_allName = "All";
  
  //=============== Histo Names ==================    
  m_etaSpectName      = "etaSpect";
  m_dPhiName          = "dPhi";
  m_effName           = "eff";

  m_recoName          = "reco";
  m_truthName         = "truth";
  m_respMatName       = "respMat";
  m_unfoldedName      = "unfolded";
    
  m_dPhiRecoName      = m_dPhiName + "_" + m_recoName;
  m_dPhiTruthName     = m_dPhiName + "_" + m_truthName;
  m_dPhiRespMatName   = m_dPhiName + "_" + m_respMatName;;

  m_dPhiUnfoldedName     = m_dPhiName     + "_" + m_unfoldedName;;
  m_dPhiRecoUnfoldedName = m_dPhiRecoName + "_" + m_unfoldedName;;
}

DiJetAnalysis::~DiJetAnalysis(){
  delete m_dPP    ; m_dPP     = NULL;
  delete anaTool  ; anaTool   = NULL;
  delete drawTool ; drawTool  = NULL;
  delete styleTool; styleTool = NULL;
}

void DiJetAnalysis::Initialize(){  
  m_labelOut = m_is_pPb ? "pPb" : "pp" ;
  m_labelOut = m_isData ? m_labelOut + "_data" : m_labelOut + "_mc";
  
  m_dirOut   = "output";
  anaTool->CheckWriteDir( m_dirOut.c_str() );
  m_dirOut   += "/output_" + m_labelOut;
  anaTool->CheckWriteDir( m_dirOut.c_str() );

  m_rootFname = m_dirOut + "/myOut_" + m_labelOut + ".root";
  
  std::cout << "fNameIn/Out: " << m_rootFname << std::endl;

  std::string unfoldingMC      = GetConfig()->GetValue( "unfoldingMC"     , "" );
  std::string unfoldingMCLabel = GetConfig()->GetValue( "unfoldingMCLabel", "" );
  
  std::string output    = "output";  
  std::string system    = m_is_pPb ? "pPb" : "pp";

  m_fNameUnfoldingMC = Form( "%s/%s_%s_mc_%s/c_myOut_%s_mc_%s.root",
			     output.c_str(), output.c_str(),
			     system.c_str(), unfoldingMC.c_str(),
			     system.c_str(), unfoldingMC.c_str() );
}

//---------------------------
//       Fill Tree
//---------------------------
void DiJetAnalysis::AddHistogram( TH1* h ){
  v_hists.push_back( h );
  h->Sumw2();
  styleTool->SetHStyle( h, 0 );

  TH1* h1 = dynamic_cast< TH1* >(h);
  if( h1 ){ h->GetXaxis()->SetNdivisions(505); }

  TH2* h2 = dynamic_cast< TH2* >(h);
  if( h2 ){ h->GetYaxis()->SetNdivisions(505); }

  TH3* h3 = dynamic_cast< TH3* >(h);
  if( h3 ){ h->GetZaxis()->SetNdivisions(505); }
}

void DiJetAnalysis::AddHistogram( THnSparse* hn ){
  v_hns.push_back( hn );
  hn->Sumw2();
  hn->GetAxis(0)->SetTitle("#it{y}_{1}*");
  hn->GetAxis(1)->SetTitle("#it{y}_{2}*");
  hn->GetAxis(2)->SetTitle("#it{p}_{T}^{1}");
  hn->GetAxis(3)->SetTitle("#it{p}_{T}^{2}");
  if( hn->GetNdimensions() == 5 ){
    hn->GetAxis(4)->SetTitle("|#Delta#phi|"); }
  if( hn->GetNdimensions() == 6 ){
    hn->GetAxis(4)->SetTitle("|#Delta#phi_{Reco}|");
    hn->GetAxis(5)->SetTitle("|#Delta#phi_{Truth}|");
  }
}


void DiJetAnalysis::SaveOutputsFromTree(){
  //----------------------------------------
  //  Close the input file, 
  //  write histos to output
  //----------------------------------------
  std::cout << "fNameOut: " << m_rootFname << std::endl;
  TFile* m_fOut = new TFile( m_rootFname.c_str(),"RECREATE");
  for( auto& h  : v_hists  ){ h-> Write(); }
  for( auto& f  : v_functs ){ f-> Write(); }
  for( auto& gr : v_graphs ){ gr->Write(); }
  for( auto& hn : v_hns    ){ hn->Write(); }
  m_fOut->Close();
}

//---------------------------
//       Analysis
//---------------------------

bool DiJetAnalysis::GetFwdCentJets( const std::vector< TLorentzVector>& v_jets, 
				    const TLorentzVector*& jetFwd,
				    const TLorentzVector*& jetCent ){
  for( auto& jet : v_jets ){
    if( !jetFwd && IsForwardDetector( jet.Eta() ) && jet.Pt() > 0 )
      { jetFwd  = &jet; }
    else if( !jetCent && IsCentralDetector( jet.Eta() ) && jet.Pt() > 0 )
      { jetCent = &jet; }
    if( jetFwd && jetCent )
      { break; }
  }
  // check we have at least one jet 
  if( !jetFwd && !jetCent ){ return false; }
  return true;
}

bool DiJetAnalysis::GetDiJets( const std::vector< TLorentzVector >& v_jets, 
			       const TLorentzVector*& jet1,
			       const TLorentzVector*& jet2 ){
  for( auto& jet : v_jets ){
    if( !jet1 && jet.Pt() > 0 )
      { jet1 = &jet; } 
    else if( jet1 && jet.Pt() > 0 )
      { jet2 = &jet;
	break; }
  }
  
  // make sure we have two jets
  if( !jet1 || !jet2 )
    { return false; }
  return true;
}

double DiJetAnalysis::AnalyzeDeltaPhi( THnSparse* hn, THnSparse* hnNent,
				       const std::vector< TLorentzVector >& v_jets,
				       WeightFcn weightFcn ){
  const TLorentzVector* jet1 = NULL; const TLorentzVector* jet2 = NULL;

  if( !GetDiJets( v_jets, jet1, jet2 ) )
    { return -1; }

  double jet1_pt    = jet1->Pt()/1000.;
  double jet1_eta   = jet1->Eta();
  double jet1_phi   = jet1->Phi();
  double jet1_ystar = GetYstar( *jet1 );

  double jet2_pt    = jet2->Pt()/1000.;
  double jet2_ystar = GetYstar( *jet2 );
  
  double deltaPhi = anaTool->DeltaPhi( *jet2, *jet1 );
  
  std::vector< double > x;
  x.resize( hn->GetNdimensions() );
    
  // wont change unless we have a weightFcn
  // and then it varys depending on eta, phi, pt.
  double weight = weightFcn ? weightFcn( jet1_eta, jet1_phi, jet1_pt ) : 1;   

  x[0] = jet1_ystar;  
  x[1] = jet2_ystar;
  x[2] = jet1_pt ;
  x[3] = jet2_pt ;
  x[4] = deltaPhi;
  hn    ->Fill( &x[0], weight );
  hnNent->Fill( &x[0], 1 );

  // for pp, fill twice. once for each side since
  // it is symmetric in pp. For pPb, continue
  if( m_is_pPb ){ return deltaPhi; }
  
  x[0] = -jet1_ystar;  
  x[1] = -jet2_ystar;
  hn    ->Fill( &x[0], weight );
  hnNent->Fill( &x[0], 1 );    
 
  return deltaPhi;
}

bool DiJetAnalysis::ApplyIsolation( std::vector<TLorentzVector>& v_jets,
				    double Rmin ){
  bool haveIsolated = false;
  std::vector<bool> isIsolated( v_jets.size(), true );
  
  for(unsigned int iTestJet = 0; iTestJet < v_jets.size(); iTestJet++){
    for(unsigned int iSecondJet = iTestJet; iSecondJet < v_jets.size(); iSecondJet++){   
      if( iSecondJet == iTestJet ) continue;

      if( anaTool->DeltaR( v_jets.at(iTestJet),
			   v_jets.at(iSecondJet)) < Rmin &&
	  v_jets.at(iSecondJet).Pt() > v_jets.at(iTestJet).Pt() * 0.5 ){       
	isIsolated[ iSecondJet ] = false;
	haveIsolated = true;
	continue;
      }
    } // end loop over iSecondJet
  } // end loop over iTestJet

  for(unsigned int iJet = 0; iJet < v_jets.size(); iJet++ ){
    if( !isIsolated.at(iJet) ) v_jets.at(iJet).SetPxPyPzE(0,0,0,-1 );
  }

  return haveIsolated; 
}

bool DiJetAnalysis::ApplyCleaning( std::vector<TLorentzVector>& v_jets, 
				   std::vector<bool>& v_isCleanJet){
  bool haveCleaned = false;
  
  for( unsigned int jn = 0; jn < v_jets.size() ; jn++ ){
    if( !v_isCleanJet.at(jn) ){
      v_jets.at(jn).SetPxPyPzE(0,0,0,-1 );
      haveCleaned = true; 
    }
  }
  
  return haveCleaned;
}

//---------------------------
//       Tools
//---------------------------
void  DiJetAnalysis::FillHistoWithJets( const TLorentzVector* jet1,
					const TLorentzVector* jet2,
					TH2* h, double weight ){
  if( jet1 )
    { h->Fill( AdjustEtaForPP( jet1->Eta() ), jet1->Pt()/1000., weight ); }
  if( jet2 )
    { h->Fill( AdjustEtaForPP( jet2->Eta() ), jet2->Pt()/1000., weight ); }
} 

double DiJetAnalysis::AdjustEtaForPP( double jetEta ){
  if( m_is_pPb ) return jetEta;
  return jetEta  > 0 ? -1 * jetEta  : jetEta;
}

double DiJetAnalysis::GetYstar( const TLorentzVector& jet )
{ return m_is_pPb ? jet.Rapidity() + constants::BETAZ : jet.Rapidity(); }

bool DiJetAnalysis::IsForwardDetector( const double& eta ){
  if( std::abs( eta ) < constants::FETAMAX &&
      std::abs( eta ) > constants::FETAMIN ){ return true; }
  return false;
}

bool DiJetAnalysis::IsCentralDetector( const double& eta )
{ return ( std::abs(eta) < constants::CETAMAX ); }

// These two are just used to see which histograms to draw
// for delta Phi plots, side doesnt really matter
// because can have forward backward correlation.
bool DiJetAnalysis::IsForwardYstar( const double& ystar )
{ return ( std::abs(ystar - constants::BETAZ) < constants::FETAMAX && 
	       std::abs(ystar - constants::BETAZ) > constants::FETAMIN ); }

bool DiJetAnalysis::IsCentralYstar( const double& ystar )
{ return ( std::abs(ystar - constants::BETAZ) < constants::CETAMAX ); }

void DiJetAnalysis::NormalizeDeltaPhi( TH1* h )
{ if( h->GetEntries() ){ h->Scale( 1./h->Integral()) ; } }

TH1* DiJetAnalysis::CombineSamples( std::vector< TH1* >& vSampleHin,
				    const std::string& name ){ return NULL; }  

TH2* DiJetAnalysis::CombineSamples( std::vector< TH2* >& vSampleHin,
				    const std::string& name ){ return NULL; }  

THnSparse* DiJetAnalysis::CombineSamples( std::vector< THnSparse* >& vSampleHin,
					  const std::string& name ){ return NULL; }  


void DiJetAnalysis::GetInfoBoth( std::string& outSuffix,
				 std::string& name_a  , std::string& name_b  ,
				 std::string& label_a , std::string& label_b ,
				 std::string& suffix_a, std::string& suffix_b ){}

void DiJetAnalysis::GetInfoUnfolding( std::string& measuredName,
				      std::string& measuredLabel ){}

//---------------------------
//   Get Quantities / Plot 
//---------------------------

void DiJetAnalysis::MakeSpectra( std::vector< TH2* >& vSampleSpect,
				 const std::vector< std::string>& vLabels,
				 const std::string& name ){
  if( !vSampleSpect.size() ){ return; }
  
  std::string yAxisTitle = "dN/d#it{p}_{T}";

  double ptSpectWidth =
    ( m_ptSpectMax - m_ptSpectMin ) / m_nPtSpectBins;

  // use this as reference because
  // it should be in every file
  TH2* hRef = vSampleSpect[0];
  int nXbins = hRef->GetNbinsX();

  uint nSamples = vSampleSpect.size();
  
  std::vector< std::vector< TH1* > > vSpect;
  vSpect.resize( nSamples );

  double max = -1;
  
  for( uint iG = 0; iG < nSamples; iG++){
    std::string label = vLabels[iG];

    for( int xBin = 1; xBin <= nXbins; xBin++ ){
      double etaMin = hRef->GetXaxis()->GetBinLowEdge( xBin );
      double etaMax = hRef->GetXaxis()->GetBinUpEdge ( xBin );

      TH1* h_etaSpect =
	vSampleSpect[iG]->
	ProjectionY( Form("h_%s_%s_%s",
			  name.c_str(),
			  label.c_str(),
			  anaTool->GetName( etaMin, etaMax, "Eta").c_str() ),
		     xBin, xBin );
      h_etaSpect->SetTitle
	( anaTool->GetEtaLabel( etaMin, etaMax, m_is_pPb ).c_str() );
      h_etaSpect->GetYaxis()->SetTitle( yAxisTitle.c_str() );
      h_etaSpect->Scale( 1./ptSpectWidth );
      vSpect[iG].push_back( h_etaSpect );
 
      // get min max from the final histograms
      if( label.compare( m_allName ) ){ continue; }
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
  //------- Draw Eta as Fucntion of Labels -------
  //------------------------------------------------
  double lX0, lY0, lX1, lY1;
  
  if( m_is_pPb ){ lX0 = 0.70; lY0 = 0.54; lX1 = 0.85; lY1 = 0.71; }
  else          { lX0 = 0.20; lY0 = 0.23; lX1 = 0.47; lY1 = 0.40; }
  
  for( uint iG = 0; iG < nSamples; iG++ ){
    std::string label = vLabels[iG];

    if( label.compare( m_allName ) ){ continue; };
    
    std::string cName  = label;
    std::string cLabel = label;

    TCanvas c( "c", cLabel.c_str(), 800, 600 );
    c.SetLogy();
    
    TLegend leg( lX0, lY0, lX1, lY1 );
    styleTool->SetLegendStyle( &leg, 0.65 );

    int style = 0;
    for( int iX = 0; iX < nXbins; iX++ ){
      int         xBin = iX + 1;
      double etaCenter = hRef->GetXaxis()->GetBinCenter ( xBin );

      // for pPb, dont draw at anything above -3.2
      if( m_is_pPb && etaCenter > -constants::FETAMIN ){ continue; }

      // for pp, dont draw central triggers below -3.2
      // or forward triggers above -3.2
      // for mb, and total draw everything
      bool isMb  = label.find("_mb_") != std::string::npos
	? true : false;
      bool isAll = !label.compare( m_allName )
	? true : false;
      
      if( !m_is_pPb && etaCenter < -constants::FETAMIN &&
	  label.find("320eta490") == std::string::npos &&
	  !isMb && !isAll ){ continue; }
      else if( !m_is_pPb && etaCenter > -constants::FETAMIN &&
	       label.find("320eta490") != std::string::npos &&
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

    DrawAtlasRight();    

    SaveAsAll( c, name, cName );
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
    double etaMin = hRef->GetXaxis()->GetBinLowEdge( xBin );
    double etaMax = hRef->GetXaxis()->GetBinUpEdge ( xBin );
    double etaCenter = hRef->GetXaxis()->GetBinCenter ( xBin );

    // for pPb, dont draw at anything above -3.2
    if( m_is_pPb && etaCenter > -constants::FETAMIN ){ continue; }
     
    std::string cName  = anaTool->GetName( etaMin, etaMax, "Eta" );
    std::string cLabel = anaTool->GetEtaLabel( etaMin, etaMax, m_is_pPb );
    
    TCanvas c( "c", cLabel.c_str(), 800, 600 );
    c.SetLogy();
    
    TLegend leg( lX0, lY0, lX1, lY1 );
    styleTool->SetLegendStyle( &leg  );
    
    int style = 1;
    for( uint iG = 0; iG < nSamples; iG++ ){
      std::string label = vLabels[iG];

      // for pp, dont draw central triggers below -3.2
      // or forward triggers above -3.2
      // for mb, and total draw everything
      bool isMb  = label.find("_mb_") != std::string::npos
	? true : false;
      bool isAll = !label.compare( m_allName )
	? true : false;
           
      if( !m_is_pPb && etaCenter < -constants::FETAMIN &&
	  label.find("320eta490") == std::string::npos &&
	  !isMb && !isAll ){ continue; }
      else if( !m_is_pPb && etaCenter > -constants::FETAMIN &&
	       label.find("320eta490") != std::string::npos &&
	       !isMb && !isAll ){ continue; }
      
      TH1* h = vSpect[iG][iX];
      if( iG == nSamples ){ style = 0; }
      styleTool->SetHStyle( h, style++ );
      h->Draw("epsame");
      h->SetMinimum( 1 );
      h->SetMaximum( max );
      leg.AddEntry( h, vLabels[iG].c_str() );
    } // end loop over iG
    
    leg.Draw("same");

    DrawAtlasRight();    
    drawTool->DrawRightLatex( 0.4, 0.2, cLabel );

    SaveAsAll( c, name, cName );
  } // end loop over iX

  // delete
  for( uint iG = 0; iG < nSamples; iG++ ){
    for( int iX = 0; iX < nXbins; iX++ )
      { delete vSpect[iG][iX]; }
  }
} 

void DiJetAnalysis::MakeDeltaPhi( std::vector< THnSparse* >& vhn,
				  const std::vector< std::string >& vLabel,
				  const std::string& name ){

  std::vector< TH1* > vDphi;
  std::vector< TF1* > vFits;
  
  // ---- loop over group  ----
  // ---- (jzn or trigger) ----
  for( uint iG = 0; iG < vhn.size(); iG++ ){      
    THnSparse* hn = vhn[iG];
  
    // truncate dphi bins with not many contents.
    // anaTool->TruncateHistoBins( hn, hnNent );

    std::string label = vLabel[iG];

    // in data only draw for all
    if( label.compare( m_allName ) ){ continue; } 
    
    TAxis* axis0 = hn->GetAxis( m_dPP->GetAxisI(0) );
    TAxis* axis1 = hn->GetAxis( m_dPP->GetAxisI(1) );
    TAxis* axis2 = hn->GetAxis( m_dPP->GetAxisI(2) );
    TAxis* axis3 = hn->GetAxis( m_dPP->GetAxisI(3) );
    
    int nAxis0Bins = axis0->GetNbins();
    int nAxis1Bins = axis1->GetNbins();
    int nAxis2Bins = axis2->GetNbins();
    int nAxis3Bins = axis3->GetNbins();

    int fAxisI     = m_dPP->GetAxisI(3);
 
    // ---- loop over ystars ----
    for( int axis0Bin = 1; axis0Bin <= nAxis0Bins; axis0Bin++ ){
      // set ranges
      axis0->SetRange( axis0Bin, axis0Bin );

      double axis0Low, axis0Up;
      anaTool->GetBinRange
	( axis0, axis0Bin, axis0Bin, axis0Low, axis0Up );

      for( int axis1Bin = 1; axis1Bin <= nAxis1Bins; axis1Bin++ ){
	// check we are in correct ystar and pt bins
	if( !m_dPP->CorrectPhaseSpace
	    ( std::vector<int>{ axis0Bin, axis1Bin, 0, 0 } ) )
	  { continue; }

	// set ranges
	axis1->SetRange( axis1Bin, axis1Bin );
	
	double axis1Low, axis1Up;
	anaTool->GetBinRange
	  ( axis1, axis1Bin, axis1Bin, axis1Low, axis1Up );
	
	std::vector< TH1* > vDphiWidthsTemp;
	TCanvas cWidths("cWidths","cWidths",800,600);
	TLegend leg(0.68, 0.64, 0.99, 0.77);
	int style = 0;
	styleTool->SetLegendStyle( &leg );

	std::string hTagC =
	  Form( "%s_%s",
	        anaTool->GetName( axis0Low, axis0Up, m_dPP->GetAxisName(0) ).c_str(),
		anaTool->GetName( axis1Low, axis1Up, m_dPP->GetAxisName(1) ).c_str() ); 
	
	// ---- loop over axis2 ----
	for( int axis2Bin = 1; axis2Bin <= nAxis2Bins; axis2Bin++ ){
	  // set ranges
	  axis2->SetRange( axis2Bin, axis2Bin );
	  // set range back to whole range otherwise histograms
	  // drawn for only one bin ( from projection ) beacuse
	  // loop over axis3 ends at last bin in.
	  axis3->SetRange( 1, -1 );
 
	  double axis2Low , axis2Up;
	  anaTool->GetBinRange
	    ( axis2, axis2Bin, axis2Bin, axis2Low, axis2Up );

	  std::string hTagW =
	  Form( "%s_%s_%s",
	     	anaTool->GetName( axis0Low, axis0Up, m_dPP->GetAxisName(0) ).c_str(),
		anaTool->GetName( axis1Low, axis1Up, m_dPP->GetAxisName(1) ).c_str(),
		anaTool->GetName( axis2Low, axis2Up, m_dPP->GetAxisName(2) ).c_str() ); 

	  
	  TH1* hDphiWidths = hn->Projection( fAxisI );
	  hDphiWidths->Reset();
	  
	  hDphiWidths->SetName
	    ( Form( "h_%s_%s_%s", name.c_str(), label.c_str(), hTagW.c_str() ) );
	  hDphiWidths->GetYaxis()->SetTitle( "|#Delta#phi| width" );
	  hDphiWidths->SetMarkerSize( hDphiWidths->GetMarkerSize() * 1.5 );
	  vDphiWidthsTemp.push_back( hDphiWidths );
	  styleTool->SetHStyle( hDphiWidths, style++ );
	  	  
	  leg.AddEntry
	    ( hDphiWidths, anaTool->GetLabel( axis2Low, axis2Up, m_dPP->GetAxisLabel(2) ) .c_str() );
	  // ---- loop over axis3 ----
	  for( int axis3Bin = 1; axis3Bin <= nAxis3Bins; axis3Bin++ ){
	    // check we are in correct ystar and pt bins
	    if( !m_dPP->CorrectPhaseSpace
		( std::vector<int>{ axis0Bin, axis1Bin, axis2Bin, axis3Bin } ) )
	      { continue; }

	    axis3->SetRange( axis3Bin, axis3Bin );

	    double axis3Low , axis3Up;
	    anaTool->GetBinRange
	      ( axis3, axis3Bin, axis3Bin, axis3Low, axis3Up );
	    
	    std::string hTag =
	      Form( "%s_%s_%s_%s",
		    anaTool->GetName( axis0Low, axis0Up, m_dPP->GetAxisName(0) ).c_str(),
		    anaTool->GetName( axis1Low, axis1Up, m_dPP->GetAxisName(1) ).c_str(),
		    anaTool->GetName( axis2Low, axis2Up, m_dPP->GetAxisName(2) ).c_str(),
		    anaTool->GetName( axis3Low, axis3Up, m_dPP->GetAxisName(3) ).c_str() ); 

	    // Take projection onto the dPhi axis
	    TH1* hDphi = hn->Projection( 4 );
	    hDphi->SetName
	      ( Form( "h_%s_%s_%s", name.c_str(), label.c_str(), hTag.c_str() ) );
	    // because variable bin width, scale by bin width
	    hDphi->Scale( 1.0, "width" );
	    styleTool->SetHStyle( hDphi, 0 );
	    vDphi.push_back( hDphi );

	    // subtract combinatoric contribution before normalizing
	    anaTool->SubtractCombinatoric( hDphi );
	    
	    // Save un-normalized result for
	    // unfolding later.
	    // NN = Not Normalized
	    TH1* hDphiNN = static_cast< TH1D* >
	      ( hDphi->Clone( Form("%s_NN", hDphi->GetName() ) ) );
	    hDphiNN->Write();

	    TCanvas c( "c", hDphi->GetName(), 800, 600 );
	    
	    // Normalize
	    NormalizeDeltaPhi( hDphi );

	    hDphi->Draw();
	    hDphi->GetYaxis()->SetTitle("Normalized Count");
	    hDphi->SetTitle("");
    
	    // now fit
	    TF1* fit = anaTool->FitDphi( hDphi );
	    styleTool->SetHStyle( fit, 0 );
	    vFits.push_back( fit );

	    fit->Draw("same");

	    double chi2NDF = fit->GetChisquare()/fit->GetNDF();

	    hDphi->GetXaxis()->SetRangeUser( m_dPhiZoomLow, m_dPhiDphiMax );

	    drawTool->DrawLeftLatex( 0.5, 0.66, Form( "#Chi^{2}/NDF=%4.2f", chi2NDF ) );

	    DrawTopLeftLabels
	      ( m_dPP, axis0Low, axis0Up, axis1Low, axis1Up,
		axis2Low, axis2Up, axis3Low, axis3Up, 0.8 );
	    
	    DrawAtlasRight();
	    
	    SaveAsROOT( c, hDphi->GetName() );
	    hDphi->Write();
	    fit->Write();
	    	    
	    if( fit->GetParameter(1) < 0 )
	      { continue; }

	    // Now, put results on histogram
	    // Check Chi2/NDF. If it is 0.5 < Chi2/NDF < 2.5
	    // We save the result. Otherwise, it is a bad fit.
	    if( chi2NDF < 0.5 || chi2NDF > 2.5 ){ continue; }
	    hDphiWidths->SetBinContent( axis3Bin, fit->GetParameter(1) );
	    hDphiWidths->SetBinError  ( axis3Bin, fit->GetParError (1) );	    
	  } // end loop over axis3
	} // end loop over axis2
	
	cWidths.cd();
	for( auto& h : vDphiWidthsTemp ){
	  h->SetMinimum( m_dPhiWidthMin );
	  h->SetMaximum( m_dPhiWidthMax );
	  h->GetYaxis()->SetNdivisions( 505 );
	  h->SetTitle("");
	  h->Draw("epsame");
	  h->Write();
	}
	leg.Draw("same");

	DrawTopLeftLabels
	  ( m_dPP, axis0Low, axis0Up, axis1Low, axis1Up,
	    0, 0, 0, 0, 0.8 );

	DrawAtlasRight();
	
        SaveAsAll( cWidths, Form("h_%s_%s_%s", name.c_str(), label.c_str(), hTagC.c_str() ) );
      } // end loop over axis1     
    } // end loop over axis0
  } // end loop over iG
}

void DiJetAnalysis::MakeDphiTogether(){

  std::string outSuffix;
  std::string name_a  , name_b  ;
  std::string label_a , label_b ;
  std::string suffix_a, suffix_b;

  GetInfoBoth( outSuffix, name_a, name_b, label_a, label_b, suffix_a, suffix_b );
  
  // Check if the directories exist.
  // If they don't, create them
  std::string output = "output";
  std::string outDir = output;
  anaTool->CheckWriteDir( outDir.c_str() );
  outDir += "/all";
  anaTool->CheckWriteDir( outDir.c_str() );
  outDir += "/" + outSuffix;
  anaTool->CheckWriteDir( outDir.c_str() );
  
  TAxis* axis0 = m_dPP->GetTAxis( 0 );
  TAxis* axis1 = m_dPP->GetTAxis( 1 );
  TAxis* axis2 = m_dPP->GetTAxis( 2 );
  TAxis* axis3 = m_dPP->GetTAxis( 3 );
  
  int nAxis0Bins = axis0->GetNbins();
  int nAxis1Bins = axis1->GetNbins();
  int nAxis2Bins = axis2->GetNbins();
  int nAxis3Bins = axis3->GetNbins();

  TFile* fIn_a = TFile::Open
    ( Form("%s/output_%s/c_myOut_%s.root",
	   output.c_str()  ,
	   suffix_a.c_str(),
	   suffix_a.c_str() ) );
  TFile* fIn_b = TFile::Open
    ( Form("%s/output_%s/c_myOut_%s.root",
	   output.c_str()  ,
	   suffix_b.c_str(),
	   suffix_b.c_str() ) );
  TFile* fOut  = new TFile
    ( Form("%s/all/%s/c_myOut_%s.root",
	   output.c_str()   ,
	   outSuffix.c_str(),
	   outSuffix.c_str() ) ,"recreate");

  for( int axis0Bin = 1; axis0Bin <= nAxis0Bins; axis0Bin++ ){
    double axis0Low, axis0Up;
    anaTool->GetBinRange
      ( axis0, axis0Bin, axis0Bin, axis0Low, axis0Up );
    
    for( int axis1Bin = 1; axis1Bin <= nAxis1Bins; axis1Bin++ ){
      // check we are in correct ystar and pt bins
      if( !m_dPP->CorrectPhaseSpace
	  ( std::vector<int>{ axis0Bin, axis1Bin, 0, 0 } ) )
	{ continue; }

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
	
	std::string hNameW_a = Form("h_%s_%s", name_a.c_str(), hTagW.c_str() );
	std::string hNameW_b = Form("h_%s_%s", name_b.c_str(), hTagW.c_str() );

	TH1* hW_a = static_cast<TH1D*>( fIn_a->Get( hNameW_a.c_str() ) );
	TH1* hW_b = static_cast<TH1D*>( fIn_b->Get( hNameW_b. c_str() ) );
	styleTool->SetHStyle( hW_a, style );
	styleTool->SetHStyle( hW_b, style + 5 );
	hW_a->SetMarkerSize( hW_a->GetMarkerSize() * 1.5 );
	hW_b->SetMarkerSize( hW_b->GetMarkerSize() * 1.5 );
	vHw.push_back( hW_a ); vHw.push_back( hW_b  );
	
	style++;

	if( hW_a->GetMean() ){
	  legW.AddEntry
	    ( hW_a,
	      Form("%s %s", label_a.c_str(),
		   anaTool->GetLabel( axis2Low, axis2Up, m_dPP->GetAxisLabel(2) ).c_str()));	
	  hW_a->Draw("ep same X0");
	}

	if( hW_b->GetMean() ){
	  legW.AddEntry
	    ( hW_b,
	      Form("%s %s", label_b.c_str(),
		   anaTool->GetLabel( axis2Low, axis2Up, m_dPP->GetAxisLabel(2) ).c_str()));	
	  hW_b->Draw("ep same X0");
	}

	for( int axis3Bin = 1; axis3Bin <= nAxis3Bins; axis3Bin++ ){
	  // check we are in correct ystar and pt bins
	  if( !m_dPP->CorrectPhaseSpace
	      ( std::vector<int>{ axis0Bin, axis1Bin, axis2Bin, axis3Bin } ) )
	    { continue; }

	  double axis3Low , axis3Up;
	  anaTool->GetBinRange
	    ( axis3, axis3Bin, axis3Bin, axis3Low, axis3Up );
	    
	  std::string hTag =
	    Form("%s_%s_%s_%s",
		 anaTool->GetName( axis0Low, axis0Up, m_dPP->GetAxisName(0) ).c_str(),
		 anaTool->GetName( axis1Low, axis1Up, m_dPP->GetAxisName(1) ).c_str(),
		 anaTool->GetName( axis2Low, axis2Up, m_dPP->GetAxisName(2) ).c_str(),
		 anaTool->GetName( axis3Low, axis3Up, m_dPP->GetAxisName(3) ).c_str() );
	  
	  std::string hName_a = Form("h_%s_%s", name_a.c_str(), hTag.c_str() );
	  std::string hName_b = Form("h_%s_%s", name_b.c_str(), hTag.c_str() );
	    
	  TH1* h_a = static_cast<TH1D*>( fIn_a->Get( hName_a.c_str() ) );
	  TH1* h_b = static_cast<TH1D*>( fIn_b->Get( hName_b.c_str() ) );
	  styleTool->SetHStyle( h_a, 0 );
	  styleTool->SetHStyle( h_b, 1 );
	  
	  TF1* f_a = static_cast<TF1*>( fIn_a->Get( Form("f_%s", hName_a.c_str())));
	  TF1* f_b = static_cast<TF1*>( fIn_b->Get( Form("f_%s", hName_b.c_str())));
	  styleTool->SetHStyle( f_a, 0 );
	  styleTool->SetHStyle( f_b , 1 );
	  f_a->SetLineColor( h_a->GetLineColor() );
	  f_b->SetLineColor( h_b->GetLineColor() );

	  double chi2NDF_a = f_a->GetChisquare()/f_a->GetNDF();
	  double chi2NDF_b = f_b->GetChisquare()/f_b->GetNDF();
	  
	  TCanvas c("c","c", 800, 600 );
	    
	  TLegend leg( 0.27, 0.41, 0.38, 0.52 );
	  styleTool->SetLegendStyle( &leg , 0.85 );

	  bool save = false;
	  
	  if( h_a->GetEntries() ){
	    h_a->SetMinimum(0);
	    h_a->GetYaxis()->SetNdivisions(505);
	    leg.AddEntry( h_a, Form("%s #Chi^{2}/NDF=%4.2f", label_a.c_str(), chi2NDF_a));
	    h_a->Draw("epsame");
	    f_a->Draw("same");
	    save = true;
	  }

	  if( h_b->GetEntries() ){
	    h_b->SetMinimum(0);
	    h_b->GetYaxis()->SetNdivisions(505);
	    leg.AddEntry( h_b, Form("%s #Chi^{2}/NDF=%4.2f", label_b.c_str(), chi2NDF_b));
	    h_b->Draw("epsame");
	    f_b->Draw("same");
	    save = true;
	  }

	  leg.Draw("same");

	  if( h_a->GetMaximum() > h_b->GetMaximum() ){
	    h_a->SetMaximum( h_a->GetMaximum() * 1.1 );
	    h_b->SetMaximum( h_a->GetMaximum() * 1.1 );
	  } else {
	    h_a->SetMaximum( h_b->GetMaximum() * 1.1 );
	    h_b->SetMaximum( h_b->GetMaximum() * 1.1 );
	  }

	  DrawTopLeftLabels
	    ( m_dPP, axis0Low, axis0Up, axis1Low, axis1Up,
	      axis2Low, axis2Up, axis3Low, axis3Up, 0.8 );

	  DrawAtlasRightBoth();

	  if( save ) c.SaveAs( Form("%s/all/%s/h_dPhi_%s_%s.png",
				    output.c_str(), outSuffix.c_str(),
				    hTag.c_str(), outSuffix.c_str() ) );
	  SaveAsROOT( c, Form("h_dPhi_%s", hTag.c_str() ) );
	  
	  delete  f_a; delete  f_b;
	  delete  h_a; delete  h_b;
	} // end loop over axis3
      } // end loop over axis2

      // back to cW canvas
      cW.cd();

      legW.Draw("same");

      DrawTopLeftLabels
	( m_dPP, axis0Low, axis0Up, axis1Low, axis1Up,
	  0, 0, 0, 0, 0.8 );

      DrawAtlasRightBoth();

      cW.SaveAs( Form("%s/all/%s/h_dPhi_%s_%s.png",
		      output.c_str(), outSuffix.c_str(),
		      hTagCW.c_str(), outSuffix.c_str() ) );
      SaveAsROOT( cW, Form("h_dPhi_%s", hTagCW.c_str() ) );

      for( auto & hW : vHw ){ delete hW; }
     } // end loop over ystar2
  } // end loop over ystar2
  
  std::cout << "DONE! Closing " << fOut->GetName() << std::endl;
  fOut->Close();
  std::cout << "......Closed  " << fOut->GetName() << std::endl;
}

THnSparse* DiJetAnalysis::UnfoldDeltaPhi( THnSparse* hn, TFile* fInMC,
					  const std::string& hnUnfoldedName ){
  std::vector< TH1* > vDphi;
  std::vector< TH2* > vDphiRespMat;

  std::vector< TF1* > vDphiFits;
  std::vector< TH1* > vDphiCI;
  std::vector< TH1* > vRatios;

  TAxis* axis0 = hn->GetAxis( m_dPP->GetAxisI(0) );
  TAxis* axis1 = hn->GetAxis( m_dPP->GetAxisI(1) );
  TAxis* axis2 = hn->GetAxis( m_dPP->GetAxisI(2) );
  TAxis* axis3 = hn->GetAxis( m_dPP->GetAxisI(3) );
      
  int nAxis0Bins = axis0->GetNbins();
  int nAxis1Bins = axis1->GetNbins();
  int nAxis2Bins = axis2->GetNbins();
  int nAxis3Bins = axis3->GetNbins();

  std::string unfoldingMCLabel = GetConfig()->GetValue( "unfoldingMCLabel", "" );
  
  std::string measuredName, measuredLabel;
  GetInfoUnfolding( measuredName, measuredLabel );

  // make a THnSparse to fill with unfolded results.
  THnSparse* m_hnUnfolded = (THnSparse* )hn->Clone
    ( Form("h_%s_%s", hnUnfoldedName.c_str(), m_allName.c_str() ) );
  m_hnUnfolded->Reset();
  
  // ---- loop over ystars ----
  for( int axis0Bin = 1; axis0Bin <= nAxis0Bins; axis0Bin++ ){
    // set ranges
    axis0->SetRange( axis0Bin, axis0Bin );
    double axis0Low, axis0Up;
    anaTool->GetBinRange
      ( axis0, axis0Bin, axis0Bin, axis0Low, axis0Up );
    for( int axis1Bin = 1; axis1Bin <= nAxis1Bins; axis1Bin++ ){
      // set ranges
      axis1->SetRange( axis1Bin, axis1Bin );
      double axis1Low, axis1Up;
      anaTool->GetBinRange
	( axis1, axis1Bin, axis1Bin, axis1Low, axis1Up );
      // ---- loop over axis2 ----
      for( int axis2Bin = 1; axis2Bin <= nAxis2Bins; axis2Bin++ ){
	// set ranges
	axis2->SetRange( axis2Bin, axis2Bin );
	double axis2Low , axis2Up;
	anaTool->GetBinRange
	  ( axis2, axis2Bin, axis2Bin, axis2Low, axis2Up );
	// ---- loop over axis3 ----
	for( int axis3Bin = 1; axis3Bin <= nAxis3Bins; axis3Bin++ ){
	  // check we are in correct ystar and pt bins
	  if( !m_dPP->CorrectPhaseSpace
	      ( std::vector<int>{ axis0Bin, axis1Bin, axis2Bin, axis3Bin } ) )
	    { continue; }

	  axis3->SetRange( axis3Bin, axis3Bin );
	  double axis3Low , axis3Up;
	  anaTool->GetBinRange
	    ( axis3, axis3Bin, axis3Bin, axis3Low, axis3Up );
	  
	  TCanvas c( "c", "c", 800, 600 );
	  TLegend leg( 0.22, 0.41, 0.33, 0.57 );
	  styleTool->SetLegendStyle( &leg , 0.85 );
	  
	  std::string hTag =
	    Form( "%s_%s_%s_%s",
		  anaTool->GetName( axis0Low, axis0Up, m_dPP->GetAxisName(0) ).c_str(),
		  anaTool->GetName( axis1Low, axis1Up, m_dPP->GetAxisName(1) ).c_str(),
		  anaTool->GetName( axis2Low, axis2Up, m_dPP->GetAxisName(2) ).c_str(),
		  anaTool->GetName( axis3Low, axis3Up, m_dPP->GetAxisName(3) ).c_str() ); 

	  // Take projection onto the dPhi axis
	  TH1* hDphiMeasured = hn->Projection( 4 );
	  hDphiMeasured->SetName
	    ( Form( "h_%s_%s_%s", measuredName.c_str(), m_allName.c_str(), hTag.c_str())) ;	  
	  // because variable bin width, scale by bin width
	  hDphiMeasured->Scale( 1.0, "width" );

	  // Get NON-NORMALIZED truth distribution
	  TH1* hDphiTruth    = (TH1D*)fInMC->Get
	    ( Form( "h_%s_%s_%s_NN", m_dPhiTruthName.c_str(), m_allName.c_str(), hTag.c_str()));
	  // Get NON-NORMALIZED response matrix
	  TH2* hDphiRespMat  = (TH2D*)fInMC->Get
	    ( Form( "h_%s_%s_%s", m_dPhiRespMatName.c_str(), m_allName.c_str(), hTag.c_str()));

	  vDphi.       push_back( hDphiMeasured );
	  vDphi.       push_back( hDphiTruth );
	  vDphiRespMat.push_back( hDphiRespMat );
	  
	  // setup unfold
	  RooUnfoldResponse response( hDphiMeasured, hDphiTruth, hDphiRespMat );
	  RooUnfoldBayes      unfold( &response, hDphiMeasured, 2, false );
	  // RooUnfoldBinByBin   unfold( &response, hDphiMeasured );
	  response.UseOverflow( kFALSE );
	  unfold.SetVerbose(1);

	  // unfold measured
	  TH1* hDphiUnfolded = (TH1D*)unfold.Hreco();

	  unfold.PrintTable( std::cout , hDphiTruth );
	  
	  hDphiUnfolded->SetName
	    ( Form( "h_%s_%s_%s",m_dPhiUnfoldedName.c_str(), m_allName.c_str(), hTag.c_str()));
	  hDphiUnfolded->SetTitle("");

	  // Now, after unfoldin, normalize before fit
	  NormalizeDeltaPhi( hDphiMeasured );
	  NormalizeDeltaPhi( hDphiUnfolded );
	  NormalizeDeltaPhi( hDphiTruth    );

	  hDphiMeasured->GetXaxis()->SetRangeUser( m_dPhiZoomLow, m_dPhiDphiMax );
	  hDphiUnfolded->GetXaxis()->SetRangeUser( m_dPhiZoomLow, m_dPhiDphiMax );
	  hDphiTruth   ->GetXaxis()->SetRangeUser( m_dPhiZoomLow, m_dPhiDphiMax );
	  
	  // fill unfolded THnSparse result
	  std::vector< int > x  = m_dPP->GetMappedBins
	    ( std::vector<int> { axis0Bin, axis1Bin, axis2Bin, axis3Bin } );
	  x.push_back(0); // dPhiBin;
	  for( int dPhiBin = 1; dPhiBin <= hDphiUnfolded->GetNbinsX(); dPhiBin++ ){
	    x[4] = dPhiBin;
	    m_hnUnfolded->SetBinContent
	      ( &x[0], hDphiUnfolded->GetBinContent( dPhiBin ) );
	    m_hnUnfolded->SetBinError
	      ( &x[0], hDphiUnfolded->GetBinError  ( dPhiBin ) );
	  }
	  
	  styleTool->SetHStyle( hDphiUnfolded, 0 );
	  styleTool->SetHStyle( hDphiMeasured, 1 );
	  styleTool->SetHStyle( hDphiTruth   , 2 );

	  hDphiUnfolded->Draw("ep same");
	  hDphiMeasured->Draw("ep same");
	  hDphiTruth   ->Draw("histo same");

	  double maximum = -1;
	  
	  maximum = hDphiUnfolded->GetMaximum() > hDphiMeasured->GetMaximum() ?
	    hDphiUnfolded->GetMaximum() : hDphiMeasured->GetMaximum();

	  maximum = maximum > hDphiTruth->GetMaximum() ?
	    maximum : hDphiTruth->GetMaximum();

	  // set ranges back to default
	  // otherwise the fitting doesnt work for combinatorics.
	  /*
	  hDphiMeasured->GetXaxis()->SetRangeUser( m_dPhiDphiMin, m_dPhiDphiMax );
	  hDphiUnfolded->GetXaxis()->SetRangeUser( m_dPhiDphiMin, m_dPhiDphiMax );
	  hDphiTruth   ->GetXaxis()->SetRangeUser( m_dPhiDphiMin, m_dPhiDphiMax );
	  */
	  
	  double xmin, xmax;
	  // fit and save the confidence intervals
	  // this is fit results + errors as histogram.
	  // first, fit ( no combinatoric subraction )
	  TF1* fitMeasured = anaTool->FitDphi( hDphiMeasured, false );
	  fitMeasured->GetRange( xmin, xmax );
	  TH1* hDphiMeasuredCI = static_cast<TH1D*>
	    ( hDphiMeasured->Clone( Form("%s_CI", hDphiMeasured->GetName() ) ) );
	  if( hDphiMeasuredCI->GetEntries() )
	    { (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hDphiMeasuredCI); }
	  
	  // fit with no combinatoric subtraction (already done);
	  TF1* fitUnfolded = anaTool->FitDphi( hDphiUnfolded, false );
	  TH1* hDphiUnfoldedCI = static_cast<TH1D*>
	    ( hDphiUnfolded->Clone( Form("%s_CI", hDphiUnfolded->GetName() ) ) );
	  if( hDphiUnfoldedCI->GetEntries() )
	    { (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hDphiUnfoldedCI); }

	  vDphiCI.push_back( hDphiMeasuredCI );
	  vDphiCI.push_back( hDphiUnfoldedCI );
	  
	  // fit with no combinatoric subtraction (already done);
	  TF1* fitTruth    = anaTool->FitDphi( hDphiTruth, false );
	  
	  styleTool->SetHStyle( fitMeasured, 0 );
	  styleTool->SetHStyle( fitUnfolded, 0 );
	  styleTool->SetHStyle( fitTruth, 0 );

	  vDphiFits.push_back( fitMeasured );
	  vDphiFits.push_back( fitUnfolded );
	  vDphiFits.push_back( fitTruth );

	  fitMeasured->SetLineColor( hDphiMeasured->GetLineColor() );
	  fitUnfolded->SetLineColor( hDphiUnfolded->GetLineColor() );
	  fitTruth   ->SetLineColor( hDphiTruth   ->GetLineColor() );
	  
	  fitMeasured->Draw("same");
	  fitUnfolded->Draw("same");
	  fitTruth   ->Draw("same");

	  double chi2NDFmeasured = fitMeasured->GetChisquare()/fitMeasured->GetNDF();
	  double chi2NDFunfolded = fitUnfolded->GetChisquare()/fitUnfolded->GetNDF();
	  double chi2NDFtruth    = fitTruth   ->GetChisquare()/fitTruth   ->GetNDF();

	  hDphiUnfolded->SetMaximum( maximum * 1.1 );
	  hDphiUnfolded->SetMaximum( maximum * 1.1 );
	  hDphiTruth   ->SetMaximum( maximum * 1.1 );

	  leg.AddEntry( hDphiTruth,
			Form( "Truth #Chi^{2}/NDF=%4.2f",
				   chi2NDFtruth ), "lf" );
	  leg.AddEntry( hDphiMeasured,
			Form("%s #Chi^{2}/NDF=%4.2f",
			     measuredLabel.c_str(), chi2NDFmeasured ) );
	  leg.AddEntry( hDphiUnfolded,
			Form("Unfolded #Chi^{2}/NDF=%4.2f",
			     chi2NDFunfolded ) );
	  
	  leg.Draw();

	  drawTool->DrawLeftLatex( 0.25, 0.39, unfoldingMCLabel.c_str() );
	  
	  DrawTopLeftLabels
	    ( m_dPP, axis0Low, axis0Up, axis1Low, axis1Up,
	      axis2Low, axis2Up, axis3Low, axis3Up, 0.8 );
	    
	  DrawAtlasRight();
	    
	  SaveAsAll( c, Form( "h_%s_%s_MUT_%s",
			       m_dPhiUnfoldedName.c_str(),
			       m_allName.c_str(),
			       hTag.c_str()));

	  TCanvas cR("cR","cR", 800, 600 );

	  TH1* hR = static_cast< TH1D* >
	    ( hDphiUnfolded->Clone( Form( "h_%s_%s_ratio_%s",
					  m_dPhiUnfoldedName.c_str(),
					  m_allName.c_str(),
					  hTag.c_str())));
	  
	  hR->Divide( hDphiTruth );
	  styleTool->SetHStyle( hR, 0 );
	  
	  hR->SetStats(kFALSE);
	  hR->SetFillColor(46);
	  hR->GetYaxis()->SetTitle("|#Delta#phi| Unfolded/Truth");
	  hR->GetXaxis()->SetTitle("|#Delta#phi|");
	  hR->Draw("e2p");
	  hR->SetMaximum( 2.0 );
	  hR->SetMinimum( 0.0 );

	  double xMin = hR->GetXaxis()->GetXmin();
	  double xMax = hR->GetXaxis()->GetXmax();
	  
	  TLine line( xMin, 1, xMax, 1 );
	  line.Draw();

	  DrawTopLeftLabels
	    ( m_dPP, axis0Low, axis0Up, axis1Low, axis1Up,
	      axis2Low, axis2Up, axis3Low, axis3Up, 0.8 );
	    
	  DrawAtlasRight();
	  
	  SaveAsROOT( cR, hR->GetName() );	  
	} // end loop over axis3
      } // end loop over axis2
    } // end loop over axis1     
  } // end loop over axis0
  for( auto f  : vDphiFits    ){ delete f;  }
  for( auto h  : vDphiCI      ){ delete h;  }

  for( auto h  : vDphi        ){ delete h ; }
  for( auto rm : vDphiRespMat ){ delete rm; }

  for( auto r  : vRatios      ){ delete r;  }
  
  return m_hnUnfolded;
}


//---------------------------
//        Drawing
//---------------------------

void DiJetAnalysis::DrawAtlasRight( double x0, double y0, double scale ){}

void DiJetAnalysis::DrawAtlasRightBoth( double x0, double y0, double scale ){} 

void DiJetAnalysis::DrawTopLeftLabels( DeltaPhiProj* dPP,
				       double axis0Low, double axis0Up,
				       double axis1Low, double axis1Up,
				       double axis2Low, double axis2Up,
				       double axis3Low, double axis3Up,
				       double scale ){
  drawTool->DrawLeftLatex
    ( 0.13, 0.86, CT::AnalysisTools::GetLabel
      ( axis0Low, axis0Up, dPP->GetAxisLabel(0) ), scale );
  // for now, modify later to be more dynamic
  // to the variables present 
  if( !( axis1Low || axis1Up ) ){ return ; }
  drawTool->DrawLeftLatex
    ( 0.13, 0.79, CT::AnalysisTools::GetLabel
      ( axis1Low, axis1Up, dPP->GetAxisLabel(1) ), scale );  
  if( !( axis2Low || axis2Up ) ){ return ; }
  drawTool->DrawLeftLatex
    ( 0.13, 0.73, CT::AnalysisTools::GetLabel
      ( axis2Low, axis2Up, dPP->GetAxisLabel(2) ), scale );
  if( !( axis3Low || axis3Up ) ){ return ; }
  drawTool->DrawLeftLatex
    ( 0.13, 0.66, CT::AnalysisTools::GetLabel
      ( axis3Low, axis3Up, dPP->GetAxisLabel(3) ), scale );
}

//---------------------------
//       Saving 
//---------------------------

void DiJetAnalysis::SaveAsPdfPng( const TCanvas& c,
				  const std::string& label1,
				  const std::string& label2,
				  const std::string& axis1,
				  double min1, double max1,
				  const std::string& axis2 ,
				  double min2, double max2 ){
  std::stringstream ss;
  ss << m_dirOut + "/";
  
  if( !label1.empty() ){ ss << boost::format("%s") % label1; }
  if( !label2.empty() ){ ss << boost::format("_%s") % label2; }
  ss << "_" + m_labelOut;
  if( !axis1.empty() ){
    ss << boost::format("_%2.0f_%s_%2.0f") % min1 % axis1 % max1;
  }
  if( !axis2.empty() ){
    ss << boost::format("_%2.0f_%s_%2.0f") % min2 % axis2 % max2;
  }

  std::string sPdf = ss.str() + ".pdf";
  std::string sPng = ss.str() + ".png";

  c.SaveAs( sPdf.c_str() );
  c.SaveAs( sPng.c_str() );
}

void DiJetAnalysis::SaveAsROOT( const TCanvas& c,
			        const std::string& label1,
				const std::string& label2,
				const std::string& axis1,
				double min1, double max1,
				const std::string& axis2 ,
				double min2, double max2 ){
  std::stringstream ss;
  ss << "c_";
  
  if( !label1.empty() ){ ss << boost::format("%s") % label1; }
  if( !label2.empty() ){ ss << boost::format("_%s") % label2; }
  ss << "_" + m_labelOut;
  if( !axis1.empty() ){
    ss << boost::format("_%2.0f_%s_%2.0f") % min1 % axis1 % max1;
  }
  if( !axis2.empty() ){
    ss << boost::format("_%2.0f_%s_%2.0f") % min2 % axis2 % max2;
  }

  c.Write( ss.str().c_str() );
}

void DiJetAnalysis::SaveAsAll( const TCanvas& c,
			       const std::string& label1,
			       const std::string& label2,
			       const std::string& axis1,
			       double min1, double max1,
			       const std::string& axis2 ,
			       double min2, double max2 ){  
  SaveAsPdfPng( c, label1, label2, axis1, min1, max1, axis2, min2, max2 );
  SaveAsROOT  ( c, label1, label2, axis1, min1, max1, axis2, min2, max2 );
}
