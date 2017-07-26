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
#include <TGaxis.h>
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

DiJetAnalysis::DiJetAnalysis() : DiJetAnalysis( false, false, 0, 0 ){}

DiJetAnalysis::DiJetAnalysis( bool is_pPb )
  : DiJetAnalysis( is_pPb, false, 0, 0 ){}

DiJetAnalysis::DiJetAnalysis( bool is_pPb, bool isData )
  : DiJetAnalysis( is_pPb, isData, 0, 0 ){}

DiJetAnalysis::DiJetAnalysis( bool is_pPb, bool isData, int mcType )
  : DiJetAnalysis( is_pPb, isData, mcType, 0 ){}


/* Constructor for DiJetAnalysis
 * 
 * Here, make histogram binning, load
 * config file. Also, define some 
 * common used names.
 * 
 * Rest of things happen in Initialize()
 * such as reading options from config file
 * and further setting names, etc.
 */ 
DiJetAnalysis::DiJetAnalysis( bool is_pPb, bool isData, int mcType, int uncertComp )
  : m_is_pPb( is_pPb ), m_isData( isData ), m_mcType(mcType), m_uncertComp( uncertComp ){
  
  //============= common strings =================
  m_s_pp  = "pp";
  m_s_pPb = "pPb";
  m_s_pt  = "pT";

  m_s_pt1  = "pT1";
  m_s_pt2  = "pT2";

  m_sOutput   = "output";  
  m_myOutName = "myOut";
  m_sMC       = "mc";
  m_sData     = "data";
  
  m_sMUT   = "MUT";
  m_sRatio = "ratio";

  // name for "All" histos
  // this is either merged JZN or Data from Triggers
  m_allName = "All";

  m_unfoldingFileSuffix = "UF";
  
  boost::assign::push_back( m_vMCtypeNames )
    ( "pythia8" )( "herwig" )( "pythia8powheg" );

  boost::assign::push_back( m_vMCtypeLabels  )
    ( "Pythia8" )( "Herwig" )( "Pythia8+Powheg" );

  //==================== Cuts ====================
  m_nMinEntriesFit = 20;
    
  m_dPhiThirdJetFraction = 0.4;

  m_dPhiUnfoldingMin = 2.0;
  m_dPhiUnfoldingMax = constants::PI; 
  
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
  m_etaForwardMin = -3.9;
  m_etaForwardMax = -2.7;
  m_nEtaForwardBinsFine = 6;
  
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
    ( 28 )( 35 )( 45 )( 90 );
  m_nVarPtBins = m_varPtBinning.size() - 1;

  m_nDphiBinsLarge  = 8 ; m_dPhiBinsLargeFactor  = 8;
  m_nDphiBinsMedium = 10; m_dPhiBinsMediumFactor = 2;
  m_nDphiBinsSmall  = 4 ; m_dPhiBinsSmallFactor  = 1;

  // --- variable dPhi binning ---
  double smallBinWidth = constants::PI /
    ( m_nDphiBinsSmall  * m_dPhiBinsSmallFactor  +
      m_nDphiBinsMedium * m_dPhiBinsMediumFactor +
      m_nDphiBinsLarge  * m_dPhiBinsLargeFactor );

  for( int largeBin = 0; largeBin < m_nDphiBinsLarge; largeBin++ )
    { m_varDphiBinning.push_back
	( largeBin * smallBinWidth * m_dPhiBinsLargeFactor ); }
  double mediumBinStart =
    m_nDphiBinsLarge * smallBinWidth * m_dPhiBinsLargeFactor;
  for( int mediumBin = 0; mediumBin < m_nDphiBinsMedium; mediumBin++ )
    { m_varDphiBinning.push_back
	( mediumBin * smallBinWidth * m_dPhiBinsMediumFactor + mediumBinStart ); }
  double smallBinStart =
    m_nDphiBinsLarge  * smallBinWidth * m_dPhiBinsLargeFactor +
    m_nDphiBinsMedium * smallBinWidth * m_dPhiBinsMediumFactor;
  for( int smallBin = 0; smallBin <= m_nDphiBinsSmall; smallBin++ )
    { m_varDphiBinning.push_back
	( smallBin * smallBinWidth * m_dPhiBinsSmallFactor + smallBinStart ); }
  
  m_nVarDphiBins = m_varDphiBinning.size() - 1;
  
  // --- dPhiBins ---  
  boost::assign::push_back( m_nDphiBins    )
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

  // default values 0->2 are in
  // SetHStyleRatio.
  // This is for "nonstandard"
  m_ratioMin = 0.0;
  m_ratioMax = 2.0;
  
  
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
  
  //========= common tools ========
  anaTool   = new CT::AnalysisTools();
  drawTool  = new CT::DrawTools();
  styleTool = new CT::StyleTools();

  //========= config file =========
  std::string configName   = "config/configJetsFwd";
  
  std::string configDataType = m_is_pPb ? "_" + m_s_pPb : "_" + m_s_pp;
  configName += configDataType;

  std::string configSuffix = m_isData ? "_" + m_sData + ".cfg" : "_" + m_sMC + ".cfg" ;
  configName += configSuffix;

  std::cout << "Reading " << configName << std::endl;

  m_config = new TEnv();
  m_config->ReadFile( configName.c_str(), EEnvLevel(0));
  
  //=============== Histo Names ==================    
  m_etaSpectName      = "etaSpect";
  m_dPhiName          = "dPhi";
  m_effName           = "eff";
  m_purityName        = "purity";

  m_recoName          = "reco";
  m_truthName         = "truth";
  m_respMatName       = "respMat";
  m_unfoldedName      = "unfolded";
    
  m_dPhiRecoName      = m_dPhiName + "_" + m_recoName;
  m_dPhiTruthName     = m_dPhiName + "_" + m_truthName;
  
  m_dPhiRespMatName   = m_dPhiName + "_" + m_respMatName;
  m_ptRespMatName     = m_s_pt     + "_" + m_respMatName;

  m_allRespMatName    = m_allName  + "_" + m_respMatName;
  
  m_dPhiCFactorsName  =
    m_dPhiName + "_" + m_recoName + "_" + m_truthName + "_" + m_sRatio;
  m_dPhiUnfoldedName     = m_dPhiName     + "_" + m_unfoldedName;
  m_dPhiRecoUnfoldedName = m_dPhiRecoName + "_" + m_unfoldedName;
}

DiJetAnalysis::~DiJetAnalysis(){
  
  delete m_dPP    ; m_dPP     = NULL;
  delete anaTool  ; anaTool   = NULL;
  delete drawTool ; drawTool  = NULL;
  delete styleTool; styleTool = NULL;
}

void DiJetAnalysis::Initialize(){
  
  m_mcTypeName  = m_vMCtypeNames [ m_mcType ];
  m_mcTypeLabel = m_vMCtypeLabels[ m_mcType ];

  std::string system = m_is_pPb ? m_s_pPb : m_s_pp;

  m_labelOut = m_isData ? system + "_" + m_sData : system + "_" + m_sMC;

  // append a suffix to. data and MC can have
  // different m_labelOut
  AdditionalSuffix( m_labelOut );
  
  // Check if the directories exist.
  // If they don't, create them
  m_dirOut   = m_sOutput;
  anaTool->CheckWriteDir( m_dirOut.c_str() );
  m_dirOut   += "/" + m_sOutput + "_" + m_labelOut;
  anaTool->CheckWriteDir( m_dirOut.c_str() );

  // need suffix for uncertainties
  std::string ucSuffix =
    m_uncertComp ? ( m_uncertComp > 0 ?
		     Form("_P%d", m_uncertComp) :
		     Form("_N%d", -1 * m_uncertComp) ) : "" ;
  
  // this is the fname of the raw histograms. I.e.
  // what is produced from reading the tree (THnMulfs, TH3, TH2 )
  // These get used to take projections, etc, later.
  m_rawHistosFname =
    m_dirOut + "/" + m_myOutName +"_" + m_labelOut + ucSuffix + ".root";
  
  std::cout << "fNameIn/Out: " << m_rawHistosFname << std::endl;
  
  std::string unfoldingMC = m_mcTypeName;

  // set names for various output files.
  // for data only, we just have one.
  // no need to add unfolding suffix.
  m_fNameOutDefault  = m_dirOut + "/c_" + m_myOutName + "_" + m_labelOut;
  m_fNameOut         = m_fNameOutDefault + ucSuffix;
  m_fNameOutUF       = m_fNameOut + "_" + m_unfoldingFileSuffix;
  m_fNameOutDefault += ".root";
  m_fNameOut        += ".root";
  m_fNameOutUF      += ".root";

  // name of file that is used for unfolding
  m_fNameUnfoldingMC = Form( "%s/%s_%s_%s_%s/c_%s_%s_%s_%s.root",
			     m_sOutput.c_str(), m_sOutput.c_str(),
			     system.c_str(), m_sMC.c_str(),
			     unfoldingMC.c_str(), m_myOutName.c_str(),
			     system.c_str(), m_sMC.c_str(),
			     unfoldingMC.c_str() );

}

//---------------------------
//       Fill Tree
//---------------------------
void DiJetAnalysis::AddHistogram( TH1* h ){
  
  v_hists.push_back( h );
  h->Sumw2();
  styleTool->SetHStyle( h, 0 );

  TH1* h1 = dynamic_cast< TH1* >(h);
  if( h1 ){ h->SetNdivisions(505, "X" ); }

  TH2* h2 = dynamic_cast< TH2* >(h);
  if( h2 ){
    h->SetNdivisions(505, "XY" );
  }

  TH3* h3 = dynamic_cast< TH3* >(h);
  if( h3 ){
    h->GetXaxis()->SetNdivisions( 505 , "XYZ" );
  }
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
}


void DiJetAnalysis::SaveOutputsFromTree(){
  
  //----------------------------------------
  //  Close the input file, 
  //  write histos to output
  //----------------------------------------
  std::cout << "fNameOut: " << m_rawHistosFname << std::endl;
  TFile* fOut = new TFile( m_rawHistosFname.c_str(),"RECREATE");
  for( auto& h  : v_hists  ){ h-> Write(); }
  for( auto& f  : v_functs ){ f-> Write(); }
  for( auto& gr : v_graphs ){ gr->Write(); }
  for( auto& hn : v_hns    ){ hn->Write(); }
  fOut->Close();
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

double DiJetAnalysis::AnalyzeDeltaPhi( THnSparse* hn,
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
  hn->Fill( &x[0], weight );
  
  // for pp, fill twice. once for each side since
  // it is symmetric in pp. For pPb, continue
  if( m_is_pPb ){ return deltaPhi; }
  
  x[0] = -jet1_ystar;  
  x[1] = -jet2_ystar;
  hn->Fill( &x[0], weight );
  
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
				      std::string& measuredLabel,
				      std::string& typeLabel ){}

// hM is measured, hC is correction factor
TH1* DiJetAnalysis::BinByBinUnfolding( TH1* hM, TH1* hC ){
  
  TH1* hUnf = static_cast<TH1D*>(hM->Clone("h_unfolded"));
  hUnf->Reset();
  
  for( int xBin = hM->FindBin( m_dPhiUnfoldingMin ); xBin <= hC->GetNbinsX(); xBin++ ){
    double vR = hM->GetBinContent( xBin  );
    double vC = hC->GetBinContent( xBin  );

    double eR = hM->GetBinError(xBin);
    double eC = hC->GetBinError(xBin);
    // correction factor;
    double newDphi = vR * vC;
    // error on correction factor    
    double newDphiError =  newDphi * 
      std::sqrt( std::pow( eR / vR, 2) +
		 std::pow( eC / vC, 2) ) ; 
     
    hUnf->SetBinContent( xBin, newDphi      );
    hUnf->SetBinError  ( xBin, newDphiError );
  }
  
  return hUnf;
}

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
      double etaMin, etaMax;
      anaTool->GetBinRange
	( hRef->GetXaxis(), xBin, xBin, etaMin, etaMax );
      
      TH1* h_etaSpect =
	vSampleSpect[iG]->
	ProjectionY( Form("h_%s_%s_%s",
			  name.c_str(),
			  label.c_str(),
			  anaTool->GetName( etaMin, etaMax, "Eta").c_str() ),
		     xBin, xBin );
      h_etaSpect->SetTitle
	( anaTool->GetEtaLabel( etaMin, etaMax, m_is_pPb ).c_str() );
      h_etaSpect->SetYTitle( yAxisTitle.c_str() );
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
    double etaMin, etaMax;
    anaTool->GetBinRange
      ( hRef->GetXaxis(), xBin, xBin, etaMin, etaMax );
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
				  const std::string& name,
				  bool isUnfolded ){

  std::vector< TH1* > vDphi;
  std::vector< TF1* > vFits;

  // ---- loop over group  ----
  // ---- (jzn or trigger) ----
  for( uint iG = 0; iG < vhn.size(); iG++ ){      
    THnSparse* hn = vhn[iG];
  
    std::string label = vLabel[iG];

    if( label.compare( m_allName ) ){ continue; } 
  
    TAxis* axis0 = hn->GetAxis( m_dPP->GetAxisI(0) ); int nAxis0Bins = axis0->GetNbins();
    TAxis* axis1 = hn->GetAxis( m_dPP->GetAxisI(1) ); int nAxis1Bins = axis1->GetNbins();
    TAxis* axis2 = hn->GetAxis( m_dPP->GetAxisI(2) ); int nAxis2Bins = axis2->GetNbins();
    TAxis* axis3 = hn->GetAxis( m_dPP->GetAxisI(3) ); int nAxis3Bins = axis3->GetNbins();
    
    int fAxisI   = m_dPP->GetAxisI(3);
 
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

	std::string hTagCW =
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
	  hDphiWidths->SetYTitle( "|#Delta#phi| width" );
	  hDphiWidths->SetMarkerSize( hDphiWidths->GetMarkerSize() * 1.5 );
	  hDphiWidths->SetNdivisions( 505, "X" );
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
	    std::string hDphiName = Form( "h_%s_%s_%s", name.c_str(), label.c_str(), hTag.c_str() );
	    TH1* hDphi = hn->Projection( 4 );
	    hDphi->SetName( hDphiName.c_str() );
	    styleTool->SetHStyle( hDphi, 0 );
	    hDphi->SetNdivisions( 505, "Y" );
	    vDphi.push_back( hDphi );

	    // if its not unfolded result, subtract combinatoric, noramlize
	    if( !isUnfolded ){
	      // because variable bin width, scale by bin width
	      hDphi->Scale( 1.0, "width" );
	      // subtract combinatoric contribution before normalizing
	      anaTool->SubtractCombinatoric( hDphi ); 
	      // Normalize
	      NormalizeDeltaPhi( hDphi );
	    }
	    
	    TCanvas c( "c", hDphi->GetName(), 800, 600 );
	    
	    hDphi->Draw();
	    hDphi->SetYTitle("Normalized Count");
	    hDphi->SetTitle("");
    
	    // now fit
	    TF1* fit = anaTool->FitDphi( hDphi, m_dPhiUnfoldingMin, m_dPhiUnfoldingMax );
	    styleTool->SetHStyle( fit, 0 );
	    vFits.push_back( fit );

	    // !!!!!!!
	    std::cout << fit->GetName() << "  " << fit->GetParameter(1) << std::endl;
	  	    
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
	    // if( chi2NDF < 0.5 || chi2NDF > 2.5 ){ continue; }
	    hDphiWidths->SetBinContent( axis3Bin, fit->GetParameter(1) );
	    hDphiWidths->SetBinError  ( axis3Bin, fit->GetParError (1) );	    
	  } // end loop over axis3
	} // end loop over axis2
	
	cWidths.cd();
	for( auto& h : vDphiWidthsTemp ){
	  h->SetMinimum( m_dPhiWidthMin );
	  h->SetMaximum( m_dPhiWidthMax );
	  h->SetNdivisions( 505, "Y" );
	  h->SetTitle("");
	  h->Draw("epsame");
	  h->Write();
	}
	leg.Draw("same");

	DrawTopLeftLabels
	  ( m_dPP, axis0Low, axis0Up, axis1Low, axis1Up,
	    0, 0, 0, 0, 0.8 );

	DrawAtlasRight();
	
        SaveAsAll( cWidths, Form("h_%s_%s_%s", name.c_str(), label.c_str(), hTagCW.c_str() ) );
      } // end loop over axis1     
    } // end loop over axis0
  } // end loop over iG

  for( auto& f : vFits ){ delete f; }
  for( auto& h : vDphi ){ delete h; }
}

THnSparse* DiJetAnalysis::UnfoldDeltaPhi( TFile* fInData, TFile* fInMC,
					  const std::string& hnUnfoldedName ){
  std::vector< TH1* > vDphi;
  std::vector< TF1* > vDphiFits;
  
  std::vector< TH1* > vCFactors;
  std::vector< TH2* > vDphiRespMat;
  
  TAxis* axis0 = m_dPP->GetTAxis(0); int nAxis0Bins = axis0->GetNbins();
  TAxis* axis1 = m_dPP->GetTAxis(1); int nAxis1Bins = axis1->GetNbins();
  TAxis* axis2 = m_dPP->GetTAxis(2); int nAxis2Bins = axis2->GetNbins();
  TAxis* axis3 = m_dPP->GetTAxis(3); int nAxis3Bins = axis3->GetNbins();

  std::string unfoldingMCLabel = m_mcTypeLabel;
  
  std::string measuredName, measuredLabel, typeMeasured;
  GetInfoUnfolding( measuredName, measuredLabel, typeMeasured );

  // make a THnSparse to fill with unfolded results.
  THnSparse* hnUnfolded =
    new THnSparseD( Form("h_%s_%s", hnUnfoldedName.c_str(), m_allName.c_str() ) , "",
		    m_nDphiDim, &m_nDphiBins[0],
		    &m_dPhiMin[0], &m_dPhiMax[0] );

  TAxis* axis0Def = m_dPP->GetDefaultTAxis( 0 );
  TAxis* axis1Def = m_dPP->GetDefaultTAxis( 1 );
  TAxis* axis2Def = m_dPP->GetDefaultTAxis( 2 );
  TAxis* axis3Def = m_dPP->GetDefaultTAxis( 3 );
  
  hnUnfolded->GetAxis(0)->
    Set( axis0Def->GetNbins(), axis0Def->GetXbins()->GetArray() );
  hnUnfolded->GetAxis(1)->
    Set( axis1Def->GetNbins(), axis1Def->GetXbins()->GetArray() );
  hnUnfolded->GetAxis(2)->
    Set( axis2Def->GetNbins(), axis2Def->GetXbins()->GetArray() );
  hnUnfolded->GetAxis(3)->
    Set( axis3Def->GetNbins(), axis3Def->GetXbins()->GetArray() );
  hnUnfolded->GetAxis(4)->Set( m_nVarDphiBins,   &( m_varDphiBinning[0]  ) );

  hnUnfolded->GetAxis(0)->SetTitle( m_dPP->GetDefaultAxisLabel(0).c_str() );
  hnUnfolded->GetAxis(1)->SetTitle( m_dPP->GetDefaultAxisLabel(1).c_str() );
  hnUnfolded->GetAxis(2)->SetTitle( m_dPP->GetDefaultAxisLabel(2).c_str() );
  hnUnfolded->GetAxis(3)->SetTitle( m_dPP->GetDefaultAxisLabel(3).c_str() );
  hnUnfolded->GetAxis(4)->SetTitle( measuredLabel.c_str() );
  
  // ---- loop over ystars ----
  for( int axis0Bin = 1; axis0Bin <= nAxis0Bins; axis0Bin++ ){
    double axis0Low, axis0Up;
    anaTool->GetBinRange
      ( axis0, axis0Bin, axis0Bin, axis0Low, axis0Up );
    for( int axis1Bin = 1; axis1Bin <= nAxis1Bins; axis1Bin++ ){
      double axis1Low, axis1Up;
      anaTool->GetBinRange
	( axis1, axis1Bin, axis1Bin, axis1Low, axis1Up );
      // ---- loop over axis2 ----
      for( int axis2Bin = 1; axis2Bin <= nAxis2Bins; axis2Bin++ ){
	double axis2Low , axis2Up;
	anaTool->GetBinRange
	  ( axis2, axis2Bin, axis2Bin, axis2Low, axis2Up );
	// ---- loop over axis3 ----
	for( int axis3Bin = 1; axis3Bin <= nAxis3Bins; axis3Bin++ ){
	  // check we are in correct ystar and pt bins
	  if( !m_dPP->CorrectPhaseSpace
	      ( std::vector<int>{ axis0Bin, axis1Bin, axis2Bin, axis3Bin } ) )
	    { continue; }

	  double axis3Low , axis3Up;
	  anaTool->GetBinRange
	    ( axis3, axis3Bin, axis3Bin, axis3Low, axis3Up );
	  	  
	  TLegend leg( 0.22, 0.41, 0.33, 0.57 );
	  styleTool->SetLegendStyle( &leg , 0.85 );
	  
	  std::string hTag =
	    Form( "%s_%s_%s_%s",
		  anaTool->GetName( axis0Low, axis0Up, m_dPP->GetAxisName(0) ).c_str(),
		  anaTool->GetName( axis1Low, axis1Up, m_dPP->GetAxisName(1) ).c_str(),
		  anaTool->GetName( axis2Low, axis2Up, m_dPP->GetAxisName(2) ).c_str(),
		  anaTool->GetName( axis3Low, axis3Up, m_dPP->GetAxisName(3) ).c_str() ); 

	  // Get measured distribution
	  TH1* hDphiMeasured = static_cast<TH1D*>
	    ( fInData->Get( Form( "h_%s_%s_%s", measuredName.c_str(), m_allName.c_str(), hTag.c_str())));
	  // Get truth distribution

	  TH1* hDphiTruth    = static_cast<TH1D*>
	    (fInMC->Get( Form( "h_%s_%s_%s", m_dPhiTruthName.c_str(), m_allName.c_str(), hTag.c_str())));
	  TH1* hCFactors     = static_cast<TH1D*>
	    (fInMC->Get( Form( "h_%s_%s_%s", m_dPhiCFactorsName.c_str(), m_allName.c_str(), hTag.c_str())));
	  TH2* hDphiRespMat  = static_cast<TH2D*>
	    (fInMC->Get( Form( "h_%s_%s_%s", m_dPhiRespMatName.c_str(), m_allName.c_str(), hTag.c_str())));

	  vDphi.push_back( hDphiMeasured );
	  vDphi.push_back( hDphiTruth    );

	  vCFactors.push_back( hCFactors     );
	  vDphiRespMat.push_back( hDphiRespMat  );
	  
       	  TF1* fitMeasured =  static_cast<TF1*>
	    ( fInData->Get( Form( "f_h_%s_%s_%s", measuredName.c_str(), m_allName.c_str(), hTag.c_str())));
	  TF1* fitTruth    =  static_cast<TF1*>
	    ( fInMC->Get( Form( "f_h_%s_%s_%s", m_dPhiTruthName.c_str(), m_allName.c_str(), hTag.c_str())));
	  
	  // ----------- Unfold -----------
	  // Unfold using bin-by-bin and the resposne factors.
	  TH1* hDphiUnfolded = BinByBinUnfolding
	    ( hDphiMeasured,  hCFactors );
	  hDphiUnfolded->SetName
	    ( Form( "h_%s_%s_%s",m_dPhiUnfoldedName.c_str(), m_allName.c_str(), hTag.c_str()));
	  
	  // fit with no combinatoric subtraction (already done);
	  TF1* fitUnfolded = anaTool->FitDphi( hDphiUnfolded, m_dPhiUnfoldingMin, m_dPhiUnfoldingMax );

	  // !!!!!!!!!
	  std::cout << fitUnfolded->GetName() << "   " << fitUnfolded->GetParameter(1) << std::endl;
	  
	  // -------- Unfold Done ---------
	  
	  // fill unfolded THnSparse result
	  std::vector< int > x  = m_dPP->GetMappedBins
	    ( std::vector<int> { axis0Bin, axis1Bin, axis2Bin, axis3Bin } );
	  x.push_back(0); // dPhiBin;
	  for( int dPhiBin = 1; dPhiBin <= hDphiUnfolded->GetNbinsX(); dPhiBin++ ){
	    x[4] = dPhiBin;
	    hnUnfolded->SetBinContent
	      ( &x[0], hDphiUnfolded->GetBinContent( dPhiBin ) );
	    hnUnfolded->SetBinError
	      ( &x[0], hDphiUnfolded->GetBinError  ( dPhiBin ) );
	  }
	  
	  styleTool->SetHStyle( hDphiUnfolded, 0 );
	  styleTool->SetHStyle( hDphiMeasured, 1 );
	  styleTool->SetHStyle( hDphiTruth   , 2 );
	  styleTool->SetHStyleRatio( hCFactors );
	  
	  hDphiMeasured->GetXaxis()->SetRangeUser( m_dPhiZoomLow, m_dPhiDphiMax );
	  hDphiUnfolded->GetXaxis()->SetRangeUser( m_dPhiZoomLow, m_dPhiDphiMax );
	  hDphiTruth   ->GetXaxis()->SetRangeUser( m_dPhiZoomLow, m_dPhiDphiMax );
	  
	  // Now Draw everything.
	  TCanvas c( "c", "c", 800, 700 );
	  TPad pad1("pad1", "", 0.0, 0.25, 1.0, 1.0 );
	  pad1.SetBottomMargin(0);
	  pad1.Draw();
	  pad1.cd();

	  hDphiMeasured->Draw("ep same");	
	  hDphiTruth   ->Draw("histo same");
	  hDphiUnfolded->Draw("ep same");
	  
	  styleTool->SetHStyle( fitMeasured, 0 );
	  styleTool->SetHStyle( fitUnfolded, 0 );
	  styleTool->SetHStyle( fitTruth   , 0 );

	  vDphiFits.push_back( fitMeasured );
	  vDphiFits.push_back( fitUnfolded );
	  vDphiFits.push_back( fitTruth    );

	  fitMeasured->SetLineColor( hDphiMeasured->GetLineColor() );
	  fitUnfolded->SetLineColor( hDphiUnfolded->GetLineColor() );
	  fitTruth   ->SetLineColor( hDphiTruth   ->GetLineColor() );

	  fitMeasured->SetLineStyle( 3 );
	  fitUnfolded->SetLineStyle( 3 );
	  fitTruth   ->SetLineStyle( 3 );

	  fitMeasured->Draw("same");
	  fitUnfolded->Draw("same");
 	  fitTruth   ->Draw("same");
	  
	  double chi2NDFmeasured = fitMeasured->GetChisquare()/fitMeasured->GetNDF();
	  double chi2NDFunfolded = fitUnfolded->GetChisquare()/fitUnfolded->GetNDF();
	  double chi2NDFtruth    = fitTruth   ->GetChisquare()/fitTruth   ->GetNDF();

	  double maximum = -1;
	  
	  maximum = hDphiUnfolded->GetMaximum() > hDphiMeasured->GetMaximum() ?
	    hDphiUnfolded->GetMaximum() : hDphiMeasured->GetMaximum();
	  maximum = maximum > hDphiTruth->GetMaximum() ?
	    maximum : hDphiTruth->GetMaximum();
	  
	  hDphiUnfolded->SetMaximum( maximum * 1.1 );
	  hDphiUnfolded->SetMaximum( maximum * 1.1 );
	  hDphiTruth   ->SetMaximum( maximum * 1.1 );

	  hDphiUnfolded->SetMinimum( 0 );
	  hDphiUnfolded->SetMinimum( 0 );
	  hDphiTruth   ->SetMinimum( 0 );

	  leg.AddEntry( hDphiTruth,
	 		Form("Truth #Chi^{2}/NDF=%4.2f",
			     chi2NDFtruth ), "lf" );
	  leg.AddEntry( hDphiMeasured,
			Form("%s #Chi^{2}/NDF=%4.2f",
			     typeMeasured.c_str(), chi2NDFmeasured ) );
	  leg.AddEntry( hDphiUnfolded,
			Form("Unfolded #Chi^{2}/NDF=%4.2f",
			     chi2NDFunfolded ) );
	  
	  leg.Draw();

	  drawTool->DrawLeftLatex( 0.25, 0.39, unfoldingMCLabel.c_str() );
	  
	  DrawTopLeftLabels
	    ( m_dPP, axis0Low, axis0Up, axis1Low, axis1Up,
	      axis2Low, axis2Up, axis3Low, axis3Up, 0.8 );
	    
	  DrawAtlasRight();

	  c.cd();
	  TPad pad2("pad2", "", 0.0, 0.0, 1.0, 0.25 );
	  pad2.SetTopMargin(0);
	  pad2.SetBottomMargin(0.2);
	  pad2.Draw();
	  pad2.cd();

	  hCFactors->SetYTitle("|#Delta#phi_{Truth}|/|#Delta#phi_{Reco}|");
	  hCFactors->GetXaxis()->SetRangeUser( m_dPhiZoomLow, m_dPhiDphiMax );
	  
	  hCFactors->Draw("e2p");

	  double xMin = m_dPhiZoomLow;
	  double xMax = hCFactors->GetXaxis()->GetXmax();
	  
	  TLine line( xMin, 1, xMax, 1 );
	  line.SetLineWidth( 2 );
	  line.Draw();

	  TLine lineP25( xMin, 1.25, xMax, 1.25 );
	  lineP25.SetLineStyle( 3  );
	  lineP25.SetLineColor( 12 );
	  lineP25.SetLineWidth( 2  );
	  lineP25.Draw();
	  
	  TLine lineN25( xMin, 0.75, xMax, 0.75 );
	  lineN25.SetLineStyle( 3 );
	  lineN25.SetLineColor( 12 );
	  lineN25.SetLineWidth( 2  );
	  lineN25.Draw();

	  TLine lineP50( xMin, 1.50, xMax, 1.50 );
	  lineP50.SetLineStyle( 3  );
	  lineP50.SetLineColor( 12 );
	  lineP50.SetLineWidth( 1  );
	  lineP50.Draw();
	  
	  TLine lineN50( xMin, 0.50, xMax, 0.50 );
	  lineN50.SetLineStyle( 3 );
	  lineN50.SetLineColor( 12 );
	  lineN50.SetLineWidth( 1  );
	  lineN50.Draw();

	  pad2.Draw();
	  c.cd();
	  
	  SaveAsAll( c, Form( "h_%s_%s_%s_%s",
			      m_dPhiUnfoldedName.c_str(),
			      m_allName.c_str(), m_sMUT.c_str(),
			      hTag.c_str()));

	} // end loop over axis3
      } // end loop over axis2
    } // end loop over axis1     
  } // end loop over axis0
  for( auto f  : vDphiFits    ){ delete f;  }
  for( auto h  : vDphi        ){ delete h ; }

  for( auto rm : vDphiRespMat ){ delete rm; }
  for( auto r  : vCFactors    ){ delete r;  }

  hnUnfolded->Write();
  
  return hnUnfolded;
}


void DiJetAnalysis::MakeDphiTogether(){

  std::vector< TF1* > vF;
  std::vector< TH1* > vH;
  std::vector< TH1* > vHw;
  std::vector< TH1* > vR;
  
  std::string outSuffix;
  std::string name_a  , name_b  ;
  std::string label_a , label_b ;
  std::string suffix_a, suffix_b;

  GetInfoBoth( outSuffix, name_a, name_b, label_a, label_b, suffix_a, suffix_b );

  std::string ratio = Form("%s/%s", label_a.c_str(), label_b.c_str() );
  
  // Check if the directories exist.
  // If they don't, create them
  std::string outDir = m_sOutput;
  anaTool->CheckWriteDir( outDir.c_str() );
  outDir += "/" + m_allName;
  anaTool->CheckWriteDir( outDir.c_str() );
  outDir += "/" + outSuffix;
  anaTool->CheckWriteDir( outDir.c_str() );
  
  TAxis* axis0 = m_dPP->GetTAxis(0); int nAxis0Bins = axis0->GetNbins();
  TAxis* axis1 = m_dPP->GetTAxis(1); int nAxis1Bins = axis1->GetNbins();
  TAxis* axis2 = m_dPP->GetTAxis(2); int nAxis2Bins = axis2->GetNbins();
  TAxis* axis3 = m_dPP->GetTAxis(3); int nAxis3Bins = axis3->GetNbins();

  // lines to be drawn along axis3. this is
  // x-axis that widths are plotted as function of 
  double xMin = axis3->GetXmin(); double xMax = axis3->GetXmax();
  
  TLine line( xMin, 1, xMax, 1 );
  line.SetLineWidth( 2 );

  TLine lineP25( xMin, 1.25, xMax, 1.25 );
  lineP25.SetLineStyle( 3  );
  lineP25.SetLineColor( 12 );
  lineP25.SetLineWidth( 2  );
	  
  TLine lineN25( xMin, 0.75, xMax, 0.75 );
  lineN25.SetLineStyle( 3 );
  lineN25.SetLineColor( 12 );
  lineN25.SetLineWidth( 2  );

  TLine lineP50( xMin, 1.50, xMax, 1.50 );
  lineP50.SetLineStyle( 3  );
  lineP50.SetLineColor( 12 );
  lineP50.SetLineWidth( 1  );
	  
  TLine lineN50( xMin, 0.50, xMax, 0.50 );
  lineN50.SetLineStyle( 3 );
  lineN50.SetLineColor( 12 );
  lineN50.SetLineWidth( 1  );
  
  TFile* fIn_a = TFile::Open
    ( Form("%s/%s_%s/c_%s_%s_UF.root",
	   m_sOutput.c_str(), m_sOutput.c_str(),suffix_a.c_str(),
	   m_myOutName.c_str(),suffix_a.c_str() ) );
  TFile* fIn_b = TFile::Open
    ( Form("%s/%s_%s/c_%s_%s_UF.root",
	   m_sOutput.c_str(), m_sOutput.c_str(),suffix_b.c_str(),
	   m_myOutName.c_str(),suffix_b.c_str() ) );
  TFile* fOut  = new TFile
    ( Form("%s/%s/%s/c_%s_%s.root",
	   m_sOutput.c_str(), m_allName.c_str(), outSuffix.c_str(),
	   m_myOutName.c_str(), outSuffix.c_str() ) ,"recreate");

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

      // Now Draw everything.
      TCanvas cW( "cW", "cW", 800, 700 );
      TPad pad1W("pad1W", "", 0.0, 0.25, 1.0, 1.0 );
      pad1W.SetBottomMargin(0);
      pad1W.Draw();
      TPad pad2W("pad2W", "", 0.0, 0.0, 1.0, 0.25 );
      pad2W.SetTopMargin(0);
      pad2W.SetBottomMargin(0.2);
      pad2W.Draw();

      // tag for widths canvas
      std::string hTagCW =
	Form ("%s_%s",
	      anaTool->GetName( axis0Low, axis0Up, m_dPP->GetAxisName(0) ).c_str(),
	      anaTool->GetName( axis1Low, axis1Up, m_dPP->GetAxisName(1) ).c_str() );
      
      TLegend legW( 0.33, 0.13, 0.87, 0.26 );
      styleTool->SetLegendStyle( &legW );
      legW.SetNColumns(2);
      
      int style = 0;

      std::vector< TH1* > vHw;
      
      for( int axis2Bin = 1; axis2Bin <= nAxis2Bins; axis2Bin++ ){
	double axis2Low , axis2Up;
	anaTool->GetBinRange
	  ( axis2, axis2Bin, axis2Bin, axis2Low, axis2Up );
	
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

	// switch back to cW canvas
	// because a new one was creted in axis3 loop
	cW.cd();
	pad1W.cd();

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

	pad2W.cd();
	// Make the ratio histogram.
	if( hW_a->GetMaximum() && hW_b->GetMean() ){
	  TH1* hW_R = static_cast< TH1D* >
	    ( hW_a->Clone( Form( "h_%s_%s_%s_%s", m_dPhiName.c_str(), m_sRatio.c_str(),
		            m_allName.c_str(), hTagW.c_str())));
	  hW_R->SetMinimum( 0.25 );
	  hW_R->SetMaximum( 1.75 );
	  hW_R->Divide( hW_b );
	  hW_R->SetYTitle( ratio.c_str() );
	  vR.push_back( hW_R );
	  hW_R->Draw("ep same");
	}
	
	for( int axis3Bin = 1; axis3Bin <= nAxis3Bins; axis3Bin++ ){
	  // check we are in correct ystar and pt bins
	  if( !m_dPP->CorrectPhaseSpace
	      ( std::vector<int>{ axis0Bin, axis1Bin, axis2Bin, axis3Bin } ) )
	    { continue; }

	  double axis3Low , axis3Up;
	  anaTool->GetBinRange
	    ( axis3, axis3Bin, axis3Bin, axis3Low, axis3Up );

	  TCanvas c( "c", "c", 800, 700 );
	  TPad pad1("pad1", "", 0.0, 0.25, 1.0, 1.0 );
	  pad1.SetBottomMargin(0);
	  pad1.Draw();
	  TPad pad2("pad2", "", 0.0, 0.0, 1.0, 0.25 );
	  pad2.SetTopMargin(0);
	  pad2.SetBottomMargin(0.2);
	  pad2.Draw();
	  
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
	  vH.push_back( h_a ); vH.push_back( h_b );

	  TH1* h_R = static_cast< TH1D* >
	    ( h_a->Clone( Form( "h_%s_%s_%s_%s", m_dPhiName.c_str(), m_sRatio.c_str(),
				 m_allName.c_str(), hTagW.c_str())));
	  styleTool->SetHStyleRatio( h_R );
	  h_R->Divide( h_b );
	  h_R->SetYTitle( ratio.c_str() );
	  vR.push_back( h_R );
	  h_R->Draw("ep same");
	  
	  TF1* f_a = static_cast<TF1*>( fIn_a->Get( Form("f_%s", hName_a.c_str())));
	  TF1* f_b = static_cast<TF1*>( fIn_b->Get( Form("f_%s", hName_b.c_str())));
	  styleTool->SetHStyle( f_a, 0 );
	  styleTool->SetHStyle( f_b, 1 );
	  f_a->SetLineColor( h_a->GetLineColor() );
	  f_b->SetLineColor( h_b->GetLineColor() );
	  vF.push_back( f_a ); vF.push_back( f_b );
	  
	  double chi2NDF_a = f_a->GetChisquare()/f_a->GetNDF();
	  double chi2NDF_b = f_b->GetChisquare()/f_b->GetNDF();
	    
	  TLegend leg( 0.27, 0.41, 0.38, 0.52 );
	  styleTool->SetLegendStyle( &leg , 0.85 );

	  bool save = false;

	  pad1.cd();
	  if( h_a->GetEntries() ){
	    h_a->SetMinimum(0);
	    h_a->SetNdivisions( 504, "Y" );
	    leg.AddEntry( h_a, Form("%s #Chi^{2}/NDF=%4.2f", label_a.c_str(), chi2NDF_a));
	    h_a->Draw("epsame");
	    f_a->Draw("same");
	    save = true;
	  }

	  if( h_b->GetEntries() ){
	    h_b->SetMinimum(0);
	    h_b->SetNdivisions( 504, "Y" );
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

	  pad2.cd(); 
	  
	  if( save ) c.SaveAs( Form("%s/%s/%s/h_%s_%s_%s.pdf",
				    m_sOutput.c_str(), m_allName.c_str(),
				    outSuffix.c_str(), m_dPhiName.c_str(),
				    hTag.c_str(),outSuffix.c_str() ) );
	  SaveAsROOT( c, Form("h_%s_%s", m_dPhiName.c_str(), hTag.c_str() ) );
	  
	} // end loop over axis3
      } // end loop over axis2

      // back to cW canvas
      cW.cd();
      pad1W.cd();
      
      legW.Draw("same");

      DrawTopLeftLabels
	( m_dPP, axis0Low, axis0Up, axis1Low, axis1Up,
	  0, 0, 0, 0, 0.8 );

      DrawAtlasRightBoth();

      pad2W.cd();
      line.Draw();
      lineP25.Draw();
      lineN25.Draw();
      lineP50.Draw();
      lineN50.Draw();
      
      cW.SaveAs( Form("%s/%s/%s/h_%s_%s_%s.pdf",
		      m_sOutput.c_str(), m_allName.c_str(),
		      outSuffix.c_str(), m_dPhiName.c_str(),
		      hTagCW.c_str(),outSuffix.c_str() ) );
      SaveAsROOT( cW, Form("h_%s_%s", m_dPhiName.c_str(), hTagCW.c_str() ) );

    } // end loop over ystar2
  } // end loop over ystar2

  for( auto & f  : vF  ){ delete f ; }
  for( auto & r  : vR  ){ delete r ; }
  for( auto & h  : vH  ){ delete h ; }
  for( auto & hW : vHw ){ delete hW; }
  
  std::cout << "DONE! Closing " << fOut->GetName() << std::endl;
  fOut->Close();
  std::cout << "......Closed  " << fOut->GetName() << std::endl;
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
  // only save png pdf if the doing nominal sample.
  // otherwise it is pretty useles.
  if( !m_uncertComp )
    { SaveAsPdfPng( c, label1, label2, axis1, min1, max1, axis2, min2, max2 ); }
  SaveAsROOT  ( c, label1, label2, axis1, min1, max1, axis2, min2, max2 );
}
