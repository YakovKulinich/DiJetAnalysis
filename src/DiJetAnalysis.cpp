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
  m_sFinal    = "final";
  
  m_sMUT      = "MUT";
  m_sRatio    = "ratio";
  m_sRebin    = "rebin";  
  
  // name for "All" histos
  // this is either merged JZN or Data from Triggers
  m_allName = "All";

  m_unfoldingFileSuffix = "UF";
  m_systematicsFileSuffix = "SYS";
  
  boost::assign::push_back( m_vMCtypeNames )
    ( "pythia8" )( "herwig" )( "pythia8powheg" );

  boost::assign::push_back( m_vMCtypeLabels  )
    ( "Pythia8" )( "Herwig" )( "Pythia8+Powheg" );

  //==================== Cuts ====================
  m_nMinEntriesFit = 20;
    
  m_dPhiThirdJetFraction = 0.4;

  // where plots are drawn from on dphi axis
  m_dPhiZoomLow      = constants::PI / 2;
  m_dPhiZoomHigh     = constants::PI;

  m_dPhiLogMin       = 1E-3;
  
  // where to fit from 
  m_dPhiFittingMin   = 2 * constants::PI / 3;
  m_dPhiFittingMax   = constants::PI;
  
  // sets range of unfolding and where
  // correction factors are derived
  m_dPhiUnfoldingMin = 0.0;
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
  m_ptSpectMin   = 20;
  m_ptSpectMax   = 100;
  m_nPtSpectBins = ( m_ptSpectMax - m_ptSpectMin ) / 2 ; 

  // -------- eff ---------
  m_effMin = 0.;
  m_effMax = 1.4;
  
  // ---- forward eta binning ---
  boost::assign::push_back( m_varFwdEtaBinning )
    ( -4.0 )( -3.3)( -3.1 )( -2.7 );
  m_nVarFwdEtaBins = m_varFwdEtaBinning.size() - 1; 

  // --- variable eta/ystar binning ---
  // ystar
  boost::assign::push_back( m_varYstarBinning )
    ( -4.0 )( -2.7 )( -1.8 )( 0 )( 1.8 )( 4.0 );
  m_nVarYstarBins = m_varYstarBinning.size() - 1;
  
  // --- variable pt binning ---
  boost::assign::push_back( m_varPtBinning )
    ( 28 )( 35 )( 45 )( 90 );
  m_nVarPtBins = m_varPtBinning.size() - 1;


  // --- variable dPhi binning ----
  int nDphiBinsLarge = 4;

  // MakeDefaultBinning( m_varDphiBinning, m_varDphiRebinnedBinning, nDphiBinsLarge );
  MakeLinearBinning( m_varDphiBinning, m_varDphiRebinnedBinning, nDphiBinsLarge );
  // m_varDphiRebinnedBinning = m_varDphiBinning;
  
  m_nVarDphiBins         = m_varDphiBinning.size()         - 1;  
  m_nVarDphiRebinnedBins = m_varDphiRebinnedBinning.size() - 1;

  m_dPhiZoomLowBin  = nDphiBinsLarge;
  m_dPhiZoomHighBin = m_nVarDphiBins;
  
  m_dPhiRebinnedZoomLowBin  = nDphiBinsLarge;
  m_dPhiRebinnedZoomHighBin = m_nVarDphiRebinnedBins;
    
  int count = 0;
  for( auto & b : m_varDphiBinning ){ std::cout << count++ << "," << b << " -> "; }
  std::cout << " --- " << std::endl;
  count = 0;
  for( auto & b : m_varDphiRebinnedBinning ){ std::cout << count++ << "," << b << " -> "; }
  std::cout << " --- " << std::endl;

  std::cout << "new amount of bins " << m_nVarDphiRebinnedBins << std::endl;

  // --- dPhiBins ---  
  boost::assign::push_back( m_vNdPhiBins )
    ( m_nVarYstarBins )( m_nVarYstarBins )
    ( m_nVarPtBins    )( m_nVarPtBins    )
    ( m_nVarDphiBins  );
    
  boost::assign::push_back( m_vDphiMin )
    ( 0 )( 0 )( 0 )( 0 )( 0 );

  boost::assign::push_back( m_vDphiMax )
    ( 1 )( 1 )( 1 )( 1 )( 1 );

  m_nDphiDim     = m_vNdPhiBins.size();
  
  m_dPhiDphiMin  = 0;
  m_dPhiDphiMax  = constants::PI;
  
  m_dPhiWidthMin = 0.00;
  m_dPhiWidthMax = 0.5;

  // default values 0->2 are in
  // SetHStyleRatio.
  // This is for "nonstandard"
  m_ratioMin = 0.0;
  m_ratioMax = 2.0;
    
  //========= Set DeltaPhi Axes Order ============
  // The DeltaPhiProj object will have the order
  // onto which to take projections. It also knows
  // what axis is what. The object knows which axis
  // has what name and range. I.e. Pt1 is actually
  // The third axis in the THnSparse. 
  m_dPP = new DeltaPhiProj( YS1, PT1, PT2, YS2 );

  m_dPP->AddTAxis( new TAxis( m_nVarYstarBins, &m_varYstarBinning[0] ) );
  m_dPP->AddTAxis( new TAxis( m_nVarYstarBins, &m_varYstarBinning[0] ) );
  m_dPP->AddTAxis( new TAxis( m_nVarPtBins   , &m_varPtBinning[0]    ) );
  m_dPP->AddTAxis( new TAxis( m_nVarPtBins   , &m_varPtBinning[0]    ) );
  
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
  m_etaSpectName       = "etaSpect";
  m_dPhiName           = "dPhi";
  m_effName            = "eff";
  m_purityName         = "purity";

  m_recoName           = "reco";
  m_truthName          = "truth";
  m_respMatName        = "respMat";
  m_unfoldedName       = "unfolded";
  m_systematicsName    = "systematics";
  
  m_dPhiRecoName       = m_dPhiName + "_" + m_recoName;
  m_dPhiTruthName      = m_dPhiName + "_" + m_truthName;
  
  m_dPhiCFactorsName   =
    m_dPhiName + "_" + m_recoName + "_" + m_truthName + "_" + m_sRatio;

  m_dPhiUnfoldedName     = m_dPhiName  + "_" + m_unfoldedName;
  m_dPhiSystematicsName  = m_dPhiName  + "_" + m_systematicsName;
  
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
  m_uncertSuffix =
    m_uncertComp ? ( m_uncertComp > 0 ?
		     Form("P%d", m_uncertComp) :
		     Form("N%d", -1 * m_uncertComp) ) : "0" ;
  
  // this is the fname of the raw histograms. I.e.
  // what is produced from reading the tree (THnMulfs, TH3, TH2 )
  // These get used to take projections, etc, later.
  m_rawHistosFname =
    m_dirOut + "/" + m_myOutName +"_" + m_labelOut + "_" + m_uncertSuffix + ".root";
  
  std::cout << "fNameIn/Out: " << m_rawHistosFname << std::endl;
  
  std::string unfoldingMC = m_mcTypeName;

  // set names for various output files.
  // for data only, we just have one.
  // no need to add unfolding suffix.
  m_fNameOutDefault =
    m_dirOut + "/c_" + m_myOutName + "_" + m_labelOut + "_0";

  m_fNameOutDefaultUF =
    m_fNameOutDefault + "_" + m_unfoldingFileSuffix;
  
  m_fNameOut =
    m_dirOut + "/c_" + m_myOutName + "_" + m_labelOut + "_" + m_uncertSuffix;

  m_fNameOutUF =
    m_fNameOut + "_" + m_unfoldingFileSuffix;

  m_fNameOutSYS =
    m_fNameOutDefault + "_" + m_systematicsFileSuffix;

  m_fNameOutDefault   += ".root";
  m_fNameOut          += ".root";
  m_fNameOutDefaultUF += ".root";
  m_fNameOutUF        += ".root";
  m_fNameOutSYS       += ".root";

  m_fNameRivetMC = "data/rivet_" + system + ".root";
  
  // name of file that is used as input for unfolding
  m_fNameUnfoldingMC = Form( "%s/%s_%s_%s_%s/c_%s_%s_%s_%s_%s.root",
			     m_sOutput.c_str(), m_sOutput.c_str(),
			     system.c_str(), m_sMC.c_str(),
			     unfoldingMC.c_str(), m_myOutName.c_str(),
			     system.c_str(), m_sMC.c_str(),
			     unfoldingMC.c_str(), m_uncertSuffix.c_str() );
}

//---------------------------
//       Fill Tree
//---------------------------
void DiJetAnalysis::AddHistogram( TH1* h ){
  
  v_hists.push_back( h );
  h->Sumw2();
  styleTool->SetHStyle( h, 0 );

  TH1* h1 = dynamic_cast< TH1* >(h);
  if( h1 ){ h->SetNdivisions( 505, "X" ); }

  TH2* h2 = dynamic_cast< TH2* >(h);
  if( h2 ){ h->SetNdivisions( 505, "XY" ); }

  TH3* h3 = dynamic_cast< TH3* >(h);
  if( h3 ){ h->SetNdivisions( 505, "XYZ" ); }
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
				       const std::vector< TLorentzVector >& v_jets ){
  
  const TLorentzVector* jet1 = NULL; const TLorentzVector* jet2 = NULL;

  if( !GetDiJets( v_jets, jet1, jet2 ) )
    { return -1; }

  double jet1_pt    = jet1->Pt()/1000.;
  double jet1_ystar = GetYstar( *jet1 );

  double jet2_pt    = jet2->Pt()/1000.;
  double jet2_ystar = GetYstar( *jet2 );
  
  double deltaPhi = anaTool->DeltaPhi( *jet2, *jet1 );
  
  std::vector< double > x;
  x.resize( hn->GetNdimensions() );
    
  double jetWeight = GetJetWeight( *jet1 );

  x[0] = jet1_ystar;  
  x[1] = jet2_ystar;
  x[2] = jet1_pt ;
  x[3] = jet2_pt ;
  x[4] = deltaPhi;
  hn->Fill( &x[0], jetWeight );
  
  // for pp, fill twice. once for each side since
  // it is symmetric in pp. For pPb, continue
  if( m_is_pPb ){ return deltaPhi; }
  
  x[0] = -jet1_ystar;  
  x[1] = -jet2_ystar;
  hn->Fill( &x[0], jetWeight );
  
  return deltaPhi;
}

double DiJetAnalysis::AnalyzeDeltaPhi( THnSparse* hn,
				       THnSparse* hnMig,
				       const std::vector< TLorentzVector >& v_rJets,
				       const std::vector< TLorentzVector >& v_tJets ){
  
  const TLorentzVector* rJet1 = NULL; const TLorentzVector* rJet2 = NULL;
  const TLorentzVector* tJet1 = NULL; const TLorentzVector* tJet2 = NULL;

  if( !GetDiJets( v_rJets, rJet1, rJet2 ) )
    { return -1; }

  if( !GetDiJets( v_tJets, tJet1, tJet2 ) )
    { return -1; }  

  double tJet2_pt    = tJet2->Pt()/1000.;
  
  double rJet1_pt    = rJet1->Pt()/1000.;
  double rJet1_ystar = GetYstar( *rJet1 );

  double rJet2_pt    = rJet2->Pt()/1000.;
  double rJet2_ystar = GetYstar( *rJet2 );

  double deltaPhi = anaTool->DeltaPhi( *rJet2, *rJet1 );

  bool differentBin = false;
  
  for( uint bin = 0; bin < m_varPtBinning.size() - 1; bin++ ){
    if( ! ( tJet2_pt >= m_varPtBinning[ bin ] && tJet2_pt < m_varPtBinning[ bin + 1 ] ) &&
	( rJet2_pt >= m_varPtBinning[ bin ] && rJet2_pt < m_varPtBinning[ bin + 1 ] ) )
      { differentBin = true;
	// std::cout << rJet2_pt << " " << tJet2_pt << std::endl;
	break;
      }
  }
  
  std::vector< double > x;
  x.resize( hn->GetNdimensions() );
    
  double rJetWeight = GetJetWeight( *rJet1 );
  
  x[0] = rJet1_ystar;  
  x[1] = rJet2_ystar;
  x[2] = rJet1_pt ;
  x[3] = rJet2_pt ;
  x[4] = deltaPhi;
  hn->Fill( &x[0], rJetWeight );

  if( differentBin ){ hnMig->Fill( &x[0], rJetWeight ); }
  
  // for pp, fill twice. once for each side since
  // it is symmetric in pp. For pPb, continue
  if( m_is_pPb ){ return deltaPhi; }
  
  x[0] = -rJet1_ystar;  
  x[1] = -rJet2_ystar;
  hn->Fill( &x[0], rJetWeight );

  if( differentBin ){ hnMig->Fill( &x[0], rJetWeight ); }
  
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
  
  if( m_is_pPb ){ return jetEta; }
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

double DiJetAnalysis::GetJetWeight( const TLorentzVector& jet )
{ return 1; }

double DiJetAnalysis::GetUncertaintyWeight( const TLorentzVector& jet1,
					    const TLorentzVector& jet2 )
{ return 1; }

TH1* DiJetAnalysis::CombineSamples( std::vector< TH1* >& vSampleHin,
				    const std::string& name ){ return NULL; }  

TH2* DiJetAnalysis::CombineSamples( std::vector< TH2* >& vSampleHin,
				    const std::string& name ){ return NULL; }  

THnSparse* DiJetAnalysis::CombineSamples( std::vector< THnSparse* >& vSampleHin,
					  const std::string& name ){ return NULL; }  


void DiJetAnalysis::GetInfoBoth( std::string& name_a  , std::string& name_b  ,
				 std::string& label_a , std::string& label_b ,
				 std::string& suffix_a, std::string& suffix_b ){}

void DiJetAnalysis::GetInfoUnfolding( std::string& measuredName,
				      std::string& measuredLabel,
				      std::string& typeLabel ){}

// hM is measured, hC is correction factor
TH1* DiJetAnalysis::BinByBinUnfolding( TH1* hM, TH1* hC ){
  
  TH1* hUnf = static_cast<TH1D*>( hM->Clone("h_unfolded") );
  hUnf->Reset();
  
  for( int xBin = hM->FindBin( m_dPhiUnfoldingMin ); xBin <= hM->GetNbinsX(); xBin++ ){
    int cBin  = hC->FindBin( hM->GetBinCenter( xBin ) );
    
    double vM = hM->GetBinContent( xBin );
    double vC = hC->GetBinContent( cBin );

    if( vM == 0 || vC == 0 ){ continue; }
    
    double eM = hM->GetBinError( xBin );

    // correction factor;
    double newDphi = vM * vC;
    // error on correction factor    
    double newDphiError = eM * vC; 
    
    hUnf->SetBinContent( xBin, newDphi      );
    hUnf->SetBinError  ( xBin, newDphiError );
  }
  
  return hUnf;
}

TFile* DiJetAnalysis::GetListOfSystUncert( std::vector< int >& v_uc,
					   std::map< int, TFile* >& mFinUC ){

  // read text file generated by
  // script to know what factors used
  std::ifstream file( "uc.txt" );
  int number;

  TFile* fInDefault = NULL;
  
  while( file >> number )
    { v_uc.push_back( number ); }
  
  for( auto uc : v_uc ){
    m_uncertComp = uc;
    Initialize();

    if( !m_fNameOutUF.compare( m_fNameOutDefaultUF ) ){
      std::cout << "Found default sample: " << m_fNameOutUF << std::endl;
      fInDefault = TFile::Open( m_fNameOutUF.c_str() );
    } else {
      std::cout << "Adding  to systematics: " << m_fNameOutUF << std::endl;
      mFinUC[ uc ] =  TFile::Open( m_fNameOutUF.c_str() ) ;
    }
  }
  
  return fInDefault;
}

void DiJetAnalysis::MakeDefaultBinning( std::vector< double >& varBinning,
					std::vector< double >& varRebinnedBinning,
					int nLargeBins ){

  std::cout << "DEFAULT BINNING" << std::endl;

  int nDphiBinsLarge  = nLargeBins ; int dPhiBinsLargeFactor  = 10; 
  int nDphiBinsMedium = 8          ; int dPhiBinsMediumFactor = 2;
  int nDphiBinsSmall  = 4          ; int dPhiBinsSmallFactor  = 1;

  double smallBinWidth =
    constants::PI / ( nDphiBinsSmall  * dPhiBinsSmallFactor  +
		      nDphiBinsMedium * dPhiBinsMediumFactor +
		      nDphiBinsLarge  * dPhiBinsLargeFactor );
  for( int largeBin = 0; largeBin < nDphiBinsLarge; largeBin++ )
    { varBinning.push_back
	( largeBin * smallBinWidth * dPhiBinsLargeFactor ); }

  double mediumBinStart =
    nDphiBinsLarge * smallBinWidth * dPhiBinsLargeFactor;
  for( int mediumBin = 0; mediumBin < nDphiBinsMedium; mediumBin++ )
    { varBinning.push_back
	( mediumBin * smallBinWidth * dPhiBinsMediumFactor + mediumBinStart ); }

  double smallBinStart =
    nDphiBinsLarge  * smallBinWidth * dPhiBinsLargeFactor +
    nDphiBinsMedium * smallBinWidth * dPhiBinsMediumFactor;
  for( int smallBin = 0; smallBin <= nDphiBinsSmall; smallBin++ )
    { varBinning.push_back
	( smallBin * smallBinWidth * dPhiBinsSmallFactor + smallBinStart ); }
  
  varRebinnedBinning.push_back( varBinning[0]  );
  varRebinnedBinning.push_back( varBinning[3]  );
  varRebinnedBinning.push_back( varBinning[5]  );
  varRebinnedBinning.push_back( varBinning[7]  );
  varRebinnedBinning.push_back( varBinning[9]  );
  varRebinnedBinning.push_back( varBinning[11] );
  varRebinnedBinning.push_back( varBinning[12] );
  varRebinnedBinning.push_back( varBinning[13] );
  varRebinnedBinning.push_back( varBinning[14] );
  varRebinnedBinning.push_back( varBinning[15] );
  varRebinnedBinning.push_back( varBinning[16] );
}
  
void DiJetAnalysis::MakeLinearBinning( std::vector< double >& varBinning,
				       std::vector< double >& varRebinnedBinning,
				       int nLargeBins ){

  std::cout << "LINEAR BINNING" << std::endl;
  
  double varBinStart = 2 * constants::PI / 3;
  double dPhiBinsLargeWidth = ( varBinStart ) / nLargeBins;
  
  for( int i = 0; i <= nLargeBins; i++ ){
    varBinning.        push_back( i * dPhiBinsLargeWidth );
    varRebinnedBinning.push_back( i * dPhiBinsLargeWidth );
  }

  int    nVarBins    = 12;
  double finalWidth  = 0.05;

  double dWidthPerBin =
    ( constants::PI - varBinStart - nVarBins * finalWidth ) /
    ( nVarBins / 2 - ( nVarBins * nVarBins ) / 2 );

  std::cout << " !!!! dWidthPerBin : " << dWidthPerBin << std::endl;
  
  for( int i = 1; i <= nVarBins; i++ ){
    double width = finalWidth - ( nVarBins - i ) * dWidthPerBin;
    std:: cout << i << " : " << width << std::endl;
    varBinning.push_back( varBinning.back() + width );
  }

  // for now, do a rebin of 2 in middle
  // be careful total bins is a multiple of 2
  for( int i = nLargeBins + 2; i < nVarBins + nLargeBins - 4; i += 2 ){
    varRebinnedBinning.push_back( varBinning[ i ] );
  }

  // leave the last 4 as is.
  for( int i = nVarBins + nLargeBins - 4; i <= nVarBins + nLargeBins; i++ ){
    varRebinnedBinning.push_back( varBinning[ i ] );
  }
}
 
void DiJetAnalysis::MakeLogBinning( std::vector< double >& varBinning,
				    std::vector< double >& varRebinnedBinning,
				    int nLargeBins ){
  
  std::cout << "LOGARITHMIC BINNING" << std::endl;

  double varBinStart = 2 * constants::PI / 3;

  int nLogBins = 12;
   
  double logMin   = TMath::Log10( varBinStart   );
  double logMax   = TMath::Log10( constants::PI );
  double logWidth = ( logMax - logMin ) / nLogBins; 

  std::vector< double > varWidths;
  
  for( int i = 0; i < nLogBins; i++ ){
    varWidths.push_back( TMath::Power( 10, logMin + logWidth * ( i + 1 ) ) -
			 TMath::Power( 10, logMin + logWidth * i ) );
  }

  double dPhiBinsLargeWidth = ( varBinStart ) / nLargeBins;
  
  for( int i = 0; i <= nLargeBins; i++ ){
    varBinning.push_back( i * dPhiBinsLargeWidth );
  }

  for( int i = 1; i <= nLogBins; i++ ){
    varBinning.push_back( varBinning.back() + varWidths[ nLogBins - i ] );
  }
  
  varRebinnedBinning  = varBinning;
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

      h_etaSpect->Write();
      // get min max from the final histograms
      if( label.compare( m_allName ) ){ continue; }
      if( max < h_etaSpect->GetMaximum() )
	{ max = h_etaSpect->GetMaximum(); }
     } // end loop over xBin
  } // end loop over iG

  // set maxima globally for all spectra hists.
  // easier to compare. Set on log scale.
  max = anaTool->GetLogMaximum( max );
  
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

    SaveAsAll( c, Form("%s_%s", name.c_str(), cName.c_str() ) );
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

    SaveAsAll( c, Form("%s_%s", name.c_str(), cName.c_str() ) );
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

  // These are "quality" check histograms
  h_mult     = new TH1D( Form("h_mult_%s", name.c_str() ),
			 ";Count;Multiplicity", 100,0,5000);
  h_tauSigma = new TH2D( Form("h_tauSigma_%s", name.c_str() ),
			 ";#tau;#sigma"   , 60, 0., 0.6, 60, 0., 0.3 );
  h_tauAmp   = new TH2D( Form("h_tauAmp_%s"  , name.c_str() ),
			 ";#tau;Amplitude", 60, 0., 0.6, 60, 0., 0.3 );
  h_tauConst = new TH2D( Form("h_tauConst_%s", name.c_str() ),
			 ";#tau;Constant" , 60, 0., 0.6, 60, 0., 0.002 );
  
  std::vector< TH1* > vDphi;
  std::vector< TF1* > vFits;
  std::vector< THnSparse* > vHnDphi;
  
  // ---- loop over group  ----
  // ---- (jzn or trigger) ----
  for( uint iG = 0; iG < vhn.size(); iG++ ){      
    THnSparse* hn = vhn[iG];
  
    std::string label = vLabel[iG];

    if( label.compare( m_allName ) ){ continue; }

    THnSparse* hnDphi = static_cast< THnSparse* >
      ( hn->Clone( Form( "h_%s_%s",  name.c_str(), label.c_str() ) ) );
    hnDphi->Reset();
    vHnDphi.push_back( hnDphi );
    
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
	  double axis2Low , axis2Up;
	  anaTool->GetBinRange
	    ( axis2, axis2Bin, axis2Bin, axis2Low, axis2Up );

	  std::string hTagW =
	  Form( "%s_%s_%s",
	     	anaTool->GetName( axis0Low, axis0Up, m_dPP->GetAxisName(0) ).c_str(),
		anaTool->GetName( axis1Low, axis1Up, m_dPP->GetAxisName(1) ).c_str(),
		anaTool->GetName( axis2Low, axis2Up, m_dPP->GetAxisName(2) ).c_str() ); 

	  // do this becaue it is set to last bin after one iteration
	  axis3->SetRange( 1, -1 );
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

	    h_mult    ->Fill( hDphi->GetEntries() );
	    
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
	    
	    // fill output THnSparse result
	    std::vector< int > x  = m_dPP->GetMappedBins
	      ( std::vector<int> { axis0Bin, axis1Bin, axis2Bin, axis3Bin } );
	    x.push_back(0); // dPhiBin;
	    for( int dPhiBin = 1; dPhiBin <= hDphi->GetNbinsX(); dPhiBin++ ){
	      x[4] = dPhiBin;
	      hnDphi->SetBinContent( &x[0], hDphi->GetBinContent( dPhiBin ) );
	      hnDphi->SetBinError  ( &x[0], hDphi->GetBinError  ( dPhiBin ) );
	    }

	    // now fit
	    TF1* fit = anaTool->FitDphi( hDphi, m_dPhiFittingMin, m_dPhiFittingMax );
	    styleTool->SetHStyle( fit, 0 );
	    vFits.push_back( fit );

	    h_tauAmp  ->Fill( fit->GetParameter(1), fit->GetParameter(0) );
	    h_tauSigma->Fill( fit->GetParameter(1), fit->GetParameter(2) );
	    h_tauConst->Fill( fit->GetParameter(1), fit->GetParameter(3) );

	    if( fit->GetParError(1) > 0.5 ) { std::cout << " !!!!!!!!!!!! " << fit->GetName()
							<< " " << fit->GetParameter(1)  
							<< " " << fit->GetParError(1)  
							<< " " << fit->GetParameter(2) 
							<< " " << fit->GetParError(2) << std::endl; }
	    
	    fit->Draw("same");

	    double chi2NDF = fit->GetChisquare()/fit->GetNDF();

	    hDphi->GetXaxis()->SetRange( m_dPhiZoomLowBin, m_dPhiZoomHighBin );

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
    // save THnSparse dPhi to file
    hnDphi->Write();
  } // end loop over iG

  
  h_mult    ->Write();
  h_tauSigma->Write();
  h_tauAmp  ->Write();
  h_tauConst->Write();
  delete h_mult;
  delete h_tauSigma;
  delete h_tauAmp;
  delete h_tauConst;
  
  for( auto& f  : vFits   ){ delete f;  }
  for( auto& h  : vDphi   ){ delete h;  }

  // Save the THnSparse to separate file.
  // This is used for unfolding later.
  // Gets put into data/ dirrectory.
  std::string fOutName =
    m_dirOut + "/" + name + "_" + m_labelOut + ".root";
  
  TFile* fOut = new TFile( fOutName.c_str(), "RECREATE" );

  for( auto& hn : vHnDphi ){ hn->Write(); delete hn; }

  fOut->Close();
}

THnSparse* DiJetAnalysis::UnfoldDeltaPhi( TFile* fInData, TFile* fInMC,
					  const std::string& hnUnfoldedName ){
  std::vector< TH1* > vDphi;
  std::vector< TF1* > vDphiFits;
  
  std::vector< TH1* > vCFactors;
  
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
		    m_nDphiDim, &m_vNdPhiBins[0],
		    &m_vDphiMin[0], &m_vDphiMax[0] );

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
  hnUnfolded->GetAxis(4)->Set( m_nVarDphiBins, &( m_varDphiBinning[0]  ) );

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
	  	  
	  TLegend leg( 0.54, 0.16, 0.66, 0.33 );
	  styleTool->SetLegendStyle( &leg , 0.85 );
	  
	  std::string hTag =
	    Form( "%s_%s_%s_%s",
		  anaTool->GetName( axis0Low, axis0Up, m_dPP->GetAxisName(0) ).c_str(),
		  anaTool->GetName( axis1Low, axis1Up, m_dPP->GetAxisName(1) ).c_str(),
		  anaTool->GetName( axis2Low, axis2Up, m_dPP->GetAxisName(2) ).c_str(),
		  anaTool->GetName( axis3Low, axis3Up, m_dPP->GetAxisName(3) ).c_str() ); 

	  // Get measured and truth distributions
	  TH1* hDphiMeasured = static_cast<TH1D*>
	    ( fInData->Get( Form( "h_%s_%s_%s", measuredName.c_str(), m_allName.c_str(), hTag.c_str())));
	  TH1* hDphiTruth    = static_cast<TH1D*>
	    (fInMC->Get( Form( "h_%s_%s_%s", m_dPhiTruthName.c_str(), m_allName.c_str(), hTag.c_str())));

	  // Get the correction factors and response matrix
	  TH1* hCFactors     = static_cast<TH1D*>
	    (fInMC->Get( Form( "h_%s_%s_%s", m_dPhiCFactorsName.c_str(), m_allName.c_str(), hTag.c_str())));

	  vDphi.push_back( hDphiMeasured );
	  vDphi.push_back( hDphiTruth    );

	  vCFactors.push_back( hCFactors     );
	  
       	  TF1* fitMeasured =  static_cast<TF1*>
	    ( fInData->Get( Form( "f_h_%s_%s_%s", measuredName.c_str(), m_allName.c_str(), hTag.c_str())));
	  TF1* fitTruth    =  static_cast<TF1*>
	    ( fInMC->Get( Form( "f_h_%s_%s_%s", m_dPhiTruthName.c_str(), m_allName.c_str(), hTag.c_str())));
	  
	  // ----------- Unfold -----------
	  // Unfold using bin-by-bin and the resposne factors.
	  TH1* hDphiUnfolded = BinByBinUnfolding
	    ( hDphiMeasured,  hCFactors );
	  hDphiUnfolded->SetName
	    ( Form( "h_%s_%s_%s", m_dPhiUnfoldedName.c_str(), m_allName.c_str(), hTag.c_str()));
	  
	  // fit with no combinatoric subtraction (already done);
	  TF1* fitUnfolded = anaTool->FitDphi( hDphiUnfolded, m_dPhiFittingMin, m_dPhiFittingMax );

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
	  
	  hDphiMeasured->GetXaxis()->SetRange( m_dPhiZoomLowBin, m_dPhiZoomHighBin );
	  hDphiUnfolded->GetXaxis()->SetRange( m_dPhiZoomLowBin, m_dPhiZoomHighBin );
	  hDphiTruth   ->GetXaxis()->SetRange( m_dPhiZoomLowBin, m_dPhiZoomHighBin );
	  
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

	  maximum = anaTool->GetLogMaximum( maximum );
	  
	  hDphiUnfolded->SetMaximum( maximum );
	  hDphiUnfolded->SetMaximum( maximum );
	  hDphiTruth   ->SetMaximum( maximum );

	  hDphiUnfolded->SetMinimum( m_dPhiLogMin );
	  hDphiUnfolded->SetMinimum( m_dPhiLogMin );
	  hDphiTruth   ->SetMinimum( m_dPhiLogMin );


	  pad1.SetLogy();
	  
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

	  drawTool->DrawLeftLatex( 0.57, 0.13, unfoldingMCLabel.c_str() );
	  
	  DrawTopLeftLabels
	    ( m_dPP, axis0Low, axis0Up, axis1Low, axis1Up,
	      axis2Low, axis2Up, axis3Low, axis3Up, 0.8 );
	    
	  DrawAtlasRight( -0.17, -0.1, 0.825 );

	  c.cd();
	  TPad pad2("pad2", "", 0.0, 0.0, 1.0, 0.25 );
	  pad2.SetTopMargin(0);
	  pad2.SetBottomMargin(0.2);
	  pad2.Draw();
	  pad2.cd();

	  hCFactors->SetYTitle( "Truth/Reco" );
	  hCFactors->GetXaxis()->SetRange
	    ( m_dPhiRebinnedZoomLowBin, m_dPhiRebinnedZoomHighBin );
	  hCFactors->GetXaxis()->SetNdivisions( 504 );	  
	  hCFactors->Draw("ep");

	  double xMin = m_dPhiZoomLow;
	  double xMax = m_dPhiZoomHigh; 
	  
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

  for( auto r  : vCFactors    ){ delete r;  }

  hnUnfolded->Write();
  
  return hnUnfolded;
}

void DiJetAnalysis::MakeDphiTogether(){

  std::vector< TF1* > vF;
  std::vector< TH1* > vH;
  std::vector< TH1* > vHw;
  std::vector< TH1* > vR;
  
  std::string name_a , name_b  ;
  std::string label_a, label_b ;
  std::string fName_a, fName_b;

  GetInfoBoth( name_a, name_b, label_a, label_b, fName_a, fName_b );

  std::string ratio = Form("%s/%s", label_a.c_str(), label_b.c_str() );
   
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
  
  TFile* fIn_a = TFile::Open( fName_a.c_str() );
  TFile* fIn_b = TFile::Open( fName_b.c_str() );
  TFile* fOut  = new TFile( m_fNameOutTogether.c_str() ,"recreate");
  
  //  bool isRivet = false;
  //  if( fName_b.find( "rivet" ) != std::string::npos ){ isRivet = true; }
  
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
      
      TLegend legW( 0.10, 0.08, 0.9, 0.22 );
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
	
	std::string hNameW_a = "h_" + name_a + "_" + hTagW;
	std::string hNameW_b = "h_" + name_b + "_" + hTagW;

	TH1* hW_a = static_cast<TH1D*>( fIn_a->Get( hNameW_a.c_str() ) );
	TH1* hW_b = static_cast<TH1D*>( fIn_b->Get( hNameW_b.c_str() ) );

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

	  std::string hName_a = "h_" + name_a + "_" + hTag;
	  std::string hName_b = "h_" + name_b + "_" + hTag;

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
	  h_R->GetXaxis()->SetNdivisions( 504 );
	  h_R->SetYTitle( ratio.c_str() );
	  vR.push_back( h_R );
	  
	  TF1* f_a = static_cast<TF1*>( fIn_a->Get( Form("f_%s", hName_a.c_str())));
	  TF1* f_b = static_cast<TF1*>( fIn_b->Get( Form("f_%s", hName_b.c_str())));
	  styleTool->SetHStyle( f_a, 0 );
	  styleTool->SetHStyle( f_b, 1 );
	  f_a->SetLineColor( h_a->GetLineColor() );
	  f_b->SetLineColor( h_b->GetLineColor() );
	  vF.push_back( f_a ); vF.push_back( f_b );
	  
	  double chi2NDF_a = f_a->GetChisquare()/f_a->GetNDF();
	  double chi2NDF_b = f_b->GetChisquare()/f_b->GetNDF();

	  TLegend leg( 0.45, 0.10, 0.66, 0.28 );
	  styleTool->SetLegendStyle( &leg , 0.85 );

	  bool save = false ;
	  
	  pad1.cd();
	  pad1.SetLogy();
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

	  double max =  h_a->GetMaximum() > h_b->GetMaximum() ?
	    h_a->GetMaximum() : h_b->GetMaximum();

	  max = anaTool->GetLogMaximum( max );
	  
	  h_a->SetMaximum( max );
	  h_a->SetMinimum( m_dPhiLogMin );

	  h_b->SetMaximum( max );
	  h_b->SetMinimum( m_dPhiLogMin );
	
	  DrawTopLeftLabels
	    ( m_dPP, axis0Low, axis0Up, axis1Low, axis1Up,
	      axis2Low, axis2Up, axis3Low, axis3Up, 0.8 );

	  DrawAtlasRightBoth();

	  pad2.cd();
	  h_R->Draw("ep same");

	  // lines to be drawn along axis3. this is
	  // x-axis that widths are plotted as function of 
	  double xMin = m_dPhiZoomLow;
	  double xMax = m_dPhiZoomHigh;
  
	  TLine dPhiLine( xMin, 1, xMax, 1 );
	  dPhiLine.SetLineWidth( 2 );
	  dPhiLine.Draw();
	  
	  TLine dPhiLineP25( xMin, 1.25, xMax, 1.25 );
	  dPhiLineP25.SetLineStyle( 3  );
	  dPhiLineP25.SetLineColor( 12 );
	  dPhiLineP25.SetLineWidth( 2  );
	  dPhiLineP25.Draw();
	  
	  TLine dPhiLineN25( xMin, 0.75, xMax, 0.75 );
	  dPhiLineN25.SetLineStyle( 3 );
	  dPhiLineN25.SetLineColor( 12 );
	  dPhiLineN25.SetLineWidth( 2  );
	  dPhiLineN25.Draw();
	  
	  TLine dPhiLineP50( xMin, 1.50, xMax, 1.50 );
	  dPhiLineP50.SetLineStyle( 3  );
	  dPhiLineP50.SetLineColor( 12 );
	  dPhiLineP50.SetLineWidth( 1  );
	  dPhiLineP50.Draw();
	  
	  TLine dPhiLineN50( xMin, 0.50, xMax, 0.50 );
	  dPhiLineN50.SetLineStyle( 3 );
	  dPhiLineN50.SetLineColor( 12 );
	  dPhiLineN50.SetLineWidth( 1  );
	  dPhiLineN50.Draw();
	  
	  if( save )
	    { SaveAsPdfPng( c, Form("h_%s_%s", m_dPhiName.c_str(), hTag.c_str() ), true ); }
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

      SaveAsPdfPng( cW, Form("h_%s_%s", m_dPhiName.c_str(), hTagCW.c_str() ), true );
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


void DiJetAnalysis::MakeFinalPlotsTogether(){
  
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

void DiJetAnalysis::SaveAsROOT( const TCanvas& c,
			        const std::string& name ){ 
  c.Write( Form( "c_%s" , name.c_str() ) );
}

void DiJetAnalysis::SaveAsPdfPng( const TCanvas& c,
				  const std::string& name,
				  bool together ){

  if( m_uncertComp && !together ){ return; }
  
  std::string outName = together ? m_dirOutTogether : m_dirOut;
  
  outName += "/" + name;
 
  std::string sPdf = outName + ".pdf";
  std::string sPng = outName + ".png";

  c.SaveAs( sPdf.c_str() );
  c.SaveAs( sPng.c_str() );
}

void DiJetAnalysis::SaveAsAll( const TCanvas& c,
			       const std::string& name,
			       bool together ){

  // only save png pdf if the doing nominal sample.
  // otherwise it is pretty useles.
  if( !m_uncertComp )
    { SaveAsPdfPng( c, name, together ); }
  SaveAsROOT  ( c, name );
}
