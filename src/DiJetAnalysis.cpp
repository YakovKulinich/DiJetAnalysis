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

#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldBinByBin.h"

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

  m_sEta    = "eta";
  m_sYstar  = "ystar";
  
  m_sOutput   = "output";  
  m_myOutName = "myOut";
  m_sMC       = "mc";
  m_sData     = "data";
  m_sFinal    = "final";

  m_sRaw      = "raw";
  m_sPerf     = "perf";
  m_sPhys     = "phys";

  m_sRuns     = "runs";
  m_sFine     = "fine";
  m_sWeights  = "weights";
  m_sCounts   = "counts";
  m_sReb      = "reb";
  m_sRatio    = "ratio";
  m_sSum      = "sum";

  m_unweightedFileSuffix  = "UW";
  m_unfoldingFileSuffix   = "UF";
  m_systematicsFileSuffix = "SYS";
  
  boost::assign::push_back( m_vMCtypeNames )
    ( "pythia8" )( "herwig" )( "pythia8powheg" );

  boost::assign::push_back( m_vMCtypeLabels  )
    ( "Pythia8" )( "Herwig" )( "Pythia8+Powheg" );

  //==================== Cuts ====================
  m_nMinEntriesFit = 20;
    
  m_dPhiThirdJetFraction = 0.4;

  // where plots are drawn from on dphi axis
  m_dPhiZoomLow      = 2 * constants::PI / 3;
  m_dPhiZoomHigh     = constants::PI;

  m_dPhiLogMin       = 5E-4;
  m_dPhiLogMax       = 1;

  // flag to fit with or without constant
  m_fitDphiWC = false;

  // where to fit from 
  m_dPhiFittingMin  = 2.48;
  m_dPhiFittingMax  = constants::PI;

  m_dPhiFittingMinB = 2 * constants::PI / 3;

  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // for this uncertainty, we change the fitting range.
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if( m_uncertComp == 22 ){
    m_dPhiFittingMin  = m_dPhiFittingMinB; 
  }
  
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
  m_ptSpectMin   = 0;
  m_ptSpectMax   = 100;
  m_nPtSpectBins = ( m_ptSpectMax - m_ptSpectMin ) / 2 ; 

  m_ptSpectYaxisMin = 1E1;
  
  // -------- eff ---------
  m_effMin = 0.;
  m_effMax = 1.5;
  
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

  // --- for resp mat ---
  // here we include underflow and overflow pt bins
  m_varPtBinningUfOf = m_varPtBinning;
  m_varPtBinningUfOf.insert( m_varPtBinningUfOf.begin(), 20  );
  m_varPtBinningUfOf.insert( m_varPtBinningUfOf.end()  , 120 );
  m_nVarPtBinsUfOf   = m_varPtBinningUfOf.size() - 1;
  
  // --- variable dPhi binning ----
  // int nDphiBinsSmall  = 5;
  // int nDphiBinsMedium = 7;
  int nDphiBinsLarge  = 4;

  /*
  MakeDefaultBinning( m_varDphiBinning, m_varDphiRebinnedBinning,
		      nDphiBinsLarge, nDphiBinsMedium, nDphiBinsSmall );
  */
  
  MakeLinearBinning ( m_varDphiBinning, m_varDphiRebinnedBinning, nDphiBinsLarge );
  // m_varDphiRebinnedBinning = m_varDphiBinning;
  
  m_nVarDphiBins         = m_varDphiBinning.size()         - 1;  
  m_nVarDphiRebinnedBins = m_varDphiRebinnedBinning.size() - 1;

  m_dPhiZoomLowBin  = nDphiBinsLarge + 1;
  m_dPhiZoomHighBin = m_nVarDphiBins;
  
  m_dPhiRebinnedZoomLowBin  = nDphiBinsLarge + 1;
  m_dPhiRebinnedZoomHighBin = m_nVarDphiRebinnedBins;

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
  
  m_dPhiWidthMin = 0.0;
  m_dPhiWidthMax = 0.6;

  m_dPhiYieldMin = 1.0e-3;
  m_dPhiYieldMax = 1.0;

  // default values 0->2 are in
  // SetHStyleRatio.
  // This is for "nonstandard"
  m_ratioMin = 0.5;
  m_ratioMax = 1.5;
  
  //========= common tools ========
  anaTool   = new CT::AnalysisTools();
  drawTool  = new CT::DrawTools();
  styleTool = new CT::StyleTools();

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
  // name for "All" histos
  // this is either merged JZN or Data from Triggers
  m_allName   = "All";

  m_spectName      = "spect";
  m_dPhiName       = "dPhi";
  m_widthName      = "width";
  m_yieldName      = "yield";
  
  m_recoName       = "reco";
  m_truthName      = "truth";
  m_pairedName     = "paired";
  m_unpairedName   = "unpaired";
    
  m_cFactorName    = "cFactor";
  m_respMatName    = "respMat";
  m_unfoldedName   = "unfolded";

  m_etaEffName     = m_sEta   + "_" + m_effName;
  m_ystarEffName   = m_sYstar + "_" + m_effName;
  
  m_etaSpectName        = m_sEta + "_" + m_spectName;
  
  m_ystarSpectName      = m_sYstar + "_" + m_spectName;
  m_ystarSpectRecoName  = m_ystarSpectName + "_" + m_recoName;
  m_ystarSpectTruthName = m_ystarSpectName + "_" + m_truthName;

  m_ystarSpectFineName      = m_ystarSpectName + "_" + m_sFine;
  m_ystarSpectFineRecoName  = m_ystarSpectFineName + "_" + m_recoName;
  m_ystarSpectFineTruthName = m_ystarSpectFineName + "_" + m_truthName;
 
  m_ystarSpectRespMatName      = m_ystarSpectName + "_" + m_respMatName;  
  m_ystarSpectCfactorsName     = m_ystarSpectName + "_" + m_cFactorName;
  m_ystarSpectUnfoldedName     = m_ystarSpectName + "_" + m_unfoldedName; 
  m_ystarSpectRecoUnfoldedName = m_ystarSpectRecoName + "_" + m_unfoldedName;
  
  m_dPhiRecoName   = m_dPhiName + "_" + m_recoName;
  m_dPhiTruthName  = m_dPhiName + "_" + m_truthName;

  m_dPhiCfactorsName    = m_dPhiName + "_" + m_cFactorName;;
  m_dPhiUnfoldedName    = m_dPhiName + "_" + m_unfoldedName;
  
  m_systematicsName     = "systematics";

  m_effName        = "eff";
  m_purityName     = "purity";
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
  
  // set names for various output files.
  // for data only, we just have one.
  // no need to add unfolding suffix.
  m_fName = m_dirOut + "/" + m_myOutName + "_" + m_labelOut;

  m_fNameRaw    = m_fName + "_" + m_sRaw;
  m_fNameRawUW  = m_fName + "_" + m_sRaw + "_" + m_unweightedFileSuffix;
  m_fNamePerf   = m_fName + "_" + m_sPerf;
  m_fNamePerfUF = m_fName + "_" + m_sPerf + "_" + m_unfoldingFileSuffix;
  m_fNamePhys   = m_fName + "_" + m_sPhys;
  m_fNamePhysUF = m_fName + "_" + m_sPhys + "_" + m_unfoldingFileSuffix;

  m_fNameDefRaw    = m_fNameRaw    + "_0.root";
  m_fNameDefRawUW  = m_fNameRawUW  + "_0.root";
  m_fNameDefPerf   = m_fNamePerf   + "_0.root";
  m_fNameDefPerfUF = m_fNamePerfUF + "_0.root";
  m_fNameDefPhys   = m_fNamePhys   + "_0.root";
  m_fNameDefPhysUF = m_fNamePhysUF + "_0.root";

  m_fNameRaw       += "_" + m_uncertSuffix + ".root";
  m_fNameRawUW     += "_" + m_uncertSuffix + ".root";
  m_fNamePerf      += "_" + m_uncertSuffix + ".root";
  m_fNamePerfUF    += "_" + m_uncertSuffix + ".root";
  m_fNamePhys      += "_" + m_uncertSuffix + ".root";
  m_fNamePhysUF    += "_" + m_uncertSuffix + ".root";

  m_fNameSYS = m_fName + "_" + m_systematicsFileSuffix + ".root";
  
  std::cout << "fNameRaw: " << m_fNameRaw << std::endl;
}

//---------------------------
// Fill Tree / Plot Controls
//---------------------------

void DiJetAnalysis::MakeResultsTogether(){

  // For Together "final"
  // Check if the directories exist.
  // If they don't, create them
  m_dirOutTogether = m_sOutput;
  anaTool->CheckWriteDir( m_dirOutTogether.c_str() );
  m_dirOutTogether += "/" + m_allName;
  anaTool->CheckWriteDir( m_dirOutTogether.c_str() );
  m_dirOutTogether += "/" + m_labelOut + "_" + m_uncertSuffix;
  anaTool->CheckWriteDir( m_dirOutTogether.c_str() );

  m_fNameTogether =
    m_dirOutTogether + "/" + m_myOutName + "_" +
    m_labelOut + "_" + m_uncertSuffix + ".root" ;

  TFile* fOut  = new TFile( m_fNameTogether.c_str() ,"recreate");

  MakeSpectTogether( fOut );
  //MakeFinalPlotsTogether( fOut, m_widthName );
  //MakeFinalPlotsTogether( fOut, m_yieldName );

  MakeDphiTogether ( fOut );
  
  std::cout << "DONE! Closing " << fOut->GetName() << std::endl;
  fOut->Close();
  std::cout << "......Closed  " << fOut->GetName() << std::endl;
}

void DiJetAnalysis::ProcessSystematics(){}

//---------------------------
//       Fill Tree
//---------------------------
void DiJetAnalysis::AddHistogram( TH1* h ){
  
  v_hists.push_back( h );
  h->Sumw2();
  styleTool->SetHStyle( h, 0 );

  TH1* h1 = dynamic_cast< TH1* >(h);
  if( h1 ){ h->SetNdivisions( 504, "X" ); }

  TH2* h2 = dynamic_cast< TH2* >(h);
  if( h2 ){
    h->SetNdivisions( 504, "X" );
    h->SetNdivisions( 504, "Y" );
  }

  TH3* h3 = dynamic_cast< TH3* >(h);
  if( h3 ){
    h->SetNdivisions( 504, "X" );
    h->SetNdivisions( 504, "YZ" );
  }
}

void DiJetAnalysis::AddHistogram( THnSparse* hn ){
  
  v_hns.push_back( hn );
  hn->Sumw2();
  hn->GetAxis(0)->SetTitle("#it{y}_{1}*");
  hn->GetAxis(1)->SetTitle("#it{y}_{2}*");
  hn->GetAxis(2)->SetTitle("#it{p}_{T}^{1}");
  hn->GetAxis(3)->SetTitle("#it{p}_{T}^{2}");
  if( hn->GetNdimensions() >= 5 ){
    hn->GetAxis(4)->SetTitle("|#Delta#phi|"); }
}


void DiJetAnalysis::SaveOutputsFromTree(){
  
  //----------------------------------------
  //  Close the input file, 
  //  write histos to output
  //----------------------------------------
  std::cout << "fName: " << m_fNameRaw << std::endl;
  TFile* fOut = new TFile( m_fNameRaw.c_str(),"RECREATE");
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
  return ( jetFwd || jetCent ) ? true : false;
}

bool DiJetAnalysis::GetDiJets( const std::vector< TLorentzVector >& v_jets, 
			       const TLorentzVector*& jet1,
			       const TLorentzVector*& jet2 ){
  
  for( auto& jet : v_jets ){
    if( !jet1 && jet.Pt() > 0 ){
      jet1 = &jet;
    } else if( jet1 && jet.Pt() > 0 ){
      jet2 = &jet;
      break;
    }
  }
  
  // make sure we have two jets
  return ( jet1 && jet2 ) ? true : false;
}

void DiJetAnalysis::AnalyzeSpectra( TH2* hSpect,
				    const std::vector< TLorentzVector >& vJets){

  // to fill a simple spectra over all jets.
  // if want to do more complicated things (such as in data)
  // do them in ProcessEvents()
  // loop over paired reco jets and fill
  for( auto& jet : vJets ){
    if( jet.Pt() <= 0 ){ continue; }
	
    double jetYstar  = GetYstar( jet );
    double jetPt     = jet.Pt()/1000.;
    
    double jetWeight = GetJetWeight( jet );

    hSpect->Fill(  jetYstar,  jetPt,  jetWeight );	

    // for pp fill both sides
    if( m_is_pPb ){ continue; }

    hSpect->Fill( -jetYstar,  jetPt,  jetWeight );	
  }
}

double DiJetAnalysis::AnalyzeDeltaPhi( THnSparse* hn,
				       const std::vector< TLorentzVector >& v_jets ){
  
  const TLorentzVector* jet1 = NULL; const TLorentzVector* jet2 = NULL;

  if( !GetDiJets( v_jets, jet1, jet2 ) ){ return -1; }

  double jet1_pt    = jet1->Pt()/1000.;
  double jet1_ystar = GetYstar( *jet1 );

  double jet2_pt    = jet2->Pt()/1000.;
  double jet2_ystar = GetYstar( *jet2 );
  
  double deltaPhi   = anaTool->DeltaPhi( *jet2, *jet1 );
  double jetWeight  = GetJetWeight( *jet1 );

  std::vector< double > p{ jet1_ystar, jet2_ystar, jet1_pt, jet2_pt, deltaPhi, 0.5 };

  hn->Fill( &p[0], jetWeight );

  // for pp, fill twice. once for each side since
  // it is symmetric in pp. For pPb, continue
  if( m_is_pPb ){ return deltaPhi; }

  p[0] = -jet1_ystar;  
  p[1] = -jet2_ystar;
  hn->Fill( &p[0], jetWeight );

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

  // fill histograms, for pp do both sides
  if( jet1 ){
    h->Fill( jet1->Eta(), jet1->Pt()/1000., weight );
    if( !m_is_pPb ){ h->Fill( -jet1->Eta(), jet1->Pt()/1000., weight ); }
  }
  if( jet2 ){
    h->Fill( jet2->Eta(), jet2->Pt()/1000., weight );
    if( !m_is_pPb ){ h->Fill( -jet2->Eta(), jet2->Pt()/1000., weight ); }
  }
} 

double DiJetAnalysis::GetYstar( const TLorentzVector& jet )
{ return m_is_pPb ? jet.Rapidity() + constants::BETAZ : jet.Rapidity(); }

bool DiJetAnalysis::IsForwardDetector( const double& eta ){

  return ( std::abs( eta ) < constants::FETAMAX &&
	   std::abs( eta ) > constants::FETAMIN ) ? true : false;
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

void DiJetAnalysis::NormalizeDeltaPhi( TH1* hIn, TH1* hNorm,
				       double xBinCenter, bool comb ){

  if( !hIn->GetEntries() ){ return; }

  // for yields, have to rescale later.
  // need this for comb subtraction, however
  hIn->Scale( 1., "width" );
  
  if( comb ){ anaTool->SubtractCombinatoric( hIn ); }

  if( !hNorm ){
    hIn->Scale( 1./hIn->Integral() ) ;
  } else {
    int      xBin = hNorm->FindBin( xBinCenter );
    double  nJets = hNorm->GetBinContent( xBin );
    std::cout << hIn->GetName() << " " << hIn->GetBinContent(6) << " " << nJets << std::endl;
    hIn->Scale( 1./nJets );
    std::cout << hIn->GetBinContent(6) << std::endl;
  }
}

double DiJetAnalysis::GetJetWeight( const TLorentzVector& jet )
{ return 1; }

TH1* DiJetAnalysis::CombineSamples( std::vector< TH1* >& vSampleHin,
				    const std::string& name ){ return NULL; }  

TH2* DiJetAnalysis::CombineSamples( std::vector< TH2* >& vSampleHin,
				    const std::string& name ){ return NULL; }  

THnSparse* DiJetAnalysis::CombineSamples( std::vector< THnSparse* >& vSampleHin,
					  const std::string& name ){ return NULL; }  


void DiJetAnalysis::GetSpectraLabels( std::string& axisLabel, std::string& axisLabelTex,
				      const std::string& name ){

  if( name.find( m_sEta ) != std::string::npos ){
    axisLabel    = "Eta";
    axisLabelTex = "#eta";
  } else if( name.find( m_sYstar ) != std::string::npos ){
    axisLabel    = "Ystar1";
    axisLabelTex = "#it{y}_{1}*";
  } else {
    axisLabel    = "";
    axisLabelTex = "";
  }
}


// hM is measured, hC is correction factor
TH1* DiJetAnalysis::BinByBinUnfolding( TH1* hM, TH1* hC ){
  
  TH1* hUnf = static_cast<TH1D*>( hM->Clone("h_unfolded") );
  hUnf->Reset();
  
  for( int xBin = 1; xBin <= hM->GetNbinsX(); xBin++ ){
    int cBin  = hC->FindBin( hM->GetBinCenter( xBin ) );
    
    double vM = hM->GetBinContent( xBin );
    double vC = hC->GetBinContent( cBin );
    
    if( vM == 0 || vC == 0 ){ continue; }
    
    double eM = hM->GetBinError( xBin );
    double eC = hC->GetBinError( cBin );
    
    // correction factor;
    double newDphi = vM * vC;
    // error on correction factor    

    double newDphiError =  newDphi *
      std::sqrt( std::pow( eM / vM, 2) +
		 std::pow( eC / vC, 2) ) ; 

    // newDphiError = eM * vC; 
    
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

  TFile* fInDef = NULL;
  
  while( file >> number )
    { v_uc.push_back( number ); }
  
  for( auto uc : v_uc ){
    m_uncertComp = uc;
    Initialize();
    
    if( !m_fNamePhysUF.compare( m_fNameDefPhysUF ) ){
      std::cout << "Found default sample: " << m_fNamePhysUF << std::endl;
      fInDef = TFile::Open( m_fNamePhysUF.c_str() );
    } else {
      std::cout << "Adding to systematics: " << m_fNamePhysUF << std::endl;
      mFinUC[ uc ] =  TFile::Open( m_fNamePhysUF.c_str() ) ;
    }
  }

  // set back to default
  m_uncertComp = 0;
  Initialize();
  
  return fInDef;
}

void DiJetAnalysis::MakeDefaultBinning( std::vector< double >& varBinning,
					std::vector< double >& varRebinnedBinning,
					int nLargeBins,
					int nMediumBins,
					int nSmallBins ){

  std::cout << "DEFAULT BINNING" << std::endl;

  int nVarBins = nLargeBins + nMediumBins + nSmallBins;
  
  int dPhiBinsSmallFactor  = 1;
  int dPhiBinsMediumFactor = 2;

  double varBinStart   = 2 * constants::PI / 3;
  double largeBinWidth = varBinStart / nLargeBins;
  
  for( int largeBin = 0; largeBin < nLargeBins; largeBin++ ){
    varBinning.push_back( largeBin * largeBinWidth );
    varRebinnedBinning.push_back( varBinning.back() );
  }

  double smallBinWidth =
    ( constants::PI - varBinStart ) /
    ( nSmallBins  * dPhiBinsSmallFactor  +
      nMediumBins * dPhiBinsMediumFactor );
  
  double mediumBinStart =
    nLargeBins * largeBinWidth;
  for( int mediumBin = 0; mediumBin < nMediumBins; mediumBin++ )
    { varBinning.push_back
	( mediumBin * smallBinWidth * dPhiBinsMediumFactor + mediumBinStart ); }

  double smallBinStart =
    mediumBinStart +  nMediumBins * smallBinWidth * dPhiBinsMediumFactor;
  for( int smallBin = 0; smallBin <= nSmallBins; smallBin++ )
    { varBinning.push_back
	( smallBin * smallBinWidth * dPhiBinsSmallFactor + smallBinStart ); }


  for( int i = nLargeBins; i <= nVarBins; i += 2 ){
    varRebinnedBinning.push_back( varBinning[ i ] );
  }
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
  
  for( int i = 1; i <= nVarBins; i++ ){
    double width = finalWidth - ( nVarBins - i ) * dWidthPerBin;
    varBinning.push_back( varBinning.back() + width );
  }
  // for now, do a rebin of 2 in middle
  // be careful total bins is a multiple of 2
  for( int i = nLargeBins + 2; i < nVarBins + nLargeBins - 4; i += 2 ){
    varRebinnedBinning.push_back( varBinning[ i ] );
  }
  
  // leave the last 4 as is => i++
  // rebin last 4 bins => i += 2
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

void DiJetAnalysis::MakeEtaPhiPtMap( std::vector< TH2* >& vSampleMaps,
				     const std::vector< std::string>& vLabels,
				     const std::string& name ){

  uint nSamples = vSampleMaps.size();
  
  TCanvas c_map( "c_map", "c_map", 800, 600 );

  for( uint iG = 0; iG < nSamples; iG++ ){
    std::string label = vLabels[iG];

    if( label.compare( m_allName ) ){ continue; }

    TH1* h = vSampleMaps[ iG ];
    
    h->Draw("col");
    styleTool->SetHStyle( h, 0 );
    DrawAtlasRight();  
    SaveAsAll( c_map, h->GetName() );
    h->Write();
  }
}

void DiJetAnalysis::MakeSpectra( std::vector< TH2* >& vSampleSpect,
				 const std::vector< std::string>& vLabels,
				 const std::string& name ){
    
  if( !vSampleSpect.size() ){ return; }

  std::string axisLabel, axisLabelTex;
  GetSpectraLabels( axisLabel, axisLabelTex, name );

  bool isEta = name.find( m_sEta ) != std::string::npos ? true : false;

  // use this as reference because
  // it should be in every file
  TH2*  hRef = vSampleSpect[0];
  int nXbins = hRef->GetNbinsX();

  uint nSamples = vSampleSpect.size();
  
  std::vector< std::vector< TH1* > > vSpect;
  std::vector< std::vector< TH1* > > vSpectCounts;
  vSpect      .resize( nSamples );
  vSpectCounts.resize( nSamples );

  std::string yAxisTitle = "dN/d#it{p}_{T} [GeV]";

  double max = -1;
  
  for( uint iG = 0; iG < nSamples; iG++){

    std::string label = vLabels[iG];

    for( int xBin = 1; xBin <= nXbins; xBin++ ){
      
      double xMin, xMax;
      anaTool->GetBinRange
	( hRef->GetXaxis(), xBin, xBin, xMin, xMax );
      
      std::string hTag = anaTool->GetName( xMin, xMax, axisLabel);
       
      TH1* hSpectCounts =
	vSampleSpect[iG]->
	ProjectionY( Form("h_%s_%s_%s_%s",
			  name.c_str(), m_sCounts.c_str(), label.c_str(), hTag.c_str() ),
		     xBin, xBin );
      
      TH1* hSpect = static_cast< TH1D* >
	( hSpectCounts->Clone
	  ( Form("h_%s_%s_%s",
		 name.c_str(), label.c_str(), hTag.c_str() ) ) );
      
      hSpect->SetTitle( anaTool->GetLabel( xMin, xMax, axisLabelTex ).c_str() );
      hSpect->SetYTitle( yAxisTitle.c_str() );
      
      vSpect      [iG].push_back( hSpect );
      vSpectCounts[iG].push_back( hSpectCounts );

      hSpect->GetXaxis()->SetRangeUser( m_varPtBinning.front(), m_varPtBinning.back() );
     
      // scale by width of bins to get dN/dpT
      hSpect->Scale( 1., "width" );
      
      hSpect->Write();
      hSpectCounts->Write();
      
      // get min max from the final histograms
      if( label.compare( m_allName ) ){ continue; }
      if( max < hSpect->GetMaximum() ){ max = hSpect->GetMaximum(); }
      
    } // end loop over xBin
  } // end loop over iG

  // set maxima globally for all spectra hists.
  // easier to compare. Set on log scale.

  max = anaTool->GetLogMaximum( max );
  
  double lX0, lY0, lX1, lY1;
  
  //------------------------------------------------
  //------- For Each xAxis Bin, draw an IGs --------
  //------------------------------------------------

  if( m_is_pPb ){ lX0 = 0.45; lY0 = 0.54; lX1 = 0.76; lY1 = 0.67; }
  else          { lX0 = 0.25; lY0 = 0.20; lX1 = 0.45; lY1 = 0.40; }
  
  for( int iX = 0; iX < nXbins; iX++ ){
    int    xBin = iX + 1;
    double xMin, xMax;
    anaTool->GetBinRange( hRef->GetXaxis(), xBin, xBin, xMin, xMax );
    double xCenter = hRef->GetXaxis()->GetBinCenter ( xBin );

    // for pPb, dont draw at anything above -3.2
    if( isEta && m_is_pPb && xCenter > -constants::FETAMIN ){ continue; }
     
    std::string cName  = anaTool->GetName ( xMin, xMax, axisLabel    );
    std::string cLabel = anaTool->GetLabel( xMin, xMax, axisLabelTex );
    
    TCanvas c( cLabel.c_str(), cLabel.c_str(), 800, 600 );
    c.SetLogy();
    
    TLegend leg( lX0, lY0, lX1, lY1 );
    styleTool->SetLegendStyle( &leg, 0.7 );
    
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

      TH1* h = vSpect[iG][iX];
           
      if( !h->GetEntries() && !isMb && !isAll ){ continue; }
      else if( !h->GetEntries() && !isMb && !isAll ){ continue; }

      h->GetXaxis()->SetRangeUser( m_varPtBinning.front(), m_varPtBinning.back() );
      
      if( iG == nSamples ){ style = 0; }
      styleTool->SetHStyle( h, style++ );
      h->Draw("epsame");
      h->SetMinimum( m_ptSpectYaxisMin );
      h->SetMaximum( max );
      leg.AddEntry( h, label.c_str() );
    } // end loop over iG
    
    leg.Draw("same");

    DrawAtlasRight();    
    drawTool->DrawRightLatex( 0.87, 0.73, cLabel );

    SaveAsAll( c, Form("%s_%s", name.c_str(), cName.c_str() ) );
  } // end loop over iX
  
  //------------------------------------------------
  //--------- For each IG, draw xAxisBins ----------
  //------------------------------------------------

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
    styleTool->SetLegendStyle( &leg, 0.7 );

    int style = 0;
    for( int iX = 0; iX < nXbins; iX++ ){
      int       xBin = iX + 1;
      double xCenter = hRef->GetXaxis()->GetBinCenter ( xBin );

      // for pPb, dont draw at anything above -3.2
      if( m_is_pPb && xCenter > -constants::FETAMIN ){ continue; }

      // for pp, dont draw central triggers below -3.2
      // or forward triggers above -3.2
      // for mb, and total draw everything
      bool isMb  = label.find("_mb_") != std::string::npos
	? true : false;
      bool isAll = !label.compare( m_allName )
	? true : false;

      TH1* h = vSpect[iG][iX];
      
      if( h->GetEntries() && !isMb && !isAll ){ continue; }
      else if( h->GetEntries() && !isMb && !isAll ){ continue; }

      h->GetXaxis()->SetRangeUser( m_varPtBinning.front(), m_varPtBinning.back() );
     
      styleTool->SetHStyle( h, style++ );
      h->Draw("epsame");
      h->SetMinimum( m_ptSpectYaxisMin );
      h->SetMaximum( max );
      leg.AddEntry( h, h->GetTitle() );
    } // end loop over iX
    leg.Draw("same");

    DrawAtlasRight();    

    SaveAsAll( c, Form("%s_%s", name.c_str(), cName.c_str() ) );
  } // end loop over iG  

  // delete
  for( uint iG = 0; iG < nSamples; iG++ ){
    for( int iX = 0; iX < nXbins; iX++ ){
      delete vSpect      [iG][iX];
      delete vSpectCounts[iG][iX];
    }
  }
} 

TH2* DiJetAnalysis::UnfoldSpectra( TFile* fInData, TFile* fInMC,
				   const std::string& hUnfoldedName ){

  std::vector< TH1* > vSpect;
  std::vector< TH1* > vCfactors;
  std::vector< TH2* > vRespMat;
  std::vector< TH1* > vUnfolded;
  
  std::string unfoldingMCLabel = m_mcTypeLabel;

  std::string axisLabel, axisLabelTex;
  GetSpectraLabels( axisLabel, axisLabelTex, hUnfoldedName );
  
  std::string measuredName, recoName, truthName, respMatName;
  std::string unfoldedLabel, typeMeasured;
  GetSpectUnfoldingInfo( measuredName, recoName, truthName, 
			 respMatName, unfoldedLabel, typeMeasured );

  // make final output histo for UnfSpect
  // leave for now, should be more general not hardcoded
  TH2* hUnfoldedCountsAll =
    new TH2D( Form( "h_%s_%s", hUnfoldedName.c_str(), m_allName.c_str() ),
	      ";#it{y}_{1}*;#it{p}_{T1} [GeV]",
	      m_nVarYstarBins, 0, 1,
	      m_nVarPtBinsUfOf, 0, 1 );
  hUnfoldedCountsAll->GetXaxis()->Set( m_nVarYstarBins    , &( m_varYstarBinning [0] ) );
  hUnfoldedCountsAll->GetYaxis()->Set( m_nVarPtBinsUfOf, &( m_varPtBinningUfOf[0] ) );

  TAxis* xAxisUnfolded = hUnfoldedCountsAll->GetXaxis();

  for( int xBin = 1; xBin <= hUnfoldedCountsAll->GetNbinsX(); xBin++ ){

    double xLow, xUp;
    anaTool->GetBinRange
      ( xAxisUnfolded, xBin, xBin, xLow, xUp );

    std::string hTag = anaTool->GetName( xLow, xUp, axisLabel );
    
    TH1* hMeasuredCounts = static_cast< TH1D* >
      ( fInData->Get( Form("h_%s_%s_%s_%s", measuredName.c_str(), m_sCounts.c_str(),
			   m_allName.c_str(), hTag.c_str() ) ) );
    TH1* hMeasured = static_cast< TH1D* >
      ( fInData->Get( Form("h_%s_%s_%s", measuredName.c_str(),
			   m_allName.c_str(), hTag.c_str() ) ) );
    TH1* hReco     = static_cast< TH1D* >
      ( fInMC->Get( Form("h_%s_%s_%s", recoName.c_str(),
			 m_allName.c_str(), hTag.c_str() ) ) );
    TH1* hTruth    = static_cast< TH1D* >
      ( fInMC->Get( Form("h_%s_%s_%s", truthName.c_str(),
			 m_allName.c_str(), hTag.c_str() ) ) );

    
    styleTool->SetHStyle( hMeasuredCounts, 1 );
    styleTool->SetHStyle( hMeasured, 1 );
    styleTool->SetHStyle( hTruth   , 2 );

    vSpect.push_back( hMeasuredCounts );
    vSpect.push_back( hMeasured );
    vSpect.push_back( hTruth    );
    
    TH1* hCfactors = static_cast<TH1D*>
      ( fInMC->Get( Form( "h_%s_%s_%s", m_ystarSpectCfactorsName.c_str(),
			  m_allName.c_str(), hTag.c_str())));
    
    TH2* hRespMat  = static_cast< TH2D* >
      ( fInMC->Get( Form("h_%s_%s_%s", respMatName.c_str(),
			 m_allName.c_str(), hTag.c_str() ) ) );

    styleTool->SetHStyle( hCfactors, 0 );
    styleTool->SetHStyle( hRespMat , 0 );

    vCfactors.push_back( hCfactors );
    vRespMat .push_back( hRespMat  );

    /*
    // setup unfold
    RooUnfoldResponse response( hMeasured, hTruth, hRespMat );
    RooUnfoldBayes      unfold( &response, hMeasured, 2, false );
    response.UseOverflow( kFALSE );
    // RooUnfoldBinByBin   unfold( &response, hMeasured );
    
    unfold.SetVerbose(0);

    // unfold measured
    TH1* hUnfolded = (TH1D*)unfold.Hreco();
	  
    */

    TH1* hUnfoldedCounts = BinByBinUnfolding( hMeasuredCounts, hCfactors );
    hUnfoldedCounts->SetName
      ( Form( "h_%s_%s_%s_%s", hUnfoldedName.c_str(), m_sCounts.c_str(),
	      m_allName.c_str(), hTag.c_str()));
   
    // Clone the counts, then normalize
    TH1* hUnfolded = static_cast< TH1D* >
      ( hUnfoldedCounts->Clone( Form( "h_%s_%s_%s", m_ystarSpectUnfoldedName.c_str(),
				      m_allName.c_str(), hTag.c_str())));
    styleTool->SetHStyle( hUnfolded, 0 );
    vSpect.push_back( hUnfolded );

    // fill final 2d unfolded spectra histo
    for( int yBin = 1; yBin <= hUnfoldedCountsAll->GetNbinsY(); yBin++ ){
      double val      = hUnfoldedCounts->GetBinContent( yBin );
      double valError = hUnfoldedCounts->GetBinError  ( yBin );

      hUnfoldedCountsAll->SetBinContent( xBin, yBin, val      );
      hUnfoldedCountsAll->SetBinError  ( xBin, yBin, valError );
    }
    
    // normalize to bin width
    hUnfolded->Scale( 1., "width" );

    hUnfolded->Write();

    // Now Draw everything.
    TCanvas c( hUnfolded->GetName(), hUnfolded->GetName(), 800, 800 );
    TPad pad1("pad1", "", 0.0, 0.35, 1.0, 1.0 );
    pad1.SetBottomMargin(0.0);
    pad1.Draw();
    TPad pad2("pad2", "", 0.0, 0.0, 1.0, 0.34 );
    pad2.SetTopMargin(0.05);
    pad2.SetBottomMargin(0.25);
    pad2.Draw();

    pad1.cd();
    pad1.SetLogy();
    
    hMeasured->SetTitle( "" );
    hTruth   ->SetTitle( "" );
    hUnfolded->SetTitle( "" );    

    double maximum =
      hMeasured->GetMaximum() > hUnfolded->GetMaximum() ? 
      hMeasured->GetMaximum() : hUnfolded->GetMaximum(); 
    
    double minimum =
      hMeasured->GetMinimum() < hUnfolded->GetMinimum() ? 
      hMeasured->GetMinimum() : hUnfolded->GetMinimum(); 

    maximum = anaTool->GetLogMaximum( maximum );
    minimum = anaTool->GetLogMinimum( minimum );
    
    hMeasured->SetMaximum( maximum );
    hTruth   ->SetMaximum( maximum );
    hUnfolded->SetMaximum( maximum );

    hMeasured->SetMinimum( minimum );
    hTruth   ->SetMinimum( minimum );
    hUnfolded->SetMinimum( minimum );
    
    TLegend leg( 0.71, 0.50, 0.90, 0.70 );
    styleTool->SetLegendStyle( &leg );

    leg.AddEntry( hMeasured, Form( "%s_{Reco}" , typeMeasured.c_str() ) );
    leg.AddEntry( hUnfolded, Form( "%s_{UF}", typeMeasured.c_str() ) );

    TH1* hRdataMC = NULL;
    
    // only draw truth when there is MC
    // otherwise, draw reco unfolded MC
    // scaled to unfolded data.
    if( !m_isData ){
      leg.AddEntry( hTruth, "MC_{Truth}" );
      hTruth->Draw( "histo same" );
    } else {
      double integralData = hUnfolded->Integral( 2, m_nVarPtBinsUfOf - 1 );
      double integralMC   = hReco    ->Integral( 2, m_nVarPtBinsUfOf - 1 );

      double scalingFactor = integralData / integralMC;

      hReco->Scale( scalingFactor );
      hReco->SetFillColor( kAzure - 9 );
      
      leg.AddEntry( hReco, "MC_{UF}", "f" );
      hReco->Draw ( "hist" );
    }

    hMeasured->Draw( "ep same" );
    hUnfolded->Draw( "ep same" );
    
    leg.Draw();
    
    DrawAtlasRight();    
    drawTool->DrawRightLatex
      ( 0.45, 0.15, anaTool->GetLabel( xLow, xUp, axisLabelTex ) );

    // make ratios and draw cfactors
    pad2.cd();

    double legXmin = m_isData ? 0.35 : 0.55;
    double legXmax = m_isData ? 0.85 : 0.90;
    
    TLegend legR( legXmin, 0.82, legXmax, 0.92 );
    styleTool->SetLegendStyle( &legR, 0.9 );
    legR.SetNColumns(3);
	  
    // make sure to set range to rebinned.
    // the c-factors are from rebinned distributions.
    styleTool->SetHStyleRatio( hCfactors );
    hCfactors->SetYTitle( "Ratio" );
    hCfactors->SetXTitle( "#it{p}_{T} [GeV]" );
    hCfactors->SetTitleOffset( 2.3, "x" );

    if( m_isData ){
      hCfactors->SetMaximum( 1.75 );
      hCfactors->SetMinimum( 0.50 );
    } else{
      hCfactors->SetMaximum( 1.15 );
      hCfactors->SetMinimum( 0.55 );
    }
    hCfactors->Draw("ep");

    TH1* hR = static_cast< TH1D* >
      ( hUnfolded->
	Clone( Form( "h_%s_%s_%s_%s_%s", m_sRatio.c_str(),
		     measuredName.c_str(), m_unfoldedName.c_str(),
		     m_allName.c_str(), hTag.c_str())));
    styleTool->SetHStyleRatio( hR, 1 );
    hR->Divide( hMeasured );

    hR->GetXaxis()->SetRangeUser( m_varPtBinning.front(), m_varPtBinning.back() );
      
    hR->Draw( "ep same" );
	  
    legR.AddEntry( hCfactors, "CF" );
    legR.AddEntry( hR , Form( "%s_{UF}/%s_{Reco}",
			      typeMeasured.c_str(), typeMeasured.c_str() ) );

    // make ratio of data to reco
    if( m_isData ){
      hRdataMC = static_cast< TH1D* >
	( hMeasured->
	  Clone( Form( "h_%s_%s_%s_%s_%s", m_sRatio.c_str(),
		       measuredName.c_str(), m_recoName.c_str(),
		       m_allName.c_str(), hTag.c_str())));
      styleTool->SetHStyleRatio( hRdataMC, 2 );
      hRdataMC->Divide( hReco );
      hRdataMC->Draw( "ep same" );
      legR.AddEntry( hRdataMC, "Data_{UF}/MC_{UF}" );
    }
    
    legR.Draw();

    styleTool->HideAxis( hMeasured, "x" );
    styleTool->HideAxis( hTruth   , "x" );

    double xMin = m_varPtBinning.front();
    double xMax = m_varPtBinning.back();
	  
    TLine line( xMin, 1, xMax, 1 );
    line.SetLineWidth( 2 );
    line.Draw();

    TLine lineP25( xMin, 1.25, xMax, 1.25 );
    lineP25.SetLineStyle( 2  );
    lineP25.SetLineColor( 12 );
    lineP25.SetLineWidth( 2  );
    lineP25.Draw();
    
    TLine lineN25( xMin, 0.75, xMax, 0.75 );
    lineN25.SetLineStyle( 2  );
    lineN25.SetLineColor( 12 );
    lineN25.SetLineWidth( 2  );
    lineN25.Draw();
    
    // MUT is Measuerd Unfolded Truth
    SaveAsAll( c, Form( "h_%s_%s_MUT_%s",
			m_ystarSpectUnfoldedName.c_str(),
			m_allName.c_str(), hTag.c_str()));
  }

  for( auto& h : vSpect    ){ delete h; }
  for( auto& r : vRespMat  ){ delete r; }
  for( auto& r : vCfactors ){ delete r; }
  for( auto& u : vUnfolded ){ delete u; }

  hUnfoldedCountsAll->Write();
  
  return hUnfoldedCountsAll;
}

void DiJetAnalysis::MakeDeltaPhi( std::vector< THnSparse* >& vhn,
				  const std::vector< std::string >& vLabels,
				  const std::string& dPhiName,
				  TFile* fInMCSpect, 
				  const std::string& spectName ){

  // lines to be drawn along axis3. this is
  // x-axis that widths are plotted as function of 
  double xMin = m_varYstarBinning.front();
  double xMax = m_varYstarBinning.back();
  
  TLine line( xMin, 1, xMax, 1 );
  line.SetLineWidth( 2 );

  TLine lineP25( xMin, 1.25, xMax, 1.25 );
  lineP25.SetLineStyle( 2  );
  lineP25.SetLineColor( 12 );
  lineP25.SetLineWidth( 2  );
	  
  TLine lineN25( xMin, 0.75, xMax, 0.75 );
  lineN25.SetLineStyle( 2  );
  lineN25.SetLineColor( 12 );
  lineN25.SetLineWidth( 2  );
  
  std::vector< std::string > v_listBadHist;
  std::vector< std::string > v_listBadFits;
 
  TH1D* h_dPhiChiS = new TH1D( Form( "h_dPhiChiS_%s", dPhiName.c_str() ),
			       ";#Chi^{2}/NDF;Count", 20, 0, 5 );
  TH1D* h_dPhiProb = new TH1D( Form( "h_dPhiProb_%s", dPhiName.c_str() ),
			       ";Probability;Count", 20, 0, 1 );

  TH1D* h_dPhiChiS2 = new TH1D( Form( "h_dPhiChiS2_%s", dPhiName.c_str() ),
			       ";#Chi^{2}/NDF;Count", 20, 0, 5 );
  TH1D* h_dPhiProb2 = new TH1D( Form( "h_dPhiProb2_%s", dPhiName.c_str() ),
			       ";Probability;Count", 20, 0, 1 );
 
  h_dPhiChiS->Sumw2();
  h_dPhiProb->Sumw2();

  styleTool->SetHStyle( h_dPhiChiS, 0 );
  styleTool->SetHStyle( h_dPhiProb, 0 );

  styleTool->SetHStyle( h_dPhiChiS2, 1 );
  styleTool->SetHStyle( h_dPhiProb2, 1 );

  TLegend legChiSProb( 0.60, 0.7, 0.74, 0.8 );
  styleTool->SetLegendStyle( &legChiSProb );
  
  legChiSProb.AddEntry( h_dPhiChiS , Form("Fit > %2.1f", m_dPhiFittingMin  ) );
  legChiSProb.AddEntry( h_dPhiChiS2, Form("Fit > %2.1f", m_dPhiFittingMinB ) );
  
  std::vector< TH1* > vHdPhi;
  std::vector< TH1* > vDphiWidths;
  std::vector< TH1* > vDphiYields;
  std::vector< TH1* > vSpect;
  std::vector< TF1* > vFits;
  std::vector< TH1* > vWR;
  
  // ---- loop over group  ----
  // ---- (jzn or trigger) ----
  for( uint iG = 0; iG < vhn.size(); iG++ ){      
    THnSparse* hn = vhn[iG];
  
    std::string label = vLabels[iG];

    if( label.compare( m_allName ) ){ continue; }
    
    TAxis* axis0 = hn->GetAxis( m_dPP->GetAxisI(0) ); int nAxis0Bins = axis0->GetNbins();
    TAxis* axis1 = hn->GetAxis( m_dPP->GetAxisI(1) ); int nAxis1Bins = axis1->GetNbins();
    TAxis* axis2 = hn->GetAxis( m_dPP->GetAxisI(2) ); int nAxis2Bins = axis2->GetNbins();
    TAxis* axis3 = hn->GetAxis( m_dPP->GetAxisI(3) ); int nAxis3Bins = axis3->GetNbins();
    
    int fAxisI   = m_dPP->GetAxisI(3);
    
    int dPhiAxisI = hn->GetNdimensions() - 1;
   
    int nCol = 5, nRow = 6;
    
    TCanvas cAll( "cAll", "cAll", 4800, 3600 );
    cAll.Divide( nCol, nRow, 0, 0 );

    int cAllPadI = 0;
    
    // ---- loop over ystars ----
    for( int axis0Bin = 1; axis0Bin <= nAxis0Bins; axis0Bin++ ){
      // set ranges
      axis0->SetRange( axis0Bin, axis0Bin );
      
      double axis0Low, axis0Up;
      anaTool->GetBinRange
	( axis0, axis0Bin, axis0Bin, axis0Low, axis0Up );
	
      // FIX!!!
      TH1* hSpectCounts = static_cast< TH1D* >
	( fInMCSpect->Get( Form( "h_%s_%s_%s_%s", spectName.c_str(), m_sCounts.c_str(),
				 m_allName.c_str(), anaTool->GetName
				 ( axis0Low, axis0Up, m_dPP->GetAxisName(0) ).c_str() ) ) );
      vSpect.push_back( hSpectCounts );
      
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
	TCanvas cWidths( "cWidths", "cWidths",800,600);
	TLegend legW( 0.55, 0.2, 0.89, 0.38 );
	styleTool->SetLegendStyle( &legW );

	std::vector< TH1* > vDphiYieldsTemp;
	TCanvas cYields( "cYields", "cYields",800,600);
	TLegend legY( 0.55, 0.55, 0.89, 0.73 );
	styleTool->SetLegendStyle( &legY );
	
	int style = 0;
	
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
	  axis3->SetRange( 0, 0 );

	  std::string hNameW = "h_" + dPhiName + "_" + m_widthName + "_" + label + "_" + hTagW;
	  std::string hNameY = "h_" + dPhiName + "_" + m_yieldName + "_" + label + "_" + hTagW;
	  
	  TH1* hDphiWidths = hn->Projection( fAxisI );
	  hDphiWidths->Reset();
	  hDphiWidths->SetName( hNameW.c_str() );
	  hDphiWidths->SetYTitle( "RMS (#pi - #Delta#phi)" );
	  hDphiWidths->SetMarkerSize( hDphiWidths->GetMarkerSize() * 1.5 );
	  styleTool->SetHStyle( hDphiWidths, style );
	  vDphiWidthsTemp.push_back( hDphiWidths );
	  vDphiWidths    .push_back( hDphiWidths );

	  TH1* hDphiWidths2 = hn->Projection( fAxisI );
	  hDphiWidths2->Reset();
	  hDphiWidths2->SetName( Form( "%s_2", hNameW.c_str() ) );
	  hDphiWidths2->SetYTitle( "RMS (#pi - #Delta#phi)" );
	  hDphiWidths2->SetMarkerSize( hDphiWidths2->GetMarkerSize() * 2 );
	  styleTool->SetHStyle( hDphiWidths2, 5 + style );
	  hDphiWidths2->SetMarkerColor( kBlue );
	  hDphiWidths2->SetLineColor  ( kBlue );
	  vDphiWidths.push_back( hDphiWidths2 );
	  
	  TH1* hDphiWidthsStat = hn->Projection( fAxisI );
	  hDphiWidthsStat->Reset();
	  hDphiWidthsStat->SetName( Form( "%s_stat", hNameW.c_str() ) );
	  hDphiWidthsStat->SetYTitle( "RMS (#pi - #Delta#phi)" );
	  hDphiWidthsStat->SetMarkerSize( hDphiWidthsStat->GetMarkerSize() * 2 );
	  styleTool->SetHStyle( hDphiWidthsStat, 5 + style );
	  vDphiWidths.push_back( hDphiWidthsStat );

	  TH1* hDphiYields = hn->Projection( fAxisI );
	  hDphiYields->Reset();
	  hDphiYields->SetName( hNameY.c_str() );
	  hDphiYields->SetYTitle( "Pair Jet Yield Per Jet_{1}" );
	  hDphiYields->SetMarkerSize( hDphiYields->GetMarkerSize() * 1.5 );
	  styleTool->SetHStyle( hDphiYields, style );
	  vDphiYieldsTemp.push_back( hDphiYields );
	  vDphiYields    .push_back( hDphiYields );

	  style++;
	  
	  legW.AddEntry
	    ( hDphiWidths, anaTool->
	      GetLabel( axis2Low, axis2Up, m_dPP->GetAxisLabel(2) ) .c_str() );

	  legY.AddEntry
	    ( hDphiWidths, anaTool->
	      GetLabel( axis2Low, axis2Up, m_dPP->GetAxisLabel(2) ) .c_str() );

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
	    std::string hDphiName =
	      Form( "h_%s_%s_%s", dPhiName.c_str(), label.c_str(), hTag.c_str() );
	    std::string hDphiCountsName =
	      Form( "h_%s_%s_%s_%s",
		    dPhiName.c_str(), m_sCounts.c_str(), label.c_str(), hTag.c_str() );
	    std::string hYieldsName =
	      Form( "h_%s_%s_%s", m_yieldName.c_str(), label.c_str(), hTag.c_str() );

	    TH1* hDphiCounts = hn->Projection( dPhiAxisI );
	    hDphiCounts->SetName( hDphiCountsName.c_str() );
	    styleTool->SetHStyle( hDphiCounts, 0 );
	    vHdPhi.push_back( hDphiCounts );
	    
	    // save the per-jet normalized counts histogram.
	    TH1* hDphi = static_cast< TH1D* >( hDphiCounts->Clone( hDphiName.c_str() ) );
	    styleTool->SetHStyle( hDphi, 0 );
	    hDphi->SetYTitle("1/N_{jet_{1}} dN_{12}/d|#Delta#phi|");
	    vHdPhi.push_back( hDphi );
	    	    
	    // Combinatoric subtraction after scaling by width
	    // then normalize by jet spectra.
	    NormalizeDeltaPhi( hDphi, hSpectCounts, 0.5 * ( axis1Up + axis1Low ), true );
	    
	    // take final dPhi histogram, and redo
	    // bin width normalization to get per-jet
	    // yields after comb subtaction
	    TH1* hYields = static_cast< TH1D* >( hDphi->Clone( hYieldsName.c_str() ) );
	    styleTool->SetHStyle( hDphi, 0 );
	    hDphi->SetXTitle("|#Delta#phi|");
	    hDphi->SetYTitle("1/N_{jet_{1}} dN_{12}/d|#Delta#phi|");
	    hDphi->SetTitle( "" );
	    vHdPhi.push_back( hYields );

	    // undo the hDphi->Scale( 1., "width" ) part
	    anaTool->UndoWidthScaling( hYields );

	    TCanvas c( hDphi->GetName(), hDphi->GetName(), 800, 600 );
	    c.SetLogy();

	    hDphi->SetMinimum( m_dPhiLogMin );
	    hDphi->SetMaximum( m_dPhiLogMax );
	    hDphi->Draw("ep same x0");

	    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
	    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
	    TF1* fit  = NULL;
	    TF1* fit2 = NULL; 
	    
	    // now fit
	    if( m_fitDphiWC ){
	      fit  = anaTool->FitDphiWC( hDphi, m_dPhiFittingMin , m_dPhiFittingMax );
	      fit2 = anaTool->FitDphiWC( hDphi, m_dPhiFittingMinB, m_dPhiFittingMax );
	    } else {
	      fit  = anaTool->FitDphi  ( hDphi, m_dPhiFittingMin , m_dPhiFittingMax );
	      fit2 = anaTool->FitDphi  ( hDphi, m_dPhiFittingMinB, m_dPhiFittingMax );
	    }
	    
	    styleTool->SetHStyle( fit, 0 );
	    vFits.push_back( fit );

	    styleTool->SetHStyle( fit2, 0 );
	    fit2->SetLineColor( kRed );
	    fit2->SetName( Form( "%s_2", fit2->GetName() ) );
	    vFits.push_back( fit2 );
	    
	    // set range back to what it should be.
	    // it can be changed in the fit funtion
	    hDphi->GetXaxis()->SetRange( m_dPhiZoomLowBin, m_dPhiZoomHighBin );

	    // draw longer fit first;
	    fit2->Draw("same");
	    fit ->Draw("same");
	    
	    double chi2NDF  = fit->GetChisquare() / fit->GetNDF();
	    double prob     = fit->GetProb();

	    double chi2NDF2 = fit2->GetChisquare() / fit2->GetNDF();
	    double prob2    = fit2->GetProb();	    
	    
	    h_dPhiChiS->Fill( chi2NDF );
	    h_dPhiProb->Fill( prob    );

	    h_dPhiChiS2->Fill( chi2NDF2 );
	    h_dPhiProb2->Fill( prob2    );

	    if( prob < 0.05 ){
	      v_listBadHist.push_back( hDphi->GetName() );
	      v_listBadFits.push_back( fit  ->GetName() );
	    }
	    
	    if( m_fitDphiWC ){
	      drawTool->DrawLeftLatex
		( 0.60, 0.41, Form( "Prob=(%4.2f, #color[2]{%4.2f})",
				    prob, prob2 ) );
	      drawTool->DrawLeftLatex
		( 0.55, 0.33, Form( "#Chi^{2}/NDF=(%4.2f, #color[2]{%4.2f})",
				    chi2NDF, chi2NDF2 ) );
	      drawTool->DrawLeftLatex
		( 0.43, 0.25, Form( "Const =(%4.2e, #color[2]{%4.2e})",
				    fit->GetParameter(3), fit2->GetParameter(3) ) );
	    } else {
	      drawTool->DrawLeftLatex
		( 0.70, 0.33, Form( "Prob = %4.2f", prob ) );
	      drawTool->DrawLeftLatex
		( 0.65, 0.25, Form( "#Chi^{2}/NDF = %4.2f", chi2NDF ) );
	    }

	    DrawTopLeftLabels
	      ( m_dPP, axis0Low, axis0Up, axis1Low, axis1Up,
		axis2Low, axis2Up, axis3Low, axis3Up );

	    DrawAtlasRight();
	    
	    SaveAsAll( c, hDphi->GetName() );
	    fit->Write();
	    
	    hDphi->Write();
	    hDphiCounts->Write();
	    
	    // draw onto common canvas
	    cAll.cd( cAllPadI + 1 );
	    gPad->SetLogy();
	    hDphi->Draw();
	    fit->Draw("same");
	    fit2->Draw("same");
	    
	    DrawTopLeftLabels
	      ( m_dPP, axis0Low, axis0Up, axis1Low, axis1Up,
		axis2Low, axis2Up, axis3Low, axis3Up );

	    double dx = 0, dy = 0;
	    
	    if( cAllPadI == nCol - 1 ){ DrawAtlasRight(); }
	    if( cAllPadI % nCol > 0 ){
	      hDphi->SetYTitle("");
	      dx = 0.1;
	    } else {
	      hDphi->GetYaxis()->SetTitleOffset(3.4);
	      dx = 0.15;
	    }

	    if( cAllPadI / nCol <  nRow - 1 ){
	      hDphi->SetXTitle("");
	      dy = - 0.1;
	    } else {
	      hDphi->GetXaxis()->SetTitleOffset(4);
	    }

	    drawTool->DrawLeftLatex
	      ( 0.60 + dx, 0.33 + dy, Form( "Prob=(%4.2f, #color[2]{%4.2f})"
					    , prob, prob2 ) );
	    drawTool->DrawLeftLatex
	      ( 0.55 + dx, 0.25 + dy, Form( "#Chi^{2}/NDF=(%4.2f, #color[2]{%4.2f})"
				  , chi2NDF, chi2NDF2 ) );
	    /*
	    drawTool->DrawLeftLatex
	      ( 0.70 + dx, 0.33 + dy, Form( "Prob = %4.2f", prob ) );
	    drawTool->DrawLeftLatex
	      ( 0.65 + dx, 0.25 + dy, Form( "#Chi^{2}/NDF = %4.2f", chi2NDF ) );
	    */

	    cAllPadI++;
	    
	    if( fit->GetParameter(1) < 0 )
	      { continue; }

	    double   tau = fit->GetParameter(1);
	    double sigma = fit->GetParameter(2);

	    double   tauError = fit->GetParError(1);
	    double sigmaError = fit->GetParError(2);
	    
	    double yieldError = 0;
	    double yield      = hYields->IntegralAndError( 1, hDphi->GetNbinsX(), yieldError );
	 
	    double width      = std::sqrt( tau * tau + sigma * sigma ); 
	    double widthError = std::sqrt( std::pow( 2 * tau   * tauError  , 2 ) +
					   std::pow( 2 * sigma * sigmaError, 2 ) );

	    double   tau2 = fit2->GetParameter(1);
	    double sigma2 = fit2->GetParameter(2);

	    double   tauError2 = fit2->GetParError(1);
	    double sigmaError2 = fit2->GetParError(2);
	    
	    double width2      = std::sqrt( tau2 * tau2 + sigma2 * sigma2 ); 
	    double widthError2 = std::sqrt( std::pow( 2 * tau2   * tauError2  , 2 ) +
					    std::pow( 2 * sigma2 * sigmaError2, 2 ) );
	    
	    // Now, put results on histogram
	    std::pair< double, double > rmsAndError =
	      anaTool->GetRMS( hDphi, m_dPhiFittingMin, m_dPhiFittingMax, constants::PI );

	    double widthStat      = rmsAndError.first;
	    double widthStatError = rmsAndError.second;

	    std::cout << hDphi->GetName() << std::endl;
	    std::cout << "Tau   = " << tau   << " TauError   = " << tauError   << std::endl;
	    std::cout << "Sigma = " << sigma << " SigmaError = " << sigmaError << std::endl;
	    std::cout << "Width = " << width << " WidthError = " << widthError << std::endl;
	    std::cout << "Yield = " << yield << " YieldError = " << yieldError << std::endl;

	    hDphiYields->SetBinContent( axis3Bin, yield      );
	    hDphiYields->SetBinError  ( axis3Bin, yieldError );

	    // !!!!!!!!!!!!!!!!!!!!!!!!!!!
	    hDphiWidths->SetBinContent( axis3Bin, width      );
	    hDphiWidths->SetBinError  ( axis3Bin, widthError );	    

	    hDphiWidths2->SetBinContent( axis3Bin, width2      );
	    hDphiWidths2->SetBinError  ( axis3Bin, widthError2 );

	    hDphiWidthsStat->SetBinContent( axis3Bin, widthStat      );
	    hDphiWidthsStat->SetBinError  ( axis3Bin, widthStatError );

	  } // end loop over axis3

	  // check we are in correct ystar and pt bins
	  if( !m_dPP->CorrectPhaseSpace
	      ( std::vector<int>{ axis0Bin, axis1Bin, axis2Bin, 0 } ) )
	    { continue; }

	  // ----------- widths -----------
	  TCanvas cWidthsCmp( "cWidthsCmp", "cWidthsCmp", 800, 800 );

	  TPad pad1("pad1", "", 0.0, 0.35, 1.0, 1.0 );
	  pad1.SetBottomMargin(0.0);
	  pad1.Draw();
	  
	  TPad pad2("pad2", "", 0.0, 0.0, 1.0, 0.34 );
	  pad2.SetTopMargin(0.05);
	  pad2.SetBottomMargin(0.25);
	  pad2.Draw();

	  TLegend legWAll( 0.50, 0.1, 0.89, 0.28 );
	  styleTool->SetLegendStyle( &legWAll );

	  pad1.cd();

	  hDphiWidths->SetMinimum( m_dPhiWidthMin );
	  hDphiWidths->SetMaximum( m_dPhiWidthMax );

	  TH1* hDphiWidthsCmp = static_cast< TH1D* >
	    ( hDphiWidths->Clone( Form( "h_%s_cmp", hNameW.c_str() ) ) );
	  hDphiWidthsCmp->SetYTitle( "RMS (#pi - #Delta#phi)" );
	  hDphiWidthsCmp->SetMarkerSize( hDphiWidthsCmp->GetMarkerSize() * 1.5 );
	  vDphiWidths.push_back( hDphiWidthsCmp );
	  
	  styleTool->SetHStyle( hDphiWidthsCmp , 0 );
	  styleTool->SetHStyle( hDphiWidths2   , 6 );
	  styleTool->SetHStyle( hDphiWidthsStat, 5 );
	   
	  hDphiWidthsCmp->Draw("epsame X0");
	  hDphiWidths2  ->Draw("epsame X0");
	
	  legWAll.AddEntry( hDphiWidthsCmp,
			    Form( "Fit %2.1f<#Delta#phi<#pi", m_dPhiFittingMin  ) );
	  legWAll.AddEntry( hDphiWidths2  ,
			    Form( "Fit %2.1f<#Delta#phi<#pi", m_dPhiFittingMinB ) );
	  
	  legWAll.Draw();

	  DrawTopLeftLabels
	    ( m_dPP, axis0Low, axis0Up, axis1Low, axis1Up,
	      axis2Low, axis2Up, 0, 0 );

	  DrawAtlasRight();
	  
	  pad2.cd();

	  std::string hNameR = hNameW + "_" + m_sRatio;
	  TH1* hWR = static_cast< TH1D* >( hDphiWidths->Clone( hNameR.c_str() ) );
	  styleTool->SetHStyle( hWR, 0 );
	  vWR.push_back( hWR );
	  
	  hWR->SetMaximum( 1.5 );
	  hWR->SetMinimum( 0.5 );
	  hWR->SetYTitle( "Ratio" );

	  hWR->Divide( hDphiWidths2 );

	  hWR->Draw("ep X0" );
	  
	  line.Draw();
	  lineP25.Draw();
	  lineN25.Draw();
	  
	  SaveAsAll( cWidthsCmp, hNameW );

	  // ----------- yields -----------
	  TCanvas cYieldsAll( "cYieldsAll", "cYieldsAll",800,600);
	  cYieldsAll.SetLogy();
	  
	  hDphiYields->SetMinimum( m_dPhiYieldMin );
	  hDphiYields->SetMaximum( m_dPhiYieldMax );
	
	  hDphiYields->Draw("ep same X0");

	  DrawTopLeftLabels
	    ( m_dPP, axis0Low, axis0Up, axis1Low, axis1Up,
	      axis2Low, axis2Up, 0, 0 );

	  DrawAtlasRight();
	  
	  SaveAsAll( cYieldsAll, hNameY );
	} // end loop over axis2
	
	cWidths.cd();
	for( auto& h : vDphiWidthsTemp ){
	  h->SetMinimum( m_dPhiWidthMin );
	  h->SetMaximum( m_dPhiWidthMax );
	  h->SetTitle("");
	  h->Draw("epsame X0");
	  h->Write();
	}
	
	legW.Draw("same");

	DrawTopLeftLabels
	  ( m_dPP, axis0Low, axis0Up, axis1Low, axis1Up, 0, 0, 0, 0 );

	DrawAtlasRight();
	
        SaveAsAll( cWidths, Form("h_%s_%s_%s_%s",
				 dPhiName.c_str(), m_widthName.c_str(),
				 label.c_str(), hTagCW.c_str() ) );

	cYields.cd();
	cYields.SetLogy();
	for( auto& h : vDphiYieldsTemp ){
	  h->SetMinimum( m_dPhiYieldMin );
	  h->SetMaximum( m_dPhiYieldMax );
	  h->SetTitle("");
	  h->Draw("epsame X0");
	  h->Write();
	}
	
	legY.Draw("same");

	DrawTopLeftLabels
	  ( m_dPP, axis0Low, axis0Up, axis1Low, axis1Up, 0, 0, 0, 0 );

	DrawAtlasRight();
	
        SaveAsAll( cYields, Form("h_%s_%s_%s_%s",
				 dPhiName.c_str(), m_yieldName.c_str(),
				 label.c_str(), hTagCW.c_str() ) );

      } // end loop over axis1     
    } // end loop over axis0
    // SaveAsAll( cAll, Form("h_%s_%s", dPhiName. c_str(),label.c_str() ) );
  } // end loop over iG
  
  TCanvas c( "c", "c", 1200, 600 );
  c.Divide( 2, 1 );

  c.cd(1);
  h_dPhiChiS2->SetMaximum
    ( h_dPhiChiS->GetMaximum() > h_dPhiChiS2->GetMaximum() ?
      h_dPhiChiS->GetMaximum() * 1.1 :  h_dPhiChiS2->GetMaximum() * 1.1 );
  h_dPhiChiS2->Draw("hist C X0 same");
  h_dPhiChiS ->Draw("hist C X0 same");
  legChiSProb.Draw();

  
  c.cd(2);
  h_dPhiProb2->SetMaximum
    ( h_dPhiProb->GetMaximum() > h_dPhiProb2->GetMaximum() ?
      h_dPhiProb->GetMaximum() * 1.1 :  h_dPhiProb2->GetMaximum() * 1.1 );
  h_dPhiProb2->Draw("hist C X0 same");
  h_dPhiProb ->Draw("hist C X0 same");
  DrawAtlasRight();

  h_dPhiChiS ->Write();
  h_dPhiProb ->Write();
  h_dPhiChiS2->Write();
  h_dPhiProb2->Write();
  SaveAsAll( c, Form( "h_chi2_prob_%s", dPhiName.c_str() ) );

  delete h_dPhiChiS;
  delete h_dPhiProb;
  delete h_dPhiChiS2;
  delete h_dPhiProb2;

  std::cout << "------------ Bad Fits - " << v_listBadFits.size() << std::endl;
  for( uint i = 0; i < v_listBadFits.size(); i++ ){
    std::cout << v_listBadHist[i] << std::endl;
    std::cout << v_listBadFits[i] << std::endl;
    std::cout << "       ----   " << std::endl;
  }

  for( auto& r : vWR    ){ delete r; }
  for( auto& f : vFits  ){ delete f; }
  for( auto& h : vHdPhi ){ delete h; }
  for( auto& h : vSpect ){ delete h; }
  for( auto& w : vDphiWidths ){ delete w; }
  for( auto& y : vDphiYields ){ delete y; }
}

THnSparse* DiJetAnalysis::UnfoldDeltaPhi( TFile* fInData, TFile* fInMC,
					  const std::string& hnUnfoldedName,
					  TFile* fInMCSpect,
					  const std::string& spectName ){
  

  std::cout << fInData->GetName() << " " << fInMC->GetName() << std::endl;
  
  std::vector< TH1* > vHdPhi;
  std::vector< TH1* > vSpect;
  std::vector< TF1* > vDphiFits;
  
  std::vector< TH1* > vCfactors;
  std::vector< TH1* > vRatios;

  std::vector< TGraph* > vGdPhi;
  
  TAxis* axis0 = m_dPP->GetTAxis(0); int nAxis0Bins = axis0->GetNbins();
  TAxis* axis1 = m_dPP->GetTAxis(1); int nAxis1Bins = axis1->GetNbins();
  TAxis* axis2 = m_dPP->GetTAxis(2); int nAxis2Bins = axis2->GetNbins();
  TAxis* axis3 = m_dPP->GetTAxis(3); int nAxis3Bins = axis3->GetNbins();

  std::string unfoldingMCLabel = m_mcTypeLabel;
  
  std::string measuredName, truthName, unfoldedLabel, typeMeasured;
  GetDphiUnfoldingInfo( measuredName, truthName, unfoldedLabel, typeMeasured );

  // make a THnSparse to fill with unfolded results.
  THnSparse* hnUnfoldedCountsAll =
    new THnSparseD( Form("h_%s_%s", hnUnfoldedName.c_str(), m_allName.c_str() ) , "",
		    m_nDphiDim, &m_vNdPhiBins[0],&m_vDphiMin[0], &m_vDphiMax[0] );

  TAxis* axis0Def = m_dPP->GetDefaultTAxis( 0 );
  TAxis* axis1Def = m_dPP->GetDefaultTAxis( 1 );
  TAxis* axis2Def = m_dPP->GetDefaultTAxis( 2 );
  TAxis* axis3Def = m_dPP->GetDefaultTAxis( 3 );
  
  hnUnfoldedCountsAll->GetAxis(0)->
    Set( axis0Def->GetNbins(), axis0Def->GetXbins()->GetArray() );
  hnUnfoldedCountsAll->GetAxis(1)->
    Set( axis1Def->GetNbins(), axis1Def->GetXbins()->GetArray() );
  hnUnfoldedCountsAll->GetAxis(2)->
    Set( axis2Def->GetNbins(), axis2Def->GetXbins()->GetArray() );
  hnUnfoldedCountsAll->GetAxis(3)->
    Set( axis3Def->GetNbins(), axis3Def->GetXbins()->GetArray() );
  hnUnfoldedCountsAll->GetAxis(4)->Set( m_nVarDphiBins, &( m_varDphiBinning[0]  ) );
  
  hnUnfoldedCountsAll->GetAxis(0)->SetTitle( m_dPP->GetDefaultAxisLabel(0).c_str() );
  hnUnfoldedCountsAll->GetAxis(1)->SetTitle( m_dPP->GetDefaultAxisLabel(1).c_str() );
  hnUnfoldedCountsAll->GetAxis(2)->SetTitle( m_dPP->GetDefaultAxisLabel(2).c_str() );
  hnUnfoldedCountsAll->GetAxis(3)->SetTitle( m_dPP->GetDefaultAxisLabel(3).c_str() );
  hnUnfoldedCountsAll->GetAxis(4)->SetTitle( unfoldedLabel.c_str() );

  // this is more for debugging...
  // prepare c-factors histograms pt2 vs pt1 in bins of ystar1 ystar2
  std::vector< std::vector< TH2* > > vHcFactorsPtMat;
  
  for( int axis0Bin = 1; axis0Bin <= nAxis0Bins; axis0Bin++ ){ 
    // check we are in correct ystar and pt bins
    if( !m_dPP->CorrectPhaseSpace
	( std::vector<int>{ axis0Bin, 0, 0, 0 } ) )
      { continue; }

    double axis0Low, axis0Up;
    anaTool->GetBinRange( axis0, axis0Bin, axis0Bin, axis0Low, axis0Up );
      
    vHcFactorsPtMat.push_back( std::vector< TH2*>() );
      
    for( int axis3Bin = 1; axis3Bin <= nAxis3Bins; axis3Bin++ ){
      double axis3Low, axis3Up;
      anaTool->GetBinRange( axis3, axis3Bin, axis3Bin, axis3Low, axis3Up );
	
      TH2* hCfactorsPtMat =
	new TH2D( Form( "h_%s_ptMat_%s_%s_%s", m_cFactorName.c_str(), m_allName.c_str(),
			anaTool->GetName( axis0Low, axis0Up, m_dPP->GetAxisName(0) ).c_str(),
			anaTool->GetName( axis3Low, axis3Up, m_dPP->GetAxisName(3) ).c_str() ), 
		  ";p_{T1};p_{T2}",
		  m_nVarPtBins, &m_varPtBinning[0], 
		  m_nVarPtBins, &m_varPtBinning[0] );

      vHcFactorsPtMat.back().push_back( hCfactorsPtMat );
      styleTool->SetHStyle( hCfactorsPtMat, 0 );
    }
  }

  double xMin = m_dPhiZoomLow;
  double xMax = m_dPhiZoomHigh; 
	  
  TLine line( xMin, 1, xMax, 1 );
  line.SetLineWidth( 2 );

  TLine lineP25( xMin, 1.25, xMax, 1.25 );
  lineP25.SetLineStyle( 2  );
  lineP25.SetLineColor( 12 );
  lineP25.SetLineWidth( 2  );
	  
  TLine lineN25( xMin, 0.75, xMax, 0.75 );
  lineN25.SetLineStyle( 2  );
  lineN25.SetLineColor( 12 );
  lineN25.SetLineWidth( 2  );

  TLine lineP50( xMin, 1.50, xMax, 1.50 );
  lineP50.SetLineStyle( 2  );
  lineP50.SetLineColor( 12 );
  lineP50.SetLineWidth( 1  );
	  
  TLine lineN50( xMin, 0.50, xMax, 0.50 );
  lineN50.SetLineStyle( 2  );
  lineN50.SetLineColor( 12 );
  lineN50.SetLineWidth( 1  );

  // ---- loop over ystars ----
  for( int axis0Bin = 1; axis0Bin <= nAxis0Bins; axis0Bin++ ){
    // check we are in correct ystar and pt bins
    if( !m_dPP->CorrectPhaseSpace
	( std::vector<int>{ axis0Bin, 0, 0, 0 } ) )
      { continue; }

    double axis0Low, axis0Up;
    anaTool->GetBinRange
      ( axis0, axis0Bin, axis0Bin, axis0Low, axis0Up );

    TH1* hSpectCounts = static_cast< TH1D* >
      ( fInMCSpect->Get( Form( "h_%s_%s_%s_%s", spectName.c_str(), m_sCounts.c_str(),
			       m_allName.c_str(), anaTool->GetName
			       ( axis0Low, axis0Up, m_dPP->GetAxisName(0) ).c_str() ) ) );
    vSpect.push_back( hSpectCounts );

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
	  	  
	  TLegend leg( 0.68, 0.04, 0.92, 0.25 );
	  styleTool->SetLegendStyle( &leg );
	  

	  std::string hTag =
	    anaTool->GetName( axis0Low, axis0Up, m_dPP->GetAxisName(0) ) + "_" + 
	    anaTool->GetName( axis1Low, axis1Up, m_dPP->GetAxisName(1) ) + "_" + 
	    anaTool->GetName( axis2Low, axis2Up, m_dPP->GetAxisName(2) ) + "_" +  
	    anaTool->GetName( axis3Low, axis3Up, m_dPP->GetAxisName(3) ); 

	  
	  // Get measured and truth distributions
	  // both normalized and not. The counts are the
	  // un-normalized histograms.
	  // I.e. hMeasured and hTruth
	  // Since we already did normalization and fitting
	  // we just retrieve them.
	  TH1* hMeasuredCounts    = static_cast<TH1D*>
	    ( fInData->Get
	      ( Form( "h_%s_%s_%s_%s", measuredName.c_str(), m_sCounts.c_str(),
		      m_allName.c_str(), hTag.c_str())));
	  TH1* hMeasured          = static_cast<TH1D*>
	    ( fInData->Get( Form( "h_%s_%s_%s", measuredName.c_str(),
				  m_allName.c_str(), hTag.c_str())));
	  TH1* hTruth             = static_cast<TH1D*>
	    ( fInMC  ->Get( Form( "h_%s_%s_%s", truthName.c_str(),
				  m_allName.c_str(), hTag.c_str())));
	  
	  styleTool->SetHStyle( hMeasured, 1 );
	  styleTool->SetHStyle( hTruth   , 2 );

	  vHdPhi.push_back( hMeasuredCounts );
	  vHdPhi.push_back( hMeasured       );
	  vHdPhi.push_back( hTruth          );

	  // Get correction factors. IMPORTANT
	  TH1* hCfactors = static_cast<TH1D*>
	    ( fInMC->Get( Form( "h_%s_%s_%s", m_dPhiCfactorsName.c_str(),
				m_allName.c_str(), hTag.c_str())));
	  vCfactors.push_back( hCfactors );

	  // for compararison later.
	  // so its all in one file.
	  hCfactors->Write();
	  
	  // ----------- Unfold -----------
	  // Unfold using bin-by-bin and the resposne factors.
	  // Do this on counts, then normalize and subtract comb.
	  TH1* hUnfoldedCounts = BinByBinUnfolding
	    ( hMeasuredCounts,  hCfactors );
	  hUnfoldedCounts->SetName
	    ( Form( "h_%s_%s_%s_%s", m_dPhiUnfoldedName.c_str(),
		    m_sCounts.c_str(), m_allName.c_str(), hTag.c_str()));
	  styleTool->SetHStyle( hUnfoldedCounts, 0 );
	  vHdPhi.push_back( hUnfoldedCounts );

	  // fill unfolded THnSparse result
	  std::vector< int > x  = m_dPP->GetMappedBins
	    ( std::vector<int> { axis0Bin, axis1Bin, axis2Bin, axis3Bin } );
	  x.push_back(0); // x[4] = dPhiBin; 
	  for( int dPhiBin = 1; dPhiBin <= hUnfoldedCounts->GetNbinsX(); dPhiBin++ ){
	    x[4] = dPhiBin;
	    hnUnfoldedCountsAll->SetBinContent
	      ( &x[0], hUnfoldedCounts->GetBinContent( dPhiBin ) );
	    hnUnfoldedCountsAll->SetBinError
	      ( &x[0], hUnfoldedCounts->GetBinError  ( dPhiBin ) );
	  }

	  // Clone the counts, then normalize
	  TH1* hUnfolded = static_cast< TH1D* >
	    ( hUnfoldedCounts->Clone( Form( "h_%s_%s_%s", m_dPhiUnfoldedName.c_str(),
						m_allName.c_str(), hTag.c_str())));
	  hUnfolded->SetYTitle( hMeasured->GetYaxis()->GetTitle() );
	  vHdPhi.push_back( hUnfolded );

	  // Normalize after subtracting combinatoric
	  NormalizeDeltaPhi( hUnfolded, hSpectCounts, 0.5 * ( axis1Up + axis1Low ), true );

	  TF1* fitUnfolded = NULL;
	  if( m_fitDphiWC ){
	    fitUnfolded = anaTool->FitDphiWC( hUnfolded, m_dPhiFittingMin, m_dPhiFittingMax ); 
	  } else {
	    fitUnfolded = anaTool->FitDphi  ( hUnfolded, m_dPhiFittingMin, m_dPhiFittingMax );
	  }
	  
	  // -------- Unfold Done ---------
	  	  	  
	  hMeasured->GetXaxis()->SetRange( m_dPhiZoomLowBin, m_dPhiZoomHighBin );
	  hUnfolded->GetXaxis()->SetRange( m_dPhiZoomLowBin, m_dPhiZoomHighBin );
	  hTruth   ->GetXaxis()->SetRange( m_dPhiZoomLowBin, m_dPhiZoomHighBin );
	  
	  TCanvas c( "c", "c", 800, 800 );
	  TPad pad1("pad1", "", 0.0, 0.35, 1.0, 1.0 );
	  pad1.SetBottomMargin(0.0);
	  pad1.Draw();
	  
	  TPad pad2("pad2", "", 0.0, 0.0, 1.0, 0.34 );
	  pad2.SetTopMargin(0.05);
	  pad2.SetBottomMargin(0.25);
	  pad2.Draw();

	  pad1.cd();
	  pad1.SetLogy();

	  hMeasured->Draw("ep x0 same");	
	  hTruth   ->Draw("ep x0 same");
	  hUnfolded->Draw("ep x0 same");

	  leg.AddEntry( hTruth   , "Truth" );
	  leg.AddEntry( hMeasured, typeMeasured.c_str() );
	  leg.AddEntry( hUnfolded, "Unfolded" );

	  double maximum = -1;
	  
	  maximum = hUnfolded->GetMaximum() > hMeasured->GetMaximum() ?
	    hUnfolded->GetMaximum() : hMeasured->GetMaximum();
	  maximum = maximum > hTruth->GetMaximum() ?
	    maximum : hTruth->GetMaximum();

	  maximum = anaTool->GetLogMaximum( maximum );

	  hMeasured->SetMinimum( m_dPhiLogMin );
	  hMeasured->SetMaximum( m_dPhiLogMax );
	  
       	  TF1* fitMeasured = static_cast<TF1*>
	    ( fInData->Get( Form( "f_h_%s_%s_%s", measuredName.c_str(),
				  m_allName.c_str(), hTag.c_str())));
	  TF1* fitTruth    = static_cast<TF1*>
	    ( fInMC  ->Get( Form( "f_h_%s_%s_%s", truthName.c_str(),
				  m_allName.c_str(), hTag.c_str())));
	  
	  styleTool->SetHStyle( fitMeasured, 0 );
	  styleTool->SetHStyle( fitUnfolded, 0 );
	  styleTool->SetHStyle( fitTruth   , 0 );
	  
	  vDphiFits.push_back( fitMeasured );
	  vDphiFits.push_back( fitUnfolded );
	  vDphiFits.push_back( fitTruth    );

	  fitMeasured->SetLineColor( hMeasured->GetLineColor() );
	  fitUnfolded->SetLineColor( hUnfolded->GetLineColor() );
	  fitTruth   ->SetLineColor( hTruth   ->GetLineColor() );

	  fitMeasured->SetLineStyle( 7 );
	  fitUnfolded->SetLineStyle( 7 );
	  fitTruth   ->SetLineStyle( 7 );

	  fitMeasured->Draw("same");
	  fitUnfolded->Draw("same");
 	  fitTruth   ->Draw("same");

	  leg.Draw();

	  DrawTopLeftLabels
	    ( m_dPP, axis0Low, axis0Up, axis1Low, axis1Up,
	      axis2Low, axis2Up, axis3Low, axis3Up );
	    
	  DrawAtlasRight();

	  // put the RMS there
	  double   tau = fitUnfolded->GetParameter(1);
	  double sigma = fitUnfolded->GetParameter(2);

	  double   tauError = fitUnfolded->GetParError(1);
	  double sigmaError = fitUnfolded->GetParError(2);

	  double width      = std::sqrt( tau * tau + sigma * sigma ); 
	  double widthError = std::sqrt( std::pow( 2 * tau   * tauError  , 2 ) +
					 std::pow( 2 * sigma * sigmaError, 2 ) );

	  drawTool->DrawLeftLatex( 0.19, 0.48, Form( "RMS = %5.3f #pm %5.3f", width, widthError ) );
	    
	  pad2.cd();

	  TLegend legR( 0.2, 0.79, 0.6, 0.89 );
	  styleTool->SetLegendStyle( &legR, 0.9 );
	  legR.SetNColumns(3);
	  
	  // make sure to set range to rebinned.
	  // the c-factors are from rebinned distributions.
	  styleTool->SetHStyleRatio( hCfactors );
	  hCfactors->SetYTitle( "Ratio" );
	  hCfactors->GetXaxis()->SetRange
	    ( m_dPhiRebinnedZoomLowBin, m_dPhiRebinnedZoomHighBin );
	  hCfactors->SetTitleOffset( 2.3, "x" );
	  hCfactors->SetMaximum( 1.5 );
	  hCfactors->SetMinimum( 0.5 );
	  hCfactors->Draw("ep");
	  
	  TH1* hR = static_cast< TH1D* >
	    ( hUnfolded->
	      Clone( Form( "h_%s_%s_%s_%s_%s", m_sRatio.c_str(),
			   measuredName.c_str(), m_unfoldedName.c_str(),
			   m_allName.c_str(), hTag.c_str())));
	  styleTool->SetHStyleRatio( hR, 1 );
	  hR->Divide( hTruth );
	  // hR->Draw( "ep same" );

	  TH1* hRfit =  static_cast< TH1D* >
	    ( hUnfolded->
	      Clone( Form( "h_%s_fit_%s_%s_%s_%s", m_sRatio.c_str(),
			   measuredName.c_str(), m_unfoldedName.c_str(),
			   m_allName.c_str(), hTag.c_str())));
	  styleTool->SetHStyleRatio( hRfit, 2 );
	  hRfit->Reset();

	  int xBinMin = hRfit->FindBin( m_dPhiFittingMin );

	  for( int xBin = xBinMin; xBin <= hRfit->GetNbinsX(); xBin++ ){

	    double binCenter = hRfit->GetBinCenter( xBin );

	    if( !hUnfolded->GetBinContent( xBin ) ){ continue; };
	    
	    double ratio = fitUnfolded->Eval( binCenter ) /
	      hUnfolded->GetBinContent( xBin );
	    
	    double ratioError = hUnfolded->GetBinError( xBin ) /
	      fitUnfolded->Eval( binCenter );
	      
	    hRfit->SetBinContent( xBin, ratio      );
	    hRfit->SetBinError  ( xBin, ratioError );
	  }

	  // hRfit->Draw( "ep same" );
	  
	  legR.AddEntry( hCfactors, "CF" );
	  // legR.AddEntry( hR    , "UF/Truth" );
	  // legR.AddEntry( hRfit    , "Fit/UF" );
	  
	  legR.Draw();

	  styleTool->HideAxis( hMeasured, "x" );
	  styleTool->HideAxis( hUnfolded, "x" );
	  styleTool->HideAxis( hTruth   , "x" );
	  
	  line.Draw();
	  lineP25.Draw();
	  lineN25.Draw();
	  // lineP50.Draw();
	  // lineN50.Draw();
	  
	  // save this for debugging
	  hUnfoldedCounts->Write();
	  
	  // MUT is Measuerd Unfolded Truth
	  SaveAsAll( c, Form( "h_%s_%s_MUT_%s",
			      m_dPhiUnfoldedName.c_str(),
			      m_allName.c_str(), hTag.c_str()));

	  // fill the pt2 vs pt1 response histos.
	  TH2*   hCfactorsPtMat   = vHcFactorsPtMat[ axis0Bin - 1 ][ axis3Bin - 1 ];	  
	  int    dPhiPtCfactorBin = hCfactors->GetNbinsX();
	  double binCfactor       = hCfactors->GetBinContent( dPhiPtCfactorBin );

	  hCfactorsPtMat->SetBinContent( axis1Bin, axis2Bin, binCfactor );
	} // end loop over axis3
      } // end loop over axis2
    } // end loop over axis1     
  } // end loop over axis0

  // now draw the pt2 pt1 matrixes with c factors
  for( int axis0Bin = 1; axis0Bin <= nAxis0Bins; axis0Bin++ ){ 
    // check we are in correct ystar and pt bins
    if( !m_dPP->CorrectPhaseSpace
	( std::vector<int>{ axis0Bin, 0, 0, 0 } ) )
      { continue; }

    double axis0Low, axis0Up;
    anaTool->GetBinRange( axis0, axis0Bin, axis0Bin, axis0Low, axis0Up );

    for( int axis3Bin = 1; axis3Bin <= nAxis3Bins; axis3Bin++ ){
      double axis3Low, axis3Up;
      anaTool->GetBinRange( axis3, axis3Bin, axis3Bin, axis3Low, axis3Up );

      TH2* hCfactorsPtMat = vHcFactorsPtMat[ axis0Bin - 1 ][ axis3Bin - 1 ];	  

      TCanvas c1( hCfactorsPtMat->GetName(), hCfactorsPtMat->GetName(), 800, 600 );
    
      gStyle->SetPaintTextFormat(".2f");
      hCfactorsPtMat->Draw("colz text");
      hCfactorsPtMat->SetTitle("");

      drawTool->DrawLeftLatex
	( 0.13, 0.86, anaTool->GetLabel
	  ( axis0Low, axis0Up, m_dPP->GetDefaultAxisLabel(0) ) );
      drawTool->DrawLeftLatex
	( 0.13, 0.79, anaTool->GetLabel
	  ( axis3Low, axis3Up, m_dPP->GetDefaultAxisLabel(1) ) );  

      DrawAtlasRight();

      SaveAsAll( c1, hCfactorsPtMat->GetName() );
      hCfactorsPtMat->Write();
    }
  }
  
  for( auto vh : vHcFactorsPtMat )
    { for( auto h : vh ){ delete h; } }

  for( auto f : vDphiFits ){ delete f;  }
  for( auto h : vHdPhi    ){ delete h ; }
  for( auto h : vSpect    ){ delete h ; }

  for( auto r : vCfactors ){ delete r;  }

  hnUnfoldedCountsAll->Write();
  
  return hnUnfoldedCountsAll;
}

void DiJetAnalysis::MakeSpectTogether( TFile* fOut ){

  // here it is only done for ystar!
  // later, tehre should be some switch
  // if this is required for other binning.
  // for example, eta.
  
  std::vector< TH1* > vH;
  std::vector< TH1* > vR;
  
  std::string name_a , name_b  ;
  std::string label_a, label_b ;
  std::string fName_a, fName_b;

  GetInfoTogether( name_a, name_b, label_a, label_b, fName_a, fName_b, 0 );

  std::string ratio = Form("%s/%s", label_a.c_str(), label_b.c_str() );

  double xMin = m_varPtBinning.front();
  double xMax = m_varPtBinning.back();
  
  TLine line( xMin, 1, xMax, 1 );
  line.SetLineWidth( 2 );

  TLine lineP15( xMin, 1.15, xMax, 1.15 );
  lineP15.SetLineStyle( 2  );
  lineP15.SetLineColor( 12 );
  lineP15.SetLineWidth( 2  );
	  
  TLine lineN15( xMin, 0.85, xMax, 0.85 );
  lineN15.SetLineStyle( 2  );
  lineN15.SetLineColor( 12 );
  lineN15.SetLineWidth( 2  );
  
  std::string axisLabel, axisLabelTex;
  GetSpectraLabels( axisLabel, axisLabelTex, name_a );
  
  TFile* fIn_a = TFile::Open( fName_a.c_str() );
  TFile* fIn_b = TFile::Open( fName_b.c_str() );

  fOut->cd();
  
  for( uint xBin = 1; xBin <= m_nVarYstarBins; xBin++ ){

    double xLow  = m_varYstarBinning[ xBin - 1 ];
    double xUp   = m_varYstarBinning[ xBin ];
   
    std::string hTag = anaTool->GetName( xLow, xUp, axisLabel);

    std::string hName_a = "h_" + name_a + "_" + m_allName + "_" + hTag;
    std::string hName_b = "h_" + name_b + "_" + m_allName + "_" + hTag;
    
    TH1* hSpect_a = static_cast< TH1D* >( fIn_a->Get( hName_a.c_str() ) );
    TH1* hSpect_b = static_cast< TH1D* >( fIn_b->Get( hName_a.c_str() ) ); 

    styleTool->SetHStyle( hSpect_a, 0 );
    styleTool->SetHStyle( hSpect_b, 1 );

    vH.push_back( hSpect_a );
    vH.push_back( hSpect_b );
    
    TLegend leg( 0.65, 0.6, 0.85, 0.8 );
    styleTool->SetLegendStyle( &leg );

    leg.AddEntry( hSpect_a, label_a.c_str() );
    leg.AddEntry( hSpect_b, label_b.c_str() );
    
    double scalingFactor =
      hSpect_a->Integral( 2, m_nVarPtBinsUfOf - 1 ) /
      hSpect_b->Integral( 2, m_nVarPtBinsUfOf - 1 );

    hSpect_b->Scale( scalingFactor );
    
    TCanvas c( "c", "c", 800, 800 );
    TPad pad1("pad1", "", 0.0, 0.35, 1.0, 1.0 );
    pad1.SetBottomMargin(0);
    pad1.Draw();
    TPad pad2("pad2", "", 0.0, 0.0, 1.0, 0.34 );
    pad2.SetTopMargin(0.05);
    pad2.SetBottomMargin(0.25);
    pad2.Draw();

    pad2.cd();
    
    TH1* hR = static_cast< TH1D* >
      ( hSpect_a->Clone( Form( "h_%s_%s_%s_%s", m_ystarSpectName.c_str(),
			      m_sRatio.c_str(), m_allName.c_str(), hTag.c_str())));
    styleTool->SetHStyleRatio( hR );
    vR.push_back( hR );

    hR->Divide( hSpect_b );
    hR->SetYTitle( ratio.c_str() );
    hR->SetTitleOffset( 2.3, "x" );
    hR->Draw( "ep same" );

    line.Draw();
    lineP15.Draw(); lineN15.Draw();
    
    pad1.cd();
    pad1.SetLogy();

    styleTool->HideAxis( hSpect_a, "x" );
    styleTool->HideAxis( hSpect_b, "x" );
    
    hSpect_a->Draw( "ep same" );
    hSpect_b->Draw( "ep same" );

    leg.Draw( "same" );

    DrawAtlasRightBoth();

    drawTool->DrawRightLatex
      ( 0.45, 0.15, anaTool->GetLabel( xLow, xUp, axisLabelTex ) );
        
    SaveAsAll( c, Form("h_%s_%s", m_ystarSpectName.c_str(), hTag.c_str() ), true );
  }

  for( auto& h : vH ){ delete h; }
  for( auto& r : vR ){ delete r; }
}

void DiJetAnalysis::MakeDphiTogether( TFile* fOut ){

  std::vector< TF1* > vF;
  std::vector< TH1* > vH;
  std::vector< TH1* > vHw;
  std::vector< TH1* > vHy;
  std::vector< TH1* > vR;
  
  std::string name_a , name_b  ;
  std::string label_a, label_b ;
  std::string fName_a, fName_b;

  GetInfoTogether( name_a, name_b, label_a, label_b, fName_a, fName_b, 1 );

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
  lineP25.SetLineStyle( 2  );
  lineP25.SetLineColor( 12 );
  lineP25.SetLineWidth( 2  );
	  
  TLine lineN25( xMin, 0.75, xMax, 0.75 );
  lineN25.SetLineStyle( 2  );
  lineN25.SetLineColor( 12 );
  lineN25.SetLineWidth( 2  );

  TLine lineP50( xMin, 1.50, xMax, 1.50 );
  lineP50.SetLineStyle( 2  );
  lineP50.SetLineColor( 12 );
  lineP50.SetLineWidth( 1  );
	  
  TLine lineN50( xMin, 0.50, xMax, 0.50 );
  lineN50.SetLineStyle( 2  );
  lineN50.SetLineColor( 12 );
  lineN50.SetLineWidth( 1  );

  TLine lineP1S( xMin, 1.34, xMax, 1.34 );
  lineP1S.SetLineStyle( 2  );
  lineP1S.SetLineColor( 12 );
  lineP1S.SetLineWidth( 1  );
	  
  TLine lineN1S( xMin, 0.66, xMax, 0.66 );
  lineN1S.SetLineStyle( 2  );
  lineN1S.SetLineColor( 12 );
  lineN1S.SetLineWidth( 1  );
  
  TFile* fIn_a = TFile::Open( fName_a.c_str() );
  TFile* fIn_b = TFile::Open( fName_b.c_str() );

  xMin = m_dPhiZoomLow;
  xMax = m_dPhiZoomHigh;
  
  TLine dPhiLine( xMin, 1, xMax, 1 );
  dPhiLine.SetLineWidth( 2 );
	  
  TLine dPhiLineP25( xMin, 1.25, xMax, 1.25 );
  dPhiLineP25.SetLineStyle( 3  );
  dPhiLineP25.SetLineColor( 12 );
  dPhiLineP25.SetLineWidth( 2  );
	  
  TLine dPhiLineN25( xMin, 0.75, xMax, 0.75 );
  dPhiLineN25.SetLineStyle( 3 );
  dPhiLineN25.SetLineColor( 12 );
  dPhiLineN25.SetLineWidth( 2  );
	  
  TLine dPhiLineP50( xMin, 1.50, xMax, 1.50 );
  dPhiLineP50.SetLineStyle( 3  );
  dPhiLineP50.SetLineColor( 12 );
  dPhiLineP50.SetLineWidth( 1  );
	  
  TLine dPhiLineN50( xMin, 0.50, xMax, 0.50 );
  dPhiLineN50.SetLineStyle( 3 );
  dPhiLineN50.SetLineColor( 12 );
  dPhiLineN50.SetLineWidth( 1  );

  fOut->cd();

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
      TCanvas cW( "cW", "cW", 800, 800 );
      TPad pad1W("pad1W", "", 0.0, 0.35, 1.0, 1.0 );
      pad1W.SetBottomMargin(0);
      pad1W.Draw();
      TPad pad2W("pad2W", "", 0.0, 0.0, 1.0, 0.34 );
      pad2W.SetTopMargin(0.05);
      pad2W.SetBottomMargin(0.25);
      pad2W.Draw();

      // Now Draw everything.
      TCanvas cY( "cY", "cY", 800, 800 );
      TPad pad1Y("pad1Y", "", 0.0, 0.35, 1.0, 1.0 );
      pad1Y.SetBottomMargin(0);
      pad1Y.Draw();
      TPad pad2Y("pad2Y", "", 0.0, 0.0, 1.0, 0.34 );
      pad2Y.SetTopMargin(0.05);
      pad2Y.SetBottomMargin(0.25);
      pad2Y.Draw();
      
      // tag for widths canvas
      std::string hTagC =
	Form ("%s_%s",
	      anaTool->GetName( axis0Low, axis0Up, m_dPP->GetAxisName(0) ).c_str(),
	      anaTool->GetName( axis1Low, axis1Up, m_dPP->GetAxisName(1) ).c_str() );

      // this is bs. works though. vary the last factor 
      double deltaYleg = ( axis1Bin - 1 ) * 0.075;
      
      TLegend legW( 0.50, 0.03, 0.85, 0.13 + deltaYleg );
      styleTool->SetLegendStyle( &legW, 0.75 );

      TLegend legY( 0.31, 0.03, 0.72, 0.13 + deltaYleg  );
      styleTool->SetLegendStyle( &legY, 0.75 );
      
      int style = 0;
      
      for( int axis2Bin = 1; axis2Bin <= nAxis2Bins; axis2Bin++ ){
	double axis2Low , axis2Up;
	anaTool->GetBinRange
	  ( axis2, axis2Bin, axis2Bin, axis2Low, axis2Up );
	
	// get widths histos
	std::string hTag =
	  Form("%s_%s_%s",
	       anaTool->GetName( axis0Low, axis0Up, m_dPP->GetAxisName(0) ).c_str(),
	       anaTool->GetName( axis1Low, axis1Up, m_dPP->GetAxisName(1) ).c_str(),
	       anaTool->GetName( axis2Low, axis2Up, m_dPP->GetAxisName(2) ).c_str() );
	
	std::string hNameW_a =
	  "h_" + name_a + "_" + m_widthName + "_" + m_allName + "_" + hTag;
	std::string hNameW_b =
	  "h_" + name_b + "_" + m_widthName + "_" + m_allName + "_" + hTag;

	TH1* hW_a = static_cast<TH1D*>( fIn_a->Get( hNameW_a.c_str() ) );
	TH1* hW_b = static_cast<TH1D*>( fIn_b->Get( hNameW_b.c_str() ) );

	styleTool->SetHStyle( hW_a, style );
	styleTool->SetHStyle( hW_b, style + 5 );
	hW_a->SetMarkerSize( hW_a->GetMarkerSize() * 1.5 );
	hW_b->SetMarkerSize( hW_b->GetMarkerSize() * 1.5 );
	vHw.push_back( hW_a ); vHw.push_back( hW_b  );

	std::string hNameY_a =
	  "h_" + name_a + "_" + m_yieldName + "_" + m_allName + "_" + hTag;
	std::string hNameY_b =
	  "h_" + name_b + "_" + m_yieldName + "_" + m_allName + "_" + hTag;

	TH1* hY_a = static_cast<TH1D*>( fIn_a->Get( hNameY_a.c_str() ) );
	TH1* hY_b = static_cast<TH1D*>( fIn_b->Get( hNameY_b.c_str() ) );

	styleTool->SetHStyle( hY_a, style );
	styleTool->SetHStyle( hY_b, style + 5 );
	hY_a->SetMarkerSize( hY_a->GetMarkerSize() * 1.5 );
	hY_b->SetMarkerSize( hY_b->GetMarkerSize() * 1.5 );
	vHy.push_back( hY_a ); vHy.push_back( hY_b  );

	// switch back to cW canvas
	// because a new one was creted in axis3 loop
	cW.cd();
	pad2W.cd();

	// Make the ratio histogram.
	TH1* hW_R = static_cast< TH1D* >
	  ( hW_a->Clone( Form( "h_%s_%s_%s_%s", m_dPhiName.c_str(), m_sRatio.c_str(),
			       m_allName.c_str(), hTag.c_str())));
	hW_R->Divide( hW_b );
	hW_R->SetYTitle( ratio.c_str() );
	hW_R->SetTitleOffset( 2.3, "x" );
	styleTool->SetHStyleRatio( hW_R, style );
	hW_R->SetMarkerSize( hW_R->GetMarkerSize() * 1.5 );
	
	vR.push_back( hW_R );
	hW_R->Draw("ep x0 same");

	pad1W.cd();
	
	styleTool->HideAxis( hW_a, "x" );
	styleTool->HideAxis( hW_b, "x" );

	if( hW_a->GetMean() ){
	  legW.AddEntry ( hW_a,Form("%s %s", label_a.c_str(), anaTool->GetLabel
				    ( axis2Low, axis2Up, m_dPP->GetAxisLabel(2) ).c_str()));
	}
	if( hW_b->GetMean() ){
	  legW.AddEntry( hW_b, Form("%s %s", label_b.c_str(), anaTool->GetLabel
				    ( axis2Low, axis2Up, m_dPP->GetAxisLabel(2) ).c_str()));	
	}
	
	hW_a->Draw("ep same X0");
	hW_b->Draw("ep same X0");

	// switch back to cY canvas
	// because a new one was creted in axis3 loop
	cY.cd();
	pad2Y.cd();

	// Make the ratio histogram.
	TH1* hY_R = static_cast< TH1D* >
	  ( hY_a->Clone( Form( "h_%s_%s_%s_%s", m_dPhiName.c_str(), m_sRatio.c_str(),
			       m_allName.c_str(), hTag.c_str())));
	hY_R->Divide( hY_b );
	hY_R->SetYTitle( ratio.c_str() );
	hY_R->SetTitleOffset( 2.3, "x" );
	styleTool->SetHStyleRatio( hY_R, style );
	hY_R->SetMarkerSize( hY_R->GetMarkerSize() * 1.5 );
	
	vR.push_back( hY_R );
	hY_R->Draw("ep x0 same");

	style++;
	
	pad1Y.cd();
	pad1Y.SetLogy();
	
	styleTool->HideAxis( hY_a, "x" );
	styleTool->HideAxis( hY_b, "x" );

	if( hY_a->GetMean() ){
	  legY.AddEntry ( hY_a,Form("%s %s", label_a.c_str(), anaTool->GetLabel
				    ( axis2Low, axis2Up, m_dPP->GetAxisLabel(2) ).c_str()));
	}
	if( hY_b->GetMean() ){
	  legY.AddEntry( hY_b, Form("%s %s", label_b.c_str(), anaTool->GetLabel
				    ( axis2Low, axis2Up, m_dPP->GetAxisLabel(2) ).c_str()));	
	}
	
	hY_a->Draw("ep same X0");
	hY_b->Draw("ep same X0");

	for( int axis3Bin = 1; axis3Bin <= nAxis3Bins; axis3Bin++ ){
	  // check we are in correct ystar and pt bins
	  if( !m_dPP->CorrectPhaseSpace
	      ( std::vector<int>{ axis0Bin, axis1Bin, axis2Bin, axis3Bin } ) )
	    { continue; }

	  double axis3Low , axis3Up;
	  anaTool->GetBinRange
	    ( axis3, axis3Bin, axis3Bin, axis3Low, axis3Up );

	  // Now Draw everything.
	  TCanvas c( "c", "c", 800, 800 );
	  TPad pad1("pad1", "", 0.0, 0.35, 1.0, 1.0 );
	  pad1.SetBottomMargin(0);
	  pad1.Draw();
	  TPad pad2("pad2", "", 0.0, 0.0, 1.0, 0.34 );
	  pad2.SetTopMargin(0.05);
	  pad2.SetBottomMargin(0.25);
	  pad2.Draw();

	  std::string hTagDphi =
	    Form("%s_%s_%s_%s",
		 anaTool->GetName( axis0Low, axis0Up, m_dPP->GetAxisName(0) ).c_str(),
		 anaTool->GetName( axis1Low, axis1Up, m_dPP->GetAxisName(1) ).c_str(),
		 anaTool->GetName( axis2Low, axis2Up, m_dPP->GetAxisName(2) ).c_str(),
		 anaTool->GetName( axis3Low, axis3Up, m_dPP->GetAxisName(3) ).c_str() );

	  std::string hName_a = "h_" + name_a + "_" + m_allName + "_" + hTagDphi;
	  std::string hName_b = "h_" + name_b + "_" + m_allName + "_" + hTagDphi;

	  std::string gName_a = "g_" + name_a + "_" + m_allName + "_" + hTagDphi;
	  std::string gName_b = "g_" + name_b + "_" + m_allName + "_" + hTagDphi;
	  
	  TH1* h_a = static_cast<TH1D*>( fIn_a->Get( hName_a.c_str() ) );
	  TH1* h_b = static_cast<TH1D*>( fIn_b->Get( hName_b.c_str() ) );
	  styleTool->SetHStyle( h_a, 0 );
	  styleTool->SetHStyle( h_b, 1 );
	  vH.push_back( h_a ); vH.push_back( h_b );
	  
	  TH1* h_R = static_cast< TH1D* >
	    ( h_a->Clone( Form( "h_%s_%s_%s_%s", m_dPhiName.c_str(), m_sRatio.c_str(),
				m_allName.c_str(), hTagDphi.c_str())));
	  styleTool->SetHStyleRatio( h_R );
	  h_R->Divide( h_b );
	  h_R->SetYTitle( ratio.c_str() );
	  h_R->SetTitleOffset( 2.3, "x" );
	  vR.push_back( h_R );

	  pad1.cd();
	  pad1.SetLogy();
	  
	  styleTool->HideAxis( h_a, "x" );
	  styleTool->HideAxis( h_b, "x" );

	  h_a->GetXaxis()->SetRange( m_dPhiZoomLowBin, m_dPhiZoomHighBin );
	  
	  h_a->Draw("ep x0 same");
	  h_b->Draw("ep x0 same");
	  
	  TF1* f_a = static_cast<TF1*>( fIn_a->Get( Form("f_%s", hName_a.c_str())));
	  TF1* f_b = static_cast<TF1*>( fIn_b->Get( Form("f_%s", hName_b.c_str())));
	  styleTool->SetHStyle( f_a, 0 );
	  styleTool->SetHStyle( f_b, 1 );
	  f_a->SetLineColor( h_a->GetLineColor() );
	  f_b->SetLineColor( h_b->GetLineColor() );
	  vF.push_back( f_a ); vF.push_back( f_b );

	  TLegend leg( 0.7, 0.03, 0.9, 0.18 );
	  styleTool->SetLegendStyle( &leg , 0.9 );

	  h_a->SetMinimum( m_dPhiLogMin );
	  h_b->SetMinimum( m_dPhiLogMin );

	  h_a->SetMaximum( m_dPhiLogMax );
	  h_b->SetMaximum( m_dPhiLogMax );
	  
	  leg.AddEntry( h_a, label_a.c_str() );
	  leg.AddEntry( h_b, label_b.c_str() );
	 
	  f_a->Draw("same");
	  f_b->Draw("same");
	  
	  leg.Draw("same");

	  double max =  h_a->GetMaximum() > h_b->GetMaximum() ?
	    h_a->GetMaximum() : h_b->GetMaximum();

	  max = anaTool->GetLogMaximum( max );
	  
	  h_a->SetMinimum( m_dPhiLogMin );
	  h_a->SetMaximum( m_dPhiLogMax );
	  
	  DrawTopLeftLabels
	    ( m_dPP, axis0Low, axis0Up, axis1Low, axis1Up,
	      axis2Low, axis2Up, axis3Low, axis3Up );

	  drawTool->DrawAtlasInternal();
	  drawTool->DrawLeftLatex( 0.75, 0.2, m_isData ? "Data" : "MC" );

	  double width_a      = hW_a->GetBinContent( axis3Bin );
	  double widthError_a = hW_a->GetBinError  ( axis3Bin );

	  double width_b      = hW_b->GetBinContent( axis3Bin );
	  double widthError_b = hW_b->GetBinError  ( axis3Bin );

	  drawTool->DrawLeftLatex
	    ( 0.19, 0.48, Form( "Width_{%s}=%4.2f#pm%4.2f",
				label_a.c_str(), width_a, widthError_a ) );
	  drawTool->DrawLeftLatex
	    ( 0.19, 0.38, Form( "Width_{%s}=%4.2f#pm%4.2f",
				label_b.c_str(), width_b, widthError_b ) );

	  double yield_a      = hY_a->GetBinContent( axis3Bin );
	  double yieldError_a = hY_a->GetBinError  ( axis3Bin );

	  double yield_b      = hY_b->GetBinContent( axis3Bin );
	  double yieldError_b = hY_b->GetBinError  ( axis3Bin );
	  
	  drawTool->DrawLeftLatex
	    ( 0.53, 0.85, Form( "Yield_{%s}=%0.3f#pm%0.3f",
				label_a.c_str(),yield_a, yieldError_a ) );
	  drawTool->DrawLeftLatex
	    ( 0.53, 0.77, Form( "Yield_{%s}=%0.3f#pm%0.3f",
				label_b.c_str(), yield_b, yieldError_b ) );
	  
	  pad2.cd();
	  h_R->GetXaxis()->SetRange( m_dPhiZoomLowBin, m_dPhiZoomHighBin );
	  h_R->Draw("ep x0 same");

	  dPhiLine.Draw();
	  dPhiLineP25.Draw(); dPhiLineN25.Draw();
	  // dPhiLineP50.Draw(); dPhiLineN50.Draw();

	  h_R->Write();
	  SaveAsAll( c, Form("h_%s_%s", m_dPhiName.c_str(), hTagDphi.c_str() ), true );
	} // end loop over axis3
      } // end loop over axis2

      // back to cW canvas
      cW.cd();
      pad1W.cd();
      
      legW.Draw("same");

      DrawTopLeftLabels
	( m_dPP, axis0Low, axis0Up, axis1Low, axis1Up,
	  0, 0, 0, 0 );

      DrawAtlasRightBoth();

      pad2W.cd();
      line.Draw();
      lineP25.Draw();
      lineN25.Draw();
      lineP50.Draw();
      lineN50.Draw();

      SaveAsPdfPng( cW, Form("h_%s_%s_%s", m_dPhiName.c_str(),
			     m_widthName.c_str(), hTagC.c_str() ), true );
      SaveAsROOT( cW, Form("h_%s_%s_%s", m_dPhiName.c_str(),
			   m_widthName.c_str(), hTagC.c_str() ) );

      // back to cW canvas
      cY.cd();
      pad1Y.cd();
      
      legY.Draw("same");

      DrawTopLeftLabels
	( m_dPP, axis0Low, axis0Up, axis1Low, axis1Up,
	  0, 0, 0, 0 );

      DrawAtlasRightBoth();

      pad2Y.cd();
      line.Draw();
      lineP25.Draw();
      lineN25.Draw();
      lineP50.Draw();
      lineN50.Draw();

      SaveAsPdfPng( cY, Form("h_%s_%s_%s", m_dPhiName.c_str(),
			     m_yieldName.c_str(), hTagC.c_str() ), true );
      SaveAsROOT( cY, Form("h_%s_%s_%s", m_dPhiName.c_str(),
			   m_yieldName.c_str(), hTagC.c_str() ) );

      
    } // end loop over ystar2
  } // end loop over ystar2

  TH1* h_chi2_a = static_cast< TH1D* >
    ( fIn_a->Get( Form( "h_dPhiChiS_%s", name_a.c_str() ) ) );
  styleTool->SetHStyle( h_chi2_a, 0 );
  
  TH1* h_chi2_b = static_cast< TH1D* >
    ( fIn_b->Get( Form( "h_dPhiChiS_%s", name_b.c_str() ) ) );
  styleTool->SetHStyle( h_chi2_b, 1 );

  TH1* h_prob_a = static_cast< TH1D* >
    ( fIn_a->Get( Form( "h_dPhiProb_%s", name_a.c_str() ) ) );
  styleTool->SetHStyle( h_prob_a, 0 );
  
  TH1* h_prob_b = static_cast< TH1D* >
    ( fIn_b->Get( Form( "h_dPhiProb_%s", name_b.c_str() ) ) );
  styleTool->SetHStyle( h_prob_b, 1 );

  TCanvas c( "c", "c", 1200, 600 );
  c.Divide( 2, 1 );

  TLegend leg( 0.6, 0.6, 0.7, 0.7 );
  styleTool->SetLegendStyle( &leg );
  
  c.cd(1);
  h_chi2_a->SetMaximum( h_chi2_a->GetMaximum() > h_chi2_b->GetMaximum() ?
			h_chi2_a->GetMaximum() * 1.1 : h_chi2_b->GetMaximum() * 1.1 );
  h_chi2_a->Draw("hist p same");
  h_chi2_b->Draw("hist p same");
  DrawAtlasRightBoth();

  c.cd(2);
  h_prob_a->SetMaximum( h_prob_a->GetMaximum() > h_prob_b->GetMaximum() ?
			h_prob_a->GetMaximum() * 1.1 : h_prob_b->GetMaximum() * 1.1 );
  h_prob_a->Draw("hist p same");
  h_prob_b->Draw("hist p same");

  leg.AddEntry( h_prob_a, label_a.c_str() );
  leg.AddEntry( h_prob_b, label_b.c_str() );

  leg.Draw();
  
  DrawAtlasRightBoth();

  SaveAsPdfPng( c, "h_chi2_prob", true );
  SaveAsROOT  ( c, "h_chi2_prob");

  /*
  delete h_chi2_a;
  delete h_chi2_b;
  delete h_prob_a;
  delete h_prob_b;
  */
   
  for( auto & f  : vF  ){ delete f ; }
  for( auto & r  : vR  ){ delete r ; }
  for( auto & h  : vH  ){ delete h ; }
  for( auto & hW : vHw ){ delete hW; }
  for( auto & hY : vHy ){ delete hY; }
}

void DiJetAnalysis::MakeFinalPlotsTogether( TFile* fOut, const std::string& name ){}

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

  double dy = 0.08 * scale;
  double ystart = 0.86 + ( 1 - scale ) * 0.1;
  
  drawTool->DrawLeftLatex
    ( 0.18, ystart, CT::AnalysisTools::GetLabel
      ( axis0Low, axis0Up, dPP->GetAxisLabel(0) ), scale );
  // for now, modify later to be more dynamic
  // to the variables present 
  if( !( axis1Low || axis1Up ) ){ return ; }
  drawTool->DrawLeftLatex
    ( 0.18, ystart - dy, CT::AnalysisTools::GetLabel
      ( axis1Low, axis1Up, dPP->GetAxisLabel(1) ), scale );  
  if( !( axis2Low || axis2Up ) ){ return ; }
  drawTool->DrawLeftLatex
    ( 0.18, ystart - 2 * dy, CT::AnalysisTools::GetLabel
      ( axis2Low, axis2Up, dPP->GetAxisLabel(2) ), scale );
  if( !( axis3Low || axis3Up ) ){ return ; }
  drawTool->DrawLeftLatex
    ( 0.18, ystart - 3 * dy, CT::AnalysisTools::GetLabel
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
  
  // only save png pdf if the doing nominal sample. (0)
  // otherwise it is too many pdf/pngs.
  if( !m_uncertComp )
    { SaveAsPdfPng( c, name, together ); }
  SaveAsROOT  ( c, name );
}
