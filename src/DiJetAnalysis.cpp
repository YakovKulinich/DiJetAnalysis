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

  m_sDphi       = "#Delta#phi";
  m_sDphiTitle  = "#it{C}_{12}";
  m_sWidthTitle = "#it{W}_{12}";
  m_sYieldTitle = "#it{I}_{12}";
  m_sWidthRatioTitle = "#rho_{#it{W}}^{pPb}";
  m_sYieldRatioTitle = "#rho_{#it{I}}^{pPb}";

  m_unweightedFileSuffix  = "UW";
  m_unfoldingFileSuffix   = "UF";
  m_systematicsFileSuffix = "SYS";
  
  boost::assign::push_back( m_vMCtypeNames )
    ( "pythia8" )( "herwig" )( "pythia8powheg" );

  boost::assign::push_back( m_vMCtypeLabels  )
    ( "Pythia8" )( "Herwig" )( "Pythia8+Powheg" );

  //==================== Cuts ====================
  m_nMinEntriesFit = 20;

  m_deltaPtCut = 0.0; // GeV

  m_jetDeltaR = 0.4;

  // increase the size of the jet or the region
  // excluded by hec for systematics.
  if( m_uncertComp == 24 ){
    m_jetDeltaR = 0.5;
  }
  
  m_hecEtaMinA = 1.5 - m_jetDeltaR;
  m_hecEtaMaxA = 3.2 + m_jetDeltaR;

  m_hecPhiMinA = -constants::PI;
  m_hecPhiMaxA = -constants::PI / 2 + m_jetDeltaR;
  
  m_hecEtaMinB = 1.5 - m_jetDeltaR;
  m_hecEtaMaxB = 3.2 + m_jetDeltaR;
  
  m_hecPhiMinB = constants::PI - m_jetDeltaR;
  m_hecPhiMaxB = constants::PI;

  // where plots are drawn from on dphi axis
  m_dPhiZoomLow      = 2 * constants::PI / 3;
  m_dPhiZoomHigh     = constants::PI;

  m_dPhiLogMin       = 5E-4;
  m_dPhiLogMax       = 1;

  // where to fit from 
  m_dPhiFittingMin  = 2.48;
  m_dPhiFittingMax  = constants::PI;

  m_dPhiFittingMinB = 2 * constants::PI / 3;
  m_dPhiFittingMinC = 1.8;
  m_dPhiFittingMinD = 2.6;
  
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // for this uncertainty, we change the fitting range.
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if( m_uncertComp == 22 ){
    double dPhiFittingTmp = m_dPhiFittingMin;
    m_dPhiFittingMin  = m_dPhiFittingMinB; 
    m_dPhiFittingMinB = dPhiFittingTmp;
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
    ( -4.4 )( -3.2)( 0 );
  m_nVarFwdEtaBins = m_varFwdEtaBinning.size() - 1; 

  // --- variable eta/ystar binning ---
  // ystar
  boost::assign::push_back( m_varYstarBinning )
    ( -4.0 )( -2.7 )( -1.8 )( 0 )( 1.8 )( 4.0 );
  m_nVarYstarBins = m_varYstarBinning.size() - 1;

  boost::assign::push_back( m_varYstarBinningFlipped )
    ( -4.0 )( -1.8 )( 0 )( 1.8 )( 2.7 )( 4.0 );

  boost::assign::push_back( m_varEtaBinning )
    ( -4.5 )( -2.7 )( -1.8 )( 0 )( 1.8 )( 4.0 );
  m_nVarEtaBins = m_varEtaBinning.size() - 1;

  boost::assign::push_back( m_varEtaBarrelBinning )
    ( -2.5 )( -1.8 )( 0 )( 1.8 )( 2.5 );
  m_nVarEtaBarrelBins = m_varEtaBarrelBinning.size() - 1;
  
  // --- variable pt binning ---
  boost::assign::push_back( m_varPtBinning )
    ( 28 )( 35 )( 45 )( 90 );
  m_nVarPtBins = m_varPtBinning.size() - 1;

  boost::assign::push_back( m_varRtrkPtBinning )
    ( 20 )( 22 )( 24 )( 26 )( 28 )( 30 )( 32 )( 34 )( 36 )( 38 )
    ( 40 )( 44 )( 48 )( 53 )( 58 )( 64 )( 70 )( 78 )( 90 );
  m_nVarRtrkPtBins = m_varRtrkPtBinning.size() - 1;

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
  m_dPhiWidthMax = 1.0;

  m_dPhiYieldMax = 1.0e-2;
  m_dPhiYieldMin = 5.0e-6;

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

  m_effName        = "eff";
  m_purityName     = "purity";

  m_etaEffName     = m_sEta   + "_" + m_effName;
  m_ystarEffName   = m_sYstar + "_" + m_effName;
  
  m_etaSpectName        = m_sEta + "_" + m_spectName;
  
  m_ystarSpectName      = m_sYstar + "_" + m_spectName;
  m_ystarSpectRecoName  = m_ystarSpectName + "_" + m_recoName;
  m_ystarSpectTruthName = m_ystarSpectName + "_" + m_truthName;

  m_ystarSpectFineName      = m_ystarSpectName + "_" + m_sFine;
  m_ystarSpectFineRecoName  = m_ystarSpectFineName + "_" + m_recoName;
  m_ystarSpectFineTruthName = m_ystarSpectFineName + "_" + m_truthName;

  m_rtrkName = "rtrk";
  
  m_ystarRespMatName           = m_sYstar         + "_" + m_respMatName;  
  
  m_ystarSpectRespMatName      = m_ystarSpectName + "_" + m_respMatName;  
  m_ystarSpectCfactorsName     = m_ystarSpectName + "_" + m_cFactorName;
  m_ystarSpectUnfoldedName     = m_ystarSpectName + "_" + m_unfoldedName; 
  m_ystarSpectRecoUnfoldedName = m_ystarSpectRecoName + "_" + m_unfoldedName;
  
  m_dPhiRecoName   = m_dPhiName + "_" + m_recoName;
  m_dPhiTruthName  = m_dPhiName + "_" + m_truthName;

  m_dPhiCfactorsName    = m_dPhiName + "_" + m_cFactorName;;
  m_dPhiUnfoldedName    = m_dPhiName + "_" + m_unfoldedName;
  
  m_systematicsName     = "systematics";

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

  // MakeSpectTogether( fOut );
  MakeFinalPlotsTogether( fOut, m_widthName );
  MakeFinalPlotsTogether( fOut, m_yieldName );

  // MakeDphiTogether( fOut );
  // CompareWeightIsoPtCuts( fOut );
  
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
    hn->GetAxis(4)->SetTitle(m_sDphi.c_str()); }
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
			       const TLorentzVector*& jet2,
			       bool requirePtCut ){
  
  for( auto& jet : v_jets ){
    if( !jet1 && jet.Pt() > 0 ){
      jet1 = &jet;
    } else if( jet1 && jet.Pt() > 0 ){
      jet2 = &jet;
      break;
    }
  }
  
  // make sure we have two jets
  bool haveBothJets = ( jet1 && jet2 ) ? true : false;

  // if we dont leave and return false
  if( !haveBothJets ){ return false; }

  // if we dont require pt Cut, we can return true now
  if( !requirePtCut ){ return true; }
  
  // otherwise, check pT cut
  double deltaPt = std::abs( ( jet1->Pt() - jet2->Pt() )/1000. ); 

  // we have both jets, return true if they pass cut
  return ( deltaPt > m_deltaPtCut ) ? true : false;
}

bool DiJetAnalysis::PassHECCuts( const TLorentzVector& jet ){

  // for pp, this isnt a problem
  if( !m_is_pPb ){ return true; }
  
  double jetEta = jet.Eta();
  double jetPhi = jet.Phi();

  // now, if second jet is in Pb direction
  // hec region, throw this event out.
  if( ( ( jetEta < m_hecEtaMaxA && jetEta > m_hecEtaMinA ) &&
	( jetPhi < m_hecPhiMaxA && jetPhi > m_hecPhiMinA ) ) ||
      ( ( jetEta < m_hecEtaMaxB && jetEta > m_hecEtaMinB ) &&
	( jetPhi < m_hecPhiMaxB && jetPhi > m_hecPhiMinB ) ) )
    { return false; }

  return true;

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

  // in pPb, quit if there is a subleading
  // jet in Pb going HEC region
  if( !PassHECCuts( *jet2 ) ){ return -1; }
  
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

  // do nothing
  return true;
  
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
      v_jets.at(jn).SetPxPyPzE( 0, 0, 0, -1);
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

TF1* DiJetAnalysis::NormalizeDeltaPhi( TH1* hIn, TH1* hNorm,
				       double xBinCenter, bool comb ){

  TF1* combFit = NULL;

  if( !hIn->GetEntries() ){ return combFit; }

  // for yields, have to rescale later.
  // need this for comb subtraction, however
  hIn->Scale( 1., "width" );
  
  if( comb ){
    combFit = anaTool->SubtractCombinatoric( hIn );
  }

  if( !hNorm ){
    hIn->Scale( 1./hIn->Integral() ) ;
  } else {
    int      xBin = hNorm->FindBin( xBinCenter );
    double  nJets = hNorm->GetBinContent( xBin );
    hIn->Scale( 1./nJets );
    if( combFit ){
      combFit->SetParameter
	( 0, combFit->GetParameter(0) * 1./nJets );
    }
  }

  return combFit;
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

    std::cout << xBin << " - " << vM << " " << eM << " " << vC
	      << " " << eC << " " << newDphi << " " << newDphiError << std::endl; 
    
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


TH1* DiJetAnalysis::FlipOverXaxis( TH1* h, std::vector< double >& binning ){

  std::string title =
    Form( "%s;%s;%s", h->GetTitle(), h->GetXaxis()->GetTitle(), h->GetYaxis()->GetTitle() );
  int nBins = binning.size() - 1;

  std::string name = h->GetName();
  
  h->SetName( Form( "%s_tmp", name.c_str() ) );
  
  TH1* hFlip = new TH1D( name.c_str(), title.c_str(), nBins, &binning[0] );

  for( int i = 1; i <= nBins; i++ ){
    hFlip->SetBinContent( i, h->GetBinContent( nBins - i + 1 ) );
    hFlip->SetBinError  ( i, h->GetBinError  ( nBins - i + 1 ) );
  }

  delete h;
  return hFlip;
}

void DiJetAnalysis::GetSignficance( TGraphAsymmErrors* gNominalR,
				    TGraphAsymmErrors* gSystematicsR ){
  //--------------------------------------------
  // Fit some bins to see the results on ratios
  //--------------------------------------------
  double posYstarFitMin = 1.5;
  double posYstarFitMax = 4.0;
  double c0Pos = 0, dc0Pos = 0, dc1Pos = 0, dc2Pos = 0;
  anaTool->FitPol0Syst
    ( gNominalR, gSystematicsR, c0Pos, dc0Pos, dc1Pos, dc2Pos,
      posYstarFitMin, posYstarFitMax );

  std::string posFitResult = Form( "%5.4f +- %5.4f + %5.4f - %5.4f",
				   c0Pos, dc0Pos, dc1Pos, dc2Pos );
	
  double negYstarFitMin = -4.0;
  double negYstarFitMax =  0.0;
  double c0Neg = 0, dc0Neg = 0, dc1Neg = 0, dc2Neg = 0;
  anaTool->FitPol0Syst
    ( gNominalR, gSystematicsR, c0Neg, dc0Neg, dc1Neg, dc2Neg,
      negYstarFitMin, negYstarFitMax );

  std::string negFitResult = Form( "%5.4f +- %5.4f + %5.4f - %5.4f",
				   c0Neg, dc0Neg, dc1Neg, dc2Neg );

  double significancePos =
    ( 1 - c0Pos ) / std::sqrt( std::pow( dc0Pos, 2 ) + std::pow( dc1Pos, 2 ) );
  double significanceNeg =
    ( 1 - c0Neg ) / std::sqrt( std::pow( dc0Neg, 2 ) + std::pow( dc1Neg, 2 ) );
	
  std::cout << "   " << gNominalR->GetName() << std::endl;
  std::cout << "-----------------Positive------------------------" << std::endl;
  std::cout << "   " << posFitResult  << std::endl;
  std::cout << "   sigma : " << significancePos << std::endl;
  std::cout << "-----------------Negative------------------------" << std::endl;
  std::cout << "   " << negFitResult  << std::endl;
  std::cout << "   sigma : " << significanceNeg << std::endl;
  std::cout << "-----------------Negative------------------------" << std::endl;


}

//---------------------------
//   Get Quantities / Plot 
//---------------------------

void DiJetAnalysis::MakeEtaPhiPtMap( std::vector< TH2* >& vSampleMaps,
				     const std::vector< std::string>& vLabels,
				     const std::string& name ){

  uint nSamples = vSampleMaps.size();
  
  TCanvas c_map( "c_map", "c_map", 800, 600 );
  c_map.SetLogz();
  
  for( uint iG = 0; iG < nSamples; iG++ ){
    std::string label = vLabels[iG];

    if( label.compare( m_allName ) ){ continue; }

    TH1* h = vSampleMaps[ iG ];
    
    h->Draw("col");
    styleTool->SetHStyle( h, 0 );

    DrawAtlasRight();
    
    if( name.find("etaPhi") != std::string::npos ){
      double yLabel = 0.66;
      if( m_isData ){ yLabel = 0.57; }
      drawTool->DrawRightLatex( 0.87, yLabel, "28<#it{p}_{T}<35 [GeV]" );

      // draw lines showing excluded regions
      TLine lA( m_hecEtaMinA, m_hecPhiMinA, m_hecEtaMinA, m_hecPhiMaxA );
      TLine rA( m_hecEtaMaxA, m_hecPhiMinA, m_hecEtaMaxA, m_hecPhiMaxA );
      TLine tA( m_hecEtaMinA, m_hecPhiMaxA, m_hecEtaMaxA, m_hecPhiMaxA );
      TLine bA( m_hecEtaMinA, m_hecPhiMinA, m_hecEtaMaxA, m_hecPhiMinA );

      TLine lB( m_hecEtaMinB, m_hecPhiMinB, m_hecEtaMinB, m_hecPhiMaxB );
      TLine rB( m_hecEtaMaxB, m_hecPhiMinB, m_hecEtaMaxB, m_hecPhiMaxB );
      TLine tB( m_hecEtaMinB, m_hecPhiMaxB, m_hecEtaMaxB, m_hecPhiMaxB );
      TLine bB( m_hecEtaMinB, m_hecPhiMinB, m_hecEtaMaxB, m_hecPhiMinB );

      TLine lAA( m_hecEtaMinA + m_jetDeltaR, m_hecPhiMinA,
		 m_hecEtaMinA + m_jetDeltaR, m_hecPhiMaxA - m_jetDeltaR );
      TLine rAA( m_hecEtaMaxA - m_jetDeltaR, m_hecPhiMinA,
		 m_hecEtaMaxA - m_jetDeltaR, m_hecPhiMaxA - m_jetDeltaR );
      TLine tAA( m_hecEtaMinA + m_jetDeltaR, m_hecPhiMaxA - m_jetDeltaR,
		 m_hecEtaMaxA - m_jetDeltaR, m_hecPhiMaxA - m_jetDeltaR );
      TLine bAA( m_hecEtaMinA + m_jetDeltaR, m_hecPhiMinA,
		 m_hecEtaMaxA - m_jetDeltaR, m_hecPhiMinA );

 
      lA.SetLineWidth(3); rA.SetLineWidth(3);
      tA.SetLineWidth(3); bA.SetLineWidth(3);

      lB.SetLineWidth(3); rB.SetLineWidth(3);
      tB.SetLineWidth(3); bB.SetLineWidth(3);

      lAA.SetLineWidth(3); rAA.SetLineWidth(3);
      tAA.SetLineWidth(3); bAA.SetLineWidth(3);

      lAA.SetLineColor(kRed); rAA.SetLineColor(kRed);
      tAA.SetLineColor(kRed); bAA.SetLineColor(kRed);

      lA.Draw(); rA.Draw();
      tA.Draw(); bA.Draw();

      lB.Draw(); rB.Draw();
      tB.Draw(); bB.Draw();

      lAA.Draw(); rAA.Draw();
      tAA.Draw(); bAA.Draw();

      // save and continue or else lines dont draw.
      SaveAsAll( c_map, h->GetName() );
      h->Write();

      continue;
    }
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
  double min = -1;
  
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
      if( min > hSpect->GetMinimum() || min < 0 ){
	min = hSpect->GetMinimum();
      }
      
    } // end loop over xBin
  } // end loop over iG

  // set maxima globally for all spectra hists.
  // easier to compare. Set on log scale.

  max = anaTool->GetLogMaximum( max );
  min = anaTool->GetLogMinimum( min ) / 10;

  if( m_is_pPb ){
    min = 1e2;
  }
  
  double lX0, lY0, lX1, lY1;
  
  //------------------------------------------------
  //------- For Each xAxis Bin, draw an IGs --------
  //------------------------------------------------

  if( m_is_pPb ){ lX0 = 0.45; lY0 = 0.54; lX1 = 0.76; lY1 = 0.67; }
  else          { lX0 = 0.20; lY0 = 0.47; lX1 = 0.40; lY1 = 0.87; }
  
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
    styleTool->SetLegendStyle( &leg, 0.75 );
    
    int style = 1;
    for( uint iG = 0; iG < nSamples; iG++ ){
      std::string label = vLabels[iG];
     
      // dont draw MB triggers. They arent used in spectra anywhere.
      // We only have them to calculate trigger efficiencies.
      if( label.find("_mb_") != std::string::npos ){ continue; }
      
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
      h->SetMinimum( min );
      leg.AddEntry( h, label.c_str() );
    } // end loop over iG
    
    leg.Draw("same");

    DrawAtlasRight();    
    drawTool->DrawRightLatex( 0.87, 0.62, cLabel );

    SaveAsAll( c, Form("%s_%s", name.c_str(), cName.c_str() ) );
  } // end loop over iX
  
  //------------------------------------------------
  //--------- For each IG, draw xAxisBins ----------
  //------------------------------------------------

  lX0 = 0.20; lY0 = 0.22; lX1 = 0.8; lY1 = 0.35; 
  
  for( uint iG = 0; iG < nSamples; iG++ ){
    std::string label = vLabels[iG];

    if( label.compare( m_allName ) ){ continue; };
    
    std::string cName  = label;
    std::string cLabel = label;

    TCanvas c( "c", cLabel.c_str(), 800, 600 );
    c.SetLogy();
    
    TLegend leg( lX0, lY0, lX1, lY1 );
    styleTool->SetLegendStyle( &leg, 0.85 );
    leg.SetNColumns( 2 );
    
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
      h->SetMinimum( min );
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

void DiJetAnalysis::MakeRtrk( std::vector< TH3* >& vSampleRtrk,
			      const std::vector< std::string>& vLabels,
			      const std::string& name ){
    
  if( !vSampleRtrk.size() ){ return; }

  // it should be in every file
  TH3*  hRef = vSampleRtrk[0];
  int nXbins = hRef->GetNbinsX();
  int nYbins = hRef->GetNbinsY();

  uint nSamples = vSampleRtrk.size();
  
  std::vector< std::vector< TH1* > > vRtrk;
  vRtrk.resize( nSamples );

  std::string yAxisTitle = "rTrk";

  TH2* hRtrk = new TH2D( Form( "h_%s", name.c_str() ), "",
			 m_nVarEtaBarrelBins, 0, 1,
			 m_nVarRtrkPtBins, 0, 1 );
  hRtrk->GetXaxis()->Set( m_nVarEtaBarrelBins, &( m_varEtaBarrelBinning[0] ) );
  hRtrk->GetYaxis()->Set( m_nVarRtrkPtBins, &( m_varRtrkPtBinning[0] ) ) ;
  hRtrk->GetXaxis()->SetTitle( hRef->GetXaxis()->GetTitle() );
  hRtrk->GetYaxis()->SetTitle( hRef->GetYaxis()->GetTitle() );
  hRtrk->GetZaxis()->SetTitle( yAxisTitle.c_str() );
  styleTool->SetHStyle( hRtrk, 0 );

  for( auto& s : vSampleRtrk ){
    std::cout << s->GetName() << std::endl;
  }
  
  for( uint iG = 0; iG < nSamples; iG++){

    std::string label = vLabels[iG];
    
    // skip for all other ones 
    if( label.compare( m_allName ) ){ continue; }

    TH3* hRtrkAll = vSampleRtrk[iG];

    // hRtrkAll->RebinY();
    
    // loop over eta
    for( int xBin = 1; xBin <= nXbins; xBin++ ){
      // loop over ept
      for( int yBin = 1; yBin <= nYbins; yBin++ ){
	TH1* hTmp = static_cast<TH1D*>
	  ( hRtrkAll->ProjectionZ( "hTmp", xBin, xBin, yBin, yBin ) );
	hRtrk->SetBinContent( xBin, yBin, hTmp->GetMean() );
	hRtrk->SetBinError  ( xBin, yBin, hTmp->GetMeanError() );
	delete hTmp;
      }
    }
  }

  hRtrk->Write();
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
    
    TLegend leg( 0.71, 0.48, 0.90, 0.68 );
    styleTool->SetLegendStyle( &leg );

    leg.AddEntry( hMeasured, Form( "%s_{Reco}" , typeMeasured.c_str() ) );
    leg.AddEntry( hUnfolded, Form( "%s_{UF}", typeMeasured.c_str() ) );

    //TH1* hRdataMC = NULL;
    
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
      
      leg.AddEntry( hReco, "MC_{Truth}", "f" );
      hReco->Draw ( "hist" );
    }

    double measuredMax = hMeasured->GetBinContent( hMeasured->GetMaximumBin() );
    double unfoldedMax = hUnfolded->GetBinContent( hUnfolded->GetMaximumBin() );

    
    double measuredMin = hMeasured->GetBinContent( hMeasured->GetMinimumBin() );
    double unfoldedMin = hUnfolded->GetBinContent( hUnfolded->GetMinimumBin() );

    
    double maximum =
      measuredMax > unfoldedMax ? 
      measuredMax : unfoldedMax;
    
    double minimum =
      measuredMin < unfoldedMin ? 
      measuredMin : unfoldedMin;

    maximum = anaTool->GetLogMaximum( maximum );
    minimum = anaTool->GetLogMinimum( minimum );

    if( !m_isData ){
      maximum = 10e4;
      minimum = 10e1;
    }
    
    hMeasured->SetMaximum( maximum );
    hTruth   ->SetMaximum( maximum );
    hUnfolded->SetMaximum( maximum );

    hMeasured->SetMinimum( minimum );
    hTruth   ->SetMinimum( minimum );
    hUnfolded->SetMinimum( minimum );

    hMeasured->Draw( "ep same" );
    hUnfolded->Draw( "ep same" );
    
    leg.Draw();
    
    DrawAtlasRight();    
    drawTool->DrawRightLatex
      ( 0.615, 0.715, anaTool->GetLabel( xLow, xUp, axisLabelTex ) );

    // make ratios and draw cfactors
    pad2.cd();

    double legXmin = m_isData ? 0.25 : 0.35;
    double legXmax = m_isData ? 0.90 : 0.90;
    
    TLegend legR( legXmin, 0.82, legXmax, 0.92 );
    styleTool->SetLegendStyle( &legR, 0.9 );
    legR.SetNColumns(3);
	  
    // make sure to set range to rebinned.
    // the c-factors are from rebinned distributions.
    styleTool->SetHStyleRatio( hCfactors );
    hCfactors->SetYTitle( "Corr. Factor" );
    hCfactors->SetXTitle( "#it{p}_{T} [GeV]" );
    hCfactors->SetTitleOffset( 2.3, "x" );

    if( m_isData ){
      hCfactors->SetMaximum( 1.75 );
      hCfactors->SetMinimum( 0.50 );
    } else{
      hCfactors->SetMaximum( 1.15 );
      hCfactors->SetMinimum( 0.55 );
    }
    TH1* hR = static_cast< TH1D* >
      ( hUnfolded->
	Clone( Form( "h_%s_%s_%s_%s_%s", m_sRatio.c_str(),
		     measuredName.c_str(), m_unfoldedName.c_str(),
		     m_allName.c_str(), hTag.c_str())));
    styleTool->SetHStyleRatio( hR, 1 );
    hR->Divide( hMeasured );

    hR->GetXaxis()->SetRangeUser
      ( m_varPtBinning.front(), m_varPtBinning.back() );
    // hR->Draw( "ep same" );

    hCfactors->GetXaxis()->SetRangeUser
      ( m_varPtBinning.front(), m_varPtBinning.back() );
    hCfactors->Draw("ep same");
    
    legR.AddEntry( hCfactors, "Correction Factor (from MC)" );

    /*
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
    */
    
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
				  TFile* fInPerf, 
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

  std::string mcLevel = "";
  if( dPhiName.find( m_recoName ) != std::string::npos ){
    mcLevel = "Reco Level";
  } else if ( dPhiName.find( m_truthName ) != std::string::npos ){
    mcLevel = "Truth Level";
  }
  
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
	( fInPerf->Get( Form( "h_%s_%s_%s_%s", spectName.c_str(), m_sCounts.c_str(),
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
	  
	  TH1* hDphiWidthsTmp = hn->Projection( fAxisI );
	  TH1* hDphiWidths = FlipOverXaxis( hDphiWidthsTmp, m_varYstarBinningFlipped );
	  hDphiWidths->Reset();
	  hDphiWidths->SetName( hNameW.c_str() );
	  hDphiWidths->SetYTitle( m_sWidthTitle.c_str() );
	  hDphiWidths->SetMarkerSize( hDphiWidths->GetMarkerSize() * 1.5 );
	  styleTool->SetHStyle( hDphiWidths, style );
	  vDphiWidthsTemp.push_back( hDphiWidths );
	  vDphiWidths    .push_back( hDphiWidths );

	  TH1* hDphiWidths2Tmp = hn->Projection( fAxisI );
	  TH1* hDphiWidths2 = FlipOverXaxis( hDphiWidths2Tmp, m_varYstarBinningFlipped );
	  hDphiWidths2->Reset();
	  hDphiWidths2->SetName( Form( "%s_2", hNameW.c_str() ) );
	  hDphiWidths2->SetYTitle( m_sWidthTitle.c_str() );
	  hDphiWidths2->SetMarkerSize( hDphiWidths2->GetMarkerSize() * 2 );
	  styleTool->SetHStyle( hDphiWidths2, 5 + style );
	  vDphiWidths.push_back( hDphiWidths2 );

	  TH1* hDphiWidths3Tmp = hn->Projection( fAxisI );
	  TH1* hDphiWidths3 = FlipOverXaxis( hDphiWidths3Tmp, m_varYstarBinningFlipped );
	  hDphiWidths3->Reset();
	  hDphiWidths3->SetName( Form( "%s_2", hNameW.c_str() ) );
	  hDphiWidths3->SetYTitle( m_sWidthTitle.c_str() );
	  hDphiWidths3->SetMarkerSize( hDphiWidths3->GetMarkerSize() * 2 );
	  styleTool->SetHStyle( hDphiWidths3, 5 + style );
	  vDphiWidths.push_back( hDphiWidths3 );

	  TH1* hDphiWidths4Tmp = hn->Projection( fAxisI );
	  TH1* hDphiWidths4 = FlipOverXaxis( hDphiWidths4Tmp, m_varYstarBinningFlipped );
	  hDphiWidths4->Reset();
	  hDphiWidths4->SetName( Form( "%s_2", hNameW.c_str() ) );
	  hDphiWidths4->SetYTitle( m_sWidthTitle.c_str() );
	  hDphiWidths4->SetMarkerSize( hDphiWidths4->GetMarkerSize() * 2 );
	  styleTool->SetHStyle( hDphiWidths4, 5 + style );
	  vDphiWidths.push_back( hDphiWidths4 );
	  
	  TH1* hDphiWidthsStatTmp = hn->Projection( fAxisI );
	  TH1* hDphiWidthsStat = FlipOverXaxis( hDphiWidthsStatTmp, m_varYstarBinningFlipped );
	  hDphiWidthsStat->Reset();
	  hDphiWidthsStat->SetName( Form( "%s_stat", hNameW.c_str() ) );
	  hDphiWidthsStat->SetYTitle( "RMS (#pi - #Delta#phi)" );
	  hDphiWidthsStat->SetMarkerSize( hDphiWidthsStat->GetMarkerSize() * 2 );
	  styleTool->SetHStyle( hDphiWidthsStat, 5 + style );
	  vDphiWidths.push_back( hDphiWidthsStat );

	  TH1* hDphiYieldsTmp = hn->Projection( fAxisI );
	  TH1* hDphiYields = FlipOverXaxis( hDphiYieldsTmp, m_varYstarBinningFlipped );
	  hDphiYields->Reset();
	  hDphiYields->SetName( hNameY.c_str() );
	  hDphiYields->SetYTitle( m_sYieldTitle.c_str() );
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
	    hDphi->SetYTitle( m_sDphiTitle.c_str() );
	    hDphi->SetXTitle(m_sDphi.c_str());
	    hDphi->SetTitle( "" );
	    vHdPhi.push_back( hDphi );

	    // save the per-jet normalized counts histogram.
	    TH1* hDphiUs = static_cast< TH1D* >( hDphi->Clone( Form( "%s_us", hDphiName.c_str() ) ) );
	    styleTool->SetHStyle( hDphiUs, 1 );
	    hDphiUs->SetYTitle( m_sDphiTitle.c_str() );
	    hDphiUs->SetXTitle(m_sDphi.c_str());
	    hDphiUs->SetTitle( "" );
	    vHdPhi.push_back( hDphiUs );
	   
	    // Combinatoric subtraction after scaling by width
	    // then normalize by jet spectra.
	    NormalizeDeltaPhi( hDphiUs, hSpectCounts, 0.5 * ( axis1Up + axis1Low ), false );
	    TF1* combFit = NormalizeDeltaPhi( hDphi  , hSpectCounts, 0.5 * ( axis1Up + axis1Low ), true  );
	    
	    // take final dPhi histogram, and redo
	    // bin width normalization to get per-jet
	    // yields after comb subtaction
	    TH1* hYields = static_cast< TH1D* >( hDphi->Clone( hYieldsName.c_str() ) );
	    vHdPhi.push_back( hYields );
	    // undo the hDphi->Scale( 1., "width" ) part
	    anaTool->UndoWidthScaling( hYields );

	    TH1* hYieldsUs = static_cast< TH1D* >( hDphiUs->Clone( Form( "%s_us" , hYieldsName.c_str() ) ) );
	    vHdPhi.push_back( hYieldsUs );
	    // undo the hDphi->Scale( 1., "width" ) part
	    anaTool->UndoWidthScaling( hYieldsUs );
	    
	    TCanvas c( hDphi->GetName(), hDphi->GetName(), 800, 600 );
	    // c.SetLogy();

	    hDphi->SetMaximum( hDphi->GetMaximum() * 0.5 );
	    hDphi->SetMinimum( -1 * hDphi->GetMaximum() * 0.05 );
	    // hDphi->SetMinimum( m_dPhiLogMin );
	    // hDphi->SetMaximum( m_dPhiLogMax );
	    hDphi  ->Draw("ep same x0");
	    hDphiUs->Draw("ep same x0");

	    // if there was something to fit, draw the fit result
	    if( combFit ){
	      styleTool->SetHStyle( combFit, 1 );
	      combFit->SetLineColor( kRed );
	      vFits.push_back( combFit );
	      combFit->Draw( "same" );
	    }

	    /*
	    TLegend legUs( 0.76, 0.36, 0.86, 0.46 );
	    styleTool->SetLegendStyle( &legUs, 0.85 );
	    */
	    
	    TLegend legUs( 0.53, 0.40, 0.63, 0.50 );
	    styleTool->SetLegendStyle( &legUs, 0.70 );
	    	    
	    legUs.AddEntry( hDphiUs, "NO CS");
	    legUs.AddEntry( hDphi  , "W/CS" );
	    
	    legUs.Draw();
	    
	    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	    TF1* fit   = NULL;
	    TF1* fitUs = NULL;
	    TF1* fit2  = NULL; 
	    TF1* fit3  = NULL; 
	    TF1* fit4  = NULL; 
	    
	    // now fit
	    fit   = anaTool->FitDphi( hDphi  , m_dPhiFittingMin , m_dPhiFittingMax );
	    fitUs = anaTool->FitDphi( hDphiUs, m_dPhiFittingMin , m_dPhiFittingMax );
	    fit2  = anaTool->FitDphi( hDphi  , m_dPhiFittingMinB, m_dPhiFittingMax );
	    fit3  = anaTool->FitDphi( hDphi  , m_dPhiFittingMinC, m_dPhiFittingMax );
	    fit4  = anaTool->FitDphi( hDphi  , m_dPhiFittingMinD, m_dPhiFittingMax );
	    
	    double   tau = fit->GetParameter(1);
	    double sigma = fit->GetParameter(2);

	    double   tauError = fit->GetParError(1);
	    double sigmaError = fit->GetParError(2);
	    
	    double yieldError = 0;
	    double yield      = hYields->IntegralAndError( 0, hDphi->GetNbinsX(), yieldError );
	 
	    double width      = std::sqrt( 2 * tau * tau + sigma * sigma ); 
	    double widthError = std::sqrt( std::pow( 4 * tau   * tauError  , 2 ) +
					   std::pow( 2 * sigma * sigmaError, 2 ) );

	    double   tauUs = fitUs->GetParameter(1);
	    double sigmaUs = fitUs->GetParameter(2);

	    double   tauErrorUs = fitUs->GetParError(1);
	    double sigmaErrorUs = fitUs->GetParError(2);
	    
	    double yieldErrorUs = 0;
	    double yieldUs      = hYieldsUs->IntegralAndError( 0, hDphiUs->GetNbinsX(), yieldErrorUs );
	 
	    double widthUs      = std::sqrt( 2 * tauUs * tauUs + sigmaUs * sigmaUs ); 
	    double widthErrorUs = std::sqrt( std::pow( 4 * tauUs   * tauErrorUs  , 2 ) +
					     std::pow( 2 * sigmaUs * sigmaErrorUs, 2 ) );

	    
	    styleTool->SetHStyle( fit, 0 );
	    vFits.push_back( fit );

	    styleTool->SetHStyle( fitUs, 0 );
	    fitUs->SetLineColor( kRed );
	    fitUs->SetName( Form( "%s_Us", fitUs->GetName() ) );
	    vFits.push_back( fitUs );

	    styleTool->SetHStyle( fit3, 0 );
	    fit2->SetName( Form( "%s_2", fit3->GetName() ) );
	    vFits.push_back( fit2 );

	    styleTool->SetHStyle( fit3, 0 );
	    fit2->SetName( Form( "%s_3", fit3->GetName() ) );
	    vFits.push_back( fit3 );

	    styleTool->SetHStyle( fit4, 0 );
	    fit2->SetName( Form( "%s_4", fit4->GetName() ) );
	    vFits.push_back( fit4 );

	    // set range back to what it should be.
	    // it can be changed in the fit funtion
	    // hDphi->GetXaxis()->SetRange( m_dPhiZoomLowBin, m_dPhiZoomHighBin );
	    hDphi->GetXaxis()->SetRange( 1, m_dPhiZoomHighBin );
 
	    // draw longer fit first;
	    // fit2->Draw("same");
	    fit  ->Draw("same");
	    fitUs->Draw("same");
	    
	    double chi2NDF  = fit->GetChisquare() / fit->GetNDF();
	    double prob     = fit->GetProb();

	    double chi2NDF2 = fit2->GetChisquare() / fit2->GetNDF();
	    double prob2    = fit2->GetProb();	    
	    
	    h_dPhiChiS->Fill( chi2NDF );
	    h_dPhiProb->Fill( prob    );

	    h_dPhiChiS2->Fill( chi2NDF2 );
	    h_dPhiProb2->Fill( prob2    );
	    
	    drawTool->DrawLeftLatex
	      ( 0.19, 0.48, Form( "%s^{CS} = %6.4f #pm %6.4f",
				  m_sYieldTitle.c_str(), yield, yieldError ), 0.7 );
	    drawTool->DrawLeftLatex
	      ( 0.19, 0.42, Form( "%s^{CS} = %6.4f #pm %6.4f",
				  m_sWidthTitle.c_str(), width, widthError ), 0.7 );

	    drawTool->DrawLeftLatex
	      ( 0.19, 0.60, Form( "#color[2]{%s^{No CS} = %6.4f #pm %6.4f}",
				  m_sYieldTitle.c_str(), yieldUs, yieldErrorUs ), 0.7 );
	    drawTool->DrawLeftLatex
	      ( 0.19, 0.54, Form( "#color[2]{%s^{No CS} = %6.4f #pm %6.4f}",
				  m_sWidthTitle.c_str(), widthUs, widthErrorUs ), 0.7 );

	    if( prob < 0.05 ){
	      v_listBadHist.push_back( hDphi->GetName() );
	      v_listBadFits.push_back( fit  ->GetName() );
	    }

	    drawTool->DrawLeftLatex
	      ( 0.53, 0.60, Form( "#Chi^{2}/NDF = %4.2f", chi2NDF ), 0.7 );
	    drawTool->DrawLeftLatex
	      ( 0.53, 0.54, Form( "Prob = %4.2f", prob ), 0.7 );
	    
	    
	    DrawTopLeftLabels
	      ( m_dPP, axis0Low, axis0Up, axis1Low, axis1Up,
		axis2Low, axis2Up, axis3Low, axis3Up, 0.85 );

	    DrawAtlasRight( CT::DrawTools::drawX0, CT::DrawTools::drawY0, 0.85 );
	    if( !m_isData ){
	      drawTool->DrawRightLatex( 0.875, 0.81, mcLevel.c_str() );
	    }
	    
	    SaveAsAll( c, hDphi->GetName() );
	    fit->Write();
	    
	    hDphi->Write();
	    hDphiCounts->Write();

	    hYields->Write();
	    
	    // draw onto common canvas
	    cAll.cd( cAllPadI + 1 );
	    gPad->SetLogy();
	    hDphi->Draw( "ep X0 same" );
	    fit->Draw("same");

	    double scale  = 1.5;
	    double dyL    = 0.08 * scale;
	    double xstart = 0.02 + ( 1 - scale ) * 0.1;
	    double ystart = 1.1  + ( 1 - scale ) * 0.1;

	    double dx = 0.0, dy = 0.;
	    
	    if( cAllPadI == nCol - 1 ){
	      DrawAtlasRight();
	      if( !m_isData ){
		drawTool->DrawRightLatex( 0.875, 0.73, mcLevel.c_str(), scale );
	      }
	    }
	    if( cAllPadI % nCol > 0 ){
	      hDphi->SetYTitle("");
	      dx = 0.1;
	      dy = -0.1;
	    } else {
	      dx = 0.25;
	    }

	    double dyLabel = 0;
	    if( cAllPadI / nCol <  nRow - 1 ){
	      hDphi->SetXTitle("");
	      dy = -0.15;
	    } else {
	      dy = 0.07;
	      dyLabel = -0.2;
	    }

	    hDphi->GetXaxis()->SetTitleOffset(4);
	    hDphi->GetYaxis()->SetTitleOffset(4);
	    
	    drawTool->DrawLeftLatex
	      ( xstart + dx, ystart + dy + dyLabel, CT::AnalysisTools::GetLabel
		( axis0Low, axis0Up, m_dPP->GetAxisLabel(0) ), scale );
	    drawTool->DrawLeftLatex
	      ( xstart + dx, ystart - dyL + dy + dyLabel, CT::AnalysisTools::GetLabel
		( axis1Low, axis1Up, m_dPP->GetAxisLabel(1) ), scale );  
	    drawTool->DrawLeftLatex
	      ( xstart + dx, ystart - 2 * dyL + dy + dyLabel, CT::AnalysisTools::GetLabel
		( axis2Low, axis2Up, m_dPP->GetAxisLabel(2) ), scale );
	    drawTool->DrawLeftLatex
	      ( xstart + dx, ystart - 3 * dyL + dy + dyLabel, CT::AnalysisTools::GetLabel
		( axis3Low, axis3Up, m_dPP->GetAxisLabel(3) ), scale );

	    hDphi->SetTitleSize( (int)(32 * scale), "xyz");
	    hDphi->SetLabelSize( (int)(32 * scale), "xyz");
	    
	    drawTool->DrawLeftLatex
	      ( 0.47 + dx, 0.36 + dy, Form( "Prob = %4.2f", prob ), scale);
	    drawTool->DrawLeftLatex
	      ( 0.42 + dx, 0.28 + dy, Form( "#Chi^{2}/NDF = %4.2f", chi2NDF ), scale );

	    cAllPadI++;
	    
	    if( fit->GetParameter(1) < 0 )
	      { continue; }

	    drawTool->DrawLeftLatex
	      ( 0.19, 0.55, Form( "Yield = %5.3f #pm %5.3f", yield, yieldError ) );
	    drawTool->DrawLeftLatex
	      ( 0.19, 0.48, Form( "RMS = %5.3f #pm %5.3f", width, widthError ) );

	    double   tau2 = fit2->GetParameter(1);
	    double sigma2 = fit2->GetParameter(2);

	    double   tauError2 = fit2->GetParError(1);
	    double sigmaError2 = fit2->GetParError(2);
	    
	    double width2      = std::sqrt( 2 * tau2 * tau2 + sigma2 * sigma2 ); 
	    double widthError2 = std::sqrt( std::pow( 4 * tau2   * tauError2  , 2 ) +
					    std::pow( 2 * sigma2 * sigmaError2, 2 ) );

	    double   tau3 = fit3->GetParameter(1);
	    double sigma3 = fit3->GetParameter(2);

	    double   tauError3 = fit3->GetParError(1);
	    double sigmaError3 = fit3->GetParError(2);
	    
	    double width3      = std::sqrt( 2 * tau3 * tau3 + sigma3 * sigma3 ); 
	    double widthError3 = std::sqrt( std::pow( 4 * tau3   * tauError3  , 2 ) +
					    std::pow( 2 * sigma3 * sigmaError3, 2 ) );

	    double   tau4 = fit4->GetParameter(1);
	    double sigma4 = fit4->GetParameter(2);

	    double   tauError4 = fit4->GetParError(1);
	    double sigmaError4 = fit4->GetParError(2);
	    
	    double width4      = std::sqrt( 2 * tau4 * tau4 + sigma4 * sigma4 ); 
	    double widthError4 = std::sqrt( std::pow( 4 * tau4   * tauError4  , 2 ) +
					    std::pow( 2 * sigma4 * sigmaError4, 2 ) );

	    // Now, put results on histogram
	    std::pair< double, double > rmsAndError =
	      anaTool->GetRMS( hDphi, 0., constants::PI, constants::PI );
	    double widthStat      = rmsAndError.first;
	    double widthStatError = rmsAndError.second;

	    widthStat      = hDphiCounts->GetRMS();
	    widthStatError = hDphiCounts->GetRMSError();
	    
	    std::cout << hDphi->GetName() << std::endl;
	    std::cout << "Tau   = " << tau   << " TauError   = " << tauError   << std::endl;
	    std::cout << "Sigma = " << sigma << " SigmaError = " << sigmaError << std::endl;
	    std::cout << "Width = " << width << " WidthError = " << widthError << std::endl;
	    std::cout << "Yield = " << yield << " YieldError = " << yieldError << std::endl;

	    std::cout << hDphi->GetName() << std::endl;
	    std::cout << "Tau2   = " << tau2   << " TauError2   = " << tauError2   << std::endl;
	    std::cout << "Sigma2 = " << sigma2 << " SigmaError2 = " << sigmaError2 << std::endl;
	    std::cout << "Width2 = " << width2 << " WidthError2 = " << widthError2 << std::endl;

	    
	    hDphiYields->SetBinContent( nAxis3Bins + 1 - axis3Bin, yield      );
	    hDphiYields->SetBinError  ( nAxis3Bins + 1 - axis3Bin, yieldError );

	    // !!!!!!!!!!!!!!!!!!!!!!!!!!!
	    hDphiWidths->SetBinContent( nAxis3Bins + 1 - axis3Bin, width      );
	    hDphiWidths->SetBinError  ( nAxis3Bins + 1 - axis3Bin, widthError );	    

	    hDphiWidths2->SetBinContent( nAxis3Bins + 1 - axis3Bin, width2      );
	    hDphiWidths2->SetBinError  ( nAxis3Bins + 1 - axis3Bin, widthError2 );

	    hDphiWidths3->SetBinContent( nAxis3Bins + 1 - axis3Bin, width3      );
	    hDphiWidths3->SetBinError  ( nAxis3Bins + 1 - axis3Bin, widthError3 );

	    hDphiWidths4->SetBinContent( nAxis3Bins + 1 - axis3Bin, width4      );
	    hDphiWidths4->SetBinError  ( nAxis3Bins + 1 - axis3Bin, widthError4 );
	    
	    hDphiWidthsStat->SetBinContent( nAxis3Bins + 1 - axis3Bin, widthStat      );
	    hDphiWidthsStat->SetBinError  ( nAxis3Bins + 1 - axis3Bin, widthStatError );
	  } // end loop over axis3

	  double scalingFactor = ( axis1Up - axis1Low ) * ( axis2Up - axis2Low );
	  hDphiYields->Scale( 1./scalingFactor, "width" );

	  
	  // check we are in correct ystar and pt bins
	  if( !m_dPP->CorrectPhaseSpace
	      ( std::vector<int>{ axis0Bin, axis1Bin, axis2Bin, 0 } ) )
	    { continue; }

	  // ----------- widths -----------
	  TCanvas cWidthsCmp( "cWidthsCmp", "cWidthsCmp", 800, 800 );

	  TPad pad1("pad1", "", 0.0, 0.40, 1.0, 1.0 );
	  pad1.SetBottomMargin(0.0);
	  pad1.Draw();
	  
	  TPad pad2("pad2", "", 0.0, 0.0, 1.0, 0.39 );
	  pad2.SetTopMargin(0.05);
	  pad2.SetBottomMargin(0.25);
	  pad2.Draw();

	  TLegend legWAll( 0.30, 0.05, 0.69, 0.20 );
	  styleTool->SetLegendStyle( &legWAll );

	  pad1.cd();

	  hDphiWidths->SetMinimum( m_dPhiWidthMin );
	  hDphiWidths->SetMaximum( m_dPhiWidthMax );

	  TH1* hDphiWidthsCmp = static_cast< TH1D* >
	    ( hDphiWidths->Clone( Form( "h_%s_cmp", hNameW.c_str() ) ) );
	  hDphiWidthsCmp->SetYTitle( m_sWidthTitle.c_str() );
	  hDphiWidthsCmp->SetMarkerSize( hDphiWidthsCmp->GetMarkerSize() * 1.5 );
	  vDphiWidths.push_back( hDphiWidthsCmp );
	  
	  styleTool->SetHStyle( hDphiWidthsCmp , 0 );
	  styleTool->SetHStyle( hDphiWidths2   , 1 );
	  styleTool->SetHStyle( hDphiWidths3   , 2 );
	  styleTool->SetHStyle( hDphiWidths4   , 3 );
	  styleTool->SetHStyle( hDphiWidthsStat, 5 );
	   
	  hDphiWidthsCmp ->Draw("epsame X0");
	  hDphiWidths2   ->Draw("epsame X0");
	  hDphiWidthsStat->Draw("epsame X0");
	  // hDphiWidths3  ->Draw("epsame X0");
	  // hDphiWidths4  ->Draw("epsame X0");
	
	  legWAll.AddEntry
	    ( hDphiWidthsCmp,Form( "Fit %2.1f<#Delta#phi<#pi", m_dPhiFittingMin  ) );
	  legWAll.AddEntry
	    ( hDphiWidths2,  Form( "Fit %2.1f<#Delta#phi<#pi", m_dPhiFittingMinB ) );
	  legWAll.AddEntry( hDphiWidthsStat, "Statistical RMS" );

	  /*
	  legWAll.AddEntry
	    ( hDphiWidths3,  Form( "Fit %2.1f<#Delta#phi<#pi", m_dPhiFittingMinC ) );
	  legWAll.AddEntry
	    ( hDphiWidths4,  Form( "Fit %2.1f<#Delta#phi<#pi", m_dPhiFittingMinD ) );
	  */
	  legWAll.Draw();

	  DrawTopLeftLabels
	    ( m_dPP, axis0Low, axis0Up, axis1Low, axis1Up,
	      axis2Low, axis2Up );

	  DrawAtlasRight();
	  
	  pad2.cd();

	  std::string hNameR = hNameW + "_" + m_sRatio;
 
	  TH1* hWR2 = static_cast< TH1D* >( hDphiWidths2->Clone( hNameR.c_str() ) );
	  styleTool->SetHStyle( hWR2, 1 );
	  vWR.push_back( hWR2 );
	  hWR2->Divide( hDphiWidthsCmp );
	  hWR2->SetYTitle( "Red/Blk" );
	  hWR2->SetMinimum( 0.5 );
	  hWR2->SetMaximum( 1.5 );
	  hWR2->Draw("ep X0 same");

	  TH1* hWR3 = static_cast< TH1D* >( hDphiWidths3->Clone( hNameR.c_str() ) );
	  styleTool->SetHStyle( hWR3, 2 );
	  vWR.push_back( hWR3 );
	  hWR3->Divide( hDphiWidthsCmp );
	  // hWR3->Draw("ep X0 same");

	  TH1* hWR4 = static_cast< TH1D* >( hDphiWidths4->Clone( hNameR.c_str() ) );
	  styleTool->SetHStyle( hWR4, 3 );
	  vWR.push_back( hWR4 );
	  hWR4->Divide( hDphiWidthsCmp );
	  // hWR4->Draw("ep X0 same");

	  TAxis* axisHwr = hWR2->GetXaxis();
    	  TF1* fWR2 = anaTool->FitPol0( hWR2, axisHwr->GetXmin(), axisHwr->GetXmax() );
	  styleTool->SetHStyle( fWR2, 0 );
	  
	  drawTool->DrawRightLatex
	    ( 0.87, 0.85, Form("%f #pm %f", fWR2->GetParameter(0), fWR2->GetParError(0) ) );
	  drawTool->DrawRightLatex
	    ( 0.87, 0.31,  Form("Prob = %4.2f", fWR2->GetProb() ) );
	  
	  fWR2->Draw("same");

	  line.Draw();
	  lineP25.Draw();
	  lineN25.Draw();

	  // for case uncertainty = 22 (fitting)
	  // now change results where fit values had large statistical
	  // error to fit + error, this is to not double count
	  // the statistical error twice.
	  // Bins with good statistics can be treated as is.
	  if( m_uncertComp == 22 ){
	    TH1* hWRF = static_cast< TH1D* >
	      ( hDphiWidthsCmp->Clone( Form( "%s_2", hNameR.c_str() ) ) );
	    styleTool->SetHStyle( hWRF, 1 );
	    vWR.push_back( hWRF );
	    hWRF->Divide( hDphiWidths2 );
	    
	    TF1* fWRF = anaTool->FitPol0( hWRF, axisHwr->GetXmin(), axisHwr->GetXmax() );
	    styleTool->SetHStyle( fWRF, 0 );
	    fWRF->SetLineColor( kRed );

	    hWRF->Draw( "ep X0 same" );
	    fWRF->Draw( "same" );

	    for( int i = 1; i <= nAxis3Bins; i++ ){
	      /*
	      double v1 = hDphiWidthsCmp->GetBinContent( i );
	      double e1 = hDphiWidthsCmp->GetBinError  ( i );
	      double v2 = hDphiWidths2  ->GetBinContent( i );
	      double e2 = hDphiWidths2  ->GetBinError  ( i );

	      // only treat bins with where one of the fits has
	      // larger error ( > 10% ) of value in this way
	      if( ( e1/v1 ) < 0.1 && ( e2/v2 ) < 0.1 ){ continue;}
	      if( std::abs( v1 - v2 ) < e1 ){ continue; }
	      */	      
	      
	      double ci = std::abs( hWRF->GetBinContent( i ) );
	      double cf = std::abs( fWRF->GetParameter(0) ) +  fWRF->GetParError(0);
	      double widthOld = hDphiWidthsCmp->GetBinContent( i );
	      double widthNew = widthOld * ( cf / ci );
	      std::cout << " ------ " << i << " -------- " << ci << " " << cf << " "
			<< widthOld << " " << widthNew << std::endl;
	      hDphiWidths->SetBinContent( i, widthNew );
	    }
	    pad1.cd();
	    styleTool->SetHStyle( hDphiWidths, 2 );
	    hDphiWidths->Draw( "ep X0 same" );
	  }

	  SaveAsAll( cWidthsCmp, hNameW );

	  styleTool->SetHStyle( hDphiWidths, 0 );

	  // ----------- yields -----------
	  TCanvas cYieldsAll( "cYieldsAll", "cYieldsAll",800,600);
	  cYieldsAll.SetLogy();
	  
	  hDphiYields->SetMinimum( m_dPhiYieldMin );
	  hDphiYields->SetMaximum( m_dPhiYieldMax );
	
	  hDphiYields->Draw("ep same X0");

	  DrawTopLeftLabels
	    ( m_dPP, axis0Low, axis0Up, axis1Low, axis1Up,
	      axis2Low, axis2Up );

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
	  ( m_dPP, axis0Low, axis0Up, axis1Low, axis1Up );

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
	  ( m_dPP, axis0Low, axis0Up, axis1Low, axis1Up );

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
  h_dPhiChiS2->Draw("hist same");
  h_dPhiChiS ->Draw("hist same");
  legChiSProb.Draw();

  
  c.cd(2);
  h_dPhiProb2->SetMaximum
    ( h_dPhiProb->GetMaximum() > h_dPhiProb2->GetMaximum() ?
      h_dPhiProb->GetMaximum() * 1.1 :  h_dPhiProb2->GetMaximum() * 1.1 );
  h_dPhiProb2->Draw("hist same");
  h_dPhiProb ->Draw("hist same");
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
					  TFile* fInPerf,
					  const std::string& spectName ){
  

  std::cout << " --- " << fInData->GetName() << " " << fInMC->GetName() << std::endl;
  
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

  TLine lineP25Shift( 2.6, 1.25, xMax, 1.25 );
  lineP25Shift.SetLineStyle( 2  );
  lineP25Shift.SetLineColor( 12 );
  lineP25Shift.SetLineWidth( 2  );

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

  int nCol = 5, nRow = 6;
    
  TCanvas cAll( "cAll", "cAll", 4800, 3600 );
  cAll.Divide( nCol, nRow, 0, 0 );

  int cAllPadI = 0;
  
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
      ( fInPerf->Get( Form( "h_%s_%s_%s_%s", spectName.c_str(), m_sCounts.c_str(),
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
	    ( fInMC->Get( Form( "h_%s_%s_%s", truthName.c_str(),
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
	  styleTool->SetHStyleRatio( hCfactors, 0 );
	  
	  // for compararison later.
	  // so its all in one file.
	  TCanvas cCFactors( hCfactors->GetName(), hCfactors->GetName(), 800, 600 );
	  cCFactors.SetLogz();
	    
	  hCfactors->Draw("ep X0");
	  hCfactors->SetTitleOffset( 1.1, "x" );
	  
	  DrawTopLeftLabels
	    ( m_dPP, axis0Low, axis0Up, axis1Low, axis1Up,
	      axis2Low, axis2Up, axis3Low, axis3Up );

	  line.Draw();
	  lineP25Shift.Draw();
	  lineN25.Draw();
	  
	  DrawAtlasRight();

	  SaveAsAll( cCFactors, hCfactors->GetName() );
  
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
	  fitUnfolded = anaTool->FitDphi( hUnfolded, m_dPhiFittingMin, m_dPhiFittingMax );
	  
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
	  hTruth->SetMinimum( m_dPhiLogMin );
	  hTruth->SetMaximum( m_dPhiLogMax );

	  
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

	  double width      = std::sqrt( 2 * tau * tau + sigma * sigma ); 
	  double widthError = std::sqrt( std::pow( 4 * tau   * tauError  , 2 ) +
					 std::pow( 2 * sigma * sigmaError, 2 ) );

	  double yieldError = 0;
	  double yield      = hUnfolded->IntegralAndError( 0, hUnfolded->GetNbinsX(), yieldError );
	  
	  drawTool->DrawLeftLatex( 0.19, 0.55, Form( "%s = %5.3f #pm %5.3f",
						     m_sYieldTitle.c_str(), yield, yieldError ) );
	  drawTool->DrawLeftLatex( 0.19, 0.48, Form( "%s = %5.3f #pm %5.3f",
						     m_sWidthTitle.c_str(), width, widthError ) );

	  pad2.cd();

	  TLegend legR( 0.2, 0.79, 0.6, 0.89 );
	  styleTool->SetLegendStyle( &legR, 0.9 );
	  legR.SetNColumns(3);
	  
	  // make sure to set range to rebinned.
	  // the c-factors are from rebinned distributions.
	  styleTool->SetHStyleRatio( hCfactors );
	  hCfactors->SetYTitle( "Corr. Factor" );
	  hCfactors->GetXaxis()->SetRange
	    ( m_dPhiRebinnedZoomLowBin, m_dPhiRebinnedZoomHighBin );
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
	  
	  legR.AddEntry( hCfactors, "Correction Factors (from MC)" );
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

	  // save to common canvas
	  
	  // draw onto common canvas
	  cAll.cd( cAllPadI + 1 );
	  hCfactors->SetMaximum( 1.8 );
	  hCfactors->Draw("ep X0");

	  double scale  = 1.5;
	  double dyL    = 0.08 * scale;
	  double xstart = 0.02 + ( 1 - scale ) * 0.1;
	  double ystart = 1.1  + ( 1 - scale ) * 0.1;

	  double dx = 0.0, dy = 0.;

	  hCfactors->GetXaxis()->SetTitleOffset(4);
	  hCfactors->GetYaxis()->SetTitleOffset(4);

	  hCfactors->SetTitleSize( (int)(32 * scale), "xyz");
	  hCfactors->SetLabelSize( (int)(32 * scale), "xyz");
	  
	  if( cAllPadI == nCol - 1 ){
	    if( m_is_pPb ){
	      drawTool->DrawRightLatex( 0.87, 0.73, "#it{p}+Pb Pythia8", scale);
	    } else {
	      drawTool->DrawRightLatex( 0.87, 0.73, "#it{pp} Pythia8", scale );
	    }
	  }
	  if( cAllPadI % nCol > 0 ){
	    hCfactors->SetYTitle("");
	    dx = 0.1;
	    dy = -0.1;
	  } else {
	    dx = 0.25;
	  }

	  double dyLabel = 0;
	  if( cAllPadI / nCol <  nRow - 1 ){
	    hCfactors->SetXTitle("");
	    dy = -0.15;
	  } else {
	    dy = 0.07;
	    dyLabel = -0.2;
	  }
	    
	  drawTool->DrawLeftLatex
	    ( xstart + dx, ystart + dy + dyLabel, CT::AnalysisTools::GetLabel
	      ( axis0Low, axis0Up, m_dPP->GetAxisLabel(0) ), scale );
	  drawTool->DrawLeftLatex
	    ( xstart + dx, ystart - dyL + dy + dyLabel, CT::AnalysisTools::GetLabel
	      ( axis1Low, axis1Up, m_dPP->GetAxisLabel(1) ), scale );  
	  drawTool->DrawLeftLatex
	    ( xstart + dx, ystart - 2 * dyL + dy + dyLabel, CT::AnalysisTools::GetLabel
	      ( axis2Low, axis2Up, m_dPP->GetAxisLabel(2) ), scale );
	  drawTool->DrawLeftLatex
	    ( xstart + dx, ystart - 3 * dyL + dy + dyLabel, CT::AnalysisTools::GetLabel
	      ( axis3Low, axis3Up, m_dPP->GetAxisLabel(3) ), scale );
	    
	  line.Draw();
	  lineP25Shift.Draw();
	  lineN25.Draw();

	  cAllPadI++;

	} // end loop over axis3
      } // end loop over axis2
    } // end loop over axis1     
  } // end loop over axis0

  // SaveAsAll( cAll, Form("h_%s_All", m_dPhiCfactorsName.c_str()) );
  
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

void DiJetAnalysis::CompareWeightIsoPtCuts( TFile* fOut ){
   
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

  TFile* fSpectData = m_is_pPb ?
    TFile::Open( "output/output_pPb_data/myOut_pPb_data_perf_0.root" ) :
    TFile::Open( "output/output_pp_data/myOut_pp_data_perf_0.root" );

  TFile* fSpectMC = m_is_pPb ?
    TFile::Open( "output/output_pPb_mc_pythia8/myOut_pPb_mc_pythia8_perf_0.root" ) :
    TFile::Open( "output/output_pp_mc_pythia8/myOut_pp_mc_pythia8_perf_0.root" );

  TFile* fSpectMCUW = m_is_pPb ?
    TFile::Open( "data/output_pPb_mc_pythia8_uw/myOut_pPb_mc_pythia8_perf_0.root" ) :
    TFile::Open( "data/output_pp_mc_pythia8_uw/myOut_pp_mc_pythia8_perf_0.root" );
  
  TFile* fData = m_is_pPb ?
    TFile::Open( "output/output_pPb_data/myOut_pPb_data_phys_0.root" ) :
    TFile::Open( "output/output_pp_data/myOut_pp_data_phys_0.root" );

  TFile* fMC = m_is_pPb ?
    TFile::Open( "output/output_pPb_mc_pythia8/myOut_pPb_mc_pythia8_phys_0.root" ) :
    TFile::Open( "output/output_pp_mc_pythia8/myOut_pp_mc_pythia8_phys_0.root" );

  TFile* fMCUW = m_is_pPb ?
    TFile::Open( "data/output_pPb_mc_pythia8_uw/myOut_pPb_mc_pythia8_phys_0.root" ) :
    TFile::Open( "data/output_pp_mc_pythia8_uw/myOut_pp_mc_pythia8_phys_0.root" );

  TFile* fDataIso = m_is_pPb ?
    TFile::Open( "data/iso/myOut_pPb_data_phys_UF_0_Iso.root" ) :
    TFile::Open( "data/iso/myOut_pp_data_phys_UF_0_Iso.root" );
  
  TFile* fDataNoIso = m_is_pPb ?
    TFile::Open( "data/iso/myOut_pPb_data_phys_UF_0_NoIso.root" ) : 
    TFile::Open( "data/iso/myOut_pp_data_phys_UF_0_NoIso.root" );
  
  TFile* fData0pT = m_is_pPb ?
    TFile::Open( "data/ptcut/myOut_pPb_data_phys_UF_0_NoIso_0pT.root" ) : 
    TFile::Open( "data/ptcut/myOut_pp_data_phys_UF_0_NoIso_0pT.root" );

  TFile* fData1pT = m_is_pPb ?
    TFile::Open( "data/ptcut/myOut_pPb_data_phys_UF_0_NoIso_1pT.root" ) : 
    TFile::Open( "data/ptcut/myOut_pp_data_phys_UF_0_NoIso_1pT.root" );

  TFile* fData2pT = m_is_pPb ?
    TFile::Open( "data/ptcut/myOut_pPb_data_phys_UF_0_NoIso_2pT.root" ) : 
    TFile::Open( "data/ptcut/myOut_pp_data_phys_UF_0_NoIso_2pT.root" );

  TFile* fData3pT = m_is_pPb ?
    TFile::Open( "data/ptcut/myOut_pPb_data_phys_UF_0_NoIso_3pT.root" ) : 
    TFile::Open( "data/ptcut/myOut_pp_data_phys_UF_0_NoIso_3pT.root" );

  TFile* fMC0pT = m_is_pPb ?
    TFile::Open( "data/ptcut/myOut_pPb_mc_pythia8_phys_UF_0_NoIso_0pT.root" ) : 
    TFile::Open( "data/ptcut/myOut_pp_mc_pythia8_phys_UF_0_NoIso_0pT.root" );

  TFile* fMC1pT = m_is_pPb ?
    TFile::Open( "data/ptcut/myOut_pPb_mc_pythia8_phys_UF_0_NoIso_1pT.root" ) : 
    TFile::Open( "data/ptcut/myOut_pp_mc_pythia8_phys_UF_0_NoIso_1pT.root" );

  TFile* fMC2pT = m_is_pPb ?
    TFile::Open( "data/ptcut/myOut_pPb_mc_pythia8_phys_UF_0_NoIso_2pT.root" ) : 
    TFile::Open( "data/ptcut/myOut_pp_mc_pythia8_phys_UF_0_NoIso_2pT.root" );

  TFile* fMC3pT = m_is_pPb ?
    TFile::Open( "data/ptcut/myOut_pPb_mc_pythia8_phys_UF_0_NoIso_3pT.root" ) : 
    TFile::Open( "data/ptcut/myOut_pp_mc_pythia8_phys_UF_0_NoIso_3pT.root" );

  bool useOverlay = true;
  std::string sOverlayOrNot = useOverlay ? "overlay" : "def"; 

  /*
  TFile* fpPbDataHEC =
    TFile::Open("data/hec/output_pPb_data_def_hec/myOut_pPb_data_phys_UF_0.root" );
  TFile* fpPbDataHEC =
    TFile::Open( Form("data/hec/output_pPb_data_%s_hec/myOut_pPb_data_phys_UF_0.root",
		      sOverlayOrNot.c_str() ) );
  TFile* fpPbDataHEC =
    TFile::Open("output/output_pPb_data/myOut_pPb_data_phys_UF_0.root" );
  TFile* fpPbDataHEC =
    TFile::Open( Form("data/hec/output_pPb_data_%s_hec/myOut_pPb_data_phys_UF_0.root",
		      sOverlayOrNot.c_str() ) );
  TFile* fpPbDataNoHEC =
    TFile::Open( Form("data/hec/output_pPb_data_%s_nohec/myOut_pPb_data_phys_UF_0.root",
		      sOverlayOrNot.c_str() ) );
  TFile* fpPbDataNoHEC =
    TFile::Open( "data/pPb_signal/myOut_pPb_data_phys_UF_0.root" );
  */
  TFile* fpPbDataNoHec =
    TFile::Open("output/output_pPb_data/myOut_pPb_data_phys_UF_0.root" );
  TFile* fpPbDataNoHecShift =
    TFile::Open( "output/output_pPb_data/myOut_pPb_data_phys_UF_P24.root" );

  TFile* fpPbDataNoJes = 
    TFile::Open( "output/output_pPb_data/myOut_pPb_data_phys_UF_0.root" );
  TFile* fpPbDataJes1 =
    TFile::Open("output/output_pPb_data/myOut_pPb_data_phys_UF_P25.root" );
  TFile* fpPbDataJes2 =
    TFile::Open("output/output_pPb_data/myOut_pPb_data_phys_UF_P26.root" );
  TFile* fpPbDataJes3 =
    TFile::Open("output/output_pPb_data/myOut_pPb_data_phys_UF_N25.root" );
  TFile* fpPbDataJes4 =
    TFile::Open("output/output_pPb_data/myOut_pPb_data_phys_UF_N26.root" );

  std::string hMCname   = "h_dPhi_reco_All";
  std::string hDataName = "h_dPhi_All";

  std::string hSpectMCname   = "h_ystar_spect_fine_reco_All";
  std::string hSpectDataName = "h_ystar_spect_fine_All";

  TH1* hProbWeight = new TH1D( "hProbWeight", ";Probability;Count", 20, 0, 1 );
  styleTool->SetHStyle( hProbWeight, 0 );
  
  fOut->cd();

  for( int axis0Bin = 1; axis0Bin <= nAxis0Bins; axis0Bin++ ){
    double axis0Low, axis0Up;
    anaTool->GetBinRange
      ( axis0, axis0Bin, axis0Bin, axis0Low, axis0Up );

    // compare spectra
    std::string hTagSpect =
      Form ( anaTool->GetName( axis0Low, axis0Up, m_dPP->GetAxisName(0) ).c_str() );

    // Now Draw everything.
    TCanvas ccSpectMC( "ccSpectMC", "ccSpectMC", 800, 700 );
    TPad padSpectMC1("padSpectMC1", "", 0.0, 0.40, 1.0, 1.0 );
    padSpectMC1.SetBottomMargin(0);
    padSpectMC1.Draw();
    TPad padSpectMC2("padSpectMC2", "", 0.0, 0.0, 1.0, 0.39 );
    padSpectMC2.SetTopMargin(0.05);
    padSpectMC2.SetBottomMargin(0.25);
    padSpectMC2.Draw();
	  
    TLegend legSpectMC( 0.2, 0.1, 0.4, 0.4 );
    styleTool->SetLegendStyle( &legSpectMC , 0.90 );
	  
    TH1* h_SpectMC    = static_cast< TH1D* >
      ( fSpectMC->Get( Form("%s_%s", hSpectMCname.c_str(), hTagSpect.c_str() ) ) );
    TH1* h_SpectMC_UW = static_cast< TH1D* >
      ( fSpectMCUW->Get( Form("%s_%s", hSpectMCname.c_str(), hTagSpect.c_str() ) ) );
    TH1* h_SpectData  = static_cast< TH1D* >
      ( fSpectData->Get( Form("%s_%s", hSpectDataName.c_str(), hTagSpect.c_str() ) ) );

    h_SpectMC->Rebin(2);
    h_SpectMC_UW->Rebin(2);
    h_SpectData->Rebin(2);
    
    h_SpectData->GetXaxis()->SetRangeUser( m_varPtBinning.front(), m_varPtBinning.back() );

    styleTool->SetHStyle( h_SpectMC, 1 );
    styleTool->SetHStyle( h_SpectMC_UW, 2 );
    styleTool->SetHStyle( h_SpectData, 0 );
	  
    padSpectMC1.cd();
    padSpectMC1.SetLogy();

    int pTbinLow = h_SpectData->FindBin( m_varPtBinning.front() ) + 1;
    int pTbinUp  = h_SpectData->FindBin( m_varPtBinning.back () ) - 1;
    
    double integralData  = h_SpectData->Integral( pTbinLow, pTbinUp );
    double integralMC    = h_SpectMC  ->Integral( pTbinLow, pTbinUp );
    
    double scalingFactor = integralMC / integralData;
    h_SpectData->Scale( scalingFactor );
    double maxSpectData = h_SpectData->GetBinContent( h_SpectData->GetMaximumBin() );
    h_SpectData->SetMaximum( maxSpectData * 2 );
    
    h_SpectData->Draw("ep X0 same");
    h_SpectMC->Draw("ep X0 same");
    h_SpectMC_UW->Draw("ep X0 same");

    legSpectMC.AddEntry( h_SpectData, "Data" );
    legSpectMC.AddEntry( h_SpectMC_UW, "Reco MC - Default" );
    legSpectMC.AddEntry( h_SpectMC, "Reco MC - Re-weighted" );

    drawTool->DrawAtlasInternal();

    DrawTopLeftLabels
      ( m_dPP, axis0Low, axis0Up );

    m_is_pPb ?
      drawTool->DrawRightLatex( 0.87, 0.77, "#it{p}+Pb Pythia8") :
      drawTool->DrawRightLatex( 0.87, 0.77, "#it{pp} Pythia8");
	  
    legSpectMC.Draw();
	  
    padSpectMC2.cd();
    TH1* hR_SpectMC = static_cast< TH1D* >
      ( h_SpectData->Clone
	( Form("%s_%s_%s", hDataName.c_str(), hTagSpect.c_str(), m_sRatio.c_str())));

    hR_SpectMC->GetXaxis()->SetRangeUser( m_varPtBinning.front(), m_varPtBinning.back() );
    
    styleTool->SetHStyleRatio( hR_SpectMC, 1 );
    hR_SpectMC->Divide( h_SpectMC );
    hR_SpectMC->Draw( "ep X0 same" );
    hR_SpectMC->SetYTitle( "Data/MC" );
    hR_SpectMC->SetMinimum( 0.6 );
	  
    TH1* hR_SpectMC_UW = static_cast< TH1D* >
      ( h_SpectData->Clone
	( Form("%s_%s_%s", hDataName.c_str(), hTagSpect.c_str(), m_sRatio.c_str() ) ) );
    styleTool->SetHStyleRatio( hR_SpectMC_UW, 2 );
    hR_SpectMC_UW->Divide( h_SpectMC_UW );
    hR_SpectMC_UW->Draw( "ep X0 same" );
    hR_SpectMC_UW->SetYTitle( "Data/SpectMC" );

    double xSr0 = 0.15;
    double xSr1 = 0.60;
    
    TLegend legSpectR( xSr0, 0.73, xSr1, 0.9 );
    styleTool->SetLegendStyle( &legSpectR, 0.85 );
    legSpectR.AddEntry( hR_SpectMC_UW, "Data/MC Reco (Spectra Weights)" );
    legSpectR.AddEntry( hR_SpectMC, "Data/MC Re-weighted" );
    
    legSpectR.Draw();

    TLine lineSpect( m_varPtBinning.front(), 1,
		     m_varPtBinning.back (), 1 );

    lineSpect.Draw();
    
    SaveAsAll( ccSpectMC, Form("hSpectMC_%s", hTagSpect.c_str() ) );
    
    for( int axis1Bin = 1; axis1Bin <= nAxis1Bins; axis1Bin++ ){
      // check we are in correct ystar and pt bins
      if( !m_dPP->CorrectPhaseSpace
	  ( std::vector<int>{ axis0Bin, axis1Bin, 0, 0 } ) )
	{ continue; }

      double axis1Low, axis1Up;
      anaTool->GetBinRange
	( axis1, axis1Bin, axis1Bin, axis1Low, axis1Up );

      // tag for widths canvas
      std::string hTagC =
	Form ("%s_%s",
	      anaTool->GetName( axis0Low, axis0Up, m_dPP->GetAxisName(0) ).c_str(),
	      anaTool->GetName( axis1Low, axis1Up, m_dPP->GetAxisName(1) ).c_str() );

      for( int axis2Bin = 1; axis2Bin <= nAxis2Bins; axis2Bin++ ){
	double axis2Low , axis2Up;
	anaTool->GetBinRange
	  ( axis2, axis2Bin, axis2Bin, axis2Low, axis2Up );

	if( !m_dPP->CorrectPhaseSpace
	    ( std::vector<int>{ axis0Bin, axis1Bin, axis2Bin, 0 } ) )
	  { continue; }
	
	// get widths histos
	std::string hTag =
	  Form("%s_%s_%s",
	       anaTool->GetName( axis0Low, axis0Up, m_dPP->GetAxisName(0) ).c_str(),
	       anaTool->GetName( axis1Low, axis1Up, m_dPP->GetAxisName(1) ).c_str(),
	       anaTool->GetName( axis2Low, axis2Up, m_dPP->GetAxisName(2) ).c_str() );

	//---------------------------------------------
	// ------------------ WIDTHS ------------------
	//---------------------------------------------

	// !!! Check Isolation
	TCanvas ccNoIsoW( "ccNoIsoW", "ccNoIsoW", 800, 700 );
	TPad padNoIsoW1("padNoIsoW1", "", 0.0, 0.35, 1.0, 1.0 );
	padNoIsoW1.SetBottomMargin(0);
	padNoIsoW1.Draw();
	TPad padNoIsoW2("padNoIsoW2", "", 0.0, 0.0, 1.0, 0.34 );
	padNoIsoW2.SetTopMargin(0.05);
	padNoIsoW2.SetBottomMargin(0.25);
	padNoIsoW2.Draw();

	TLegend legNoIsoW( 0.4, 0.05, 0.6, 0.25 );
	styleTool->SetLegendStyle( &legNoIsoW );

	TH1* h_DataIsoW = static_cast< TH1D* >
	  ( fDataIso->Get( Form("h_dPhi_unfolded_width_All_%s", hTag.c_str() ) ) );
	styleTool->SetHStyle( h_DataIsoW, 0 );
	
	TH1* h_DataNoIsoW = static_cast< TH1D* >
	  ( fDataNoIso->Get( Form("h_dPhi_unfolded_width_All_%s", hTag.c_str() ) ) );
	styleTool->SetHStyle( h_DataNoIsoW, 5 );
	
	padNoIsoW1.cd();
	h_DataIsoW->Draw("ep X0 same");
	h_DataNoIsoW->Draw("ep X0 same");
	legNoIsoW.AddEntry( h_DataIsoW, "Iso Cut" );
	legNoIsoW.AddEntry( h_DataNoIsoW, "No Iso Cut" );
	legNoIsoW.Draw();

	m_is_pPb ?
	  drawTool->DrawRightLatex( 0.8, 0.7, "#it{p}+Pb" ) :
	  drawTool->DrawRightLatex( 0.8, 0.7, "#it{pp}" );

	drawTool->DrawAtlasInternal();

	DrawTopLeftLabels
	  ( m_dPP, axis0Low, axis0Up, axis1Low, axis1Up,
	    axis2Low, axis2Up );
	
	TH1D* hRNoIsoW = static_cast< TH1D* >
	  ( h_DataIsoW->Clone( Form("%s_NoIsoR", h_DataIsoW->GetName() ) ) );
	styleTool->SetHStyleRatio( hRNoIsoW, 0 );
	hRNoIsoW->GetYaxis()->SetTitle("Ratio");
	hRNoIsoW->Divide( h_DataNoIsoW );
	hRNoIsoW->SetMaximum( 1.1 );
	hRNoIsoW->SetMinimum( 0.9 );

	padNoIsoW2.cd();
	hRNoIsoW->Draw("p histo X0 same");

	line.Draw();
	
	SaveAsAll( ccNoIsoW, hRNoIsoW->GetName() );

	// !!! Check pT cuts
	TCanvas ccPtCutW( "ccPtCutW", "ccPtCutW", 800, 700 );
	TPad padPtCutW1("padPtCutW1", "", 0.0, 0.35, 1.0, 1.0 );
	padPtCutW1.SetBottomMargin(0);
	padPtCutW1.Draw();
	TPad padPtCutW2("padPtCutW2", "", 0.0, 0.0, 1.0, 0.34 );
	padPtCutW2.SetTopMargin(0.05);
	padPtCutW2.SetBottomMargin(0.25);
	padPtCutW2.Draw();

	TLegend legPtCutW( 0.4, 0.05, 0.6, 0.25 );
	styleTool->SetLegendStyle( &legPtCutW );

	TH1* h_DataNoPtCutW = static_cast< TH1D* >
	  ( fData0pT->Get( Form("h_dPhi_unfolded_width_All_%s", hTag.c_str() ) ) );
	styleTool->SetHStyle( h_DataNoPtCutW, 0 );

	std::cout << fData1pT->GetName() << std::endl;
	std::cout <<  Form("h_dPhi_unfolded_width_All_%s", hTag.c_str() ) << std::endl;
	
	TH1* h_Data1PtCutW = static_cast< TH1D* >
	  ( fData1pT->Get( Form("h_dPhi_unfolded_width_All_%s", hTag.c_str() ) ) );
	styleTool->SetHStyle( h_Data1PtCutW, 1 );
	
	TH1* h_Data2PtCutW = static_cast< TH1D* >
	  ( fData2pT->Get( Form("h_dPhi_unfolded_width_All_%s", hTag.c_str() ) ) );
	styleTool->SetHStyle( h_Data2PtCutW, 2 );

	TH1* h_Data3PtCutW = static_cast< TH1D* >
	  ( fData3pT->Get( Form("h_dPhi_unfolded_width_All_%s", hTag.c_str() ) ) );
	styleTool->SetHStyle( h_Data3PtCutW, 3 );

	padPtCutW1.cd();
	h_DataNoPtCutW->Draw("ep X0 same");
	h_Data1PtCutW ->Draw("ep X0 same");
	h_Data2PtCutW ->Draw("ep X0 same");
	h_Data3PtCutW ->Draw("ep X0 same");
	legPtCutW.AddEntry( h_DataNoPtCutW, "#Delta#it{p}_{T} = 0 GeV" );
	legPtCutW.AddEntry( h_Data1PtCutW,  "#Delta#it{p}_{T} = 1 GeV" );
	legPtCutW.AddEntry( h_Data2PtCutW,  "#Delta#it{p}_{T} = 2 GeV" );
	legPtCutW.AddEntry( h_Data3PtCutW,  "#Delta#it{p}_{T} = 3 GeV" );
	legPtCutW.Draw();

	m_is_pPb ?
	  drawTool->DrawRightLatex( 0.8, 0.7, "#it{p}+Pb" ) :
	  drawTool->DrawRightLatex( 0.8, 0.7, "#it{pp}" );

	drawTool->DrawAtlasInternal();

	DrawTopLeftLabels( m_dPP, axis0Low, axis0Up, axis1Low, axis1Up, axis2Low, axis2Up );
	
	TH1D* hR1PtCutW = static_cast< TH1D* >
	  ( h_Data1PtCutW->Clone( Form("%s_1PtCutR", h_Data1PtCutW->GetName() ) ) );
	hR1PtCutW->GetYaxis()->SetTitle("Ratio");
	hR1PtCutW->Divide( h_DataNoPtCutW );
	hR1PtCutW->SetMaximum( 1.5 );
	hR1PtCutW->SetMinimum( 0.5 );

	TH1D* hR2PtCutW = static_cast< TH1D* >
	  ( h_Data2PtCutW->Clone( Form("%s_2PtCutR", h_Data2PtCutW->GetName() ) ) );
	hR2PtCutW->Divide( h_DataNoPtCutW );

	TH1D* hR3PtCutW = static_cast< TH1D* >
	  ( h_Data3PtCutW->Clone( Form("%s_3PtCutR", h_Data3PtCutW->GetName() ) ) );
	hR3PtCutW->Divide( h_DataNoPtCutW );

	padPtCutW2.cd();
	hR1PtCutW->Draw("ep X0 same");
	hR2PtCutW->Draw("ep X0 same");
	hR3PtCutW->Draw("ep X0 same");

	line.Draw();
	
	SaveAsAll( ccPtCutW, hR1PtCutW->GetName() );

	// compare the results when there is no hec default, and a shifted region
	if( m_is_pPb ){
	  TFile* fNoHec   = fpPbDataNoHec;
	  TFile* fNoHecShift = fpPbDataNoHecShift;

	  TH1* hW_NoHec = static_cast< TH1D* >
	    ( fNoHec->Get( Form("h_dPhi_unfolded_width_All_%s", hTag.c_str() ) ) );
	  styleTool->SetHStyle( hW_NoHec, 0 );
	  TH1* hW_noHecShift = static_cast< TH1D* >
	    ( fNoHecShift->Get( Form("h_dPhi_unfolded_width_All_%s", hTag.c_str() ) ) );
	  styleTool->SetHStyle( hW_noHecShift, 1 );

	  TH1* hWR = static_cast< TH1D* >
	    ( hW_NoHec->Clone( Form("h_dPhi_unfolded_width_All_hec_%s", hTag.c_str() ) ) );
	  styleTool->SetHStyleRatio( hWR, 0 );
	  // hWR->GetYaxis()->SetTitle("Include/Exclude");
	  hWR->GetYaxis()->SetTitle("Ratio");
	  // hWR->GetYaxis()->SetTitle("Overlay/Signal");
	  hWR->Divide( hW_noHecShift );
	  hWR->SetMaximum( 1.2 );
	  hWR->SetMinimum( 0.8 );

	  TCanvas ccHecW( "ccHecW", "ccHecW", 800, 700 );
	  TPad padHecW1("padHecW1", "", 0.0, 0.35, 1.0, 1.0 );
	  padHecW1.SetBottomMargin(0);
	  padHecW1.Draw();
	  TPad padHecW2("padHecW2", "", 0.0, 0.0, 1.0, 0.34 );
	  padHecW2.SetTopMargin(0.05);
	  padHecW2.SetBottomMargin(0.25);
	  padHecW2.Draw();

	  padHecW1.cd();
	  hW_NoHec  ->Draw("ep X0 same");
	  hW_noHecShift->Draw("ep X0 same");

	  DrawTopLeftLabels
	    ( m_dPP, axis0Low, axis0Up, axis1Low, axis1Up,
	      axis2Low, axis2Up );

	  DrawAtlasRight();
	  
	  TLegend legWhec( 0.2, 0.12, 0.3, 0.27 );
	  styleTool->SetLegendStyle( &legWhec );
	  // legWhec.AddEntry( hW_NoHec, "Include NoHec" );
	  // legWhec.AddEntry( hW_noHecShift, "Exclude NoHec" );
	  legWhec.AddEntry( hW_NoHec, "Default HEC Region" );
	  legWhec.AddEntry( hW_noHecShift, "Increased HEC Region" );
	  // legWhec.AddEntry( hW_NoHec, "Overlay" );
	  // legWhec.AddEntry( hW_noHecShift, "Signal" );

	  legWhec.Draw();
	  drawTool->DrawLeftLatex( 0.2, 0.07, sOverlayOrNot );
	  
	  padHecW2.cd();
	  hWR->Draw("ep X0");

	  line.Draw();
	  
	  SaveAsAll( ccHecW, Form("h_dPhi_unfolded_width_All_hec_%s", hTag.c_str() ) );
	}
	// compare the results with and without new JES uncertainty
	if( m_is_pPb ){
	  TFile* fNoJes = fpPbDataNoJes;
	  TFile* fJes1   = fpPbDataJes1;
	  TFile* fJes2   = fpPbDataJes2;
	  TFile* fJes3   = fpPbDataJes3;
	  TFile* fJes4   = fpPbDataJes4;

	  TH1* hW_noJes = static_cast< TH1D* >
	    ( fNoJes->Get( Form("h_dPhi_unfolded_width_All_%s", hTag.c_str() ) ) );
	  styleTool->SetHStyle( hW_noJes, 0 );
	  TH1* hW_jes1 = static_cast< TH1D* >
	    ( fJes1->Get( Form("h_dPhi_unfolded_width_All_%s", hTag.c_str() ) ) );
	  styleTool->SetHStyle( hW_jes1, 1 );
	  TH1* hW_jes2 = static_cast< TH1D* >
	    ( fJes2->Get( Form("h_dPhi_unfolded_width_All_%s", hTag.c_str() ) ) );
	  styleTool->SetHStyle( hW_jes2, 2 );
	  TH1* hW_jes3 = static_cast< TH1D* >
	    ( fJes3->Get( Form("h_dPhi_unfolded_width_All_%s", hTag.c_str() ) ) );
	  styleTool->SetHStyle( hW_jes3, 3 );
	  TH1* hW_jes4 = static_cast< TH1D* >
	    ( fJes4->Get( Form("h_dPhi_unfolded_width_All_%s", hTag.c_str() ) ) );
	  styleTool->SetHStyle( hW_jes4, 4 );
	  
	  TH1* hWR1 = static_cast< TH1D* >
	    ( hW_noJes->Clone( Form("h_dPhi_unfolded_width_All_jes_%s", hTag.c_str() ) ) );
	  styleTool->SetHStyleRatio( hWR1, 1 );
	  hWR1->GetYaxis()->SetTitle("Def/Shifted");
	  hWR1->Divide( hW_jes1 );
	  hWR1->SetMaximum( 1.1 );
	  hWR1->SetMinimum( 0.9 );

	  TH1* hWR2 = static_cast< TH1D* >
	    ( hW_noJes->Clone( Form("h_dPhi_unfolded_width_All_jes_%s", hTag.c_str() ) ) );
	  styleTool->SetHStyleRatio( hWR2, 2 );
	  hWR2->Divide( hW_jes2 );
	  
	  TH1* hWR3 = static_cast< TH1D* >
	    ( hW_noJes->Clone( Form("h_dPhi_unfolded_width_All_jes_%s", hTag.c_str() ) ) );
	  styleTool->SetHStyleRatio( hWR3, 3 );
	  hWR3->Divide( hW_jes3 );
	  
	  TH1* hWR4 = static_cast< TH1D* >
	    ( hW_noJes->Clone( Form("h_dPhi_unfolded_width_All_jes_%s", hTag.c_str() ) ) );
	  styleTool->SetHStyleRatio( hWR4, 4 );
	  hWR4->Divide( hW_jes4 );
	  
	  TCanvas ccJesW( "ccJesW", "ccJesW", 800, 700 );
	  TPad padJesW1("padJesW1", "", 0.0, 0.35, 1.0, 1.0 );
	  padJesW1.SetBottomMargin(0);
	  padJesW1.Draw();
	  TPad padJesW2("padJesW2", "", 0.0, 0.0, 1.0, 0.34 );
	  padJesW2.SetTopMargin(0.05);
	  padJesW2.SetBottomMargin(0.25);
	  padJesW2.Draw();

	  padJesW1.cd();
	  hW_noJes->Draw("ep X0 same");
	  hW_jes1 ->Draw("ep X0 same");
	  hW_jes2 ->Draw("ep X0 same");
	  hW_jes3 ->Draw("ep X0 same");
	  hW_jes4 ->Draw("ep X0 same");
	 
	  DrawTopLeftLabels
	    ( m_dPP, axis0Low, axis0Up, axis1Low, axis1Up,
	      axis2Low, axis2Up );

	  DrawAtlasRight();
	  
	  TLegend legWhec( 0.4, 0.10, 0.5, 0.30 );
	  styleTool->SetLegendStyle( &legWhec, 0.8 );
	  legWhec.AddEntry( hW_noJes, "No Correction"   );
	  legWhec.AddEntry( hW_jes1, "JES Correction 1 +" );
	  legWhec.AddEntry( hW_jes2, "JES Correction 2 +" );
	  legWhec.AddEntry( hW_jes3, "JES Correction 1 -" );
	  legWhec.AddEntry( hW_jes4, "JES Correction 2 -" );
	  
	  legWhec.Draw();
	  drawTool->DrawLeftLatex( 0.4, 0.07, sOverlayOrNot );

	  padJesW2.cd();
	  hWR1->Draw("hist p X0 same");
	  hWR2->Draw("hist p X0 same ");
	  hWR3->Draw("hist p X0 same ");
	  hWR4->Draw("hist p X0 same ");
	   
	  line.Draw();
	  
	  SaveAsAll( ccJesW, Form("h_dPhi_unfolded_width_All_jes_%s", hTag.c_str() ) );
	}
	
	//---------------------------------------------
	// ------------------ YIELDS ------------------
	//---------------------------------------------

	// !!! Check Isolation
	TCanvas ccNoIsoY( "ccNoIsoY", "ccNoIsoY", 800, 700 );
	TPad padNoIsoY1("padNoIsoY1", "", 0.0, 0.35, 1.0, 1.0 );
	padNoIsoY1.SetBottomMargin(0);
	padNoIsoY1.Draw();
	TPad padNoIsoY2("padNoIsoY2", "", 0.0, 0.0, 1.0, 0.34 );
	padNoIsoY2.SetTopMargin(0.05);
	padNoIsoY2.SetBottomMargin(0.25);
	padNoIsoY2.Draw();
	
	TLegend legNoIsoY( 0.5, 0.05, 0.7, 0.25 );
	styleTool->SetLegendStyle( &legNoIsoY );

	TH1* h_DataIsoY = static_cast< TH1D* >
	  ( fDataIso->Get( Form("h_dPhi_unfolded_yield_All_%s", hTag.c_str() ) ) );
	styleTool->SetHStyle( h_DataIsoY, 0 );
	
	TH1* h_DataNoIsoY = static_cast< TH1D* >
	  ( fDataNoIso->Get( Form("h_dPhi_unfolded_yield_All_%s", hTag.c_str() ) ) );
	styleTool->SetHStyle( h_DataNoIsoY, 5 );
	
	padNoIsoY1.cd();
	padNoIsoY1.SetLogy();
	h_DataIsoY->Draw("ep X0 same");
	h_DataNoIsoY->Draw("ep X0 same");
	legNoIsoY.AddEntry( h_DataIsoY, "Iso Cut" );
	legNoIsoY.AddEntry( h_DataNoIsoY, "No Iso Cut" );
	legNoIsoY.Draw();

	m_is_pPb ?
	  drawTool->DrawRightLatex( 0.8, 0.7, "#it{p}+Pb" ) :
	  drawTool->DrawRightLatex( 0.8, 0.7, "#it{pp}" );

	drawTool->DrawAtlasInternal();

	DrawTopLeftLabels
	  ( m_dPP, axis0Low, axis0Up, axis1Low, axis1Up,
	    axis2Low, axis2Up );

	TH1D* hRNoIsoY = static_cast< TH1D* >
	  ( h_DataIsoY->Clone( Form("%s_NoIsoR", h_DataIsoY->GetName() ) ) );
	styleTool->SetHStyleRatio( hRNoIsoY, 0 );
	hRNoIsoY->GetYaxis()->SetTitle("Ratio");
	hRNoIsoY->Divide( h_DataNoIsoY );

	padNoIsoY2.cd();
	hRNoIsoY->Draw("ep X0 same");
	hRNoIsoY->SetMaximum( 1.1 );
	hRNoIsoY->SetMinimum( 0.9 );
	
	line.Draw();
	
	SaveAsAll( ccNoIsoY, hRNoIsoY->GetName() );

	// !!! Check pT cuts
	TCanvas ccPtCutY( "ccPtCutY", "ccPtCutY", 800, 700 );
	TPad padPtCutY1("padPtCutY1", "", 0.0, 0.35, 1.0, 1.0 );
	padPtCutY1.SetBottomMargin(0);
	padPtCutY1.Draw();
	TPad padPtCutY2("padPtCutY2", "", 0.0, 0.0, 1.0, 0.34 );
	padPtCutY2.SetTopMargin(0.05);
	padPtCutY2.SetBottomMargin(0.25);
	padPtCutY2.Draw();

	TLegend legPtCutY( 0.4, 0.05, 0.6, 0.25 );
	styleTool->SetLegendStyle( &legPtCutY );

	TH1* h_DataNoPtCutY = static_cast< TH1D* >
	  ( fData0pT->Get( Form("h_dPhi_unfolded_yield_All_%s", hTag.c_str() ) ) );
	styleTool->SetHStyle( h_DataNoPtCutY, 0 );
	
	TH1* h_Data1PtCutY = static_cast< TH1D* >
	  ( fData1pT->Get( Form("h_dPhi_unfolded_yield_All_%s", hTag.c_str() ) ) );
	styleTool->SetHStyle( h_Data1PtCutY, 1 );

	TH1* h_Data2PtCutY = static_cast< TH1D* >
	  ( fData2pT->Get( Form("h_dPhi_unfolded_yield_All_%s", hTag.c_str() ) ) );
	styleTool->SetHStyle( h_Data2PtCutY, 2 );

	TH1* h_Data3PtCutY = static_cast< TH1D* >
	  ( fData3pT->Get( Form("h_dPhi_unfolded_yield_All_%s", hTag.c_str() ) ) );
	styleTool->SetHStyle( h_Data3PtCutY, 3 );

	padPtCutY1.cd();
	padPtCutY1.SetLogy();
	h_DataNoPtCutY->Draw("ep X0 same");
	h_Data1PtCutY ->Draw("ep X0 same");
	h_Data2PtCutY ->Draw("ep X0 same");
	h_Data3PtCutY ->Draw("ep X0 same");
	legPtCutY.AddEntry( h_DataNoPtCutY, "#Delta#it{p}_{T} = 0 GeV" );
	legPtCutY.AddEntry( h_Data1PtCutY,  "#Delta#it{p}_{T} = 1 GeV" );
	legPtCutY.AddEntry( h_Data2PtCutY,  "#Delta#it{p}_{T} = 2 GeV" );
	legPtCutY.AddEntry( h_Data3PtCutY,  "#Delta#it{p}_{T} = 3 GeV" );
	legPtCutY.Draw();

	m_is_pPb ?
	  drawTool->DrawRightLatex( 0.8, 0.7, "#it{p}+Pb" ) :
	  drawTool->DrawRightLatex( 0.8, 0.7, "#it{pp}"   );

	drawTool->DrawAtlasInternal();

	DrawTopLeftLabels( m_dPP, axis0Low, axis0Up, axis1Low, axis1Up, axis2Low, axis2Up );
	
	TH1D* hR1PtCutY = static_cast< TH1D* >
	  ( h_Data1PtCutY->Clone( Form("%s_1PtCutR", h_Data1PtCutY->GetName() ) ) );
	hR1PtCutY->GetYaxis()->SetTitle("Ratio");
	hR1PtCutY->Divide( h_DataNoPtCutY );
	hR1PtCutY->SetMaximum( 1.1 );
	hR1PtCutY->SetMinimum( 0.0 );

	TH1D* hR2PtCutY = static_cast< TH1D* >
	  ( h_Data2PtCutY->Clone( Form("%s_2PtCutR", h_Data2PtCutY->GetName() ) ) );
	hR2PtCutY->Divide( h_DataNoPtCutY );
	hR1PtCutY->SetMaximum( 1.1 );
	hR1PtCutY->SetMinimum( 0.0 );

	TH1D* hR3PtCutY = static_cast< TH1D* >
	  ( h_Data3PtCutY->Clone( Form("%s_3PtCutR", h_Data3PtCutY->GetName() ) ) );
	hR3PtCutY->Divide( h_DataNoPtCutY );
	hR1PtCutY->SetMaximum( 1.1 );
	hR1PtCutY->SetMinimum( 0.0 );

	padPtCutY2.cd();
	hR1PtCutY->Draw("ep X0 same");
	hR2PtCutY->Draw("ep X0 same");
	hR3PtCutY->Draw("ep X0 same");

	line.Draw();
	
	SaveAsAll( ccPtCutY, hR1PtCutY->GetName() );

	// compare the results when there is no hec default, and a shifted region
	if( m_is_pPb ){
	  TFile* fNoHec   = fpPbDataNoHec;
	  TFile* fNoHecShift = fpPbDataNoHecShift;

	  TH1* hY_NoHec = static_cast< TH1D* >
	    ( fNoHec->Get( Form("h_dPhi_unfolded_yield_All_%s", hTag.c_str() ) ) );
	  styleTool->SetHStyle( hY_NoHec, 0 );
	  TH1* hY_noHecShift = static_cast< TH1D* >
	    ( fNoHecShift->Get( Form("h_dPhi_unfolded_yield_All_%s", hTag.c_str() ) ) );
	  styleTool->SetHStyle( hY_noHecShift, 1 );

	  TH1* hYR = static_cast< TH1D* >
	    ( hY_NoHec->Clone( Form("h_dPhi_unfolded_yield_All_hec_%s", hTag.c_str() ) ) );
	  styleTool->SetHStyleRatio( hYR, 0 );
	  // hYR->GetYaxis()->SetTitle("Include/Exclude");
	  hYR->GetYaxis()->SetTitle("Ratio");
	  // hYR->GetYaxis()->SetTitle("Overlay/Signal");
	  hYR->Divide( hY_noHecShift );
	  hYR->SetMaximum( 1.2 );
	  hYR->SetMinimum( 0.8 );

	  TCanvas ccHecY( "ccHecY", "ccHecY", 800, 700 );
	  TPad padHecY1("padHecY1", "", 0.0, 0.35, 1.0, 1.0 );
	  padHecY1.SetBottomMargin(0);
	  padHecY1.Draw();
	  padHecY1.SetLogy();
	  TPad padHecY2("padHecY2", "", 0.0, 0.0, 1.0, 0.34 );
	  padHecY2.SetTopMargin(0.05);
	  padHecY2.SetBottomMargin(0.25);
	  padHecY2.Draw();

	  padHecY1.cd();
	  hY_NoHec  ->Draw("ep X0 same");
	  hY_noHecShift->Draw("ep X0 same");

	  DrawTopLeftLabels
	    ( m_dPP, axis0Low, axis0Up, axis1Low, axis1Up,
	      axis2Low, axis2Up );

	  DrawAtlasRight();

	  TLegend legYhec( 0.4, 0.15, 0.5, 0.30 );
	  styleTool->SetLegendStyle( &legYhec );
	  // legYhec.AddEntry( hY_NoHec, "Include NoHec" );
	  // legYhec.AddEntry( hY_noHecShift, "Exclude NoHec" );
	  legYhec.AddEntry( hY_NoHec, "Default HEC Region" );
	  legYhec.AddEntry( hY_noHecShift, "Increased HEC  Region" );
	  // legYhec.AddEntry( hY_NoHec,  "Overlay" );
	  // legYhec.AddEntry( hY_noHecShift, "Signal" );

	  legYhec.Draw();
	  drawTool->DrawLeftLatex( 0.4, 0.07, sOverlayOrNot );

	  padHecY2.cd();
	  hYR->Draw("ep X0");

	  line.Draw();
	  
	  SaveAsAll( ccHecY, Form("h_dPhi_unfolded_yield_All_hec_%s", hTag.c_str() ) );
	}
	// compare the results with and without new JES uncertainty
	if( m_is_pPb ){
	  TFile* fNoJes = fpPbDataNoJes;
	  TFile* fJes1   = fpPbDataJes1;
	  TFile* fJes2   = fpPbDataJes2;
	  TFile* fJes3   = fpPbDataJes3;
	  TFile* fJes4   = fpPbDataJes4;

	  TH1* hY_noJes = static_cast< TH1D* >
	    ( fNoJes->Get( Form("h_dPhi_unfolded_yield_All_%s", hTag.c_str() ) ) );
	  styleTool->SetHStyle( hY_noJes, 0 );
	  TH1* hY_jes1 = static_cast< TH1D* >
	    ( fJes1->Get( Form("h_dPhi_unfolded_yield_All_%s", hTag.c_str() ) ) );
	  styleTool->SetHStyle( hY_jes1, 1 );
	  TH1* hY_jes2 = static_cast< TH1D* >
	    ( fJes2->Get( Form("h_dPhi_unfolded_yield_All_%s", hTag.c_str() ) ) );
	  styleTool->SetHStyle( hY_jes2, 2 );
	  TH1* hY_jes3 = static_cast< TH1D* >
	    ( fJes3->Get( Form("h_dPhi_unfolded_yield_All_%s", hTag.c_str() ) ) );
	  styleTool->SetHStyle( hY_jes3, 3 );
	  TH1* hY_jes4 = static_cast< TH1D* >
	    ( fJes4->Get( Form("h_dPhi_unfolded_yield_All_%s", hTag.c_str() ) ) );
	  styleTool->SetHStyle( hY_jes4, 4 );

	  TH1* hYR1 = static_cast< TH1D* >
	    ( hY_noJes->Clone( Form("h_dPhi_unfolded_yield_All_jes_%s", hTag.c_str() ) ) );
	  styleTool->SetHStyleRatio( hYR1, 1 );
	  hYR1->GetYaxis()->SetTitle("Def/Shifted");
	  hYR1->Divide( hY_jes1 );
	  hYR1->SetMaximum( 1.1 );
	  hYR1->SetMinimum( 0.9 );

	  TH1* hYR2 = static_cast< TH1D* >
	    ( hY_noJes->Clone( Form("h_dPhi_unfolded_yield_All_jes_%s", hTag.c_str() ) ) );
	  styleTool->SetHStyleRatio( hYR2, 2 );
	  hYR2->Divide( hY_jes2 );

	  TH1* hYR3 = static_cast< TH1D* >
	    ( hY_noJes->Clone( Form("h_dPhi_unfolded_yield_All_jes_%s", hTag.c_str() ) ) );
	  styleTool->SetHStyleRatio( hYR3, 3 );
	  hYR3->Divide( hY_jes3 );
	  
	  TH1* hYR4 = static_cast< TH1D* >
	    ( hY_noJes->Clone( Form("h_dPhi_unfolded_yield_All_jes_%s", hTag.c_str() ) ) );
	  styleTool->SetHStyleRatio( hYR4, 4 );
	  hYR4->Divide( hY_jes4 );

	  TCanvas ccJesY( "ccJesY", "ccJesY", 800, 700 );
	  TPad padJesY1("padJesY1", "", 0.0, 0.35, 1.0, 1.0 );
	  padJesY1.SetBottomMargin(0);
	  padJesY1.Draw();
	  TPad padJesY2("padJesY2", "", 0.0, 0.0, 1.0, 0.34 );
	  padJesY2.SetTopMargin(0.05);
	  padJesY2.SetBottomMargin(0.25);
	  padJesY2.Draw();

	  padJesY1.cd();
	  padJesY1.SetLogy();
	  hY_noJes->Draw("ep X0 same");
	  hY_jes1 ->Draw("ep X0 same");
	  hY_jes2 ->Draw("ep X0 same");
	  hY_jes3 ->Draw("ep X0 same");
	  hY_jes4 ->Draw("ep X0 same");
	 
	  DrawTopLeftLabels
	    ( m_dPP, axis0Low, axis0Up, axis1Low, axis1Up,
	      axis2Low, axis2Up );

	  DrawAtlasRight();
	  
	  TLegend legYhec( 0.4, 0.10, 0.5, 0.30 );
	  styleTool->SetLegendStyle( &legYhec, 0.8 );
	  legYhec.AddEntry( hY_noJes, "No Correction"   );
	  legYhec.AddEntry( hY_jes1, "JES Correction 1 +" );
	  legYhec.AddEntry( hY_jes2, "JES Correction 2 +" );
	  legYhec.AddEntry( hY_jes3, "JES Correction 1 -" );
	  legYhec.AddEntry( hY_jes4, "JES Correction 2 -" );
	  
	  legYhec.Draw();
	  drawTool->DrawLeftLatex( 0.4, 0.07, sOverlayOrNot );

	  padJesY2.cd();
	  hYR1->Draw("hist p X0 same");
	  hYR2->Draw("hist p X0 same ");
	  hYR3->Draw("hist p X0 same ");
	  hYR4->Draw("hist p X0 same ");
	   
	  line.Draw();
	  
	  SaveAsAll( ccJesY, Form("h_dPhi_unfolded_yield_All_jes_%s", hTag.c_str() ) );
	}

	// ----- weighting vs y* -----
	TH1* h_MC_ys    =
	  new TH1D( Form("h_MC_ys_%s", hTag.c_str() )   , ";#it{y}_{2}*Integral C_{12} Near #pi;",
		    m_nVarYstarBins, &m_varYstarBinning[0] );
	TH1* h_MC_UW_ys =
	  new TH1D( Form("h_MC_UW_ys_%s", hTag.c_str() ), ";#it{y}_{2}*Integral C_{12} Near #pi;",
		    m_nVarYstarBins, &m_varYstarBinning[0] );
	TH1* h_Data_ys  =
	  new TH1D( Form("h_Data_ys_%s", hTag.c_str() ) , ";#it{y}_{2}*Integral C_{12} Near #pi;",
		    m_nVarYstarBins, &m_varYstarBinning[0] );

	styleTool->SetHStyle( h_MC_ys   , 1 );
	styleTool->SetHStyle( h_MC_UW_ys, 2 );
	styleTool->SetHStyle( h_Data_ys , 0 );

	TCanvas cYS( "cYS", "cYS", 800, 600 );
	TPad padYS1("padYS1", "", 0.0, 0.35, 1.0, 1.0 );
	padYS1.SetBottomMargin(0);
	padYS1.Draw();
	TPad padYS2("padYS2", "", 0.0, 0.0, 1.0, 0.34 );
	padYS2.SetTopMargin(0.05);
	padYS2.SetBottomMargin(0.25);
	padYS2.Draw();

	for( int axis3Bin = 1; axis3Bin <= nAxis3Bins; axis3Bin++ ){
	  // check we are in correct ystar and pt bins
	  if( !m_dPP->CorrectPhaseSpace
	      ( std::vector<int>{ axis0Bin, axis1Bin, axis2Bin, axis3Bin } ) )
	    { continue; }

	  double axis3Low , axis3Up;
	  anaTool->GetBinRange
	    ( axis3, axis3Bin, axis3Bin, axis3Low, axis3Up );

 	  // Now Draw everything.
	  TCanvas ccMC( "ccMC", "ccMC", 800, 700 );
	  TPad padMC1("padMC1", "", 0.0, 0.35, 1.0, 1.0 );
	  padMC1.SetBottomMargin(0);
	  padMC1.Draw();
	  TPad padMC2("padMC2", "", 0.0, 0.0, 1.0, 0.34 );
	  padMC2.SetTopMargin(0.05);
	  padMC2.SetBottomMargin(0.25);
	  padMC2.Draw();
	  
	  std::string hTagDphi =
	    Form("%s_%s_%s_%s",
		 anaTool->GetName( axis0Low, axis0Up, m_dPP->GetAxisName(0) ).c_str(),
		 anaTool->GetName( axis1Low, axis1Up, m_dPP->GetAxisName(1) ).c_str(),
		 anaTool->GetName( axis2Low, axis2Up, m_dPP->GetAxisName(2) ).c_str(),
		 anaTool->GetName( axis3Low, axis3Up, m_dPP->GetAxisName(3) ).c_str() );
	  
	  TLegend legMC( 0.2, 0.4, 0.4, 0.55 );
	  styleTool->SetLegendStyle( &legMC , 0.90 );
	  
	  TH1* h_MC    = static_cast< TH1D* >
	    ( fMC->Get( Form("%s_%s", hMCname.c_str(), hTagDphi.c_str() ) ) );
	  TH1* h_MC_UW = static_cast< TH1D* >
	    ( fMCUW->Get( Form("%s_%s", hMCname.c_str(), hTagDphi.c_str() ) ) );
	  TH1* h_Data  = static_cast< TH1D* >
	    ( fData->Get( Form("%s_%s", hDataName.c_str(), hTagDphi.c_str() ) ) );

	  styleTool->SetHStyle( h_MC, 1 );
	  styleTool->SetHStyle( h_MC_UW, 2 );
	  styleTool->SetHStyle( h_Data, 0 );

	  double maxMC = h_MC->GetBinContent( h_MC->GetMaximumBin() );
	  
	  padMC1.cd();
	  h_Data->SetMaximum( maxMC * 1.5 );
	  h_Data->Draw("ep X0 same");
	  h_MC->Draw("ep X0 same");
	  h_MC_UW->Draw("ep X0 same");

	  legMC.AddEntry( h_Data, "Data" );
	  legMC.AddEntry( h_MC, "MC - Re-Weighted" );
	  legMC.AddEntry( h_MC_UW, "MC - Default" );

	  drawTool->DrawAtlasInternal();

	  DrawTopLeftLabels
	    ( m_dPP, axis0Low, axis0Up, axis1Low, axis1Up,
	      axis2Low, axis2Up, axis3Low, axis3Up );

	  m_is_pPb ?
	    drawTool->DrawRightLatex( 0.87, 0.77, "#it{p}+Pb Pythia8") :
	    drawTool->DrawRightLatex( 0.87, 0.77, "#it{pp} Pythia8");
	  
	  legMC.Draw();
	  
	  padMC2.cd();
	  TH1* hR_MC = static_cast< TH1D* >
	    ( h_Data->Clone
	      ( Form("%s_%s_%s", hDataName.c_str(), hTagDphi.c_str(), m_sRatio.c_str())));
	  
	  styleTool->SetHStyleRatio( hR_MC, 1 );
	  hR_MC->Divide( h_MC );
	  hR_MC->Draw( "ep X0 same" );
	  hR_MC->SetYTitle( "Data/MC" );
	  hR_MC->SetMinimum( 0.25 );
	  
	  TH1* hR_MC_UW = static_cast< TH1D* >
	    ( h_Data->Clone
	      ( Form("%s_%s_%s", hDataName.c_str(), hTagDphi.c_str(), m_sRatio.c_str() ) ) );
	  styleTool->SetHStyleRatio( hR_MC_UW, 2 );
	  hR_MC_UW->Divide( h_MC_UW );
	  hR_MC_UW->Draw( "ep X0 same" );
	  hR_MC_UW->SetYTitle( "Data/MC" );

	  TLegend legR( 0.6, 0.75, 0.8, 0.9 );
	  styleTool->SetLegendStyle( &legR, 0.85 );
	  legR.AddEntry( hR_MC, "Data/MC Reweighted" );
	  legR.AddEntry( hR_MC_UW, "Data/MC Default" );

	  TF1* fR = anaTool->FitPol0( hR_MC, m_dPhiFittingMin , m_dPhiFittingMax );
	  styleTool->SetHStyle( fR , 0 );

	  std::cout << " ------------- " << fR->GetProb() << std::endl; 
	  hProbWeight->Fill( fR->GetProb() );
	  
	  legR.Draw();

	  fR->Draw("same");
		       
	  SaveAsAll( ccMC, Form("hMC_%s_%s", m_dPhiName.c_str(), hTagDphi.c_str() ) );

	  // fill the y* histograms with integral near peak (four bins near peak)

	  double mc_ys_integralError;
	  double mc_uw_ys_integralError;
	  double data_integralError;
	  
	  double mc_ys_integral    = h_MC   ->IntegralAndError
	    ( m_nVarDphiBins - 2, m_nVarDphiBins, mc_ys_integralError    );
	  double mc_uw_ys_integral = h_MC_UW->IntegralAndError
	    ( m_nVarDphiBins - 2, m_nVarDphiBins, mc_uw_ys_integralError );
	  double data_integral     = h_Data ->IntegralAndError
	    ( m_nVarDphiBins - 2, m_nVarDphiBins, data_integralError     );
	  
	  h_MC_ys   ->SetBinContent( axis3Bin, mc_ys_integral    );
	  h_MC_UW_ys->SetBinContent( axis3Bin, mc_uw_ys_integral );
	  h_Data_ys ->SetBinContent( axis3Bin, data_integral     );

	  h_MC_ys   ->SetBinError( axis3Bin, mc_ys_integralError    );
	  h_MC_UW_ys->SetBinError( axis3Bin, mc_uw_ys_integralError );
	  h_Data_ys ->SetBinError( axis3Bin, data_integralError     );
	
	  
	  // -----------------------------------------
	  //    compare mc with different pT cuts
	  // -----------------------------------------

	  // TRUTH
	  TCanvas cTruthMCpT( "cTruthMCpT" ,"cTruthMCpT", 800, 600 );

	  TH1* h_dPhiTruthMCNoPtCut = static_cast< TH1D* >
	    ( fMC0pT->Get( Form("h_dPhi_truth_All_%s", hTagDphi.c_str() ) ) );
	  styleTool->SetHStyle( h_dPhiTruthMCNoPtCut, 0 );
	
	  TH1* h_dPhiTruthMC1PtCut = static_cast< TH1D* >
	    ( fMC1pT->Get( Form("h_dPhi_truth_All_%s", hTagDphi.c_str() ) ) );
	  styleTool->SetHStyle( h_dPhiTruthMC1PtCut, 1 );

	  TH1* h_dPhiTruthMC2PtCut = static_cast< TH1D* >
	    ( fMC2pT->Get( Form("h_dPhi_truth_All_%s", hTagDphi.c_str() ) ) );
	  styleTool->SetHStyle( h_dPhiTruthMC2PtCut, 2 );

	  TH1* h_dPhiTruthMC3PtCut = static_cast< TH1D* >
	    ( fMC3pT->Get( Form("h_dPhi_truth_All_%s", hTagDphi.c_str() ) ) );
	  styleTool->SetHStyle( h_dPhiTruthMC3PtCut, 3 );

	  TF1* f_dPhiTruthMCNoPtCut = static_cast< TF1* >
	    ( fMC0pT->Get( Form("f_h_dPhi_truth_All_%s", hTagDphi.c_str() ) ) );
	  f_dPhiTruthMCNoPtCut->SetLineColor( h_dPhiTruthMCNoPtCut->GetLineColor() );
	  
	  TF1* f_dPhiTruthMC1PtCut = static_cast< TF1* >
	    ( fMC1pT->Get( Form("f_h_dPhi_truth_All_%s", hTagDphi.c_str() ) ) );
	  f_dPhiTruthMC1PtCut->SetLineColor( h_dPhiTruthMC1PtCut->GetLineColor() );
	  
	  TF1* f_dPhiTruthMC2PtCut = static_cast< TF1* >
	    ( fMC2pT->Get( Form("f_h_dPhi_truth_All_%s", hTagDphi.c_str() ) ) );
	  f_dPhiTruthMC2PtCut->SetLineColor( h_dPhiTruthMC2PtCut->GetLineColor() );
	  
	  TF1* f_dPhiTruthMC3PtCut = static_cast< TF1* >
	    ( fMC3pT->Get( Form("f_h_dPhi_truth_All_%s", hTagDphi.c_str() ) ) );
	  f_dPhiTruthMC3PtCut->SetLineColor( h_dPhiTruthMC3PtCut->GetLineColor() );
	  
	  TLegend legDphiTruthMC( 0.18, 0.35, 0.45, 0.55 );
	  styleTool->SetLegendStyle( &legDphiTruthMC );
	  legDphiTruthMC.AddEntry( h_dPhiTruthMCNoPtCut, "#Delta#it{p}_{T} = 0 GeV" );
	  legDphiTruthMC.AddEntry( h_dPhiTruthMC1PtCut , "#Delta#it{p}_{T} = 1 GeV" );
	  legDphiTruthMC.AddEntry( h_dPhiTruthMC2PtCut , "#Delta#it{p}_{T} = 2 GeV" );
	  legDphiTruthMC.AddEntry( h_dPhiTruthMC3PtCut , "#Delta#it{p}_{T} = 3 GeV" );

	  double dPhiTruthMCmax =
	    h_dPhiTruthMCNoPtCut->GetBinContent( h_dPhiTruthMCNoPtCut->GetMaximumBin() );

	  h_dPhiTruthMCNoPtCut->SetMaximum( dPhiTruthMCmax * 1.5 );
	  
	  h_dPhiTruthMCNoPtCut->Draw( "ep X0 same" );
	  h_dPhiTruthMC1PtCut ->Draw( "ep X0 same" );
	  h_dPhiTruthMC2PtCut ->Draw( "ep X0 same" );
	  h_dPhiTruthMC3PtCut ->Draw( "ep X0 same" );

	  f_dPhiTruthMCNoPtCut->Draw( "same" );
	  f_dPhiTruthMC1PtCut ->Draw( "same" );
	  f_dPhiTruthMC2PtCut ->Draw( "same" );
	  f_dPhiTruthMC3PtCut ->Draw( "same" );
 
	  legDphiTruthMC.Draw();
	  
	  drawTool->DrawAtlasInternal();

	  DrawTopLeftLabels
	    ( m_dPP, axis0Low, axis0Up, axis1Low, axis1Up,
	      axis2Low, axis2Up, axis3Low, axis3Up );

	  m_is_pPb ?
	    drawTool->DrawRightLatex( 0.87, 0.79, "Truth #it{p}+Pb Pythia8") :
	    drawTool->DrawRightLatex( 0.87, 0.79, "Truth #it{pp} Pythia8");

	  SaveAsAll( cTruthMCpT, Form("hTruthMCpT_%s_%s", m_dPhiName.c_str(), hTagDphi.c_str() ) );

	  // RECO
	  TCanvas cRecoMCpT( "cRecoMCpT" ,"cRecoMCpT", 800, 600 );

	  TH1* h_dPhiRecoMCNoPtCut = static_cast< TH1D* >
	    ( fMC0pT->Get( Form("h_dPhi_reco_All_%s", hTagDphi.c_str() ) ) );
	  styleTool->SetHStyle( h_dPhiRecoMCNoPtCut, 0 );
	
	  TH1* h_dPhiRecoMC1PtCut = static_cast< TH1D* >
	    ( fMC1pT->Get( Form("h_dPhi_reco_All_%s", hTagDphi.c_str() ) ) );
	  styleTool->SetHStyle( h_dPhiRecoMC1PtCut, 1 );

	  TH1* h_dPhiRecoMC2PtCut = static_cast< TH1D* >
	    ( fMC2pT->Get( Form("h_dPhi_reco_All_%s", hTagDphi.c_str() ) ) );
	  styleTool->SetHStyle( h_dPhiRecoMC2PtCut, 2 );

	  TH1* h_dPhiRecoMC3PtCut = static_cast< TH1D* >
	    ( fMC3pT->Get( Form("h_dPhi_reco_All_%s", hTagDphi.c_str() ) ) );
	  styleTool->SetHStyle( h_dPhiRecoMC3PtCut, 3 );

	  TF1* f_dPhiRecoMCNoPtCut = static_cast< TF1* >
	    ( fMC0pT->Get( Form("f_h_dPhi_reco_All_%s", hTagDphi.c_str() ) ) );
	  f_dPhiRecoMCNoPtCut->SetLineColor( h_dPhiRecoMCNoPtCut->GetLineColor() );
	  
	  TF1* f_dPhiRecoMC1PtCut = static_cast< TF1* >
	    ( fMC1pT->Get( Form("f_h_dPhi_reco_All_%s", hTagDphi.c_str() ) ) );
	  f_dPhiRecoMC1PtCut->SetLineColor( h_dPhiRecoMC1PtCut->GetLineColor() );
	  
	  TF1* f_dPhiRecoMC2PtCut = static_cast< TF1* >
	    ( fMC2pT->Get( Form("f_h_dPhi_reco_All_%s", hTagDphi.c_str() ) ) );
	  f_dPhiRecoMC2PtCut->SetLineColor( h_dPhiRecoMC2PtCut->GetLineColor() );
	  
	  TF1* f_dPhiRecoMC3PtCut = static_cast< TF1* >
	    ( fMC3pT->Get( Form("f_h_dPhi_reco_All_%s", hTagDphi.c_str() ) ) );
	  f_dPhiRecoMC3PtCut->SetLineColor( h_dPhiRecoMC3PtCut->GetLineColor() );
	  
	  TLegend legDphiRecoMC( 0.18, 0.35, 0.45, 0.55 );
	  styleTool->SetLegendStyle( &legDphiRecoMC );
	  legDphiRecoMC.AddEntry( h_dPhiRecoMCNoPtCut, "#Delta#it{p}_{T} = 0 GeV" );
	  legDphiRecoMC.AddEntry( h_dPhiRecoMC1PtCut , "#Delta#it{p}_{T} = 1 GeV" );
	  legDphiRecoMC.AddEntry( h_dPhiRecoMC2PtCut , "#Delta#it{p}_{T} = 2 GeV" );
	  legDphiRecoMC.AddEntry( h_dPhiRecoMC3PtCut , "#Delta#it{p}_{T} = 3 GeV" );

	  double dPhiRecoMCmax =
	    h_dPhiRecoMCNoPtCut->GetBinContent( h_dPhiRecoMCNoPtCut->GetMaximumBin() );

	  h_dPhiRecoMCNoPtCut->SetMaximum( dPhiRecoMCmax * 1.5 );
	  
	  h_dPhiRecoMCNoPtCut->Draw( "ep X0 same" );
	  h_dPhiRecoMC1PtCut ->Draw( "ep X0 same" );
	  h_dPhiRecoMC2PtCut ->Draw( "ep X0 same" );
	  h_dPhiRecoMC3PtCut ->Draw( "ep X0 same" );

	  f_dPhiRecoMCNoPtCut->Draw( "same" );
	  f_dPhiRecoMC1PtCut ->Draw( "same" );
	  f_dPhiRecoMC2PtCut ->Draw( "same" );
	  f_dPhiRecoMC3PtCut ->Draw( "same" );
 
	  legDphiRecoMC.Draw();
	  
	  drawTool->DrawAtlasInternal();

	  DrawTopLeftLabels
	    ( m_dPP, axis0Low, axis0Up, axis1Low, axis1Up,
	      axis2Low, axis2Up, axis3Low, axis3Up );

	  m_is_pPb ?
	    drawTool->DrawRightLatex( 0.87, 0.79, "Reco #it{p}+Pb Pythia8") :
	    drawTool->DrawRightLatex( 0.87, 0.79, "Reco #it{pp} Pythia8");

	  SaveAsAll( cRecoMCpT, Form("hRecoMCpT_%s_%s", m_dPhiName.c_str(), hTagDphi.c_str() ) );

	  // Data RECO
	  TCanvas cRecoDatapT( "cRecoDatapT" ,"cRecoDatapT", 800, 600 );

	  TH1* h_dPhiRecoDataNoPtCut = static_cast< TH1D* >
	    ( fData0pT->Get( Form("h_dPhi_All_%s", hTagDphi.c_str() ) ) );
	  styleTool->SetHStyle( h_dPhiRecoDataNoPtCut, 0 );
	
	  TH1* h_dPhiRecoData1PtCut = static_cast< TH1D* >
	    ( fData1pT->Get( Form("h_dPhi_All_%s", hTagDphi.c_str() ) ) );
	  styleTool->SetHStyle( h_dPhiRecoData1PtCut, 1 );

	  TH1* h_dPhiRecoData2PtCut = static_cast< TH1D* >
	    ( fData2pT->Get( Form("h_dPhi_All_%s", hTagDphi.c_str() ) ) );
	  styleTool->SetHStyle( h_dPhiRecoData2PtCut, 2 );

	  TH1* h_dPhiRecoData3PtCut = static_cast< TH1D* >
	    ( fData3pT->Get( Form("h_dPhi_All_%s", hTagDphi.c_str() ) ) );
	  styleTool->SetHStyle( h_dPhiRecoData3PtCut, 3 );

	  TF1* f_dPhiRecoDataNoPtCut = static_cast< TF1* >
	    ( fData0pT->Get( Form("f_h_dPhi_All_%s", hTagDphi.c_str() ) ) );
	  f_dPhiRecoDataNoPtCut->SetLineColor( h_dPhiRecoDataNoPtCut->GetLineColor() );
	  
	  TF1* f_dPhiRecoData1PtCut = static_cast< TF1* >
	    ( fData1pT->Get( Form("f_h_dPhi_All_%s", hTagDphi.c_str() ) ) );
	  f_dPhiRecoData1PtCut->SetLineColor( h_dPhiRecoData1PtCut->GetLineColor() );
	  
	  TF1* f_dPhiRecoData2PtCut = static_cast< TF1* >
	    ( fData2pT->Get( Form("f_h_dPhi_All_%s", hTagDphi.c_str() ) ) );
	  f_dPhiRecoData2PtCut->SetLineColor( h_dPhiRecoData2PtCut->GetLineColor() );
	  
	  TF1* f_dPhiRecoData3PtCut = static_cast< TF1* >
	    ( fData3pT->Get( Form("f_h_dPhi_All_%s", hTagDphi.c_str() ) ) );
	  f_dPhiRecoData3PtCut->SetLineColor( h_dPhiRecoData3PtCut->GetLineColor() );
	  
	  TLegend legDphiRecoData( 0.18, 0.35, 0.45, 0.55 );
	  styleTool->SetLegendStyle( &legDphiRecoData );
	  legDphiRecoData.AddEntry( h_dPhiRecoDataNoPtCut, "#Delta#it{p}_{T} = 0 GeV" );
	  legDphiRecoData.AddEntry( h_dPhiRecoData1PtCut , "#Delta#it{p}_{T} = 1 GeV" );
	  legDphiRecoData.AddEntry( h_dPhiRecoData2PtCut , "#Delta#it{p}_{T} = 2 GeV" );
	  legDphiRecoData.AddEntry( h_dPhiRecoData3PtCut , "#Delta#it{p}_{T} = 3 GeV" );

	  double dPhiRecoDatamax =
	    h_dPhiRecoDataNoPtCut->GetBinContent( h_dPhiRecoDataNoPtCut->GetMaximumBin() );

	  h_dPhiRecoDataNoPtCut->SetMaximum( dPhiRecoDatamax * 1.5 );
	  
	  h_dPhiRecoDataNoPtCut->Draw( "ep X0 same" );
	  h_dPhiRecoData1PtCut ->Draw( "ep X0 same" );
	  h_dPhiRecoData2PtCut ->Draw( "ep X0 same" );
	  h_dPhiRecoData3PtCut ->Draw( "ep X0 same" );

	  f_dPhiRecoDataNoPtCut->Draw( "same" );
	  f_dPhiRecoData1PtCut ->Draw( "same" );
	  f_dPhiRecoData2PtCut ->Draw( "same" );
	  f_dPhiRecoData3PtCut ->Draw( "same" );
 
	  legDphiRecoData.Draw();
	  
	  drawTool->DrawAtlasInternal();

	  DrawTopLeftLabels
	    ( m_dPP, axis0Low, axis0Up, axis1Low, axis1Up,
	      axis2Low, axis2Up, axis3Low, axis3Up );

	  m_is_pPb ?
	    drawTool->DrawRightLatex( 0.87, 0.79, "Reco Data #it{p}+Pb" ) :
	    drawTool->DrawRightLatex( 0.87, 0.79, "Reco Data #it{pp}"   );

	  SaveAsAll( cRecoDatapT, Form("hRecoDatapT_%s_%s", m_dPhiName.c_str(), hTagDphi.c_str() ) );

	  // Data UNFOLDED
	  TCanvas cUnfoldedDatapT( "cUnfoldedDatapT" ,"cUnfoldedDatapT", 800, 600 );

	  TH1* h_dPhiUnfoldedDataNoPtCut = static_cast< TH1D* >
	    ( fData0pT->Get( Form("h_dPhi_unfolded_All_%s", hTagDphi.c_str() ) ) );
	  styleTool->SetHStyle( h_dPhiUnfoldedDataNoPtCut, 0 );
	
	  TH1* h_dPhiUnfoldedData1PtCut = static_cast< TH1D* >
	    ( fData1pT->Get( Form("h_dPhi_unfolded_All_%s", hTagDphi.c_str() ) ) );
	  styleTool->SetHStyle( h_dPhiUnfoldedData1PtCut, 1 );

	  TH1* h_dPhiUnfoldedData2PtCut = static_cast< TH1D* >
	    ( fData2pT->Get( Form("h_dPhi_unfolded_All_%s", hTagDphi.c_str() ) ) );
	  styleTool->SetHStyle( h_dPhiUnfoldedData2PtCut, 2 );

	  TH1* h_dPhiUnfoldedData3PtCut = static_cast< TH1D* >
	    ( fData3pT->Get( Form("h_dPhi_unfolded_All_%s", hTagDphi.c_str() ) ) );
	  styleTool->SetHStyle( h_dPhiUnfoldedData3PtCut, 3 );

	  TF1* f_dPhiUnfoldedDataNoPtCut = static_cast< TF1* >
	    ( fData0pT->Get( Form("f_h_dPhi_unfolded_All_%s", hTagDphi.c_str() ) ) );
	  f_dPhiUnfoldedDataNoPtCut->SetLineColor( h_dPhiUnfoldedDataNoPtCut->GetLineColor() );
	  
	  TF1* f_dPhiUnfoldedData1PtCut = static_cast< TF1* >
	    ( fData1pT->Get( Form("f_h_dPhi_unfolded_All_%s", hTagDphi.c_str() ) ) );
	  f_dPhiUnfoldedData1PtCut->SetLineColor( h_dPhiUnfoldedData1PtCut->GetLineColor() );
	  
	  TF1* f_dPhiUnfoldedData2PtCut = static_cast< TF1* >
	    ( fData2pT->Get( Form("f_h_dPhi_unfolded_All_%s", hTagDphi.c_str() ) ) );
	  f_dPhiUnfoldedData2PtCut->SetLineColor( h_dPhiUnfoldedData2PtCut->GetLineColor() );
	  
	  TF1* f_dPhiUnfoldedData3PtCut = static_cast< TF1* >
	    ( fData3pT->Get( Form("f_h_dPhi_unfolded_All_%s", hTagDphi.c_str() ) ) );
	  f_dPhiUnfoldedData3PtCut->SetLineColor( h_dPhiUnfoldedData3PtCut->GetLineColor() );
	  
	  TLegend legDphiUnfoldedData( 0.18, 0.35, 0.45, 0.55 );
	  styleTool->SetLegendStyle( &legDphiUnfoldedData );
	  legDphiUnfoldedData.AddEntry( h_dPhiUnfoldedDataNoPtCut, "#Delta#it{p}_{T} = 0 GeV" );
	  legDphiUnfoldedData.AddEntry( h_dPhiUnfoldedData1PtCut , "#Delta#it{p}_{T} = 1 GeV" );
	  legDphiUnfoldedData.AddEntry( h_dPhiUnfoldedData2PtCut , "#Delta#it{p}_{T} = 2 GeV" );
	  legDphiUnfoldedData.AddEntry( h_dPhiUnfoldedData3PtCut , "#Delta#it{p}_{T} = 3 GeV" );

	  double dPhiUnfoldedDatamax =
	    h_dPhiUnfoldedDataNoPtCut->GetBinContent( h_dPhiUnfoldedDataNoPtCut->GetMaximumBin() );

	  h_dPhiUnfoldedDataNoPtCut->SetMaximum( dPhiUnfoldedDatamax * 1.5 );
	  
	  h_dPhiUnfoldedDataNoPtCut->Draw( "ep X0 same" );
	  h_dPhiUnfoldedData1PtCut ->Draw( "ep X0 same" );
	  h_dPhiUnfoldedData2PtCut ->Draw( "ep X0 same" );
	  h_dPhiUnfoldedData3PtCut ->Draw( "ep X0 same" );

	  f_dPhiUnfoldedDataNoPtCut->Draw( "same" );
	  f_dPhiUnfoldedData1PtCut ->Draw( "same" );
	  f_dPhiUnfoldedData2PtCut ->Draw( "same" );
	  f_dPhiUnfoldedData3PtCut ->Draw( "same" );
 
	  legDphiUnfoldedData.Draw();
	  
	  drawTool->DrawAtlasInternal();

	  DrawTopLeftLabels
	    ( m_dPP, axis0Low, axis0Up, axis1Low, axis1Up,
	      axis2Low, axis2Up, axis3Low, axis3Up );

	  m_is_pPb ?
	    drawTool->DrawRightLatex( 0.87, 0.79, "Unfolded Data #it{p}+Pb" ) :
	    drawTool->DrawRightLatex( 0.87, 0.79, "Unfolded Data #it{pp}"   );
                     
	  SaveAsAll( cUnfoldedDatapT, Form("hUnfoldedDatapT_%s_%s", m_dPhiName.c_str(), hTagDphi.c_str() ) );

	} // end loop over axis3

	TLegend legYS( 0.6, 0.4, 0.8, 0.6 );
	styleTool->SetLegendStyle( &legYS, 0.9 );
	
	legYS.AddEntry( h_Data_ys , "Data" );
	legYS.AddEntry( h_MC_ys   , "MC - Re-Weighted" );
	legYS.AddEntry( h_MC_UW_ys, "MC - Default" );
	
	cYS.cd();
	padYS1.cd();

	h_MC_ys->SetMaximum( h_MC_ys->GetMaximum() * 1.5 );
	
	h_MC_ys   ->Draw("ep X0 same");
	h_MC_UW_ys->Draw("ep X0 same");
	h_Data_ys ->Draw("ep X0 same");

	drawTool->DrawAtlasInternal();

	DrawTopLeftLabels
	  ( m_dPP, axis0Low, axis0Up, axis1Low, axis1Up,
	    axis2Low, axis2Up );

	m_is_pPb ?
	  drawTool->DrawRightLatex( 0.87, 0.77, "#it{p}+Pb Pythia8") :
	  drawTool->DrawRightLatex( 0.87, 0.77, "#it{pp} Pythia8");
	  
	legYS.Draw();

	padYS2.cd();
	TH1* hR_DataMC_ys = static_cast< TH1D* >
	  ( h_Data_ys->Clone( Form("hR_MC_ys_%s", hTag.c_str() ) ) );
	styleTool->SetHStyleRatio( hR_DataMC_ys, 1 );
	hR_DataMC_ys->Divide( h_MC_ys );
	hR_DataMC_ys->SetTitle( "Data/MC" );

	TH1* hR_DataMC_UW_ys = static_cast< TH1D* >
	  ( h_Data_ys->Clone( Form("hR_MC_UW_ys_%s", hTag.c_str() ) ) );
	styleTool->SetHStyleRatio( hR_DataMC_UW_ys, 2 );
	hR_DataMC_UW_ys->Divide( h_MC_UW_ys );

	hR_DataMC_ys->SetMaximum(1.3);
	
	hR_DataMC_ys   ->Draw("ep X0 same");
	hR_DataMC_UW_ys->Draw("ep X0 same");

	hR_DataMC_ys   ->Write();
	hR_DataMC_UW_ys->Write();
	
	TLegend legR_ys( 0.6, 0.75, 0.8, 0.9 );
	styleTool->SetLegendStyle( &legR_ys, 0.85 );
	legR_ys.AddEntry( hR_DataMC_ys   , "Data/MC Reweighted" );
	legR_ys.AddEntry( hR_DataMC_UW_ys, "Data/MC Default"  );

	legR_ys.Draw();
		
	SaveAsAll( cYS, Form("hMC_YS_%s_%s", m_dPhiName.c_str(), hTag.c_str() ) );

      } // end loop over axis2
    } // end loop over ystar2
  } // end loop over ystar2

  TCanvas cProb( "cProb", "cProb", 800, 600 );
  hProbWeight->Draw("hist");

  DrawAtlasRightBoth();
  if( m_is_pPb ){
    drawTool->DrawRightLatex( 0.8, 0.8, "#it{p}+Pb" );
  } else {
    drawTool->DrawRightLatex( 0.8, 0.8, "#it{pp}" );
  }

  std::cout << m_is_pPb << std::endl;
  std::string probName = "h_probWeights_";
  probName +=  m_is_pPb ? "pPb" : "pp" ; 
  
  SaveAsAll( cProb, probName );
  hProbWeight->Write();
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

  TFile* fIn_a = TFile::Open( fName_a.c_str() );
  TFile* fIn_b = TFile::Open( fName_b.c_str() );

  /*
  int pTcut = 0;

  std::string fName_pPb =
    Form( "data/myOut_pPb_data_phys_UF_0_NoIso_%dpT.root", pTcut );
  std::string fName_pp =
    Form( "data/myOut_pp_data_phys_UF_0_NoIso_%dpT.root", pTcut );
  
  TFile* fIn_a = TFile::Open( fName_pPb.c_str() );
  TFile* fIn_b = TFile::Open( fName_pp.c_str() );
  */
  
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

  bool isPythia8Closure = !m_isData && m_mcType == 0 ? true : false; 
  
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
      
      TLegend legW( 0.18, 0.03, 0.53, 0.13 + deltaYleg );
      styleTool->SetLegendStyle( &legW, 0.75 );

      TLegend legY( 0.31, 0.03, 0.72, 0.13 + deltaYleg  );
      styleTool->SetLegendStyle( &legY, 0.75 );
      
      int style = 0;
      
      for( int axis2Bin = 1; axis2Bin <= nAxis2Bins; axis2Bin++ ){
	double axis2Low , axis2Up;
	anaTool->GetBinRange
	  ( axis2, axis2Bin, axis2Bin, axis2Low, axis2Up );

	if( !m_dPP->CorrectPhaseSpace
	    ( std::vector<int>{ axis0Bin, axis1Bin, axis2Bin, 0 } ) )
	  { continue; }
	
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

	// for pythia8 MC closure comparisons, zoom in on the ratios
	if( isPythia8Closure ){
	  hW_R->SetMaximum( 1.2 );
	  hW_R->SetMinimum( 0.8 );
	}
	
	vR.push_back( hW_R );
	hW_R->Draw("ep x0 same");

	pad1W.cd();
	
	styleTool->HideAxis( hW_a, "x" );
	styleTool->HideAxis( hW_b, "x" );

	if( hW_a->GetMean() ){
	  legW.AddEntry ( hW_a, Form("%s %s", label_a.c_str(), anaTool->GetLabel
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

	// for pythia8 MC closure comparisons, zoom in on the ratios
	if( isPythia8Closure ){
	  hY_R->SetMaximum( 1.2 );
	  hY_R->SetMinimum( 0.8 );
	}
	
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

	  TCanvas cc( "cc", "cc", 800, 600 );
	  
	  // Now Draw everything.
	  TCanvas ccMC( "ccMC", "ccMC", 800, 700 );
	  TPad padMC1("padMC1", "", 0.0, 0.35, 1.0, 1.0 );
	  padMC1.SetBottomMargin(0);
	  padMC1.Draw();
	  TPad padMC2("padMC2", "", 0.0, 0.0, 1.0, 0.34 );
	  padMC2.SetTopMargin(0.05);
	  padMC2.SetBottomMargin(0.25);
	  padMC2.Draw();
	  
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

	  h_a->GetXaxis()->SetRange( m_dPhiZoomLowBin, m_dPhiZoomHighBin );

	  h_a->GetXaxis()->SetTitle(m_sDphi.c_str());
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

	  TLegend legF( 0.75, 0.22, 0.85, 0.32 );
	  styleTool->SetLegendStyle( &legF , 0.95 );

	  h_a->SetMinimum( m_dPhiLogMin );
	  h_b->SetMinimum( m_dPhiLogMin );

	  h_a->SetMaximum( m_dPhiLogMax );
	  h_b->SetMaximum( m_dPhiLogMax );
	  
	  leg.AddEntry( h_a, label_a.c_str() );
	  leg.AddEntry( h_b, label_b.c_str() );
	  legF.AddEntry( h_a, label_a.c_str() );
	  legF.AddEntry( h_b, label_b.c_str() );
	 
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
	  drawTool->DrawLeftLatex( 0.76, 0.2, m_isData ? "Data" : "MC" );

	  double width_a      = hW_a->GetBinContent( axis3Bin );
	  double widthError_a = hW_a->GetBinError  ( axis3Bin );

	  double width_b      = hW_b->GetBinContent( axis3Bin );
	  double widthError_b = hW_b->GetBinError  ( axis3Bin );

	  drawTool->DrawLeftLatex
	    ( 0.19, 0.54, Form( "%s_{%s}=%4.2f#pm%4.2f", m_sWidthTitle.c_str(),
				label_a.c_str(), width_a, widthError_a ) );
	  drawTool->DrawLeftLatex
	    ( 0.19, 0.46, Form( "%s_{%s}=%4.2f#pm%4.2f", m_sWidthTitle.c_str(),
				label_b.c_str(), width_b, widthError_b ) );

	  double yield_a      = hY_a->GetBinContent( axis3Bin );
	  double yieldError_a = hY_a->GetBinError  ( axis3Bin );

	  double yield_b      = hY_b->GetBinContent( axis3Bin );
	  double yieldError_b = hY_b->GetBinError  ( axis3Bin );

	  std::cout << yield_a << std::endl;
	  
	  drawTool->DrawLeftLatex
	    ( 0.53, 0.79, Form( "%s_{%s}=%e#pm%e", m_sYieldTitle.c_str(),
				label_a.c_str(), yield_a, yieldError_a ) );
	  drawTool->DrawLeftLatex
	    ( 0.53, 0.70, Form( "%s_{%s}=%e#pm%e", m_sYieldTitle.c_str(),
				label_b.c_str(), yield_b, yieldError_b ) );

	  cc.cd();

	  double maxBinA = h_a->GetMaximumBin();
	  double maxBinB = h_b->GetMaximumBin(); 

	  double maxA = h_a->GetBinContent( maxBinA ) * 1.5;
	  double maxB = h_b->GetBinContent( maxBinB ) * 1.5;

	  double maximum = maxA > maxB ? maxA : maxB;
	  
	  h_a->SetMaximum( maximum );
	  
	  h_a->Draw("ep x0 same");
	  h_b->Draw("ep x0 same");
	  
	  f_a->Draw("same");
	  f_b->Draw("same");
	  
	  legF.Draw("same");
	  	  
	  DrawTopLeftLabels
	    ( m_dPP, axis0Low, axis0Up, axis1Low, axis1Up,
	      axis2Low, axis2Up, axis3Low, axis3Up );

	  double labelX0 = m_isData ? 0.80 : 0.78;
	  
	  drawTool->DrawAtlasInternal();
	  drawTool->DrawLeftLatex( labelX0, 0.345, m_isData ? "Data" : "MC" );

	  drawTool->DrawLeftLatex
	    ( 0.19, 0.55, Form( "%s_{%s}=%4.2f #pm %4.2f", m_sWidthTitle.c_str(),
				label_a.c_str(), width_a, widthError_a ) );
	  drawTool->DrawLeftLatex
	    ( 0.19, 0.47, Form( "%s_{%s}=%4.2f #pm %4.2f", m_sWidthTitle.c_str(),
				label_b.c_str(), width_b, widthError_b ) );

	  drawTool->DrawLeftLatex
	    ( 0.50, 0.78, Form( "%s_{%s}=%3.1e #pm %3.1e", m_sYieldTitle.c_str(),
				label_a.c_str(),yield_a, yieldError_a ) );
	  drawTool->DrawLeftLatex
	    ( 0.50, 0.70, Form( "%s_{%s}=%3.1e #pm %3.1e", m_sYieldTitle.c_str(),
				label_b.c_str(), yield_b, yieldError_b ) );

	  SaveAsAll( cc, Form("hf_%s_%s", m_dPhiName.c_str(), hTagDphi.c_str() ), true );
	  
	  styleTool->HideAxis( h_a, "x" );
	  styleTool->HideAxis( h_b, "x" );

	  pad2.cd();
	  h_R->GetXaxis()->SetRange( m_dPhiZoomLowBin, m_dPhiZoomHighBin );
	  h_R->Draw("ep x0 same");

	  dPhiLine.Draw();
	  dPhiLineP25.Draw(); dPhiLineN25.Draw();

	  h_R->Write();
	  SaveAsAll( c , Form("hr_%s_%s", m_dPhiName.c_str(), hTagDphi.c_str() ), true );
	} // end loop over axis3
      } // end loop over axis2

      // back to cW canvas
      cW.cd();
      pad1W.cd();
      
      legW.Draw("same");

      DrawTopLeftLabels( m_dPP, axis0Low, axis0Up, axis1Low, axis1Up );
      DrawAtlasRightBoth();      

      // drawTool->DrawRightLatex( 0.8, 0.8, Form("#Delta#it{p}_{T}=%d GeV", pTcut ) );
      
      pad2W.cd();
      line.Draw();

      if( !isPythia8Closure ){
	lineP25.Draw();
	lineN25.Draw();
      }

      SaveAsPdfPng( cW, Form("h_%s_%s_%s", m_dPhiName.c_str(),
			     m_widthName.c_str(), hTagC.c_str() ), true );
      SaveAsROOT( cW, Form("h_%s_%s_%s", m_dPhiName.c_str(),
			   m_widthName.c_str(), hTagC.c_str() ) );

      // back to cY canvas
      cY.cd();
      pad1Y.cd();
      
      legY.Draw("same");

      DrawTopLeftLabels( m_dPP, axis0Low, axis0Up, axis1Low, axis1Up );
      DrawAtlasRightBoth();

      // drawTool->DrawRightLatex( 0.8, 0.8, Form("#Delta#it{p}_{T}=%d GeV", pTcut ) );

      pad2Y.cd();
      line.Draw();

      if( !isPythia8Closure ){
	lineP25.Draw();
	lineN25.Draw();
      }

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
  h_chi2_a->Draw("hist same");
  h_chi2_b->Draw("hist same");
  DrawAtlasRightBoth();

  c.cd(2);
  h_prob_a->SetMaximum( h_prob_a->GetMaximum() > h_prob_b->GetMaximum() ?
			h_prob_a->GetMaximum() * 1.1 : h_prob_b->GetMaximum() * 1.1 );
  h_prob_a->Draw("hist same");
  h_prob_b->Draw("hist same");

  leg.AddEntry( h_prob_a, label_a.c_str() );
  leg.AddEntry( h_prob_b, label_b.c_str() );

  leg.Draw();
  
  DrawAtlasRightBoth();

  SaveAsPdfPng( c, "h_chi2_prob", true );
  SaveAsROOT  ( c, "h_chi2_prob");
  
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
    ( 0.19, ystart, CT::AnalysisTools::GetLabel
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
    ( 0.19, ystart - 3 * dy, CT::AnalysisTools::GetLabel
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
