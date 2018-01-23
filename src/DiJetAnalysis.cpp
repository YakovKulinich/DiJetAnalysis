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
  
  m_sCounts   = "counts";
  m_sReb      = "reb";
  m_sRatio    = "ratio";

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

  m_dPhiLogMin       = 5E-4;
  m_dPhiLogMax       = 1E-0;
  
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

  m_ptSpectYaxisMin = 1E0;
  
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

  // --- for resp mat ---
  // here we include underflow and overflow pt bins
  m_varPtBinningRespMat = m_varPtBinning;
  m_varPtBinningRespMat.insert( m_varPtBinningRespMat.begin(), 10  );
  m_varPtBinningRespMat.insert( m_varPtBinningRespMat.end()  , 120 );
  m_nVarPtBinsRespMat = m_varPtBinningRespMat.size() - 1;
  
  // --- variable dPhi binning ----
  int nDphiBinsLarge = 4;

  MakeLinearBinning ( m_varDphiBinning, m_varDphiRebinnedBinning, nDphiBinsLarge );
  // m_varDphiRebinnedBinning = m_varDphiBinning;
  // MakeDefaultBinning( m_varDphiBinning, m_varDphiRebinnedBinning, nDphiBinsLarge );
  

  /*
  // for testing rightmost binst
  m_varDphiBinning.push_back( 2.9 );
  m_varDphiBinning.push_back( constants::PI );
  m_varDphiRebinnedBinning = m_varDphiBinning;
  */
  
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
    ( m_nVarDphiBins  )(       2         );
    
  boost::assign::push_back( m_vDphiMin )
    ( 0 )( 0 )( 0 )( 0 )( 0 )( 0 );

  boost::assign::push_back( m_vDphiMax )
    ( 1 )( 1 )( 1 )( 1 )( 1 )( 2 );

  m_nDphiDim     = m_vNdPhiBins.size();
  
  m_dPhiDphiMin  = 0;
  m_dPhiDphiMax  = constants::PI;
  
  m_dPhiWidthMin = 0.0;
  m_dPhiWidthMax = 0.6;

  // default values 0->2 are in
  // SetHStyleRatio.
  // This is for "nonstandard"
  m_ratioMin = 0.0;
  m_ratioMax = 2.0;
  
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
  
  m_recoName       = "reco";
  m_truthName      = "truth";
  m_pairedName     = "paired";
    
  m_cFactorName    = "cFactor";
  m_respMatName    = "respMat";
  m_unfoldedName   = "unfolded";

  m_etaSpectName   = m_sEta  + "_" + m_spectName;
  m_ystarSpectName = m_sYstar+ "_" + m_spectName;

  m_etaSpectRecoName       = m_etaSpectName   + "_" + m_recoName;
  m_ystarSpectRecoName     = m_ystarSpectName + "_" + m_recoName;

  m_etaSpectTruthName      = m_etaSpectName   + "_" + m_truthName;
  m_ystarSpectTruthName    = m_ystarSpectName + "_" + m_truthName;

  m_ystarSpectRespMatName  = m_ystarSpectName + "_" + m_respMatName;
  
  m_etaSpectUnfoldedName   = m_etaSpectName   + "_" + m_unfoldedName;
  m_ystarSpectUnfoldedName = m_ystarSpectName + "_" + m_unfoldedName; 

  m_etaSpectRecoUnfoldedName   = m_etaSpectRecoName   + "_" + m_unfoldedName;
  m_ystarSpectRecoUnfoldedName = m_ystarSpectRecoName + "_" + m_unfoldedName;

  m_ystarSpectCfactorsName = m_ystarSpectName + "_" + m_cFactorName;
  
  m_dPhiRecoName   = m_dPhiName + "_" + m_recoName;
  m_dPhiTruthName  = m_dPhiName + "_" + m_truthName;

  m_dPhiCfactorsName     = m_dPhiName + "_" + m_cFactorName;;
  m_dPhiUnfoldedName     = m_dPhiName + "_" + m_unfoldedName;
  
  m_systematicsName      = "systematics";
  m_dPhiSystematicsName  = m_dPhiName + "_" + m_systematicsName;

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
  m_fNamePerf   = m_fName + "_" + m_sPerf;
  m_fNamePerfUF = m_fName + "_" + m_sPerf + "_" + m_unfoldingFileSuffix;
  m_fNamePhys   = m_fName + "_" + m_sPhys;
  m_fNamePhysUF = m_fName + "_" + m_sPhys + "_" + m_unfoldingFileSuffix;

  m_fNameDefRaw    = m_fNameRaw    + "_0.root";
  m_fNameDefPerf   = m_fNamePerf   + "_0.root";
  m_fNameDefPerfUF = m_fNamePerfUF + "_0.root";
  m_fNameDefPhys   = m_fNamePhys   + "_0.root";
  m_fNameDefPhysUF = m_fNamePhysUF + "_0.root";

  m_fNameRaw       += "_" + m_uncertSuffix + ".root";
  m_fNamePerf      += "_" + m_uncertSuffix + ".root";
  m_fNamePerfUF    += "_" + m_uncertSuffix + ".root";
  m_fNamePhys      += "_" + m_uncertSuffix + ".root";
  m_fNamePhysUF    += "_" + m_uncertSuffix + ".root";

  m_fNameSYS = m_fName + "_" + m_systematicsFileSuffix;
  
  std::cout << "fNameRaw: " << m_fNameRaw << std::endl;

  m_fNameRivetMC = "data/rivet_" + system + ".root";

  // directory of where unfolding MC files are
  std::string unfoldingMCdir =
    Form( "%s/%s_%s_%s_%s",
	  m_sOutput.c_str(), m_sOutput.c_str(), system.c_str(),
	  m_sMC.c_str(), m_mcTypeName.c_str());
  
  // name of file that is used as input for performance unfolding
  m_fNamePerfUnfoldingMC
    = Form( "%s/%s_%s_%s_%s_%s_%s.root",
	    unfoldingMCdir.c_str(), m_myOutName.c_str(),
	    system.c_str(), m_sMC.c_str(), m_mcTypeName.c_str(),
	    m_sPerf.c_str(), m_uncertSuffix.c_str() );

  // name of file that is used as input for physics unfolding
  m_fNamePhysUnfoldingMC
    = Form( "%s/%s_%s_%s_%s_%s_%s.root",
	    unfoldingMCdir.c_str(), m_myOutName.c_str(),
	    system.c_str(), m_sMC.c_str(), m_mcTypeName.c_str(),
	    m_sPhys.c_str(), m_uncertSuffix.c_str() );  
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
  if( hn->GetNdimensions() == 5 ){
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

void DiJetAnalysis::AnalyzeSpectra( TH2* hSpect,
				    const std::vector< TLorentzVector >& vJets){

  // to fill a simple spectra over all jets.
  // if want to do more complicated things (such as in data)
  // do them in ProcessEvents()
  // loop over paired reco jets and fill
  for( auto& jet : vJets ){
    if( jet.Pt() <= 0 ){ continue; }
	
    double jetEta    = jet.Eta();
    double jetPt     = jet.Pt()/1000.;

    double jetWeight = GetJetWeight( jet );

    hSpect->Fill(  jetEta,  jetPt,  jetWeight );	

    // for pp fill both sides
    if( m_is_pPb ){ continue; }

    hSpect->Fill( -jetEta,  jetPt,  jetWeight );	
  }
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
  x[5] = 0;
  hn->Fill( &x[0], jetWeight );

  x[5] = 1;
  hn->Fill( &x[0], jetWeight * deltaPhi );
  
  // for pp, fill twice. once for each side since
  // it is symmetric in pp. For pPb, continue
  if( m_is_pPb ){ return deltaPhi; }
  
  x[0] = -jet1_ystar;  
  x[1] = -jet2_ystar;
  hn->Fill( &x[0], jetWeight );

  x[5] = 1;
  hn->Fill( &x[0], jetWeight * deltaPhi );
 
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

void DiJetAnalysis::NormalizeDeltaPhi( TH1* hIn, TH1* hNorm,
				       double xBinCenter, bool comb ){

  if( !hIn->GetEntries() ){ return; }

  // this is true for all, do
  // before comb subtraction
  hIn->Scale( 1., "width" );
  
  if( comb ){ anaTool->SubtractCombinatoric( hIn ); }

  if( !hNorm ){
    hIn->Scale( 1./hIn->Integral() ) ;
  } else {
    int xBin = hNorm->FindBin( xBinCenter );
    std::cout << "  " << hNorm->GetName() << std::endl;
    std::cout << "  " << hIn->Integral() << " " << " " << xBin << " "
	      << hNorm->GetBinContent( xBin ) << " "
	      << hIn->Integral()/hNorm->GetBinContent( xBin ) << std::endl;
    hIn->Scale( 1./hNorm->GetBinContent( xBin ) );
  }
}

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
    // double eC = hC->GetBinError( cBin );
    
    // correction factor;
    double newDphi = vM * vC;
    // error on correction factor    
    double newDphiError = eM * vC; 

    /*
    double newDphiError =  newDphi * 
      std::sqrt( std::pow( eM / vM, 2) +
		 std::pow( eC / vC, 2) ) ; 
    */
    
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
      std::cout << "Adding  to systematics: " << m_fNamePhysUF << std::endl;
      mFinUC[ uc ] =  TFile::Open( m_fNamePhysUF.c_str() ) ;
    }
  }
  
  return fInDef;
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

void DiJetAnalysis::MakeEtaPhiPtMap( std::vector< TH2* >& vSampleMaps,
				     const std::vector< std::string>& vLabels,
				     const std::string& name ){

  uint nSamples = vSampleMaps.size();
  
  TCanvas c_map("c_map","c_map",800,600);

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

  std::string axisLabel, axisLabelTex;
  GetSpectraLabels( axisLabel, axisLabelTex, name );

  bool isEta = name.find( m_sEta ) != std::string::npos ? true : false;
    
  if( !vSampleSpect.size() ){ return; }
  
  std::string yAxisTitle = "dN/d#it{p}_{T} [GeV]";
  // std::string yAxisTitle = "N_{jets}";

  // use this as reference because
  // it should be in every file
  TH2*  hRef = vSampleSpect[0];
  int nXbins = hRef->GetNbinsX();

  uint nSamples = vSampleSpect.size();
  
  std::vector< std::vector< TH1* > > vSpect;
  vSpect.resize( nSamples );

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
      
      vSpect[iG].push_back( hSpectCounts );
      vSpect[iG].push_back( hSpect );
      
      // in case var binning, scale by width
      hSpect->Scale( 1., "width" );
      
      hSpect->Write();
      hSpectCounts->Write();
      
      // get min max from the final histograms
      if( label.compare( m_allName )  ){ continue; }
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
  else          { lX0 = 0.35; lY0 = 0.67; lX1 = 0.65; lY1 = 0.88; }
  
  for( int iX = 0; iX < nXbins; iX++ ){
    int    xBin = iX + 1;
    double xMin, xMax;
    anaTool->GetBinRange
      ( hRef->GetXaxis(), xBin, xBin, xMin, xMax );
    double xCenter = hRef->GetXaxis()->GetBinCenter ( xBin );

    // for pPb, dont draw at anything above -3.2
    if( isEta && m_is_pPb && xCenter > -constants::FETAMIN ){ continue; }
     
    std::string cName  = anaTool->GetName ( xMin, xMax, axisLabel    );
    std::string cLabel = anaTool->GetLabel( xMin, xMax, axisLabelTex );
    
    TCanvas c( cLabel.c_str(), cLabel.c_str(), 800, 600 );
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

      TH1* h = vSpect[iG][iX];
           
      if( !h->GetEntries() && !isMb && !isAll ){ continue; }
      else if( !h->GetEntries() && !isMb && !isAll ){ continue; }
      
      if( iG == nSamples ){ style = 0; }
      styleTool->SetHStyle( h, style++ );
      h->Draw("epsame");
      h->SetMinimum( m_ptSpectYaxisMin );
      h->SetMaximum( max );
      leg.AddEntry( h, label.c_str() );
    } // end loop over iG
    
    leg.Draw("same");

    DrawAtlasRight();    
    drawTool->DrawRightLatex( 0.4, 0.2, cLabel );

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
    styleTool->SetLegendStyle( &leg, 0.65 );

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
    for( int iX = 0; iX < nXbins; iX++ )
      { delete vSpect[iG][iX]; }
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
	      m_nVarPtBinsRespMat, 0, 1 );
  hUnfoldedCountsAll->GetXaxis()->Set( m_nVarYstarBins    , &( m_varYstarBinning[0]     ) );
  hUnfoldedCountsAll->GetYaxis()->Set( m_nVarPtBinsRespMat, &( m_varPtBinningRespMat[0] ) );

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
    
    TH1* hTruthCounts    = static_cast< TH1D* >
      ( fInMC->Get( Form("h_%s_%s_%s_%s", truthName.c_str(), m_sCounts.c_str(),
			 m_allName.c_str(), hTag.c_str() ) ) );

    TH1* hTruth    = static_cast< TH1D* >
      ( fInMC->Get( Form("h_%s_%s_%s", truthName.c_str(),
			 m_allName.c_str(), hTag.c_str() ) ) );

    
    styleTool->SetHStyle( hMeasuredCounts, 1 );
    styleTool->SetHStyle( hTruthCounts   , 2 );
    styleTool->SetHStyle( hMeasured, 1 );
    styleTool->SetHStyle( hTruth   , 2 );

    vSpect.push_back( hMeasuredCounts );
    vSpect.push_back( hTruthCounts    );
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
      double integralData = hUnfolded->Integral( 2, m_nVarPtBinsRespMat - 1 );
      double integralMC   = hReco    ->Integral( 2, m_nVarPtBinsRespMat - 1 );

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

    double xMin = hR->GetXaxis()->GetXmin();
    double xMax = hR->GetXaxis()->GetXmax();
	  
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

  // These are "quality" check histograms
  TH1D* h_mult     = new TH1D( Form("h_mult_%s", dPhiName.c_str() ),
			 ";Count;Multiplicity", 100,0,5000);
  TH2D* h_tauSigma = new TH2D( Form("h_tauSigma_%s", dPhiName.c_str() ),
			 ";#tau;#sigma"   , 60, 0., 0.6, 60, 0., 0.3 );
  TH2D* h_tauAmp   = new TH2D( Form("h_tauAmp_%s"  , dPhiName.c_str() ),
			 ";#tau;Amplitude", 60, 0., 0.6, 60, 0., 0.3 );


  std::vector< std::string > v_listBadHist;
  std::vector< std::string > v_listBadFits;
  
  TH1D* h_dPhiChi2 = new TH1D( Form( "h_dPhiChi2_%s", dPhiName.c_str() ),
			       ";#Chi^{2}/NDF;Count", 20, 0, 5 );
  TH1D* h_dPhiProb = new TH1D( Form( "h_dPhiProb_%s", dPhiName.c_str() ),
			       ";Probability;Count", 20, 0, 1 );

  h_dPhiChi2->Sumw2();
  h_dPhiProb->Sumw2();
  
  styleTool->SetHStyle( h_dPhiChi2, 0 );
  styleTool->SetHStyle( h_dPhiProb, 0 );
  
  std::vector< TH1*    > vHdPhi;
  std::vector< TGraph* > vGdPhi;
  std::vector< TH1* > vDphi;
  std::vector< TH1* > vSpect;
  std::vector< TF1* > vFits;
  std::vector< THnSparse* > vHnDphi;
  
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
    
    int    dPhiAxisI = hn->GetNdimensions() - 2;
    int   entryAxisI = hn->GetNdimensions() - 1;
    TAxis* entryAxis = hn->GetAxis( entryAxisI );
    
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
	  TH1* hDphiWidths = hn->Projection( fAxisI );
	  hDphiWidths->Reset();
	  
	  hDphiWidths->SetName
	    ( Form( "h_%s_%s_%s", dPhiName.c_str(), label.c_str(), hTagW.c_str() ) );
	  hDphiWidths->SetYTitle( "RMS (#pi - #Delta#phi)" );
	  hDphiWidths->SetMarkerSize( hDphiWidths->GetMarkerSize() * 1.5 );
	  styleTool->SetHStyle( hDphiWidths, style );
	  vDphiWidthsTemp.push_back( hDphiWidths );

	  TH1* hDphiWidthsStat = hn->Projection( fAxisI );
	  hDphiWidthsStat->Reset();
	  
	  hDphiWidthsStat->SetName
	    ( Form( "h_%s_stat_%s_%s", dPhiName.c_str(), label.c_str(), hTagW.c_str() ) );
	  hDphiWidthsStat->SetYTitle( "RMS (#pi - #Delta#phi)" );
	  hDphiWidthsStat->SetMarkerSize( hDphiWidthsStat->GetMarkerSize() * 2 );
	  styleTool->SetHStyle( hDphiWidthsStat, 5 + style );

	  TH1* hDphiWidthsStat2 = hn->Projection( fAxisI );
	  hDphiWidthsStat2->Reset();
	  
	  hDphiWidthsStat2->SetName
	    ( Form( "h_%s_stat_%s_%s", dPhiName.c_str(), label.c_str(), hTagW.c_str() ) );
	  hDphiWidthsStat2->SetYTitle( "RMS (#pi - #Delta#phi)" );
	  hDphiWidthsStat2->SetMarkerSize( hDphiWidthsStat2->GetMarkerSize() * 2 );
	  styleTool->SetHStyle( hDphiWidthsStat2, 5 + style );
	  hDphiWidthsStat2->SetMarkerColor( kGreen );
	  hDphiWidthsStat2->SetLineColor  ( kGreen );

	  style++;
	  
	  legW.AddEntry
	    ( hDphiWidths,
	      anaTool->GetLabel( axis2Low, axis2Up, m_dPP->GetAxisLabel(2) ) .c_str() );

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
	    std::string gDphiName =
	      Form( "g_%s_%s_%s", dPhiName.c_str(), label.c_str(), hTag.c_str() );
	    std::string hDphiCountsName =
	      Form( "h_%s_%s_%s_%s",
		    dPhiName.c_str(), m_sCounts.c_str(), label.c_str(), hTag.c_str() );
	    std::string hDphiCountsSumName =
	      Form( "h_%s_%s_sum_%s_%s",
		    dPhiName.c_str(), m_sCounts.c_str(), label.c_str(), hTag.c_str() );
	    
	    entryAxis->SetRange( 1, 1 );
	    
	    TH1* hDphiCounts = hn->Projection( dPhiAxisI );
	    hDphiCounts->SetName( hDphiCountsName.c_str() );
	    styleTool->SetHStyle( hDphiCounts, 0 );
	    vHdPhi.push_back( hDphiCounts );

	    entryAxis->SetRange( 2, 2 );
	    
	    TH1* hDphiCountsSum = hn->Projection( dPhiAxisI );
	    hDphiCountsSum->SetName( hDphiCountsSumName.c_str() );
	    styleTool->SetHStyle( hDphiCountsSum, 0 );
	    vHdPhi.push_back( hDphiCountsSum );

	    TGraph* gDphi = anaTool->Barycenters( hDphiCounts, hDphiCountsSum );
	    vGdPhi.push_back( gDphi );
	    gDphi->SetName( gDphiName.c_str() );
	    styleTool->SetHStyle( gDphi, 0 );

	    h_mult->Fill( hDphiCounts->GetEntries() );
	    
	    // save the normalized counts histogram.
	    TH1* hDphi = static_cast< TH1D* >
	      ( hDphiCounts->Clone( hDphiName.c_str() ) );
	    styleTool->SetHStyle( hDphi, 0 );
	    vHdPhi.push_back( hDphi );
	      
	    // Normalize with combinatoric subtraction after
	    // scaling by width (the last parameter true)
	    NormalizeDeltaPhi
	      ( hDphi, hSpectCounts, 0.5 * ( axis1Up + axis1Low ), true );
	    
	    anaTool->MatchGraphHistoX( gDphi, hDphi );

	    // for canvas, since its tgraphasymmerrors.
	    double x0 = m_dPhiZoomLow;
	    double x1 = m_dPhiZoomHigh;
	    double y0 = m_dPhiLogMin;
	    double y1 = m_dPhiLogMax;
	    
	    y1 = anaTool->GetGraphMax( gDphi, m_dPhiFittingMin, m_dPhiFittingMax ) * 1.1;
	    
	    TCanvas c( hDphi->GetName(), hDphi->GetName(), 800, 600 );
	    styleTool->SetCStyleGraph
	      ( c, x0, y0, x1, y1, ";#Delta#phi;1/N_{jet_{1}} dN_{pair}/d|#Delta#phi|" );
	    
	    gDphi->Draw( "p" );
	    // hDphi->Draw();
	    // hDphi->SetYTitle("1/N_{jet_{1}} dN_{pair}/d|#Delta#phi|");
	    gDphi->SetTitle("");
	    
	    // now fit
	    TF1* fit = anaTool->FitDphi( gDphi, m_dPhiFittingMin, m_dPhiFittingMax );
	    styleTool->SetHStyle( fit, 0 );
	    vFits.push_back( fit );

	    double   amp = fit->GetParameter(0);
	    double   tau = fit->GetParameter(1);
	    double sigma = fit->GetParameter(2);

	    double   tauError = fit->GetParError(1);
	    double sigmaError = fit->GetParError(2);
	    
	    h_tauAmp  ->Fill( tau, amp   );
	    h_tauSigma->Fill( tau, sigma );
	  
	    fit->Draw("same");

	    double chi2NDF = fit->GetChisquare()/fit->GetNDF();
	    double prob    = fit->GetProb();

	    h_dPhiChi2->Fill( chi2NDF );
	    h_dPhiProb->Fill( prob    );

	    if( prob < 0.05 ){
	      v_listBadHist.push_back( hDphi->GetName() );
	      v_listBadFits.push_back( fit  ->GetName() );
	    }
	    
	    hDphi->GetXaxis()->SetRange( m_dPhiZoomLowBin, m_dPhiZoomHighBin );

	    drawTool->DrawLeftLatex( 0.47, 0.62, Form( "#Chi^{2}/NDF=%4.2f", chi2NDF ) );
	    drawTool->DrawLeftLatex( 0.47, 0.54, Form( "Prob=%4.2f"        , prob    ) );

	    DrawTopLeftLabels
	      ( m_dPP, axis0Low, axis0Up, axis1Low, axis1Up,
		axis2Low, axis2Up, axis3Low, axis3Up );
	    
	    DrawAtlasRight();
	    
	    SaveAsROOT( c, hDphi->GetName() );
	    fit->Write();

	    gDphi->Write();
	    hDphi->Write();
	    hDphiCounts->Write();
	    hDphiCountsSum->Write();
	    
	    if( fit->GetParameter(1) < 0 )
	      { continue; }

	    double width      = std::sqrt( tau * tau + sigma * sigma ); 
	    double widthError = std::sqrt( std::pow( 2 * tau   * tauError  , 2 ) +
					   std::pow( 2 * sigma * sigmaError, 2 ) );

	    std::cout << hDphi->GetName() << " " << width << " " << widthError << std::endl;
	    
	    // Now, put results on histogram
	    std::pair< double, double > rmsAndError =
	      anaTool->GetRMS( hDphi, m_dPhiFittingMin, m_dPhiFittingMax, constants::PI );

	    double widthStat      = rmsAndError.first;
	    double widthStatError = rmsAndError.second;

	    // Now, put results on histogram
	    std::pair< double, double > rmsAndError2 =
	      anaTool->GetRMS( hDphi, 1, m_dPhiFittingMax, constants::PI );

	    double widthStat2      = rmsAndError2.first;
	    double widthStat2Error = rmsAndError2.second;

	    std::cout << "Tau   = " << tau   << " TauError   = " << tauError   << std::endl;
	    std::cout << "Sigma = " << sigma << " SigmaError = " << sigmaError << std::endl;
	    std::cout << "Width = " << width << " WidthError = " << widthError << std::endl;
	    
	    hDphiWidths->SetBinContent( axis3Bin, width      );
	    hDphiWidths->SetBinError  ( axis3Bin, widthError );	    

	    hDphiWidthsStat->SetBinContent( axis3Bin, widthStat      );
	    hDphiWidthsStat->SetBinError  ( axis3Bin, widthStatError );

	    hDphiWidthsStat2->SetBinContent( axis3Bin, widthStat2      );
	    hDphiWidthsStat2->SetBinError  ( axis3Bin, widthStat2Error );
	  } // end loop over axis3

	  TCanvas cWidthsAll( "cWidthsAll", "cWidthsAll",800,600);
	  TLegend legWAll( 0.50, 0.2, 0.89, 0.38 );
	  styleTool->SetLegendStyle( &legWAll );

	  hDphiWidths->SetMinimum( m_dPhiWidthMin );
	  hDphiWidths->SetMaximum( m_dPhiWidthMax );
	
	  hDphiWidths     ->Draw("epsame");
	  hDphiWidthsStat ->Draw("epsame");
	  hDphiWidthsStat2->Draw("epsame");

	  legWAll.AddEntry( hDphiWidths     , "RMS Fit" );
	  legWAll.AddEntry( hDphiWidthsStat , "Statistics 2#pi/3<#Delta#phi<#pi" );
	  legWAll.AddEntry( hDphiWidthsStat2, "Statistics 1<#Delta#phi<#pi" );

	  legWAll.Draw();

	  DrawTopLeftLabels
	    ( m_dPP, axis0Low, axis0Up, axis1Low, axis1Up,
	      axis2Low, axis2Up, 0, 0 );

	  DrawAtlasRight();
	  
	  SaveAsAll( cWidthsAll, Form("h_%s_%s_All_%s", dPhiName.c_str(),
				      label.c_str(), hTagW.c_str() ) );
	  
	} // end loop over axis2
	
	cWidths.cd();
	for( auto& h : vDphiWidthsTemp ){
	  h->SetMinimum( m_dPhiWidthMin );
	  h->SetMaximum( m_dPhiWidthMax );
	  h->SetTitle("");
	  h->Draw("epsame");
	  h->Write();
	}
	legW.Draw("same");

	DrawTopLeftLabels
	  ( m_dPP, axis0Low, axis0Up, axis1Low, axis1Up, 0, 0, 0, 0 );

	DrawAtlasRight();
	
        SaveAsAll( cWidths, Form("h_%s_%s_%s", dPhiName.c_str(),
				 label.c_str(), hTagCW.c_str() ) );
      } // end loop over axis1     
    } // end loop over axis0
  } // end loop over iG

  TCanvas c( "c", "c", 1200, 600 );
  c.Divide( 2, 1 );

  c.cd(1);
  h_dPhiChi2->Draw("ep");
  DrawAtlasRight();

  c.cd(2);
  h_dPhiProb->Draw("ep");
  DrawAtlasRight();

  SaveAsAll( c, Form( "h_chi2_prob_%s", dPhiName.c_str() ) );

  delete h_dPhiChi2;
  delete h_dPhiProb;

  std::cout << " ------------ Bad Fits - " << v_listBadFits.size() << std::endl;
  for( uint i = 0; i < v_listBadFits.size(); i++ ){
    std::cout << v_listBadHist[i] << std::endl;
    std::cout << v_listBadFits[i] << std::endl;
    std::cout << "    --" << std::endl;
  }
  
  h_mult    ->Write();
  h_tauSigma->Write();
  h_tauAmp  ->Write();
  delete h_mult;
  delete h_tauSigma;
  delete h_tauAmp;
  
  for( auto& f : vFits  ){ delete f; }
  for( auto& h : vHdPhi ){ delete h; }
  for( auto& h : vGdPhi ){ delete h; }
  for( auto& h : vSpect ){ delete h; }
  
  // Save the THnSparse to separate file.
  // This is used for unfolding later.
  // Gets put into data/ dirrectory.
  std::string fOutName =
    m_dirOut + "/" + dPhiName + "_" + m_labelOut + ".root";
  
  TFile* fOut = new TFile( fOutName.c_str(), "RECREATE" );

  for( auto& hn : vHnDphi ){ hn->Write(); delete hn; }

  fOut->Close();
}

THnSparse* DiJetAnalysis::UnfoldDeltaPhi( TFile* fInData, TFile* fInMC,
					  const std::string& hnUnfoldedName,
					  TFile* fInMCSpect,
					  const std::string& spectName ){

  std::vector< TH1* > vDphi;
  std::vector< TH1* > vSpect;
  std::vector< TF1* > vDphiFits;
  
  std::vector< TH1* > vCfactors;
  std::vector< TH1* > vRatios;
  
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
		    m_nDphiDim, &m_vNdPhiBins[0],
		    &m_vDphiMin[0], &m_vDphiMax[0] );

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
  
  // ---- loop over ystars ----
  for( int axis0Bin = 1; axis0Bin <= nAxis0Bins; axis0Bin++ ){
    // check we are in correct ystar and pt bins
    if( !m_dPP->CorrectPhaseSpace
	( std::vector<int>{ axis0Bin, 0, 0, 0 } ) )
      { continue; }

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
	    Form( "%s_%s_%s_%s",
		  anaTool->GetName( axis0Low, axis0Up, m_dPP->GetAxisName(0) ).c_str(),
		  anaTool->GetName( axis1Low, axis1Up, m_dPP->GetAxisName(1) ).c_str(),
		  anaTool->GetName( axis2Low, axis2Up, m_dPP->GetAxisName(2) ).c_str(),
		  anaTool->GetName( axis3Low, axis3Up, m_dPP->GetAxisName(3) ).c_str() ); 

	  // Get measured and truth distributions
	  // both normalized and not. The counts are the
	  // un-normalized histograms.
	  // I.e. hMeasured and hTruth
	  // Since we already did normalization and fitting
	  // we just retrieve them.
	  TH1* hMeasuredCounts = static_cast<TH1D*>
	    ( fInData->Get
	      ( Form( "h_%s_%s_%s_%s", measuredName.c_str(), m_sCounts.c_str(),
		      m_allName.c_str(), hTag.c_str())));

	  TH1* hTruthCounts    = static_cast<TH1D*>
	    ( fInMC  ->Get
	      ( Form( "h_%s_%s_%s_%s", truthName.c_str(), m_sCounts.c_str(),
		      m_allName.c_str(), hTag.c_str())));
	  TH1* hMeasured       = static_cast<TH1D*>
	    ( fInData->Get( Form( "h_%s_%s_%s", measuredName.c_str(),
				  m_allName.c_str(), hTag.c_str())));
	  TH1* hTruth          = static_cast<TH1D*>
	    ( fInMC  ->Get( Form( "h_%s_%s_%s", truthName.c_str(),
				  m_allName.c_str(), hTag.c_str())));
	  styleTool->SetHStyle( hMeasured, 1 );
	  styleTool->SetHStyle( hTruth   , 2 );
	  
	  vDphi.push_back( hMeasuredCounts );
	  vDphi.push_back( hTruthCounts    );
	  vDphi.push_back( hMeasured       );
	  vDphi.push_back( hTruth          );

	  // Get correction factors. IMPORTANT
	  TH1* hCfactors = static_cast<TH1D*>
	    ( fInMC->Get( Form( "h_%s_%s_%s", m_dPhiCfactorsName.c_str(),
				m_allName.c_str(), hTag.c_str())));
	  vCfactors.push_back( hCfactors );
	  
	  // ----------- Unfold -----------
	  // Unfold using bin-by-bin and the resposne factors.
	  // Do this on counts, then normalize and subtract comb.
	  TH1* hUnfoldedCounts = BinByBinUnfolding
	    ( hMeasuredCounts,  hCfactors );
	  hUnfoldedCounts->SetName
	    ( Form( "h_%s_%s_%s_%s", m_dPhiUnfoldedName.c_str(),
		    m_sCounts.c_str(), m_allName.c_str(), hTag.c_str()));
	  styleTool->SetHStyle( hUnfoldedCounts, 0 );
	  vDphi.push_back( hUnfoldedCounts );

	  // fill unfolded THnSparse result
	  std::vector< int > x  = m_dPP->GetMappedBins
	    ( std::vector<int> { axis0Bin, axis1Bin, axis2Bin, axis3Bin } );
	  x.push_back(0); // dPhiBin;
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
	  vDphi.push_back( hUnfolded );
	  
	  // Normalize
	  NormalizeDeltaPhi( hUnfolded, hSpectCounts, 0.5 * ( axis1Up + axis1Low ) );
	  TF1* fitUnfolded = anaTool->FitDphi( hUnfolded, m_dPhiFittingMin, m_dPhiFittingMax );
	  
	  // -------- Unfold Done ---------
	  	  	  
	  hMeasured->GetXaxis()->SetRange( m_dPhiZoomLowBin, m_dPhiZoomHighBin );
	  hUnfolded->GetXaxis()->SetRange( m_dPhiZoomLowBin, m_dPhiZoomHighBin );
	  hTruth   ->GetXaxis()->SetRange( m_dPhiZoomLowBin, m_dPhiZoomHighBin );
	  
	  // Now Draw everything.
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
	  
	  hMeasured->Draw("ep same");	
	  hTruth   ->Draw("histo same");
	  hUnfolded->Draw("ep same");
	  
       	  TF1* fitMeasured = static_cast<TF1*>
	    ( fInData->Get( Form( "f_g_%s_%s_%s", measuredName.c_str(),
				  m_allName.c_str(), hTag.c_str())));
	  TF1* fitTruth    = static_cast<TF1*>
	    ( fInMC  ->Get( Form( "f_g_%s_%s_%s", truthName.c_str(),
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
	  
	  double maximum = -1;
	  
	  maximum = hUnfolded->GetMaximum() > hMeasured->GetMaximum() ?
	    hUnfolded->GetMaximum() : hMeasured->GetMaximum();
	  maximum = maximum > hTruth->GetMaximum() ?
	    maximum : hTruth->GetMaximum();

	  maximum = anaTool->GetLogMaximum( maximum );

	  hMeasured->SetMinimum( m_dPhiLogMin );
	  hMeasured->SetMaximum( m_dPhiLogMax );

	  leg.AddEntry( hTruth   , "Truth", "lf" );
	  leg.AddEntry( hMeasured, typeMeasured.c_str() );
	  leg.AddEntry( hUnfolded, "Unfolded" );
	  
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

	  hRfit->Draw( "ep same" );
	  
	  legR.AddEntry( hCfactors,       "CF" );
	  // legR.AddEntry( hR       , "UF/Truth" );
	  legR.AddEntry( hRfit    ,   "Fit/UF" );
	  
	  legR.Draw();

	  styleTool->HideAxis( hMeasured, "x" );
	  styleTool->HideAxis( hUnfolded, "x" );
	  styleTool->HideAxis( hTruth   , "x" );
	  
	  double xMin = m_dPhiZoomLow;
	  double xMax = m_dPhiZoomHigh; 
	  
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

	  /*
	  TLine lineP50( xMin, 1.50, xMax, 1.50 );
	  lineP50.SetLineStyle( 2  );
	  lineP50.SetLineColor( 12 );
	  lineP50.SetLineWidth( 1  );
	  lineP50.Draw();
	  
	  TLine lineN50( xMin, 0.50, xMax, 0.50 );
	  lineN50.SetLineStyle( 2  );
	  lineN50.SetLineColor( 12 );
	  lineN50.SetLineWidth( 1  );
	  lineN50.Draw();
	  */
	  
	  // save this for debugging
	  // hUnfoldedCounts->Write();
	  
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
	( 0.13, 0.86, CT::AnalysisTools::GetLabel
	  ( axis0Low, axis0Up, m_dPP->GetDefaultAxisLabel(0) ) );
      drawTool->DrawLeftLatex
	( 0.13, 0.79, CT::AnalysisTools::GetLabel
	  ( axis3Low, axis3Up, m_dPP->GetDefaultAxisLabel(1) ) );  

      DrawAtlasRight();

      SaveAsAll( c1, hCfactorsPtMat->GetName() );
      hCfactorsPtMat->Write();
    }
  }
  
  for( auto vh : vHcFactorsPtMat )
    { for( auto h : vh ){ delete h; } }

  for( auto f : vDphiFits    ){ delete f;  }
  for( auto h : vDphi        ){ delete h ; }
  for( auto h : vSpect       ){ delete h ; }

  for( auto r : vCfactors    ){ delete r;  }

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

  double xMin = m_varPtBinningRespMat.front();
  double xMax = m_varPtBinningRespMat.back();
  
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

    std::string hName_a = "h_" + name_a + "_" + hTag;
    std::string hName_b = "h_" + name_b + "_" + hTag;
    
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
      hSpect_a->Integral( 2, m_nVarPtBinsRespMat - 1 ) /
      hSpect_b->Integral( 2, m_nVarPtBinsRespMat - 1 );

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

      // tag for widths canvas
      std::string hTagCW =
	Form ("%s_%s",
	      anaTool->GetName( axis0Low, axis0Up, m_dPP->GetAxisName(0) ).c_str(),
	      anaTool->GetName( axis1Low, axis1Up, m_dPP->GetAxisName(1) ).c_str() );
      
      TLegend legW( 0.47, 0.05, 0.88, 0.45 );
      styleTool->SetLegendStyle( &legW, 0.75 );
      // legW.SetNColumns(2);
      
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

	// switch back to cW canvas
	// because a new one was creted in axis3 loop
	cW.cd();
	pad2W.cd();

	// Make the ratio histogram.
	TH1* hW_R = static_cast< TH1D* >
	  ( hW_a->Clone( Form( "h_%s_%s_%s_%s", m_dPhiName.c_str(), m_sRatio.c_str(),
			       m_allName.c_str(), hTagW.c_str())));
	hW_R->Divide( hW_b );
	hW_R->SetYTitle( ratio.c_str() );
	hW_R->SetTitleOffset( 2.3, "x" );
	styleTool->SetHStyleRatio( hW_R, style );
	hW_R->SetMarkerSize( hW_R->GetMarkerSize() * 1.5 );
	
	vR.push_back( hW_R );
	hW_R->Draw("ep same");

	style++;
	
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
	  
	  TH1* hR = static_cast< TH1D* >
	    ( h_a->Clone( Form( "h_%s_%s_%s_%s", m_dPhiName.c_str(), m_sRatio.c_str(),
				m_allName.c_str(), hTagW.c_str())));
	  styleTool->SetHStyleRatio( hR );
	  hR->Divide( h_b );
	  hR->SetYTitle( ratio.c_str() );
	  hR->SetTitleOffset( 2.3, "x" );
	  vR.push_back( hR );

	  pad1.cd();
	  pad1.SetLogy();
	  
	  styleTool->HideAxis( h_a, "x" );
	  styleTool->HideAxis( h_b, "x" );
	  
	  h_a->Draw("epsame");
	  h_b->Draw("epsame");
	  
	  TF1* f_a = static_cast<TF1*>( fIn_a->Get( Form("f_g_%s", hName_a.c_str())));
	  TF1* f_b = static_cast<TF1*>( fIn_b->Get( Form("f_g_%s", hName_b.c_str())));
	  styleTool->SetHStyle( f_a, 0 );
	  styleTool->SetHStyle( f_b, 1 );
	  f_a->SetLineColor( h_a->GetLineColor() );
	  f_b->SetLineColor( h_b->GetLineColor() );
	  vF.push_back( f_a ); vF.push_back( f_b );

	  // TLegend leg( 0.45, 0.10, 0.66, 0.28 );
	  TLegend leg( 0.7, 0.06, 0.9, 0.23 );
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

	  DrawAtlasRightBoth();

	  pad2.cd();
	  hR->Draw("ep same");

	  dPhiLine.Draw();
	  dPhiLineP25.Draw(); dPhiLineN25.Draw();
	  // dPhiLineP50.Draw(); dPhiLineN50.Draw();

	  SaveAsAll( c, Form("h_%s_%s", m_dPhiName.c_str(), hTag.c_str() ), true );
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

      SaveAsPdfPng( cW, Form("h_%s_%s", m_dPhiName.c_str(), hTagCW.c_str() ), true );
      SaveAsROOT( cW, Form("h_%s_%s", m_dPhiName.c_str(), hTagCW.c_str() ) );

    } // end loop over ystar2
  } // end loop over ystar2

  for( auto & f  : vF  ){ delete f ; }
  for( auto & r  : vR  ){ delete r ; }
  for( auto & h  : vH  ){ delete h ; }
  for( auto & hW : vHw ){ delete hW; }
}

void DiJetAnalysis::MakeFinalPlotsTogether(){}

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
    ( 0.18, 0.86, CT::AnalysisTools::GetLabel
      ( axis0Low, axis0Up, dPP->GetAxisLabel(0) ), scale );
  // for now, modify later to be more dynamic
  // to the variables present 
  if( !( axis1Low || axis1Up ) ){ return ; }
  drawTool->DrawLeftLatex
    ( 0.18, 0.76, CT::AnalysisTools::GetLabel
      ( axis1Low, axis1Up, dPP->GetAxisLabel(1) ), scale );  
  if( !( axis2Low || axis2Up ) ){ return ; }
  drawTool->DrawLeftLatex
    ( 0.18, 0.67, CT::AnalysisTools::GetLabel
      ( axis2Low, axis2Up, dPP->GetAxisLabel(2) ), scale );
  if( !( axis3Low || axis3Up ) ){ return ; }
  drawTool->DrawLeftLatex
    ( 0.18, 0.57, CT::AnalysisTools::GetLabel
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
