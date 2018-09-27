#include <TROOT.h>
#include <TEnv.h>
#include <TGraphAsymmErrors.h>
#include <TLorentzVector.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TGaxis.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TF1.h>
#include <THnSparse.h>

#include <iostream>
#include <sstream>
#include <cmath>
#include <cfloat>
#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <boost/assign.hpp>

#include "DiJetAnalysisMC.h"
#include "DeltaPhiProj.h"
#include "UncertaintyProvider.h"

static const bool compToPythia = true;

DiJetAnalysisMC::DiJetAnalysisMC()
  : DiJetAnalysisMC( false, 0, 0 ) {}

DiJetAnalysisMC::DiJetAnalysisMC( bool is_pPb )
  : DiJetAnalysisMC( is_pPb, 0 , 0 ) {}

DiJetAnalysisMC::DiJetAnalysisMC( bool is_pPb, int mcType )
  : DiJetAnalysisMC( is_pPb, mcType, 0 ){}

DiJetAnalysisMC::DiJetAnalysisMC( bool is_pPb, int mcType, int uncertComp )
  : DiJetAnalysis( is_pPb, false, mcType, uncertComp ){
  
  //========== Set Histogram Binning =============
  // ------ truth binning --------
  m_ptTruthWidth  = 5;
  m_ptTruthMin    = 25;
  m_ptTruthMax    = 90;
  m_nPtTruthBins  = (m_ptTruthMax - m_ptTruthMin) / m_ptTruthWidth;

  // ---- JES/PRes/Etc ----- 
  m_scaleResUseEta = true;
  
  // --- variable eta/ystar binning ---
  m_nRPtRecoTruthBins = 100;
  m_rPtRecoTruthMin   = 0;
  m_rPtRecoTruthMax   = 2;

  m_nDAngleRecoTruthBins = 100;
  m_dAngleRecoTruthMin   = -0.5;
  m_dAngleRecoTruthMax   = 0.5;

  // --- Dphi Resp Matrix binning ---  

  boost::assign::push_back( m_vNdPhiRespMatBins )
    ( m_nVarYstarBins  )( m_nVarYstarBins  )
    ( m_nVarPtBinsUfOf )( m_nVarPtBinsUfOf )
    ( m_nVarPtBinsUfOf )( m_nVarPtBinsUfOf )
    ( m_nVarDphiBins   )( m_nVarDphiBins   );

  boost::assign::push_back( m_vNdPhiRespMatRebBins )
    ( m_nVarYstarBins  )( m_nVarYstarBins  )
    ( m_nVarPtBinsUfOf )( m_nVarPtBinsUfOf )
    ( m_nVarPtBinsUfOf )( m_nVarPtBinsUfOf )
    ( m_nVarDphiRebinnedBins  )( m_nVarDphiRebinnedBins  );
  
  boost::assign::push_back( m_vDphiRespMatMin  )
    ( 0 )( 0 )( 0 )( 0 )( 0 )( 0 )( 0 )( 0 );

  boost::assign::push_back( m_vDphiRespMatMax  )
    ( 1 )( 1 )( 1 )( 1 )( 1 )( 1 )( 1 )( 1 );

  m_nDphiRespMatDim = m_vNdPhiRespMatBins.size();
  
  //==================== Cuts ====================    
  m_dRmax = 0.2;

  //=============== Histo Names ==================
  m_ystarSpectFineTruthUPName =
    m_ystarSpectFineName + "_" + m_unpairedName + "_" + m_truthName;

  
  m_dPhiRespMatName        = m_dPhiName + "_" + m_respMatName;
  m_dPhiRespMatRebName     = m_dPhiName + "_" + m_respMatName + "_" + m_sReb;
  m_ptRespMatName          = m_s_pt     + "_" + m_respMatName;

  m_dPhiRecoUnfoldedName   = m_dPhiRecoName + "_" + m_unfoldedName;
  
  //============== tool ===============
  m_uncertaintyProvider = new UncertaintyProvider( m_uncertComp, m_is_pPb );

}

DiJetAnalysisMC::~DiJetAnalysisMC(){
  delete m_uncertaintyProvider;
}

void DiJetAnalysisMC::Initialize()
{
  // Initalize things common to everything
  DiJetAnalysis::Initialize();
  
  std::string mcMenu = GetMCMenu();

  // Get Info On jzn slices
  std::vector< int > vJznUsedD = anaTool->vectoriseI
    ( GetConfig()->GetValue( Form("jznUsed.%s", mcMenu.c_str()),"" )," ");
  // double -> int 
  for( auto d : vJznUsedD ){ m_vJznUsed.push_back(d); }
  m_vJznLabels  =  anaTool->vectorise
    ( GetConfig()->GetValue( Form("jznLabel.%s"  ,mcMenu.c_str()),"" )," ");
  m_vJznFnameIn =  anaTool->vectorise
    ( GetConfig()->GetValue( Form("jznFnameIn.%s",mcMenu.c_str()),"" )," ");
  m_vJznSigma   = anaTool->vectoriseD
    ( GetConfig()->GetValue( Form("jznSigma.%s"  ,mcMenu.c_str()),"" )," ");
  m_vJznEff     = anaTool->vectoriseD
    ( GetConfig()->GetValue( Form("jznEff.%s"    ,mcMenu.c_str()),"" )," ");
  m_vJznSumOverlayWeights = anaTool->vectoriseD
    ( GetConfig()->GetValue( Form("jznSumOverlayWeights.%s",mcMenu.c_str()),"" )," ");
  m_vJznPtThreshold       = anaTool->vectoriseD
    ( GetConfig()->GetValue( Form("jznPtThreshold.%s"      ,mcMenu.c_str()),"" )," ");

  std::string weightFile =
    GetConfig()->GetValue( Form("weightFile.%s",mcMenu.c_str()), "");
  if( !weightFile.empty() ){
    // Get Associated weights for powheg.
    // These need to be used at filling time (per event basis).
    TFile* file =
      TFile::Open( weightFile.c_str() );
    if( !m_hPowhegWeights ) // because its a static member 
      { m_hPowhegWeights = static_cast< TH3D* >( file->Get("h3_pT_y_phi_rw") ); }
    m_hPowhegWeights->SetDirectory(0);
    file->Close();
  }

  m_nJzn = m_vJznUsed.size();  

  //========== Cuts and jzn slices =============
  
  // calculate sum of sigma and eff
  m_sumSigmaEff = 0;

  for( uint iG = 0; iG < m_nJzn; iG++ )
    { m_sumSigmaEff += ( m_vJznSigma[iG] * m_vJznEff[iG] ); }

  std::cout << "SumSigmaEff = " << m_sumSigmaEff << std::endl;

  for( uint iG = 0; iG < m_nJzn; iG++ ){
    TFile* fIn  = TFile::Open( m_vJznFnameIn[iG].c_str() );
    TTree* tree = static_cast< TTree* >( fIn->Get( "tree" ) );

    int nEventsTotal = tree->GetEntries();
    double     sigma = m_vJznSigma[iG];      
    double       eff = m_vJznEff  [iG];
  
    m_vJznNev.push_back( nEventsTotal );
    m_vJznWeights.push_back
      ( ( 1./nEventsTotal * m_sumSigmaEff ) * ( sigma * eff ) );
  
    std::cout << iG << "   weight = "
	      << m_vJznWeights.back() << std::endl;
    
    fIn->Close();
  }

  std::string system = m_is_pPb ? m_s_pPb : m_s_pp;

  std::string weightDir =
    Form( "%s/%s_%s_%s", m_sOutput.c_str(), m_sOutput.c_str(),
	  system.c_str(), m_sData.c_str() );
  
  // name of file that is used as input for performance weights
  m_fNamePerfWeightData
    = Form( "%s/%s_%s_%s_%s_%s.root",
	    weightDir.c_str(), m_myOutName.c_str(), system.c_str(),
	    m_sData.c_str(), m_sPerf.c_str(), m_uncertSuffix.c_str() );

  // name of file that is used as input for physics weights
  m_fNamePhysWeightData
    = Form( "%s/%s_%s_%s_%s_%s.root",
	    weightDir.c_str(), m_myOutName.c_str(), system.c_str(),
	    m_sData.c_str(), m_sPhys.c_str(), m_uncertSuffix.c_str() );

  // Get FCal Weights histogram
  TFile* fInFCalWeight = TFile::Open("data/pPbFCalWeight.root");
  h_pPbFCalWeights = static_cast< TH1D* >( fInFCalWeight->Get("hR") );
}

void DiJetAnalysisMC::AdditionalSuffix( std::string& label ){
  label += "_" + m_mcTypeName;
}

std::string DiJetAnalysisMC::GetMCMenu(){
  return m_is_pPb ?
    m_s_pPb + "." + m_mcTypeName :
    m_s_pp  + "." + m_mcTypeName;  
}

//---------------------------
// Fill Tree / Plot Controls
//---------------------------

void DiJetAnalysisMC::RunOverTreeFillSpect( int nEvents, 
					    int startEvent ){

  SetupHistograms();
  ProcessEventsForWeights( nEvents, startEvent, 0 );
  SaveOutputsFromTree();

  // copy to the unweighted file.
  // needed for uncertainties later
  std::cout << "Copy " << m_fNameRaw << " -> " << m_fNameRawUW << std::endl;
  TFile::Cp( m_fNameRaw.c_str(), m_fNameRawUW.c_str() );
}

void DiJetAnalysisMC::ProcessSpectWeights(){

  LoadHistograms( 1 );
  // add a slice "all" to collection
  // rest of plots include combined slices
  m_vJznLabels.push_back( m_allName );

  // open MC and Data files.
  // make spectra, and dPhi weights.
  // put into a data file (data/pp_mc_weights.root)

  TFile* fOutMCPerf = new TFile( m_fNamePerf.c_str(),"RECREATE");

  m_hAllYstarSpectFineReco =
    CombineSamples( m_vHjznYstarSpectFineReco , m_ystarSpectFineRecoName  );
  MakeSpectra( m_vHjznYstarSpectFineReco, m_vJznLabels, m_ystarSpectFineRecoName );
  fOutMCPerf->Close();

  
  // Open a TFile ( data/pp_mc_weights.root ) for writing )
  // Have function that makes weights for spect, for dPhi.
  std::string fNameOut = "data/" + m_spectName + "_" + m_sWeights + "_" + m_labelOut + ".root";
  TFile* fOut = new TFile( fNameOut.c_str(), "RECREATE" );

  // avoid memory leak.
  delete MakeSpectWeights( fOut );
  
  std::cout << "DONE! Closing " << fOut->GetName() << std::endl;
  fOut->Close();
  delete fOut;
  std::cout << "......Closed  " << std::endl;
}

void DiJetAnalysisMC::RunOverTreeFillDphi( int nEvents, 
					    int startEvent ){

  SetupHistograms();
  LoadSpectWeights();
  ProcessEventsForWeights( nEvents, startEvent, 1 );
  SaveOutputsFromTree();

  // copy to the unweighted file.
  // needed for uncertainties later
  std::cout << "Copy " << m_fNameRaw << " -> " << m_fNameRawUW << std::endl;
  TFile::Cp( m_fNameRaw.c_str(), m_fNameRawUW.c_str() );
}

void DiJetAnalysisMC::ProcessDphiWeights(){

  LoadHistograms( 1 );
  // add a slice "all" to collection
  // rest of plots include combined slices
  m_vJznLabels.push_back( m_allName );

  TFile* fInMCPerf  = TFile::Open( m_fNamePerf.c_str() );
  TFile* fOutMCPhys = new TFile( m_fNamePhys.c_str(),"RECREATE");
  m_hAllDphiReco = CombineSamples( m_vHjznDphiReco, m_dPhiRecoName   );
  MakeDeltaPhi( m_vHjznDphiReco , m_vJznLabels, m_dPhiRecoName,
		fInMCPerf, m_ystarSpectRecoName );
  fOutMCPhys->Close();
  fInMCPerf ->Close();
  
  // Open a TFile ( data/pp_mc_weights.root ) for writing )
  // Have function that makes weights for spect, for dPhi.
  std::string fNameOut = "data/" + m_dPhiName + "_" + m_sWeights + "_" + m_labelOut + ".root";
  TFile* fOut = new TFile( fNameOut.c_str(), "RECREATE" );

  // avoid memory leak.
  delete MakeDphiWeights( fOut );
  
  std::cout << "DONE! Closing " << fOut->GetName() << std::endl;
  fOut->Close();
  delete fOut;
  std::cout << "......Closed  " << std::endl;
}

void DiJetAnalysisMC::RunOverTreeFillHistos( int nEvents, 
					     int startEvent ){

  SetupHistograms();
  LoadSpectWeights();
  LoadDphiWeights();
  ProcessEvents( nEvents, startEvent );
  SaveOutputsFromTree();
}

void DiJetAnalysisMC::ProcessPerformance(){

  LoadHistograms();

  TFile* fOut = new TFile( m_fNamePerf.c_str(),"RECREATE");
  
  // add a slice "all" to collection
  // rest of plots include combined slices
  m_vJznLabels.push_back( m_allName );

  m_hAllEtaPtMap = CombineSamples( m_vHjznEtaPtMap, "etaPtMap" );
  MakeEtaPhiPtMap( m_vHjznEtaPtMap , m_vJznLabels, "etaPtMap" );
  m_hAllEtaPhiMap = CombineSamples( m_vHjznEtaPhiMap, "etaPhiMap" );
  MakeEtaPhiPtMap( m_vHjznEtaPhiMap , m_vJznLabels, "etaPhiMap" );
  
  // only do this for the default sample
  if( !m_uncertComp ){
    // need to cd into the original fout, because scale and res files are created 
    MakeScaleRes( m_vHjznRecoTruthRpt , m_vHjznRecoTruthRptNent , "recoTruthRpt"  );
    fOut->cd();
    MakeScaleRes( m_vHjznRecoTruthDeta, m_vHjznRecoTruthDetaNent, "recoTruthDeta" );
    fOut->cd();
    MakeScaleRes( m_vHjznRecoTruthDphi, m_vHjznRecoTruthDphiNent, "recoTruthDphi" );
    fOut->cd();
    MakeEfficiencies
    ( m_vHjznYstarSpectFineTruth, m_vHjznYstarSpectFineTruthUP, m_ystarEffName );
  }


  m_hAllYstarSpectReco  = CombineSamples( m_vHjznYstarSpectReco , m_ystarSpectRecoName  );
  m_hAllYstarSpectTruth = CombineSamples( m_vHjznYstarSpectTruth, m_ystarSpectTruthName );
  m_hAllYstarSpectFineReco  = CombineSamples( m_vHjznYstarSpectFineReco , m_ystarSpectFineRecoName  );
  m_hAllYstarSpectFineTruth = CombineSamples( m_vHjznYstarSpectFineTruth, m_ystarSpectFineTruthName );

  MakeSpectra( m_vHjznYstarSpectReco , m_vJznLabels, m_ystarSpectRecoName  );
  MakeSpectra( m_vHjznYstarSpectTruth, m_vJznLabels, m_ystarSpectTruthName );
  MakeSpectra( m_vHjznYstarSpectFineReco  , m_vJznLabels, m_ystarSpectFineRecoName  );
  MakeSpectra( m_vHjznYstarSpectFineTruth , m_vJznLabels, m_ystarSpectFineTruthName );

  // make ystar spectra response matrix
  m_hAllYstarSpectRespMat = CombineSamples( m_vHjznYstarSpectRespMat, m_ystarSpectRespMatName );

  MakeSpectCFactorsRespMat( m_vHjznYstarSpectReco, m_vHjznYstarSpectTruth, m_vHjznYstarSpectRespMat,
			    m_vJznLabels, m_ystarSpectCfactorsName, m_ystarSpectRespMatName );

  // make ystar response matrix
  m_hAllYstarRespMat = CombineSamples( m_vHjznYstarRespMat, m_ystarRespMatName );
  MakeYstarRespMat( m_vHjznYstarRespMat, m_vJznLabels, m_ystarRespMatName );

  if( m_is_pPb && doRtrk ){
    m_hAllJznRtrk1 = CombineSamples( m_vHjznRtrk1, Form("%s1", m_rtrkName.c_str() ) );
    m_hAllJznRtrk2 = CombineSamples( m_vHjznRtrk2, Form("%s2", m_rtrkName.c_str() ) );
    m_hAllJznRtrk4 = CombineSamples( m_vHjznRtrk4, Form("%s4", m_rtrkName.c_str() ) );

    m_hAllJznRtrk1->Write();
    m_hAllJznRtrk2->Write();
    m_hAllJznRtrk4->Write();
    
    MakeRtrk( m_vHjznRtrk1, m_vJznLabels, Form("%s1", m_rtrkName.c_str() ) );
    MakeRtrk( m_vHjznRtrk2, m_vJznLabels, Form("%s2", m_rtrkName.c_str() ) );
    MakeRtrk( m_vHjznRtrk4, m_vJznLabels, Form("%s4", m_rtrkName.c_str() ) );
  }
  
  
  std::cout << "DONE! Closing " << fOut->GetName() << std::endl;
  fOut->Close();
  delete fOut;
  std::cout << "......Closed  " << std::endl;

  fOut = new TFile( m_fNamePerf.c_str(), "UPDATE");
  // CompareAngularRes( fOut );
  CompareScaleRes( fOut, "recoTruthRpt"  );
  // CompareScaleRes( fOut, "recoTruthDeta" );
  // CompareScaleRes( fOut, "recoTruthDphi" );
  fOut->Close(); delete fOut;

  CompareRtrk( fOut );
}

void DiJetAnalysisMC::UnfoldPerformance(){

  // for now
  // Copy File with original dPhi, spectra, etc,
  // into file where histos with corrections are
  // going to be appended. 
  std::cout << "Copy " << m_fNamePerf << " -> " << m_fNamePerfUF << std::endl;
  TFile::Cp( m_fNamePerf.c_str(), m_fNamePerfUF.c_str() );

  // Open two for reading one for updating.
  // open the MC file used for unfolding info.
  // open the data file used for measured info.
  // passed to unfolding function.
  std::string fNameMC = m_fNamePerf;
  if( compToPythia && !m_is_pPb ){
    fNameMC = "output/output_pp_mc_pythia8/myOut_pp_mc_pythia8_perf_0.root";
  }

  TFile* fInData = TFile::Open( m_fNamePerf.c_str() );
  TFile* fInMC   = TFile::Open(     fNameMC.c_str() );
  TFile* fOut    = new TFile( m_fNamePerfUF.c_str(),"UPDATE");

  std::cout << "----- Unfolding MC ------" << std::endl;
  // make a vector with just the unfolded result.
  // this is to send it to MakeDeltaPhi(..) to have
  // unfolded results plotted separately
  std::vector< TH2* >        m_vHspectUnfolded;
  std::vector< std::string > m_vLabelsUnfolded;

  // make unfolded THnSparse with similar naming convention
  // as the other histograms. At this point, don't care about
  // doing this for all triggers. Altohugh, this can be
  // repeated in a loop with m_allName subsitituted for trigger,
  // and subsequently added to the vectors above.  
  TH2* m_hAllspectRecoUnfolded =
    UnfoldSpectra( fInData, fInMC, m_ystarSpectRecoUnfoldedName );
  m_vHspectUnfolded.push_back( m_hAllspectRecoUnfolded );
  m_vLabelsUnfolded.push_back( m_allName );

  // unfold on MC, just used for testing purpose.
  // so there is no comb subt or normalization or scaling
  MakeSpectra( m_vHspectUnfolded, m_vLabelsUnfolded, m_ystarSpectRecoUnfoldedName );

  std::cout << "DONE! Closing " << fOut->GetName() << std::endl;
  fOut->Close(); delete fOut;
  std::cout << "......Closed  " << std::endl;
}

void DiJetAnalysisMC::ProcessPhysics(){

  LoadHistograms();

  // THIS HAS TO BE CHANGED
  TFile* fInMCPerf  = TFile::Open( m_fNamePerf.c_str() );
  
  TFile* fOut = new TFile( m_fNamePhys.c_str(),"RECREATE");

  // add a slice "all" to collection
  // rest of plots include combined slices
  m_vJznLabels.push_back( m_allName );
  
  m_hAllDphiReco  = CombineSamples( m_vHjznDphiReco, m_dPhiRecoName   );
  MakeDeltaPhi( m_vHjznDphiReco , m_vJznLabels, m_dPhiRecoName,
		fInMCPerf, m_ystarSpectRecoName );

  m_hAllDphiTruth = CombineSamples( m_vHjznDphiTruth, m_dPhiTruthName );
  MakeDeltaPhi( m_vHjznDphiTruth, m_vJznLabels, m_dPhiTruthName,
		fInMCPerf, m_ystarSpectTruthName );

  
  m_hAllDphiRespMat    = CombineSamples( m_vHjznDphiRespMat   , m_dPhiRespMatName    );
  m_hAllDphiRespMatReb = CombineSamples( m_vHjznDphiRespMatReb, m_dPhiRespMatRebName );

  MakeDphiCFactorsRespMat( m_vHjznDphiTruth, m_vHjznDphiReco, m_vHjznDphiRespMatReb,
			   m_vJznLabels, m_dPhiCfactorsName, m_dPhiRespMatName );

  
  // MakePtRespMat( m_vHjznDphiRespMatReb, m_vJznLabels, m_ptRespMatName );

  std::cout << "DONE! Closing " << fOut->GetName() << std::endl;
  fOut->Close();
  delete fOut;
  std::cout << "......Closed  " << std::endl;
}

void DiJetAnalysisMC::UnfoldPhysics(){
  
  // for now
  // Copy File with original dPhi, spectra, etc,
  // into file where histos with corrections are
  // going to be appended. 
  std::cout << "Copy " << m_fNamePhys << " -> " << m_fNamePhysUF << std::endl;
  TFile::Cp( m_fNamePhys.c_str(), m_fNamePhysUF.c_str() );

  std::string fNameMC = m_fNamePhys;
  if( compToPythia && !m_is_pPb ){
    fNameMC = "output/output_pp_mc_pythia8/myOut_pp_mc_pythia8_phys_0.root";
  }
  
  // Open two for reading one for updating.
  // open the MC file used for unfolding info.
  // open the data file used for measured info.
  // passed to unfolding function.
  TFile* fInData  = TFile::Open(   m_fNamePhys.c_str() );
  TFile* fInMC    = TFile::Open(       fNameMC.c_str() );
  TFile* fInPerf  = TFile::Open( m_fNamePerfUF.c_str() );
  TFile* fOut     = new TFile( m_fNamePhysUF.c_str(),"UPDATE");
 
  std::cout << "----- Unfolding MC ------" << std::endl;
  // make a vector with just the unfolded result.
  // this is to send it to MakeDeltaPhi(..) to have
  // unfolded results plotted separately
  std::vector< THnSparse*  > m_vHDphiUnfolded;
  std::vector< std::string > m_vLabelsUnfolded;

  // make unfolded THnSparse with similar naming convention
  // as the other histograms. At this point, don't care about
  // doing this for all triggers. Altohugh, this can be
  // repeated in a loop with m_allName subsitituted for trigger,
  // and subsequently added to the vectors above.  
  THnSparse* m_hAllDphiRecoUnfolded =
    UnfoldDeltaPhi( fInData, fInMC, m_dPhiRecoUnfoldedName,
		    fInPerf, m_ystarSpectRecoUnfoldedName );
  m_vHDphiUnfolded .push_back( m_hAllDphiRecoUnfolded );
  m_vLabelsUnfolded.push_back( m_allName );

  // unfold on MC, just used for testing purpose.
  // make deltaPhi, give flag (true) that its unfolded response
  // so there is no comb subt or normalization or scaling
  MakeDeltaPhi( m_vHDphiUnfolded, m_vLabelsUnfolded, m_dPhiRecoUnfoldedName,
	        fInPerf, m_ystarSpectRecoUnfoldedName );
  
  std::cout << "DONE! Closing " << fOut->GetName() << std::endl;
  fOut->Close(); delete fOut;
  std::cout << "......Closed  " << std::endl;


  // fOut = new TFile( m_fNamePhysUF.c_str(),"UPDATE");
  // CompareCfactorsWUW( fOut );
  // fOut->Close(); delete fOut;
}

//---------------------------------
//            Read Data
//---------------------------------
 
void DiJetAnalysisMC::SetupHistograms(){
  
  for( uint iG = 0; iG < m_nJzn; iG++){
    std::string jzn = m_vJznLabels[iG];
    
    std::cout << "Making - " << jzn << " histograms " << std::endl;
    
    // -------- maps ---------
    m_vHjznEtaPhiMap.
      push_back( new TH2D( Form("h_etaPhiMap_%s", jzn.c_str() ),
			  ";#eta_{reco};#phi_{reco}",
			  m_nEtaMapBins, m_etaMapMin, m_etaMapMax,
			  m_nPhiMapBins, m_phiMapMin, m_phiMapMax ) );
    AddHistogram( m_vHjznEtaPhiMap.back() );
      
    m_vHjznEtaPtMap.
      push_back( new TH2D( Form("h_etaPtMap_%s", jzn.c_str() ),
			   ";#eta_{reco};#it{p}_{T}^{reco}",
			   m_nEtaMapBins, m_etaMapMin, m_etaMapMax,
			   m_nPtMapBins , m_ptMapMin , m_ptMapMax ) );
    AddHistogram( m_vHjznEtaPtMap.back() );
    
    // -------- spect --------
    m_vHjznYstarSpectReco.push_back
      ( new TH2D( Form("h_%s_%s", m_ystarSpectRecoName.c_str(), jzn.c_str() ), 
		  ";#it{y}_{1}*;#it{p}_{T1}^{reco} [GeV]",
		  m_nVarYstarBins, 0, 1,
		  m_nVarPtBinsUfOf, 0, 1 ) );
    m_vHjznYstarSpectReco.back()->GetXaxis()->
      Set( m_nVarYstarBins, &( m_varYstarBinning[0] ) );
    m_vHjznYstarSpectReco.back()->GetYaxis()->
      Set( m_nVarPtBinsUfOf, &( m_varPtBinningUfOf[0] ) );
    AddHistogram( m_vHjznYstarSpectReco.back() );
    
    m_vHjznYstarSpectTruth.push_back
      ( new TH2D( Form("h_%s_%s", m_ystarSpectTruthName.c_str(), jzn.c_str() ), 
		  ";#it{y}_{1}*;#it{p}_{T1}^{truth} [GeV]",
		  m_nVarYstarBins, 0, 1,
		  m_nVarPtBinsUfOf, 0, 1 ) );
    m_vHjznYstarSpectTruth.back()->GetXaxis()->
      Set( m_nVarYstarBins, &( m_varYstarBinning[0] ) );
    m_vHjznYstarSpectTruth.back()->GetYaxis()->
      Set( m_nVarPtBinsUfOf, &( m_varPtBinningUfOf[0] ) );
    AddHistogram( m_vHjznYstarSpectTruth.back() );

    m_vHjznYstarSpectFineReco.push_back
      ( new TH2D( Form("h_%s_%s", m_ystarSpectFineRecoName.c_str(), jzn.c_str() ), 
		  ";#eta_{reco};#it{p}_{T}^{reco} [GeV]",
		  m_nVarYstarBins, 0, 1,
		  m_nPtSpectBins, m_ptSpectMin, m_ptSpectMax ) );
    m_vHjznYstarSpectFineReco.back()->GetXaxis()->
      Set( m_nVarYstarBins, &( m_varYstarBinning[0] ) );
    AddHistogram( m_vHjznYstarSpectFineReco.back() );
    
    m_vHjznYstarSpectFineTruth.push_back
      ( new TH2D( Form("h_%s_%s", m_ystarSpectFineTruthName.c_str(), jzn.c_str() ), 
		  ";#eta_{truth};#it{p}_{T}^{truth} [GeV]",
		  m_nVarYstarBins, 0, 1,
		  m_nPtSpectBins, m_ptSpectMin, m_ptSpectMax ) );
    m_vHjznYstarSpectFineTruth.back()->GetXaxis()->
      Set( m_nVarYstarBins, &( m_varYstarBinning[0] ) );
    AddHistogram( m_vHjznYstarSpectFineTruth.back() );

    m_vHjznYstarSpectFineTruthUP.push_back
      ( new TH2D( Form("h_%s_%s", m_ystarSpectFineTruthUPName.c_str(), jzn.c_str() ), 
		  ";#eta_{truth};#it{p}_{T}^{truth} [GeV]",
		  m_nVarYstarBins, 0, 1,
		  m_nPtSpectBins, m_ptSpectMin, m_ptSpectMax ) );
    m_vHjznYstarSpectFineTruthUP.back()->GetXaxis()->
      Set( m_nVarYstarBins, &( m_varYstarBinning[0] ) );
    AddHistogram( m_vHjznYstarSpectFineTruthUP.back() );
    
    // --- spectra response matrix ----
    m_vHjznYstarSpectRespMat.push_back
      ( new TH3D( Form("h_%s_%s", m_ystarSpectRespMatName.c_str(), jzn.c_str() ), 
		  ";#it{y}_{1}*;#it{p}_{T1}^{reco} [GeV];#it{p}_{T1}^{truth} [GeV]",
		  m_nVarYstarBins , 0, 1,
		  m_nVarPtBinsUfOf, 0, 1,
		  m_nVarPtBinsUfOf, 0, 1 ) );
    m_vHjznYstarSpectRespMat.back()->GetXaxis()->
      Set( m_nVarYstarBins, &( m_varYstarBinning[0] ) );
    m_vHjznYstarSpectRespMat.back()->GetYaxis()->
      Set( m_nVarPtBinsUfOf, &( m_varPtBinningUfOf[0] ) );
    m_vHjznYstarSpectRespMat.back()->GetZaxis()->
      Set( m_nVarPtBinsUfOf, &( m_varPtBinningUfOf[0] ) );
    AddHistogram( m_vHjznYstarSpectRespMat.back() );

    // --- ystar response matrix ----
    m_vHjznYstarRespMat.push_back
      ( new TH3D( Form("h_%s_%s", m_ystarRespMatName.c_str(), jzn.c_str() ), 
		  ";#it{y}_{reco}*;#it{y}_{truth}*;#it{p}_{T}^{Truth}",
		  m_nVarYstarBins, 0, 1,
		  m_nVarYstarBins, 0, 1,
		  m_nVarPtBins   , 0, 1 ) );
    m_vHjznYstarRespMat.back()->GetXaxis()->
      Set( m_nVarYstarBins, &( m_varYstarBinning[0] ) );
    m_vHjznYstarRespMat.back()->GetYaxis()->
      Set( m_nVarYstarBins, &( m_varYstarBinning[0] ) );
    m_vHjznYstarRespMat.back()->GetZaxis()->
      Set( m_nVarPtBins, &( m_varPtBinning[0] ) );
    AddHistogram( m_vHjznYstarRespMat.back() );
    
    // --------- recoTruthRpt ---------
    int nScaleResBins = m_scaleResUseEta ?
      m_nVarEtaBins : m_nVarYstarBins;
    double* scaleResBinning = m_scaleResUseEta ?
      &m_varEtaBinning[0] : &m_varYstarBinning[0];
    
    m_vHjznRecoTruthRpt.push_back
      ( new TH3D( Form("h_recoTruthRpt_%s", jzn.c_str() ),
		  ";#it{y}_{truth}*;#it{p}_{T}^{truth} [GeV];#it{p}_{T}^{reco}/#it{p}_{T}^{truth}",
		  nScaleResBins, 0, 1,
		  m_nPtTruthBins,m_ptTruthMin, m_ptTruthMax,
		  m_nRPtRecoTruthBins, m_rPtRecoTruthMin, m_rPtRecoTruthMax) );
    m_vHjznRecoTruthRpt.back()->GetXaxis()->
      Set( nScaleResBins, scaleResBinning );
    AddHistogram( m_vHjznRecoTruthRpt.back() );
    
    m_vHjznRecoTruthRptNent.push_back
      ( new TH2D( Form("h_recoTruthRptNent_%s", jzn.c_str() ),
		  ";#it{y}_{truth}*;#it{p}_{T}^{truth} [GeV]",
		  nScaleResBins, 0, 1,
		  m_nPtTruthBins, m_ptTruthMin, m_ptTruthMax ) );
    m_vHjznRecoTruthRptNent.back()->GetXaxis()->
      Set( nScaleResBins, scaleResBinning );
    AddHistogram( m_vHjznRecoTruthRptNent.back() );    
    
    // --------- recoTruthDeta ---------
    m_vHjznRecoTruthDeta.push_back
      ( new TH3D( Form("h_recoTruthDeta_%s", jzn.c_str() ),
		  ";#it{y}_{truth}*;#it{p}_{T}^{truth};#eta^{reco}-#eta^{truth}",
		  nScaleResBins, 0, 1,
		  m_nPtTruthBins,m_ptTruthMin, m_ptTruthMax,
		  m_nDAngleRecoTruthBins, m_dAngleRecoTruthMin, m_dAngleRecoTruthMax ) );
    m_vHjznRecoTruthDeta.back()->GetXaxis()->
      Set( nScaleResBins, scaleResBinning );
    AddHistogram( m_vHjznRecoTruthDeta.back() );    
    
    m_vHjznRecoTruthDetaNent.push_back
      ( new TH2D( Form("h_recoTruthDetaNent_%s", jzn.c_str() ),
		  ";#it{y}_{truth};#it{p}_{T}^{truth}",
		  nScaleResBins, 0, 1,
		  m_nPtTruthBins, m_ptTruthMin, m_ptTruthMax ) );
    m_vHjznRecoTruthDetaNent.back()->GetXaxis()->
      Set( nScaleResBins, scaleResBinning );
    AddHistogram( m_vHjznRecoTruthDetaNent.back() );    

    // --------- recoTruthDphi ---------
    m_vHjznRecoTruthDphi.push_back
      ( new TH3D( Form("h_recoTruthDphi_%s", jzn.c_str() ),
		  ";#it{y}_{truth}*;#it{p}_{T}^{truth};#phi^{reco}-#phi^{truth}",
		  nScaleResBins, 0, 1,
		  m_nPtTruthBins, m_ptTruthMin, m_ptTruthMax,
		  m_nDAngleRecoTruthBins, m_dAngleRecoTruthMin, m_dAngleRecoTruthMax ) );
    m_vHjznRecoTruthDphi.back()->GetXaxis()->
      Set( nScaleResBins, scaleResBinning );
    AddHistogram( m_vHjznRecoTruthDphi.back() );
    
    m_vHjznRecoTruthDphiNent.push_back
      ( new TH2D( Form("h_recoTruthDphiNent_%s", jzn.c_str() ),
		  ";#it{y}_{truth}*;#it{p}_{T}^{truth}",
		  nScaleResBins, 0, 1,
		  m_nPtTruthBins, m_ptTruthMin, m_ptTruthMax ) );
    m_vHjznRecoTruthDphiNent.back()->GetXaxis()->
      Set( nScaleResBins, scaleResBinning );
    AddHistogram( m_vHjznRecoTruthDphiNent.back() );

    // rtrk
    m_vHjznRtrk1.push_back
      ( new TH3D( Form("h_%s1_%s", m_rtrkName.c_str(), jzn.c_str() ), 
		  ";#eta;#it{p}_{T} [GeV];Rtrk",
		  m_nVarEtaBarrelBins, 0, 1,
		  m_nVarRtrkPtBins, 0, 1,
		  50, 0, 1 ) ) ;
    m_vHjznRtrk1.back()->GetXaxis()->
      Set( m_nVarEtaBarrelBins, &( m_varEtaBarrelBinning[0] ) );
    m_vHjznRtrk1.back()->GetYaxis()->
      Set( m_nVarRtrkPtBins, &( m_varRtrkPtBinning[0] ) );
    AddHistogram( m_vHjznRtrk1.back() );
   
    m_vHjznRtrk2.push_back
      ( new TH3D( Form("h_%s2_%s", m_rtrkName.c_str(), jzn.c_str() ), 
		  ";#eta;#it{p}_{T} [GeV];Rtrk",
		  m_nVarEtaBarrelBins, 0, 1,
		  m_nVarRtrkPtBins, 0, 1,
		  50, 0, 1 ) ) ;
    m_vHjznRtrk2.back()->GetXaxis()->
      Set( m_nVarEtaBarrelBins, &( m_varEtaBarrelBinning[0] ) );
    m_vHjznRtrk2.back()->GetYaxis()->
      Set( m_nVarRtrkPtBins, &( m_varRtrkPtBinning[0] ) );
    AddHistogram( m_vHjznRtrk2.back() );

    m_vHjznRtrk4.push_back
      ( new TH3D( Form("h_%s4_%s", m_rtrkName.c_str(), jzn.c_str() ), 
		  ";#eta;#it{p}_{T} [GeV];Rtrk",
		  m_nVarEtaBarrelBins, 0, 1,
		  m_nVarRtrkPtBins, 0, 1,
		  50, 0, 1 ) ) ;
    m_vHjznRtrk4.back()->GetXaxis()->
      Set( m_nVarEtaBarrelBins, &( m_varEtaBarrelBinning[0] ) );
    m_vHjznRtrk4.back()->GetYaxis()->
      Set( m_nVarRtrkPtBins, &( m_varRtrkPtBinning[0] ) );
    AddHistogram( m_vHjznRtrk4.back() );

    // -------- dPhi --------
    THnSparse* hnDphiReco =
      new THnSparseD( Form("h_%s_%s", m_dPhiRecoName.c_str(), jzn.c_str() ), "",
		      m_nDphiDim, &m_vNdPhiBins[0], &m_vDphiMin[0], &m_vDphiMax[0] );
    hnDphiReco->GetAxis(0)->Set( m_nVarYstarBins, &( m_varYstarBinning[0] ) );
    hnDphiReco->GetAxis(1)->Set( m_nVarYstarBins, &( m_varYstarBinning[0] ) );
    hnDphiReco->GetAxis(2)->Set( m_nVarPtBins   , &( m_varPtBinning[0]    ) );
    hnDphiReco->GetAxis(3)->Set( m_nVarPtBins   , &( m_varPtBinning[0]    ) );
    hnDphiReco->GetAxis(4)->Set( m_nVarDphiBins , &( m_varDphiBinning[0]  ) );
    m_vHjznDphiReco.push_back( hnDphiReco );
    AddHistogram( hnDphiReco );
    
    THnSparse* hnDphiTruth =
      new THnSparseD( Form("h_%s_%s", m_dPhiTruthName.c_str(), jzn.c_str() ), "",
		      m_nDphiDim, &m_vNdPhiBins[0], &m_vDphiMin[0], &m_vDphiMax[0] );
    hnDphiTruth->GetAxis(0)->Set( m_nVarYstarBins, &( m_varYstarBinning[0] ) );
    hnDphiTruth->GetAxis(1)->Set( m_nVarYstarBins, &( m_varYstarBinning[0] ) );
    hnDphiTruth->GetAxis(2)->Set( m_nVarPtBins   , &( m_varPtBinning[0]    ) );
    hnDphiTruth->GetAxis(3)->Set( m_nVarPtBins   , &( m_varPtBinning[0]    ) );
    hnDphiTruth->GetAxis(4)->Set( m_nVarDphiBins , &( m_varDphiBinning[0]  ) );
    m_vHjznDphiTruth.push_back( hnDphiTruth );
    AddHistogram( hnDphiTruth );

    // -------- Dphi Response Matrix --------    
    THnSparse* hnDphiRespMat =
      new THnSparseD( Form("h_%s_%s", m_dPhiRespMatName.c_str(), jzn.c_str() ), "",
		      m_nDphiRespMatDim, &m_vNdPhiRespMatBins[0],
		      &m_vDphiRespMatMin[0], &m_vDphiRespMatMax[0] );
    hnDphiRespMat->GetAxis(0)->Set( m_nVarYstarBins , &( m_varYstarBinning[0]  ) );
    hnDphiRespMat->GetAxis(1)->Set( m_nVarYstarBins , &( m_varYstarBinning[0]  ) );
    hnDphiRespMat->GetAxis(2)->Set( m_nVarPtBinsUfOf, &( m_varPtBinningUfOf[0] ) );
    hnDphiRespMat->GetAxis(3)->Set( m_nVarPtBinsUfOf, &( m_varPtBinningUfOf[0] ) );
    hnDphiRespMat->GetAxis(4)->Set( m_nVarPtBinsUfOf, &( m_varPtBinningUfOf[0] ) );
    hnDphiRespMat->GetAxis(5)->Set( m_nVarPtBinsUfOf, &( m_varPtBinningUfOf[0] ) );
    hnDphiRespMat->GetAxis(6)->Set( m_nVarDphiBins  , &( m_varDphiBinning[0]   ) );
    hnDphiRespMat->GetAxis(7)->Set( m_nVarDphiBins  , &( m_varDphiBinning[0]   ) );

    m_vHjznDphiRespMat.push_back( hnDphiRespMat );
    AddHistogram( hnDphiRespMat );
    
    // change some of the axis names
    hnDphiRespMat->GetAxis(2)->SetTitle( "Reco #it{p}_{T,1}"  );
    hnDphiRespMat->GetAxis(3)->SetTitle( "Truth #it{p}_{T,1}" );
    hnDphiRespMat->GetAxis(4)->SetTitle( "Reco #it{p}_{T,2}"  );
    hnDphiRespMat->GetAxis(5)->SetTitle( "Truth #it{p}_{T,2}" );
    hnDphiRespMat->GetAxis(6)->SetTitle( "#Delta#phi_{reco}"  );
    hnDphiRespMat->GetAxis(7)->SetTitle( "#Delta#phi_{truth}" );
    
    // --- Rebinned Dphi Response Matrix -----    
    THnSparse* hnDphiRespMatReb =
      new THnSparseD( Form("h_%s_%s", m_dPhiRespMatRebName.c_str(), jzn.c_str() ), "",
		      m_nDphiRespMatDim, &m_vNdPhiRespMatRebBins[0],
		      &m_vDphiRespMatMin[0], &m_vDphiRespMatMax[0] );
    hnDphiRespMatReb->GetAxis(0)->Set( m_nVarYstarBins , &( m_varYstarBinning[0]  ) );
    hnDphiRespMatReb->GetAxis(1)->Set( m_nVarYstarBins , &( m_varYstarBinning[0]  ) );
    hnDphiRespMatReb->GetAxis(2)->Set( m_nVarPtBinsUfOf, &( m_varPtBinningUfOf[0] ) );
    hnDphiRespMatReb->GetAxis(3)->Set( m_nVarPtBinsUfOf, &( m_varPtBinningUfOf[0] ) );
    hnDphiRespMatReb->GetAxis(4)->Set( m_nVarPtBinsUfOf, &( m_varPtBinningUfOf[0] ) );
    hnDphiRespMatReb->GetAxis(5)->Set( m_nVarPtBinsUfOf, &( m_varPtBinningUfOf[0] ) );
    hnDphiRespMatReb->GetAxis(6)->Set( m_nVarDphiRebinnedBins, &( m_varDphiRebinnedBinning[0] ) );
    hnDphiRespMatReb->GetAxis(7)->Set( m_nVarDphiRebinnedBins, &( m_varDphiRebinnedBinning[0] ) );

    m_vHjznDphiRespMatReb.push_back( hnDphiRespMatReb );
    AddHistogram( hnDphiRespMatReb );
    
    // change some of the axis names
    hnDphiRespMatReb->GetAxis(2)->SetTitle( "Reco #it{p}_{T,1}"  );
    hnDphiRespMatReb->GetAxis(3)->SetTitle( "Truth #it{p}_{T,1}" );
    hnDphiRespMatReb->GetAxis(4)->SetTitle( "Reco #it{p}_{T,2}"  );
    hnDphiRespMatReb->GetAxis(5)->SetTitle( "Truth #it{p}_{T,2}" );
    hnDphiRespMatReb->GetAxis(6)->SetTitle( "#Delta#phi_{reco}"  );
    hnDphiRespMatReb->GetAxis(7)->SetTitle( "#Delta#phi_{truth}" );    
  } 
}

void DiJetAnalysisMC::LoadSpectWeights(){

  std::string fNameWeights = "data/" + m_spectName + "_" + m_sWeights + "_" + m_labelOut + ".root";
  TFile* fWeights = TFile::Open( fNameWeights.c_str() );

  std::string hName = "h_" + m_spectName + "_" + m_sWeights + "_" + m_allName;

  std::string axisLabel, axisLabelTex;
  GetSpectraLabels( axisLabel, axisLabelTex, m_sYstar );
  
  for( uint xBin = 1; xBin <= m_nVarYstarBins; xBin++ ){
    double xLow  = m_varYstarBinning[ xBin - 1 ];
    double xUp   = m_varYstarBinning[ xBin ];

    std::string hTag    = anaTool->GetName( xLow, xUp, axisLabel);
    std::string fitName = Form( "f_h_%s_%s_%s_%s",
				m_spectName.c_str(), m_sWeights.c_str(),
				m_allName.c_str(), hTag.c_str() );

    m_vSpectWeightFits.push_back
      ( static_cast< TF1* >( fWeights->Get( fitName.c_str() ) ) );
  }

  fWeights->Close();
}

void DiJetAnalysisMC::LoadDphiWeights(){

  std::string fNameWeights = "data/" + m_dPhiName + "_" + m_sWeights + "_" + m_labelOut + ".root";
  TFile* fWeights = TFile::Open( fNameWeights.c_str() );

  m_dPhiWeight  = static_cast< TH3* >
    ( fWeights->Get( Form( "h_%s_%s_%s", m_dPhiName.c_str(), m_sWeights.c_str(), m_allName.c_str() ) ) );

  m_dPhiWeight->SetDirectory(0);
  
  fWeights->Close();
}

TH2* hW = new TH2D( "hW", "hW", 16, 0, 16, 100, 0, 10 );

void DiJetAnalysisMC::ProcessEventsForWeights( int nEventsIn, int startEventIn, int mode ){ 

  // collections and variables
  std::vector< TLorentzVector >    vT_jets;
  std::vector< TLorentzVector >* p_vT_jets = &vT_jets;
  
  std::vector< TLorentzVector >    vR_jets;
  std::vector< TLorentzVector >* p_vR_jets = &vR_jets;
  
  std::vector< bool >    v_isCleanJet;
  std::vector< bool >* p_v_isCleanJet = &v_isCleanJet;

  for( uint iG = 0; iG < m_nJzn; iG++){

    std::cout << "fNameIn: " << m_vJznFnameIn[iG] << std::endl;
  
    TFile* fIn  = TFile::Open( m_vJznFnameIn[iG].c_str() );
    TTree* tree = static_cast< TTree*>( fIn->Get( "tree" ) );

    // Connect to tree
    tree->SetBranchAddress( "vT_jets"      , &p_vT_jets      );
    tree->SetBranchAddress( "vR_C_jets"    , &p_vR_jets      );
    tree->SetBranchAddress( "v_isCleanJet" , &p_v_isCleanJet );

    if( m_is_pPb ){
      m_FCalEt = tree->SetBranchAddress( "FCalEtA"  , &m_FCalEt );
    } else {
      m_FCalEt = 0;
    }

    // n events
    int nEvents, startEvent, nEventsTotal, endEvent;
    
    nEventsTotal = tree->GetEntries();
    nEvents      = nEventsIn > 0 ? nEventsIn : nEventsTotal;
    startEvent   = startEventIn < nEventsTotal ?
				  startEventIn : nEventsTotal - 1;
    endEvent     = startEvent + nEvents < nEventsTotal ?
					  startEvent + nEvents : nEventsTotal;
    
    std::cout << startEvent << " " << endEvent << " "
	      << nEventsTotal << " " << nEvents << std::endl;
    
    // -------- EVENT LOOP ---------
    for( m_ev = startEvent; m_ev < endEvent; m_ev++ ){
      tree->GetEntry( m_ev );
      
      if( anaTool->DoPrint( m_ev ) ) {
	std::cout << "\nEvent : " << m_ev 
		  << "    has : " << vR_jets.size() << " reco jets"
		  << "    and : " << vT_jets.size() << " truth jets"
		  << std::endl;
      }

      ApplyCleaning ( vR_jets, v_isCleanJet );
      // ApplyIsolation( vR_jets, 1.0 );
      // ApplyIsolation( vT_jets, 1.0 );
      
      std::vector< TLorentzVector > vTR_paired_jets;
      std::vector< TLorentzVector > vTT_paired_jets;
      PairJets( vT_jets, vR_jets, vTT_paired_jets, vTR_paired_jets );

      std::sort( vTR_paired_jets.begin(), vTR_paired_jets.end(),
		 anaTool->sortByDecendingPt );
      std::sort( vTT_paired_jets.begin(), vTT_paired_jets.end(),
		 anaTool->sortByDecendingPt );
      
      // Do Dphi analysis
      if( mode == 0 ){
	AnalyzeSpectra( m_vHjznYstarSpectFineReco[iG], vTR_paired_jets );	
      } else if( mode == 1 ){
	AnalyzeDeltaPhiWithWeight
	  ( m_vHjznDphiReco[iG], vTR_paired_jets, vTT_paired_jets, mode );
      }      
    } // -------- END EVENT LOOP ---------
   
    std::cout << "DONE WITH " << m_vJznLabels[iG] << std::endl;
    
    fIn->Close(); delete fIn;
  } // end loop over a JZ sample
}

void DiJetAnalysisMC::ProcessEvents( int nEventsIn, int startEventIn ){ 

  // collections and variables
  std::vector< TLorentzVector >    vT_jets;
  std::vector< TLorentzVector >* p_vT_jets = &vT_jets;
  
  std::vector< TLorentzVector >    vR_jets;
  std::vector< TLorentzVector >* p_vR_jets = &vR_jets;
  
  std::vector< bool >    v_isCleanJet;
  std::vector< bool >* p_v_isCleanJet = &v_isCleanJet;

  std::vector< std::vector< float > >    v_sysUncert;
  std::vector< std::vector< float > >* p_v_sysUncert = & v_sysUncert;

  std::vector< double > vRtrk1;
  std::vector< double > vRtrk2;
  std::vector< double > vRtrk4;

  std::vector< double >* p_vRtrk1 = &vRtrk1;
  std::vector< double >* p_vRtrk2 = &vRtrk2;
  std::vector< double >* p_vRtrk4 = &vRtrk4;
  
  TH2* hJerComp = new TH2D( "hJerComp", ";#it{p}_T;After/Before", 124, 28, 90, 50, 0, 2 );
  styleTool->SetHStyle( hJerComp, 0 );
  
  for( uint iG = 0; iG < m_nJzn; iG++){

    std::cout << "fNameIn: " << m_vJznFnameIn[iG] << std::endl;
  
    TFile* fIn  = TFile::Open( m_vJznFnameIn[iG].c_str() );
    TTree* tree = static_cast< TTree*>( fIn->Get( "tree" ) );

    // Connect to tree
    tree->SetBranchAddress( "vT_jets"      , &p_vT_jets      );
    tree->SetBranchAddress( "vR_C_jets"    , &p_vR_jets      );
    tree->SetBranchAddress( "v_isCleanJet" , &p_v_isCleanJet );
    tree->SetBranchAddress( "v_sysUncert"  , &p_v_sysUncert  );

    if( m_is_pPb && doRtrk ){
      tree->SetBranchAddress( "vRtrk1", &p_vRtrk1 );
      tree->SetBranchAddress( "vRtrk2", &p_vRtrk2 );
      tree->SetBranchAddress( "vRtrk4", &p_vRtrk4 );
    }
    
    if( m_is_pPb ){
      m_FCalEt = tree->SetBranchAddress( "FCalEtA"  , &m_FCalEt );
    } else {
      m_FCalEt = 0;
    }

    // n events
    int nEvents, startEvent, nEventsTotal, endEvent;
    
    nEventsTotal = tree->GetEntries();
    nEvents      = nEventsIn > 0 ? nEventsIn : nEventsTotal;
    startEvent   = startEventIn < nEventsTotal ?
				  startEventIn : nEventsTotal - 1;
    endEvent     = startEvent + nEvents < nEventsTotal ?
					  startEvent + nEvents : nEventsTotal;
    
    std::cout << startEvent << " " << endEvent << " "
	      << nEventsTotal << " " << nEvents << std::endl;
    
    // -------- EVENT LOOP ---------
    for( m_ev = startEvent; m_ev < endEvent; m_ev++ ){
      tree->GetEntry( m_ev );
            
      if( anaTool->DoPrint( m_ev ) ) {
	std::cout << "\nEvent : " << m_ev 
		  << "    has : " << vR_jets.size() << " reco jets"
		  << "    and : " << vT_jets.size() << " truth jets"
		  << std::endl;
      }
      
      ApplyCleaning ( vR_jets, v_isCleanJet );
      // ApplyIsolation( vR_jets, 1.0 );
      // ApplyIsolation( vT_jets, 1.0 );
            
      std::vector< TLorentzVector > vTR_paired_jets;
      std::vector< TLorentzVector > vTT_paired_jets;
      PairJets( vT_jets, vR_jets, vTT_paired_jets, vTR_paired_jets );

      // If not running on default sample.
      // Apply uncertainties to all reco jets.
      if( m_uncertComp ){
	m_uncertaintyProvider->RegisterUFactors  ( &v_sysUncert );
	m_uncertaintyProvider->ApplyUncertainties( vTR_paired_jets, vTT_paired_jets );
      }
      
      // analyze ystar resp mat BEFORE sorting
      // done to not include JER effects.
      AnalyzeYstarRespMat( m_vHjznYstarRespMat[iG],
			   vTR_paired_jets, vTT_paired_jets );

      // do rtrk stuff in pPb before sorting
      if( m_is_pPb && doRtrk ){
	
	for( uint iJet = 0; iJet < vR_jets.size(); iJet++ ){

	  if( !PassHECCuts( vR_jets[iJet] ) ){ continue; }

	  double jetPt   = vR_jets[iJet].Pt() / 1000.;
	  double jetEta  = vR_jets[iJet].Eta();

	  double rTrkPt1 = vRtrk1[iJet] / 1000.;
	  double rTrkPt2 = vRtrk2[iJet] / 1000.;
	  double rTrkPt4 = vRtrk4[iJet] / 1000.;

	  if( std::abs( jetEta ) > 2.5 || !jetPt ){ continue; }

	  if( rTrkPt1 ){
	    m_vHjznRtrk1[iG]->Fill( jetEta, jetPt, jetPt/rTrkPt1 );
	  } if( rTrkPt2 ){
	    m_vHjznRtrk2[iG]->Fill( jetEta, jetPt, jetPt/rTrkPt2 );
	  } if( rTrkPt4 ){
	    m_vHjznRtrk4[iG]->Fill( jetEta, jetPt, jetPt/rTrkPt4 );
	  } 
	}
      }
    
      
      std::sort( vTR_paired_jets.begin(), vTR_paired_jets.end(),
		 anaTool->sortByDecendingPt );
      std::sort( vTT_paired_jets.begin(), vTT_paired_jets.end(),
		 anaTool->sortByDecendingPt );

      // last parameter is mode, we use 2 to get dphi weight.
      int mode = 2;
      AnalyzeDeltaPhiWithWeight( m_vHjznDphiReco [iG], vTR_paired_jets, vTT_paired_jets, mode );
      AnalyzeDeltaPhiWithWeight( m_vHjznDphiTruth[iG], vTT_paired_jets, vTT_paired_jets, mode, false );
      
      AnalyzeDphiRespMat
        ( m_vHjznDphiRespMat[iG], m_vHjznDphiRespMatReb[iG], vTR_paired_jets, vTT_paired_jets );
      	          
      // fill single jet spectra 
      if( vTT_paired_jets.size() ){
	
	TLorentzVector& rJetFront = vTR_paired_jets.front();
	TLorentzVector& tJetFront = vTT_paired_jets.front();
	
	double recoJetWeight  = GetJetWeight( rJetFront ) * GetSpectWeight( tJetFront );
	double truthJetWeight = GetJetWeight( tJetFront ) * GetSpectWeight( tJetFront );
	
	m_vHjznYstarSpectReco[iG] ->Fill
	  ( GetYstar( rJetFront ), rJetFront.Pt()/1000., recoJetWeight  );
	m_vHjznYstarSpectTruth[iG]->Fill
	  ( GetYstar( tJetFront ), tJetFront.Pt()/1000., truthJetWeight );

	// fill both sides for pp
	if( !m_is_pPb ){
	  m_vHjznYstarSpectReco[iG] ->Fill
	    ( -GetYstar( rJetFront ), rJetFront.Pt()/1000., recoJetWeight  );
	  m_vHjznYstarSpectTruth[iG]->Fill
	    ( -GetYstar( tJetFront ), tJetFront.Pt()/1000., truthJetWeight );
	}
      }

      // fill single speectra response matrix
      AnalyzeSpectRespMat( m_vHjznYstarSpectRespMat[iG],
			   vTR_paired_jets, vTT_paired_jets );

      // make spectra
      AnalyzeSpectra( m_vHjznYstarSpectFineTruth[iG], vTT_paired_jets );

      // for reco fine, use weights (for comparison later)
      for( uint iJet = 0; iJet < vTR_paired_jets.size(); iJet++ ){
	TLorentzVector& rJet = vTR_paired_jets[ iJet ];
	TLorentzVector& tJet = vTT_paired_jets[ iJet ];

	double jetYstar  = GetYstar( rJet );
	double jetPt     = rJet.Pt()/1000.;
    
	double jetWeight = GetJetWeight( rJet ) * GetSpectWeight( tJet );

	m_vHjznYstarSpectFineReco[iG]->Fill( jetYstar, jetPt, jetWeight );

	// for pp fill both sides
	if( !m_is_pPb ){
	  m_vHjznYstarSpectFineReco[iG]->Fill( -jetYstar, jetPt, jetWeight );
	}
      }
            
      // for efficiencies. These have different binning
      AnalyzeSpectra( m_vHjznYstarSpectFineTruthUP[iG], vT_jets );
      
      // do JER/JES, angular scales and resolution.
      AnalyzeScaleResolution( vTR_paired_jets, vTT_paired_jets, iG );
    } // -------- END EVENT LOOP ---------
   
    std::cout << "DONE WITH " << m_vJznLabels[iG] << std::endl;

    fIn->Close(); delete fIn;
  } // end loop over a JZ sample
  TFile* f = new TFile("myOut.root", "RECREATE" );
  hW->Write();
  hJerComp->Write();
  f->Close();   
}

//---------------------------
//          Analysis
//---------------------------

double DiJetAnalysisMC::AnalyzeDeltaPhiWithWeight
( THnSparse* hn,
  const std::vector< TLorentzVector >& v_jets,
  const std::vector< TLorentzVector >& vW_jets,
  int mode,
  bool isReco){

  const TLorentzVector* jet1  = NULL; const TLorentzVector* jet2  = NULL;
  const TLorentzVector* wJet1 = NULL; const TLorentzVector* wJet2 = NULL;

  if( !GetDiJets(  v_jets,  jet1,  jet2 ) ||
      !GetDiJets( vW_jets, wJet1, wJet2, false ) )
    { return -1; }

  // in pPb, quit if there is a subleading
  // jet in Pb going HEC region
  if( isReco && !PassHECCuts( *jet2 ) ){ return -1; }
  
  double jetPt1    = jet1->Pt()/1000.;
  double jetYstar1 = GetYstar( *jet1 );
  double jetPt2    = jet2->Pt()/1000.;
  double jetYstar2 = GetYstar( *jet2 );

  double dPhi      = anaTool->DeltaPhi(  *jet2,  *jet1 );
  // double jetWeight = GetJetWeight( *jet1 ) * GetSpectWeight( *wJet1 );
  double jetWeight = GetJetWeight( *jet1 ) * GetSpectWeight( *wJet1 );
  if( mode == 2 ){
    jetWeight *= GetDphiWeight( *wJet1, *wJet2 );
  }
    
  std::vector< double > p{ jetYstar1, jetYstar2, jetPt1, jetPt2, dPhi, 0.5 };
  
  hn->Fill( &p[0], jetWeight );

  // for pp, fill twice. once for each side since
  // it is symmetric in pp. For pPb, continue
  if( m_is_pPb ){ return dPhi; }

  p[0] = -jetYstar1;  
  p[1] = -jetYstar2;
  hn->Fill( &p[0], jetWeight );
 
  return dPhi;
}


void DiJetAnalysisMC::AnalyzeScaleResolution( const std::vector< TLorentzVector >& vR_jets,
					      const std::vector< TLorentzVector >& vT_jets,
					      const int iG ){
  // loop over jets
  for( uint iJet = 0; iJet < vR_jets.size(); iJet ++ ){ 
    const TLorentzVector& rJet = vR_jets[iJet];
    const TLorentzVector& tJet = vT_jets[iJet];

    // reco jet
    double    jetEtaReco  = rJet.Eta();
    double    jetPhiReco  = rJet.Phi();
    double     jetPtReco  = rJet.Pt()/1000.;
    //  double jetWeightReco  = GetJetWeight( rJet );
    double jetWeightReco  = 1;

    // truth jet
    double    jetEtaTruth = tJet.Eta();
    double    jetPhiTruth = tJet.Phi();
    double     jetPtTruth = tJet.Pt()/1000.;
    double  jetYstarTruth = GetYstar( tJet );
    // double jetWeightTruth = GetJetWeight( tJet );
    double jetWeightTruth = 1;

    if( jetPtReco > 28 && jetPtReco < 35 ){
      // fill negative because the coordinate system is
      // with p going in positive eta
      m_vHjznEtaPhiMap[iG]->Fill( jetEtaReco, jetPhiReco, jetWeightReco );
    }
    
    m_vHjznEtaPtMap [iG]->Fill( jetEtaReco, jetPtReco , jetWeightReco );

    double angle = m_scaleResUseEta ? jetEtaTruth : jetYstarTruth;
    
    if( PassHECCuts( rJet) ){
      m_vHjznRecoTruthRpt       [iG]->
	Fill( angle, jetPtTruth, jetPtReco/jetPtTruth, jetWeightTruth );
      m_vHjznRecoTruthRptNent   [iG]->
	Fill( angle, jetPtTruth );

      m_vHjznRecoTruthDeta      [iG]->
	Fill( angle, jetPtTruth, jetEtaReco - jetEtaTruth, jetWeightTruth );
      m_vHjznRecoTruthDetaNent  [iG]->
	Fill( angle, jetPtTruth );

      m_vHjznRecoTruthDphi      [iG]->
	Fill( angle, jetPtTruth, jetPhiReco - jetPhiTruth, jetWeightTruth );
      m_vHjznRecoTruthDphiNent  [iG]->
	Fill( angle, jetPtTruth  );
    }

    /*
    // for pp, fill both
    if( m_is_pPb ){ continue; }

    if( jetEtaReco > 28 && jetEtaReco < 35 ){
      m_vHjznEtaPhiMap[iG]->Fill( -jetEtaReco, jetPhiReco, jetWeightReco );
    }
    m_vHjznEtaPtMap [iG]->Fill( -jetEtaReco, jetPtReco , jetWeightReco );

    if( PassHECCuts( rJet) ){
      m_vHjznRecoTruthRpt       [iG]->
	Fill( -angle, jetPtTruth, jetPtReco/jetPtTruth, jetWeightTruth );
      m_vHjznRecoTruthRptNent   [iG]->
	Fill( -angle, jetPtTruth );

      m_vHjznRecoTruthDeta      [iG]->
	Fill( -angle, jetPtTruth, jetEtaReco - jetEtaTruth, jetWeightTruth );
      m_vHjznRecoTruthDetaNent  [iG]->
	Fill( -angle, jetPtTruth );

      m_vHjznRecoTruthDphi      [iG]->
	Fill( -angle, jetPtTruth, jetPhiReco - jetPhiTruth, jetWeightTruth );
      m_vHjznRecoTruthDphiNent  [iG]->
	Fill( -angle, jetPtTruth  );
    }
    */
  } // end loop over pairs
}

void DiJetAnalysisMC::AnalyzeYstarRespMat( TH3* hRespMatYstar,
					   const std::vector< TLorentzVector >& vR_jets,
					   const std::vector< TLorentzVector >& vT_jets ){
  
  if( vR_jets.size() && vT_jets.size() ){

    const TLorentzVector* rJet = &vR_jets.front();
    const TLorentzVector* tJet = &vT_jets.front();

    // reco jet
    double rJet_ystar = GetYstar( *rJet );

    // truth jet
    double tJet_pt    = tJet->Pt()/1000.;
    double tJet_ystar = GetYstar( *tJet );
    double jetWeight = GetJetWeight( *tJet );
    
    hRespMatYstar->Fill(  tJet_ystar, rJet_ystar, tJet_pt, jetWeight );
    
    // for pp fill plus minus ystar
    if( m_is_pPb ){ return; }

    hRespMatYstar->Fill( -tJet_ystar, -rJet_ystar, tJet_pt, jetWeight );
  }
}


void DiJetAnalysisMC::AnalyzeSpectRespMat( TH3* hRespMatPt,
					   const std::vector< TLorentzVector >& vR_jets,
					   const std::vector< TLorentzVector >& vT_jets ){
  
  if( vR_jets.size() && vT_jets.size() ){

    const TLorentzVector* rJet = &vR_jets.front();
    const TLorentzVector* tJet = &vT_jets.front();

    // reco jet
    double rJet_pt    = rJet->Pt()/1000.;

    // truth jet
    double tJet_pt    = tJet->Pt()/1000.;
    double tJet_ystar = GetYstar( *tJet );
    double jetWeight = GetJetWeight( *tJet );
    
    hRespMatPt->Fill( tJet_ystar, rJet_pt, tJet_pt, jetWeight );
    
    // for pp fill plus minus ystar
    if( m_is_pPb ){ return; }

    hRespMatPt->Fill( -tJet_ystar, rJet_pt, tJet_pt, jetWeight );
  }
}

void DiJetAnalysisMC::AnalyzeDphiRespMat( THnSparse* hnDphi,
					  THnSparse* hnDphiReb,
					  const std::vector< TLorentzVector >& vR_jets,
					  const std::vector< TLorentzVector >& vT_jets ){

  const TLorentzVector* rJet1 = NULL; const TLorentzVector* rJet2 = NULL;
  const TLorentzVector* tJet1 = NULL; const TLorentzVector* tJet2 = NULL;

  if( !GetDiJets( vR_jets, rJet1, rJet2 ) ||
      !GetDiJets( vT_jets, tJet1, tJet2 ) ) { return; }

  // reco jet 1
  double rjetPt1    = rJet1->Pt()/1000.;

  // reco jet 2
  double rjetPt2    = rJet2->Pt()/1000.;

  // in pPb, quit if there is a subleading
  // reco jet in Pb going HEC region
  if( !PassHECCuts( *rJet2 ) ){ return; }
  
  // truth jet 1
  double tjetPt1    = tJet1->Pt()/1000.;
  double tjetYstar1 = GetYstar( *tJet1 );

  // truth jet 2
  double tjetPt2    = tJet2->Pt()/1000.;
  double tjetYstar2 = GetYstar( *tJet2 );

  // dphi reco and truth
  double recoDphi    = anaTool->DeltaPhi( *rJet2, *rJet1 );
  double truthDphi   = anaTool->DeltaPhi( *tJet2, *tJet1 );
  
  std::vector< double > xDphi( hnDphi->GetNdimensions(), 0 );
  
  double jetWeight = GetJetWeight( *tJet1 ) *
    GetSpectWeight( *tJet1 ) * GetDphiWeight( *tJet1, *tJet2 );
 
  // fill for Dphi resp mat
  xDphi[0] = tjetYstar1;  
  xDphi[1] = tjetYstar2;
  xDphi[2] = rjetPt1 ;
  xDphi[3] = tjetPt1;
  xDphi[4] = rjetPt2 ;
  xDphi[5] = tjetPt2;
  xDphi[6] = recoDphi;
  xDphi[7] = truthDphi;
  hnDphi   ->Fill( &xDphi[0], jetWeight );
  hnDphiReb->Fill( &xDphi[0], jetWeight );
  
  // for pp, fill twice. once for each side since
  // it is symmetric in pp. For pPb, continue
  if( m_is_pPb ){ return; }
  
  // fill for Dphi resp mat
  xDphi[0] = -tjetYstar1;  
  xDphi[1] = -tjetYstar2;
  hnDphi   ->Fill( &xDphi[0], jetWeight );
  hnDphiReb->Fill( &xDphi[0], jetWeight );
}

//---------------------------
//          Tools 
//---------------------------

// pair reco and truth jets with deltaR parameter 
void DiJetAnalysisMC::PairJets( std::vector< TLorentzVector >& vA_jets,
			        std::vector< TLorentzVector >& vB_jets,
				std::vector< TLorentzVector >& vA_paired_jets,
				std::vector< TLorentzVector >& vB_paired_jets ){
  // exit if either is empty
  if( !vA_jets.size() || !vB_jets.size() ){ return; }
  
  for( auto& aJet : vA_jets){
    
    // also ignore jets outside of a certain ystar range.
    double aJetYstar = GetYstar( aJet );
    if( aJetYstar < m_varYstarBinning.front() ||
	aJetYstar > m_varYstarBinning.back() )
      { continue; }
      
    // for each A jet, need to find closest B jet
    // set deltaRmin to something large, find jet with smallest
    // deltaRmin less than Rmax and pair that to the a jet
    double           deltaRmin = 4 * constants::ETAMAX;
    TLorentzVector* pairedBJet = NULL;
    for( auto& bJet : vB_jets ){
      // "cleaned" jets had their px, py, pz set to 0
      // to avoid changing the vector size, skip if
      // its pt is zero.
      if( bJet.Pt() == 0 ){ continue; }
      // take the closest bjet to ajet
      double deltaR = anaTool->DeltaR( bJet, aJet );
      if( deltaR <= deltaRmin ) {	
	pairedBJet = &bJet;
	deltaRmin  = deltaR;
      }
    }  // end loop over b jets
    // just for safety. if there is at least
    // one b jet, it should pair to all a jets
    if( !pairedBJet ){ continue; }

    // Require jets to be within some radius
    // If they aren't, dont add a pair.
    if( deltaRmin > m_dRmax ){ continue; }
    
    vA_paired_jets.push_back( aJet );
    vB_paired_jets.push_back( *pairedBJet );
  } // end loop over truth jets
}

double DiJetAnalysisMC::GetSpectWeight( const TLorentzVector& jet ){

  // use just the -4.0 < y*1 < -2.7 histogram
  // return m_vSpectWeightFits[ bin - 1 ]->Eval( jet.Pt()/1000. );
  // if it is 23rd uncertainty, do not reweight
  if( m_uncertComp == 23 ){ return 1; }
  return m_vSpectWeightFits[ 0 ]->Eval( jet.Pt()/1000. );
}

double DiJetAnalysisMC::GetDphiWeight( const TLorentzVector& jet1,
				       const TLorentzVector& jet2 ){

  // if it is 23rd uncertainty, do not reweight
  if( m_uncertComp == 23 ){ return 1; }
  
  double jetYstar1 = GetYstar( jet1 );
  double jetYstar2 = GetYstar( jet2 );
  double dPhi      = anaTool->DeltaPhi( jet1, jet2 );

  // fix this later.
  // weights only given for negative jetYstar1
  // only do this for pp case.
  if( !m_is_pPb && jetYstar1 > 0 ){ jetYstar1 *= -1; jetYstar2 *= -1; }

  int bin = m_dPhiWeight->FindBin( jetYstar1, jetYstar2, dPhi );

  return m_dPhiWeight->GetBinContent( bin );
}

TH1* DiJetAnalysisMC::CombineSamples( std::vector< TH1* >& vSampleHin,
				      const std::string& name ){
  if( !vSampleHin.size() ){ return NULL; }

  TH1* h_res = static_cast< TH1D* >
    ( vSampleHin[0]->Clone( Form("h_%s_%s", name.c_str() , m_allName.c_str() ) ) );
  h_res->Reset();
  
  for( uint iG = 0; iG < vSampleHin.size(); iG++ ){
    double scale = m_vJznEff[iG] * m_vJznSigma[iG] /
      m_vJznSumOverlayWeights[iG];
    h_res->Add( vSampleHin[iG], scale );
  }
  h_res->Scale( 1./m_sumSigmaEff );

  vSampleHin.push_back( h_res );
  
  return h_res;
}

TH2* DiJetAnalysisMC::CombineSamples( std::vector< TH2* >& vSampleHin,
				      const std::string& name ){
  if( !vSampleHin.size() ){ return NULL; }

  TH2* h_res = static_cast< TH2D* >
    ( vSampleHin[0]->Clone( Form("h_%s_%s", name.c_str() , m_allName.c_str() ) ) );
  h_res->Reset();
  
  for( uint iG = 0; iG < vSampleHin.size(); iG++ ){
    double scale = m_vJznEff[iG] * m_vJznSigma[iG] /
      m_vJznSumOverlayWeights[iG];
    h_res->Add( vSampleHin[iG], scale );
  }
  h_res->Scale( 1./m_sumSigmaEff );

  vSampleHin.push_back( h_res );
  
  return h_res;
}

TH3* DiJetAnalysisMC::CombineSamples( std::vector< TH3* >& vSampleHin,
				      const std::string& name ){
  if( !vSampleHin.size() ){ return NULL; }

  TH3* h_res = static_cast< TH3D* >
    ( vSampleHin[0]->Clone( Form("h_%s_%s", name.c_str() , m_allName.c_str() ) ) );
  h_res->Reset();
  
  for( uint iG = 0; iG < vSampleHin.size(); iG++ ){
    double scale = m_vJznEff[iG] * m_vJznSigma[iG] /
      m_vJznSumOverlayWeights[iG];
    h_res->Add( vSampleHin[iG], scale );
  }
  h_res->Scale( 1./m_sumSigmaEff );

  vSampleHin.push_back( h_res );
  
  return h_res;
}

THnSparse* DiJetAnalysisMC::CombineSamples( std::vector< THnSparse* >& vSampleHin,
				      const std::string& name ){
  if( !vSampleHin.size() ){ return NULL; }

  THnSparse* h_res = static_cast< THnSparseD* >
    ( vSampleHin[0]->Clone( Form("h_%s_%s", name.c_str() , m_allName.c_str() ) ) );
  h_res->Reset();
  
  for( uint iG = 0; iG < vSampleHin.size(); iG++ ){
    double scale = m_vJznEff[iG] * m_vJznSigma[iG] /
      m_vJznSumOverlayWeights[iG];
    h_res->Add( vSampleHin[iG], scale );
  }
  h_res->Scale( 1./m_sumSigmaEff );

  vSampleHin.push_back( h_res );
  
  return h_res;
}

void DiJetAnalysisMC::CombineSamples( TH1* h_res,
				      std::vector< TH1* >& vSampleHin,
				      std::vector< TH1* >& vSampleNentIn,
				      const std::string& name ){

  // loop over xBins
  for( int xBin = 1; xBin <= h_res->GetNbinsX(); xBin++ ){

    double valTot      = 0;
    double valErrorTot = 0;
    double denomTot    = 0;

    for( uint iG = 0; iG < vSampleHin.size(); iG++ ){
      if( vSampleHin[iG]->GetEntries() == 0 ){ continue; }
      
      double valueBin    = vSampleHin   [iG]->GetBinContent( xBin );
      double valueBinErr = vSampleHin   [iG]->GetBinError( xBin );
      double nEntriesBin = vSampleNentIn[iG]->GetBinContent( xBin );
      double weight      = m_vJznWeights[iG] / m_vJznSumOverlayWeights[iG];

      if( nEntriesBin < m_nMinEntriesFit || valueBin == 0 ){ continue; }
      
      double val    = weight * nEntriesBin * valueBin;
      double valErr = weight * nEntriesBin * valueBinErr;
      valTot       += val;
      valErrorTot  += valErr * valErr;

      double denom  = weight * nEntriesBin;
      denomTot     += denom;
    }

    double valFinal      = valTot / denomTot;
    double valErrorFinal = valErrorTot / ( denomTot * denomTot );
    valErrorFinal        = std::sqrt( valErrorFinal );

    // check if we have NaN from
    // dividing valTot by zero (denomTot)
    if( std::isnan( valFinal ) ){ continue; }
    
    h_res->SetBinContent( xBin, valFinal );
    h_res->SetBinError  ( xBin, valErrorFinal );
  } // end loop over xBins
}

void DiJetAnalysisMC::CombineSamples( TH1* h_res,
				      std::vector< TH1* >& vSampleHin,
				      std::vector< TH1* >& vSampleNentInT,
				      std::vector< TH1* >& vSampleNentInR,
				      std::vector< TH2* >& vSampleRespMat,
				      const std::string& name ){

  // loop over xBins
  for( int xBin = 1; xBin <= h_res->GetNbinsX(); xBin++ ){

    double valTot      = 0;
    double valErrorTot = 0;
    double denomTot    = 0;

    for( uint iG = 0; iG < vSampleHin.size(); iG++ ){
      if( vSampleHin[iG]->GetEntries() == 0 ){ continue; }

      double valueBin = vSampleHin    [iG]->GetBinContent( xBin );
      double nEntT    = vSampleNentInT[iG]->GetBinContent( xBin );
      double nEntR    = vSampleNentInR[iG]->GetBinContent( xBin );
      double nEntMat  = vSampleRespMat[iG]->GetBinContent( xBin );
      
      double weight   = m_vJznWeights[iG] / m_vJznSumOverlayWeights[iG];

      double newError =
	std::sqrt( std::pow( nEntT, 2 ) / std::pow( nEntR, 3 ) *
		   ( 1 - std::pow( nEntMat, 2 ) / ( nEntT * nEntR ) ) );

      // if one of them is zero, leave it alone
      if( nEntT == 0 || nEntR == 0 ){
	newError = 0;
      }

      double val    = weight * nEntR * valueBin;
      double valErr = weight * nEntR * newError;

      valTot       += val;
      valErrorTot  += valErr * valErr;

      double denom  = weight * nEntR;
      denomTot     += denom;
    }

    double valFinal      = valTot / denomTot;
    double valErrorFinal = valErrorTot / ( denomTot * denomTot );
    valErrorFinal        = std::sqrt( valErrorFinal );
    
    // check if we have NaN from
    // dividing valTot by zero (denomTot)
    if( std::isnan( valFinal ) ){ continue; }
    
    h_res->SetBinContent( xBin, valFinal );
    h_res->SetBinError  ( xBin, valErrorFinal );
  } // end loop over xBins
}

TGraphAsymmErrors* DiJetAnalysisMC::CombineSamples( std::vector< TGraphAsymmErrors* >& vSampleVIN,
						    std::vector< TH1* >& vSampleNentIN ){

  // old code that used maps. in the efficiency part, construct map from vector.
  // no time to rewrite right now. (04.06.18)
  
  // check if we have one
  if( !vSampleNentIN.size() ){ return NULL; }

  // temp histos for this scope
  // fill values from TGraph into these histos
  std::vector< TH1* > vSampleVal;
  std::vector< TH1* > vSampleValErrorLow;
  std::vector< TH1* > vSampleValErrorHigh;

  TH1* hRef = vSampleNentIN[0];
  
  int  nXbins = hRef->GetNbinsX();
  double xMin = hRef->GetXaxis()->GetXmin();
  double xMax = hRef->GetXaxis()->GetXmax();

  uint nSamples = vSampleVIN.size();
  
  for( uint iG = 0; iG < nSamples; iG++ ){
    std::string label = m_vJznLabels[ iG ];
        
    vSampleVal.push_back
      ( new TH1D( Form("h_jznVal_%s", label.c_str() ),
		  Form("h_jznVal_%s", label.c_str() ),
		  nXbins, xMin, xMax ) );
			 
    vSampleValErrorLow.push_back
      (new TH1D( Form("h_jznValErrorLow_%s", label.c_str() ),
		 Form("h_jznValErrorLow_%s", label.c_str() ),
		 nXbins, xMin, xMax ) );

    vSampleValErrorHigh.push_back
      (new TH1D( Form("h_jznValErrorHigh_%s", label.c_str() ),
		 Form("h_jznValErrorHigh_%s", label.c_str() ),
		 nXbins, xMin, xMax ) );
    
    TGraphAsymmErrors* gPts = vSampleVIN[ iG ];
    double x, y;
    
    for( int i = 0; i < gPts->GetN(); i++ ){
      gPts->GetPoint(i, x, y);
      int bin = vSampleVal[ iG ]->FindBin( x );

      double eYlow  = gPts->GetErrorYlow ( i );
      double eYhigh = gPts->GetErrorYhigh( i );
      
      vSampleVal         [ iG ]->SetBinContent( bin, y      );
      vSampleValErrorLow [ iG ]->SetBinContent( bin, eYlow  );
      vSampleValErrorHigh[ iG ]->SetBinContent( bin, eYhigh );
    }    
  }

  std::vector< double > vX;
  std::vector< double > vY;
  std::vector< double > vEx;
  std::vector< double > vEyLow;
  std::vector< double > vEyHigh;
  
  for( int xBin = 1; xBin <= nXbins; xBin++ ){

    double valTot   = 0;
    double denomTot = 0;
    double valErrorLowTot  = 0;
    double valErrorHighTot = 0;

    double xCent  = hRef->GetBinCenter( xBin );
    double xWidth = hRef->GetBinWidth ( xBin );
    
    for( uint iG = 0; iG < nSamples; iG++ ){
      double nEntriesBin = vSampleNentIN   [ iG ]->GetBinContent( xBin );
      double weight      = m_vJznWeights[ iG ] / m_vJznSumOverlayWeights[ iG ];

      double valueBin    = vSampleVal[ iG ]->GetBinContent( xBin );
     
      double valueBinErrLow  = vSampleValErrorLow [ iG ]->GetBinContent( xBin );
      double valueBinErrHigh = vSampleValErrorHigh[ iG ]->GetBinContent( xBin );

      if( valueBin == 0 && valueBinErrLow == 0 &&  valueBinErrHigh == 0 )
	{ continue; }
      
      if( nEntriesBin <  5 )
	{ continue; }
      
      double val         = weight * nEntriesBin * valueBin;
      double valErrLow   = weight * nEntriesBin * valueBinErrLow;
      double valErrHigh  = weight * nEntriesBin * valueBinErrHigh;
      valTot            += val;
      valErrorLowTot    += valErrLow  * valErrLow;
      valErrorHighTot   += valErrHigh * valErrHigh;

      double denom  = weight * nEntriesBin;
      denomTot     += denom;
    }

    double valFinal      = valTot / denomTot;
    
    double valErrorLowFinal  = valErrorLowTot  / ( denomTot * denomTot );
    double valErrorHighFinal = valErrorHighTot / ( denomTot * denomTot );

    valErrorLowFinal  = std::sqrt( valErrorLowFinal );
    valErrorHighFinal = std::sqrt( valErrorHighFinal );
    
    // check if we have NaN from
    // dividing valTot by zero (denomTot)
    if( std::isnan( valFinal ) ){ continue; }

    vX.push_back ( xCent  ); vY.push_back( valFinal );
    vEx.push_back( xWidth * 0.5 );
    vEyLow.push_back ( valErrorLowFinal  );
    vEyHigh.push_back( valErrorHighFinal );
  }

  for( auto& h : vSampleVal          ){ delete h; }
  for( auto& h : vSampleValErrorLow  ){ delete h; }
  for( auto& h : vSampleValErrorHigh ){ delete h; }
    
  return new TGraphAsymmErrors( vX.size(), &vX[0] , &vY[0],&vEx[0], &vEx[0],
				&vEyLow[0], &vEyHigh[0] );
}

void DiJetAnalysisMC::SetCfactorsErrors( TH1* hR, TH1* hT, TH2* hM, TH1* hC ){

  std::cout << hT->GetName() << " +++ " << hT->GetNbinsX() << std::endl;
  // set factors
  for( int xBin = 1; xBin <= hT->GetNbinsX(); xBin++ ){
    int cBin = hC->FindBin( hT->GetBinCenter(xBin) );

    // if one of them is zero, leave it alone
    if( hT->GetBinContent( xBin ) == 0 ||
	hR->GetBinContent( xBin ) == 0 ){
      hC->SetBinContent( cBin, 1 );
      hC->SetBinError  ( cBin, 0 );
      continue;
    }
	      		
    double vR = hR->GetBinContent( xBin  );
    double vT = hT->GetBinContent( xBin  );

    double vM = hM->GetBinContent( cBin, cBin );
	     
    // error on correction factor    
    double newDphiError =
      std::sqrt( std::pow( vT, 2 ) / std::pow( vR, 3 ) *
		 ( 1 - std::pow( vM, 2 ) / ( vT * vR ) ) );

    std::cout << xBin << "   T = " << vT << "   R = " << vR << "   M = "
	      << vM << "    C = " << hC->GetBinContent( xBin )
	      << "   err = " << newDphiError << std::endl;
    
    hC->SetBinError( cBin, newDphiError );
  }
}

double DiJetAnalysisMC::GetJetWeight( const TLorentzVector& jet ){

  // !!ONLY FOR pPb!!
  
  if( m_is_pPb ){
    if( m_FCalEt > 0.09 ){
      return h_pPbFCalWeights->GetBinContent
	( h_pPbFCalWeights->FindBin( 0.09 ) );
    }
    return h_pPbFCalWeights->GetBinContent
      ( h_pPbFCalWeights->FindBin( m_FCalEt ) );
  }
  
  return 1;
}

void DiJetAnalysisMC::GetTypeTitle( const std::string& type,
				    std::string& yTitleMean,
				    std::string& yTitleSigma ){ 
  if( type.find("Rpt") != std::string::npos ){
    yTitleMean  = "<#it{p}_{T}^{reco}/#it{p}_{T}^{truth}>";
    yTitleSigma = "#sigma(#it{p}_{T}^{reco}/#it{p}_{T}^{truth})"; 
  } else if( type.find("Deta") != std::string::npos ){
    yTitleMean  = "<#Delta#eta>";
    yTitleSigma = "#sigma(#Delta#eta)"; 
  } else if( type.find("Dphi") != std::string::npos ){
    yTitleMean  = "<#Delta#phi>";
    yTitleSigma = "#sigma(#Delta#phi)"; 
  } 
}

void DiJetAnalysisMC::GetSpectWeightInfo( std::string& name_a , std::string& name_b,
					  std::string& fName_a, std::string& fName_b ){

  
  // a is Data, b is MC
  name_a   = m_ystarSpectFineName;
  name_b   = m_ystarSpectFineRecoName;
  fName_a  = m_fNamePerfWeightData;
  fName_b  = m_fNamePerf;
}

void DiJetAnalysisMC::GetDphiWeightInfo ( std::string& name_a , std::string& name_b,
					  std::string& fName_a, std::string& fName_b ){
 
  // a is Data, b is MC
  name_a   = m_dPhiName;
  name_b   = m_dPhiRecoName;
  fName_a  = m_fNamePhysWeightData;
  fName_b  = m_fNamePhys;
}

void DiJetAnalysisMC::GetSpectUnfoldingInfo( std::string& measuredName,
					     std::string& recoName,
					     std::string& truthName,
					     std::string& respMatName,
					     std::string& unfoldedLabel,
					     std::string& typeLabel ){

  measuredName  = m_ystarSpectRecoName;
  recoName      = m_ystarSpectRecoName;
  truthName     = m_ystarSpectTruthName;
  respMatName   = m_ystarSpectRespMatName;
  unfoldedLabel = "dN/d#it{p}_{T} [GeV]";;
  typeLabel     = "MC";
}

void DiJetAnalysisMC::GetDphiUnfoldingInfo( std::string& measuredName,
					    std::string& truthName,
					    std::string& unfoldedLabel,
					    std::string& typeLabel ){

  measuredName  = m_dPhiRecoName;
  truthName     = m_dPhiTruthName;
  unfoldedLabel = m_sDphi;
  typeLabel     = "MC";
}

void DiJetAnalysisMC::GetInfoTogether( std::string& name_a , std::string& name_b ,
				       std::string& label_a, std::string& label_b,
				       std::string& fName_a, std::string& fName_b,
				       int option ){

  std::string* pFname = NULL;
  
  switch( option ){
  case 0:
    pFname = &m_fNamePerfUF;
    name_a = m_ystarSpectName;
    name_b = m_ystarSpectName;
    break;
  case 1:
    pFname = &m_fNamePhysUF;
    name_a = m_dPhiName;
    name_b = m_dPhiName;
    break;
  }
  
  int combinationBoth = GetConfig()->GetValue( "combinationBoth", 0 );

  if( combinationBoth == 0 ){
    name_a   +=  "_" + m_recoName;
    name_b   +=  "_" + m_truthName;
    label_a  = "Reco_{MC}";
    label_b  = "Truth_{MC}";
    fName_a  = *pFname;
    fName_b  = *pFname;
  } else if ( combinationBoth == 1 ){
    name_a   +=  "_" + m_truthName;
    name_b   +=  "_" + m_recoName + "_" + m_unfoldedName;
    label_a  = "Truth_{MC}";
    label_b  = "UF_{MC}";
    fName_a  = *pFname;
    fName_b  = *pFname;
  } else if ( combinationBoth == 2 ){
    name_a   +=  "_" + m_recoName + "_" + m_unfoldedName;
    name_b   +=  "_" + m_recoName + "_" + m_unfoldedName;
    label_a  = "#it{p}+Pb";
    label_b  = "#it{pp}";
    m_is_pPb = true;  DiJetAnalysis::Initialize();
    fName_a  = *pFname;
    m_is_pPb = false; DiJetAnalysis::Initialize();
    fName_b  = *pFname;
  } else if ( combinationBoth == 3 ){
    name_a   +=  "_" + m_truthName;
    name_b   +=  "_" + m_truthName;
    label_a  = "#it{p}+Pb";
    label_b  = "#it{pp}";
    m_is_pPb = true;  DiJetAnalysis::Initialize();
    fName_a  = *pFname;
    m_is_pPb = false; DiJetAnalysis::Initialize();
    fName_b  = *pFname;
  } else if ( combinationBoth == 4 ){
    name_a   +=  "_" + m_truthName;
    name_b   +=   "_" + m_recoName + "_" + m_unfoldedName;
    label_a  = "Truth_{MC}";
    label_b  = "Unfolded{MC}";
    m_is_pPb = false;  DiJetAnalysis::Initialize();
    fName_a  = *pFname + ".1";
    m_is_pPb = false; DiJetAnalysis::Initialize();
    fName_b  = *pFname + ".2";
  } 
}

void DiJetAnalysisMC::GetPurityEff( TH2* hRespMat, TH1* hPurity, TH1* hEff ){

  for( int xBin = 1 ; xBin <= hRespMat->GetNbinsX(); xBin++ ){
    double sum = 0;
    for( int yBin = 1 ; yBin <= hRespMat->GetNbinsY(); yBin++ ){
      // if( xBin == yBin ){ continue; }
      sum += hRespMat->GetBinContent( xBin, yBin );
    }
    double eff = sum ? hRespMat->GetBinContent( xBin, xBin )/sum : 0 ;    
    hEff->SetBinContent( xBin, eff );
  }

  for( int yBin = 1 ; yBin <= hRespMat->GetNbinsY(); yBin++ ){
    double sum = 0;
    for( int xBin = 1 ; xBin <= hRespMat->GetNbinsX(); xBin++ ){
      // if( xBin == yBin ){ continue; }
      sum += hRespMat->GetBinContent( xBin, yBin );
    }
    double purity = sum ? hRespMat->GetBinContent( yBin, yBin )/sum : 0 ;    
    hPurity->SetBinContent( yBin, purity );
  }
}

//---------------------------------
//     Get Quantities / Plot 
//---------------------------------
void DiJetAnalysisMC::LoadHistograms( int opt ){

  std::string fNameRaw = "";
  switch( opt ){
  case 0:
    fNameRaw = m_fNameRaw;
    break;
  case 1:
    fNameRaw = m_fNameRawUW;
    break;
  default:
    fNameRaw = m_fNameRaw;
  }
  
  TFile* fIn = TFile::Open( fNameRaw.c_str() ); 
  
  for( uint iG = 0; iG < m_nJzn; iG++){
    std::string jzn = m_vJznLabels[iG];
    
    // -------- maps ---------
    m_vHjznEtaPhiMap.push_back
      ( static_cast< TH2D* >
	( fIn->Get( Form("h_etaPhiMap_%s", jzn.c_str() ))));
    m_vHjznEtaPhiMap.back()->SetDirectory(0);

    m_vHjznEtaPtMap.push_back
      ( static_cast< TH2D* >
	( fIn->Get( Form("h_etaPtMap_%s", jzn.c_str() ))));
    m_vHjznEtaPtMap.back()->SetDirectory(0);

    // -------- spect --------
    m_vHjznYstarSpectReco.push_back
      ( static_cast< TH2D* >
	( fIn->Get
	  ( Form("h_%s_%s", m_ystarSpectRecoName.c_str(), jzn.c_str() ))));
    m_vHjznYstarSpectReco.back()->SetDirectory(0);

    m_vHjznYstarSpectTruth.push_back 
      ( static_cast< TH2D* >
	( fIn->Get( Form("h_%s_%s", m_ystarSpectTruthName.c_str(), jzn.c_str() ))));
    m_vHjznYstarSpectTruth.back()->SetDirectory(0);

    m_vHjznYstarSpectFineReco.push_back
      ( static_cast< TH2D* >
	( fIn->Get( Form("h_%s_%s", m_ystarSpectFineRecoName.c_str(), jzn.c_str() ))));
    m_vHjznYstarSpectFineReco.back()->SetDirectory(0);

    m_vHjznYstarSpectFineTruth.push_back 
      ( static_cast< TH2D* >
	( fIn->Get( Form("h_%s_%s", m_ystarSpectFineTruthName.c_str(), jzn.c_str() ))));
    m_vHjznYstarSpectFineTruth.back()->SetDirectory(0);

    m_vHjznYstarSpectFineTruthUP.push_back 
      ( static_cast< TH2D* >
	( fIn->Get( Form("h_%s_%s", m_ystarSpectFineTruthUPName.c_str(), jzn.c_str() ))));
    m_vHjznYstarSpectFineTruthUP.back()->SetDirectory(0);

    // --- spectra response matrix ----
    m_vHjznYstarSpectRespMat.push_back 
      ( static_cast< TH3D* >
	( fIn->Get( Form("h_%s_%s", m_ystarSpectRespMatName.c_str(), jzn.c_str() ))));
    m_vHjznYstarSpectRespMat.back()->SetDirectory(0);

    // --- ystar response matrix ----
    m_vHjznYstarRespMat.push_back 
      ( static_cast< TH3D* >
	( fIn->Get( Form("h_%s_%s", m_ystarRespMatName.c_str(), jzn.c_str() ))));
    m_vHjznYstarRespMat.back()->SetDirectory(0);

    // --------- recoTruthRpt ---------
    m_vHjznRecoTruthRpt.push_back
      ( static_cast< TH3D* >
	( fIn->Get( Form("h_recoTruthRpt_%s", jzn.c_str() ))));
    m_vHjznRecoTruthRpt.back()->SetDirectory(0);

    m_vHjznRecoTruthRptNent.push_back
      ( static_cast< TH2D* >
	( fIn->Get( Form("h_recoTruthRptNent_%s", jzn.c_str() ))));
    m_vHjznRecoTruthRptNent.back()->SetDirectory(0);
    
    // --------- recoTruthDeta ---------
    m_vHjznRecoTruthDeta.push_back
      ( static_cast< TH3D* >
	( fIn->Get( Form("h_recoTruthDeta_%s", jzn.c_str() ))));
    m_vHjznRecoTruthDeta.back()->SetDirectory(0);

    m_vHjznRecoTruthDetaNent.push_back
      ( static_cast< TH2D* >
	( fIn->Get( Form("h_recoTruthDetaNent_%s", jzn.c_str() ))));
    m_vHjznRecoTruthDetaNent.back()->SetDirectory(0);
    
    // --------- recoTruthDphi ---------
    m_vHjznRecoTruthDphi.push_back
      ( static_cast< TH3D* >
	( fIn->Get( Form("h_recoTruthDphi_%s", jzn.c_str() )))); 
    m_vHjznRecoTruthDphi.back()->SetDirectory(0);

    m_vHjznRecoTruthDphiNent.push_back
      ( static_cast< TH2D* >
	( fIn->Get( Form("h_recoTruthDphiNent_%s", jzn.c_str() ))));
    m_vHjznRecoTruthDphiNent.back()->SetDirectory(0);

    // ----- rtrk ----
    if( m_is_pPb && doRtrk ){
      m_vHjznRtrk1.push_back
	( static_cast< TH3D* >
	  ( fIn->Get( Form("h_%s1_%s", m_rtrkName.c_str(), jzn.c_str() ))));
      m_vHjznRtrk1.back()->SetDirectory(0);

      m_vHjznRtrk2.push_back
	( static_cast< TH3D* >
	  ( fIn->Get( Form("h_%s2_%s", m_rtrkName.c_str(), jzn.c_str() ))));
      m_vHjznRtrk2.back()->SetDirectory(0);

      m_vHjznRtrk4.push_back
	( static_cast< TH3D* >
	  ( fIn->Get( Form("h_%s4_%s", m_rtrkName.c_str(), jzn.c_str() ))));
      m_vHjznRtrk4.back()->SetDirectory(0);
    }
    
    // -------- dPhi- --------
    m_vHjznDphiReco.push_back
      ( static_cast< THnSparse *>
	( fIn->Get( Form("h_%s_%s", m_dPhiRecoName.c_str(), jzn.c_str() ))));  
    	
    m_vHjznDphiTruth.push_back
      ( static_cast< THnSparse *>
	( fIn->Get( Form("h_%s_%s", m_dPhiTruthName.c_str(), jzn.c_str() ))));  
    
    // -------- Dphi Response Matrix --------    
    m_vHjznDphiRespMat.push_back
      ( static_cast< THnSparse *>
	( fIn->Get( Form("h_%s_%s", m_dPhiRespMatName.c_str(), jzn.c_str() ))));  
      
    // -------- Dphi Response Matrix --------    
    m_vHjznDphiRespMatReb.push_back
      ( static_cast< THnSparse *>
	( fIn->Get( Form("h_%s_%s", m_dPhiRespMatRebName.c_str(), jzn.c_str() ))));  

  }
  
  fIn->Close(); delete fIn;
}

TH2* DiJetAnalysisMC::MakeSpectWeights( TFile* fOut ){

  std::vector< TH1* > vH;
  std::vector< TH1* > vR;
  std::vector< TF1* > vF;

  std::string name_d, name_mc, fName_d, fName_mc;
  GetSpectWeightInfo( name_d, name_mc, fName_d, fName_mc );

  std::cout << fName_d  << std::endl;
  std::cout << fName_mc << std::endl;

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
  GetSpectraLabels( axisLabel, axisLabelTex, name_d );

  TFile* fIn_d  = TFile::Open( fName_d .c_str() );
  TFile* fIn_mc = TFile::Open( fName_mc.c_str() );

  fOut->cd();

  std::string hName = "h_" + m_spectName + "_" + m_sWeights + "_" + m_allName;

  TH2* hAll =
    new TH2D( hName.c_str(), ";#it{p}_{T1};#eta_{1}",
	      m_nVarYstarBins, 0, 1, m_nPtSpectBins, m_ptSpectMin, m_ptSpectMax );
  hAll->GetXaxis()->Set( m_nVarYstarBins , &( m_varYstarBinning [0] ) );
  
  TCanvas cAll( "cAll", "cAll", 800, 600 );

  double y0 = 0.6, y1 = 0.8;
  if( m_is_pPb ){ y0 = 0.3; y1 = 0.4; }
  
  TLegend leg( 0.35, y0, 0.7, y1 );

  styleTool->SetLegendStyle( &leg );
  
  for( uint xBin = 1; xBin <= m_nVarYstarBins; xBin++ ){
    
    double xLow  = m_varYstarBinning[ xBin - 1 ];
    double xUp   = m_varYstarBinning[ xBin ];
   
    std::string hTag = anaTool->GetName( xLow, xUp, axisLabel);

    std::string hName_d  = "h_" + name_d  + "_" + m_allName + "_" + hTag;
    std::string hName_mc = "h_" + name_mc + "_" + m_allName + "_" + hTag;
    
    TH1* h_d  = static_cast< TH1D* >( fIn_d ->Get( hName_d. c_str() ) );
    TH1* h_mc = static_cast< TH1D* >( fIn_mc->Get( hName_mc.c_str() ) ); 
    vH.push_back( h_d ); vH.push_back( h_mc );

    int pTbinLow = h_d->FindBin( m_varPtBinning.front() ) + 1;
    int pTbinUp  = h_d->FindBin( m_varPtBinning.back () ) - 1;
    
    double integralD  = h_d ->Integral( pTbinLow, pTbinUp );
    double integralMC = h_mc->Integral( pTbinLow, pTbinUp );

    double scalingFactor = integralD / integralMC;

    h_mc->Scale( scalingFactor );

    h_d ->Rebin(2);
    h_mc->Rebin(2);
    
    h_d ->Write();
    h_mc->Write();
    
    TH1* hR = static_cast< TH1D* >
      ( h_d->Clone( Form( "h_%s_%s_%s_%s", m_spectName.c_str(), m_sWeights.c_str(),
			  m_allName.c_str(), hTag.c_str() ) ) );
    styleTool->SetHStyleRatio( hR, 0 );
    vR.push_back( hR );

    hR->Divide( h_mc );

    for( int yBin = 1; yBin <= hAll->GetYaxis()->GetNbins(); yBin++ ){
      double val      = hR->GetBinContent( yBin );
      double valError = hR->GetBinError  ( yBin );
      hAll->SetBinContent( xBin, yBin, val      );
      hAll->SetBinError  ( xBin, yBin, valError );
    }


    double fitMin = m_varPtBinning.front() + 1;
    double fitMax = m_varPtBinning.back () - 1;
    
    TF1* fitR = anaTool->FitPol2( hR, fitMin, fitMax );
    styleTool->SetHStyle( fitR, 1 );
    vF.push_back( fitR );
    
    TCanvas c( "c", "c", 800, 600 );

    hR->SetYTitle( "Data/MC" );

    hR->GetXaxis()->SetRangeUser( m_varPtBinning.front(), m_varPtBinning.back() );
    
    hR  ->Draw();
    fitR->Draw("same");
    
    line.Draw();
    lineP15.Draw();
    lineN15.Draw();

    DrawAtlasRight();
    
    drawTool->DrawRightLatex
      ( 0.45, 0.85, anaTool->GetLabel( xLow, xUp, axisLabelTex ) );
    
    SaveAsAll( c, hR->GetName() );
    
    hR  ->Write();
    fitR->Write();

    // in pPb skip the other bins.
    if( m_is_pPb && xBin != 1 ){ continue; }

    cAll.cd();
    styleTool->SetHStyle( hR, xBin - 1 );
    hR->Draw( "ep X0 same" );
    leg.AddEntry( hR,  anaTool->GetLabel( xLow, xUp, axisLabelTex ).c_str() ); 
  }

  hAll->Write();

  cAll.cd();

  leg.Draw();

  line.Draw();
  lineP15.Draw();
  lineN15.Draw();

  DrawAtlasRight();

  SaveAsAll( cAll, hName.c_str() );

  // for( auto& f : vF ){ delete f; }
  for( auto& r : vR ){ delete r; }
  for( auto& h : vH ){ delete h; }
  
  return hAll;
}

TH3* DiJetAnalysisMC::MakeDphiWeights( TFile* fOut ){
  
  std::vector< TF1* > vF;
  std::vector< TH1* > vH;
  std::vector< TH1* > vR;

  std::string name_d, name_mc, fName_d, fName_mc;
  GetDphiWeightInfo( name_d, name_mc, fName_d, fName_mc );

  TAxis* y1Axis  = m_dPP->GetDefaultTAxis( 0 );
  TAxis* y2Axis  = m_dPP->GetDefaultTAxis( 1 ); int y2AxisBins  = y2Axis->GetNbins();
  TAxis* pt1Axis = m_dPP->GetDefaultTAxis( 2 ); int pt1AxisBins = pt1Axis->GetNbins();
  TAxis* pt2Axis = m_dPP->GetDefaultTAxis( 3 ); int pt2AxisBins = pt2Axis->GetNbins();

  // make a TH3 to fill with unfolded results.
  TH3* hnAll = new TH3D
    ( Form( "h_%s_%s_%s", m_dPhiName.c_str(), m_sWeights.c_str(), m_allName.c_str() ), "",
      y1Axis->GetNbins(), y1Axis->GetXbins()->GetArray(),
      y2Axis->GetNbins(), y2Axis->GetXbins()->GetArray(),
      m_nVarDphiBins, & m_varDphiBinning[0]);

  hnAll->GetXaxis()->SetTitle( m_dPP->GetDefaultAxisLabel(0).c_str() );
  hnAll->GetYaxis()->SetTitle( m_dPP->GetDefaultAxisLabel(1).c_str() );
  hnAll->GetZaxis()->SetTitle( m_sDphi.c_str() );

  TAxis* axis0 = m_dPP->GetTAxis(0); int nAxis0Bins = axis0->GetNbins();
  TAxis* axis1 = m_dPP->GetTAxis(1); int nAxis1Bins = axis1->GetNbins();
  TAxis* axis2 = m_dPP->GetTAxis(2); int nAxis2Bins = axis2->GetNbins();
  TAxis* axis3 = m_dPP->GetTAxis(3); int nAxis3Bins = axis3->GetNbins();

  // lines to be drawn along axis3. this is
  // x-axis that widths are plotted as function of 

  TFile* fIn_d  = TFile::Open( fName_d .c_str() );
  TFile* fIn_mc = TFile::Open( fName_mc.c_str() );

  double xMin = m_dPhiZoomLow;
  double xMax = m_dPhiZoomHigh;
  
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

  fOut->cd();

  std::vector< TCanvas* > vC;
  std::vector< std::vector< TH1* > > vData;
  std::vector< std::vector< TH1* > > vMC;
  vData.resize( y2AxisBins );
  vMC  .resize( y2AxisBins );
  
  const int y1Bin = 1;
  for( int y2Bin = 1; y2Bin <= y2AxisBins; y2Bin++ ){
    double y1Low, y1Up;
    anaTool->GetBinRange
      ( y1Axis, y1Bin, y1Bin, y1Low, y1Up );
    double y2Low, y2Up;
    anaTool->GetBinRange
      ( y2Axis, y2Bin, y2Bin, y2Low, y2Up );

    std::string cName = "cw_" +
      anaTool->GetName( y1Low, y1Up, m_dPP->GetDefaultAxisName(0) ) + "_" +
      anaTool->GetName( y2Low, y2Up, m_dPP->GetDefaultAxisName(1) );
    
    vC.push_back( new TCanvas( cName.c_str(), cName.c_str(), 800, 600 ) );
  }

  TLegend leg( 0.18, 0.22, 0.5, 0.42 );
  styleTool->SetLegendStyle( &leg, 0.75 );

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
      
      for( int axis2Bin = 1; axis2Bin <= nAxis2Bins; axis2Bin++ ){
	double axis2Low , axis2Up;
	anaTool->GetBinRange
	  ( axis2, axis2Bin, axis2Bin, axis2Low, axis2Up );

	std::string legLabel =
	  anaTool->GetLabel( axis1Low, axis2Up, m_dPP->GetDefaultAxisLabel(0) ) + " " +
	  anaTool->GetLabel( axis2Low, axis2Up, m_dPP->GetDefaultAxisLabel(0) );

	for( int axis3Bin = 1; axis3Bin <= nAxis3Bins; axis3Bin++ ){
	  // check we are in correct ystar and pt bins
	  if( !m_dPP->CorrectPhaseSpace
	      ( std::vector<int>{ axis0Bin, axis1Bin, axis2Bin, axis3Bin } ) )
	    { continue; }

	  double axis3Low , axis3Up;
	  anaTool->GetBinRange
	    ( axis3, axis3Bin, axis3Bin, axis3Low, axis3Up );

	  std::string hTagDphi =
	    Form("%s_%s_%s_%s",
		 anaTool->GetName( axis0Low, axis0Up, m_dPP->GetAxisName(0) ).c_str(),
		 anaTool->GetName( axis1Low, axis1Up, m_dPP->GetAxisName(1) ).c_str(),
		 anaTool->GetName( axis2Low, axis2Up, m_dPP->GetAxisName(2) ).c_str(),
		 anaTool->GetName( axis3Low, axis3Up, m_dPP->GetAxisName(3) ).c_str() );

	  std::string hName_d  = "h_" + name_d  + "_" + m_allName + "_" + hTagDphi;
	  std::string hName_mc = "h_" + name_mc + "_" + m_allName + "_" + hTagDphi;
	  std::string hName    = "h_" + m_dPhiName + "_both_" + m_allName + "_" + hTagDphi;
	  
	  TH1* h_d  = static_cast<TH1D*>( fIn_d ->Get( hName_d .c_str() ) );
	  TH1* h_mc = static_cast<TH1D*>( fIn_mc->Get( hName_mc.c_str() ) );
	  vH.push_back( h_d ); vH.push_back( h_mc );
    
	  double integralD  = h_d ->Integral();
	  double integralMC = h_mc->Integral();

	  double scalingFactor = integralD / integralMC;

	  h_mc->Scale( scalingFactor );
	  h_d ->Write();
	  h_mc->Write();

	  if( axis1Bin == 1 ){
	    TCanvas c( "c", "c", 800, 600 );
	    c.SetLogy();
	    
	    h_mc->SetMarkerColor( kRed );
	    h_mc->SetLineColor  ( kRed );

	    h_d ->SetMarkerSize( h_d ->GetMarkerSize() * 1.5 );
	    h_mc->SetMarkerSize( h_mc->GetMarkerSize() * 1.5 );

	    h_d ->Draw( "ep X0" );
	    h_mc->Draw( "ep X0 same" );

	    TLegend leg( 0.7, 0.3, 0.85, 0.4 );
	    styleTool->SetLegendStyle( &leg );
	    leg.AddEntry( h_d , "Data" );
	    leg.AddEntry( h_mc, "MC_{reco}" );

	    leg.Draw();
	    
	    DrawAtlasRight();
	  
	    DrawTopLeftLabels
	      ( m_dPP, axis0Low, axis0Up, axis1Low, axis1Up,
		axis2Low, axis2Up, axis3Low, axis3Up );

	    SaveAsAll( c, hName );
	  }

	  vData[ axis3Bin - 1 ].push_back( h_d  );
	  vMC  [ axis3Bin - 1 ].push_back( h_mc );
	  
	  TH1* hR = static_cast< TH1D* >
	    ( h_d->Clone( Form( "h_%s_%s_%s_%s", m_dPhiName.c_str(), m_sWeights.c_str(),
				m_allName.c_str(), hTagDphi.c_str() ) ) );
	  styleTool->SetHStyleRatio( hR );
	  hR->Divide( h_mc );
	  vR.push_back( hR );

	  hR->SetYTitle( "Data/MC" );
	  // hR->Smooth( 1, "R" );
	  
	  TCanvas c1( "c1", "c1", 800, 600 );
	  
	  hR->GetXaxis()->SetRange( m_dPhiZoomLowBin, m_dPhiZoomHighBin );
	  
	  hR->Draw( "ep X0" );
	  line.Draw();
	  lineP25.Draw();
	  lineN25.Draw();
	  lineP50.Draw();
	  lineN50.Draw();

	  DrawAtlasRight();
	  
	  DrawTopLeftLabels
	    ( m_dPP, axis0Low, axis0Up, axis1Low, axis1Up,
	      axis2Low, axis2Up, axis3Low, axis3Up );

	  SaveAsAll( c1, hR->GetName() );
	  
	  hR->Write();

	  vC[ axis3Bin - 1 ]->cd();
	  styleTool->SetHStyle( hR, ( axis1Bin - 1) * nAxis2Bins + axis2Bin - 1 );
	  if( axis3Bin == 1 ){
	    std::string legend =
	      anaTool->GetLabel( axis1Low, axis1Up, m_dPP->GetAxisLabel(1) ) + " , " + 
	      anaTool->GetLabel( axis2Low, axis2Up, m_dPP->GetAxisLabel(2) );
	    leg.AddEntry( hR, legend.c_str() );
	  }
	  hR->Draw("ep same X0");
	}  // end loop over axis3
      }  // end loop over axis2
    }  // end loop over axis1
  } // end loop over axis0
  
  // set pt axis range to all.
  // sum up over the pt.
  pt1Axis->SetRange( 1, pt1AxisBins );
  pt2Axis->SetRange( 1, pt2AxisBins );

  std::string hName = "h_" + m_dPhiName + "_" + m_sWeights + "_" + m_allName;
  
  TCanvas cAll( "cAll", "cAll", 800, 600 );
  TLegend legAll( 0.6, 0.5, 0.8, 0.75 );
  styleTool->SetLegendStyle( &legAll );
  
  for( int y2Bin = 1; y2Bin <= y2Axis->GetNbins(); y2Bin++ ){
        
    double y1Low, y1Up;
    anaTool->GetBinRange
      ( y1Axis, y1Bin, y1Bin, y1Low, y1Up );
    double y2Low, y2Up;
    anaTool->GetBinRange
      ( y2Axis, y2Bin, y2Bin, y2Low, y2Up );

    std::string hTag = 
      anaTool->GetName( y1Low, y1Up, m_dPP->GetDefaultAxisName(0) ) + "_" +
      anaTool->GetName( y2Low, y2Up, m_dPP->GetDefaultAxisName(1) );

    if( !vData[ y2Bin - 1 ].size() && !vMC[ y2Bin - 1 ].size() )
      { break; }

    std::string hName_d  = "h_" + name_d  + "_" + m_allName + "_" + hTag;
    std::string hName_mc = "h_" + name_mc + "_" + m_allName + "_" + hTag;
	 
    TH1* hData = static_cast< TH1D* >
      ( vData[ y2Bin - 1 ].front()->Clone( hName_d .c_str() ) );
    TH1* hMC   = static_cast< TH1D* >
      ( vMC  [ y2Bin - 1 ].front()->Clone( hName_mc.c_str() ) );

    for( uint i = 1; i < vData[ y2Bin - 1 ].size(); i++ ){
      hData->Add( vData[ y2Bin - 1 ][i] );
      hMC  ->Add( vMC  [ y2Bin - 1 ][i] );
    }

    double integralD  = hData->Integral();
    double integralMC = hMC  ->Integral();

    double scalingFactor = integralD / integralMC;
    hMC->Scale( scalingFactor );

    hData->Write();
    hMC  ->Write();
    
    TH1* hR = static_cast< TH1D* >
      ( hData->Clone( Form( "h_%s_%s_%s_%s", m_dPhiName.c_str(), m_sWeights.c_str(),
			    m_allName.c_str(), hTag.c_str() ) ) );
    styleTool->SetHStyle( hR, 0 );
    hR->Divide( hMC );
    vR.push_back( hR ) ;    

    hR->Smooth();

    // make the 3D histogram with all the weights
    for( int dPhiBin = 1; dPhiBin <= hR->GetNbinsX(); dPhiBin++ ){
      for( uint y1BinTmp = 1; y1BinTmp <= m_nVarYstarBins; y1BinTmp++ ){
	hnAll->SetBinContent( y1BinTmp, y2Bin, dPhiBin, hR->GetBinContent( dPhiBin ) );
	hnAll->SetBinError  ( y1BinTmp, y2Bin, dPhiBin, hR->GetBinError  ( dPhiBin ) );
      }
    }

    TCanvas cc( "cc", "cc", 800, 600 );
    hR->Draw( "ep X0" );
    
    DrawAtlasRight();

    drawTool->DrawLeftLatex
      ( 0.18, 0.86, anaTool->GetLabel
	( y1Low, y1Up, m_dPP->GetDefaultAxisLabel(0) ) );
    drawTool->DrawLeftLatex
      ( 0.18, 0.79, anaTool->GetLabel
	( y2Low, y2Up, m_dPP->GetDefaultAxisLabel(1) ) );  
    DrawAtlasRight();

    SaveAsAll( cc, hR->GetName() );
    
    TCanvas* c = vC[ y2Bin - 1 ];
    c->cd();
    
    line.Draw();
    lineP25.Draw();
    lineN25.Draw();
    lineP50.Draw();
    lineN50.Draw();

    leg.Draw();
    
    DrawAtlasRight();

    drawTool->DrawLeftLatex
      ( 0.18, 0.86, anaTool->GetLabel
	( y1Low, y1Up, m_dPP->GetDefaultAxisLabel(0) ) );
    drawTool->DrawLeftLatex
      ( 0.18, 0.79, anaTool->GetLabel
	( y2Low, y2Up, m_dPP->GetDefaultAxisLabel(1) ) );  
    DrawAtlasRight();

    SaveAsAll( *c, c->GetName() );

    hR->SetYTitle( "Weight" );
    
    cAll.cd();
    styleTool->SetHStyle( hR, y2Bin - 1 );
    hR->Draw( "ep X0 same" );
    hR->GetXaxis()->SetRangeUser( m_dPhiZoomLow, m_dPhiZoomHigh );
    if( m_is_pPb ){
      hR->SetMinimum( 0.77 );
    } else {
      hR->SetMinimum( 0.79 );
    }
    legAll.AddEntry( hR, anaTool->GetLabel
		  ( y2Low, y2Up, m_dPP->GetDefaultAxisLabel(1) ).c_str() ); 

    hR->Write();
  }
  
  hnAll->Write();

  cAll.cd();

  legAll.Draw();
  
  DrawAtlasRight();

  // hardcoded...
  drawTool->DrawLeftLatex( 0.42, 0.78, "2.7<#it{y}_{1}*<4.0");
  
  SaveAsAll( cAll, hName.c_str() );
  
  for( auto& f : vF ){ delete f; }
  for( auto& r : vR ){ delete r; }
  for( auto& h : vH ){ delete h; }

  return hnAll;
}

void DiJetAnalysisMC::MakeScaleRes( std::vector< TH3* >& vJznHin,
				    std::vector< TH2* >& vJznNentIn,
				    const std::string& type ){

  // check if we have one
  if( !vJznHin.size() ){ return; }

  // is same for all the other histos
  TH3* hRef = vJznHin.front();
  
  int  nBinsX = hRef->GetNbinsX();
  double xMin = hRef->GetXaxis()->GetXmin();
  double xMax = hRef->GetXaxis()->GetXmax();

  int  nBinsY = hRef->GetNbinsY();
  double yMin = hRef->GetYaxis()->GetXmin();
  double yMax = hRef->GetYaxis()->GetXmax();

  // to not make spelling errors later
  std::string sMean  = "mean";
  std::string sSigma = "sigma";

  std::string yTitleMean;
  std::string yTitleSigma;
  GetTypeTitle( type, yTitleMean, yTitleSigma );

  std::vector< std::vector< std::vector< TH1* > > > vProj;
  std::vector< std::vector< std::vector< TF1* > > > vFit;
  std::vector< std::vector< TH1* > > vMeans;
  std::vector< std::vector< TH1* > > vSigmas;
  std::vector< std::vector< TH1* > > vNent;
  vProj  .resize( m_nJzn );
  vFit   .resize( m_nJzn );
  vMeans .resize( m_nJzn );
  vSigmas.resize( m_nJzn );
  vNent  .resize( m_nJzn );
  
  // Get The mean, sigma histograms ready.
  for( uint iG = 0; iG < m_nJzn; iG++ ){
    std::string jzn = m_vJznLabels[iG];
    TH3* hJznHin    = vJznHin[iG];
   
    for( int xBin = 1; xBin <= nBinsX; xBin++ ){
      double xLow, xUp;
      anaTool->GetBinRange
	( hJznHin->GetXaxis(), xBin, xBin, xLow, xUp );
      
      //std::string xAxisTitle = hJznHin->GetYaxis()->GetTitle();
      std::string xAxisTitle = "#it{p}_{T}^{truth} [GeV]";

      std::string angleLabel = m_scaleResUseEta ?
	anaTool->GetEtaLabel  ( xLow, xUp, 1 ) :
	anaTool->GetYstarLabel( xLow, xUp, 1 );
      
      // build mean, sigma, project nev
      TH1* hMean = new TH1D
	( Form("h_%s_%s_%s_%s",
	       type.c_str(), sMean.c_str(), jzn.c_str(),
	       anaTool->GetName(xLow, xUp, "Ystar").c_str() ),
	  Form("%s;%s;%s", angleLabel.c_str(),
	       xAxisTitle.c_str(), yTitleMean.c_str() ),
	  nBinsY, yMin, yMax );
      vMeans[iG].push_back( hMean );
      
      TH1* hSigma = new TH1D
	( Form("h_%s_%s_%s_%s",
	       type.c_str(), sSigma.c_str(), jzn.c_str(),
	       anaTool->GetName(xLow, xUp, "Ystar").c_str()),
	  Form("%s;%s;%s", angleLabel.c_str(),
	       xAxisTitle.c_str(), yTitleSigma.c_str() ),
	  nBinsY, yMin, yMax );
      vSigmas[iG].push_back( hSigma );
      
      TH1* hNent = vJznNentIn[iG]->ProjectionY
	( Form("h_%s_%s_N_%s", type.c_str(), jzn.c_str(),
	       anaTool->GetName(xLow, xUp, "Ystar").c_str() )
	  , xBin, xBin );
      vNent[iG].push_back( hNent );

      anaTool->TruncateHistoBins( hNent );
      
      // now loop over the second axis and project onto z axis
      // to get the gaussian distributions. fit them.
      vProj[iG].resize( nBinsX );
      vFit [iG].resize( nBinsX );
      
      for( int yBin = 1; yBin <= nBinsY; yBin++ ){
	double yLow, yUp;
	anaTool->GetBinRange
	  ( hJznHin->GetYaxis(), yBin, yBin, yLow, yUp );

	int iX = xBin - 1;

	TH1* hProj = hJznHin->ProjectionZ
	  ( Form("h_%s_%s_%s_%s",  type.c_str(), jzn.c_str(), 
		 anaTool->GetName(xLow, xUp, "Ystar").c_str(),
		 anaTool->GetName(yLow, yUp, "Pt").c_str() ),
	    xBin, xBin, yBin, yBin );

	styleTool->SetHStyle( hProj, 0 );
	vProj[iG][iX].push_back( hProj );

 	TF1* fit = anaTool->FitGaussian( hProj );
	styleTool->SetHStyle( fit, 0 );
	vFit[iG][iX].push_back( fit );

	TCanvas c( hProj->GetName(), hProj->GetName(), 800, 600 );
	
	hProj->SetTitle("");
	hProj->Draw();

	fit->Draw("same");

	hProj->Write();
	
	DrawAtlasRight();
	drawTool->DrawLeftLatex( 0.18, 0.74, Form("%s", jzn.c_str() ) );

	std::string angleLabel = m_scaleResUseEta ?
	  anaTool->GetEtaLabel  ( xLow, xUp, 1 ) :
	  anaTool->GetYstarLabel( xLow, xUp, 1 );
	
	drawTool->DrawLeftLatex( 0.18, 0.88, angleLabel.c_str() );
	drawTool->DrawLeftLatex
	  ( 0.18, 0.81, anaTool->GetLabel( yLow, yUp, "#it{p}_{T}^{truth} [GeV]") );

	SaveAsROOT( c, hProj->GetName() );

	if( fit->GetParError(1) > 0.20 ){ continue; }
	if( fit->GetParError(2) > 0.02 ){ continue; }
 
	hMean->SetBinContent( yBin, fit->GetParameter(1) );
	hMean->SetBinError  ( yBin, fit->GetParError (1) );

	hSigma->SetBinContent ( yBin, fit->GetParameter(2) );
	hSigma->SetBinError   ( yBin, fit->GetParError (2) );
      }
    } // end loop over xBin (eta)
  }// end loop over jzn

  // Now, loop over the sigmas and means, and plot
  // as function of jz sample and eta
  std::vector< TH1* > vMeansFinal;
  std::vector< TH1* > vSigmasFinal;

  std::string xAxisTitle = vJznHin[0]->GetXaxis()->GetTitle();
  std::string yAxisTitle = "#it{p}_{T}^{truth} [GeV]";
  
  // make 2D histos with scale and resolution (mean, sigma)
  TH2D* hMeanAll = new TH2D
    ( Form("h_%s_%s", type.c_str(), sMean.c_str() ),
      Form(";%s;%s;%s", xAxisTitle.c_str(),
	   yAxisTitle.c_str(), yTitleMean.c_str() ),
      nBinsX, xMin, xMax,
      nBinsY, yMin, yMax );

  TH2D* hSigmaAll = new TH2D
    ( Form("h_%s_%s",
	   type.c_str(), sSigma.c_str() ),
      Form(";%s;%s;%s",
	   xAxisTitle.c_str(), yAxisTitle.c_str(),
	   yTitleSigma.c_str() ),
      nBinsX, xMin, xMax,
      nBinsY, yMin, yMax );

  hMeanAll ->GetXaxis()->Set( nBinsX, hRef->GetXaxis()->GetXbins()->GetArray() );
  hSigmaAll->GetXaxis()->Set( nBinsX, hRef->GetXaxis()->GetXbins()->GetArray() );
  
  for( int iX = 0; iX < nBinsX; iX++ ){
    int xBin = iX + 1;

    double xLow, xUp;
    anaTool->GetBinRange
      ( vJznHin[0]->GetXaxis(), xBin, xBin, xLow, xUp );
        
    std::vector< TH1* > vYstarMeanTemp;
    std::vector< TH1* > vYstarSigmaTemp;
    std::vector< TH1* > vYstarNentTemp;
    
    for( uint iG = 0; iG < m_nJzn; iG++ ){
      vYstarMeanTemp .push_back( vMeans [iG][iX] );
      vYstarSigmaTemp.push_back( vSigmas[iG][iX] );
      vYstarNentTemp .push_back( vNent  [iG][iX] );
    }

    std::string angleLabel = m_scaleResUseEta ?
      anaTool->GetEtaLabel  ( xLow, xUp, 1 ) :
      anaTool->GetYstarLabel( xLow, xUp, 1 );
    
    // Make the final, combined mean, sigma
    // Note - 
    // we use yAxisTitle from 3D histo
    // as xAxis on the final plots.
    TH1* hMeanFinal = new TH1D
      ( Form("h_%s_%s_%s",
	     type.c_str(), sMean.c_str(),
	     anaTool->GetName(xLow, xUp, "Ystar").c_str() ),
	Form("%s;%s;%s", angleLabel.c_str(),
	     yAxisTitle.c_str(), yTitleMean.c_str() ),
	nBinsY, yMin, yMax );
    vMeansFinal.push_back( hMeanFinal );
    
    TH1* hSigmaFinal = new TH1D
      ( Form("h_%s_%s_%s",
	     type.c_str(), sSigma.c_str(),
	     anaTool->GetName( xLow, xUp, "Ystar" ).c_str() ),
	Form("%s;%s;%s", angleLabel.c_str(),
	     yAxisTitle.c_str(), yTitleSigma.c_str() ),
	nBinsY, yMin, yMax );
    vSigmasFinal.push_back( hSigmaFinal );

    hMeanFinal ->GetXaxis()->
      SetRangeUser( m_varPtBinning.front(), m_varPtBinning.back() );
    hSigmaFinal->GetXaxis()->
      SetRangeUser( m_varPtBinning.front(), m_varPtBinning.back() );
    
    CombineSamples( hMeanFinal , vYstarMeanTemp , vYstarNentTemp );
    CombineSamples( hSigmaFinal, vYstarSigmaTemp, vYstarNentTemp );

    for( int yBin = 1; yBin <= nBinsY; yBin++ ){
      hMeanAll ->SetBinContent( xBin, yBin, hMeanFinal ->GetBinContent(yBin) );
      hSigmaAll->SetBinContent( xBin, yBin, hSigmaFinal->GetBinContent(yBin) );
      hMeanAll ->SetBinError( xBin, yBin, hMeanFinal ->GetBinError(yBin) );
      hSigmaAll->SetBinError( xBin, yBin, hSigmaFinal->GetBinError(yBin) );
    }
  } // end loop over xBin

  std::string nameFinal = "h_" + type;

  nameFinal += m_scaleResUseEta ? "_eta" : "_ystar";
  
  DrawCanvas( vMeansFinal , nameFinal.c_str(), sMean  );
  DrawCanvas( vSigmasFinal, nameFinal.c_str(), sSigma );

  std::string fOutName = "data/" + type + "_" + m_labelOut + ".root";
  TFile* fOut = new TFile( fOutName.c_str(), "RECREATE" );
  
  hMeanAll->Write();
  hSigmaAll->Write();

  fOut->Close();
  
  for( uint iG = 0; iG < m_nJzn; iG++){
    delete vMeansFinal [iG];
    delete vSigmasFinal[iG];
    for( int iX = 0; iX < nBinsX; iX++ ){
      delete vMeans[iG][iX];
      delete vSigmas[iG][iX];
      delete vNent[iG][iX];
      for( int iY = 0; iY < nBinsY; iY++ ){
	delete vProj[iG][iX][iY];
	delete vFit [iG][iX][iY];
      }
    }
  }
}

void DiJetAnalysisMC::MakeEfficiencies( std::vector< TH2* >& vSampleSpect,
					std::vector< TH2* >& vSampleSpectUP,
					const std::string& name ){  

  // if there are none, return
  if( !vSampleSpect.size() ){ return; }  

  std::string axisLabel, axisLabelTex;
  GetSpectraLabels( axisLabel, axisLabelTex, name );

  std::string gTitle = ";#it{p}_{T}^{truth} [GeV];#it{#varepsilon}_{reco}";

  // use this as reference because
  // it should be in every file
  TH2*  hRef = vSampleSpect[0];
  int nXbins = hRef->GetNbinsX();

  uint nSamples = vSampleSpect.size();
  
  std::vector< TH1* > vSpect;
  std::vector< TGraphAsymmErrors* > vEffGrf;
 
  std::vector< TGraphAsymmErrors* > vEffGrfFinal;
  
  double xMin = 0.; 
  double xMax = 70.; 

  TLine line( xMin, 1, xMax, 1 );
  line.SetLineWidth( 2 );

  TCanvas cEff("cEff","cEff",800,600);
  styleTool->SetCStyleGraph( cEff, xMin, m_effMin, xMax, 1.4, gTitle );

  double lX0 = 0.60;
  double lY0 = 0.25;
  double lX1 = 0.90;
  double lY1 = 0.55;

  TLegend leg( lX0, lY0, lX1, lY1 );
  styleTool->SetLegendStyle( &leg );
  
  // make individual jzn histograms
  for( int xBin = 1; xBin <= nXbins; xBin++ ){

    double xMin, xMax;
    anaTool->GetBinRange
      ( hRef->GetXaxis(), xBin, xBin, xMin, xMax );
      
    std::string hTag = anaTool->GetName( xMin, xMax, axisLabel );

    std::vector< TH1* > vSpectUPTmp;
    std::vector< TGraphAsymmErrors* > vEffGrfTmp;
    
    // loop over JZN samples
    for( uint iG = 0; iG < nSamples; iG++){

      std::string label = m_vJznLabels[ iG ];
      
      std::string hNameSpect
	= "h_" + name + "_" + label + "_" + hTag;
      std::string hNameSpectUP
	= "h_" + name + "_" + label + "_" + m_unpairedName + "_" + hTag;
      
      TH1* hSpect = static_cast< TH1D* >
	( vSampleSpect[ iG ]->ProjectionY( hNameSpect.c_str(), xBin, xBin ) );
      TH1* hSpectUP = static_cast< TH1D* >
	( vSampleSpectUP[ iG ]->ProjectionY( hNameSpectUP.c_str(), xBin, xBin ) );
      vSpectUPTmp.push_back( hSpectUP );
      vSpect.push_back( hSpect );
      vSpect.push_back( hSpectUP );
      
      std::string gEffName = "g_" + name + "_" + label + "_" + m_effName + "_" + hTag;
      
      TGraphAsymmErrors* g_etaEff = new TGraphAsymmErrors();
      g_etaEff->SetName( gEffName.c_str() ); 
      styleTool->SetHStyle( g_etaEff, xBin - 1 );
      vEffGrfTmp.push_back( g_etaEff );
      vEffGrf.push_back( g_etaEff );
      
      g_etaEff->Divide( hSpect, hSpectUP, "cl=0.683 b(1,1) mode" );
    } // end loop over JZN samples

    std::string gEffName = "g_" + name + "_" + m_allName + "_" + m_effName + "_" + hTag;

    TGraphAsymmErrors* gEffAll = CombineSamples( vEffGrfTmp, vSpectUPTmp );
    gEffAll->SetName( gEffName.c_str() );
    vEffGrfFinal.push_back( gEffAll );
    styleTool->SetHStyle( gEffAll, xBin - 1 );
    
    leg.AddEntry( gEffAll, anaTool->GetLabel( xMin, xMax, axisLabelTex ).c_str(), "lp" );
    
    gEffAll->Draw( "p" );
  }

  leg.Draw();

  line.Draw();
  
  DrawAtlasRight();
  
  std::string effName = "h_" + m_effName + "_" + m_allName;
  
  SaveAsAll( cEff, effName.c_str() ); 
  
  for( auto& h : vSpect       ){ delete h; }
  for( auto& g : vEffGrf      ){ delete g; }
  for( auto& g : vEffGrfFinal ){ delete g; }
}

void DiJetAnalysisMC::MakeSpectCFactorsRespMat( std::vector< TH2* >& vHspectReco,
						std::vector< TH2* >& vHspectTruth,
					        std::vector< TH3* >& vHspectRespMat,
						const std::vector< std::string >& vLabels,
						const std::string& nameCFactors,
						const std::string& nameRespMat ){

  std::string axisLabel, axisLabelTex;
  GetSpectraLabels( axisLabel, axisLabelTex, nameRespMat );  

  std::vector< TH1* > vSpect;
  std::vector< TH1* > vCfactors;
  std::vector< TH2* > vRespMat;
  std::vector< TF1* > vCfits;

  // these are needed to later combine the samples
  std::vector< std::vector< TH1* > > vCfactorsJzn;
  std::vector< std::vector< TH1* > > vNentTJzn;
  std::vector< std::vector< TH1* > > vNentRJzn;
  std::vector< std::vector< TH2* > > vRespMatJzn;

  std::vector< TH1* > vCfactorsAll;
  
  vCfactorsJzn.resize( m_nJzn );
  vNentTJzn   .resize( m_nJzn );
  vNentRJzn   .resize( m_nJzn );
  vRespMatJzn .resize( m_nJzn );

  std::string system = m_is_pPb ? "#it{p}+Pb" : "#it{pp}";
 
  for( uint iG = 0; iG < vHspectRespMat.size(); iG++ ){

    std::string label = vLabels[iG];

    TH2* hReco    = vHspectReco   [ iG ];
    TH2* hTruth   = vHspectTruth  [ iG ];
    TH3* hRespMatFull = vHspectRespMat[ iG ];
    
    for( int xBin = 1; xBin <= hRespMatFull->GetNbinsX(); xBin++ ){

      double xMin, xMax;
      anaTool->GetBinRange
	( hRespMatFull->GetXaxis(), xBin, xBin, xMin, xMax );

      hRespMatFull->GetXaxis()->SetRange( xBin, xBin );

      std::string hTag = anaTool->GetName( xMin, xMax, axisLabel);
      
      TH2* hRespMat = static_cast< TH2D* >( hRespMatFull->Project3D("zy") );
      hRespMat->SetName( Form("h_%s_%s_%s", nameRespMat.c_str(), label.c_str(), hTag.c_str() ) ) ;
      styleTool->SetHStyle( hRespMat, 0 );
      vRespMat.push_back( hRespMat );

      TCanvas cSpectRespMat( hRespMat->GetName(), hRespMat->GetName(), 800, 600 );
      cSpectRespMat.SetLogz();

      gStyle->SetPaintTextFormat(".1e");
      // anaTool->AverageOver( hRespMat, "row" );
      
      hRespMat->Draw("col text25");
      hRespMat->SetTitle("");

      drawTool->DrawLeftLatex
	( 0.2, 0.85, anaTool->GetLabel( xMin, xMax, axisLabelTex ).c_str() );

      drawTool->DrawAtlasInternalMCRight( CT::DrawTools::drawX0, CT::DrawTools::drawY0,
					  m_mcTypeLabel, 3, CT::StyleTools::lSS );

      drawTool->DrawRightLatex( 0.87, 0.25, Form( "%s %s", system.c_str(), m_mcTypeLabel.c_str() ) );
      
      hRespMat->Write();

      if( !label.compare( m_allName ) ){
	SaveAsAll( cSpectRespMat, hRespMat->GetName() );
      }

      TH1* hR = static_cast< TH1D* >
	( hReco->ProjectionY( Form("h_R_%s_%s_%s_%s",
				   m_ystarSpectName.c_str(), m_sCounts.c_str(),
				   label.c_str(), hTag.c_str() ) , xBin, xBin ) );

      TH1* hT = static_cast< TH1D* >
	( hTruth->ProjectionY( Form("h_T_%s_%s_%s_%s",
				    m_ystarSpectName.c_str(), m_sCounts.c_str(),
				    label.c_str(), hTag.c_str() ) , xBin, xBin ) );

      TH1* hC = static_cast< TH1D* >
	  ( hT->Clone( Form("h_%s_%s_%s",
			    nameCFactors.c_str(), label.c_str(), hTag.c_str() ) ) );
      styleTool->SetHStyleRatio( hC );	  
      hC->Divide( hR );
      vCfactors.push_back( hC );

      hC->SetYTitle( "Correction Factor" );

      SetCfactorsErrors( hR, hT, hRespMat, hC );

      std::cout << hT->GetName()  << " ---- " << hT->GetNbinsX() << std::endl; 
      for( int i = 1; i <= hT->GetNbinsX(); i++ ){
	std::cout << hT->GetBinContent(i) << " " << hR->GetBinContent(i) << " " << hC->GetBinContent(i) << std::endl;
      }
      
      // Do the fitting and drawing only for
      // the jzn individual samples. Do this later
      // for the combined sample.
      if( label.compare( m_allName ) ){
	      
	TF1* cFit = anaTool->FitPol2
	  ( hC, hC->GetXaxis()->GetXmin(), hC->GetXaxis()->GetXmax() );
	vCfits.push_back( cFit );
	    
	cFit->SetLineStyle( 2 );
	cFit->SetLineColor( kBlue );

	TCanvas cCFactors( hC->GetName(), hC->GetName(), 800, 600 );
	cCFactors.SetLogz();
	    
	hC->Draw("ep");
	cFit->Draw("same");

	drawTool->DrawLeftLatex
	  ( 0.2, 0.85, anaTool->GetLabel( xMin, xMax, axisLabelTex ).c_str() );
      	    
	DrawAtlasRight();

	hC->Write();
	cFit->Write();

	// keep track of only JZ slice counts and response matrixes.
	// these vectors are the ones that will be used to make the
	// the final correction factors.
	// For CFactor histograms, only keep track of ones with
	// "All". These factors will be changed on CombineSamples.
	vCfactorsJzn[ iG ].push_back( hC );
	vNentTJzn   [ iG ].push_back( hT );
	vNentRJzn   [ iG ].push_back( hR );
	vRespMatJzn [ iG ].push_back( hRespMat );

	// SaveAsAll( cCFactors, hC->GetName() );
      } else {
	vCfactorsAll.push_back( hC );
      }
    }
  }

  // now combine to make the final samples
  for( uint spectBin = 0; spectBin < vCfactorsAll.size(); spectBin++ ){

    // prepare all jzn samples for each ystar1,2 pt1,2 bin
    std::vector< TH1* > vCfactorsTmp;
    std::vector< TH1* > vNentTTmp;
    std::vector< TH1* > vNentRTmp;
    std::vector< TH2* > vRespMatTmp;

    for( uint iG = 0; iG < m_nJzn; iG++ ){
      vCfactorsTmp.push_back( vCfactorsJzn[ iG ][ spectBin ] );
      vNentTTmp   .push_back( vNentTJzn   [ iG ][ spectBin ] );
      vNentRTmp   .push_back( vNentRJzn   [ iG ][ spectBin ] );
      vRespMatTmp .push_back( vRespMatJzn [ iG ][ spectBin ] );
    }

    TH1* hCAll = vCfactorsAll[ spectBin ]; 
    CombineSamples( hCAll, vCfactorsTmp, vNentTTmp, vNentRTmp, vRespMatTmp );

    hCAll->Write();
  }
  

  
  for( auto& h : vSpect    ){ delete h; }
  for( auto& h : vCfactors ){ delete h; }
  for( auto& h : vRespMat  ){ delete h; }
}

void DiJetAnalysisMC::MakeDphiCFactorsRespMat( std::vector< THnSparse* >& vHnT,
					       std::vector< THnSparse* >& vHnR,
					       std::vector< THnSparse* >& vHnDphiRespMat,
					       const std::vector< std::string >& vLabels,
					       const std::string& nameCFactors,
					       const std::string& nameRespMat ){

  // these are just for garbage collection
  std::vector< TH1* > vProj;
  std::vector< TH1* > vCfactors;
  std::vector< TH2* > vRespMat;
  std::vector< TF1* > vCfits;
  std::vector< TH1* > vCfitsQuality;
  
  // these are needed to later combine the samples
  std::vector< std::vector< TH1* > > vCfactorsJzn;
  std::vector< std::vector< TH1* > > vNentTJzn;
  std::vector< std::vector< TH1* > > vNentRJzn;
  std::vector< std::vector< TH2* > > vRespMatJzn;

  std::vector< TH1* > vCfactorsAll;
  
  vCfactorsJzn.resize( m_nJzn );
  vNentTJzn   .resize( m_nJzn );
  vNentRJzn   .resize( m_nJzn );
  vRespMatJzn .resize( m_nJzn );

  std::string system = m_is_pPb ? "#it{p}+Pb" : "#it{pp}";
  
  // ---- loop over group  ----
  // ---- (jzn or trigger) ----
  for( uint iG = 0; iG < vHnT.size(); iG++ ){      
    std::string label = vLabels[iG];

    TH1D* hCfactorsProb =
      new TH1D( Form( "hCfactorsProb_%s", label.c_str() ), ";Fit Probability", 25, 0, 1 );
    vCfitsQuality.push_back( hCfactorsProb );
	
    TH1D* hCfactorsChi2 =
      new TH1D( Form( "hCfactorsChi2_%s", label.c_str() ), ";#Chi^{2}/NDF", 25, 0, 10 );
    vCfitsQuality.push_back( hCfactorsChi2 );

    TH2D* hMiiTiRi = new TH2D( Form( "hMiiTiRi_%s", label.c_str() ),
			       ";M_{ii}/R_{i};M_{ii}/T_{i}",
			       50, 0, 1,
			       50, 0, 1 );
    vRespMat.push_back( hMiiTiRi );

    THnSparse* hnT    = vHnT[iG];
    THnSparse* hnR    = vHnR[iG];
    THnSparse* hnDphi = vHnDphiRespMat[iG];
    
    // just in case, reset the histogram's ranges
    // in case they were set in previous functions
    anaTool->ResetAxisRanges( hnT    );
    anaTool->ResetAxisRanges( hnR    );
    anaTool->ResetAxisRanges( hnDphi );
    
    TAxis* axisT0 = hnT->GetAxis( m_dPP->GetAxisI(0) ); int nAxis0Bins = axisT0->GetNbins();
    TAxis* axisT1 = hnT->GetAxis( m_dPP->GetAxisI(1) ); int nAxis1Bins = axisT1->GetNbins();
    TAxis* axisT2 = hnT->GetAxis( m_dPP->GetAxisI(2) ); int nAxis2Bins = axisT2->GetNbins();
    TAxis* axisT3 = hnT->GetAxis( m_dPP->GetAxisI(3) ); int nAxis3Bins = axisT3->GetNbins();
    
    TAxis* axisR0 = hnR->GetAxis( m_dPP->GetAxisI(0) ); 
    TAxis* axisR1 = hnR->GetAxis( m_dPP->GetAxisI(1) ); 
    TAxis* axisR2 = hnR->GetAxis( m_dPP->GetAxisI(2) ); 
    TAxis* axisR3 = hnR->GetAxis( m_dPP->GetAxisI(3) );

    
    for( int axis0Bin = 1; axis0Bin <= nAxis0Bins; axis0Bin++ ){ 
      // check we are in correct ystar and pt bins
      if( !m_dPP->CorrectPhaseSpace
	  ( std::vector<int>{ axis0Bin, 0, 0, 0 } ) )
	{ continue; }

      axisT0->SetRange( axis0Bin, axis0Bin );
      axisR0->SetRange( axis0Bin, axis0Bin );
      double axis0Low, axis0Up;
      anaTool->GetBinRange( axisT0, axis0Bin, axis0Bin, axis0Low, axis0Up );
	   
      for( int axis1Bin = 1; axis1Bin <= nAxis1Bins; axis1Bin++ ){	
	axisT1->SetRange( axis1Bin, axis1Bin ); 
	axisR1->SetRange( axis1Bin, axis1Bin );
	double axis1Low, axis1Up;
	anaTool->GetBinRange( axisT1, axis1Bin, axis1Bin, axis1Low, axis1Up );

	for( int axis2Bin = 1; axis2Bin <= nAxis2Bins; axis2Bin++ ){
	  axisT2->SetRange( axis2Bin, axis2Bin ); 
	  axisR2->SetRange( axis2Bin, axis2Bin );
	  double axis2Low, axis2Up;
	  anaTool->GetBinRange( axisT2, axis2Bin, axis2Bin, axis2Low, axis2Up );

	  for( int axis3Bin = 1; axis3Bin <= nAxis3Bins; axis3Bin++ ){
	    // check we are in correct ystar and pt bins
	    if( !m_dPP->CorrectPhaseSpace
		( std::vector<int>{ axis0Bin, axis1Bin, axis2Bin, axis3Bin } ) )
	      { continue; }
	    
	    axisT3->SetRange( axis3Bin, axis3Bin ); 
	    axisR3->SetRange( axis3Bin, axis3Bin ); 	    
	    double axis3Low, axis3Up;
	    anaTool->GetBinRange( axisT3, axis3Bin, axis3Bin, axis3Low, axis3Up );
	    
	    std::string hTag =
	      Form( "%s_%s_%s_%s",
		    anaTool->GetName( axis0Low, axis0Up, m_dPP->GetAxisName(0) ).c_str(),
		    anaTool->GetName( axis1Low, axis1Up, m_dPP->GetAxisName(1) ).c_str(),
		    anaTool->GetName( axis2Low, axis2Up, m_dPP->GetAxisName(2) ).c_str(),
		    anaTool->GetName( axis3Low, axis3Up, m_dPP->GetAxisName(3) ).c_str() ); 

	    //---------------------------------------------------
	    //                 Response Matrixes
	    //---------------------------------------------------
	    // -------------- make dphi resp matrix -------------
	    // Here the axis are in order. ystar1 = axis0, etc.
	    // we need to get the mapped axis for this to be true.
	    // because axis0Bin as of now could be a pT axis.
	    // The mapped is back to original. Once we get that
	    // vector, the [0] is ystar1, [1] is ystar2, etc.
	    // Before, axis0, axis1 are arbitrary (chosen from m_dPP)
	    std::vector< int > vMappedaAxisBins = m_dPP->GetMappedBins
	      ( std::vector<int>{ axis0Bin, axis1Bin, axis2Bin, axis3Bin } );

	    int ystar1Bin =  vMappedaAxisBins[0];
	    int ystar2Bin =  vMappedaAxisBins[1];
	    // because of underflow bin we add 1 to whatever 
  	    // bin it actually is first actual bin is 2 not 1.
	    int pt1Bin    =  vMappedaAxisBins[2] + 1; 
	    int pt2Bin    =  vMappedaAxisBins[3] + 1;

	    // set the ranges, project onto dPhi x dPhi axis.
	    // for 
	    hnDphi->GetAxis(0)->SetRange( ystar1Bin, ystar1Bin  );
	    hnDphi->GetAxis(1)->SetRange( ystar2Bin, ystar2Bin  );
	    hnDphi->GetAxis(2)->SetRange( pt1Bin   , pt1Bin     );
	    hnDphi->GetAxis(3)->SetRange( pt1Bin   , pt1Bin     );
	    hnDphi->GetAxis(4)->SetRange( pt2Bin   , pt2Bin     );
	    hnDphi->GetAxis(5)->SetRange( pt2Bin   , pt2Bin     );

	    uint dPhiProjDimX = hnDphi->GetNdimensions() - 2;
	    uint dPhiProjDimY = hnDphi->GetNdimensions() - 1;
	    
	    TH2* hDphiRespMat = static_cast<TH2D*>
	      ( hnDphi->Projection( dPhiProjDimY, dPhiProjDimX ) );
	    hDphiRespMat->SetName
	      ( Form( "h_%s_%s_%s", nameRespMat.c_str(), label.c_str(), hTag.c_str() ) );
	    styleTool->SetHStyle( hDphiRespMat, 0 );
	    vRespMat.push_back( hDphiRespMat );
	    
	    TCanvas cDphiRespMat( hDphiRespMat->GetName(), hDphiRespMat->GetName(), 800, 600 );
	    cDphiRespMat.SetLogz();

	    hDphiRespMat->Draw("col");
	    hDphiRespMat->SetTitle("");

	    hDphiRespMat->GetXaxis()->SetRange
	      ( m_dPhiRebinnedZoomLowBin, m_dPhiRebinnedZoomHighBin );
	    hDphiRespMat->GetYaxis()->SetRange
	      ( m_dPhiRebinnedZoomLowBin, m_dPhiRebinnedZoomHighBin );
	    
	    DrawTopLeftLabels
	      ( m_dPP, axis0Low, axis0Up, axis1Low, axis1Up,
		axis2Low, axis2Up, axis3Low, axis3Up );
	    
	    drawTool->DrawAtlasInternalMCRight( CT::DrawTools::drawX0, CT::DrawTools::drawY0,
						m_mcTypeLabel, 3, CT::StyleTools::lSS);
	
	    drawTool->DrawRightLatex
	      ( 0.87, 0.25, Form( "%s %s", system.c_str(), m_mcTypeLabel.c_str() ) );
    
	    hDphiRespMat->Write();

	    if( !label.compare( m_allName ) )
	      { SaveAsAll( cDphiRespMat, hDphiRespMat->GetName() ); }
	
	    //---------------------------------------------------
	    //   Counts, Rebinned, Normalized dPhi Histograms
	    //---------------------------------------------------
	    TH1* hT = hnT->Projection( 4 );
	    TH1* hR = hnR->Projection( 4 );

	    // TRUNCATE HISTOGRAMS
	    for( int i = 0; i <= hT->GetNbinsX(); i++ ){

	      if( hT->GetBinContent(i) <= 2 ){
		hT->SetBinContent( i, 0 );
	      }
	      if( hR->GetBinContent(i) <= 2 ){
		hR->SetBinContent( i, 0 );
	      }
	    }
	   
	    
	    hT->SetName( Form( "h_T_%s_%s_%s_%s",
			       m_dPhiName.c_str(), m_sCounts.c_str(),
			       label.c_str(), hTag.c_str() ) );
	    hR->SetName( Form( "h_R_%s_%s_%s_%s",
			       m_dPhiName.c_str(), m_sCounts.c_str(),
			       label.c_str(), hTag.c_str() ) );
	    
	    vProj.push_back( hT );
	    vProj.push_back( hR );

	    hT->Write();
	    hR->Write();

	    TH1* hTreb = static_cast< TH1D* >
	      ( hT->Rebin
		( m_nVarDphiRebinnedBins,
		  Form( "h_T_%s_%s_%s_%s_%s",
			m_dPhiName.c_str(), m_sCounts.c_str(),
			m_sReb.c_str(), label.c_str(), hTag.c_str() ),
		  &m_varDphiRebinnedBinning[0] ) );   
	    TH1* hRreb = static_cast< TH1D* >
	      ( hR->Rebin
		( m_nVarDphiRebinnedBins,
		  Form( "h_R_%s_%s_%s_%s_%s",
			m_dPhiName.c_str(), m_sCounts.c_str(),
			m_sReb.c_str(), label.c_str(), hTag.c_str() ),
		  &m_varDphiRebinnedBinning[0] ) );
   
	    vProj.push_back( hTreb );
	    vProj.push_back( hRreb );
   
	    hTreb->Write();
	    hRreb->Write();
	    
	    //---------------------------------------------------
	    //               Correction Factors
	    //---------------------------------------------------
	    
	    // Make the ratio histogram. This is used to unfold bin-by-bin.
	    // later these are combined on a per-jz sample basis.
	    TH1* hC = static_cast< TH1D* >
	      ( hTreb->Clone
		( Form( "h_%s_%s_%s", nameCFactors.c_str(), label.c_str(), hTag.c_str())));
	    styleTool->SetHStyleRatio( hC );	  
	    hC->Divide( hRreb );
	    vCfactors.push_back( hC );

	    hC->SetYTitle( "Correction Factor" );

	    SetCfactorsErrors( hRreb, hTreb, hDphiRespMat, hC );
	    	    
	    // Smooth the cfactors
	    hC->GetXaxis()->SetRange
	      ( m_dPhiRebinnedZoomLowBin, m_dPhiRebinnedZoomHighBin );
	    // hC->Smooth( 1, "R" );

	    // Do the fitting and drawing only for
	    // the jzn individual samples. Do this later
	    // for the combined sample.
	    if( label.compare( m_allName ) ){
	      // keep track of only JZ slice counts and response matrixes.
	      // these vectors are the ones that will be used to make the
	      // the final correction factors.
	      // For CFactor histograms, only keep track of ones with
	      // "All". These factors will be changed on CombineSamples.
	      vCfactorsJzn[ iG ].push_back( hC );
	      vNentTJzn   [ iG ].push_back( hTreb );
	      vNentRJzn   [ iG ].push_back( hRreb );
	      vRespMatJzn [ iG ].push_back( hDphiRespMat );
	    } else {
	      vCfactorsAll.push_back( hC );
	    }
	  } // end loop over axis4
	} // end loop over axis2
      } // end loop over axis1
    } // end loop over axis0
    
    hCfactorsProb->Write();
    hCfactorsChi2->Write();
    hMiiTiRi     ->Write();
  } // end loop over jzn sample

  // now combine to make the final samples
  for( uint ystarPtBin = 0; ystarPtBin < vCfactorsAll.size(); ystarPtBin++ ){

    // prepare all jzn samples for each ystar1,2 pt1,2 bin
    std::vector< TH1* > vCfactorsTmp;
    std::vector< TH1* > vNentTTmp;
    std::vector< TH1* > vNentRTmp;
    std::vector< TH2* > vRespMatTmp;

    for( uint iG = 0; iG < m_nJzn; iG++ ){
      vCfactorsTmp.push_back( vCfactorsJzn[ iG ][ ystarPtBin ] );
      vNentTTmp   .push_back( vNentTJzn   [ iG ][ ystarPtBin ] );
      vNentRTmp   .push_back( vNentRJzn   [ iG ][ ystarPtBin ] );
      vRespMatTmp .push_back( vRespMatJzn [ iG ][ ystarPtBin ] );
    }

    TH1* hCAll = vCfactorsAll[ ystarPtBin ]; 
    CombineSamples( hCAll, vCfactorsTmp, vNentTTmp, vNentRTmp, vRespMatTmp );

    hCAll->Write();
  }
  
  for( auto& f : vCfits        ){ delete f; }
  for( auto& h : vProj         ){ delete h; }
  for( auto& f : vCfactors     ){ delete f; }
  for( auto& r : vRespMat      ){ delete r; }
  for( auto& h : vCfitsQuality ){ delete h; }
}


void DiJetAnalysisMC::MakeYstarRespMat( std::vector< TH3* >& vhnYstar,
					const std::vector< std::string >& vLabels,
					const std::string& nameYstar ){

  std::string system = m_is_pPb ? "#it{p}+Pb" : "#it{pp}";
  
  // ---- loop over group  ----
  // ---- (jzn or trigger) ----
  for( uint iG = 0; iG < vhnYstar.size(); iG++ ){      
    std::string label = vLabels[iG];

    // only draw for all
    if( label.compare( m_allName ) ){ continue; } 
    
    TH3* hYstarPt = vhnYstar[iG];
    
    TAxis* ptAxis = hYstarPt->GetZaxis();
    
    for( uint ptBin = 1; ptBin <= m_nVarPtBins; ptBin++ ){

      ptAxis->SetRange( ptBin, ptBin );
      
      double ptLow = ptAxis->GetBinLowEdge( ptBin );
      double ptUp  = ptAxis->GetBinUpEdge ( ptBin );

      TH2* hYstarRespMat = static_cast< TH2D* >
	( hYstarPt->Project3D("yx") );
      styleTool->SetHStyle( hYstarRespMat, 0 );
      
      std::string hName = "h_yStarRespMat_" + anaTool->GetName( ptLow, ptUp, "Pt" );
      
      hYstarRespMat->SetName( hName.c_str() );
      hYstarRespMat->SetTitle("");

      // normalize rows to 1 to get purity
      gStyle->SetPaintTextFormat(".3f");
      anaTool->AverageOver( hYstarRespMat, "row" );
 
      TCanvas c("c","c", 800, 600 );
      c.SetLogz();

      hYstarRespMat->Draw("col");
      hYstarRespMat->Draw("text same");

      drawTool->DrawAtlasInternalMCRight( CT::DrawTools::drawX0, CT::DrawTools::drawY0,
					  m_mcTypeLabel, 3, CT::StyleTools::lSS );
      drawTool->DrawLeftLatex( 0.19, 0.87, Form( "%s %s", system.c_str(), m_mcTypeLabel.c_str() ) );
      drawTool->DrawRightLatex( 0.87, 0.25, anaTool->GetLabel( ptLow, ptUp, ptAxis->GetTitle() ) );  
      
      hYstarRespMat->Write();
      SaveAsAll( c, hName );

      delete hYstarRespMat;
    } // end loop over pT bin
  }
}


void DiJetAnalysisMC::MakePtRespMat( std::vector< THnSparse* >& vhnPt,
				     const std::vector< std::string >& vLabels,
				     const std::string& namePt ){  
  
  // ----------- pT response matrix -------------
  std::vector< TH2* > vRespMat;
  std::vector< TH1* > vPurityEff;
 
  // ---- loop over group  ----
  // ---- (jzn or trigger) ----
  for( uint iG = 0; iG < vhnPt.size(); iG++ ){      
    // only draw for all
    std::string label = vLabels[iG];

    if( label.compare( m_allName ) ){ continue; } 

    THnSparse* hnPt   = vhnPt  [iG];

    // just in case, reset the axis ranges if this
    // histogram was used in previous functions
    anaTool->ResetAxisRanges( hnPt );
    
    // want the ystar axes (0, 1);
    // dont need to do anything fancy like
    // get mapped axis.
    TAxis* axis0yStar = hnPt->GetAxis(0);
    TAxis* axis1yStar = hnPt->GetAxis(1);

    TAxis* dPhiRaxis  = hnPt->GetAxis( 6 );
    TAxis* dPhiTaxis  = hnPt->GetAxis( 7 );

    int nAxis0BinsYstar = axis0yStar->GetNbins();
    int nAxis1BinsYstar = axis1yStar->GetNbins();

    int nDphiRbins      = dPhiRaxis->GetNbins();
    int nDphiTbins      = dPhiTaxis->GetNbins();

    // -------------- do resp mat for Pt    
    for( int axis0Bin = 1; axis0Bin <= nAxis0BinsYstar; axis0Bin++ ){ 
      // check we are in correct ystar and pt bins
      if( !m_dPP->CorrectPhaseSpace
	  ( std::vector<int>{ axis0Bin, 0, 0, 0 } ) )
	{ continue; }

      axis0yStar->SetRange( axis0Bin, axis0Bin );
      double axis0Low, axis0Up;
      anaTool->GetBinRange
	( axis0yStar, axis0Bin, axis0Bin, axis0Low, axis0Up );

      for( int axis1Bin = 1; axis1Bin <= nAxis1BinsYstar; axis1Bin++ ){

	axis1yStar->SetRange( axis1Bin, axis1Bin );
	double axis1Low, axis1Up;
	anaTool->GetBinRange
	  ( axis1yStar, axis1Bin, axis1Bin, axis1Low, axis1Up );
	
	// Make pT resp matrix
	std::string hTag =
	  Form( "%s_%s",
		anaTool->GetName( axis0Low, axis0Up, m_dPP->GetDefaultAxisName(0) ).c_str(),
		anaTool->GetName( axis1Low, axis1Up, m_dPP->GetDefaultAxisName(1) ).c_str() ); 

	std::vector< std::string > vPtLabel{ m_s_pt1, m_s_pt2 };
	std::vector< int > vPtProjDimY{ 3, 5 };
	std::vector< int > vPtProjDimX{ 2, 4 };

	// Reset the dphi bin ranges so we can still get
	// total pT migration. Do not call ResetAxisRanges
	// because we are in some range in ystar atm.
	dPhiRaxis->SetRange( 1, nDphiRbins );
	dPhiTaxis->SetRange( 1, nDphiTbins );
	
	// do the same for pt1 and pt2
	for( uint ptN = 0; ptN < vPtLabel.size(); ptN++ ){
	  // Take projection onto the two pT axis
	  TH2* hPtRespMat = static_cast< TH2D* >
	    ( hnPt->Projection( vPtProjDimY[ ptN ], vPtProjDimX[ ptN ] ) );
	  hPtRespMat->SetName
	    ( Form( "h_%s_%s_%s_%s", namePt.c_str(), vPtLabel[ ptN ].c_str(),
		    label.c_str(), hTag.c_str() ) );
	  styleTool->SetHStyle( hPtRespMat, 0 );
	  vRespMat.push_back( hPtRespMat );

	  // kind of a pointless bin
	  hPtRespMat->SetBinContent( hPtRespMat->FindBin( 20., 20. ), 0 );
     	    
	  TCanvas c1( hPtRespMat->GetName(), hPtRespMat->GetName(), 800, 600 );
	  c1.SetLogz();

	  gStyle->SetPaintTextFormat(".3f");
	  anaTool->AverageOver( hPtRespMat, "row" );
	  hPtRespMat->Draw("col text");
	  hPtRespMat->SetTitle("");

	  drawTool->DrawLeftLatex
	    ( 0.13, 0.86, CT::AnalysisTools::GetLabel
	      ( axis0Low, axis0Up, m_dPP->GetDefaultAxisLabel(0) ) );
	  drawTool->DrawLeftLatex
	    ( 0.13, 0.79, CT::AnalysisTools::GetLabel
	      ( axis1Low, axis1Up, m_dPP->GetDefaultAxisLabel(1) ) );  
	
	  DrawAtlasRight();
	
	  SaveAsAll( c1, hPtRespMat->GetName() );
	  hPtRespMat->Write();

	  TH1* hPtPurity = static_cast<TH1D*>
	    (hPtRespMat->ProjectionX
	     ( Form( "h_%s_%s_%s_%s_%s", namePt.c_str(), vPtLabel[ptN].c_str(),
		     m_purityName.c_str(), label.c_str(), hTag.c_str())) );
	  TH1* hPtEff    = static_cast<TH1D*>
	    ( hPtRespMat->ProjectionY
	      ( Form( "h_%s_%s_%s_%s_%s", namePt.c_str(), vPtLabel[ptN].c_str(),
		      m_effName.c_str(), label.c_str(), hTag.c_str())) );
	  vPurityEff.push_back( hPtPurity );
	  vPurityEff.push_back( hPtEff    );
	  
	  hPtPurity->Reset();
	  hPtEff->Reset();

	  styleTool->SetHStyle( hPtPurity, 0 );
	  styleTool->SetHStyle( hPtEff   , 1 );
	  
	  GetPurityEff( hPtRespMat, hPtPurity, hPtEff );
	  
	  TCanvas c2( hPtPurity->GetName(), hPtPurity->GetName(), 800, 600 );
	  TLegend leg( 0.60, 0.22, 0.90, 0.33 );
	  styleTool->SetLegendStyle( &leg );

	  leg.AddEntry( hPtPurity,  "Purity" );
	  leg.AddEntry( hPtEff, "Efficiency" );

	  hPtPurity->SetMaximum(1.0);
	  hPtPurity->SetMinimum(0.0);
	  hPtEff   ->SetMaximum(1.0);
	  hPtEff   ->SetMinimum(0.0);
	  
	  hPtPurity->Draw("ep same");
	  hPtEff   ->Draw("ep same");
	  leg.Draw();

	  drawTool->DrawLeftLatex
	    ( 0.13, 0.86, CT::AnalysisTools::GetLabel
	      ( axis0Low, axis0Up, m_dPP->GetDefaultAxisLabel(0) ) );
	  drawTool->DrawLeftLatex
	    ( 0.13, 0.79, CT::AnalysisTools::GetLabel
	      ( axis1Low, axis1Up, m_dPP->GetDefaultAxisLabel(1) ) );  
	  
	  DrawAtlasRight();
	
	  SaveAsROOT( c2, hPtPurity->GetName() );

	  // Now loop over delta Phi bins, and make these pt projections
	  // over some combination (top 4x4). Dont forget to reset
	  // the delta Phi ranges after this!

	  // keep this same for both of them for now.
	  // assume the have > 4 bins.
	  int dPhiBinRange = dPhiRaxis->GetNbins() > 4 ? 4 : dPhiRaxis->GetNbins();
	  int dPhiBinBegin = nDphiRbins - dPhiBinRange + 1;
	  for( int dPhiRbin = dPhiBinBegin; dPhiRbin <= nDphiRbins; dPhiRbin++ ){
	    dPhiRaxis->SetRange( dPhiRbin, dPhiRbin );
	    	  
	    for( int dPhiTbin = dPhiBinBegin; dPhiTbin <= nDphiTbins; dPhiTbin++ ){
	      dPhiTaxis->SetRange( dPhiTbin, dPhiTbin );
	      
	      for( uint ptN = 0; ptN < vPtLabel.size(); ptN++ ){
		// Make pT dPhi resp matrix
		std::string hTagDphi =
		  Form( "%s_dPhiR_%d_dPhiT_%d", hTag.c_str(), dPhiRbin, dPhiTbin ); 
		
		TH2* hPtRespMatDphi = static_cast< TH2D* >
		  ( hnPt->Projection( vPtProjDimY[ ptN ], vPtProjDimX[ ptN ] ) );	      	      
		hPtRespMatDphi->SetName
		  ( Form( "h_%s_%s_%s_%s_%s", namePt.c_str(), m_dPhiName.c_str(),
			  vPtLabel[ ptN ].c_str(), label.c_str(), hTagDphi.c_str() ) );
		styleTool->SetHStyle( hPtRespMatDphi, 0 );
		vRespMat.push_back( hPtRespMatDphi );

		// kind of a pointless bin
		hPtRespMatDphi->SetBinContent( hPtRespMatDphi->FindBin( 20., 20. ), 0 );
	      
		TCanvas c1( hPtRespMatDphi->GetName(), hPtRespMatDphi->GetName(), 800, 600 );
		c1.SetLogz();

		gStyle->SetPaintTextFormat(".3f");
		anaTool->AverageOver( hPtRespMatDphi, "row" );
		hPtRespMatDphi->Draw("col text");
		hPtRespMatDphi->SetTitle("");

		drawTool->DrawLeftLatex
		  ( 0.13, 0.86, CT::AnalysisTools::GetLabel
		    ( axis0Low, axis0Up, m_dPP->GetDefaultAxisLabel(0) ) );
		drawTool->DrawLeftLatex
		  ( 0.13, 0.79, CT::AnalysisTools::GetLabel
		    ( axis1Low, axis1Up, m_dPP->GetDefaultAxisLabel(1) ) );  
		drawTool->DrawLeftLatex
		  ( 0.4, 0.86, Form( "#delta#phi_{reco} Bin %d" , dPhiRbin ) );
		drawTool->DrawLeftLatex
		  ( 0.4, 0.79, Form( "#delta#phi_{truth} Bin %d", dPhiTbin ) );  
		DrawAtlasRight();
	
		SaveAsAll( c1, hPtRespMatDphi->GetName() );
		hPtRespMatDphi->Write();
	      } // end loop over ptN
	    } // end loop over dPhiTbin
	  } // end loop over dPhiRbin
	  
	} // end loop over ptN
      } // end loop over ystar2
    } // end loop over ystar1
  } // end loop over iG
  for( auto rm : vRespMat     ){ delete rm; }
  for( auto pe : vPurityEff   ){ delete pe; }
}

// CLEAN THIS UP! ITS NOT WORTH THE HEADACHE NOW
// BUT STILL, NOT GOOD (04.07.18)
void DiJetAnalysisMC::CompareAngularRes( TFile* fOut ){

  std::vector< TH1* > vAngRes;
  std::vector< TH1* > vD;

   // get the pythia8 eta sigm+mean file 
  TFile* f_eta_p8 = new TFile( "data/recoTruthDeta_pp_mc_pythia8.root", "read" );
  TH2D* hAngularResEta = static_cast< TH2D* >( f_eta_p8->Get( "h_recoTruthDeta_sigma" ) );
  hAngularResEta->SetName( "hAngularResEta" );
  hAngularResEta->SetDirectory(0);
  f_eta_p8->Close();
  delete f_eta_p8;
  
  // clone this before modifying it.
  hAngularResEta = static_cast< TH2D* >( hAngularResEta->Clone( "hAngularResEta" ) );

  // get the reference herwig eta sigm+mean file 
  TFile* f_eta_ref = new TFile( "data/recoTruthDeta_pp_mc_herwig.root", "read" );
  TH2D* hAngularResEtaTmp = static_cast< TH2D* >( f_eta_ref->Get( "h_recoTruthDeta_sigma" ) );
  hAngularResEta->SetName( "hAngularResEtaTmp" );
  hAngularResEtaTmp->SetDirectory(0);
  f_eta_ref->Close();
  delete f_eta_ref;

  // get the pythia8 phi sigm+mean file 
  TFile* f_phi_p8 = new TFile( "data/recoTruthDphi_pp_mc_pythia8.root", "read" );
  TH2D* hAngularResPhi = static_cast< TH2D* >( f_phi_p8->Get( "h_recoTruthDphi_sigma" ) );
  hAngularResPhi->SetName( "hAngularResPhi" );
  hAngularResPhi->SetDirectory(0);
  f_phi_p8->Close();
  delete f_phi_p8;

  // clone this before modifying it.
  hAngularResPhi = static_cast< TH2D* >( hAngularResPhi->Clone( "hAngularResPhi" ) );

    // get the reference herwig phi sigm+mean file 
  TFile* f_phi_ref = new TFile( "data/recoTruthDphi_pp_mc_herwig.root", "read" );
  TH2D* hAngularResPhiTmp = static_cast< TH2D* >( f_phi_ref->Get( "h_recoTruthDphi_sigma" ) );
   hAngularResPhi->SetName( "hAngularResPhiTmp" );
  hAngularResPhiTmp->SetDirectory(0);
  f_phi_ref->Close();
  delete f_phi_ref;

  fOut->cd();
  
  TH2D* hDiffEta = static_cast< TH2D* >
    ( hAngularResEta->Clone( "hDiffEta" ) );
  styleTool->SetHStyle( hDiffEta, 0 );

  TH2D* hDiffPhi = static_cast< TH2D* >
    ( hAngularResPhi->Clone( "hDiffPhi" ) );
  styleTool->SetHStyle( hDiffPhi, 0 );

  // now subtract the pythia8 and the reference (herwig)
  hDiffEta->Add( hAngularResEtaTmp, -1 );
  hDiffPhi->Add( hAngularResPhiTmp, -1 );

  std::string yTitleMean;
  std::string yTitleSigma;
  GetTypeTitle( "Deta", yTitleMean, yTitleSigma );
  
  TLegend leg( 0.5, 0.45, 0.85, 0.75 );
  styleTool->SetLegendStyle( &leg );

  // ------------ Eta ------------
  TCanvas cEta( "cEta", "cEta", 800, 800 );
  TPad padEta1("padEta1", "", 0.0, 0.35, 1.0, 1.0 );
  padEta1.SetBottomMargin(0.0);
  padEta1.Draw();
	  
  TPad padEta2("padEta2", "", 0.0, 0.0, 1.0, 0.34 );
  padEta2.SetTopMargin(0.05);
  padEta2.SetBottomMargin(0.25);
  padEta2.Draw();

  padEta1.cd();

  std::string label = m_scaleResUseEta ? "#eta" : "#it{y}*";
  std::string name  = m_scaleResUseEta ? "Eta"  : "Ystar";
  
  std::vector< int > vXbins{ 1, 4 };
  TLine line( m_ptTruthMin, 0, m_ptTruthMax, 0 );
  line.SetLineWidth( 2 );

  int style = 0;
  for( auto& xBin: vXbins ){

    double xMin, xMax;
    anaTool->GetBinRange
      ( hAngularResEta->GetXaxis(), xBin, xBin, xMin, xMax );
    
    TH1* hProjAngResEta    = static_cast< TH1D* >
      ( hAngularResEta->ProjectionY( Form( "hAngularResEta_%d", xBin ), xBin, xBin ) );
    TH1* hProjAngResEtaTmp = static_cast< TH1D* >
      ( hAngularResEtaTmp->ProjectionY( Form( "hAngularResEtaTmp_%d", xBin ), xBin, xBin ) );
    vAngRes.push_back( hProjAngResEta );
    vAngRes.push_back( hProjAngResEtaTmp );

    styleTool->SetHStyle( hProjAngResEta, style );
    styleTool->SetHStyle( hProjAngResEtaTmp, style + 5 );

    leg.AddEntry( hProjAngResEta,
		  Form("Pythia8 %s", anaTool->GetLabel( xMin, xMax, label ).c_str() ) );
    leg.AddEntry( hProjAngResEtaTmp,
		  Form("Herwig  %s", anaTool->GetLabel( xMin, xMax, label ).c_str() ) );

    padEta1.cd();
    hProjAngResEta->SetYTitle( yTitleSigma.c_str() );
    hProjAngResEta->SetTitleOffset( 2.0, "y" );
      
    SetMinMax( hProjAngResEta, "Deta", "sigma" );
    
    hProjAngResEta->GetXaxis()->SetRangeUser( m_varPtBinning.front(), m_varPtBinning.back() );

    hProjAngResEta   ->Draw("ep same");
    hProjAngResEtaTmp->Draw("ep same");
    
    TH1* hDiff = static_cast< TH1D* >
      ( hDiffEta->ProjectionY( Form( "hDiffEta_%d", xBin ), xBin, xBin ) );
    vD.push_back( hDiff );
    styleTool->SetHStyleRatio( hDiff, style );

    hDiff->GetXaxis()->SetRangeUser( m_varPtBinning.front(), m_varPtBinning.back() );
      
    padEta2.cd();

    std::string yTitle = hAngularResEta->GetZaxis()->GetTitle();
    
    hDiff->SetYTitle( Form( "#Delta (%s)", yTitle.c_str() ) );
    hDiff->SetMaximum( 0.005 );
    hDiff->SetMinimum( -0.005 );
    hDiff->SetTitleOffset( 2.0, "y" );
    hDiff->Draw( "ep same" );

    style++;
  }

  padEta1.cd();

  DrawAtlasRight();
  leg.Draw();

  padEta2.cd();
  line.Draw();    
  
  SaveAsAll( cEta, Form( "h%s_angularResEta", name.c_str() ) );

  // ------------ Phi ------------
  GetTypeTitle( "Dphi", yTitleMean, yTitleSigma );

  TCanvas cPhi( "cPhi", "cPhi", 800, 800 );
  TPad padPhi1("padPhi1", "", 0.0, 0.35, 1.0, 1.0 );
  padPhi1.SetBottomMargin(0.0);
  padPhi1.Draw();
	  
  TPad padPhi2("padPhi2", "", 0.0, 0.0, 1.0, 0.34 );
  padPhi2.SetTopMargin(0.05);
  padPhi2.SetBottomMargin(0.25);
  padPhi2.Draw();

  padPhi1.cd();

  style = 0;
  
  for( auto& xBin: vXbins ){

    double xMin, xMax;
    anaTool->GetBinRange
      ( hAngularResPhi->GetXaxis(), xBin, xBin, xMin, xMax );
    
    TH1* hProjAngResPhi    = static_cast< TH1D* >
      ( hAngularResPhi->ProjectionY( Form( "hAngularResPhi_%d", xBin ), xBin, xBin ) );
    TH1* hProjAngResPhiTmp = static_cast< TH1D* >
      ( hAngularResPhiTmp->ProjectionY( Form( "hAngularResPhiTmp_%d", xBin ), xBin, xBin ) );
    vAngRes.push_back( hProjAngResPhi );
    vAngRes.push_back( hProjAngResPhiTmp );

    styleTool->SetHStyle( hProjAngResPhi, style );
    styleTool->SetHStyle( hProjAngResPhiTmp, style + 5 );

    padPhi1.cd();
    hProjAngResPhi->SetYTitle( yTitleSigma.c_str() );
    hProjAngResPhi->SetTitleOffset( 2.0, "y" );
      
    SetMinMax( hProjAngResPhi, "Dphi", "sigma" );

    hProjAngResPhi->GetXaxis()->SetRangeUser( m_varPtBinning.front(), m_varPtBinning.back() );
    
    hProjAngResPhi   ->Draw("ep same");
    hProjAngResPhiTmp->Draw("ep same");
    
    TH1* hDiff = static_cast< TH1D* >
      ( hDiffPhi->ProjectionY( Form( "hDiffPhi_%d", xBin ), xBin, xBin ) );
    vD.push_back( hDiff );
    styleTool->SetHStyleRatio( hDiff, style );

    hDiff->GetXaxis()->SetRangeUser( m_varPtBinning.front(), m_varPtBinning.back() );
    
    padPhi2.cd();

    std::string yTitle = hAngularResPhi->GetZaxis()->GetTitle();
    
    hDiff->SetYTitle( Form( "#Delta (%s)", yTitle.c_str() ) );
    hDiff->SetMaximum( 0.005 );
    hDiff->SetMinimum( -0.005 );
    hDiff->SetTitleOffset( 2.0, "y" );
    hDiff->Draw( "ep same" );

    style++;
  }

  padPhi1.cd();

  DrawAtlasRight();
  leg.Draw();

  padPhi2.cd();
  line.Draw();    
	     
  SaveAsAll( cPhi, Form( "h%s_angularResPhi", name.c_str() ) );
}

void DiJetAnalysisMC::CompareScaleRes( TFile* fOut, const std::string& type ){

  std::vector< TH1* > vH;
  
  // get the pPb mean/res file
  TFile* f_pPb = new TFile( Form( "data/%s_pPb_mc_pythia8.root", type.c_str() ), "read" );
  TH2D* hSigma_pPb = static_cast< TH2D* >( f_pPb->Get( Form( "h_%s_sigma", type.c_str() ) ) );
  hSigma_pPb->SetName( "hSigma_pPb" );
  hSigma_pPb->SetDirectory(0);
  TH2D* hMean_pPb = static_cast< TH2D* >( f_pPb->Get( Form( "h_%s_mean", type.c_str() ) ) );
  hMean_pPb->SetName( "hMean_pPb" );
  hMean_pPb->SetDirectory(0);
  f_pPb->Close();
  delete f_pPb;

  vH.push_back( hSigma_pPb );
  vH.push_back( hMean_pPb );

  // get the pPb_signal only mean/res file
  TFile* f_pPb_sig = new TFile( Form( "data/pPb_signal/%s_pPb_mc_pythia8.root", type.c_str() ), "read" );
  TH2D* hSigma_pPb_sig = static_cast< TH2D* >( f_pPb_sig->Get( Form( "h_%s_sigma", type.c_str() ) ) );
  hSigma_pPb_sig->SetName( "hSigma_pPb_sig" );
  hSigma_pPb_sig->SetDirectory(0);
  TH2D* hMean_pPb_sig = static_cast< TH2D* >( f_pPb_sig->Get( Form( "h_%s_mean", type.c_str() ) ) );
  hMean_pPb_sig->SetName( "hMean_pPb_sig" );
  hMean_pPb_sig->SetDirectory(0);
  f_pPb_sig->Close();
  delete f_pPb_sig;

  vH.push_back( hSigma_pPb_sig );
  vH.push_back( hMean_pPb_sig );

  // get the pPb_signal only mean/res file
  TFile* f_pPb_sig_new = new TFile( Form( "data/pPb_new_signal/%s_pPb_mc_pythia8.root", type.c_str() ), "read" );
  TH2D* hSigma_pPb_sig_new = static_cast< TH2D* >( f_pPb_sig_new->Get( Form( "h_%s_sigma", type.c_str() ) ) );
  hSigma_pPb_sig_new->SetName( "hSigma_pPb_sig_new" );
  hSigma_pPb_sig_new->SetDirectory(0);
  TH2D* hMean_pPb_sig_new = static_cast< TH2D* >( f_pPb_sig_new->Get( Form( "h_%s_mean", type.c_str() ) ) );
  hMean_pPb_sig_new->SetName( "hMean_pPb_sig_new" );
  hMean_pPb_sig_new->SetDirectory(0);
  f_pPb_sig_new->Close();
  delete f_pPb_sig_new;

  vH.push_back( hSigma_pPb_sig_new );
  vH.push_back( hMean_pPb_sig_new );
  
  // get the pp mean/res file
  TFile* f_pp = new TFile( Form( "data/%s_pp_mc_pythia8.root", type.c_str() ), "read" );
  TH2D* hSigma_pp = static_cast< TH2D* >( f_pp->Get( Form( "h_%s_sigma", type.c_str() ) ) );
  hSigma_pp->SetName( "hSigma_pp" );
  hSigma_pp->SetDirectory(0);
  TH2D* hMean_pp = static_cast< TH2D* >( f_pp->Get( Form( "h_%s_mean", type.c_str() ) ) );
  hMean_pp->SetName( "hMean_pp" );
  hMean_pp->SetDirectory(0);
  f_pp->Close();
  delete f_pp;

  vH.push_back( hSigma_pp );
  vH.push_back( hMean_pp );
  
  fOut->cd();
  
  std::string yTitleMean;
  std::string yTitleSigma;
  GetTypeTitle( type, yTitleMean, yTitleSigma );
  
  std::vector< int > vXbins{ 1, 2, 3, 4, 5 };
  TLine line0( m_ptTruthMin, 0, m_ptTruthMax, 0 );
  line0.SetLineWidth( 2 );
  TLine line1( m_ptTruthMin, 1, m_ptTruthMax, 1 );
  line1.SetLineWidth( 2 );

  std::string label = m_scaleResUseEta ? "#eta" : "#it{y}*";
  std::string name  = m_scaleResUseEta ? "Eta"  : "Ystar";

  // ------------------- Mean --------------------
  TH2* hDeltaMeanAll1 = static_cast< TH2D* >
    ( hMean_pPb->Clone( Form("h_delta_mean_%s_1", type.c_str() ) ) );
  TH2* hDeltaMeanAll2 = static_cast< TH2D* >
    ( hMean_pPb->Clone( Form("h_delta_mean_%s_2", type.c_str() ) ) );

  hDeltaMeanAll1->Reset();
  hDeltaMeanAll1->Reset();

  hDeltaMeanAll1->SetDirectory(0);
  hDeltaMeanAll1->SetDirectory(0);
  
  for( auto& xBin: vXbins ){
    TLegend leg( 0.5, 0.58, 0.85, 0.75 );
    styleTool->SetLegendStyle( &leg, 0.9 );

    TCanvas cMean( "cMean", "cMean", 800, 800 );
    TPad padMean1("padMean1", "", 0.0, 0.35, 1.0, 1.0 );
    padMean1.SetBottomMargin(0.0);
    padMean1.Draw();
	  
    TPad padMean2("padMean2", "", 0.0, 0.0, 1.0, 0.34 );
    padMean2.SetTopMargin(0.05);
    padMean2.SetBottomMargin(0.25);
    padMean2.Draw();

    padMean1.cd();

    double xMin, xMax;
    anaTool->GetBinRange
      ( hMean_pPb->GetXaxis(), xBin, xBin, xMin, xMax );

    std::string hTag = anaTool->GetName( xMin, xMax, name );
    
    TH1* hProjMean_pPb     = static_cast< TH1D* >
      ( hMean_pPb->ProjectionY( Form( "hMean_pPb_%d", xBin ), xBin, xBin ) );
    TH1* hProjMean_pPb_sig = static_cast< TH1D* >
      ( hMean_pPb_sig->ProjectionY( Form( "hMean_pPb_sig_%d", xBin ), xBin, xBin ) );
    TH1* hProjMean_pPb_sig_new = static_cast< TH1D* >
      ( hMean_pPb_sig_new->ProjectionY( Form( "hMean_pPb_sig_new_%d", xBin ), xBin, xBin ) );
    TH1* hProjMean_pp = static_cast< TH1D* >
      ( hMean_pp->ProjectionY ( Form( "hMean_pp_%d" , xBin ), xBin, xBin ) );

    vH.push_back( hProjMean_pPb );
    vH.push_back( hProjMean_pPb_sig );
    vH.push_back( hProjMean_pPb_sig_new );
    vH.push_back( hProjMean_pp  );
    
    styleTool->SetHStyle( hProjMean_pPb, 0 );
    styleTool->SetHStyle( hProjMean_pPb_sig, 1 );
    styleTool->SetHStyle( hProjMean_pPb_sig_new, 3 );
    styleTool->SetHStyle( hProjMean_pp , 2 );

    leg.AddEntry( hProjMean_pPb, "#it{p}+Pb w/Overlay" );
    leg.AddEntry( hProjMean_pPb_sig, "#it{p}+Pb Signal Old" );
    leg.AddEntry( hProjMean_pPb_sig_new, "#it{p}+Pb Signal New" );
    leg.AddEntry( hProjMean_pp, "#it{pp}" );

    padMean1.cd();
    hProjMean_pPb->SetYTitle( yTitleMean.c_str() );
    hProjMean_pPb->SetTitleOffset( 2.0, "y" );
      
    SetMinMax( hProjMean_pPb, type, "mean" );

    hProjMean_pPb->SetMaximum( 1.07 );
    hProjMean_pPb->SetMinimum( 0.96 );
    
    hProjMean_pPb->GetXaxis()->
      SetRangeUser( m_varPtBinning.front(), m_varPtBinning.back() );
    
    hProjMean_pPb->Draw("ep same");
    hProjMean_pPb_sig->Draw("ep same");
    hProjMean_pPb_sig_new->Draw("ep same");
    hProjMean_pp ->Draw("ep same");

    drawTool->DrawLeftLatex( 0.18, 0.85, anaTool->GetLabel( xMin, xMax, label ) );
  
    DrawAtlasRight();
    leg.Draw();

    if( type.find("Rpt") == std::string::npos ){
      line0.Draw();
    } else {
      line1.Draw();
    }
    
    padMean2.cd();

    TH1* hD1 = static_cast< TH1D* >
      ( hProjMean_pPb_sig->Clone( Form( "hMeanRatio1_%d", xBin ) ) );
    styleTool->SetHStyleRatio( hD1, 5 );
    hD1->Add( hProjMean_pp, -1);
    
    TH1* hD2 = static_cast< TH1D* >
      ( hProjMean_pPb_sig->Clone( Form( "hMeanRatio2_%d", xBin ) ) );
    styleTool->SetHStyleRatio( hD2, 6 );
    hD2->Add( hProjMean_pPb_sig_new, -1);

    TH1* hD3 = static_cast< TH1D* >
      ( hProjMean_pPb->Clone( Form( "hMeanRatio3_%d", xBin ) ) );
    styleTool->SetHStyleRatio( hD3, 7 );
    hD3->Add( hProjMean_pPb_sig, -1 );

    vH.push_back( hD1 );
    vH.push_back( hD2 );
    vH.push_back( hD3 );

    hD1->GetXaxis()->SetRangeUser( m_varPtBinning.front(), m_varPtBinning.back() );    

    std::string yTitle = hMean_pPb->GetZaxis()->GetTitle();

    hD1->SetYTitle( "Difference" );
    hD1->SetMaximum( 0.065 );
    hD1->SetMinimum( -0.06 );

    // for the angles, Deta Dphi, set different range
    if( type.find("Rpt") == std::string::npos ){
      hD1->SetMaximum( 0.01 );
      hD1->SetMinimum( -0.005 );
    }

    hD1->SetTitleOffset( 2.0, "y" );
    hD1->Draw( "ep same" );
    hD2->Draw( "ep same" );
    hD3->Draw( "ep same" );

    TLegend legR( 0.2, 0.85, 0.9, 0.95 );
    styleTool->SetLegendStyle( &legR, 0.7 );
    legR.SetNColumns(3);
    legR.AddEntry( hD3, "#color[4]{overlay - signal}" );
    legR.AddEntry( hD2, "signal_{old} - signal_{new}" );
    legR.AddEntry( hD1, "signal - pp" );
 
    legR.Draw();
    
    line0.Draw();   
    
    SaveAsAll( cMean, Form( "h_meanComp_%s_%s", type.c_str(), hTag.c_str() ) );

    // go through and fill all histos
    for( int yBin = 1; yBin <= hDeltaMeanAll1->GetNbinsY(); yBin++ ){
      hDeltaMeanAll1->SetBinContent( xBin, yBin, hD3->GetBinContent(yBin) );
      hDeltaMeanAll2->SetBinContent( xBin, yBin, hD1->GetBinContent(yBin) );
      hDeltaMeanAll1->SetBinError( xBin, yBin, hD3->GetBinError(yBin) );
      hDeltaMeanAll2->SetBinError( xBin, yBin, hD1->GetBinError(yBin) );
    }
  }

  // ------------------- Sigma --------------------
  for( auto& xBin: vXbins ){
    TLegend leg( 0.5, 0.60, 0.85, 0.75 );
    styleTool->SetLegendStyle( &leg );

    TCanvas cSigma( "cSigma", "cSigma", 800, 800 );
    TPad padSigma1("padSigma1", "", 0.0, 0.35, 1.0, 1.0 );
    padSigma1.SetBottomMargin(0.0);
    padSigma1.Draw();
	  
    TPad padSigma2("padSigma2", "", 0.0, 0.0, 1.0, 0.34 );
    padSigma2.SetTopMargin(0.05);
    padSigma2.SetBottomMargin(0.25);
    padSigma2.Draw();

    padSigma1.cd();

    double xMin, xMax;
    anaTool->GetBinRange
      ( hSigma_pPb->GetXaxis(), xBin, xBin, xMin, xMax );

    std::string hTag = anaTool->GetName( xMin, xMax, name );
    
    TH1* hProjSigma_pPb     = static_cast< TH1D* >
      ( hSigma_pPb->ProjectionY( Form( "hSigma_pPb_%d", xBin ), xBin, xBin ) );
    TH1* hProjSigma_pPb_sig = static_cast< TH1D* >
      ( hSigma_pPb_sig->ProjectionY( Form( "hSigma_pPb_sig_%d", xBin ), xBin, xBin ) );
    TH1* hProjSigma_pp = static_cast< TH1D* >
      ( hSigma_pp->ProjectionY ( Form( "hSigma_pp_%d" , xBin ), xBin, xBin ) );

    vH.push_back( hProjSigma_pPb );
    vH.push_back( hProjSigma_pPb_sig );
    vH.push_back( hProjSigma_pp  );
    
    styleTool->SetHStyle( hProjSigma_pPb, 0 );
    styleTool->SetHStyle( hProjSigma_pPb_sig, 1 );
    styleTool->SetHStyle( hProjSigma_pp , 2 );

    leg.AddEntry( hProjSigma_pPb, "#it{p}+Pb w/Overlay" );
    leg.AddEntry( hProjSigma_pPb_sig, "#it{p}+Pb Signal" );
    leg.AddEntry( hProjSigma_pp, "#it{pp}" );

    padSigma1.cd();
    hProjSigma_pPb->SetYTitle( yTitleSigma.c_str() );
    hProjSigma_pPb->SetTitleOffset( 2.0, "y" );
      
    SetMinMax( hProjSigma_pPb, type, "sigma" );

    hProjSigma_pPb->GetXaxis()->
      SetRangeUser( m_varPtBinning.front(), m_varPtBinning.back() );
    
    hProjSigma_pPb->Draw("ep same");
    hProjSigma_pPb_sig->Draw("ep same");
    hProjSigma_pp ->Draw("ep same");

    drawTool->DrawLeftLatex( 0.18, 0.85, anaTool->GetLabel( xMin, xMax, label ) );
  
    DrawAtlasRight();
    leg.Draw();

    padSigma2.cd();
 
    TH1* hD1 = static_cast< TH1D* >
      ( hProjSigma_pPb_sig->Clone( Form( "hSigmaRatio1_%d", xBin ) ) );
    styleTool->SetHStyleRatio( hD1, 5 );
    hD1->Add( hProjSigma_pp, -1);

    TH1* hD2 = static_cast< TH1D* >
      ( hProjSigma_pPb->Clone( Form( "hSigmaRatio2_%d", xBin ) ) );
    styleTool->SetHStyleRatio( hD2, 6 );
    hD2->Add( hProjSigma_pp, -1);
    
    TH1* hD3 = static_cast< TH1D* >
      ( hProjSigma_pPb->Clone( Form( "hSigmaRatio3_%d", xBin ) ) );
    styleTool->SetHStyleRatio( hD3, 7 );
    hD3->Add( hProjSigma_pPb_sig, -1 );


    vH.push_back( hD1 );
    // vH.push_back( hD2 );
    vH.push_back( hD3 );

    hD1->GetXaxis()->SetRangeUser( m_varPtBinning.front(), m_varPtBinning.back() );    

    std::string yTitle = hSigma_pPb->GetZaxis()->GetTitle();

    hD1->SetYTitle( "Difference" );
    hD1->SetMaximum( 0.05 );
    hD1->SetMinimum( -0.026 );

    // for the angles, Deta Dphi, set different range
    if( type.find("Rpt") == std::string::npos ){
      hD1->SetMaximum( 0.01 );
      hD1->SetMinimum( -0.01 );
    }

    hD1->SetTitleOffset( 2.0, "y" );
    hD1->Draw( "ep same" );
    hD2->Draw( "ep same" );
    hD3->Draw( "ep same" );

    TLegend legR( 0.2, 0.85, 0.9, 0.95 );
    styleTool->SetLegendStyle( &legR, 0.7 );
    legR.SetNColumns(3);

    legR.AddEntry( hD3, "#color[4]{overlay - signal}" );
    legR.AddEntry( hD2, "#color[2]{overlay - signal}" );
    legR.AddEntry( hD1, "signal - pp" );

    legR.Draw();
    
    line0.Draw();   
    
    SaveAsAll( cSigma, Form( "h_sigmaComp_%s_%s", type.c_str(), hTag.c_str() ) );
  }

  for( auto& h : vH ){ delete h; }

  // dont write this for pp
  // dont write if its in y* in pPb either
  if( !m_is_pPb || !m_scaleResUseEta ){ return; }
  
  TFile* fMeanOut = new TFile
    ( Form("data/pPb_delta_mean%s_%s.root", name.c_str(), type.c_str() ), "recreate");
  hDeltaMeanAll1->Write();
  hDeltaMeanAll2->Write();
  fMeanOut->Close();

  delete hDeltaMeanAll1;
  delete hDeltaMeanAll2;
}

void DiJetAnalysisMC::CompareRtrk( TFile* fOut ){

  TFile* fMCov = TFile::Open("data/rtrk/myOut_pPb_mc_pythia8_perf_0.overlay.root");
  TFile* fMCsi = TFile::Open("data/rtrk/myOut_pPb_mc_pythia8_perf_0.signal.root");
  TFile* fData = TFile::Open("data/rtrk/myOut_pPb_data_perf_0.root");

  TH1* hMCov  = static_cast< TH2D* >( fMCov->Get("h_rtrk1") );
  TH1* hMCsi  = static_cast< TH2D* >( fMCsi->Get("h_rtrk4") );
  TH1* hData1 = static_cast< TH2D* >( fData->Get("h_rtrk1") );
  TH1* hData2 = static_cast< TH2D* >( fData->Get("h_rtrk4") );

  TH2* hRov = static_cast< TH2D* >( hMCov->Clone("hRov") );
  TH2* hRsi = static_cast< TH2D* >( hMCsi->Clone("hRov") );

  hRov->Divide( hData1 );
  hRsi->Divide( hData2 );

  for( int xBin = 1; xBin <= hRov->GetNbinsX(); xBin++ ){

    double xLow  = hRov->GetXaxis()->GetBinLowEdge( xBin );
    double xHigh = hRov->GetXaxis()->GetBinUpEdge ( xBin );
    
    TH1* hOv = static_cast<TH1D*>
      ( hRov->ProjectionY( Form( "h_rtrk_ov_%d", xBin ), xBin, xBin ) );    
    styleTool->SetHStyleRatio( hOv, 0 );
    TH1* hSi = static_cast<TH1D*>
      ( hRsi->ProjectionY( Form( "h_rtrk_si_%d", xBin ), xBin, xBin ) );    
    styleTool->SetHStyleRatio( hSi, 1 );

    hSi->GetYaxis()->SetTitle( "#it{r}_{trk}^{MC}/#it{r}_{trk}^{Data}" );
    hSi->GetXaxis()->SetTitleOffset( 1.5 );
    
    TCanvas c("c","c", 800, 600 );

    double maximum = 1.15;
    double minimum = 0.91;

    if( xBin == 1 ){
      maximum = 1.30;
      minimum = 0.75;
    }
    
    hSi->SetMaximum( maximum );
    hSi->SetMinimum( minimum );
    
    hSi->Draw("ep same");
    hOv->Draw("ep same");

    double scale = 0.9;

    TLegend leg( 0.2, 0.2, 0.4, 0.4 );
    styleTool->SetLegendStyle( &leg );

    leg.AddEntry( hOv, "MC w/Overlay #it{p}_{T}^{trk} > 1 GeV " );
    leg.AddEntry( hSi, "MC Signal #it{p}_{T}^{trk} > 4 GeV " );

    leg.Draw();
    
    drawTool->DrawAtlasInternal( CT::DrawTools::drawX0, 0.87, 1.0 );
    drawTool->DrawRightLatex( CT::DrawTools::drawX0, 0.80, drawTool->GetLumipPb(), scale );
    drawTool->DrawAtlasOverlayInfo ( CT::DrawTools::drawX0, 0.725, scale );
    drawTool->DrawAtlasJetInfo( 0.46, 0.87, 2, scale);
    drawTool->DrawAtlasEnergy ( 0.46, 0.80, m_is_pPb, scale );

    TLine line( 20, 1, 90, 1 );
    line.Draw();
    
    drawTool->DrawLeftLatex( 0.20, 0.72, anaTool->GetLabel( xLow, xHigh, "#eta" ) );

    SaveAsPdfPng( c, Form( "h_rtrk_%d", xBin ) );
  }  
}


// CLEAN THIS UP! ITS NOT WORTH THE HEADACHE NOW
// BUT STILL, NOT GOOD (03.22.18)
void DiJetAnalysisMC::CompareCfactorsWUW( TFile* fOut ){

  std::vector< TH1* > vC;
  std::vector< TH1* > vR;
  
  TFile* fW  =
    TFile::Open( m_fNamePhysUF.c_str() );
  TFile* fUW = m_is_pPb ?
    TFile::Open( "data/output_pPb_mc_pythia8_uw/myOut_pPb_mc_pythia8_phys_UF_0.root" ) : 
    TFile::Open( "data/output_pp_mc_pythia8_uw/myOut_pp_mc_pythia8_phys_UF_0.root" );
  
  fOut->cd();
  
  TAxis* axis0 = m_dPP->GetTAxis(0); int nAxis0Bins = axis0->GetNbins();
  TAxis* axis1 = m_dPP->GetTAxis(1); int nAxis1Bins = axis1->GetNbins();
  TAxis* axis2 = m_dPP->GetTAxis(2); int nAxis2Bins = axis2->GetNbins();
  TAxis* axis3 = m_dPP->GetTAxis(3); int nAxis3Bins = axis3->GetNbins();

  // lines to be drawn along axis3. this is
  // x-axis that widths are plotted as function of 
  double xMin = m_dPhiZoomLow; double xMax = m_dPhiZoomHigh;
  
  TLine line( xMin, 1, xMax, 1 );
  line.SetLineWidth( 2 );

  TLine lineP05( xMin, 1.05, xMax, 1.05 );
  lineP05.SetLineStyle( 2  );
  lineP05.SetLineColor( 12 );
  lineP05.SetLineWidth( 1  );
	  
  TLine lineN05( xMin, 0.95, xMax, 0.95 );
  lineN05.SetLineStyle( 2  );
  lineN05.SetLineColor( 12 );
  lineN05.SetLineWidth( 1  );
	
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

  xMin = m_varYstarBinning.front();
  xMax = m_varYstarBinning.back();
  
  TLine lineY( xMin, 1, xMax, 1 );
  lineY.SetLineWidth( 2 );

  TLine lineP05Y( xMin, 1.05, xMax, 1.05 );
  lineP05Y.SetLineStyle( 2  );
  lineP05Y.SetLineColor( 12 );
  lineP05Y.SetLineWidth( 1  );
	  
  TLine lineN05Y( xMin, 0.95, xMax, 0.95 );
  lineN05Y.SetLineStyle( 2  );
  lineN05Y.SetLineColor( 12 );
  lineN05Y.SetLineWidth( 1  );
	
  TLine lineP25Y( xMin, 1.25, xMax, 1.25 );
  lineP25Y.SetLineStyle( 2  );
  lineP25Y.SetLineColor( 12 );
  lineP25Y.SetLineWidth( 2  );
  
  TLine lineN25Y( xMin, 0.75, xMax, 0.75 );
  lineN25Y.SetLineStyle( 2  );
  lineN25Y.SetLineColor( 12 );
  lineN25Y.SetLineWidth( 2  );

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

      for( int axis2Bin = 1; axis2Bin <= nAxis2Bins; axis2Bin++ ){

	// check we are in correct ystar and pt bins
	if( !m_dPP->CorrectPhaseSpace
	    ( std::vector<int>{ axis0Bin, axis1Bin, axis2Bin, 0 } ) )
	  { continue; }
	
	double axis2Low , axis2Up;
	anaTool->GetBinRange
	  ( axis2, axis2Bin, axis2Bin, axis2Low, axis2Up );
	
	std::string hTag =
	  anaTool->GetName( axis0Low, axis0Up, m_dPP->GetAxisName(0) ) + "_" + 
	  anaTool->GetName( axis1Low, axis1Up, m_dPP->GetAxisName(1) ) + "_" + 
	  anaTool->GetName( axis2Low, axis2Up, m_dPP->GetAxisName(2) ); 
	  
	std::string hName =
	  "h_" + m_dPhiRecoUnfoldedName + "_" + m_widthName + "_" + m_allName + "_" + hTag;

	TH1* hCw  = static_cast< TH1D* >( fW ->Get( hName.c_str() ) );
	TH1* hCuw = static_cast< TH1D* >( fUW->Get( hName.c_str() ) ) ;
	styleTool->SetHStyle( hCw , 0 );
	styleTool->SetHStyle( hCuw, 1 );
	vC.push_back( hCw  );
	vC.push_back( hCuw );
	  
	TLegend leg( 0.2, 0.1, 0.70, 0.2 );
	styleTool->SetLegendStyle( &leg );
	leg.SetNColumns(2);
	leg.AddEntry( hCw , "Weighted"  );
	leg.AddEntry( hCuw, "UnWeighted");
	  
	TCanvas c( "c", "c", 800, 800 );
	TPad pad1("pad1", "", 0.0, 0.35, 1.0, 1.0 );
	pad1.SetBottomMargin(0.0);
	pad1.Draw();
	  
	TPad pad2("pad2", "", 0.0, 0.0, 1.0, 0.34 );
	pad2.SetTopMargin(0.05);
	pad2.SetBottomMargin(0.25);
	pad2.Draw();

	pad1.cd();
	hCw ->Draw("ep X0 same" );
	hCuw->Draw("ep X0 same" );

	DrawTopLeftLabels
	  ( m_dPP, axis0Low, axis0Up, axis1Low, axis1Up,
	    axis2Low, axis2Up, 0, 0 );
	    
	DrawAtlasRight();

	leg.Draw();

	pad2.cd();

	std::string hNameR = hName + "_" + m_sRatio;

	TH1* hR = static_cast< TH1D* >( hCw->Clone( hNameR.c_str() ) );
	styleTool->SetHStyleRatio( hR );
	vR.push_back( hR );

	hR->Divide( hCuw );
	hR->SetYTitle( "Ratio" );

	hR->Draw( "hist p X0" );

	lineY.Draw();
	lineP25Y.Draw();
	lineN25Y.Draw();
	  
	SaveAsAll( c, Form("%s_WUW", hNameR.c_str() ) );
	
	for( int axis3Bin = 1; axis3Bin <= nAxis3Bins; axis3Bin++ ){
	  // check we are in correct ystar and pt bins
	  if( !m_dPP->CorrectPhaseSpace
	      ( std::vector<int>{ axis0Bin, axis1Bin, axis2Bin, axis3Bin } ) )
	    { continue; }

	  double axis3Low , axis3Up;
	  anaTool->GetBinRange
	    ( axis3, axis3Bin, axis3Bin, axis3Low, axis3Up );

	  std::string hTag =
	    anaTool->GetName( axis0Low, axis0Up, m_dPP->GetAxisName(0) ) + "_" + 
	    anaTool->GetName( axis1Low, axis1Up, m_dPP->GetAxisName(1) ) + "_" + 
	    anaTool->GetName( axis2Low, axis2Up, m_dPP->GetAxisName(2) ) + "_" +  
	    anaTool->GetName( axis3Low, axis3Up, m_dPP->GetAxisName(3) ); 
	  
	  std::string hName = "h_" + m_dPhiCfactorsName + "_" + m_allName + "_" + hTag;

	  TH1* hCw  = static_cast< TH1D* >( fW ->Get( hName.c_str() ) );
	  TH1* hCuw = static_cast< TH1D* >( fUW->Get( hName.c_str() ) ) ;
	  styleTool->SetHStyle( hCw , 0 );
	  styleTool->SetHStyle( hCuw, 1 );
	  vC.push_back( hCw  );
	  vC.push_back( hCuw );
	  
	  TLegend leg( 0.2, 0.1, 0.70, 0.2 );
	  styleTool->SetLegendStyle( &leg );
	  leg.SetNColumns(2);
	  leg.AddEntry( hCw , "Weighted"  );
	  leg.AddEntry( hCuw, "UnWeighted");
	  
	  TCanvas c( "c", "c", 800, 800 );
	  TPad pad1("pad1", "", 0.0, 0.35, 1.0, 1.0 );
	  pad1.SetBottomMargin(0.0);
	  pad1.Draw();
	  
	  TPad pad2("pad2", "", 0.0, 0.0, 1.0, 0.34 );
	  pad2.SetTopMargin(0.05);
	  pad2.SetBottomMargin(0.25);
	  pad2.Draw();

	  pad1.cd();
	  hCw ->Draw("ep X0 same" );
	  hCuw->Draw("ep X0 same" );

	  leg.Draw();

	  line.Draw();
	  lineP25Shift.Draw();
	  lineN25.Draw();
	  
	  DrawTopLeftLabels
	    ( m_dPP, axis0Low, axis0Up, axis1Low, axis1Up,
	      axis2Low, axis2Up, axis3Low, axis3Up );
	    
	  DrawAtlasRight();
	  
	  pad2.cd();

	  std::string hNameR = hName + "_" + m_sRatio;

	  TH1* hR = static_cast< TH1D* >( hCw->Clone( hNameR.c_str() ) );
	  styleTool->SetHStyleRatio( hR );
	  vR.push_back( hR );
	  
	  hR->SetMaximum( 1.2 );
	  hR->SetMinimum( 0.8 );

	  hR->Divide( hCuw );
	  hR->SetYTitle( "Ratio" );

	  hR->Draw( "hist p X0" );

	  line.Draw();
	  lineP05.Draw();
	  lineN05.Draw();
	  
	  SaveAsAll( c, hNameR );
	}
      }
    }
  }

  for( auto& c : vC ){ delete c; }
  for( auto& r : vR ){ delete r; }
}

//---------------------------
//        Drawing
//---------------------------
void DiJetAnalysisMC::DrawCanvas( std::vector< TH1* >& vHIN,
				  const std::string& type1,
				  const std::string& type2 ){
  TCanvas c("c","c",800,600);

  double lx0 = 0.20;
  double ly0 = 0.65;
  double lx1 = 0.40;
  double ly1 = 0.88;
  
  if( type2.find("sigma") != std::string::npos ){
    lx0 = 0.36;
    ly0 = 0.61;
    lx1 = 0.85;
    ly1 = 0.82;

    if( m_is_pPb ){
      lx0 = 0.34;
      lx1 = 0.83;
    }
  
  }

  // trying to finish my phd.
  // just hard code this part.
  // final paper only has 4 MC plots
  if( finalPlots  && type2.find("sigma") != std::string::npos ){
    lx0 = 0.64;
    ly0 = 0.575;
    lx1 = 0.84;
    ly1 = 0.885;
  } else if( finalPlots  && type2.find("mean") != std::string::npos ){
    lx0 = 0.19;
    ly0 = 0.575;
    lx1 = 0.39;
    ly1 = 0.885;
  }
  
  TLegend leg( lx0, ly0, lx1, ly1 );
  styleTool->SetLegendStyle( &leg );

  int style = 0;

  // for situations where dont want to
  // plot every single bin 
  // plot every n on canvas
  for( uint xRange = 0; xRange < vHIN.size(); xRange++ ){
    TH1* h = vHIN[ xRange ];
    styleTool->SetHStyle( h, style++ );
    // hardcoded shit
    if( xRange <= 2){
      leg.AddEntry( h, h->GetTitle() );
    } else {
      leg.AddEntry( h, Form( " %s", h->GetTitle() ) );
    }
    h->SetTitle("");
    // h->SetMarkerSize( h->GetMarkerSize() * 1.3 );
    h->Draw("epsame");
    SetMinMax( h, type1, type2 );
    h->Write();
  }

  leg.Draw();
  
  double y0 = GetLineHeight( type1 );
  
  TLine line( m_ptTruthMin, y0, m_ptTruthMax, y0);
  line.SetLineWidth(2);
  if( type2.find("mean") != std::string::npos ){ line.Draw(); }

  double scale = 0.9;
  
  if( finalPlots && type2.find("sigma") != std::string::npos ){
    drawTool->DrawAtlasSimulationInternal( 0.64, 0.27, 1.0 );
    if( m_is_pPb ){
      drawTool->DrawAtlasJetInfo     ( 0.60, 0.85, m_is_pPb, scale );
      drawTool->DrawAtlasOverlayInfo ( 0.60, 0.785, scale );
      drawTool->DrawAtlasEnergy      ( 0.43, 0.36, m_is_pPb, scale );
    } else {
      drawTool->DrawAtlasJetInfo     ( 0.60, 0.85, m_is_pPb, scale );
      drawTool->DrawAtlasEnergy      ( 0.43, 0.36, m_is_pPb,  scale );
    }
  } else if( finalPlots && type2.find("mean") != std::string::npos ){
    if( m_is_pPb ){
      drawTool->DrawAtlasSimulationInternal( 0.875, 0.855 );
      drawTool->DrawAtlasJetInfo     ( 0.875, 0.78, m_is_pPb, scale );
      drawTool->DrawAtlasOverlayInfo ( 0.875, 0.71, scale );
      drawTool->DrawAtlasEnergy      ( 0.875, 0.64, m_is_pPb, scale );
    } else {
      drawTool->DrawAtlasSimulationInternal( 0.875, 0.855 );
      drawTool->DrawAtlasJetInfo     ( 0.875, 0.78, m_is_pPb, scale );
      drawTool->DrawAtlasEnergy      ( 0.875, 0.71, m_is_pPb, scale );
    }
  } else {
    DrawAtlasRight();
  }
  
  SaveAsAll( c, Form("%s_%s", type1.c_str(), type2.c_str() ) );
}

//===== MinMax and line drawing =====

void DiJetAnalysisMC::
SetMinMax( TH1* h1, const std::string& type1, const std::string& type2 ){
  // JES JER
  if( type1.find("Rpt") != std::string::npos ){ 
    if( type2.find("mean") != std::string::npos ){ // sigma
      h1->SetMaximum(1.10);
      h1->SetMinimum(0.97);
    } else if( !type2.compare("sigma") ){ // sigma
      h1->SetMaximum(0.34);
      if( finalPlots ){ h1->SetMaximum( 0.20 ); }
      h1->SetMinimum(0.);
      if( finalPlots ){ h1->SetMinimum( 0.051 ); }
    }
  }
  // ANGLES
  else if( type1.find("Deta") != std::string::npos ||
	   type1.find("Dphi") != std::string::npos ) { 
    if( !type2.compare("mean") ){ // mean
      h1->SetMaximum(0.030);      
      h1->SetMinimum(-0.020);
      if( m_is_pPb && type1.find("Deta") != std::string::npos){
	h1->SetMinimum(-0.010);
      }
    } else if( !type2.compare("sigma") ){ // sigma
      h1->SetMaximum(0.07);
      if( finalPlots ){ h1->SetMaximum( 0.050 ); }
      h1->SetMinimum(0.);
      if( finalPlots ){ h1->SetMinimum( 0.005 ); }
      if( m_is_pPb && type1.find("Deta") != std::string::npos){
	h1->SetMaximum( 0.065 );
      }
    } 
  } 
}

double DiJetAnalysisMC::GetLineHeight( const std::string& type ){
  double y0 = 0;
  
  if( type.find("Rpt") != std::string::npos ){ // JES/JER
    y0 = 1;
    y0 = 1;
  } else if( type.find("Deta") != std::string::npos ||
	     type.find("Dphi") != std::string::npos ) { // ANGLES
    y0 = 0;
    y0 = 0;
  }

  return y0;
}

void DiJetAnalysisMC::DrawAtlasRight( double x0, double y0, double scale )
{ drawTool->DrawAtlasInternalMCRight( x0, y0, m_mcTypeLabel, m_is_pPb, scale ); } 


void DiJetAnalysisMC::DrawAtlasRightBoth( double x0, double y0, double scale )
{ drawTool->DrawAtlasInternalMCRight( x0, y0, m_mcTypeLabel, m_is_pPb, scale); }
