#ifndef DIJETANALYSISMC_H
#define DIJETANALYSISMC_H

#include "DiJetAnalysis.h"

class TLegend;
class UncertaintyProvider;

class DiJetAnalysisMC : public DiJetAnalysis{
 public:
  DiJetAnalysisMC();

  DiJetAnalysisMC( bool );

  DiJetAnalysisMC( bool, int );

  DiJetAnalysisMC( bool, int, int );

  ~DiJetAnalysisMC();

  //---------------------------
  // Initialization Methods
  //---------------------------
  void Initialize();

  void AdditionalSuffix( std::string& );

  std::string GetMCMenu();
  
  //---------------------------
  // Fill Tree / Plot Controls
  //---------------------------
  void RunOverTreeFillSpect( int, int );

  void ProcessSpectWeights();
  
  void RunOverTreeFillDphi( int, int );

  void ProcessDphiWeights();
  
  void RunOverTreeFillHistos( int, int );

  void ProcessPerformance   ();

  void UnfoldPerformance    ();

  void ProcessPhysics       ();

  void UnfoldPhysics        ();

  //---------------------------
  //       Fill Tree
  //---------------------------
  void SetupHistograms();

  void LoadSpectWeights();

  void LoadDphiWeights();
  
  void ProcessEventsForWeights( int, int, int = 1 );
   
  void ProcessEvents( int, int );

 protected:

  //---------------------------
  //       Analysis
  //---------------------------
  double AnalyzeDeltaPhiWithWeight
    ( THnSparse*,
      const std::vector <TLorentzVector >&,
      const std::vector <TLorentzVector >&,
      int = 0,
      bool = true);

  double AnalyzeDeltaPhiWithDphiWeight
    ( THnSparse*,
      const std::vector <TLorentzVector >&,
      const std::vector <TLorentzVector >& );
  
  void AnalyzeScaleResolution( const std::vector< TLorentzVector >&,
			       const std::vector< TLorentzVector >&,
			       const int );

  void AnalyzeYstarRespMat( TH3*,
			    const std::vector< TLorentzVector >&,
			    const std::vector< TLorentzVector >& );

  void AnalyzeSpectRespMat( TH3*,
			    const std::vector< TLorentzVector >&,
			    const std::vector< TLorentzVector >& );

  void AnalyzeDphiRespMat( THnSparse*, THnSparse*,
			   const std::vector< TLorentzVector >&,
			   const std::vector< TLorentzVector >& );
  
  //---------------------------
  //          Tools 
  //---------------------------  
  void PairJets( std::vector< TLorentzVector >&,
		 std::vector< TLorentzVector >&,
		 std::vector< TLorentzVector >&,
		 std::vector< TLorentzVector >& );

  double GetSpectWeight( const TLorentzVector& );

  double GetDphiWeight( const TLorentzVector&, const TLorentzVector& );
  
  TH1*       CombineSamples( std::vector< TH1* >&,
			     const std::string& = "" );

  TH2*       CombineSamples( std::vector< TH2* >&,
			     const std::string& = "" );

  TH3*       CombineSamples( std::vector< TH3* >&,
			     const std::string& = "" );
  
  THnSparse* CombineSamples( std::vector< THnSparse* >&,
			     const std::string& = "" );    
  
  void CombineSamples( TH1*,
		       std::vector< TH1* >&,
		       std::vector< TH1* >&,
		       const std::string& = "" );

  void CombineSamples( TH1*,
		       std::vector< TH1* >&,
		       std::vector< TH1* >&,
		       std::vector< TH1* >&,
		       std::vector< TH2* >&,
		       const std::string& = "" );

  TGraphAsymmErrors* CombineSamples( std::vector< TGraphAsymmErrors* >&,
				     std::vector< TH1* >&  );
  
  void SetCfactorsErrors( TH1*, TH1*, TH2*, TH1* );
  
  double GetJetWeight( const TLorentzVector& );
  
  void GetTypeTitle( const std::string&,
		     std::string&, std::string& );

  void GetSpectWeightInfo( std::string&, std::string&,
			   std::string&, std::string& );

  void GetDphiWeightInfo ( std::string&, std::string&,
			   std::string&, std::string& );
  
  void GetSpectUnfoldingInfo( std::string&, std::string&, std::string&,
			      std::string&, std::string&, std::string& );;
  
  void GetDphiUnfoldingInfo ( std::string&, std::string&,
			      std::string&, std::string& );

  void GetInfoTogether( std::string&, std::string&, std::string&,
			std::string&, std::string&, std::string&,
			int = 0 );

  void GetPurityEff( TH2*, TH1*, TH1* );
  
  //---------------------------
  //  Get Quantities / Plot 
  //---------------------------
  void LoadHistograms( int = 0 );

  TH2* MakeSpectWeights( TFile* = NULL );

  TH3* MakeDphiWeights( TFile* = NULL );

  void MakeScaleRes( std::vector< TH3* >&,
		     std::vector< TH2* >&,
		     const std::string& );

  void MakeEfficiencies( std::vector< TH2* >&,
			 std::vector< TH2* >&,
			 const std::string& );
  
  void MakeSpectCFactorsRespMat( std::vector< TH2* >&,
				 std::vector< TH2* >&,
				 std::vector< TH3* >&,
				 const std::vector< std::string >&,
				 const std::string& = "",
				 const std::string& = "" );
  
  void MakeDphiCFactorsRespMat( std::vector< THnSparse* >&,
				std::vector< THnSparse* >&,
				std::vector< THnSparse* >&,
				const std::vector< std::string >&,
				const std::string& = "" ,
				const std::string& = "" ); 

  void MakePtRespMat( std::vector< THnSparse* >&,
		      const std::vector< std::string >&,
		      const std::string& = "" );

  void MakeYstarRespMat( std::vector< TH3* >&,
			 const std::vector< std::string >&,
			 const std::string& = "" );

  void CompareAngularRes ( TFile* = NULL );

  void CompareRtrk       ( TFile* = NULL ); 
  
  void CompareScaleRes   ( TFile* = NULL, const std::string& = "" );

  void CompareCfactorsWUW( TFile* = NULL );
  
  //---------------------------
  //        Drawing
  //---------------------------
  void DrawCanvas( std::vector< TH1* >&,
		   const std::string& = "",
		   const std::string& = "" );
		 
  void DrawAtlasRight    ( double = CT::DrawTools::drawX0,
			   double = CT::DrawTools::drawY0,
			   double = CT::StyleTools::lSS );

  void DrawAtlasRightBoth( double = CT::DrawTools::drawX0,
			   double = CT::DrawTools::drawY0,
			   double = CT::StyleTools::lSS ); 
  
  //===== MinMax and line drawing =====
  void SetMinMax( TH1*,
		  const std::string&,
		  const std::string& );

  double GetLineHeight( const std::string& );
    
 private:
  //============== cuts ===============
  double m_dRmax;

  //============== tool ===============
  UncertaintyProvider* m_uncertaintyProvider;

  TH1* h_pPbFCalWeights;
  double m_FCalEt;
  //============ settings ============= 
  std::vector< int > m_vJznUsed;
  
  std::vector< std::string > m_vJznLabels;
  std::vector< std::string > m_vJznFnameIn;

  std::vector< double >  m_vJznSigma;
  std::vector< double >  m_vJznEff;
  std::vector< double >  m_vJznSumOverlayWeights;
  std::vector< double >  m_vJznPtThreshold;
 
  std::vector< int    >  m_vJznNev;
  std::vector< double >  m_vJznWeights;

  uint m_nJzn;
    
  double m_sumSigmaEff;

  TH3* m_hPowhegWeights;

  std::string m_fNamePerfWeightData;
  std::string m_fNamePhysWeightData;

  TH2* m_spectWeight;
  TH3* m_dPhiWeight;
  std::vector< TF1* > m_vSpectWeightFits;
  
  //============ data =============
  // -------- maps ---------
  std::vector< TH2* > m_vHjznEtaPhiMap;
  std::vector< TH2* > m_vHjznEtaPtMap;

  TH2* m_hAllEtaPtMap;
  TH2* m_hAllEtaPhiMap;
  
  // -------- spect --------
  std::string m_ystarSpectFineTruthUPName;
 
  std::vector< TH2* > m_vHjznYstarSpectReco;
  std::vector< TH2* > m_vHjznYstarSpectTruth;
  std::vector< TH2* > m_vHjznYstarSpectFineReco;
  std::vector< TH2* > m_vHjznYstarSpectFineTruth;
  std::vector< TH2* > m_vHjznYstarSpectFineTruthUP;
  
  TH2* m_hAllYstarSpectReco;
  TH2* m_hAllYstarSpectTruth;
  TH2* m_hAllYstarSpectFineReco;
  TH2* m_hAllYstarSpectFineTruth;

  // --- spectra response matrix ----
  std::vector< TH3* > m_vHjznYstarSpectRespMat;
  std::vector< TH3* > m_vHjznYstarRespMat;

  TH3* m_hAllYstarSpectRespMat;
  TH3* m_hAllYstarRespMat;
  
  // --------- recoTruthRpt ---------
  std::vector< TH3* > m_vHjznRecoTruthRpt;
  std::vector< TH2* > m_vHjznRecoTruthRptNent;

  // --------- recoTruthDeta ---------
  std::vector< TH3* > m_vHjznRecoTruthDeta;
  std::vector< TH2* > m_vHjznRecoTruthDetaNent;

  // --------- recoTruthDphi ---------
  std::vector< TH3* > m_vHjznRecoTruthDphi;
  std::vector< TH2* > m_vHjznRecoTruthDphiNent;

  // ----- rtrk ----
  std::vector< TH3* > m_vHjznRtrk1;
  std::vector< TH3* > m_vHjznRtrk2;
  std::vector< TH3* > m_vHjznRtrk4;
  TH3* m_hAllJznRtrk1;
  TH3* m_hAllJznRtrk2;
  TH3* m_hAllJznRtrk4;
  
  // -------------- dPhi -------------
  std::vector< THnSparse* > m_vHjznDphiReco;
  std::vector< THnSparse* > m_vHjznDphiTruth;
  
  THnSparse* m_hAllDphiReco;
  THnSparse* m_hAllDphiTruth;
  
  // ----- dPhi response matrix ------
  std::string m_dPhiRespMatName;
  std::string m_dPhiRespMatRebName;

  std::string m_ptRespMatName;

  std::vector< THnSparse* > m_vHjznDphiRespMat;  
  THnSparse* m_hAllDphiRespMat;

  // -------- Rebinned Response Matrix --------     
  std::vector< THnSparse* > m_vHjznDphiRespMatReb;  
  THnSparse* m_hAllDphiRespMatReb;

  // ------- dPhi reco unfolded ------
  std::string m_dPhiRecoUnfoldedName;

  //========= histos binning ========
  
  // ------ truth binning --------
  double m_ptTruthWidth;
  double m_ptTruthMin;
  double m_ptTruthMax;
  int    m_nPtTruthBins;

  // ---- JES/PRes/Etc -----
  bool m_scaleResUseEta;

  // --- var range ----
  int    m_nRPtRecoTruthBins;
  double m_rPtRecoTruthMin;
  double m_rPtRecoTruthMax;

  int    m_nDAngleRecoTruthBins;
  double m_dAngleRecoTruthMin;
  double m_dAngleRecoTruthMax;
  
  // ---- dPhi Reponse Matrix ----- 
  uint m_nDphiRespMatDim;

  std::vector< int >    m_vNdPhiRespMatBins;
  std::vector< int >    m_vNdPhiRespMatRebBins;
  std::vector< double > m_vDphiRespMatMin;
  std::vector< double > m_vDphiRespMatMax; 
};

#endif
