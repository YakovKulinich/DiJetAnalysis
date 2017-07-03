#ifndef DIJETANALYSISMC_H
#define DIJETANALYSISMC_H

#include "DiJetAnalysis.h"

class TLegend;

class DiJetAnalysisMC : public DiJetAnalysis{
 public:
  DiJetAnalysisMC();
  DiJetAnalysisMC( bool, int );
  ~DiJetAnalysisMC();

  void Initialize();

  //---------------------------
  // Fill Tree / Plot Controls
  //---------------------------
  void RunOverTreeFillHistos( int, int );

  void ProcessPlotHistos();

  void PlotHistosTogether();
  
  //---------------------------
  //       Fill Tree
  //---------------------------
  void SetupHistograms();

  void ProcessEvents( int, int );

  //---------------------------
  //       Analysis
  //---------------------------
  void AnalyzeScaleResolution( const std::vector< TLorentzVector >&,
			       const std::vector< TLorentzVector >&,
			       const int );

  void AnalyzeResponseMatrix( THnSparse*,
			      const std::vector<TLorentzVector>&,
			      const std::vector<TLorentzVector>&,
			      WeightFcn = NULL );

  std::pair<double,double> AnalyzeDeltaPhiTruthReco
    ( THnSparse*, THnSparse*,
      const std::vector <TLorentzVector>&,
      const std::vector <TLorentzVector>&,
      WeightFcn = NULL );

  //---------------------------
  //          Tools 
  //---------------------------  
  void PairJets( std::vector< TLorentzVector >&,
		 std::vector< TLorentzVector >&,
		 std::vector< TLorentzVector >&,
		 std::vector< TLorentzVector >& );

  TH1*       CombineSamples( std::vector< TH1* >&,
			     const std::string& = "" );

  TH2*       CombineSamples( std::vector< TH2* >&,
			     const std::string& = "" );
  
  THnSparse* CombineSamples( std::vector< THnSparse* >&,
			     const std::string& = "" );    
  
  void CombineSamples( TH1*,
		       std::vector< TH1*>&,
		       std::vector< TH1*>&,
		       const std::string& = "" );

  TGraphAsymmErrors* CombineJZN
    ( std::vector< TGraphAsymmErrors* >&,
      std::vector< TH1*>& );
  
  static double GetJetWeight( double, double, double );

  void GetTypeTitle( const std::string&,
		     std::string&, std::string& );

  void GetInfoBoth( std::string&, std::string&, std::string&, std::string&,
		    std::string&, std::string&, std::string& );

  void GetInfoBothRecoTruth( std::string&, std::string&,
			     std::string&, std::string& );
  
  void GetInfoUnfolding( std::string&, std::string& );

  void GetInfoBothRecoTruth();
  
  //---------------------------
  //  Get Quantities / Plot 
  //---------------------------
  void LoadHistograms();

  void MakeScaleRes( std::vector< TH3* >&,
		     std::vector< TH2* >&,
		     const std::string& );
  
  void MakeDphiRecoTruth();
  
  void MakeResponseMatrix( std::vector<THnSparse*>&,
			   const std::vector< std::string >&,
			   const std::string& = "" );
  
  //---------------------------
  //        Drawing
  //---------------------------
  void DrawCanvas( std::vector< TH1* >&,
		   const std::string& = "",
		   const std::string& = "",
		   int = 1);

  void DrawAtlasRight( double = 0, double = 0, double = CT::StyleTools::lSS );

  void DrawAtlasRightBoth( double = 0, double = 0, double = CT::StyleTools::lSS ); 
  
  //===== MinMax and line drawing =====
  void SetMinMax( TH1*,
		  const std::string&,
		  const std::string& );

  double GetLineHeight( const std::string& );
    
 private:
  //============== cuts ===============
  double m_dRmax;
  
  //============ settings ============= 
  std::string m_mcTypeLabel;

  std::vector< int > m_vJznUsed;
  
  std::vector< std::string > m_vJznLabel;
  std::vector< std::string > m_vJznFnameIn;

  std::vector< double >  m_vJznSigma;
  std::vector< double >  m_vJznEff;
  std::vector< double >  m_vJznSumOverlayWeights;
  std::vector< double >  m_vJznPtThreshold;
 
  std::vector< int    >  m_vJznNev;
  std::vector< double >  m_vJznWeights;

  uint m_nJzn;
    
  double m_sumSigmaEff;

  static TH3* m_hPowhegWeights;

  //============ data =============
  // -------- maps ---------
  std::vector< TH2* > m_vHjznEtaPhiMap;
  std::vector< TH2* > m_vHjznEtaPtMap;

  // -------- spect --------
  std::string m_etaSpectRecoName;
  std::string m_etaSpectTruthName; 

  std::vector< TH2* > m_vHjznEtaSpectReco;
  std::vector< TH2* > m_vHjznEtaSpectTruth;
  std::vector< TH2* > m_vHjznEtaSpectTruthPaired;

  TH2* m_hAllEtaSpectReco;
  TH2* m_hAllEtaSpectTruth;
  
  // --------- recoTruthRpt ---------
  std::vector< TH3* > m_vHjznRecoTruthRpt;
  std::vector< TH2* > m_vHjznRecoTruthRptNent;

  // --------- recoTruthDeta ---------
  std::vector< TH3* > m_vHjznRecoTruthDeta;
  std::vector< TH2* > m_vHjznRecoTruthDetaNent;

  // --------- recoTruthDphi ---------
  std::vector< TH3* > m_vHjznRecoTruthDphi;
  std::vector< TH2* > m_vHjznRecoTruthDphiNent;

  // -------------- dPhi -------------
  std::vector< THnSparse* > m_vHjznDphiReco;
  std::vector< THnSparse* > m_vHjznDphiTruth;
  
  THnSparse* m_hAllDphiReco;
  THnSparse* m_hAllDphiTruth;

  // --- dPhi truth reco together ----
  std::string m_dPhiRecoPtTruthName;
  std::string m_dPhiTruthPtRecoName;
  
  std::vector< THnSparse* > m_vHjznDphiRecoPtTruth;
  std::vector< THnSparse* > m_vHjznDphiTruthPtReco;

  THnSparse* m_hAllDphiRecoPtTruth;
  THnSparse* m_hAllDphiTruthPtReco;
  
  // -------- Dphi Response Matrix --------     
  std::vector< THnSparse* > m_vHjznDphiRespMat;
  
  THnSparse* m_hAllDphiRespMat;
  
  //========= histos binning ========
  // ------ truth binning --------
  double m_ptTruthWidth;
  double m_ptTruthMin;
  double m_ptTruthMax;
  int    m_nPtTruthBins;

  // ---- JES/PRes/Etc ----- 
  int    m_nRPtRecoTruthBins;
  double m_rPtRecoTruthMin;
  double m_rPtRecoTruthMax;

  int    m_nDAngleRecoTruthBins;
  double m_dAngleRecoTruthMin;
  double m_dAngleRecoTruthMax;

  // ---- dPhi Reponse Matrix ----- 
  uint m_nDphiRespMatDim;

  std::vector< int >    m_nDphiRespMatBins;
  std::vector< double > m_dPhiRespMatMin;
  std::vector< double > m_dPhiRespMatMax;
  
};

#endif
