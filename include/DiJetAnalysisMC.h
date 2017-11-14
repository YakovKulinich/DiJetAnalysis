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
  void RunOverTreeFillHistos( int, int );

  void ProcessPlotHistos();

  void DataMCCorrections();

  void PlotHistosTogether();

  void ProcessSystematics(){};
  
  //---------------------------
  //       Fill Tree
  //---------------------------
  void SetupHistograms();

  void ProcessEvents( int, int );

 protected:
  //---------------------------
  //       Analysis
  //---------------------------
  void AnalyzeScaleResolution( const std::vector< TLorentzVector >&,
			       const std::vector< TLorentzVector >&,
			       const int );

  void AnalyzeResponseMatrix( THnSparse*, THnSparse*,
			      const std::vector<TLorentzVector>&,
			      const std::vector<TLorentzVector>& );

  std::pair<double,double> AnalyzeDeltaPhiTruthReco
    ( THnSparse*, THnSparse*,
      const std::vector <TLorentzVector>&,
      const std::vector <TLorentzVector>& );

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
		       std::vector< TH1* >&,
		       std::vector< TH1* >&,
		       const std::string& = "" );

  void CombineSamples( TH1*,
		       std::vector< TH1* >&,
		       std::vector< TH1* >&,
		       std::vector< TH1* >&,
		       std::vector< TH2* >&,
		       const std::string& = "" );

  double GetJetWeight( const TLorentzVector& );
  
  double GetUncertaintyWeight( const TLorentzVector&,
			       const TLorentzVector& );

  void GetTypeTitle( const std::string&,
		     std::string&, std::string& );

  void GetInfoBoth( std::string&, std::string&, std::string&,
		    std::string&, std::string&, std::string& );

  void GetInfoBothRecoTruth( std::string&, std::string&,
			     std::string&, std::string& );
  
  void GetInfoUnfolding( std::string&, std::string&,
			 std::string&, std::string& );

  void AnalyzePurityEff( TH2*, TH1*, TH1* );
  
  //---------------------------
  //  Get Quantities / Plot 
  //---------------------------
  void LoadHistograms();

  void MakeScaleRes( std::vector< TH3* >&,
		     std::vector< TH2* >&,
		     const std::string& );
  
  void MakeDphiRecoTruth();
  
  void MakeDphiCFactorsRespMat( std::vector<THnSparse*>&,
				std::vector<THnSparse*>&,
				std::vector<THnSparse*>&,
				const std::vector< std::string >&,
				const std::string& = "" ,
				const std::string& = "" ); 

  void MakePtResponseMatrix( std::vector<THnSparse*>&,
			     const std::vector< std::string >&,
			     const std::string& = "" );
  
  //---------------------------
  //        Drawing
  //---------------------------
  void DrawCanvas( std::vector< TH1* >&,
		   const std::string& = "",
		   const std::string& = "" );
		 

  void DrawAtlasRight    ( double = 0, double = 0, double = CT::StyleTools::lSS );

  void DrawAtlasRightBoth( double = 0, double = 0, double = CT::StyleTools::lSS ); 
  
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
  
  //============ settings ============= 
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

  TH3* m_hPowhegWeights;

  //============ data =============
  // -------- maps ---------
  std::vector< TH2* > m_vHjznEtaPhiMap;
  std::vector< TH2* > m_vHjznEtaPtMap;

  // -------- spect --------
  std::string m_etaSpectRecoName;
  std::string m_etaSpectTruthName; 
  std::string m_spectTruthName;
  
  std::vector< TH2* > m_vHjznEtaSpectReco;
  std::vector< TH2* > m_vHjznEtaSpectTruth;
  std::vector< TH1* > m_vHjznSpectTruth;
  
  TH2* m_hAllEtaSpectReco;
  TH2* m_hAllEtaSpectTruth;
  TH1* m_hAllSpectTruth;
  
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
  std::vector< THnSparse* > m_vHjznDphiRecoPairedTruth;
  
  THnSparse* m_hAllDphiReco;
  THnSparse* m_hAllDphiTruth;
  THnSparse* m_hAllDphiRecoPairedTruth;

  // ------- response matrix ---------
  std::string m_dPhiRespMatName;
  std::string m_dPhiRespMatRebName;

  std::string m_ptRespMatName;

  // ----- dPhi reco paired truth ----
  std::string m_dPhiRecoPairedTruthName;
  
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

  // -------- Rebinned Response Matrix --------     
  std::vector< THnSparse* > m_vHjznDphiRespMatReb;  
  THnSparse* m_hAllDphiRespMatReb;

  //========= histos binning ========
  
  // ------ truth binning --------
  double m_ptTruthWidth;
  double m_ptTruthMin;
  double m_ptTruthMax;
  int    m_nPtTruthBins;

  // ---- JES/PRes/Etc -----
  // --- var range ----
  int    m_nRPtRecoTruthBins;
  double m_rPtRecoTruthMin;
  double m_rPtRecoTruthMax;

  int    m_nDAngleRecoTruthBins;
  double m_dAngleRecoTruthMin;
  double m_dAngleRecoTruthMax;

  // --- variable pt binning w/over+underflow ---
  std::vector< double > m_varPtBinningRespMat;
  uint m_nVarPtBinsRespMat;

  
  // ---- dPhi Reponse Matrix ----- 
  uint m_nDphiRespMatDim;

  std::vector< int >    m_vNdPhiRespMatBins;
  std::vector< int >    m_vNdPhiRespMatRebBins;
  std::vector< double > m_vDphiRespMatMin;
  std::vector< double > m_vDphiRespMatMax; 
};

#endif
