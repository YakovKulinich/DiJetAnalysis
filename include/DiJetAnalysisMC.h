#ifndef DIJETANALYSISMC_H
#define DIJETANALYSISMC_H

#include "DiJetAnalysis.h"

class JetPair;
class TLegend;

class DiJetAnalysisMC : public DiJetAnalysis{
 public:
  DiJetAnalysisMC();
  DiJetAnalysisMC( bool, bool, int );
  ~DiJetAnalysisMC();

  void Initialize();

  //---------------------------
  // Fill Tree / Plot Controls
  //---------------------------
  void RunOverTreeFillHistos( int, int );

  void ProcessPlotHistos();

  //---------------------------
  //       Fill Tree
  //---------------------------
  void SetupHistograms();

  void ProcessEvents( int, int );  
  
  //---------------------------
  //       Plot Data 
  //---------------------------
  void LoadHistograms();
  
  void PlotSpectra( std::vector< TH2* >&,
		    const std::string&, 
		    const std::string& );
  
  void PlotVsEtaPt( std::vector< TH3* >&,
		    std::vector< TH2* >&,
		    const std::string& );

  void PlotEfficiencies( std::vector< TH2* >&,
			 std::vector< TH2* >&,
			 std::vector< TH2* >&,
			 const std::string& );

  void PlotDeltaPhi( std::vector< THnSparse*>&,
		     std::vector< THnSparse*>&,
		     const std::vector< std::string >&,
		     const std::string& = "" ,
		     const std::string& = "" );

  void PlotDphiTogether();

  void PlotCombinedDphiWidthsTogether();
  
  void PlotEtaPhiPtMap( std::vector< TH2* >& );
  
  //---------------------------
  //          Tools 
  //---------------------------  
  void PairJets( std::vector< TLorentzVector >&,
		 std::vector< TLorentzVector >&,
		 std::vector< JetPair >& );

  void CombineJZN( TH1*,
		   std::vector< TH1*>& );
  
  void CombineJZN( TH1*,
		   std::vector< TH1*>&,
		   std::vector< TH1*>& );

  TGraphAsymmErrors* CombineJZN
    ( std::vector< TGraphAsymmErrors* >&,
      std::vector< TH1*>& );
  
  static double GetJetWeight( double, double, double );

  void GetTypeTitle( const std::string&,
		     std::string&, std::string& );
  
 private:
  //============== cuts ===============
  double m_dRmax;
  
  //============ settings ============= 
  std::vector< int > m_vJznUsed;
  uint m_nJzn;
  
  std::vector< std::string > m_vJznLabel;
  std::vector< std::string > m_vJznFnameIn;

  std::vector< double >  m_vJznSigma;
  std::vector< double >  m_vJznEff;
  std::vector< double >  m_vJznSumOverlayWeights;
  std::vector< double >  m_vJznPtThreshold;
 
  std::vector< int    >  m_vJznNev;
  std::vector< double >  m_vJznWeights;

  double m_sumSigmaEff;

  static TH3* m_hPowhegWeights;
  //============ data =============
  // -------- maps ---------
  std::vector< TH2* > m_vHjznEtaPhiMap;
  std::vector< TH2* > m_vHjznEtaPtMap;

  // -------- spect --------
  std::vector< TH2* > m_vHjznEtaSpectReco;
  std::vector< TH2* > m_vHjznEtaSpectTruth;
  std::vector< TH2* > m_vHjznEtaSpectTruthNent; 
  std::vector< TH2* > m_vHjznEtaSpectTruthPaired;

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
  std::vector< THnSparse* > m_vHjznDphiRecoNent;
  std::vector< THnSparse* > m_vHjznDphiTruth;
  std::vector< THnSparse* > m_vHjznDphiTruthNent;
  
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
};

#endif
