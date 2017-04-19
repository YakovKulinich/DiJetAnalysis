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

  void PairJets( std::vector< TLorentzVector >&,
		 std::vector< TLorentzVector >&,
		 std::vector< JetPair >& );
  
  //---------------------------
  //       Plot Data 
  //---------------------------
  void LoadHistograms();
  
  void PlotSpectra( std::map< std::string, TH2* >&,
		    const std::string&, 
		    const std::string& );
  
  void PlotVsEtaPt( std::map< std::string, TH3* >&,
		    std::map< std::string, TH2* >&,
		    const std::string& );

  void PlotEfficiencies( std::map< std::string, TH2* >&,
			 std::map< std::string, TH2* >&,
			 std::map< std::string, TH2* >&,
			 const std::string& );

  void PlotEtaPhiPtMap( std::map< std::string, TH2* >& );
  
  //---------------------------
  //          Tools 
  //---------------------------  
  void
    CombineJZN( TH1*,
		std::map< std::string, TH1*>& );
  
  void
    CombineJZN( TH1*,
		std::map< std::string, TH1*>&,
		std::map< std::string, TH1*>& );

  TGraphAsymmErrors*
    CombineJZN( std::map< std::string, TGraphAsymmErrors* >&,
		std::map< std::string, TH1*>& );
  
  double GetJetWeight( double, double, double );
  
 private:
  //============== cuts ===============
  double m_dRmax;
  
  //============ settings ============= 
  std::vector< std::string > m_vUsedJZN;
  std::map< std::string, std::string > m_mJznFnameIn;

  std::map< std::string, double >  m_mJznSigma;
  std::map< std::string, double >  m_mJznEff;
  std::map< std::string, double >  m_mJznSumPowhegWeights;
  std::map< std::string, double >  m_mJznPtThreshold;
 
  std::map< std::string, int    >  m_mJznNev;
  std::map< std::string, double >  m_mJznWeights;

  double m_sumSigmaEff;

  TH3* m_hPowhegWeights;
  //============ data =============
  // -------- maps ---------
  std::map< std::string, TH2* > m_mJznEtaPhiMap;
  std::map< std::string, TH2* > m_mJznEtaPtMap;

  // -------- spect --------
  std::map< std::string, TH2* > m_mJznEtaSpectReco;
  std::map< std::string, TH2* > m_mJznEtaSpectTruth;
  std::map< std::string, TH2* > m_mJznEtaSpectTruthNent; 
  std::map< std::string, TH2* > m_mJznEtaSpectTruthPaired;

  // --------- recoTruthRpt ---------
  std::map< std::string, TH3* > m_mJznRecoTruthRpt;
  std::map< std::string, TH2* > m_mJznRecoTruthRptNent;

  // --------- recoTruthDeta ---------
  std::map< std::string, TH3* > m_mJznRecoTruthDeta;
  std::map< std::string, TH2* > m_mJznRecoTruthDetaNent;

  // --------- recoTruthDphi ---------
  std::map< std::string, TH3* > m_mJznRecoTruthDphi;
  std::map< std::string, TH2* > m_mJznRecoTruthDphiNent;

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
