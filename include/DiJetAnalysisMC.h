#ifndef DIJETANALYSISMC_H
#define DIJETANALYSISMC_H

#include "THmulf.h"

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
  //       Fill Tree
  //---------------------------
  void RunOverTreeFillHistos( int, int );

  void SetupHistograms();

  void ProcessEvents( int, int );  

  void PairJets( std::vector< TLorentzVector >&,
		 std::vector< TLorentzVector >&,
		 std::vector< JetPair >& );
  
  //---------------------------
  //       Plotting 
  //---------------------------
  void ProcessPlotHistos();

  void LoadHistograms();
  
  void PlotEtaPhiPtMap( std::map< std::string, TH2* >& );

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
  
  //---------------------------
  //          Drawing
  //---------------------------
  void DrawCanvas( std::map< std::string, TH1* >&, TH1*,
		   double, double, 
		   const std::string&,
		   const std::string& );
  
  void DrawCanvas( std::vector< TH1* >&,
		   const std::string&,
		   const std::string&,
		   int = 1);

  void DrawCanvas( std::vector< TH1* >&,
		   const std::string&,
		   const std::string&,
		   bool );

  void DrawCanvas( std::vector< TGraphAsymmErrors* >&,
		   const std::string&,
		   const std::string&,
		   double, double );


  
  //===== MinMax and line drawing =====
  void SetMinMax( TH1*,
		  const std::string&,
		  const std::string& );

  double GetLineHeight( const std::string& );
  
 private:
  //============== cuts ===============
  double m_dRmax;
  
  //============ settings =============
  int m_mcType;
  std::string m_mcTypeLabel;
  
  std::vector< std::string > m_vUsedJZN;
  std::map< std::string, std::string > m_mJznFnameIn;

  std::map< std::string, double >  m_mJznPtThreshold;
  std::map< std::string, double >  m_mJznSigma;
  std::map< std::string, double >  m_mJznEff;
  std::map< std::string, double >  m_mJznSumPowhegWeights;
 
  std::map< std::string, int    >  m_mJznNev;
  std::map< std::string, double >  m_mJznWeights;

  double m_sumSigmaEff;

  TH3* m_hPowhegWeights;
  //============ data =============
  // These are for all ETA, PHI
  std::map< std::string, TH2* > m_mJznEtaPhiMap;
  std::map< std::string, TH2* > m_mJznEtaPtMap;

  std::map< std::string, TH2* > m_mJznEtaSpectReco;
  std::map< std::string, TH2* > m_mJznEtaSpectTruth;
  std::map< std::string, TH2* > m_mJznEtaSpectTruthNent; 
  std::map< std::string, TH2* > m_mJznEtaSpectTruthPaired;
  
  std::map< std::string, TH3* > m_mJznRecoTruthRpt;
  std::map< std::string, TH2* > m_mJznRecoTruthRptNent;

  std::map< std::string, TH3* > m_mJznRecoTruthDeta;
  std::map< std::string, TH2* > m_mJznRecoTruthDetaNent;
  
  std::map< std::string, TH3* > m_mJznRecoTruthDphi;
  std::map< std::string, TH2* > m_mJznRecoTruthDphiNent;

  //========= histos binning ========
  // truth bins
  double m_ptTruthWidth;
  double m_ptTruthMin;
  double m_ptTruthMax;
  int    m_nPtTruthBins;

  // Jes Jer
  int    m_nRPtBins;
  double m_rPtMin;
  double m_rPtMax;

  // angular bins
  int    m_nDAngleBins;
  double m_dAngleMin;
  double m_dAngleMax;
};


#endif
