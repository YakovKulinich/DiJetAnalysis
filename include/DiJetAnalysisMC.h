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

  void PlotSpectra( std::map< int, TH2* >&,
		    const std::string&, 
		    const std::string& );
  
  void PlotEtaPhiPtMap( std::map< int, TH2* >& );

  void PlotVsEtaPt( std::map< int, TH3* >&,
		    std::map< int, TH2* >&,
		    const std::string& );

  //---------------------------
  //          Tools 
  //---------------------------
  void ProjectAndFit( TH3*, TH1*, TH1*, int, int, int );

  void CombineJZN( TH1*,
		   std::map< int, TH1*>& );

  void CombineJZN( TH1*,
		   std::map< int, TH1*>&,
		   std::map< int, TH1*>& );
  
  double GetJetWeight( double, double, double );
  
  //---------------------------
  //          Drawing
  //---------------------------
  void DrawCanvas( std::map< int, TH1* >&, TH1*,
		   double, double, 
		   const std::string&,
		   const std::string& );

  void DrawCanvas( std::vector< TH1* >&,
		   const std::string&,
		   const std::string&,
		   bool );
  
  void DrawCanvas( std::vector< TH1* >&,
		   const std::string&,
		   const std::string& );


  //===== MinMax and line drawing =====
  void SetMinMax( TH1*,
		  const std::string&,
		  const std::string& );

  double GetLineHeight( const std::string& );
  
 private:
  //============== cuts ===============
  double m_dRmax;
  double m_ptFitMin;
  int    m_nMinEntriesGausFit;
  
  //============ settings =============
  int m_mcType;
  std::string m_mcTypeLabel;
  
  std::vector< int > m_vUsedJZN;
  std::map< int, std::string > m_mJznFnameIn;

  std::map< int, double >      m_mJznSigma;
  std::map< int, double >      m_mJznEff;
  std::map< int, double >      m_mJznSumPowhegWeights;
 
  std::map< int, int    >      m_mJznNev;
  std::map< int, double >      m_mJznWeights;

  double m_sumSigmaEff;

  TH3* m_hPowhegWeights;
  //============ data =============
  // These are for all ETA, PHI
  std::map< int, TH2* > m_mJznEtaPhiMap;
  std::map< int, TH2* > m_mJznEtaPtMap;

  std::map< int, TH2* > m_mJznEtaSpectReco;
  std::map< int, TH2* > m_mJznEtaSpectTruth;
  
  std::map< int, TH3* > m_mJznRpt;
  std::map< int, TH3* > m_mJznDeta;
  std::map< int, TH3* > m_mJznDphi;
  std::map< int, TH2* > m_mJznNentries;

  //============ histos =============
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
