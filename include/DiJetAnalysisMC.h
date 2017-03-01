#ifndef DIJETANALYSISMC_H
#define DIJETANALYSISMC_H

#include "THmulf.h"

#include "DiJetAnalysis.h"

class JetPair;

class DiJetAnalysisMC : public DiJetAnalysis{
 public:
  DiJetAnalysisMC();
  DiJetAnalysisMC( bool, bool );
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

  void PlotSpectra( int, int );
  
  void PlotEtaPhiPtMap( std::map< int, TH2* >& );

  void PlotVsEtaPt( int, int, std::map< int, TH3* >&, int );
  
  void ProjectEtaPtAndFit( TH3*, TH1*, TH1*, int, int );
  
  void FitGaussian( TH1*, TF1* );
  
 private:
  //============   cuts   =============
  double m_dRmax;
  
  //============ settings =============
  std::vector< int > v_usedJZN;
  std::map< int, std::string > m_jznFnameIn;
  std::map< int, double >      m_jznSigma;
  std::map< int, double >      m_jznEff;

  //============ data =============
  std::map< int, TH2* > m_jznEtaSpect;
  std::map< int, TH2* > m_jznEtaPhi;
  std::map< int, TH2* > m_jznEtaPt;

  std::map< int, TH3* > m_jznRPt;
  std::map< int, TH3* > m_jznDeta;
  std::map< int, TH3* > m_jznDphi;

};

#endif
