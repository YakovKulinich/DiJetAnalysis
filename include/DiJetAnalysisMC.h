#ifndef DIJETANALYSISMC_H
#define DIJETANALYSISMC_H

#include "THmulf.h"

#include "DiJetAnalysis.h"

class JetPair;
class TLegend;

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

  void PlotSpectra( std::map< int, TH2* >& );
  
  void PlotEtaPhiPtMap( std::map< int, TH2* >& );

  void PlotVsEtaPt( std::map< int, TH3* >&,
		    std::map< int, TH2* >&,
		    int );
  
  void ProjectEtaPtAndFit( TH3*, TH1*, TH1*, int, int );

  //---------------------------
  //          Tools 
  //---------------------------
  void CombineJZN( TH1*,
		   std::map< int, TH1*>& );

  void CombineJZN( TH1*,
		   std::map< int, TH1*>&,
		   std::map< int, TH1*>& );
    
  void DrawCanvas( std::map< int, TH3* >&,
		   TCanvas&, TLegend&,
		   double, double,
		   int, bool );

  void DrawCanvas( std::vector< TH1* >&,
		   TCanvas&, TLegend&,
		   int, bool );
  
 private:
  //============== cuts ===============
  double m_dRmax;
  
  //============ settings =============
  std::vector< int > m_vUsedJZN;
  std::map< int, std::string > m_mJznFnameIn;
  std::map< int, double >      m_mJznSigma;
  std::map< int, double >      m_mJznEff;
  std::map< int, int    >      m_mJznNev;
  std::map< int, double >      m_mJznWeights;

  double m_sumSigmaEff;
  //============ data =============
  std::map< int, TH2* > m_mJznEtaSpect;
  std::map< int, TH2* > m_mJznEtaPhi;
  std::map< int, TH2* > m_mJznEtaPt;

  std::map< int, TH3* > m_mJznRpt;
  std::map< int, TH3* > m_mJznDeta;
  std::map< int, TH3* > m_mJznDphi;
  std::map< int, TH3* > m_mJznRecoEff;
  std::map< int, TH2* > m_mJznNentries;  
};

#endif
