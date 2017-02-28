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
  
  //---------------------------
  //       Fill Tree
  //---------------------------
  void RunOverTreeFillHistos( int, int );

  void SetupHistograms();

  void ProcessEvents( int, int );  

  void PairJets( std::vector< TLorentzVector >&,
		 std::vector< TLorentzVector >&,
		 std::vector< JetPair >&,
		 std::vector< double >& );
  
  //---------------------------
  //       Plotting 
  //---------------------------
  void ProcessPlotHistos();

  void LoadHistograms();

  void PlotSpectra();

  void PlotEtaPhi();

  void PlotEtaPt();
  
 private:
  //============ data =============
  std::map< int, std::string > m_jznFname;
  std::map< int, double >      m_jznSigma;
  std::map< int, double >      m_jznEff;

  std::map< int, TH1* > m_jznRPt;
  std::map< int, TH1* > m_jznDeta;
  std::map< int, TH1* > m_jznDphi;

  // eta-phi or eta-pt maps of inclusive jets
  std::map< int, TH1* > m_jznEtaPhi;
  std::map< int, TH1* > m_jznEtaPt;
};

#endif
