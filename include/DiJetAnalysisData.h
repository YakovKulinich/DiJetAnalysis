#ifndef DIJETANALYSISDATA_H
#define DIJETANALYSISDATA_H

#include "DiJetAnalysis.h"

class DiJetAnalysisData : public DiJetAnalysis{
 public:
  DiJetAnalysisData();
  DiJetAnalysisData( bool, bool );
  ~DiJetAnalysisData();
  
  //---------------------------
  //       Fill Tree
  //---------------------------
  void RunOverTreeFillHistos( int, int );

  void loadTriggers( bool );

  void setupHistograms();

  void processEvents( int, int );  

  void processEfficiencies( std::vector< TLorentzVector >&,
			    std::vector< TLorentzVector >&,
			    std::map< std::string, bool >&);

  //---------------------------
  //       Plotting 
  //---------------------------
  void PlotExistingHistos();

  void loadHistograms();

  void plotSpectra();

  void plotEfficiencies();

  void plotEtaPhi();

  void plotPtEta();
  
 private:
  //============ data =============
  std::vector< std::string > v_triggers;

  std::map< std::string, TH1* > m_triggerSpect;
  std::map< std::string, TH1* > m_triggerEff;
  std::map< std::string, TGraphAsymmErrors* > m_triggerEffGrf;

  std::map< std::string, TH1* > m_triggerRunPrescale;
  std::map< std::string, TH1* > m_triggerEtaPhi;
  std::map< std::string, TH1* > m_triggerPtEta;

  std::vector< int > v_tJetPt;
  std::multimap< int , std::string > m_tJetPtTrigger;

  std::string  mbTrigger;
};

#endif
