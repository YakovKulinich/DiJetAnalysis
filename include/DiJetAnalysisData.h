#ifndef DIJETANALYSISDATA_H
#define DIJETANALYSISDATA_H

#include <TH2.h>

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

  void loadTriggers();

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

  void plotSpectra( int, int );

  void plotEfficiencies( int, int );

  void plotEtaPhi();

  void plotEtaPt();
  
 private:
  //============ data =============
  std::vector< std::string > v_triggers;

  std::map< std::string, TH2* > m_triggerSpectEta;
  std::map< std::string, TH2* > m_triggerEffEta;
  
  std::map< std::string, TH2* > m_triggerRunPrescale;
  std::map< std::string, TH2* > m_triggerEtaPhi;
  std::map< std::string, TH2* > m_triggerEtaPt;

  std::vector< int > v_tJetPt;
  std::map   < int , std::string > m_tJetPtTrigger;

  std::string  mbTrigger;
};

#endif
