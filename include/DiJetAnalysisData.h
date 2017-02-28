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

  void LoadTriggers();

  void SetupHistograms();

  void ProcessEvents( int, int );  

  void ProcessEfficiencies( std::vector< TLorentzVector >&,
			    std::vector< TLorentzVector >&,
			    std::map< std::string, bool >&);

  //---------------------------
  //       Plotting 
  //---------------------------
  void ProcessPlotHistos();

  void LoadHistograms();

  void PlotSpectra( int, int );

  void PlotEfficiencies( int, int );

  void PlotEtaPhi();

  void PlotEtaPt();
  
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
