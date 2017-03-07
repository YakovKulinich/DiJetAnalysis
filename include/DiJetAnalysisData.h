#ifndef DIJETANALYSISDATA_H
#define DIJETANALYSISDATA_H

#include "DiJetAnalysis.h"

class DiJetAnalysisData : public DiJetAnalysis{
 public:
  DiJetAnalysisData();
  DiJetAnalysisData( bool, bool );
  ~DiJetAnalysisData();

  void Initialize();
  
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

  void PlotSpectra( std::map< std::string, TH2* >& );

  void PlotEfficiencies( std::map< std::string, TH2* >& );

  void PlotEtaPhiPtMap( std::map< std::string, TH2* >& );
  
 private:
  //============ settings =============
  std::string m_fNameIn;
  
  //============ data =============
  std::vector< std::string > m_vTriggers;

  std::map< std::string, TH2* > m_mTriggerEtaSpect;
  std::map< std::string, TH2* > m_mTriggerEtaEff;
  
  std::map< std::string, TH2* > m_mTriggerRunPrescale;
  std::map< std::string, TH2* > m_mTriggerEtaPhi;
  std::map< std::string, TH2* > m_mTriggerEtaPt;

  std::vector< int > v_tJetPt;
  std::map   < int , std::string > m_tJetPtTrigger;

  std::string  mbTrigger;
};

#endif
