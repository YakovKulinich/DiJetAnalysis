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

  void PlotEtaPhiPtMap( std::map< std::string, TH2* >& );
  
  void PlotSpectra( std::map< std::string, TH2* >&,
		    const std::string& );

  void PlotEfficiencies( std::map< std::string, TH2* >&,
			 const std::string& );
  
 private:
  //============ settings =============
  std::string m_fNameIn;

  std::vector< int > v_tJetPt;
  std::map   < int , std::string > m_tJetPtTrigger;

  std::vector< std::string > m_vTriggers;
  std::string  m_mbTrigger;  
  //============ data =============
  // -------- maps ---------
  std::map< std::string, TH2* > m_mTriggerEtaPhiMap;
  std::map< std::string, TH2* > m_mTriggerEtaPtMap;

  // -------- spect --------
  std::map< std::string, TH2* > m_mTriggerEtaSpect;
  std::map< std::string, TH2* > m_mTriggerEtaSpectEff;

  // ------- runInfo ------
  std::map< std::string, TH2* > m_mTriggerRunPrescale;

  //========= histos binning ========
};

#endif
