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
  // Fill Tree / Plot Controls
  //---------------------------
  void RunOverTreeFillHistos( int, int );

  void ProcessPlotHistos();

  //---------------------------
  //       Fill Tree
  //---------------------------
  void LoadTriggers();

  void SetupHistograms();

  void ProcessEvents( int, int );  

  //---------------------------
  //       Analysis
  //---------------------------
  void AnalyzeEff( std::vector< TLorentzVector >&,
		   std::vector< TLorentzVector >&,
		   std::map< int, bool >&);

  //---------------------------
  //       Plot Data 
  //---------------------------
  void LoadHistograms();

  void PlotSpectra( std::vector< TH2* >&,
		    const std::string& );

  void PlotEfficiencies( std::vector< TH2* >&,
			 std::vector< TH2* >&,
			 const std::string& );

  virtual void PlotDeltaPhi
    ( std::vector< THnSparse*>&,
      const std::string& );

  void PlotEtaPhiPtMap( std::vector< TH2* >& );
  
  virtual void PlotDataTogether();

   
 private:
  //============ settings =============
  std::string m_fNameIn;

  std::vector< std::string > m_vTriggers;
  std::vector< int    >      m_vRefTriggerIndex;
  std::vector< double >      m_vTholdPtTriggers;
  std::vector< double >      m_vEffPtTriggers;

  std::string m_mbTriggerName;
  int m_mbTriggerI;
  
  uint m_nTriggers;
 
  //============ data =============
  // -------- maps ---------
  std::vector< TH2* > m_vTriggerEtaPhiMap;
  std::vector< TH2* > m_vTriggerEtaPtMap;

  // -------- spect --------
  std::vector< TH2* > m_vTriggerEtaSpect;
  std::vector< TH2* > m_vTriggerEtaSpectSim;

  // -------- dPhi ---------
  std::vector< THnSparse* > m_vDphi;
  
  //========= histos binning ========
};

#endif
