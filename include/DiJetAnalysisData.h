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
  void SetupHistograms();

  void ProcessEvents( int, int );  

  void LoadTriggers();
  
  //---------------------------
  //       Analysis
  //---------------------------
  bool IsInTriggerRange( double, uint );
  
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

  void PlotEtaPhiPtMap( std::vector< TH2* >& );
  
  virtual void PlotDphiTogether();

   
 private:
  //============ settings =============
  std::string m_fNameIn;

  std::vector< std::string > m_vTriggers;
  std::vector< int    >      m_vRefTriggerIndex;
  std::vector< double >      m_vTholdPtTriggers;
  std::vector< double >      m_vEffPtTriggers;

  std::string m_mbTriggerName;
  std::string m_allName;

  int m_mbTriggerI;
  
  uint m_nTriggers;
 
  //============ data =============
  // -------- maps ---------
  std::vector< TH2* > m_vHtriggerEtaPhiMap;
  std::vector< TH2* > m_vHtriggerEtaPtMap;

  // -------- spect --------
  std::vector< TH2* > m_vHtriggerEtaSpect;
  std::vector< TH2* > m_vHtriggerEtaSpectSim;

  TH2* m_hAllEtaSpect;
  // -------- dPhi ---------
  std::vector< THnSparse* > m_vHtriggerDphi;

  THnSparse* m_hAllDphi;
  //========= histos binning ========
};

#endif
