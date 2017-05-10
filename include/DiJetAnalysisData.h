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
  bool IsInTriggerPtRange( double, int, double = 0);
  
  bool IsInTriggerEtaRange( double, int );
  
  bool IsInTriggerRange( TLorentzVector&, int );

  bool TriggerJetAboveThreshold( TLorentzVector&, int );
  
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

  void PlotDphiTogether();
  
  void PlotEtaPhiPtMap( std::vector< TH2* >& );
   
 private:
  //============ settings =============
  std::string m_fNameIn;

  std::vector< std::string > m_vTriggers;
  std::vector< int    >      m_vTriggersRefIndex;
  std::vector< double >      m_vTriggersTholdPt;
  std::vector< double >      m_vTriggersEffPtLow;
  std::vector< double >      m_vTriggersEffPtHigh;
  std::vector< double >      m_vTriggersEtaMin;
  std::vector< double >      m_vTriggersEtaMax;

  std::string m_mbTriggerName;
  std::string m_allName;

  int m_mbTriggerI;
  int m_lowestCentTriggerI;  

  double m_centMbCorrection;
  
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
  std::vector< THnSparse* > m_vHtriggerDphiNent;

  THnSparse* m_hAllDphi;
  THnSparse* m_hAllDphiNent;
  //========= histos binning ========
};

#endif
