#ifndef DIJETANALYSISDATA_H
#define DIJETANALYSISDATA_H

#include "DiJetAnalysis.h"

class DiJetAnalysisData : public DiJetAnalysis{
 public:
  DiJetAnalysisData();

  DiJetAnalysisData( bool );

  DiJetAnalysisData( bool, int );

  DiJetAnalysisData( bool, int, int );

  ~DiJetAnalysisData();

  //---------------------------
  // Initialization Methods
  //---------------------------
  void Initialize();

  void AdditionalSuffix( std::string& ){}
  
  //---------------------------
  // Fill Tree / Plot Controls
  //---------------------------
  void RunOverTreeFillHistos( int, int );

  void ProcessPlotHistos();

  void DataMCCorrections();
  
  void PlotHistosTogether();

  void ProcessSystematics();
  
  //---------------------------
  //       Fill Tree
  //---------------------------
  void SetupHistograms();

  void ProcessEvents( int, int );  

  void LoadTriggers();
  
  //---------------------------
  //       Analysis
  //---------------------------
  bool JetInTrigPtRange  ( const TLorentzVector&, int,
			   double = 0);
  
  bool JetInTrigEtaRange ( const TLorentzVector&, int );
  
  bool JetInTrigRange    ( const TLorentzVector&, int );

  bool TrigJetAboveThold ( const TLorentzVector&, int );
  
  bool TrigJetInTrigRange( const TLorentzVector&, int );
  
  void AnalyzeEff( std::vector< TLorentzVector >&,
		   std::vector< TLorentzVector >&,
		   std::map< int, bool >&);
  
  //---------------------------
  //       Tools
  //---------------------------  
  void CleanEfficiency( TGraphAsymmErrors*, int );

  TH1*       CombineSamples( std::vector< TH1* >&,
			     const std::string& = "" );

  TH2*       CombineSamples( std::vector< TH2* >&,
			     const std::string& = "" );
  
  THnSparse* CombineSamples( std::vector< THnSparse* >&,
			     const std::string& = "" );   

  void GetInfoBoth( std::string&, std::string&, std::string&,
		    std::string&, std::string&, std::string& );

  void GetInfoUnfolding( std::string&, std::string&,
			 std::string&, std::string& );
  
  //---------------------------
  //  Get Quantities / Plot 
  //---------------------------
  void LoadHistograms();

  void MakeEfficiencies( std::vector< TH2* >&,
			 std::vector< TH2* >&,
			 const std::string& );

  void MakeSystematicsGraphs();
  
  //---------------------------
  //        Drawing
  //---------------------------
  void DrawAtlasRight( double = 0, double = 0, double = CT::StyleTools::lSS );
  
  void DrawAtlasRightBoth( double = 0, double = 0, double = CT::StyleTools::lSS ); 

 private:
  //============ settings =============
  std::string m_fNameIn;

  std::vector< std::string > m_vTriggers;
  std::vector< int    >      m_vTriggersRefIndex;
  std::vector< double >      m_vTriggersPrescale;
  std::vector< double >      m_vTriggersTholdPt;
  std::vector< double >      m_vTriggersEffPtLow;
  std::vector< double >      m_vTriggersEffPtHigh;
  std::vector< double >      m_vTriggersEtaMin;
  std::vector< double >      m_vTriggersEtaMax;

  uint m_nTriggers;
  
  std::string m_mbTriggerName;

  int m_mbTriggerI;
  int m_lowestCentTriggerI;  

  double m_centMbCorrection;

  //============ data =============
  // -------- maps ---------
  std::vector< TH2* > m_vHtriggerEtaPhiMap;
  std::vector< TH2* > m_vHtriggerEtaPtMap;

  // -------- spect --------  
  std::vector< TH2* > m_vHtriggerEtaSpect;

  TH2* m_hAllEtaSpect;
  // ----- efficiencies ----
  std::vector< TH2* > m_vHtriggerEtaSpectSim;
  std::vector< TH2* > m_vHtriggerEtaSpectDenom;
  
  // -------- dPhi ---------
  
  std::vector< THnSparse* > m_vHtriggerDphi;
  std::vector< THnSparse* > m_vHtriggerDphiNent;

  THnSparse* m_hAllDphi;
  THnSparse* m_hAllDphiNent;

  THnSparse* m_hAllDphiUnfolded;
  //========= histos binning ========
};

#endif
