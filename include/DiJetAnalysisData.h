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

  void ProcessPerformance   ();

  void UnfoldPerformance    ();

  void ProcessPhysics       ();

  void UnfoldPhysics        ();

  void ProcessSystematics   ();
  
  //---------------------------
  //       Fill Tree
  //---------------------------
  void SetupHistograms();

  void ProcessEvents( int, int );  

  void LoadTriggers();
  
  //---------------------------
  //       Analysis
  //---------------------------
  bool JetInTrigPtRange  ( const TLorentzVector*, int, double = 0 );
  
  bool JetInTrigEtaRange ( const TLorentzVector*, int );
  
  bool JetInTrigRange    ( const TLorentzVector*, int );

  bool TrigJetAboveThold ( const TLorentzVector*, int );
  
  bool TrigJetInTrigRange( const TLorentzVector*, int );
  
  void AnalyzeEff( std::vector< TLorentzVector >&,
		   std::vector< TLorentzVector >&,
		   std::map< int, bool >& );
  
  //---------------------------
  //       Tools
  //---------------------------  
  void CleanEfficiency( TGraphAsymmErrors*, int );

  TH1*       CombineSamples( std::vector< TH1* >&,
			     const std::string& = "" );

  TH2*       CombineSamples( std::vector< TH2* >&,
			     const std::string& = "" );

  TH3*       CombineSamples( std::vector< TH3* >&,
			     const std::string& = "",
			     bool = true );
  
  THnSparse* CombineSamples( std::vector< THnSparse* >&,
			     const std::string& = "" );   

  void GetSpectUnfoldingInfo( std::string&, std::string&, std::string&,
			      std::string&, std::string&, std::string& );

  void GetDphiUnfoldingInfo( std::string&, std::string&,
			     std::string&, std::string& );
  
  void GetInfoTogether( std::string&, std::string&, std::string&,
			std::string&, std::string&, std::string&,
			int = 0 );

  //---------------------------
  //  Get Quantities / Plot 
  //---------------------------
  void LoadHistograms( int = 0 );

  void MakeEfficiencies( std::vector< TH2* >&,
			 std::vector< TH2* >&,
			 const std::string& );

  void MakeNjetsRun( std::vector< TH3* >&,
		     const std::vector< std::string >&,
		     const std::string& );
    
  void MakeSystematicsGraphs( TFile* = NULL, const std::string& = "" );

  void MakeFinalPlotsTogether( TFile* fOut, const std::string& = "" );
  
  void CompareCfactorsRBnRB( TFile* = NULL );
  
  //---------------------------
  //        Drawing
  //---------------------------
  void DrawAtlasRight    ( double = CT::DrawTools::drawX0,
			   double = CT::DrawTools::drawY0,
			   double = CT::StyleTools::lSS );
  
  void DrawAtlasRightBoth( double = CT::DrawTools::drawX0,
			   double = CT::DrawTools::drawY0,
			   double = CT::StyleTools::lSS, bool = false );

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

  std::string m_nJetsRunName;
  
  std::map< int, int    > m_mRunBin;
  std::map< int, double > m_mRunLumi;

  uint m_nRuns;
  
  std::string m_fNamePerfUnfoldingMC;
  std::string m_fNamePhysUnfoldingMC;
  
  //============ data =============
  // -------- maps ---------
  std::vector< TH2* > m_vHtriggerEtaPhiMap;
  std::vector< TH2* > m_vHtriggerEtaPtMap;
  std::vector< TH3* > m_vHtriggerEtaPhiPtMap;

  TH2* m_hAllEtaPhiMap;
  TH2* m_hAllEtaPtMap;
  TH3* m_hAllEtaPhiPtMap;
  // -------- spect --------
  std::vector< TH2* > m_vHtriggerYstarSpect;
  std::vector< TH2* > m_vHtriggerYstarSpectFine;

  TH2* m_hAllYstarSpect;
  TH2* m_hAllYstarSpectFine;

  std::vector< TH3* > m_vHtriggerNjetsRun;
  TH3* m_hAllNjetsRun;
  std::vector< TH3* > m_vHtriggerNjetsRunFine;
  TH3* m_hAllNjetsRunFine;

  // ----- rtrk ----
  std::vector< TH3* > m_vHtriggerRtrk1;
  std::vector< TH3* > m_vHtriggerRtrk2;
  std::vector< TH3* > m_vHtriggerRtrk4;
  TH3* m_hAllTriggerRtrk1;
  TH3* m_hAllTriggerRtrk2;
  TH3* m_hAllTriggerRtrk4;
  
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
