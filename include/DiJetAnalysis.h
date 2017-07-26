#ifndef DIJETANALYSIS_H
#define DIJETANALYSIS_H

#include <string>
#include <vector>
#include <map>

#include "MyRoot.h"

class TH1;
class TH2;
class TH3;
class TF1;
class TGraphAsymmErrors;
class THnSparse;
class TLorentzVector;
class TEnv;
class TFile;
class TTree;

class DeltaPhiProj;

typedef double (*WeightFcn)( double, double, double );
typedef std::vector<std::vector<std::vector<std::vector<TH1*>>>> FourDTH1vector;

class DiJetAnalysis{
 public:
  DiJetAnalysis();
  DiJetAnalysis( bool );
  DiJetAnalysis( bool, bool );
  DiJetAnalysis( bool, bool, int );
  DiJetAnalysis( bool, bool, int, int );
  virtual ~DiJetAnalysis();

  //---------------------------
  // Initialization Methods
  //---------------------------
  virtual void Initialize();

  virtual void AdditionalSuffix( std::string& ) = 0;
  
  //---------------------------
  // Fill Tree / Plot Controls
  //---------------------------
  virtual void RunOverTreeFillHistos( int, int ) = 0;

  virtual void ProcessPlotHistos() = 0;

  virtual void DataMCCorrections() = 0;

  virtual void PlotHistosTogether() = 0;
  
  //---------------------------
  //       Fill Tree
  //---------------------------
  virtual void SetupHistograms() = 0;

  virtual void AddHistogram( TH1* );

  virtual void AddHistogram( THnSparse* );
  
  virtual void ProcessEvents( int, int ) = 0;  
  
  virtual void SaveOutputsFromTree();

  //---------------------------
  //       Analysis
  //---------------------------
  virtual bool GetFwdCentJets( const std::vector<TLorentzVector>&,
			       const TLorentzVector*&,
			       const TLorentzVector*& );


  virtual bool GetDiJets( const std::vector<TLorentzVector>&,
			  const TLorentzVector*&,
			  const TLorentzVector*& );
  
  virtual double AnalyzeDeltaPhi( THnSparse*,
				  const std::vector <TLorentzVector>&,
				  WeightFcn = NULL );

  virtual bool  ApplyIsolation( std::vector<TLorentzVector>&,
				double );
    
  virtual bool  ApplyCleaning ( std::vector<TLorentzVector>&, 
			        std::vector<bool>& );
  //---------------------------
  //       Tools
  //---------------------------  
  void FillHistoWithJets( const TLorentzVector*,
			  const TLorentzVector*,
			  TH2*, double = 1.);
  
  double AdjustEtaForPP( double );  

  double GetYstar( const TLorentzVector& );

  bool IsForwardDetector( const double& eta );

  bool IsCentralDetector( const double& eta );

  bool IsForwardYstar( const double& ystar );

  bool IsCentralYstar( const double& ystar );
  
  void NormalizeDeltaPhi( TH1* );

  virtual TH1*       CombineSamples( std::vector< TH1* >&,
				     const std::string& = "" );

  virtual TH2*       CombineSamples( std::vector< TH2* >&,
				     const std::string& = "" );
  
  virtual THnSparse* CombineSamples( std::vector< THnSparse* >&,
				     const std::string& = "" );  

  virtual void GetInfoBoth( std::string&, std::string&, std::string&, std::string&,
			    std::string&, std::string&, std::string& );

  virtual void GetInfoUnfolding( std::string&, std::string&, std::string& );

  TH1* BinByBinUnfolding( TH1*, TH1* );
  
  //---------------------------
  //   Get Quantities / Plot 
  //---------------------------
  virtual void LoadHistograms() = 0;

  virtual void MakeSpectra( std::vector< TH2* >&,
			    const std::vector< std::string >&, 
			    const std::string& = "" );
  
  virtual void MakeDeltaPhi( std::vector<THnSparse*>&,
			     const std::vector< std::string >&,
			     const std::string& = "",
			     bool = false );
  
  virtual THnSparse* UnfoldDeltaPhi( TFile*, TFile*,
				     const std::string& = "" );

  virtual void MakeDphiTogether();

  //---------------------------
  //        Drawing
  //---------------------------
  virtual void DrawAtlasRight( double = 0, double = 0, double = CT::StyleTools::lSS );

  virtual void DrawAtlasRightBoth( double = 0, double = 0, double = CT::StyleTools::lSS );

  virtual void DrawTopLeftLabels( DeltaPhiProj*, 
				  double = 0, double = 0,
				  double = 0, double = 0,
				  double = 0, double = 0,
				  double = 0, double = 0,
				  double = CT::StyleTools::lSS );
  
  //---------------------------
  //       Saving 
  //---------------------------
  virtual void SaveAsPdfPng( const TCanvas&,
			     const std::string& = "", const std::string& = "",
			     const std::string& = "", double = 0, double = 0,
			     const std::string& = "", double = 0, double = 0);
  
  virtual void SaveAsROOT( const TCanvas&,
			   const std::string& = "", const std::string& = "",
			   const std::string& = "", double = 0, double = 0,
			   const std::string& = "", double = 0, double = 0);

  virtual void SaveAsAll( const TCanvas&,
			  const std::string& = "", const std::string& = "",
			  const std::string& = "", double = 0, double = 0,
			  const std::string& = "", double = 0, double = 0);
  

 private:
  //========== config file  ===========
  TEnv* m_config;

 protected:
  //=========== common strings ========
  std::string m_s_pp;
  std::string m_s_pPb;
  std::string m_s_pt;
  
  std::string m_s_pt1;
  std::string m_s_pt2;
  
  std::string m_sOutput;
  std::string m_myOutName;
  std::string m_sMC;
  std::string m_sData;

  std::string m_recoName;
  std::string m_truthName;
  
  std::string m_sMUT;
  std::string m_sRatio;

  std::string m_allName;

  std::string m_unfoldingFileSuffix;

  std::vector< std::string > m_vMCtypeLabels;
  std::vector< std::string > m_vMCtypeNames;
  //============ cuts =============
  int    m_nMinEntriesFit;

  double m_dPhiThirdJetFraction;

  double m_dPhiZoomLow;
  
  double m_dPhiUnfoldingMin;
  double m_dPhiUnfoldingMax; 
  //===== settings and names ======
  TEnv* GetConfig(){ return m_config; }
  
  bool m_is_pPb;
  bool m_isData;
  int  m_mcType;
  int  m_uncertComp;
  
  std::string m_mcTypeName;
  std::string m_mcTypeLabel;
  
  std::string m_labelOut;
  std::string m_dirOut;
  
  std::string m_rawHistosFname;

  std::string m_fNameOutDefault;
  std::string m_fNameOut;
  std::string m_fNameOutUF;
  
  std::string m_fNameUnfoldingMC;
  
  DeltaPhiProj* m_dPP;
  
  // -------- eventCounter --------
  int m_ev;
  
 private:
  //============ data =============
  std::vector< TH1*       > v_hists;  // for writing
  std::vector< THnSparse* > v_hns;    // for writing
  std::vector< TF1*       > v_functs; // for writing
  std::vector< TGraphAsymmErrors* > v_graphs; // for writing
  
 protected:
  //========= common tools ========
  CT::AnalysisTools* anaTool;
  CT::DrawTools*     drawTool;
  CT::StyleTools*    styleTool;

 protected:
  //======== histos binning =======
  // -------- maps ---------
  int    m_nEtaMapBins;
  double m_etaMapMin;
  double m_etaMapMax;
  
  int    m_nPhiMapBins;
  double m_phiMapMin;
  double m_phiMapMax;

  int    m_nPtMapBins;
  double m_ptMapMin;
  double m_ptMapMax;

  // -------- spect --------
  int    m_nPtSpectBins;
  double m_ptSpectMin;
  double m_ptSpectMax;
  
  // ---- JES/PRes/Etc ----- 
  int    m_nEtaForwardBinsFine;
  double m_etaForwardMin;
  double m_etaForwardMax;

  // ---- forward eta binning ---
  std::vector<double> m_varFwdEtaBinning;
  unsigned int m_nVarFwdEtaBins;

  // --- whole range variable ----
  // ------- eta/ystar binning --------
  std::vector<double> m_varYstarBinningA;
  unsigned int m_nVarYstarBinsA;

  std::vector<double> m_varYstarBinningB;
  unsigned int m_nVarYstarBinsB;
  
  // ---- variable pt binning ----
  std::vector<double> m_varPtBinning;
  unsigned int m_nVarPtBins;

  // --- variable dphi binnign ---
  std::vector<double> m_varDphiBinning;
  unsigned int m_nVarDphiBins;

  int m_nDphiBinsLarge ; 
  int m_nDphiBinsMedium; 
  int m_nDphiBinsSmall ;
  
  int m_dPhiBinsLargeFactor ;
  int m_dPhiBinsMediumFactor;
  int m_dPhiBinsSmallFactor ;

  // -------- eff ---------
  double m_effMin;
  double m_effMax;

  // -------- dphi --------
  uint   m_nDphiDim;

  double m_dPhiDphiMin;
  double m_dPhiDphiMax;
  
  std::vector< int >    m_nDphiBins;
  std::vector< double > m_dPhiMin;
  std::vector< double > m_dPhiMax;
    
  double m_dPhiWidthMin;
  double m_dPhiWidthMax;

  // -------- ratios ------
  double m_ratioMax;
  double m_ratioMin;
  
  //===== common histo names =======
  std::string m_etaSpectName;
  std::string m_dPhiName;
  std::string m_effName;
  std::string m_purityName;
  
  std::string m_respMatName;
  std::string m_unfoldedName;
  
  std::string m_dPhiRecoName;
  std::string m_dPhiTruthName;
  
  std::string m_dPhiRespMatName;
  std::string m_ptRespMatName;
  
  std::string m_allRespMatName;  

  std::string m_dPhiCFactorsName;
  std::string m_dPhiUnfoldedName;
  std::string m_dPhiRecoUnfoldedName;
};

#endif
