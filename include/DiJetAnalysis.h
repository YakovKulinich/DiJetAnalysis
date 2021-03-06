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

// this variable is to change labels and or axis ranges
// on the final plots that go into the paper.
// this is for ther performance plots, dPhi w/syst,
// overall systematics, and final plots.
static const bool finalPlots = true;

static const bool doRtrk = false;

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
  virtual void RunOverTreeFillSpect( int, int ){}

  virtual void ProcessSpectWeights(){}

  virtual void RunOverTreeFillDphi( int, int ){}
    
  virtual void ProcessDphiWeights(){}

  virtual void RunOverTreeFillHistos( int, int ) = 0;

  virtual void ProcessPerformance   () = 0;

  virtual void UnfoldPerformance    () = 0;

  virtual void ProcessPhysics       () = 0;

  virtual void UnfoldPhysics        () = 0;

  virtual void MakeResultsTogether  ();
  
  virtual void ProcessSystematics   ();
  
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
			  const TLorentzVector*&,
			  bool = true );

  virtual bool PassHECCuts( const TLorentzVector& );
  
  virtual void AnalyzeSpectra( TH2*, const std::vector< TLorentzVector >& );
  
  virtual double AnalyzeDeltaPhi( THnSparse*,
				  const std::vector <TLorentzVector >& );
  
  virtual bool  ApplyIsolation( std::vector< TLorentzVector >&,
				double );
    
  virtual bool  ApplyCleaning ( std::vector< TLorentzVector >&, 
			        std::vector< bool >& );
  //---------------------------
  //       Tools
  //---------------------------  
  void FillHistoWithJets( const TLorentzVector*,
			  const TLorentzVector*,
			  TH2*, double = 1.);
  
  double GetYstar( const TLorentzVector& );

  bool IsForwardDetector( const double& eta );

  bool IsCentralDetector( const double& eta );

  bool IsForwardYstar( const double& ystar );

  bool IsCentralYstar( const double& ystar );
  
  TF1* NormalizeDeltaPhi( TH1*, TH1* = NULL, double = 0, bool = true );
 
  virtual double GetJetWeight( const TLorentzVector& );

  virtual TH1*       CombineSamples( std::vector< TH1* >&,
				     const std::string& = "" );

  virtual TH2*       CombineSamples( std::vector< TH2* >&,
				     const std::string& = "" );
  
  virtual THnSparse* CombineSamples( std::vector< THnSparse* >&,
				     const std::string& = "" );  

  virtual void GetSpectraLabels( std::string&, std::string&,
				 const std::string& = "" );

  virtual void GetSpectUnfoldingInfo( std::string&, std::string&, std::string&,
				      std::string&, std::string&, std::string& ){}
  
  virtual void GetDphiUnfoldingInfo( std::string&, std::string&,
				     std::string&, std::string& ){}

  virtual void GetInfoTogether( std::string&, std::string&, std::string&,
				std::string&, std::string&, std::string&,
				int = 0 ){}

  TH1*   BinByBinUnfolding( TH1*, TH1* );

  TFile* GetListOfSystUncert( std::vector< int >&, std::map< int, TFile* >& );

  void MakeDefaultBinning( std::vector< double >&, std::vector< double >&,
			   int = 4, int = 8, int = 4);

  void MakeLinearBinning( std::vector< double >&, std::vector< double >&, int );

  void MakeLogBinning( std::vector< double >&, std::vector< double >&, int );

  TH1* FlipOverXaxis( TH1*, std::vector< double >& );

  void GetSignficance(TGraphAsymmErrors*, TGraphAsymmErrors*);
  
  //---------------------------
  //   Get Quantities / Plot 
  //---------------------------
  virtual void LoadHistograms( int = 0 ) = 0;

  virtual void MakeEtaPhiPtMap( std::vector< TH2* >&,
				const std::vector< std::string >&, 
				const std::string& = "" );
  
  virtual void MakeSpectra( std::vector< TH2* >&,
			    const std::vector< std::string >&, 
			    const std::string& = "" );

  virtual void MakeRtrk( std::vector< TH3* >&,
			 const std::vector< std::string >&, 
			 const std::string& = "" );
  
  virtual TH2* UnfoldSpectra( TFile*, TFile*,
			      const std::string& = "" );
  
  virtual void MakeDeltaPhi( std::vector<THnSparse*>&,
			     const std::vector< std::string >&,
			     const std::string& = "",
			     TFile* = NULL, 
			     const std::string& = "" );
  
  virtual THnSparse* UnfoldDeltaPhi( TFile*, TFile*,
				     const std::string& = "",
				     TFile* = NULL,
				     const std::string& = "" );

  virtual void MakeSpectTogether( TFile* = NULL );
  
  virtual void MakeDphiTogether( TFile* = NULL );

  virtual void CompareWeightIsoPtCuts( TFile* = NULL );

  virtual void MakeFinalPlotsTogether( TFile* fOut, const std::string& = "" );

  //---------------------------
  //        Drawing
  //---------------------------
  virtual void DrawAtlasRight( double = CT::DrawTools::drawX0,
			       double = CT::DrawTools::drawY0,
			       double = CT::StyleTools::lSS );

  virtual void DrawAtlasRightBoth(  double = CT::DrawTools::drawX0,
				    double = CT::DrawTools::drawY0,
				    double = CT::StyleTools::lSS );
  
  virtual void DrawTopLeftLabels( DeltaPhiProj*, 
				  double = 0, double = 0,
				  double = 0, double = 0,
				  double = 0, double = 0,
				  double = 0, double = 0,
				  double = CT::StyleTools::lSS );
  
  //---------------------------
  //       Saving 
  //---------------------------
  virtual void SaveAsROOT( const TCanvas&,
			   const std::string& = "" );


  virtual void SaveAsPdfPng( const TCanvas&,
			     const std::string& = "",
			     bool = false );
  
  virtual void SaveAsAll( const TCanvas&,
			  const std::string& = "",
			  bool = false );
  

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

  std::string m_sEta;
  std::string m_sYstar;
  
  std::string m_sOutput;
  std::string m_myOutName;
  std::string m_sMC;
  std::string m_sData;
  std::string m_sFinal;

  std::string m_sRaw;
  std::string m_sPhys;
  std::string m_sPerf;

  std::string m_sRuns;
  std::string m_sFine;
  std::string m_sWeights;
  std::string m_sCounts;
  std::string m_sReb;
  std::string m_sRatio;
  std::string m_sSum;

  std::string m_sDphi;
  std::string m_sDphiTitle;
  std::string m_sWidthTitle;
  std::string m_sYieldTitle;
  std::string m_sWidthRatioTitle;
  std::string m_sYieldRatioTitle;
  
  std::string m_rawFileSuffix;
  std::string m_performanceFileSuffix;
  std::string m_physicsFileSuffix;

  std::string m_unweightedFileSuffix;
  std::string m_unfoldingFileSuffix;
  std::string m_systematicsFileSuffix;

  std::vector< std::string > m_vMCtypeLabels;
  std::vector< std::string > m_vMCtypeNames;
  //============ cuts =============
  int    m_nMinEntriesFit;

  double m_deltaPtCut;
  
  double m_jetDeltaR;
  
  double m_hecEtaMinA;
  double m_hecEtaMaxA;

  double m_hecPhiMinA;
  double m_hecPhiMaxA;

  double m_hecEtaMinB;
  double m_hecEtaMaxB;

  double m_hecPhiMinB;
  double m_hecPhiMaxB;
  
  double m_dPhiZoomLow;
  double m_dPhiZoomHigh;
  
  double m_dPhiZoomLowBin;
  double m_dPhiZoomHighBin;

  double m_dPhiRebinnedZoomLowBin;
  double m_dPhiRebinnedZoomHighBin;

  double m_dPhiLogMin;
  double m_dPhiLogMax;

  double m_dPhiFittingMin;
  double m_dPhiFittingMax;

  double m_dPhiFittingMinB;
  double m_dPhiFittingMinC;
  double m_dPhiFittingMinD;
  
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
  
  std::string m_uncertSuffix;
  
  std::string m_fName;

  std::string m_fNameRaw;
  std::string m_fNameRawUW;
  std::string m_fNamePerf;
  std::string m_fNamePerfUF;
  std::string m_fNamePhys;
  std::string m_fNamePhysUF;
  
  std::string m_fNameDefRaw;
  std::string m_fNameDefRawUW;
  std::string m_fNameDefPerf;
  std::string m_fNameDefPerfUF;
  std::string m_fNameDefPhys;
  std::string m_fNameDefPhysUF;

  std::string m_fNameSYS;

  std::string m_dirOutTogether; 
  std::string m_fNameTogether;

  //===== common histo names =======  
  std::string m_allName;

  std::string m_spectName;
  std::string m_dPhiName;
  std::string m_widthName;
  std::string m_yieldName;

  std::string m_recoName;
  std::string m_truthName;
  std::string m_pairedName;
  std::string m_unpairedName;

  std::string m_cFactorName;
  std::string m_respMatName;
  std::string m_unfoldedName;

  // -------- effic --------
  std::string m_etaEffName;
  std::string m_ystarEffName;
  
  // -------- spect --------
  std::string m_etaSpectName;
  
  std::string m_ystarSpectName;
  std::string m_ystarSpectRecoName;
  std::string m_ystarSpectTruthName;

  std::string m_ystarSpectFineName;
  std::string m_ystarSpectFineRecoName;
  std::string m_ystarSpectFineTruthName;

  // --- rtrk ----
  std::string m_rtrkName;
  
  // --- ystar response matrix ----
  std::string m_ystarRespMatName;

  // --- spectra response matrix ----
  std::string m_ystarSpectRespMatName;
  
  // ------- spectra cfactors -------
  std::string m_ystarSpectCfactorsName;

  // ------- unfolded spectra -------
  std::string m_ystarSpectUnfoldedName;
  std::string m_ystarSpectRecoUnfoldedName;

  // ------------- dPhi ---------------
  std::string m_dPhiRecoName;
  std::string m_dPhiTruthName;

  // ------- dPhi cfactors ------------
  std::string m_dPhiCfactorsName;

  // ------- dPhi unfolded ------------
  std::string m_dPhiUnfoldedName;

  // -------- systematics plots -------
  std::string m_systematicsName;
  std::string m_dPhiSystematicsName;
    
  std::string m_effName;
  std::string m_purityName;

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

  //=== Set DeltaPhi Axes Order ===
  DeltaPhiProj* m_dPP;

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

  double m_ptSpectYaxisMin;
  
  // -------- eff ---------
  double m_effMin;
  double m_effMax;
  
  // ---- forward eta binning ---
  std::vector< double > m_varFwdEtaBinning;
  uint m_nVarFwdEtaBins;

  // --- whole range variable ----
  std::vector<double> m_varYstarBinning;
  std::vector<double> m_varYstarBinningFlipped;
  uint m_nVarYstarBins;

  std::vector<double> m_varEtaBinning;
  uint m_nVarEtaBins;

  // ---- barrel eta binning ----
  std::vector<double> m_varEtaBarrelBinning;
  uint m_nVarEtaBarrelBins;
  
  // ---- variable pt binning ----
  std::vector<double> m_varPtBinning;
  uint m_nVarPtBins;

  std::vector<double> m_varRtrkPtBinning;
  uint m_nVarRtrkPtBins;
  
  // --- variable pt binning w/over+underflow ---
  std::vector< double > m_varPtBinningUfOf;
  uint m_nVarPtBinsUfOf;
  
  // --- variable dphi binning ---
  std::vector<double> m_varDphiBinning;
  uint m_nVarDphiBins;

  // --- rebinned variable dphi binning ---
  std::vector<double> m_varDphiRebinnedBinning;
  uint m_nVarDphiRebinnedBins;

  // -------- dphi --------
  uint   m_nDphiDim;

  double m_dPhiDphiMin;
  double m_dPhiDphiMax;
  
  std::vector< int >    m_vNdPhiBins;
  std::vector< double > m_vDphiMin;
  std::vector< double > m_vDphiMax;
    
  double m_dPhiWidthMin;
  double m_dPhiWidthMax;

  double m_dPhiYieldMin;
  double m_dPhiYieldMax;
  
  // -------- ratios ------
  double m_ratioMax;
  double m_ratioMin;
};

#endif
