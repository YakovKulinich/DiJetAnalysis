#ifndef DIJETANALYSIS_H
#define DIJETANALYSIS_H

#include <string>
#include <vector>
#include <map>

#include "MyRoot.h"

namespace CS{
  class AnalysisTools;
  class DrawTools;
  class StyleTools;
}

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
  virtual ~DiJetAnalysis();

  virtual void Initialize();

  //---------------------------
  // Fill Tree / Plot Controls
  //---------------------------
  virtual void RunOverTreeFillHistos( int, int ) = 0;

  virtual void ProcessPlotHistos() = 0;
  
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
  
  virtual double AnalyzeDeltaPhi( THnSparse*, THnSparse*,
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

  // Make These Templates Later...
  virtual TH1*       CombineSamples( std::vector< TH1* >&,
				     const std::string& = "" );

  virtual TH2*       CombineSamples( std::vector< TH2* >&,
				     const std::string& = "" );
  
  virtual THnSparse* CombineSamples( std::vector< THnSparse* >&,
				     const std::string& = "" );  

  virtual void GetInfoBoth( std::string&, std::string&, std::string&, std::string&,
			    std::string&, std::string&, std::string& );
  
  //---------------------------
  //   Get Quantities / Plot 
  //---------------------------
  virtual void LoadHistograms() = 0;

  virtual void MakeSpectra( std::vector< TH2* >&,
			    const std::vector< std::string >&, 
			    const std::string& = "" );
  
  virtual void MakeDeltaPhi( std::vector<THnSparse*>&,
			     const std::vector< std::string >&,
			     const std::string& = "" );

  virtual void MakeDphiTogether();
  
  //---------------------------
  //        Drawing
  //---------------------------
  virtual void DrawAtlasRight( double = 0, double = 0, double = CT::StyleTools::lSS );

  virtual void DrawAtlasRightBoth( double = 0, double = 0, double = CT::StyleTools::lSS );
  
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
  //========== settings ===========
  TEnv* m_config;

 protected:
  //============ cuts =============
  int    m_nMinEntriesFit;

  double m_dPhiThirdJetFraction;
  //========== settings ===========
  TEnv* GetConfig(){ return m_config; }
  
  bool m_is_pPb;
  bool m_isData;
  int  m_mcType;

  std::string m_allName;
    
  std::string m_labelOut;
  std::string m_dirOut;
  std::string m_rootFname;
  
  DeltaPhiProj* m_dPP;

  // -------- eventCounter --------
  int m_ev;
  
  //============ files =============
  TFile* m_fIn;
  TFile* m_fOut;
  TTree* m_tree;
  
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

  // --- whole range variable ---
  // ------- eta/ystar binning --------
  std::vector<double> m_varYstarBinningA;
  unsigned int m_nVarYstarBinsA;

  std::vector<double> m_varYstarBinningB;
  unsigned int m_nVarYstarBinsB;
  
  // ---- variable pt binning ---
  std::vector<double> m_varPtBinning;
  unsigned int m_nVarPtBins;
  
  // -------- eff ---------
  double m_effMin;
  double m_effMax;

  // -------- dphi --------
  uint   m_nDphiDim;
  
  int    m_nDphiDphiBins;
  double m_dPhiDphiMin;
  double m_dPhiDphiMax;

  std::vector< int >    m_nDphiBins;
  std::vector< double > m_dPhiMin;
  std::vector< double > m_dPhiMax;
    
  double m_dPhiWidthMin;
  double m_dPhiWidthMax;

  //===== common histo names =======
  std::string m_etaSpectName;
  std::string m_dPhiName;
  std::string m_effName;
  
  std::string m_recoName;
  std::string m_truthName;
  std::string m_respMatName;
  std::string m_unfoldedName;
  
  std::string m_dPhiRecoName;
  std::string m_dPhiTruthName;
  std::string m_dPhiRespMatName;
  std::string m_dPhiUnfoldedName;
};

#endif
