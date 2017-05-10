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
class THnSparse;
class TLorentzVector;
class TEnv;
class TFile;
class TTree;

typedef double (*WeightFcn)( double, double, double );
typedef std::vector<std::vector<std::vector<std::vector<TH1*>>>> fourDTH1vector;

class DiJetAnalysis{
 public:
  DiJetAnalysis();
  DiJetAnalysis( bool, bool, int = 0 );
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
  virtual bool GetDiJets( const std::vector
			  <TLorentzVector>&,
			  TLorentzVector&,
			  TLorentzVector& );
  
  virtual double AnalyzeDeltaPhi( THnSparse*, THnSparse*,
				  const std::vector <TLorentzVector>&,
				  double = 1,
				  WeightFcn = NULL );
  
  virtual void  ApplyIsolation( double, std::vector<TLorentzVector>& );
    
  virtual void  ApplyCleaning ( std::vector<TLorentzVector>&, 
			        std::vector<bool>& );
  //---------------------------
  //       Tools
  //---------------------------  
  double AdjustEtaForPP( double );  

  double GetYstar( TLorentzVector& );

  bool IsForwardDetector( const double& eta );

  bool IsCentralDetector( const double& eta );

  bool IsForwardYstar( const double& ystar );

  bool IsCentralYstar( const double& ystar );

  virtual void CombineJZN( TH1*,
			   std::vector< TH1*>& ){};
  
  virtual void CombineJZN( TH1*,
			   std::vector< TH1*>&,
			   std::vector< TH1*>& ){};

  
  //---------------------------
  //       Plotting 
  //---------------------------
  virtual void LoadHistograms() = 0;
 
  virtual void PlotDphiTogether();
  
  virtual void PlotDeltaPhi( std::vector<THnSparse*>&,
			     std::vector<THnSparse*>&,
			     fourDTH1vector&,
			     fourDTH1vector&,
			     const std::vector< std::string >&,
			     const std::string& = "" ,
			     const std::string& = "" );

  virtual void PlotDeltaPhi( std::vector<THnSparse*>&,
			     std::vector<THnSparse*>&,
			     const std::vector< std::string >&,
			     const std::string& = "" ,
			     const std::string& = "" );

  
  //---------------------------
  //        Drawing
  //---------------------------
  void DrawCanvas( std::vector< TH1* >&,
		   const std::string& = "",
		   const std::string& = "",
		   int = 1);

  void DrawCanvas( std::vector< TH1* >&,
		   const std::string& = "",
		   const std::string& = "",
		   bool = true );

  void DrawCanvas( std::vector< TGraphAsymmErrors* >&,
		   const std::string& = "",
		   const std::string& = "",
		   double = 0, double = 0);
  
  //===== MinMax and line drawing =====
  void SetMinMax( TH1*,
		  const std::string&,
		  const std::string& );

  double GetLineHeight( const std::string& );

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
  
  bool m_isData;
  bool m_is_pPb;

  int m_mcType;
  std::string m_mcTypeLabel;
  
  std::string m_labelOut;
  std::string m_dirOut;
  std::string m_rootFname;
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

  // -------- dphi- --------
  uint   m_nDphiDim;
  uint   m_nDphiNentDim;
  
  int    m_nDphiDphiBins;
  double m_nDphiDphiMin;
  double m_nDphiDphiMax;

  std::vector< int >    m_nDphiBins;
  std::vector< double > m_dPhiMin;
  std::vector< double > m_dPhiMax;
    
  double m_dPhiWidthMin;
  double m_dPhiWidthMax;
};

#endif
