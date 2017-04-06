#ifndef DIJETANALYSIS_H
#define DIJETANALYSIS_H

#include <TLorentzVector.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TF1.h>
#include <TGraphAsymmErrors.h>
#include <TFile.h>
#include <TTree.h>
#include <THnSparse.h>

#include <string>
#include <vector>
#include <map>

class DiJetAnalysis{
 public:
  DiJetAnalysis();
  DiJetAnalysis( bool, bool );
  virtual ~DiJetAnalysis();

  virtual void Initialize();
  
  //---------------------------
  //       Fill Tree
  //---------------------------
  virtual void RunOverTreeFillHistos( int, int ) = 0;

  virtual void SetupHistograms() = 0;

  virtual void AddHistogram( TH1* );

  virtual void AddHistogram( THnSparse* );
  
  virtual void ProcessEvents( int, int ) = 0;  
  
  virtual void SaveOutputsFromTree();

  //---------------------------
  //       Analysis
  //---------------------------
  virtual double AnalyzeDeltaPhi( THnSparse*,
				  const std::vector
				  <TLorentzVector>& );
  
  virtual void  ApplyIsolation( double, std::vector<TLorentzVector>& );
    
  virtual void  ApplyCleaning ( std::vector<TLorentzVector>&, 
			        std::vector<bool>& );
  //---------------------------
  //       Tools
  //---------------------------
  void ProjectAndFitGaus( TH3*, TH1*, TH1*,
			  int, int,
			  const std::string& = "" ,
			  const std::string& = "" );
  
  void FitGaussian( TH1*, TF1* );
  
  std::string GetEtaLabel( double, double );
  
  std::string GetName( double, double, const std::string& );
  
  double AdjustEtaForPP( double );

  //---------------------------
  //       Plotting 
  //---------------------------
  virtual void ProcessPlotHistos() = 0;

  virtual void LoadHistograms() = 0;
  
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
  
  
 protected:
  //============ cuts =============
  int    m_nMinEntriesGausFit;
  double m_ptFitMin;

  double m_dPhiThirdJetFraction;
  //========== settings ===========
  bool m_isData;
  bool m_is_pPb;
  
  std::string m_labelOut;
  std::string m_dirOut;
  std::string m_rootFname;
  
  //============ files =============
  TFile* m_fIn;
  TFile* m_fOut;
  TTree* m_tree;

  // -------- eventCounter --------
  int m_ev;
  
 private:
  //============ data =============
  std::vector< TH1*       > v_hists;  // for writing
  std::vector< THnSparse* > v_hns;    // for writing
  std::vector< TF1*       > v_functs; // for writing
  std::vector< TGraphAsymmErrors* > v_graphs; // for writing

 protected:
  // -------- dPhi --------
  std::map< std::string, THnSparse* > m_mDphi;

  
  
 protected:
  //========= histos binning ========
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

  std::vector<double> m_varEtaBinning;
  unsigned int m_nVarEtaBins;
  
  // -------- eff ---------
  double m_effMin;
  double m_effMax;

  // -------- dphi- --------
  int    m_nDphiPtBins;
  double m_nDphiPtMin;
  double m_nDphiPtMax;
  
  int    m_nDphiDphiBins;
  double m_nDphiDphiMin;
  double m_nDphiDphiMax;

};

#endif
