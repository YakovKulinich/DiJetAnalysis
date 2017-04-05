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
#include <THmulf.h>

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
  
  virtual void ProcessEvents( int, int ) = 0;  
  
  virtual void SaveOutputsFromTree();

  //---------------------------
  //       Analysis
  //---------------------------
  virtual void ApplyIsolation( double, std::vector<TLorentzVector>& );
    
  virtual void ApplyCleaning( std::vector<TLorentzVector>&, 
			      std::vector<bool>& );

  //---------------------------
  //       Plotting 
  //---------------------------
  virtual void ProcessPlotHistos() = 0;

  virtual void LoadHistograms() = 0;

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
  
  //---------------------------
  //       Tools
  //---------------------------
  void ProjectAndFitGaus( TH3*, TH1*, TH1*,
			  int, int,
			  const std::string& = "" ,
			  const std::string& = "" );
  
  void FitGaussian( TH1*, TF1* );
  
  std::string GetEtaLabel( double, double );

  double AdjustEtaForPP( double );
  
 protected:
  //============ cuts =============
  int    m_nMinEntriesGausFit;
  double m_ptFitMin;
  
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

  //============ data =============
 private:
  std::vector< TH1*    > v_hists;  // for writing
  std::vector< TF1*    > v_functs; // for writing
  std::vector< TGraphAsymmErrors* > v_graphs; // for writing

 protected:
  //========= histos binning ========
  int    m_nPtSpectBins;
  double m_ptSpectMin;
  double m_ptSpectMax;
  
  // Eta-Phi Maps
  int    m_nEtaBins;
  double m_etaMin;
  double m_etaMax;
  
  int    m_nPhiBins;
  double m_phiMin;
  double m_phiMax;
  
  double m_ptWidth;
  double m_ptMin;
  double m_ptMax;
  int    m_nPtBins;

  // for jes, jer, efficiency, etc.
  int    m_nEtaForwardBinsFine;
  double m_etaForwardMin;
  double m_etaForwardMax;

  std::vector<double> m_varEtaBinning;
  int m_nVarEtaBins;

  // efficiency
  double m_effMin;
  double m_effMax;
 public:
};

#endif
