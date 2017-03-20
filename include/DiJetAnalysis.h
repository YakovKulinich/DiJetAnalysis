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

  virtual void ProcessEvents( int, int ) = 0;  

  virtual bool ApplyIsolation( double, std::vector<TLorentzVector>& );
    
  virtual void ApplyCleaning( std::vector<TLorentzVector>&, 
			      std::vector<bool>& );
      
  //---------------------------
  //       Plotting 
  //---------------------------
  virtual void ProcessPlotHistos() = 0;

  virtual void LoadHistograms() = 0;

  virtual void SaveOutputs();

  //---------------------------
  //       Tools
  //---------------------------
  void AddHistogram( TH1* );

  void FitGaussian( TH1*, TF1* );
  
  std::string GetEtaLabel( double, double );

  double AdjustEtaForPP( double );
  
 protected:
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
  //============ histos =============
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
 public:
 

};

#endif
