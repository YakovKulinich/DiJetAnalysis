#ifndef DIJETANALYSIS_H
#define DIJETANALYSIS_H

#include <TLorentzVector.h>
#include <TH1.h>
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

  void Initialize();
  
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
  
 protected:
  //========== settings ===========
  bool m_isData;
  bool m_is_pPb;

  std::string m_labelOut;
  std::string m_dirOut;
  std::string m_fNameOut;
  std::string m_fNameIn;

  //============ files =============
  TFile* m_fIn;
  TFile* m_fOut;
  TTree* m_tree;
  
  //============ data =============
  std::vector< TH1*    > v_hists;  // for writing
  std::vector< TF1*    > v_functs; // for writing
  std::vector< TGraphAsymmErrors* > v_graphs; // for writing
};

#endif
