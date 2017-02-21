#ifndef DIJETANALYSIS_H
#define DIJETANALYSIS_H

#include <TLorentzVector.h>
#include <TLine.h>
#include <TLegend.h>
#include <TH1.h>

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

  virtual void setupHistograms() = 0;

  virtual void processEvents( int, int ) = 0;  

  bool applyIsolation( double, std::vector<TLorentzVector>& );
    
  void applyCleaning( std::vector<TLorentzVector>&, 
		      std::vector<bool>& );
      
  //---------------------------
  //       Plotting 
  //---------------------------
  virtual void PlotExistingHistos() = 0;

  virtual void loadHistograms() = 0;

  virtual void plotSpectra() = 0;

  virtual void plotEtaPhi() = 0;

  virtual void plotPtEta() = 0;
  
 protected:
  //========== settings ===========
  bool m_isData;
  bool m_is_pPb;

  std::string m_labelOut;
  std::string m_fNameOut;
  
  //============ data =============
  std::vector< TH1* > v_hists; // for writing

  TCanvas* c_spect = NULL;
  TLegend* l_spect = NULL;

  TCanvas* c_etaPhi = NULL;
  TCanvas* c_ptEta  = NULL;
};

#endif
