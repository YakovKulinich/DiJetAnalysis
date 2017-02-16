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
  ~DiJetAnalysis();
  
  //---------------------------
  //       Fill Tree
  //---------------------------
  void RunOverTreeFillHistos( int, int, int, std::string&, std::string& );

  void loadTriggers();

  void setupHistograms();

  void processEventsInData( int, int, std::string&, std::string& );  

  void processEfficiencies( std::vector< TLorentzVector >&,
			    std::vector< TLorentzVector >&,
			    std::map< std::string, bool >&);

  bool applyIsolation( double, std::vector<TLorentzVector>& );
    
  void applyCleaning( std::vector<TLorentzVector>&, 
		      std::vector<bool>& );
      
  //---------------------------
  //       Plotting 
  //---------------------------
  void PlotExistingHistos( int, std::string&, std::string& );

  void loadFinalTriggers();

  void loadHistograms( std::string& );

  void plotSpectra( std::string& );

  void plotEfficiencies( std::string& );

  void plotEtaPhi( std::string& );

  void plotPtEta( std::string& );

 private:
  std::vector< std::string > v_triggers;

  std::map< std::string, TH1* > m_triggerSpect;
  std::map< std::string, TH1* > m_triggerEff;
  std::map< std::string, TH1* > m_triggerRunPrescale;

  std::map< std::string, TH1* > m_triggerEtaPhi;
  std::map< std::string, TH1* > m_triggerPtEta;

  std::vector< TH1* > v_hists; // for writing

  std::vector< int > v_tJetPt;
  std::multimap< int , std::string > m_tJetPtTrigger;

  std::string  mbTrigger;

  TH1* h_EtaPhi;

  TCanvas* c_spect = NULL;
  TLegend* l_spect = NULL;

  TCanvas* c_eff = NULL;
  TLegend* l_eff = NULL;
  TLine*   line  = NULL;

  TCanvas* c_etaPhi = NULL;
  TCanvas* c_ptEta  = NULL;
};

#endif
