#ifndef DIJETANALYSISMC_H
#define DIJETANALYSISMC_H

#include "THmulf.h"

#include "DiJetAnalysis.h"

class DiJetAnalysisMC : public DiJetAnalysis{
 public:
  DiJetAnalysisMC();
  DiJetAnalysisMC( bool, bool );
  ~DiJetAnalysisMC();
  
  //---------------------------
  //       Fill Tree
  //---------------------------
  void RunOverTreeFillHistos( int, int );

  void setupHistograms();

  void processEvents( int, int );  

  //---------------------------
  //       Plotting 
  //---------------------------
  void PlotExistingHistos();

  void loadHistograms();

  void plotSpectra();

  void plotEtaPhi();

  void plotEtaPt();
  
 private:
  const int nJZNmax;
  
  //============ data =============
  std::vector< std::string > m_jznFname;
  std::vector< TH1* > m_etaPhi;
  std::vector< TH1* > m_ptEta;
};

#endif
