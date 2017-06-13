#ifndef DIJETANALYSISBOTH_H
#define DIJETANALYSISBOTH_H

#include "DiJetAnalysis.h"

class DiJetAnalysisBoth : public DiJetAnalysis{

 public:
  DiJetAnalysisBoth();
  DiJetAnalysisBoth( bool, bool );
  ~DiJetAnalysisBoth();

  // pure virtual function from DiJetAnalysis
  void RunOverTreeFillHistos( int, int ){}
  void ProcessPlotHistos(){}

  void SetupHistograms(){}
  void ProcessEvents( int, int ){}

  void LoadHistograms(){}
  
  void Initialize();
  
  void PlotDphiTogether();

 private:
  bool m_isReco;

  std::string m_system;
  std::string m_mcLevel;
  
  std::vector< std::string > m_vMC;
  std::vector< std::string > m_vMClabels;  
};

#endif
