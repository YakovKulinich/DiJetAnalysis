#ifndef DIJETANALYSISBOTH_H
#define DIJETANALYSISBOTH_H

#include "DiJetAnalysis.h"

class DiJetAnalysisBoth : public DiJetAnalysis{

 public:
  DiJetAnalysisBoth();
  DiJetAnalysisBoth( bool, bool );
  ~DiJetAnalysisBoth();

  void Initialize();

  void AdditionalSuffix( std::string& ){}

  //---------------------------
  // Fill Tree / Plot Controls
  //---------------------------
  void RunOverTreeFillHistos( int, int ){}
  void ProcessPerformance   (){};
  void UnfoldPerformance    (){};
  void ProcessPhysics       (){};
  void UnfoldPhysics        (){};
  void MakeResultsTogether  ();

  //---------------------------
  //       Fill Tree
  //---------------------------
  void SetupHistograms(){}

  void ProcessEvents( int, int ){}

  //---------------------------
  //   Get Quantities / Plot 
  //---------------------------
  void LoadHistograms( int = 0 ){}
  
  void MakeDphiTogether();

 private:
  bool m_isReco;

  std::string m_system;
  std::string m_mcLevel;

  std::string m_sBoth;
  
  std::vector< std::string > m_vMC;
  std::vector< std::string > m_vMClabels;  
};

#endif
