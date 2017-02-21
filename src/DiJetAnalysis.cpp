#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TROOT.h>
#include <TEnv.h>
#include <TCanvas.h>

#include <iostream>

#include "MyRoot.h"

#include "DiJetAnalysis.h"

DiJetAnalysis::DiJetAnalysis()
  : m_isData(true), m_is_pPb( true )
{}

DiJetAnalysis::DiJetAnalysis( bool isData, bool is_pPb )
  : m_isData( isData ), m_is_pPb( is_pPb )
{}

DiJetAnalysis::~DiJetAnalysis(){}

void DiJetAnalysis::Initialize(){
  m_labelOut = m_isData ? "_data" : "_mc" ;
  m_labelOut = m_is_pPb ? m_labelOut + "_pPb" : m_labelOut + "_pp";
  
  m_fNameOut = "output/myOut" + m_labelOut + ".root";
  
  std::cout << "fNameOut: " << m_fNameOut << std::endl;
}

bool DiJetAnalysis::applyIsolation( double Rmin, std::vector<TLorentzVector>& v_jets ){
  
  std::vector<bool> isIsolated;

  for(unsigned int iTestJet = 0; iTestJet < v_jets.size(); iTestJet++){
    for(unsigned int iSecondJet = 0; iSecondJet < v_jets.size(); iSecondJet++){   
      if( iSecondJet == iTestJet ) continue;

      if( DeltaR( v_jets.at(iTestJet), v_jets.at(iSecondJet)) < Rmin &&
	  v_jets.at(iSecondJet).Pt() > v_jets.at(iTestJet).Pt() * 0.5 ){       
	isIsolated.push_back(false);
	continue;
      }
    } // end loop over iSecondJet
    isIsolated.push_back(true);
  } // end loop over iTestJet

  for(unsigned int iJet = 0; iJet < v_jets.size(); iJet++ ){
    if( !isIsolated.at(iJet) ) v_jets.at(iJet).SetPxPyPzE(0,0,0,-1 );
  }

  return true;
}

void DiJetAnalysis::applyCleaning( std::vector<TLorentzVector>& v_jets, 
		    std::vector<bool>& v_isCleanJet){
  for( unsigned int jn = 0; jn < v_jets.size() ; jn++ ){
    if( !v_isCleanJet.at(jn) ) v_jets.at(jn).SetPxPyPzE(0,0,0,-1 );
  }
}
