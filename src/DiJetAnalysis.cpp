#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TROOT.h>
#include <TEnv.h>
#include <TCanvas.h>

#include <iostream>
#include <boost/filesystem.hpp>

#include "MyRoot.h"

#include "DiJetAnalysis.h"

DiJetAnalysis::DiJetAnalysis()
  : m_isData(true), m_is_pPb( true ),
    m_fIn(NULL), m_fOut(NULL), m_tree(NULL)
{}

DiJetAnalysis::DiJetAnalysis( bool isData, bool is_pPb )
  : m_isData( isData ), m_is_pPb( is_pPb ),
    m_fIn(NULL), m_fOut(NULL), m_tree(NULL)
{}

DiJetAnalysis::~DiJetAnalysis(){}

void DiJetAnalysis::Initialize(){
  m_labelOut = m_isData ? "_data" : "_mc" ;
  m_labelOut = m_is_pPb ? m_labelOut + "_pPb" : m_labelOut + "_pp";

  auto checkWriteDir = []( const char* c_dirOut ){
    boost::filesystem::path dir( c_dirOut );  
    if(!(boost::filesystem::exists(dir))){
      std::cout<< c_dirOut << " doesn't Exist."<<std::endl;
      if (boost::filesystem::create_directory(dir))
	std::cout << "....Successfully Created !" << std::endl;
    }
  };

  // Check if the directories exist.
  // If they don't, create them
  m_dirOut   = "output";
  checkWriteDir( m_dirOut.c_str() );
  m_dirOut   += "/output" + m_labelOut;
  checkWriteDir( m_dirOut.c_str() );

  m_rootFname = m_dirOut + "/myOut" + m_labelOut + ".root";
  
  std::cout << "fNameIn/Out: " << m_rootFname << std::endl;
}

bool DiJetAnalysis::ApplyIsolation( double Rmin, std::vector<TLorentzVector>& v_jets ){
  
  std::vector<bool> isIsolated;

  for(unsigned int iTestJet = 0; iTestJet < v_jets.size(); iTestJet++){
    for(unsigned int iSecondJet = 0; iSecondJet < v_jets.size(); iSecondJet++){   
      if( iSecondJet == iTestJet ) continue;

      if( AnalysisTools::DeltaR( v_jets.at(iTestJet),
				 v_jets.at(iSecondJet)) < Rmin &&
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

void DiJetAnalysis::ApplyCleaning( std::vector<TLorentzVector>& v_jets, 
		    std::vector<bool>& v_isCleanJet){
  for( unsigned int jn = 0; jn < v_jets.size() ; jn++ ){
    if( !v_isCleanJet.at(jn) ) v_jets.at(jn).SetPxPyPzE(0,0,0,-1 );
  }
}

void DiJetAnalysis::SaveOutputs(){
  //----------------------------------------
  //  Close the input file, 
  //  write histos to output
  //----------------------------------------
  std::cout << "fNameOut: " << m_rootFname << std::endl;
  TFile* m_fOut = new TFile( m_rootFname.c_str(),"RECREATE");
  for( auto& h  : v_hists  ) { h-> Write(); }
  for( auto& f  : v_functs ) { f-> Write(); }
  for( auto& gr : v_graphs ) { gr->Write(); }
  m_fOut->Close();
}

void DiJetAnalysis::AddHistogram( TH1* h ){
  v_hists.push_back( h );
  h->Sumw2();
  h->GetXaxis()->SetNdivisions(505);  
  h->GetYaxis()->SetNdivisions(505);  
  StyleTools::SetHStyle( h, 0, 0.6 );
}
