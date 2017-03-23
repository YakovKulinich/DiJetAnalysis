#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TROOT.h>
#include <TEnv.h>
#include <TCanvas.h>

#include <cmath>
#include <iostream>
#include <sstream>
#include <boost/format.hpp>
#include <boost/filesystem.hpp>

#include "MyRoot.h"

#include "DiJetAnalysis.h"

DiJetAnalysis::DiJetAnalysis() : DiJetAnalysis( true, true )
{}

DiJetAnalysis::DiJetAnalysis( bool isData, bool is_pPb )
  : m_isData( isData ), m_is_pPb( is_pPb ),
    m_fIn(NULL), m_fOut(NULL), m_tree(NULL)
{
  //========== Set Histogram Binning =============
  // Common for all analysis
  
  // Triggers and Spectra
  m_nPtSpectBins = 50; 
  m_ptSpectMin   = 10;
  m_ptSpectMax   = 2 * m_nPtSpectBins + m_ptSpectMin;
  
  // Eta-Phi Maps
  m_nEtaBins = 100; 
  m_etaMin   = constants::ETAMIN;
  m_etaMax   = constants::ETAMAX;
  
  m_nPhiBins = 64; 
  m_phiMin   = -constants::PI;
  m_phiMax   = constants::PI; 
  
  m_ptWidth  = 2;
  m_ptMin    = 10;
  m_ptMax    = 100;
  m_nPtBins  = (m_ptMax - m_ptMin)/m_ptWidth;
  
  m_nEtaForwardBinsFine   = 12;
  m_etaForwardMin   = -constants::FETAMAX;
  m_etaForwardMax   = -constants::FETAMIN;

  m_varEtaBinning.push_back( -4.4 );
  m_varEtaBinning.push_back( -4.0 );
  m_varEtaBinning.push_back( -3.6 );
  m_varEtaBinning.push_back( -3.3 );
  m_nVarEtaBins = static_cast<int>( m_varEtaBinning.size() ) - 1;
  
  //==================== Cuts ====================
}

DiJetAnalysis::~DiJetAnalysis(){}

void DiJetAnalysis::Initialize(){
  m_labelOut = m_isData ? "data" : "mc" ;
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
  m_dirOut   += "/output_" + m_labelOut;
  checkWriteDir( m_dirOut.c_str() );

  m_rootFname = m_dirOut + "/myOut_" + m_labelOut + ".root";
  
  std::cout << "fNameIn/Out: " << m_rootFname << std::endl;
}

//---------------------------
//       Fill Tree
//---------------------------
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

//---------------------------
//       Plotting 
//---------------------------
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

void DiJetAnalysis::SaveAsPdfPng( const TCanvas& c,
				  const std::string& label1,
				  const std::string& label2,
				  const std::string& axis1,
				  double min1, double max1,
				  const std::string& axis2 ,
				  double min2, double max2 ){
  std::stringstream ss;
  ss << m_dirOut + "/";
  
  if( !label1.empty() ){ ss << boost::format("%s") % label1; }
  if( !label2.empty() ){ ss << boost::format("_%s") % label2; }
  ss << "_" + m_labelOut;
  if( !axis1.empty() ){
    ss << boost::format("_%2.0f_%s_%2.0f") % min1 % axis1 % max1;
  }
  if( !axis2.empty() ){
    ss << boost::format("_%2.0f_%s_%2.0f") % min2 % axis2 % max2;
  }

  std::string sPdf = ss.str() + ".pdf";
  std::string sPng = ss.str() + ".png";

  c.SaveAs( sPdf.c_str() );
  c.SaveAs( sPng.c_str() );
}

void DiJetAnalysis::SaveAsROOT( const TCanvas& c,
			        const std::string& label1,
				const std::string& label2,
				const std::string& axis1,
				double min1, double max1,
				const std::string& axis2 ,
				double min2, double max2 ){
  std::stringstream ss;
  ss << "c_";
  
  if( !label1.empty() ){ ss << boost::format("%s") % label1; }
  if( !label2.empty() ){ ss << boost::format("_%s") % label2; }
  ss << "_" + m_labelOut;
  if( !axis1.empty() ){
    ss << boost::format("_%2.0f_%s_%2.0f") % min1 % axis1 % max1;
  }
  if( !axis2.empty() ){
    ss << boost::format("_%2.0f_%s_%2.0f") % min2 % axis2 % max2;
  }

  c.Write( ss.str().c_str() );
}

void DiJetAnalysis::SaveAsAll( const TCanvas& c,
			       const std::string& label1,
			       const std::string& label2,
			       const std::string& axis1,
			       double min1, double max1,
			       const std::string& axis2 ,
			       double min2, double max2 ){  
  SaveAsPdfPng( c, label1, label2, axis1, min1, max1, axis2, min2, max2 );
  SaveAsROOT  ( c, label1, label2, axis1, min1, max1, axis2, min2, max2 );
}

//---------------------------
//       Tools
//---------------------------
void DiJetAnalysis::AddHistogram( TH1* h ){
  v_hists.push_back( h );
  h->Sumw2();
  h->GetXaxis()->SetNdivisions(505);  
  h->GetYaxis()->SetNdivisions(505);  
  StyleTools::SetHStyle( h, 0, 0.6 );
}

void DiJetAnalysis::FitGaussian( TH1* hProj, TF1* fit ){
  // fit once 
  double mean   = hProj->GetMean();
  double rms    = hProj->GetRMS();

  double fitmin = mean - 2.0 * rms;
  double fitmax = mean + 2.0 * rms;
      
  hProj->Fit( fit->GetName(), "NQR", "", fitmin, fitmax );

  // fit second time with better parameters
  mean   = fit->GetParameter(1);
  rms    = fit->GetParameter(2);

  fitmin = mean - 1.8 * rms;
  fitmax = mean + 2.2 * rms;

  fit->SetRange( fitmin, fitmax );
  
  hProj->Fit( fit->GetName(), "NQR", "", fitmin, fitmax );
}

std::string DiJetAnalysis::GetEtaLabel( double etaMin, double etaMax ){

  std::stringstream ss;
  
  if( m_is_pPb ){
    ss << boost::format("%3.1f<#eta<%3.1f") % etaMin % etaMax;
  } else {
    if( abs(etaMin) > abs(etaMax) )
      ss << boost::format("%3.1f<|#eta|<%3.1f") % abs(etaMax) % abs(etaMin);
    else
      ss << boost::format("%3.1f<|#eta|<%3.1f") % abs(etaMin) % abs(etaMax);
  }
  
  return ss.str();
}

double DiJetAnalysis::AdjustEtaForPP( double jetEta ){
  if( m_is_pPb ) return jetEta;
  return jetEta  > 0 ? -1 * jetEta  : jetEta;
}
