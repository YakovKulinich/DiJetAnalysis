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

  m_effMin = 0.;
  m_effMax = 1.3;
  
  //==================== Cuts ====================
  m_ptFitMin = 20;
  m_nMinEntriesGausFit = 20;;
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
void DiJetAnalysis::AddHistogram( TH1* h ){
  v_hists.push_back( h );
  h->Sumw2();
  h->GetXaxis()->SetNdivisions(505);  
  h->GetYaxis()->SetNdivisions(505);  
  StyleTools::SetHStyle( h, 0, StyleTools::hSS );
}

void DiJetAnalysis::SaveOutputsFromTree(){
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

//---------------------------
//       Analysis
//---------------------------

void DiJetAnalysis::ApplyIsolation( double Rmin, std::vector<TLorentzVector>& v_jets ){
  
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

void DiJetAnalysis::
ProjectAndFitGaus( TH3* h3,
		   TH1* h1Mean, TH1* h1Sigma,
		   int etaBinLow, int etaBinUp,
		   int jzn,
		   const std::string& mcLabel){
  
  TCanvas c_proj("c_proj","c_proj",800,600);

  double etaMin = h3->GetXaxis()->GetBinLowEdge( etaBinLow );
  double etaMax = h3->GetXaxis()->GetBinUpEdge( etaBinUp );

  int  nPtBins  = h3->GetNbinsY();

  std::vector< TH1* > v_hProj;
  std::vector< TF1* > v_fit;
  
  // loop over ptBins
  for( int ptBin = 1; ptBin <= nPtBins; ptBin++ ){
    double ptMin = h3->GetYaxis()->GetBinLowEdge(ptBin);
    double ptMax = h3->GetYaxis()->GetBinUpEdge(ptBin);
      
    TH1* hProj = h3->
      ProjectionZ( Form("%s_%2.0f_Eta_%2.0f_%2.0f_Pt_%2.0f",
			h3->GetName(),
			10*std::abs(etaMin),
			10*std::abs(etaMax),
			ptMin, ptMax ),
		   etaBinLow, etaBinUp, ptBin, ptBin );
    StyleTools::SetHStyle( hProj, 0, StyleTools::hSS);
    v_hProj.push_back( hProj );
    hProj->SetTitle("");
    
    TF1* fit  = new TF1( Form("f_%s_%2.0f_Eta_%2.0f_%2.0f_Pt_%2.0f",
			      h3->GetName(),
			      10*std::abs(etaMin),
			      10*std::abs(etaMax),
			      ptMin, ptMax ),
			 "gaus(0)" );
    StyleTools::SetHStyle( fit, 0, StyleTools::hSS);
    v_fit.push_back( fit );

    if( hProj->GetEntries() < m_nMinEntriesGausFit ){ continue; }
    
    FitGaussian( hProj, fit );
        
    h1Mean->SetBinContent( ptBin, fit->GetParameter(1) );
    h1Mean->SetBinError  ( ptBin, fit->GetParError (1) );

    h1Sigma->SetBinContent ( ptBin, fit->GetParameter(2) );
    h1Sigma->SetBinError   ( ptBin, fit->GetParError (2) );
    
    hProj->Draw();
    fit->Draw("same");

    if( !m_isData ){
      DrawTools::DrawAtlasInternalMCRight
	( 0, 0, StyleTools::lSS, mcLabel  ); 
    } else if( m_isData ){
      DrawTools::DrawAtlasInternalDataRight
	( 0, 0, StyleTools::lSS, m_is_pPb ); 
    }
    
    DrawTools::DrawRightLatex( 0.88, 0.81,
			      GetEtaLabel( etaMin, etaMax ).c_str(),
			      StyleTools::lSS, 1 );
    DrawTools::DrawRightLatex( 0.88, 0.74,
			      Form("%3.0f<%s<%3.1f",
				   ptMin,
				   h1Mean->GetXaxis()->GetTitle(),
				   ptMax ),
			      StyleTools::lSS, 1 );
    if( jzn >= 0 ){
      DrawTools::DrawRightLatex( 0.88, 0.68,
				 Form("JZ%i", jzn ),
				 StyleTools::lSS, 1 );
    }
    
    c_proj.Write( Form("c_%s_%s_%2.0f_Eta_%2.0f_%2.0f_Pt_%2.0f",
		       h3->GetName(),
		       m_labelOut.c_str(),
		       std::abs(etaMin)*10,
		       std::abs(etaMax)*10,
		       ptMin, ptMax ) );
    
    SaveAsROOT( c_proj,
	        h3->GetName(), "",
		"Eta", std::abs(etaMin)*10, std::abs(etaMax)*10,
		"Pt" , ptMin, ptMax );    
  } // end loop over ptBins
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

  fitmin = mean - 2.0 * rms;
  fitmax = mean + 2.0 * rms;

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
