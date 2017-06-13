#include <TROOT.h>
#include <TEnv.h>
#include <TGraphAsymmErrors.h>
#include <TLorentzVector.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TF1.h>
#include <THnSparse.h>

#include <iostream>
#include <boost/assign.hpp>

#include "DiJetAnalysisBoth.h"
#include "DeltaPhiProj.h"

DiJetAnalysisBoth::DiJetAnalysisBoth() : DiJetAnalysisBoth( false, false ) {}

DiJetAnalysisBoth::DiJetAnalysisBoth( bool is_pPb, bool isReco )
  : DiJetAnalysis( is_pPb ), m_isReco( isReco ) {}

DiJetAnalysisBoth::~DiJetAnalysisBoth(){}

void DiJetAnalysisBoth::Initialize(){
  m_system  = m_is_pPb ? "pPb"  : "pp" ;
  m_mcLevel = m_isReco ? "reco" : "truth";

  // get list of mc used
  m_vMC = anaTool->vectorise( GetConfig()->GetValue( "usedMCs", "" ) , " " );

  m_dirOut   = "output/all/both";
  m_labelOut = "both_" + m_system + "_" + m_mcLevel; 
  
  // labels for the various mc
  for( auto & mc : m_vMC ){
    m_vMClabels.push_back
      ( GetConfig()->GetValue(Form("mcTypeLabel.2015.pp.%s",mc.c_str()),"" ));
  }
}

void DiJetAnalysisBoth::PlotDphiTogether(){
  // Check if the directories exist.
  // If they don't, create them
  std::string outDir = "output";
  anaTool->CheckWriteDir( outDir.c_str() );
  outDir += "/all";
  anaTool->CheckWriteDir( outDir.c_str() );
  outDir += "/both";
  anaTool->CheckWriteDir( outDir.c_str() );

  TAxis* axis0 = m_dPP->GetTAxis( 0 );
  TAxis* axis1 = m_dPP->GetTAxis( 1 );
  TAxis* axis2 = m_dPP->GetTAxis( 2 );
  
  int nAxis0Bins = axis0->GetNbins();
  int nAxis1Bins = axis1->GetNbins();
  int nAxis2Bins = axis2->GetNbins();

  std::string trigger = "All";
  
  TFile* fIn_data  = TFile::Open( Form("output/output_%s_data/c_myOut_%s_data.root",
				      m_system.c_str(), m_system.c_str() ) );

  std::vector< TFile* > vFinMC;  
  for( auto & mc : m_vMC ){
    std::string labelOut = m_system + "_mc_" + mc;
    std::string fNameIn  = Form( "output/output_%s/c_myOut_%s.root", labelOut.c_str(), labelOut.c_str() );
    vFinMC.push_back( new TFile( fNameIn.c_str(), "read" ) );
  }

  uint nMC = m_vMC.size();
  
  TFile* fOut = new TFile( Form("%s/c_myOut_both.root", m_dirOut.c_str() ) ,"recreate");
  
  for( int axis0Bin = 1; axis0Bin <= nAxis0Bins; axis0Bin++ ){
    // set ranges
    double axis0Low, axis0Up;
    anaTool->GetBinRange
      ( axis0, axis0Bin, axis0Bin, axis0Low, axis0Up );
    
    for( int axis1Bin = 1; axis1Bin <= nAxis1Bins; axis1Bin++ ){
      double axis1Low, axis1Up;
      anaTool->GetBinRange
	( axis1, axis1Bin, axis1Bin, axis1Low, axis1Up );

      std::vector< TH1* > vHw;
      for( int axis2Bin = 1; axis2Bin <= nAxis2Bins; axis2Bin++ ){
	double axis2Low, axis2Up;
	anaTool->GetBinRange
	  ( axis2, axis2Bin, axis2Bin, axis2Low, axis2Up );

	std::string hTagW =
	Form ("%s_%s_%s",
	      anaTool->GetName( axis0Low, axis0Up, m_dPP->GetAxisName(0) ).c_str(),
	      anaTool->GetName( axis1Low, axis1Up, m_dPP->GetAxisName(1) ).c_str(),
	      anaTool->GetName( axis2Low, axis2Up, m_dPP->GetAxisName(2) ).c_str() );
	
	TCanvas c( "c", "", 800, 600 );
	
	TLegend legW( 0.33, 0.13, 0.87, 0.26 );
	styleTool->SetLegendStyle( &legW );
	legW.SetNColumns(2);	

	int style = 0;

	std::string hNameW_data = Form( "h_dPhi_%s_%s", hTagW.c_str(), trigger.c_str() );
	std::string hNameW_mc   = Form("h_dPhi_%s_%s", m_mcLevel.c_str(), hTagW.c_str() );
	  
	TH1* hW_data = static_cast< TH1D* >
	  ( fIn_data->Get( hNameW_data.c_str() ) );
	vHw.push_back( hW_data );
	styleTool->SetHStyle( hW_data, style++ );
	legW.AddEntry( hW_data, "Data" );
	hW_data->Draw("ep same x0");
	hW_data->SetMarkerSize( hW_data-> GetMarkerSize() * 1.5 );
	
	for( uint iMC = 0; iMC < nMC; iMC++ ){
	  TH1* hW_mc = static_cast< TH1D* >( vFinMC[iMC]->Get( hNameW_mc.c_str() ) );
	  styleTool->SetHStyle( hW_mc, style++ );
	  legW.AddEntry( hW_mc, m_vMClabels[iMC].c_str() );
	  hW_mc->Draw("ep same x0");
	  hW_mc->SetMarkerSize( hW_mc-> GetMarkerSize() * 1.5 );
	  vHw.push_back( hW_mc );
	}
	
	legW.Draw("same");
	
	drawTool->DrawTopLeftLabels
	( m_dPP, axis0Low, axis0Up, axis1Low, axis1Up,
	  axis2Low, axis2Up, 0, 0, 0.8 );

	drawTool->DrawRightLatex
	  ( 0.88, 0.76, Form("MC %s level", m_mcLevel.c_str() ) );
	
	drawTool->DrawAtlasInternalDataRight( 0, 0, m_is_pPb );

	SaveAsAll( c, Form("h_dPhi_%s", hTagW.c_str() ) );
      } // end loop over axis2

      for( auto& hW: vHw ){ delete hW; }
    } // end loop over axis1
  } // end loop over axis0

  std::cout << "DONE! Closing " << fOut->GetName() << std::endl;
  fOut->Close();
  std::cout << "......Closed  " << fOut->GetName() << std::endl;
}
