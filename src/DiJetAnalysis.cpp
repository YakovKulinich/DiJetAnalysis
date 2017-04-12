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
#include <boost/assign.hpp>

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
  
  // -------- maps ---------
  m_nEtaMapBins = 100; 
  m_etaMapMin   = constants::ETAMIN;
  m_etaMapMax   = constants::ETAMAX;
  
  m_nPhiMapBins = 64; 
  m_phiMapMin   = -constants::PI;
  m_phiMapMax   = constants::PI; 

  m_nPtMapBins  = 45;
  m_ptMapMin    = 10;
  m_ptMapMax    = 100;
  
  // -------- spect --------
  m_nPtSpectBins = 50; 
  m_ptSpectMin   = 10;
  m_ptSpectMax   = 2 * m_nPtSpectBins + m_ptSpectMin;

  // ---- JES/PRes/Etc ----- 
  m_nEtaForwardBinsFine   = 12;
  m_etaForwardMin   = -constants::FETAMAX;
  m_etaForwardMax   = -constants::FETAMIN;

  // ---- forward eta binning ---
  boost::assign::push_back( m_varFwdEtaBinning )
    ( -4.4 )( -4.0 )( -3.6 )( -3.3 );
  m_nVarFwdEtaBins = m_varFwdEtaBinning.size() - 1;

  // -------- eff ---------
  m_effMin = 0.;
  m_effMax = 1.3;

  // -------- dphi- --------
  m_nDphiPtBins = 5; 
  m_nDphiPtMin  = 10;
  m_nDphiPtMax  = 60;
  
  m_nDphiDphiBins = 60;
  m_nDphiDphiMin  = 0;
  m_nDphiDphiMax  = constants::PI;

  boost::assign::push_back( m_varEtaBinning )
    ( -4.4 )( -3.3 )( -2.8 )( 0 )( 2.8 )( 3.3 )( 4.4 );
  m_nVarEtaBins = m_varEtaBinning.size() - 1;

  // --- variable pt binning ---
  boost::assign::push_back( m_varPtBinning )
    ( 20 )( 30 )(40)( 200 );
  m_nVarPtBins = m_varPtBinning.size() - 1;
  
  //==================== Cuts ====================
  m_nMinEntriesGausFit = 20;
  m_ptFitMin           = 20.;
  
  m_dPhiThirdJetFraction = 0.4;

  //========= common tools ========
  anaTool   = new CT::AnalysisTools();
  drawTool  = new CT::DrawTools();
  styleTool = new CT::StyleTools();
}

DiJetAnalysis::~DiJetAnalysis(){
  delete anaTool; delete drawTool; delete styleTool;
}

void DiJetAnalysis::Initialize(){
  m_labelOut = m_isData ? "data" : "mc" ;
  m_labelOut = m_is_pPb ? m_labelOut + "_pPb" : m_labelOut + "_pp";

  // Check if the directories exist.
  // If they don't, create them
  auto checkWriteDir = []( const char* c_dirOut ){
    boost::filesystem::path dir( c_dirOut );  
    if(!(boost::filesystem::exists(dir))){
      std::cout<< c_dirOut << " doesn't Exist."<<std::endl;
      if (boost::filesystem::create_directory(dir))
	std::cout << "....Successfully Created !" << std::endl;
    }
  };

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
  styleTool->SetHStyle( h, 0, CT::StyleTools::hSS );

  TH1* h1 = dynamic_cast< TH1* >(h);
  if( h1 ){ h->GetXaxis()->SetNdivisions(505); }

  TH2* h2 = dynamic_cast< TH2* >(h);
  if( h2 ){ h->GetYaxis()->SetNdivisions(505); }

  TH3* h3 = dynamic_cast< TH3* >(h);
  if( h3 ){ h->GetZaxis()->SetNdivisions(505); }
}

void DiJetAnalysis::AddHistogram( THnSparse* hn ){
  v_hns.push_back( hn );
  hn->Sumw2();
  hn->GetAxis(0)->SetTitle("#eta_{1}");
  hn->GetAxis(1)->SetTitle("#eta_{2}");
  hn->GetAxis(2)->SetTitle("#it{p}_{T}^{1}");
  hn->GetAxis(3)->SetTitle("#it{p}_{T}^{2}");
  hn->GetAxis(4)->SetTitle("#Delta#phi");  
}


void DiJetAnalysis::SaveOutputsFromTree(){
  //----------------------------------------
  //  Close the input file, 
  //  write histos to output
  //----------------------------------------
  std::cout << "fNameOut: " << m_rootFname << std::endl;
  TFile* m_fOut = new TFile( m_rootFname.c_str(),"RECREATE");
  for( auto& h  : v_hists  ){ h-> Write(); }
  for( auto& f  : v_functs ){ f-> Write(); }
  for( auto& gr : v_graphs ){ gr->Write(); }
  for( auto& hn : v_hns    ){ hn->Write(); }
  m_fOut->Close();
}

//---------------------------
//       Analysis
//---------------------------

double DiJetAnalysis::AnalyzeDeltaPhi
( THnSparse* hn, const std::vector <TLorentzVector>& v_jets ){
  
  const TLorentzVector* jet1 = &v_jets.at(0);
  const TLorentzVector* jet2 = &v_jets.at(1);

  // because removing from vectors is expensive, I set
  // "bad" or unmached jets to have (0,0,0,-1) 4vector.
  if( jet2->Pt() == 0 || jet1->Pt() == 0 ) return -1;

  unsigned int nJets = v_jets.size();
  
  // if we have only 2 jets, we consider the 
  // "third" to have Pt=0
  double followingJetPt = 0;

  if( nJets > 2 )
    followingJetPt = v_jets.at(2).Pt()/1000.;

  double jet1_pt  = jet1->Pt()/1000.;
  double jet1_phi = jet1->Phi();
  double jet1_eta = jet1->Eta();
  
  double jet2_pt  = jet2->Pt()/1000.;
  double jet2_phi = jet2->Phi();
  double jet2_eta = jet2->Eta();      

  double ptAvg    = ( jet2_pt + jet1_pt )/2;

  double deltaPhi = anaTool->DeltaPhi( jet2_phi, jet1_phi );
  
  // some cuts for phi inter calibration
  if( followingJetPt > ptAvg * m_dPhiThirdJetFraction ) return -1;

  std::vector< double > x;
  x.resize( hn->GetNdimensions() );

  // some min pt cut
  for( int pt2Bin = 1; pt2Bin < hn->GetAxis(3)->GetNbins(); pt2Bin++ ){  

    if( jet2_pt < hn->GetAxis(3)->GetBinUpEdge ( pt2Bin ) &&
	jet2_pt < hn->GetAxis(3)->GetBinLowEdge( pt2Bin ) ){ break; }
    
    x[0] = jet1_eta;  
    x[1] = jet2_eta;
    x[2] = jet1_pt ;
    x[3] = hn->GetAxis(3)->GetBinCenter( pt2Bin );
    x[4] = deltaPhi;
      
    hn->Fill( &x[0], 1 );
  }
  
  return deltaPhi;
}

void DiJetAnalysis::ApplyIsolation( double Rmin, std::vector<TLorentzVector>& v_jets ){
  
  std::vector<bool> isIsolated;

  for(unsigned int iTestJet = 0; iTestJet < v_jets.size(); iTestJet++){
    for(unsigned int iSecondJet = 0; iSecondJet < v_jets.size(); iSecondJet++){   
      if( iSecondJet == iTestJet ) continue;

      if( anaTool->DeltaR( v_jets.at(iTestJet),
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
//       Tools
//---------------------------

void DiJetAnalysis::
ProjectAndFitGaus( TH3* h3,
		   TH1* h1Mean, TH1* h1Sigma,
		   int etaBinLow, int etaBinUp,
		   const std::string& jzn,
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
    styleTool->SetHStyle( hProj, 0, CT::StyleTools::hSS);
    v_hProj.push_back( hProj );
    hProj->SetTitle("");
    
    TF1* fit  = new TF1( Form("f_%s_%2.0f_Eta_%2.0f_%2.0f_Pt_%2.0f",
			      h3->GetName(),
			      10*std::abs(etaMin),
			      10*std::abs(etaMax),
			      ptMin, ptMax ),
			 "gaus(0)" );
    styleTool->SetHStyle( fit, 0, CT::StyleTools::hSS);
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
      drawTool->DrawAtlasInternalMCRight
	( 0, 0, CT::StyleTools::lSS, mcLabel  );
      drawTool->DrawLeftLatex
	( 0.18, 0.74, jzn.c_str(), CT::StyleTools::lSS, 1 );
    } else if( m_isData ){
      drawTool->DrawAtlasInternalDataRight
	( 0, 0, CT::StyleTools::lSS, m_is_pPb ); 
    }
    
    drawTool->DrawLeftLatex( 0.18, 0.88,
			     GetEtaLabel( etaMin, etaMax ).c_str(),
			     CT::StyleTools::lSS, 1 );
    drawTool->DrawLeftLatex( 0.18, 0.81,
			     Form("%3.0f<%s<%3.1f",
				  ptMin,
				  h1Mean->GetXaxis()->GetTitle(),
				  ptMax ),
			     CT::StyleTools::lSS, 1 );
    
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

  for( auto& h : v_hProj ){ delete h; }
  for( auto& f : v_fit   ){ delete f; }
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

double DiJetAnalysis::AdjustEtaForPP( double jetEta ){
  if( m_is_pPb ) return jetEta;
  return jetEta  > 0 ? -1 * jetEta  : jetEta;
}

std::string DiJetAnalysis::GetLabel
( double vMin,double vMax,
  const std::string& var, const std::string& unit ){
  
  std::stringstream ss;

  if( vMin != vMax ){
    ss << boost::format("%3.1f<%s<%3.1f")
      % vMin % var % vMax;
  }
  else{
    ss << boost::format("%s>%3.1f")
      % var % vMin ;  
  }
  if( !unit.empty() ){ ss << boost::format(" %s") % unit; }
  
  return ss.str();
}

std::string DiJetAnalysis::GetEtaLabel( double etaMin,
					double etaMax ){

  std::stringstream ss;
  
  if( m_is_pPb ){
    ss << boost::format("%3.1f<#eta<%3.1f") % etaMin % etaMax;
  } else {
    if( std::abs(etaMin) > std::abs(etaMax) )
      ss << boost::format("%3.1f<|#eta|<%3.1f")
	% std::abs(etaMax) % std::abs(etaMin);
    else
      ss << boost::format("%3.1f<|#eta|<%3.1f")
	% std::abs(etaMin) % std::abs(etaMax);
  }
  
  return ss.str();
}

//---------------------------
//       Plotting
//---------------------------

//---------------------------
//       Saving 
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
