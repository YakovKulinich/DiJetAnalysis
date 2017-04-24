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

DiJetAnalysis::DiJetAnalysis() : DiJetAnalysis( true, true, 0)
{}

DiJetAnalysis::DiJetAnalysis( bool isData, bool is_pPb, int mcType )
  : m_isData( isData ), m_is_pPb( is_pPb ), m_mcType( mcType ),
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
  m_nPtSpectBins = 40; 
  m_ptSpectMin   = 10;
  m_ptSpectMax   = 110;

  // ---- JES/PRes/Etc ----- 
  m_nEtaForwardBinsFine   = 11;
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
  m_nMinEntriesExpFit  = 20;

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
  styleTool->SetHStyle( h, 0 );

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

bool DiJetAnalysis::GetDiJets( const std::vector <TLorentzVector>& v_jets, 
			       TLorentzVector& jet1,
			       TLorentzVector& jet2 ){
  bool haveFirstJet   = false;
  bool haveSecondJet  = false;

  for( const auto& jet : v_jets ){
    if( !haveFirstJet && jet.Pt() > 0 ){
      jet1 = jet;
      haveFirstJet  = true;
    } else if( haveFirstJet && jet.Pt() > 0 ){
      jet2 = jet;
      haveSecondJet = true;
      break;
    }
  }

  // make sure we have two jets
  if( !haveFirstJet || !haveSecondJet )
    { return false; }

  // Require one of the jets to be forward
  if( !anaTool->IsForward( jet1.Eta() ) &&
      !anaTool->IsForward( jet2.Eta() ) )
    { return false; }
  
  // met all cuts, have two jets
  return true;
}

double DiJetAnalysis::AnalyzeDeltaPhi
( THnSparse* hn,
  const std::vector <TLorentzVector>& v_jets,
  double weight  ){

  TLorentzVector jet1, jet2;

  if( !GetDiJets( v_jets, jet1, jet2 ) )
    { return -1; }

  double jet1_pt  = jet1.Pt()/1000.;
  double jet1_phi = jet1.Phi();
  double jet1_eta = jet1.Eta();
  
  double jet2_pt  = jet2.Pt()/1000.;
  double jet2_phi = jet2.Phi();
  double jet2_eta = jet2.Eta();      
  
  double deltaPhi = anaTool->DeltaPhi( jet2_phi, jet1_phi );
  
  std::vector< double > x;
  x.resize( hn->GetNdimensions() );

  // some min pt cut
  for( int pt2Bin = 1; pt2Bin < hn->GetAxis(3)->GetNbins(); pt2Bin++ ){  

    if( jet2_pt < hn->GetAxis(3)->GetBinLowEdge ( pt2Bin ) ){ break; }
    
    x[0] = jet1_eta;  
    x[1] = jet2_eta;
    x[2] = jet1_pt ;
    x[3] = hn->GetAxis(3)->GetBinCenter( pt2Bin );
    x[4] = deltaPhi;
      
    hn->Fill( &x[0], weight );
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

//---------------------------
//       Plotting
//---------------------------
void DiJetAnalysis::PlotDeltaPhi( std::vector< THnSparse* >& vhn,
				  const std::vector< std::string >& vLabel,
				  const std::string& type1,
				  const std::string& type2 ){
  std::vector< TH1* > vDphi;
  std::vector< TF1* > vFits;

  std::string mcType = !type1.empty() ? "_" + type1 : "" ;
  
  for( uint iG = 0; iG < vLabel.size(); iG++ ){
    THnSparse* hn = vhn[iG];
    std::string label = vLabel[iG];
    
    TAxis* eta1Axis = hn->GetAxis(0); TAxis* pt1Axis = hn->GetAxis(2);
    TAxis* eta2Axis = hn->GetAxis(1); TAxis* pt2Axis = hn->GetAxis(3); 

    int nEta1Bins = eta1Axis->GetNbins();
    int nEta2Bins = eta2Axis->GetNbins();
    int nPt1Bins  =  pt1Axis->GetNbins();
    int nPt2Bins  =  pt2Axis->GetNbins();
    
    for( int eta1Bin = 1; eta1Bin <= nEta1Bins; eta1Bin++ ){
      eta1Axis->SetRange( eta1Bin, eta1Bin );
      for( int eta2Bin = 1; eta2Bin <= nEta2Bins; eta2Bin++ ){
	eta2Axis->SetRange( eta2Bin, eta2Bin );
	
	for( int pt1Bin = 1; pt1Bin <= nPt1Bins; pt1Bin++ ){
	  pt1Axis->SetRange( pt1Bin, pt1Bin );
	  for( int pt2Bin = 1; pt2Bin <= nPt2Bins; pt2Bin++ ){
	    pt2Axis->SetRange( pt2Bin, nPt2Bins );

	    if( !anaTool->IsForward( eta1Axis->GetBinCenter( eta1Bin ) ) &&
		!anaTool->IsForward( eta2Axis->GetBinCenter( eta2Bin ) ) )
	      { continue; }
	    if( eta1Axis->GetBinLowEdge( pt1Bin ) < 
		eta2Axis->GetBinLowEdge( pt2Bin ) )
	      { continue; }
	    
	    // Take projection onto the dPhi axis
	    TH1* hDphi = hn->Projection( 4 );
	    styleTool->SetHStyle( hDphi, 0 );
	    vDphi.push_back( hDphi );

	    double eta1Low, eta1Up, eta2Low, eta2Up;
	    double pt1Low , pt1Up , pt2Low , pt2Up;

	    anaTool->GetBinRange
	      ( eta1Axis, eta1Bin, eta1Bin, eta1Low, eta1Up );
	    anaTool->GetBinRange
	      ( eta2Axis, eta2Bin, eta2Bin, eta2Low, eta2Up );
	    anaTool->GetBinRange
	      ( pt1Axis, pt1Bin, pt1Bin, pt1Low, pt1Up );
	    anaTool->GetBinRange
	      ( pt2Axis, pt2Bin, pt2Bin, pt2Low, pt2Up );

	    hDphi->SetName
	      ( Form( "h_dPhi%s_%s_%s_%s_%s_%s",
		      mcType.c_str(),
		      anaTool->GetName( eta1Low, eta1Up, "Eta1").c_str(),
		      anaTool->GetName( eta2Low, eta2Up, "Eta2").c_str(),
		      anaTool->GetName( pt1Low , pt1Up , "Pt1" ).c_str(),
		      anaTool->GetName( pt2Low , pt2Low, "Pt2" ).c_str(),
		      label.c_str() ) );

	    TCanvas c( "c", hDphi->GetName(), 800, 600 );

	    hDphi->Draw();
	    if( hDphi->GetEntries() )
	      { hDphi->Scale( 1./hDphi->Integral() ); }
	    hDphi->GetYaxis()->SetTitle("Normalized Count");
	    hDphi->SetTitle("");
	    
	    drawTool->DrawLeftLatex
	      ( 0.13, 0.87,anaTool->GetLabel( eta1Low, eta1Up, "#eta_{1}" ) );
	    drawTool->DrawLeftLatex
	      ( 0.13, 0.82,anaTool->GetLabel( eta2Low, eta2Up, "#eta_{2}" ) );
	    drawTool->DrawLeftLatex
	      ( 0.13, 0.76,anaTool->GetLabel( pt1Low, pt1Up, "#it{p}_{T}^{1}" ) );
	    drawTool->DrawLeftLatex
	      ( 0.13, 0.69,anaTool->GetLabel( pt2Low, pt2Low, "#it{p}_{T}^{2}" ) );
	    drawTool->DrawLeftLatex( 0.13, 0.62, label );

	    if( m_isData ){
	      drawTool->DrawAtlasInternalDataRight( 0, 0, m_is_pPb );
	    } else {
	      drawTool->DrawRightLatex( 0.88, 0.82, type1 );
	      drawTool->DrawAtlasInternalMCRight( 0, 0, type2 );
	    }
	    
	    TF1* fit = anaTool->FitDphi( hDphi );
	    styleTool->SetHStyle( fit, 0 );
	    vFits.push_back( fit );
	    
	    fit->Draw("same");

	    SaveAsROOT( c, hDphi->GetName() );
	  }
      	}
      }      
    }
  } // end loop over iG
}

//---------------------------
//        Drawing
//---------------------------
void DiJetAnalysis::DrawCanvas( std::vector< TH1* >& vHIN,
				const std::string& type1,
				const std::string& type2,
				int spacing ){
  TCanvas c("c","c",800,600);
  
  TLegend leg(0.68, 0.64, 0.99, 0.77);
  styleTool->SetLegendStyle( &leg );
  leg.SetFillStyle(0);

  int style = 0;

  // for situations where dont want to
  // plot every single bin 
  // plot every n on canvas
  int dX = vHIN.size()/spacing; // plot every n
  for( unsigned int xRange = 0;
       xRange < vHIN.size();
       xRange += dX){
    styleTool->SetHStyle( vHIN[ xRange], style++ );
    leg.AddEntry( vHIN[ xRange ], vHIN[ xRange ]->GetTitle() );
    vHIN[ xRange ]->SetTitle("");
    vHIN[ xRange ]->Draw("epsame");
    SetMinMax( vHIN[ xRange ], type1, type2 );
  }

  leg.Draw();
  
  double y0 = GetLineHeight( type1 );
  
  double xMin = vHIN.front()->GetXaxis()->GetXmin();
  double xMax = vHIN.front()->GetXaxis()->GetXmax();
  
  TLine line( xMin, y0, xMax, y0);
  line.Draw();

  if( !m_isData)
    { drawTool->DrawAtlasInternalMCRight( 0, 0, m_mcTypeLabel ); } 
  else if( m_isData)
    { drawTool->DrawAtlasInternalDataRight( 0, 0, m_is_pPb ); } 

  SaveAsAll( c, type1, type2 );
}

void DiJetAnalysis::DrawCanvas( std::vector< TH1* >& vHIN,
				const std::string& type1,
				const std::string& type2,
				bool logY ){
  TCanvas c("c","c",800,600);
  if( logY ) { c.SetLogy(); }
  
  TLegend leg(0.68, 0.64, 0.99, 0.77);
  styleTool->SetLegendStyle( &leg );
  leg.SetFillStyle(0);

  double max = -1;
  for( auto& h : vHIN ){
    max = h->GetMaximum() > max ? h->GetMaximum() : max;
  }

  if( logY ){
    double power = log10(max);
    power = std::ceil(power); 
    max   = pow( 10, power );
  }
  
  int style = 0;  
  for( auto& h : vHIN ){
    styleTool->SetHStyle( h, style++ );
    leg.AddEntry( h, h->GetTitle() );
    h->SetTitle("");
    h->Draw("epsame");
    h->SetMaximum( max );
    h->SetMinimum( 1. );
  }

  leg.Draw();
  
  if( !m_isData)
    { drawTool->DrawAtlasInternalMCRight( 0, 0, m_mcTypeLabel ); } 
  else if( m_isData)
    { drawTool->DrawAtlasInternalDataRight( 0, 0, m_is_pPb ); } 

  SaveAsAll( c, type1, type2 ); 
}


void DiJetAnalysis::DrawCanvas( std::vector< TGraphAsymmErrors* >& vGIN,
				const std::string& type,
				const std::string& title,
				double xMin, double xMax ){
  TCanvas c("c","c",800,600);
  styleTool->SetCStyleEff( c, xMin, m_effMin, xMax, m_effMax, title );
   
  TLegend leg(0.64, 0.20, 0.95, 0.34);
  styleTool->SetLegendStyle( &leg );
  leg.SetFillStyle(0);
  
  int style = 0;  
  for( auto& gr : vGIN ){
    styleTool->SetHStyle( gr, style++ );
    leg.AddEntry( gr, gr->GetTitle() );
    gr->SetTitle("");
    gr->Draw("epsame");
  }

  leg.Draw();
  
  TLine line( xMin, 1, xMax, 1);
  line.Draw();

  if( !m_isData)
    { drawTool->DrawAtlasInternalMCRight( 0, 0, m_mcTypeLabel ); } 
  else if( m_isData)
    { drawTool->DrawAtlasInternalDataRight( 0, 0, m_is_pPb ); } 

  SaveAsAll( c, type ); 
}


//===== MinMax and line drawing =====

void DiJetAnalysis::
SetMinMax( TH1* h1, const std::string& type1, const std::string& type2 ){
  // JES JER
  if( !type1.compare("recoTruthRpt") ){ 
    if( !type2.compare("mean") ){ // sigma
      h1->SetMaximum(1.25);
      h1->SetMinimum(0.75);
    } else if( !type2.compare("sigma") ){ // sigma
      h1->SetMaximum(0.34);
      h1->SetMinimum(0.);
    }
  }
  // ANGLES
  else if( !type1.compare("recoTruthDeta") ||
	   !type1.compare("recoTruthDphi") ) { 
    if( !type2.compare("mean") ){ // mean
      h1->SetMaximum(0.075);      
      h1->SetMinimum(-0.075);
    } else if( !type2.compare("sigma") ){ // sigma
      h1->SetMaximum(0.056);
      h1->SetMinimum(0.);
    } 
  } 
}

double DiJetAnalysis::GetLineHeight( const std::string& type ){
  double y0 = 0;
  
  if( !type.compare("recoTruthRpt") ){ // JES/JER
    y0 = 1;
    y0 = 1;
  } else if( !type.compare("rEta") ||
	     !type.compare("recoTruthDphi") ) { // ANGLES
    y0 = 0;
    y0 = 0;
  }

  return y0;
}


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
