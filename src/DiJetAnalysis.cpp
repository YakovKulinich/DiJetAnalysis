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
#include <TVirtualFitter.h>

#include <cmath>
#include <iostream>
#include <sstream>
#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <boost/assign.hpp>

#include "MyRoot.h"

#include "DiJetAnalysis.h"
#include "DeltaPhiProj.h"

DiJetAnalysis::DiJetAnalysis() : DiJetAnalysis( true, true, 0 )
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
  m_ptSpectMin   = 10;
  m_ptSpectMax   = 110;
  m_nPtSpectBins = ( m_ptSpectMax - m_ptSpectMin ) / 2 ; 
  
  // ---- JES/PRes/Etc ----- 
  m_etaForwardMin = -constants::FETAMAX;
  m_etaForwardMax = -constants::FYSTARMIN;
  m_nEtaForwardBinsFine = ( m_etaForwardMax - m_etaForwardMin ) / 0.2 ;
  
  // ---- forward eta binning ---
  boost::assign::push_back( m_varFwdEtaBinning )
    ( -constants::FYSTARMAX )( -3.9 )( -3.5 )
    ( -constants::FETAMIN )( -3.1 )( -constants::FYSTARMIN );
  m_nVarFwdEtaBins = m_varFwdEtaBinning.size() - 1;
  
  // -------- eff ---------
  m_effMin = 0.;
  m_effMax = 1.4;

  // -------- dphi- --------
  m_nDphiDphiBins = 60;
  m_dPhiDphiMin   = 0;
  m_dPhiDphiMax   = constants::PI;

  // --- variable eta/ystar binning ---
  // ystarB
  boost::assign::push_back( m_varYstarBinningA )
    ( -4.0 )( -2.7 )( -1.8 )( 0 )
    ( 1.8 )( 2.7 )( 4.0 );
  m_nVarYstarBinsA = m_varYstarBinningA.size() - 1;
  
  // ystarA
  boost::assign::push_back( m_varYstarBinningB )
    ( -4.0 )( -2.7 )( -1.8 )( 0 )
    ( 1.8 )( 2.7 )( 4.0 );
  m_nVarYstarBinsB = m_varYstarBinningB.size() - 1;
  
  // --- variable pt binning ---
  boost::assign::push_back( m_varPtBinning )
    ( 25 )( 35 )( 45 )( 90 );
  m_nVarPtBins = m_varPtBinning.size() - 1;

  // --- dPhiBins ---  
  boost::assign::push_back( m_nDphiBins )
    ( m_nVarYstarBinsA )( m_nVarYstarBinsB )
    ( m_nVarPtBins     )( m_nVarPtBins     )
    ( m_nDphiDphiBins  );
    
  boost::assign::push_back( m_dPhiMin  )
    ( 0 )( 0 )( 0 )( 0 )( m_dPhiDphiMin );

  boost::assign::push_back( m_dPhiMax  )
    ( 1 )( 1 ) ( 1 )( 1 )( m_dPhiDphiMax );
    
  m_dPhiWidthMin = 0.00;
  m_dPhiWidthMax = 0.5;

  //========= Set DeltaPhi Axes Order ============
  // The DeltaPhiProj object will have the order
  // onto which to take projections. It also knows
  // what axis is what. The object knows which axis
  // has what name and range. I.e. Pt1 is actually
  // The third axis in the THnSparse. 
  m_dPP = new DeltaPhiProj( YS1, PT1, PT2, YS2 );

  m_dPP->AddTAxis( new TAxis( m_nVarYstarBinsA, &m_varYstarBinningA[0] ) );
  m_dPP->AddTAxis( new TAxis( m_nVarYstarBinsB, &m_varYstarBinningB[0] ) );
  m_dPP->AddTAxis( new TAxis( m_nVarPtBins    , &m_varPtBinning[0]     ) );
  m_dPP->AddTAxis( new TAxis( m_nVarPtBins    , &m_varPtBinning[0]     ) );
  
  //==================== Cuts ====================
  m_nMinEntriesFit = 20;

  m_dPhiThirdJetFraction = 0.4;

  //========= common tools ========
  anaTool   = new CT::AnalysisTools();
  drawTool  = new CT::DrawTools();
  styleTool = new CT::StyleTools();

  //========= config file =========
  std::string configName   = "config/configJetsFwd";
  
  std::string configDataType = m_is_pPb ? "_pPb" : "_pp";
  configName += configDataType;

  std::string configSuffix = m_isData ? "_data.cfg" : "_mc.cfg" ;
  configName += configSuffix;

  std::cout << "Reading " << configName << std::endl;

  m_config = new TEnv();
  m_config->ReadFile( configName.c_str(), EEnvLevel(0));
}

DiJetAnalysis::~DiJetAnalysis(){
  delete m_dPP;
  delete anaTool;
  delete drawTool;
  delete styleTool;
}

void DiJetAnalysis::Initialize(){  
  m_labelOut = m_is_pPb ? "pPb" : "pp" ;
  m_labelOut = m_isData ? m_labelOut + "_data" : m_labelOut + "_mc";
  
  m_dirOut   = "output";
  anaTool->CheckWriteDir( m_dirOut.c_str() );
  m_dirOut   += "/output_" + m_labelOut;
  anaTool->CheckWriteDir( m_dirOut.c_str() );

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
  hn->GetAxis(0)->SetTitle("#it{y}_{1}*");
  hn->GetAxis(1)->SetTitle("#it{y}_{2}*");
  hn->GetAxis(2)->SetTitle("#it{p}_{T}^{1}");
  hn->GetAxis(3)->SetTitle("#it{p}_{T}^{2}");
  if( hn->GetNdimensions() == 5 )
    { hn->GetAxis(4)->SetTitle("|#Delta#phi|"); }
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

bool DiJetAnalysis::GetFwdCentJets( const std::vector< TLorentzVector>& v_jets, 
				    const TLorentzVector*& jetFwd,
				    const TLorentzVector*& jetCent ){
  for( auto& jet : v_jets ){
    if( !jetFwd && IsForwardDetector( jet.Eta() ) && jet.Pt() > 0 )
      { jetFwd  = &jet; }
    else if( !jetCent && IsCentralDetector( jet.Eta() ) && jet.Pt() > 0 )
      { jetCent = &jet; }
    if( jetFwd && jetCent )
      { break; }
  }
  // check we have at least one jet 
  if( !jetFwd && !jetCent ){ return false; }
  return true;
}

bool DiJetAnalysis::GetDiJets( const std::vector< TLorentzVector >& v_jets, 
			       const TLorentzVector*& jet1,
			       const TLorentzVector*& jet2 ){
  for( auto& jet : v_jets ){
    if( !jet1 && jet.Pt() > 0 )
      { jet1 = &jet; } 
    else if( jet1 && jet.Pt() > 0 )
      { jet2 = &jet;
	break; }
  }
  
  // make sure we have two jets
  if( !jet1 || !jet2 )
    { return false; }
  return true;
}

double DiJetAnalysis::AnalyzeDeltaPhi( THnSparse* hn, THnSparse* hnNent,
				       const std::vector< TLorentzVector >& v_jets,
				       double weightIn, WeightFcn weightFcn ){
  const TLorentzVector* jet1 = NULL; const TLorentzVector* jet2 = NULL;

  if( !GetDiJets( v_jets, jet1, jet2 ) )
    { return -1; }

  double jet1_pt    = jet1->Pt()/1000.;
  double jet1_eta   = jet1->Eta();
  double jet1_phi   = jet1->Phi();
  double jet1_ystar = GetYstar( *jet1 );

  double jet2_pt    = jet2->Pt()/1000.;
  double jet2_ystar = GetYstar( *jet2 );
  
  double deltaPhi = anaTool->DeltaPhi( *jet2, *jet1 );
  
  std::vector< double > x;
  x.resize( hn->GetNdimensions() );
    
  // wont change unless we have a weightFcn
  // and then it varys depending on eta, phi, pt.
  double weight = weightFcn ? weightFcn( jet1_eta, jet1_phi, jet1_pt ) : weightIn;   

  x[0] = jet1_ystar;  
  x[1] = jet2_ystar;
  x[2] = jet1_pt ;
  x[3] = jet2_pt ;
  x[4] = deltaPhi;
  hn    ->Fill( &x[0], weight );
  hnNent->Fill( &x[0], 1 );

  // just to check everything was sorted.
  // can get rid of later
  if( jet1_pt < jet2_pt ){
    for( auto& jet : v_jets ){
      std::cout << jet.Pt()/1000. << " " << jet.Eta() << std::endl;
    }
  }
  
  // for pp, fill twice. once for each side since
  // it is symmetric in pp. For pPb, continue
  if( m_is_pPb ){ return deltaPhi; }
  
  x[0] = -jet1_ystar;  
  x[1] = -jet2_ystar;
  hn    ->Fill( &x[0], weight );
  hnNent->Fill( &x[0], 1 );    
 
  return deltaPhi;
}

bool DiJetAnalysis::ApplyIsolation( std::vector<TLorentzVector>& v_jets,
				    double Rmin ){
  bool haveIsolated = false;
  std::vector<bool> isIsolated( v_jets.size(), true );
  
  for(unsigned int iTestJet = 0; iTestJet < v_jets.size(); iTestJet++){
    for(unsigned int iSecondJet = iTestJet; iSecondJet < v_jets.size(); iSecondJet++){   
      if( iSecondJet == iTestJet ) continue;

      if( anaTool->DeltaR( v_jets.at(iTestJet),
			   v_jets.at(iSecondJet)) < Rmin &&
	  v_jets.at(iSecondJet).Pt() > v_jets.at(iTestJet).Pt() * 0.5 ){       
	isIsolated[ iSecondJet ] = false;
	haveIsolated = true;
	continue;
      }
    } // end loop over iSecondJet
  } // end loop over iTestJet

  for(unsigned int iJet = 0; iJet < v_jets.size(); iJet++ ){
    if( !isIsolated.at(iJet) ) v_jets.at(iJet).SetPxPyPzE(0,0,0,-1 );
  }

  return haveIsolated; 
}

bool DiJetAnalysis::ApplyCleaning( std::vector<TLorentzVector>& v_jets, 
				   std::vector<bool>& v_isCleanJet){
  bool haveCleaned = false;
  
  for( unsigned int jn = 0; jn < v_jets.size() ; jn++ ){
    if( !v_isCleanJet.at(jn) ){
      v_jets.at(jn).SetPxPyPzE(0,0,0,-1 );
      haveCleaned = true; 
    }
  }
  
  return haveCleaned;
}

//---------------------------
//       Tools
//---------------------------
void  DiJetAnalysis::FillHistoWithJets( const TLorentzVector* jet1,
					const TLorentzVector* jet2,
					TH2* h, double weight ){
  if( jet1 )
    { h->Fill( AdjustEtaForPP( jet1->Eta() ), jet1->Pt()/1000., weight ); }
  if( jet2 )
    { h->Fill( AdjustEtaForPP( jet2->Eta() ), jet2->Pt()/1000., weight ); }
} 

double DiJetAnalysis::AdjustEtaForPP( double jetEta ){
  if( m_is_pPb ) return jetEta;
  return jetEta  > 0 ? -1 * jetEta  : jetEta;
}

double DiJetAnalysis::GetYstar( const TLorentzVector& jet )
{ return m_is_pPb ? jet.Rapidity() + constants::BETAZ : jet.Rapidity(); }

bool DiJetAnalysis::IsForwardDetector( const double& eta ){
  if( std::abs( eta ) < constants::FETAMAX &&
      std::abs( eta ) > constants::FETAMIN ){ return true; }
  return false;
}

bool DiJetAnalysis::IsCentralDetector( const double& eta )
{ return ( std::abs(eta) < constants::CETAMAX ); }

// These two are just used to see which histograms to draw
// for delta Phi plots, side doesnt really matter
// because can have forward backward correlation.
bool DiJetAnalysis::IsForwardYstar( const double& ystar )
{ return ( std::abs(ystar - constants::BETAZ) < constants::FETAMAX && 
	       std::abs(ystar - constants::BETAZ) > constants::FETAMIN ); }

bool DiJetAnalysis::IsCentralYstar( const double& ystar )
{ return ( std::abs(ystar - constants::BETAZ) < constants::CETAMAX ); }

//---------------------------
//       Plotting
//---------------------------

void DiJetAnalysis::PlotDeltaPhi( std::vector< THnSparse* >& vhn,
				  std::vector< THnSparse* >& vhnNent,
				  const std::vector< std::string >& vLabel,
				  const std::string& type1,
				  const std::string& type2 ){
  FourDTH1vector vDphiWidths;
  FourDTH1vector vDphiNent;
  
  PlotDeltaPhi( vhn, vhnNent, vDphiWidths, vDphiNent, vLabel, type1, type2 );
}

void DiJetAnalysis::PlotDeltaPhi( std::vector< THnSparse* >& vhn,
				  std::vector< THnSparse* >& vhnNent,
				  FourDTH1vector& vDphiWidths,
				  FourDTH1vector& vDphiNent,
				  const std::vector< std::string >& vLabel,
				  const std::string& type1,
				  const std::string& type2 ){
  std::vector< TH1* > vDphi;
  std::vector< TF1* > vFits;

  std::string mcType = !type1.empty() ? "_" + type1 : "" ;
  
  // ---- loop over group  ----
  // ---- (jzn or trigger) ----
  vDphiWidths.resize( vLabel.size() );
  vDphiNent  .resize( vLabel.size() );
  for( uint iG = 0; iG < vLabel.size(); iG++ ){      
    THnSparse* hn     = vhn[iG];
    THnSparse* hnNent = vhnNent[iG];

    // truncate dphi bins with not many contents.
    // anaTool->TruncateHistoBins( hn, hnNent );

    std::string label = vLabel[iG];

    // in data only draw for all
    // if( m_isData && label.compare("All") ){ continue; } 

    TAxis* axis0 = hn->GetAxis( m_dPP->GetAxisI(0) );
    TAxis* axis1 = hn->GetAxis( m_dPP->GetAxisI(1) );
    TAxis* axis2 = hn->GetAxis( m_dPP->GetAxisI(2) );
    TAxis* axis3 = hn->GetAxis( m_dPP->GetAxisI(3) );

    TAxis* axis0Nent = hnNent->GetAxis( m_dPP->GetAxisI(0) );
    TAxis* axis1Nent = hnNent->GetAxis( m_dPP->GetAxisI(1) );
    TAxis* axis2Nent = hnNent->GetAxis( m_dPP->GetAxisI(2) );
    TAxis* axis3Nent = hnNent->GetAxis( m_dPP->GetAxisI(3) );
    
    int nAxis0Bins = axis0->GetNbins();
    int nAxis1Bins = axis1->GetNbins();
    int nAxis2Bins = axis2->GetNbins();
    int nAxis3Bins = axis3->GetNbins();

    int fAxisI     = m_dPP->GetAxisI(3);
 
    // ---- loop over ystars ----
    // resize collections
    vDphiWidths[iG].resize( nAxis0Bins );
    vDphiNent  [iG].resize( nAxis0Bins );

    for( int axis0Bin = 1; axis0Bin <= nAxis0Bins; axis0Bin++ ){
      // set ranges
      axis0    ->SetRange( axis0Bin, axis0Bin );
      axis0Nent->SetRange( axis0Bin, axis0Bin );

      double axis0Low, axis0Up;
      anaTool->GetBinRange
	( axis0, axis0Bin, axis0Bin, axis0Low, axis0Up );

      // resize collections
      vDphiWidths[iG][ axis0Bin - 1 ].resize( nAxis1Bins );
      vDphiNent  [iG][ axis0Bin - 1 ].resize( nAxis1Bins );

      for( int axis1Bin = 1; axis1Bin <= nAxis1Bins; axis1Bin++ ){
	// set ranges
	axis1    ->SetRange( axis1Bin, axis1Bin );
	axis1Nent->SetRange( axis1Bin, axis1Bin ); 
	
	double axis1Low, axis1Up;
	anaTool->GetBinRange
	  ( axis1, axis1Bin, axis1Bin, axis1Low, axis1Up );
	
	std::vector< TH1* > vDphiWidthsTemp;
	TCanvas cWidths("cWidths","cWidths",800,600);
	TLegend leg(0.68, 0.64, 0.99, 0.77);
	int style = 0;
	styleTool->SetLegendStyle( &leg );

	// general name for final canvas/histo. This is used in
	// naming other histograms also.
	std::string hTag =
	  Form( "dPhi%s_%s_%s",
		mcType.c_str(),
		anaTool->GetName( axis0Low, axis0Up, m_dPP->GetAxisName(0) ).c_str(),
		anaTool->GetName( axis1Low, axis1Up, m_dPP->GetAxisName(1) ).c_str() ); 
	
	// ---- loop over axis2 ----
	for( int axis2Bin = 1; axis2Bin <= nAxis2Bins; axis2Bin++ ){
	  // set ranges
	  axis2    ->SetRange( axis2Bin, axis2Bin );
	  axis2Nent->SetRange( axis2Bin, axis2Bin );
	  
	  double axis2Low , axis2Up;
	  anaTool->GetBinRange
	    ( axis2, axis2Bin, axis2Bin, axis2Low, axis2Up );
	  
	  TH1* hDphiWidths = hn->Projection( fAxisI );
	  hDphiWidths->Reset();
	  
	  hDphiWidths->SetName
	    ( Form( "h_%s_%s_%s", hTag.c_str(),
		    anaTool->GetName( axis2Low , axis2Up, m_dPP->GetAxisName(2) ).c_str(),
		    label.c_str() ) );
	  hDphiWidths->GetYaxis()->SetTitle( "|#Delta#phi| width" );
	  vDphiWidths[iG][ axis0Bin - 1 ][ axis1Bin - 1 ].push_back( hDphiWidths );
	  hDphiWidths->SetMarkerSize( hDphiWidths->GetMarkerSize() * 1.5 );
	  vDphiWidthsTemp.push_back( hDphiWidths );
	  styleTool->SetHStyle( hDphiWidths, style++ );
	  
	  TH1* hNent = hnNent->Projection( fAxisI );
	  hNent->SetName
	    ( Form( "h_dPhiNent%s_%s_%s_%s",
		    mcType.c_str(), hTag.c_str(),
		    anaTool->GetName( axis2Low, axis2Up, m_dPP->GetAxisName(2) ).c_str(),
		    label.c_str() ) );
	  vDphiNent[iG][ axis0Bin - 1 ][ axis1Bin - 1 ].push_back( hNent );
	  
	  leg.AddEntry
	    ( hDphiWidths, anaTool->GetLabel( axis2Low, axis2Up, m_dPP->GetAxisLabel(2) ) .c_str() );

	  // ---- loop over axis3 ----
	  for( int axis3Bin = 1; axis3Bin <= nAxis3Bins; axis3Bin++ ){
	    axis3->SetRange( axis3Bin, axis3Bin );
	    double axis3Low , axis3Up;
	    anaTool->GetBinRange
	      ( axis3, axis3Bin, axis3Bin, axis3Low, axis3Up );

	    // Take projection onto the dPhi axis
	    TH1* hDphi = hn->Projection( 4 );
	    hDphi->SetName
	      ( Form( "h_%s_%s_%s_%s", hTag.c_str(),
		      anaTool->GetName( axis2Low , axis2Up, m_dPP->GetAxisName(2) ).c_str(),
		      anaTool->GetName( axis3Low , axis3Up, m_dPP->GetAxisName(3) ).c_str(),
		      label.c_str() ) );
	    styleTool->SetHStyle( hDphi, 0 );
	    vDphi.push_back( hDphi );

	    TCanvas c( "c", hDphi->GetName(), 800, 600 );

	    hDphi->Draw();
	    if( hDphi->GetEntries() )
	      { hDphi->Scale( 1./hDphi->Integral() ); }
	    hDphi->GetYaxis()->SetTitle("Normalized Count");
	    hDphi->SetTitle("");
	    
	    drawTool->DrawTopLeftLabels
	      ( m_dPP, axis0Low, axis0Up, axis1Low, axis1Up,
		axis2Low, axis2Up, axis3Low, axis3Up, 0.8 );
	    
	    if( m_isData ){
	      drawTool->DrawAtlasInternalDataRight( 0, 0, m_is_pPb );
	    } else {
	      drawTool->DrawRightLatex( 0.88, 0.82, type1 );
	      drawTool->DrawAtlasInternalMCRight( 0, 0, type2 );
	    }

	    // now fit
	    TF1* fit = anaTool->FitDphi( hDphi );
	    styleTool->SetHStyle( fit, 0 );
	    vFits.push_back( fit );

	    // save the confidence intervals
	    // this is fit results + errors as histo.
	    // TH1* hDphiCI = static_cast<TH1D*>
	    //   ( hDphi->Clone( Form("%s_CI", hDphi->GetName() ) ) );
	    
	    TH1* hDphiCI = 
	      new TH1D( Form("%s_CI", hDphi->GetName() ), "", m_nDphiDphiBins, m_dPhiDphiMin, m_dPhiDphiMax );

	    if( !m_isData ){
	      (TVirtualFitter::GetFitter())->GetConfidenceIntervals(hDphiCI);
	      hDphiCI->Write();
	    }
	    
	    fit->Draw("same");

	    SaveAsROOT( c, hDphi->GetName() );

	    if( fit->GetParameter(1) < 0 )
	      { continue; }

	    // dont add entries for bins that have small counts
	    // aka truncate
	    if( hNent->GetBinContent( axis3Bin ) < 20 ){ continue; }
	    // Now, put results on histogram
	    hDphiWidths->SetBinContent( axis3Bin, fit->GetParameter(1) );
	    hDphiWidths->SetBinError  ( axis3Bin, fit->GetParError (1) );	    
	  } // end loop over axis3

	  // set range back to whole range otherwise
	  // histograms drawn for only one bin
	  axis3    ->SetRange( 1, -1 );
	  axis3Nent->SetRange( 1, -1 );

	} // end loop over axis2
	
	cWidths.cd();
	for( auto& h : vDphiWidthsTemp ){
	  h->SetMinimum( m_dPhiWidthMin );
	  h->SetMaximum( m_dPhiWidthMax );
	  h->GetYaxis()->SetNdivisions( 505 );
	  h->SetTitle("");
	  h->Draw("epsame");
	}
	leg.Draw("same");

	drawTool->DrawLeftLatex
	  ( 0.13, 0.87,anaTool->GetLabel
	    ( axis0Low, axis0Up, m_dPP->GetAxisLabel(0) ) );
	drawTool->DrawLeftLatex
	  ( 0.13, 0.82,anaTool->GetLabel
	    ( axis1Low, axis1Up, m_dPP->GetAxisLabel(1) ) );

	if( m_isData ){
	  drawTool->DrawAtlasInternalDataRight( 0, 0, m_is_pPb );
	}
	else {
	  drawTool->DrawRightLatex    ( 0.88, 0.82, type1 );
	  drawTool->DrawAtlasInternalMCRight( 0, 0, type2 );
	}

	
	// SaveAsAll( cWidths, Form("h_%s_%s", hTag.c_str(), label.c_str() ) );
      } // end loop over axis1     
    } // end loop over axis0
  } // end loop over iG
}

void DiJetAnalysis::PlotDphiTogether(){}

//---------------------------
//        Drawing
//---------------------------
void DiJetAnalysis::DrawCanvas( std::vector< TH1* >& vHIN,
				const std::string& type1,
				const std::string& type2,
				int spacing ){
  TCanvas c("c","c",800,600);
  
  TLegend leg(0.64, 0.61, 0.99, 0.82);
  styleTool->SetLegendStyle( &leg, 0.8 );
  leg.SetFillStyle(0);

  int style = 0;

  // for situations where dont want to
  // plot every single bin 
  // plot every n on canvas
  int dX = vHIN.size()/spacing; // plot every n
  for( uint xRange = 2; xRange < vHIN.size(); xRange += dX){
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
  
  TLegend leg(0.64, 0.63, 0.99, 0.8);
  styleTool->SetLegendStyle( &leg, 0.7 );
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
