#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <boost/assign.hpp>

#include "MyRoot.h"

const int pPbLumi2016 = 437; // ub
const int ppLumi2015  = 26;  // pb

//===================================
//         COMMON FUNCTIONS
//===================================
bool CT::AnalysisTools::EpsilonEqual( double a, double b )
{ return fabs( a - b ) < constants::EPSILON; }

// returns dphi in range 0<dphi<2pi
double CT::AnalysisTools::DPhiFC( double phi1, double phi2 ){
  double deltaPhi = phi1 - phi2;
  while( deltaPhi < 0 ) deltaPhi += 2*constants::PI;
  return deltaPhi;
}

double CT::AnalysisTools::DeltaPhi
( double phi1, double phi2 ){
  double deltaPhi = TMath::Abs(phi1 - phi2);
  if( deltaPhi > constants::PI ){ deltaPhi = 2*constants::PI - deltaPhi; };
  return deltaPhi;
}

double CT::AnalysisTools::DeltaR
(  const TLorentzVector& jet1, const TLorentzVector& jet2 ){  
  double deltaEta = jet1.Eta() - jet2.Eta();
  double deltaPhi = TMath::Abs( jet1.Phi() - jet2.Phi() );
  if(deltaPhi > TMath::Pi())
    deltaPhi = 2*TMath::Pi() - deltaPhi;
  return TMath::Sqrt( deltaEta*deltaEta + deltaPhi*deltaPhi );
}

bool CT::AnalysisTools::sortByDecendingPt
( const TLorentzVector& jet1, const TLorentzVector& jet2 ){
  return ( jet1.Pt() > jet2.Pt() );
}


bool CT::AnalysisTools::TruncateHistoBins( THnSparse* hn,
					   THnSparse* hnNent){
  TAxis* aAxis = hn->GetAxis(0); TAxis* cAxis = hn->GetAxis(2);
  TAxis* bAxis = hn->GetAxis(1); TAxis* dAxis = hn->GetAxis(3); 

  int nAbins = aAxis->GetNbins();
  int nBbins = bAxis->GetNbins();
  int nCbins = cAxis->GetNbins();
  int nDbins = dAxis->GetNbins();

  std::vector< int > x;
  x.resize( hn->GetNdimensions() );
  
  for( int aBin = 1; aBin <= nAbins; aBin++ ){
    for( int bBin = 1; bBin <= nBbins; bBin++ ){
      for( int cBin = 1; cBin <= nCbins; cBin++ ){
	for( int dBin = 1; dBin <= nDbins; dBin++ ){
	  x[0] = aBin; x[1] = bBin; x[2] = cBin; x[3] = dBin;
	  if( hnNent->GetBinContent( &x[0] )  < 5. )
	    { hn->SetBinContent( &x[0], 0. ); }
	}
      }
    }
  }
  
  return true;
}


bool CT::AnalysisTools::TruncateHistoBins( TH3* h3 ){
  for( int z = 1; z <= h3->GetZaxis()->GetNbins(); z++ ){
    for( int y = 1; y <= h3->GetYaxis()->GetNbins(); y++ ){    
      for( int x = 1; x <= h3->GetXaxis()->GetNbins(); x++ ){
	if( h3->GetBinContent( x, y, z ) < 3 )
	  { h3->SetBinContent( x, y, z, 0 ); }
      }
    }
  }
  return true;
}

bool CT::AnalysisTools::DoPrint( int ev ) {
  int statSize=1;
  if( ev != 0){
    double power=std::floor(log10(ev));
    statSize=(int)std::pow(10.,power);
  }
  if( ev%statSize == 0 ) return true;
  return false;
}

std::vector<std::string> CT::AnalysisTools::vectorise
(TString str, TString sep) {
  std::vector<std::string> result;
  TObjArray *strings = str.Tokenize(sep.Data());
  if (strings->GetEntries()==0) { delete strings; return result; }
  TIter istr(strings);
  while (TObjString* os=(TObjString*)istr())
    {result.push_back(std::string(os->GetString()));}
  delete strings; return result;
}   

std::vector<double> CT::AnalysisTools::vectoriseD
(TString str, TString sep) {
  std::vector<double> result;
  std::vector<std::string> vecS = vectorise(str,sep);
  for (uint i=0;i<vecS.size();++i)
    {result.push_back(std::stod(vecS[i]));}
  return result;
}   


TF1* CT::AnalysisTools::FitDphi( TH1* hProj, double xLow, double xHigh ){
  auto expoFunction = [&]( double* x, double* par){
    return par[0]*std::exp((x[0]-constants::PI)/par[1])+par[2];
  };
  
  TF1* fit = new TF1( Form("f_%s", hProj->GetName()),
		      expoFunction, 0, constants::PI, 3);
  
  if( !hProj->GetEntries() )
    { return fit; }

  fit->SetParameters( 0.3, 0.3, 0 );
  
  hProj->Fit( fit->GetName(), "NQR", "" , 0, constants::PI);

  return fit;
}

TF1* CT::AnalysisTools::FitGaussian( TH1* hProj, double xLow, double xHigh ){

  TF1* fit  = new TF1( Form("f_%s", hProj->GetName()), "gaus(0)" );

  double hXmin = hProj->GetXaxis()->GetXmin();
  double hXmax = hProj->GetXaxis()->GetXmax();
  
  // fit once 
  double mean   = hProj->GetMean();
  double rms    = hProj->GetRMS();

  double fitMin = mean - 2.0 * rms;
  double fitMax = mean + 2.0 * rms;

  if( hProj->GetEntries() < 5 || fitMin < hXmin || fitMax > hXmax )
    { return fit; }
  
  hProj->Fit( fit->GetName(), "NQR", "", fitMin, fitMax );

  // fit second time with better parameters
  mean   = fit->GetParameter(1);
  rms    = fit->GetParameter(2);

  fitMin = mean - 2.0 * rms;
  fitMax = mean + 2.0 * rms;

  if( hProj->GetEntries() < 5 || fitMin < hXmin || fitMax > hXmax )
    { return fit; }
  
  fit->SetRange( fitMin, fitMax );
  
  hProj->Fit( fit->GetName(), "NQR", "", fitMin, fitMax );

  return fit;
}

bool CT::AnalysisTools::GetRatioFromFits( TH1* hR, TF1* f1, TF1* f2 ){
  /*
  // histos need to have same binning...
  for( int xBin = 1; xBin <= h1->GetNbinsX(); xBin++ ){
    double binCenter = hR->GetBinCenter( xBin );
    double ratio = f1->Eval( binCenter )/f2->Eval( binCenter );
    hR->SetBinContent( xBin, ratio );
  }
  */
  return true;
}

void CT::AnalysisTools::GetBinRange( TAxis* a,
				 int b1, int b2,
				 double& x1, double& x2){
  x1 = a->GetBinLowEdge( b1 );
  x2 = a->GetBinUpEdge ( b2 );
}

std::string CT::AnalysisTools::GetName( double v1, double v2,
					const std::string& var){
  // for now, edit this later
  int multiplier = std::abs(v1) < 10 ? 10 : 1;
  
  std::stringstream ss;
  if( v1 != v2 ){
    ss <<  boost::format("%d_%s_%d")
      % static_cast<int>(std::abs(multiplier * v1))
      % var
      % static_cast<int>(std::abs(multiplier * v2)) ;
  } else {
    ss <<  boost::format("%d_%s")
      % static_cast<int>(std::abs(multiplier * v1))
      % var ;
  }
  return ss.str();
}

std::string CT::AnalysisTools::GetEtaLabel( double etaMin,
					    double etaMax,
					    bool is_pPb){
  std::stringstream ss;
  if( is_pPb )
    { ss << boost::format("%3.1f<#eta<%3.1f") % etaMin % etaMax;}
  else {
    if( std::abs(etaMin) > std::abs(etaMax) )
      { ss << boost::format("%3.1f<|#eta|<%3.1f")
	  % std::abs(etaMax) % std::abs(etaMin); }
    else
      { ss << boost::format("%3.1f<|#eta|<%3.1f")
	  % std::abs(etaMin) % std::abs(etaMax); }
  }
  return ss.str();
}

std::string CT::AnalysisTools::GetYstarLabel( double ystarMin,
					      double ystarMax,
					      bool is_pPb,
					      std::string label ){
  std::stringstream ss;
  if( is_pPb )
    { ss << boost::format("%3.1f<%s<%3.1f") % ystarMin % label % ystarMax;}
  else {
    if( std::abs(ystarMin) > std::abs(ystarMax) )
      { ss << boost::format("%3.1f<|%s|<%3.1f")
	  % std::abs(ystarMax) % label % std::abs(ystarMin); }
    else
      { ss << boost::format("%3.1f<|%s|<%3.1f")
	  % std::abs(ystarMin) % label % std::abs(ystarMax); }
  }
  return ss.str();
}

std::string CT::AnalysisTools::GetLabel
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

// should be in miscallaneous
void CT::AnalysisTools::CheckWriteDir( const char* c_dirOut ){
  boost::filesystem::path dir( c_dirOut );  
  if(!(boost::filesystem::exists(dir))){
    std::cout<< c_dirOut << " doesn't Exist."<<std::endl;
    if (boost::filesystem::create_directory(dir))
      std::cout << "....Successfully Created !" << std::endl;
  }
}

//===================================
//          STYLE STUFF
//===================================
// ======= Styles for Stuff ======

void CT::StyleTools::SetCustomMarkerStyle( TH1* his , int iflag ){	
  //Set Color
  his->SetLineWidth(2);
  if( iflag == 0 ){
    his->SetLineColor(kBlack);
    his->SetMarkerColor(kBlack);
    his->SetMarkerStyle(20);
    his->SetMarkerSize(1.2);
  } 
  else if(iflag == 1 ){
    his->SetLineColor(kRed);
    his->SetMarkerColor(kRed);
    his->SetMarkerStyle(21);
    his->SetMarkerSize(1.1);
  }
  else if(iflag == 2 ){
    his->SetLineColor(kAzure-3);
    his->SetMarkerColor(kAzure-3);
    his->SetMarkerStyle(33);
    his->SetMarkerSize(1.8);
  }
  else if(iflag == 3 ){
    his->SetLineColor(kSpring-6);
    his->SetMarkerColor(kSpring-6);
    his->SetMarkerStyle(34);
    his->SetMarkerSize(1.5);
  }
  else if(iflag == 4 ){
    his->SetLineColor(kOrange+1);
    his->SetLineWidth(2);
    his->SetMarkerColor(kOrange+1);
    his->SetMarkerStyle(29);
    his->SetMarkerSize(1.6);
  }
  else if(iflag == 5 ){
    his->SetLineColor(kBlack);
    his->SetMarkerColor(kBlack);
    his->SetMarkerStyle(24);
    his->SetMarkerSize(1.3);
  }
  else if(iflag == 6 ){
    his->SetLineColor(kRed);
    his->SetMarkerColor(kRed);
    his->SetMarkerStyle(25);
    his->SetMarkerSize(1.2);
  }  
  else if(iflag == 7 ){
    his->SetLineColor(kAzure-3);
    his->SetMarkerColor(kAzure-3);
    his->SetMarkerStyle(27);
    his->SetMarkerSize(1.9);
  }
  else if(iflag == 8 ){
    his->SetLineColor(kSpring-6);
    his->SetMarkerColor(kSpring-6);
    his->SetMarkerStyle(28);
    his->SetMarkerSize(1.6);
  }
  else if(iflag == 9 ){
    his->SetLineColor(kOrange+1);
    his->SetLineWidth(2);
    his->SetMarkerColor(kOrange+1);
    his->SetMarkerStyle(30);
    his->SetMarkerSize(1.7);
  }

}

void CT::StyleTools::SetCustomMarkerStyle( TGraph* graph , int iflag ){
  //Set Color
  graph->SetFillColor(0);
  graph->SetLineWidth(2);
  if( iflag == 0 ){
    graph->SetLineColor(kBlack);
    graph->SetMarkerColor(kBlack);
    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(1.2);
  } 
  else if(iflag == 1 ){
    graph->SetLineColor(kRed);
    graph->SetMarkerColor(kRed);
    graph->SetMarkerStyle(21);
    graph->SetMarkerSize(1.1);
  }
  else if(iflag == 2 ){
    graph->SetLineColor(kAzure-3);
    graph->SetMarkerColor(kAzure-3);
    graph->SetMarkerStyle(33);
    graph->SetMarkerSize(1.8);
  }
  else if(iflag == 3 ){
    graph->SetLineColor(kSpring-6);
    graph->SetMarkerColor(kSpring-6);
    graph->SetMarkerStyle(34);
    graph->SetMarkerSize(1.5);
  }
  else if(iflag == 4 ){
    graph->SetLineColor(kBlack);
    graph->SetMarkerColor(kBlack);
    graph->SetMarkerStyle(24);
    graph->SetMarkerSize(1.1);
  }
  else if(iflag == 5 ){
    graph->SetLineColor(kRed);
    graph->SetMarkerColor(kRed);
    graph->SetMarkerStyle(25);
    graph->SetMarkerSize(1.4);
  }  
  else if(iflag == 6 ){
    graph->SetLineColor(kAzure-3);
    graph->SetMarkerColor(kAzure-3);
    graph->SetMarkerStyle(27);
    graph->SetMarkerSize(1.8);
  }
  else if(iflag == 7 ){
    graph->SetLineColor(kSpring-6);
    graph->SetMarkerColor(kSpring-6);
    graph->SetMarkerStyle(28);
    graph->SetMarkerSize(1.5);
  }
}

void CT::StyleTools::SetHStyle( TH1* his, int iflag, double scale)
{
  his->SetLineWidth(2);
  his->SetStats(0);

  his->SetTitleFont( 43, "t");
  his->SetTitleSize( (int)(32 * scale), "t"); 

  his->SetTitleFont( 43, "xyz");
  his->SetTitleSize( (int)(32 * scale), "xyz");
  his->SetTitleOffset(1.1, "xyz");
  
  his->SetLabelFont( 43, "xyz");
  his->SetLabelSize( (int)(30 * scale), "xyz");
  
  SetCustomMarkerStyle( his, iflag );  
}

void CT::StyleTools::SetHStyle( TGraph* graph, int iflag, double scale)
{
  graph->SetLineWidth(2);

  graph->GetXaxis()->SetTitleFont( 43 );
  graph->GetXaxis()->SetTitleSize( (int)(32 * scale) );  
  graph->GetXaxis()->SetTitleOffset(1.1);

  graph->GetXaxis()->SetLabelFont( 43 );
  graph->GetXaxis()->SetLabelSize( (int)(30 * scale) );

  graph->GetYaxis()->SetTitleFont( 43 );
  graph->GetYaxis()->SetTitleSize( (int)(32 * scale) );
  graph->GetYaxis()->SetTitleOffset(1.1);

  graph->GetYaxis()->SetLabelFont( 43 );
  graph->GetYaxis()->SetLabelSize( (int)(30 * scale) );
  
  SetCustomMarkerStyle( graph, iflag );
}


void CT::StyleTools::SetHStyle( TF1* funct, int iflag, double scale)
{
  funct->SetLineWidth(2);
  
  funct->GetXaxis()->SetTitleFont( 43 );
  funct->GetXaxis()->SetTitleSize( (int)(32 * scale) );  
  funct->GetXaxis()->SetTitleOffset(1.1);

  funct->GetXaxis()->SetLabelFont( 43 );
  funct->GetXaxis()->SetLabelSize( (int)(30 * scale) );

  funct->GetYaxis()->SetTitleFont( 43 );
  funct->GetYaxis()->SetTitleSize( (int)(32 * scale) );
  funct->GetYaxis()->SetTitleOffset(1.1);

  funct->GetYaxis()->SetLabelFont( 43 );
  funct->GetYaxis()->SetLabelSize( (int)(30 * scale) );
  
  // SetCustomMarkerStyle( funct, iflag );
}

TH1F* CT::StyleTools::SetCStyleEff ( TCanvas& c,
				     double x0, double y0,
				     double x1, double y1,
				     const std::string& title ){
  c.DrawFrame( x0, y0, x1, y1, title.c_str() );
  TH1F* hF = c.DrawFrame( x0, y0, x1, y1, title.c_str() );
  hF->GetXaxis()->SetNdivisions(505);  
  hF->GetYaxis()->SetNdivisions(510);  
  StyleTools::SetHStyle( hF, 0, StyleTools::hSS );
  return hF;
}

void CT::StyleTools::SetLegendStyle(TLegend * legend, double scale)
{
  legend->SetBorderSize(0);
  legend->SetTextFont(43);
  legend->SetTextSize(28 * scale);
  legend->SetLineColor(1);
  legend->SetLineStyle(1);
  legend->SetLineWidth(1);
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
}

const double CT::StyleTools::lSS = 0.60;

const double CT::StyleTools::hSS = 0.80;

//===================================
//          DRAWING STUFF
//===================================

// ============ GENERAL ================

void CT::DrawTools::DrawRightLatex
( double x, double y , const std::string& s, double scale, int color ){
  TLatex tltx; 
  tltx.SetTextFont(43);
  tltx.SetTextSize((int)(32*scale));
  tltx.SetTextColor(color);
  tltx.SetNDC(kTRUE);
  tltx.SetTextAlign(32);
  tltx.DrawLatex( x, y, s.c_str() );
}

void CT::DrawTools::DrawLeftLatex
( double x, double y , const std::string& s, double scale, int color  ){
  TLatex tltx; 
  tltx.SetTextFont(43);
  tltx.SetTextSize((int)(32*scale));
  tltx.SetTextColor(color);
  tltx.SetNDC(kTRUE);
  tltx.SetTextAlign(12);
  tltx.DrawLatex( x, y, s.c_str() );
}

void CT::DrawTools::DrawCenterLatex
( double x, double y , const std::string& s, double scale, int color ){
  TLatex tltx; 
  tltx.SetTextFont(43);
  tltx.SetTextSize((int)(32*scale));
  tltx.SetTextColor(color);
  tltx.SetNDC(kTRUE);
  tltx.SetTextAlign(22);
  tltx.DrawLatex( x, y, s.c_str() );
}

// ============ DATA ================

void CT::DrawTools::DrawAtlasInternal( double scale ){
  DrawRightLatex
    (0.88 , 0.93, "#bf{#font[72]{ATLAS}} Internal", scale, 1 );
}

void CT::DrawTools::DrawAtlasInternalDataRight
( double x0, double y0, bool is_pPb, double scale ){
  DrawRightLatex
    (0.88 , 0.93,"#bf{#font[72]{ATLAS}} Internal", scale, 1 );
  if( is_pPb ){
    DrawRightLatex
      (0.88, 0.87 + y0, Form("#it{p}+Pb 2016, %i #mub^{-1}",
			     pPbLumi2016), scale, 1 );
  } else {
    DrawRightLatex
      (0.88 + x0, 0.87 + y0,Form("#it{pp} 2015, %i pb^{-1}",
			ppLumi2015), scale, 1 );
  }
  DrawRightLatex(0.88 + x0, 0.81 + y0, 
		 "#sqrt{s_{NN}}=5.02 TeV", scale, 1 );
}

void CT::DrawTools::DrawAtlasInternalDataLeft
( double x0, double y0, bool is_pPb, double scale ){
  DrawRightLatex(0.875, 0.93, 
		 "#bf{#font[72]{ATLAS}} Internal", scale, 1 );
  if( is_pPb ){
    DrawLeftLatex(0.18 + x0, 0.87 + y0, 
		   Form("#it{p}+Pb 2016, %i #mub^{-1}",
			pPbLumi2016), scale, 1 );
  } else {
    DrawLeftLatex(0.18 + x0, 0.87 + y0, 
		   Form("#it{pp} 2015, %i pb^{-1}",
			ppLumi2015), scale, 1 );
  }
  DrawLeftLatex(0.18 + x0, 0.81 + y0, 
		"#sqrt{s_{NN}}=5.02 TeV", scale, 1 ); 
}

// ============ MC ================

void CT::DrawTools::DrawAtlasInternalMCRight
( double x0, double y0, const std::string& mcType, double scale ){ 
  DrawRightLatex(0.88, 0.93, 
		 "#bf{#font[72]{ATLAS}} Simulation Internal", scale, 1 );
  DrawRightLatex(0.88 + x0, 0.87, mcType, scale, 1 );
 }


void CT::DrawTools::DrawAtlasInternalMCLeft
( double x0, double y0, const std::string& mcType, double scale ){ 
  DrawRightLatex(0.88, 0.93, 
		 "#bf{#font[72]{ATLAS}} Simulation Internal", scale, 1 );

  DrawLeftLatex(0.18 + x0, 0.87, mcType, scale, 1 );
}
