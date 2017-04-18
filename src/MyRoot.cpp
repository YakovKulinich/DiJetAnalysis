#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <boost/assign.hpp>

#include "MyRoot.h"

const int pPbLumi2016 = 437; // ub
const int ppLumi2015  = 26;  // pb

//===================================
//         COMMON FUNCTIONS
//===================================
bool CT::AnalysisTools::IsForward( const double& eta )
{ return ( TMath::Abs(eta) < constants::FETAMAX && 
	   TMath::Abs(eta) > constants::FETAMIN ); }

bool CT::AnalysisTools::IsCentral( const double& eta )
{ return ( TMath::Abs(eta) < constants::CETAMAX ); }

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

bool CT::AnalysisTools::TruncateHistoBins( TH3* h3 ){
  for( int z = 1; z <= h3->GetZaxis()->GetNbins(); z++ ){
    for( int y = 1; y <= h3->GetYaxis()->GetNbins(); y++ ){    
      for( int x = 1; x <= h3->GetXaxis()->GetNbins(); x++ ){
	if( h3->GetBinContent( x, y, z ) < 3 ) h3->SetBinContent( x, y, z, 0 );
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


void CT::AnalysisTools::FitGaussian( TH1* hProj, TF1* fit ){
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

void CT::AnalysisTools::GetBinRange( TAxis* a,
				 int b1, int b2,
				 double& x1, double& x2){
  x1 = a->GetBinLowEdge( b1 );
  x2 = a->GetBinUpEdge ( b2 );
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
    his->SetLineColor(kViolet);
    his->SetMarkerColor(kViolet);
    his->SetMarkerStyle(22);
    his->SetMarkerSize(1.6);
  }
  else if(iflag == 6 ){
    his->SetLineColor(kCyan+1);
    his->SetMarkerColor(kCyan+1);
    his->SetMarkerStyle(23);
    his->SetMarkerSize(1.6);
  }
  else if(iflag == 7 ){
    his->SetLineColor(46);
    his->SetMarkerColor(46);
    his->SetMarkerStyle(21);
    his->SetMarkerSize(1.4);
  }
  else if(iflag == 8 ){
    his->SetLineColor(kRed+3);
    his->SetMarkerColor(kRed+3);
    his->SetMarkerStyle(25);
    his->SetMarkerSize(1.4);
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
    graph->SetLineColor(kOrange+1);
    graph->SetMarkerColor(kOrange+1);
    graph->SetMarkerStyle(29);
    graph->SetMarkerSize(1.6);
  }
  else if(iflag == 5 ){
    graph->SetLineColor(kViolet);
    graph->SetMarkerColor(kViolet);
    graph->SetMarkerStyle(22);
    graph->SetMarkerSize(1.6);
  }
  else if(iflag == 6 ){
    graph->SetLineColor(kCyan+1);
    graph->SetMarkerColor(kCyan+1);
    graph->SetMarkerStyle(23);
    graph->SetMarkerSize(1.6);
  }
  else if(iflag == 7 ){
    graph->SetLineColor(46);
    graph->SetMarkerColor(46);
    graph->SetMarkerStyle(21);
    graph->SetMarkerSize(1.4);
  }
  else if(iflag == 8 ){
    graph->SetLineColor(kRed+3);
    graph->SetMarkerColor(kRed+3);
    graph->SetMarkerStyle(25);
    graph->SetMarkerSize(1.4);
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
  his->SetTitleOffset(1.2, "xyz");
  
  his->SetLabelFont( 43, "xyz");
  his->SetLabelSize( (int)(30 * scale), "xyz");
  
  SetCustomMarkerStyle( his, iflag );  
}

void CT::StyleTools::SetHStyle( TGraph* graph, int iflag, double scale)
{
  graph->SetLineWidth(2);

  graph->GetXaxis()->SetTitleFont( 43 );
  graph->GetXaxis()->SetTitleSize( (int)(32 * scale) );  
  graph->GetXaxis()->SetTitleOffset(1.2);

  graph->GetXaxis()->SetLabelFont( 43 );
  graph->GetXaxis()->SetLabelSize( (int)(30 * scale) );

  graph->GetYaxis()->SetTitleFont( 43 );
  graph->GetYaxis()->SetTitleSize( (int)(32 * scale) );
  graph->GetYaxis()->SetTitleOffset(1.2);

  graph->GetYaxis()->SetLabelFont( 43 );
  graph->GetYaxis()->SetLabelSize( (int)(30 * scale) );
  
  SetCustomMarkerStyle( graph, iflag );
}


void CT::StyleTools::SetHStyle( TF1* funct, int iflag, double scale)
{
  funct->SetLineWidth(2);

  funct->GetXaxis()->SetTitleFont( 43 );
  funct->GetXaxis()->SetTitleSize( (int)(32 * scale) );  
  funct->GetXaxis()->SetTitleOffset(1.2);

  funct->GetXaxis()->SetLabelFont( 43 );
  funct->GetXaxis()->SetLabelSize( (int)(30 * scale) );

  funct->GetYaxis()->SetTitleFont( 43 );
  funct->GetYaxis()->SetTitleSize( (int)(32 * scale) );
  funct->GetYaxis()->SetTitleOffset(1.2);

  funct->GetYaxis()->SetLabelFont( 43 );
  funct->GetYaxis()->SetLabelSize( (int)(30 * scale) );
  
  // SetCustomMarkerStyle( funct, iflag );
}

TH1F* CT::StyleTools::SetCStyleEff
( TCanvas& c,
  double x0, double y0, double x1, double y1,
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

const double CT::StyleTools::hSS = 0.75;

//===================================
//          DRAWING STUFF
//===================================

// ============ GENERAL ================

void CT::DrawTools::DrawRightLatex
( double x, double y , const char* s, int color, double scale){
  TLatex tltx; 
  tltx.SetTextFont(43);
  tltx.SetTextSize((int)(32*scale));
  tltx.SetTextColor(color);
  tltx.SetNDC(kTRUE);
  tltx.SetTextAlign(32);
  tltx.DrawLatex( x, y, s );
}

void CT::DrawTools::DrawLeftLatex
( double x, double y , const char* s, int color, double scale ){
  TLatex tltx; 
  tltx.SetTextFont(43);
  tltx.SetTextSize((int)(32*scale));
  tltx.SetTextColor(color);
  tltx.SetNDC(kTRUE);
  tltx.SetTextAlign(12);
  tltx.DrawLatex( x, y, s );
}

void CT::DrawTools::DrawCenterLatex
( double x, double y , const char* s, int color, double scale){
  TLatex tltx; 
  tltx.SetTextFont(43);
  tltx.SetTextSize((int)(32*scale));
  tltx.SetTextColor(color);
  tltx.SetNDC(kTRUE);
  tltx.SetTextAlign(22);
  tltx.DrawLatex( x, y, s );
}

// ============ DATA ================

void CT::DrawTools::DrawAtlasInternal( double scale ){
  DrawRightLatex
    (0.88 , 0.93, "#bf{#font[72]{ATLAS}} Internal", scale);
}

void CT::DrawTools::DrawAtlasInternalDataRight
( double x0, double y0, bool is_pPb, double scale ){
  DrawRightLatex
    (0.88 , 0.93,"#bf{#font[72]{ATLAS}} Internal", scale);
  if( is_pPb ){
    DrawRightLatex
      (0.88, 0.87 + y0, Form("#it{p}+Pb 2016, %i #mub^{-1}",
			     pPbLumi2016), scale);
  } else {
    DrawRightLatex
      (0.88 + x0, 0.87 + y0,Form("#it{pp} 2015, %i pb^{-1}",
			ppLumi2015), scale);
  }
  DrawRightLatex(0.88 + x0, 0.81 + y0, 
		 "#sqrt{s_{NN}}=5.02 TeV", scale);
}

void CT::DrawTools::DrawAtlasInternalDataLeft
( double x0, double y0, bool is_pPb, double scale ){
  DrawRightLatex(0.875, 0.93, 
		 "#bf{#font[72]{ATLAS}} Internal", scale);
  if( is_pPb ){
    DrawLeftLatex(0.18 + x0, 0.87 + y0, 
		   Form("#it{p}+Pb 2016, %i #mub^{-1}",
			pPbLumi2016), scale);
  } else {
    DrawLeftLatex(0.18 + x0, 0.87 + y0, 
		   Form("#it{pp} 2015, %i pb^{-1}",
			ppLumi2015), scale);
  }
  DrawLeftLatex(0.18 + x0, 0.81 + y0, 
		"#sqrt{s_{NN}}=5.02 TeV", scale); 
}

// ============ MC ================

void CT::DrawTools::DrawAtlasInternalMCRight
( double x0, double y0, const std::string& mcType, double scale ){ 
  DrawRightLatex(0.88, 0.93, 
		 "#bf{#font[72]{ATLAS}} Simulation Internal", scale);
  DrawRightLatex(0.88 + x0, 0.87,
		 mcType.c_str(), scale);
 }


void CT::DrawTools::DrawAtlasInternalMCLeft
( double x0, double y0, const std::string& mcType, double scale ){ 
  DrawRightLatex(0.88, 0.93, 
		 "#bf{#font[72]{ATLAS}} Simulation Internal", scale);

  DrawLeftLatex(0.18 + x0, 0.87,
		mcType.c_str(), scale);
}
