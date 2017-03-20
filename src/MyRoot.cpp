#include "MyRoot.h"

const int pPbLumi2016 = 437; // ub
const int ppLumi2015  = 26;  // pb

//===================================
//         COMMON FUNCTIONS
//===================================
bool AnalysisTools::isForward( const double& eta ){
  return ( TMath::Abs(eta) < constants::FETAMAX && 
	   TMath::Abs(eta) > constants::FETAMIN );
}

bool AnalysisTools::isCentral( const double& eta ){
  return ( TMath::Abs(eta) < constants::CETAMAX );
}

bool AnalysisTools::EpsilonEqual( double a, double b ){
  return fabs( a - b ) < constants::EPSILON;
}

// returns dphi in range 0<dphi<2pi
double AnalysisTools::DPhiFC( double phi1, double phi2 ){
  double deltaPhi = phi1 - phi2;

  while( deltaPhi < 0 ) deltaPhi += 2*constants::PI;

  return deltaPhi;
}

double AnalysisTools::DeltaPhi
( double phi1, double phi2 ){
  double deltaPhi = TMath::Abs(phi1 - phi2);

  if( deltaPhi > constants::PI ){ deltaPhi = 2*constants::PI - deltaPhi; };

  return deltaPhi;
}

double AnalysisTools::DeltaR
(  const TLorentzVector& jet1, const TLorentzVector& jet2 ){  
  double deltaEta = jet1.Eta() - jet2.Eta();
  double deltaPhi = TMath::Abs( jet1.Phi() - jet2.Phi() );
  if(deltaPhi > TMath::Pi())
    deltaPhi = 2*TMath::Pi() - deltaPhi;
  return TMath::Sqrt( deltaEta*deltaEta + deltaPhi*deltaPhi );
}

bool AnalysisTools::sortByDecendingPt
( const TLorentzVector& jet1, const TLorentzVector& jet2 ){
  return ( jet1.Pt() > jet2.Pt() );
}

bool AnalysisTools::TruncateHistoBins( TH3* h3 ){
  for( int z = 1; z <= h3->GetZaxis()->GetNbins(); z++ ){
    for( int y = 1; y <= h3->GetYaxis()->GetNbins(); y++ ){    
      for( int x = 1; x <= h3->GetXaxis()->GetNbins(); x++ ){
	if( h3->GetBinContent( x, y, z ) < 3 ) h3->SetBinContent( x, y, z, 0 );
      }
    }
  }
  return true;
}

bool AnalysisTools::DoPrint( int ev ) {
  int statSize=1;
  if( ev != 0){
    double power=std::floor(log10(ev));
    statSize=(int)std::pow(10.,power);
  }
  if( ev%statSize == 0 ) return true;
  return false;
}

std::vector<std::string> AnalysisTools::vectorise(TString str, TString sep) {
  std::vector<std::string> result; TObjArray *strings = str.Tokenize(sep.Data());
  if (strings->GetEntries()==0) { delete strings; return result; }
  TIter istr(strings);
  while (TObjString* os=(TObjString*)istr()) result.push_back(std::string(os->GetString()));
  delete strings; return result;
}   


void AnalysisTools::FitGaussian( TH1* hProj, TF1* fit ){
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

//===================================
//          DRAWING STUFF
//===================================

// ============ GENERAL ================

void DrawTools::DrawRightLatex
( double x, double y , const char* s, float scale = 1, int color = 1){
  TLatex tltx; 
  tltx.SetTextFont(43); tltx.SetTextSize((int)(32*scale)); tltx.SetTextColor(color);
  tltx.SetNDC(kTRUE);
  tltx.SetTextAlign(32);
  tltx.DrawLatex( x, y, s );
}

void DrawTools::DrawLeftLatex
( double x, double y , const char* s, float scale = 1, int color = 1){
  TLatex tltx; 
  tltx.SetTextFont(43); tltx.SetTextSize((int)(32*scale)); tltx.SetTextColor(color);
  tltx.SetNDC(kTRUE);
  tltx.SetTextAlign(12);
  tltx.DrawLatex( x, y, s );
}

void DrawTools::DrawCenterLatex
( double x, double y , const char* s, float scale = 1, int color = 1){
  TLatex tltx; 
  tltx.SetTextFont(43); tltx.SetTextSize((int)(32*scale)); tltx.SetTextColor(color);
  tltx.SetNDC(kTRUE);
  tltx.SetTextAlign(22);
  tltx.DrawLatex( x, y, s );
}

// ============ DATA ================

void DrawTools::DrawAtlasInternalDataRight
( double x0, double y0, double scale, bool is_pPb ){
  DrawRightLatex(0.88 , 0.93,
		 "#bf{#font[72]{ATLAS}} Internal", scale);
  if( is_pPb ){
    DrawRightLatex(0.88, 0.87 + y0, 
		   Form("#it{p}+Pb 2016, %i #mub^{-1}",
			pPbLumi2016), scale);
  } else {
    DrawRightLatex(0.88 + x0, 0.87 + y0, 
		   Form("#it{pp} 2015, %i pb^{-1}",
			ppLumi2015), scale);
  }
  DrawRightLatex(0.88 + x0, 0.81 + y0, 
		 "#sqrt{s_{NN}}=5.02 TeV", scale);
}

void DrawTools::DrawAtlasInternalDataLeft
( double x0, double y0, double scale, bool is_pPb ){
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

void DrawTools::DrawAtlasInternalMCRight
( double x0, double y0, double scale,
  const std::string& mcType ){ 
  DrawRightLatex(0.88, 0.93, 
		 "#bf{#font[72]{ATLAS}} Simulation Internal", scale);
  DrawRightLatex(0.88 + x0, 0.87,
		 mcType.c_str(), scale);
 }


void DrawTools::DrawAtlasInternalMCLeft
( double x0, double y0, double scale,
  const std::string& mcType ){ 
  DrawRightLatex(0.88, 0.93, 
		 "#bf{#font[72]{ATLAS}} Simulation Internal", scale);

  DrawLeftLatex(0.18 + x0, 0.87,
		mcType.c_str(), scale);
}

// ======= Styles for Stuff ======

void StyleTools::SetCustomMarkerStyle( TH1* his , int iflag ){	
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

void StyleTools::SetCustomMarkerStyle( TGraph* graph , int iflag ){
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

void StyleTools::SetHStyle( TH1* his, int iflag, float scale)
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

void StyleTools::SetHStyle( TGraph* graph, int iflag, float scale)
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


void StyleTools::SetHStyle( TF1* funct, int iflag, float scale)
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

void StyleTools::SetLegendStyle(TLegend * legend, float scale)
{
  legend->SetBorderSize(0);
  legend->SetTextFont(43);
  legend->SetTextSize(28 * scale);
  legend->SetLineColor(1);
  legend->SetLineStyle(1);
  legend->SetLineWidth(1);
  legend->SetFillColor(0);
  legend->SetFillStyle(1001);
}
