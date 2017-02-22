#include "MyRoot.h"

const int pPbLumi2016 = 437; // ub
const int ppLumi2015  = 26;  // pb

//===================================
//         COMMON FUNCTIONS
//===================================

bool isForward( const double& eta ){
  return ( TMath::Abs(eta) < constants::FETAMAX && 
	   TMath::Abs(eta) > constants::FETAMIN );
}

bool isCentral( const double& eta ){
  return ( TMath::Abs(eta) < constants::CETAMAX );
}


bool isInTriggerEtaRange( const double& eta ){
  return ( TMath::Abs(eta) < constants::TETAMAX && 
	   TMath::Abs(eta) > constants::TETAMIN );
}


bool EpsilonEqual( double a, double b ){
  return fabs( a - b ) < constants::EPSILON;
}

// returns dphi in range 0<dphi<2pi
double DPhiFC( double phi1, double phi2 ){
  double deltaPhi = phi1 - phi2;

  while( deltaPhi < 0 ) deltaPhi += 2*constants::PI;

  return deltaPhi;
}

double DeltaPhi( double phi1, double phi2 ){
  double deltaPhi = TMath::Abs(phi1 - phi2);

  if( deltaPhi > constants::PI ){ deltaPhi = 2*constants::PI - deltaPhi; };

  return deltaPhi;
}

double DeltaR(  const TLorentzVector& jet1, const TLorentzVector& jet2 ){  
  double deltaEta = jet1.Eta() - jet2.Eta();
  double deltaPhi = TMath::Abs( jet1.Phi() - jet2.Phi() );
  if(deltaPhi > TMath::Pi())
    deltaPhi = 2*TMath::Pi() - deltaPhi;
  return TMath::Sqrt( deltaEta*deltaEta + deltaPhi*deltaPhi );
}

bool sortByDecendingPt( const TLorentzVector& jet1, const TLorentzVector& jet2 ){
  return ( jet1.Pt() > jet2.Pt() );
}

bool TruncateHistoBins( TH3* h3 ){
  for( int z = 1; z <= h3->GetZaxis()->GetNbins(); z++ ){
    for( int y = 1; y <= h3->GetYaxis()->GetNbins(); y++ ){    
      for( int x = 1; x <= h3->GetXaxis()->GetNbins(); x++ ){
	if( h3->GetBinContent( x, y, z ) < 3 ) h3->SetBinContent( x, y, z, 0 );
      }
    }
  }
  return true;
}

bool DoPrint( int ev ) {
  int statSize=1;
  if( ev != 0){
    double power=std::floor(log10(ev));
    statSize=(int)std::pow(10.,power);
  }
  if( ev%statSize == 0 ) return true;
  return false;
}

std::vector<std::string> vectorise(TString str, TString sep) {
  std::vector<std::string> result; TObjArray *strings = str.Tokenize(sep.Data());
  if (strings->GetEntries()==0) { delete strings; return result; }
  TIter istr(strings);
  while (TObjString* os=(TObjString*)istr()) result.push_back(std::string(os->GetString()));
  delete strings; return result;
}   

//===================================
//          DRAWING STUFF
//===================================


void DrawRightLatex( double x, double y , const char* s, float scale = 1, int color = 1){
  TLatex tltx; 
  tltx.SetTextFont(43); tltx.SetTextSize((int)(32*scale)); tltx.SetTextColor(color);
  tltx.SetNDC(kTRUE);
  tltx.SetTextAlign(32);
  tltx.DrawLatex( x, y, s );
}

void DrawLeftLatex( double x, double y , const char* s, float scale = 1, int color = 1){
  TLatex tltx; 
  tltx.SetTextFont(43); tltx.SetTextSize((int)(32*scale)); tltx.SetTextColor(color);
  tltx.SetNDC(kTRUE);
  tltx.SetTextAlign(12);
  tltx.DrawLatex( x, y, s );
}

void DrawCenterLatex( double x, double y , const char* s, float scale = 1, int color = 1){
  TLatex tltx; 
  tltx.SetTextFont(43); tltx.SetTextSize((int)(32*scale)); tltx.SetTextColor(color);
  tltx.SetNDC(kTRUE);
  tltx.SetTextAlign(22);
  tltx.DrawLatex( x, y, s );
}

// ============ DATA ================

void DrawAtlasInternalDataRight( double x0, double y0,
				 double scale, bool is_pPb ){
  DrawRightLatex(0.875 + x0, 0.95,
		 "#bf{#font[72]{ATLAS}} Internal", scale);
  if( is_pPb ){
    DrawRightLatex(0.875 + x0, 0.88 + y0, 
		   Form("#it{p}+Pb 2016, %i #mub^{-1}",
			pPbLumi2016), scale);
  } else {
    DrawRightLatex(0.875 + x0, 0.88 + y0, 
		   Form("#it{pp} 2015, %i pb^{-1}",
			ppLumi2015), scale);
  }
  DrawRightLatex(0.875 + x0, 0.81 + y0, 
		 "#sqrt{s_{NN}}=5.02 TeV", scale);
}

void DrawAtlasInternalDataLeft( double x0, double y0,
				double scale, bool is_pPb ){
  DrawRightLatex(0.875 + x0, 0.95, 
		 "#bf{#font[72]{ATLAS}} Internal", scale);
  if( is_pPb ){
    DrawRightLatex(0.875 + x0, 0.88 + y0, 
		   Form("#it{p}+Pb 2016, %i #mub^{-1}",
			pPbLumi2016), scale);
  } else {
    DrawRightLatex(0.875 + x0, 0.88 + y0, 
		   Form("#it{pp} 2015, %i pb^{-1}",
			ppLumi2015), scale);
  }
  DrawLeftLatex(0.18 + x0, 0.81 + y0, 
		"#sqrt{s_{NN}}=5.02 TeV", scale); 
}

// ============ MC ================

void DrawAtlasInternalMC( double x0, double y0, bool isReco, double scale ){ 
  DrawRightLatex(0.875 + x0, 0.95, 
		 "#bf{#font[72]{ATLAS}} Simulation Internal", scale);
  if( isReco  ) DrawRightLatex(0.875 + x0, 0.815,
			       "Reco Level", scale);
  if( !isReco ) DrawRightLatex(0.875 + x0, 0.815, 
			       "Truth Level", scale);
}

void SetCustomMarkerStyle( TH1* his , int iflag ){
  	
  //Set Color
  if( iflag == 0 ){
    his->SetLineColor(kBlack);
    his->SetLineWidth(2);
    his->SetMarkerColor(kBlack);
    his->SetMarkerStyle(20);
    his->SetMarkerSize(1.2);
  } 
  else if(iflag == 1 ){
    his->SetLineColor(kRed);
    his->SetLineWidth(2);
    his->SetMarkerColor(kRed);
    his->SetMarkerStyle(21);
    his->SetMarkerSize(1.1);
  }
  else if(iflag == 2 ){
    his->SetLineColor(kAzure-3);
    his->SetLineWidth(2);
    his->SetMarkerColor(kAzure-3);
    his->SetMarkerStyle(33);
    his->SetMarkerSize(1.8);
  }
  else if(iflag == 3 ){
    his->SetLineColor(kSpring-6);
    his->SetLineWidth(2);
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
    his->SetLineWidth(2);
    his->SetMarkerColor(kViolet);
    his->SetMarkerStyle(22);
    his->SetMarkerSize(1.6);
  }
  else if(iflag == 6 ){
    his->SetLineColor(kCyan+1);
    his->SetLineWidth(2);
    his->SetMarkerColor(kCyan+1);
    his->SetMarkerStyle(23);
    his->SetMarkerSize(1.6);
  }
  else if(iflag == 7 ){
    his->SetLineColor(46);
    his->SetLineWidth(2);
    his->SetMarkerColor(46);
    his->SetMarkerStyle(21);
    his->SetMarkerSize(1.4);
  }
  else if(iflag == 8 ){
    his->SetLineColor(kRed+3);
    his->SetLineWidth(2);
    his->SetMarkerColor(kRed+3);
    his->SetMarkerStyle(25);
    his->SetMarkerSize(1.4);
  } 

}

void SetHStyle( TH1* his, int iflag, float scale)
{
  his->SetLineWidth(2);
  his->SetStats(0);

  his->SetTitleFont( 43, "t");
  his->SetTitleSize( (int)(32 * scale), "t"); 

  his->SetLabelFont( 43, "xyz");
  his->SetLabelSize( (int)(30 * scale), "xyz");
  his->SetTitleFont( 43, "xyz");
  his->SetTitleSize( (int)(32 * scale), "xyz");
  his->SetTitleOffset(1.2, "xy");
  
  SetCustomMarkerStyle( his, iflag );  
}

void SetHStyle( TF1* his, int iflag)
{
  his->SetLineWidth(1);
  his->GetXaxis()->SetLabelFont(43);
  his->GetXaxis()->SetLabelSize(30);
  his->GetXaxis()->SetTitleSize(32);
  his->GetXaxis()->SetTitleOffset(1.2);
  his->GetXaxis()->SetTitleFont(43);
  his->GetYaxis()->SetLabelFont(43);
  his->GetYaxis()->SetLabelSize(30);
  his->GetYaxis()->SetTitleSize(32);
  his->GetYaxis()->SetTitleOffset(1.2);
  his->GetYaxis()->SetTitleFont(43);
  
  // SetCustomMarkerStyle( his, iflag ); 
}


void SetHStyle_open( TH1* his, int iflag)
{
  his->SetLineWidth(2);
  his->SetStats(0);
  his->GetXaxis()->SetLabelFont(43);
  his->GetXaxis()->SetLabelSize(30);
  his->GetXaxis()->SetTitleSize(32);
  his->GetXaxis()->SetTitleOffset(1.2);
  his->GetXaxis()->SetTitleFont(43);
  his->GetYaxis()->SetLabelFont(43);
  his->GetYaxis()->SetLabelSize(30);
  his->GetYaxis()->SetTitleSize(32);
  his->GetYaxis()->SetTitleOffset(1.2);
  his->GetYaxis()->SetTitleFont(43);
	
  SetCustomMarkerStyle( his, iflag );
}

void SetHStyle_graph( TGraphErrors* his, int iflag)
{
  his->SetLineWidth(2);
  //his->SetStats(0);
  his->GetXaxis()->SetLabelFont(43);
  his->GetXaxis()->SetLabelSize(30);
  his->GetXaxis()->SetTitleSize(32);
  his->GetXaxis()->SetTitleOffset(1.2);
  his->GetXaxis()->SetTitleFont(43);
  his->GetYaxis()->SetLabelFont(43);
  his->GetYaxis()->SetLabelSize(30);
  his->GetYaxis()->SetTitleSize(32);
  his->GetYaxis()->SetTitleOffset(1.2);
  his->GetYaxis()->SetTitleFont(43);
	
  //  SetCustomMarkerStyle( his, iflag );
}

void SetHStyle_TF1( TF1* his, int iflag)
{	
  //Set Color
  if( iflag == 0 ){
    his->SetLineColor(kBlack);
    his->SetLineWidth(2);
  } 
  else if(iflag == 1 ){
    his->SetLineColor(kRed+1);
    his->SetLineWidth(2);
  }
  else if(iflag == 2 ){
    his->SetLineColor(kBlue+1);
  }
  else if(iflag == 3 ){
    his->SetLineColor(kGreen+2);
    his->SetLineWidth(2);

  }
  else if(iflag == 4 ){
    his->SetLineColor(kOrange+1);
    his->SetLineWidth(2);
  }
  else if(iflag == 5 ){
    his->SetLineColor(kViolet);
  }
  else if(iflag == 6 ){
    his->SetLineColor(kCyan+1);
    his->SetLineWidth(2);
  }
  else if(iflag == 7 ){
    his->SetLineColor(46);
    his->SetLineWidth(2);
  }
  else if(iflag == 8 ){
    his->SetLineColor(kRed+3);
    his->SetLineWidth(2);
  } 
 
  
}

void SetHStyle_graph( TGraphAsymmErrors* his, int iflag)
{
  his->SetLineWidth(2);
  //his->SetStats(0);
  his->GetXaxis()->SetLabelFont(43);
  his->GetXaxis()->SetLabelSize(30);
  his->GetXaxis()->SetTitleSize(32);
  his->GetXaxis()->SetTitleOffset(1.2);
  his->GetXaxis()->SetTitleFont(43);
  his->GetYaxis()->SetLabelFont(43);
  his->GetYaxis()->SetLabelSize(30);
  his->GetYaxis()->SetTitleSize(32);
  his->GetYaxis()->SetTitleOffset(1.2);
  his->GetYaxis()->SetTitleFont(43);
	
  //Set Color
  if( iflag == 0 ){
    his->SetLineColor(kBlack);
    his->SetFillColor(kGray);
    his->SetLineWidth(2);
    his->SetMarkerColor(kBlack);
    his->SetMarkerStyle(20);
    his->SetMarkerSize(1.2);
    //his->SetFillStyle(3002);
  } 
  else if(iflag == 1 ){
    his->SetLineColor(kRed+1);
    his->SetFillColor(kRed-10);
    his->SetLineWidth(2);
    his->SetMarkerColor(kRed+1);
    his->SetMarkerStyle(21);
    his->SetMarkerSize(1.2);
    //his->SetFillStyle(3002);
  }
  else if(iflag == 2 ){
    his->SetLineColor(kAzure-3);
    his->SetFillColor(kAzure-9);
    his->SetLineWidth(2);
    his->SetMarkerColor(kAzure-3);
    his->SetMarkerStyle(33);
    his->SetMarkerSize(1.5);
    //his->SetFillStyle(3002);
  }
  else if(iflag == 3 ){
    his->SetLineColor(kSpring-6);
    //his->SetFillColor(kSpring+1);
    his->SetFillColor(kGreen-10);
    his->SetLineWidth(2);
    his->SetMarkerColor(kSpring-6);
    his->SetMarkerStyle(34);
    his->SetMarkerSize(1.5);
    //his->SetFillStyle(3002);
  }
  else if(iflag == 4 ){
    his->SetLineColor(kOrange+1);
    his->SetFillColor(kOrange+1);
    his->SetLineWidth(2);
    his->SetMarkerColor(kOrange+1);
    his->SetMarkerStyle(29);
    his->SetMarkerSize(1.6);
    his->SetFillStyle(3002);
  }
  else if(iflag == 5 ){
    his->SetLineColor(kViolet);
    his->SetFillColor(kViolet);
    his->SetLineWidth(2);
    his->SetMarkerColor(kViolet);
    his->SetMarkerStyle(22);
    his->SetMarkerSize(1.6);
    his->SetFillStyle(3002);
  }
  else if(iflag == 6 ){
    his->SetLineColor(kCyan+1);
    his->SetFillColor(kCyan+1);
    his->SetLineWidth(2);
    his->SetMarkerColor(kCyan+1);
    his->SetMarkerStyle(23);
    his->SetMarkerSize(1.6);
    his->SetFillStyle(3002);
  }
  else if(iflag == 7 ){
    his->SetLineColor(46);
    his->SetLineWidth(2);
    his->SetMarkerColor(46);
    his->SetMarkerStyle(21);
    his->SetMarkerSize(1.4);
  }
  else if(iflag == 8 ){
    his->SetLineColor(kRed+3);
    his->SetLineWidth(2);
    his->SetMarkerColor(kRed+3);
    his->SetMarkerStyle(25);
    his->SetMarkerSize(1.4);
  }  
 
  
}

void SetHStyle_graph_open( TGraphErrors* his, int iflag)
{
  his->SetLineWidth(2);
  //his->SetStats(0);
  his->GetXaxis()->SetLabelFont(43);
  his->GetXaxis()->SetLabelSize(30);
  his->GetXaxis()->SetTitleSize(32);
  his->GetXaxis()->SetTitleOffset(1.2);
  his->GetXaxis()->SetTitleFont(43);
  his->GetYaxis()->SetLabelFont(43);
  his->GetYaxis()->SetLabelSize(30);
  his->GetYaxis()->SetTitleSize(32);
  his->GetYaxis()->SetTitleOffset(1.2);
  his->GetYaxis()->SetTitleFont(43);
	
  //Set Color
  if( iflag == 0 ){
    his->SetLineColor(kBlack);
    his->SetLineWidth(2);
    his->SetMarkerColor(kBlack);
    his->SetMarkerStyle(24);
    his->SetMarkerSize(1.2);
  } 
  else if(iflag == 1 ){
    his->SetLineColor(kRed+1);
    his->SetLineWidth(2);
    his->SetMarkerColor(kRed+1);
    his->SetMarkerStyle(25);
    his->SetMarkerSize(1.2);
  }
  else if(iflag == 2 ){
    his->SetLineColor(kBlue+1);
    his->SetLineWidth(2);
    his->SetMarkerColor(kBlue+1);
    his->SetMarkerStyle(26);
    his->SetMarkerSize(1.5);
  }
  else if(iflag == 3 ){
    his->SetLineColor(kGreen+2);
    his->SetLineWidth(2);
    his->SetMarkerColor(kGreen+2);
    his->SetMarkerStyle(32);
    his->SetMarkerSize(1.5);
  }
  else if(iflag == 4 ){
    his->SetLineColor(kOrange+1);
    his->SetLineWidth(2);
    his->SetMarkerColor(kOrange+1);
    his->SetMarkerStyle(30);
    his->SetMarkerSize(1.6);
  }
  else if(iflag == 5 ){
    his->SetLineColor(kViolet);
    his->SetLineWidth(2);
    his->SetMarkerColor(kViolet);
    his->SetMarkerStyle(27);
    his->SetMarkerSize(1.6);
  }
  else if(iflag == 6 ){
    his->SetLineColor(kCyan+1);
    his->SetLineWidth(2);
    his->SetMarkerColor(kCyan+1);
    his->SetMarkerStyle(28);
    his->SetMarkerSize(1.6);
  }
  
}


void SetLegendStyle(TLegend * legend, float scale)
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
