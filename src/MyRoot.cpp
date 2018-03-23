#include <boost/format.hpp>
#include <boost/filesystem.hpp>
#include <boost/assign.hpp>

#include <TVirtualFitter.h>

#include "MyRoot.h"

#include <iostream>
#include <sstream>

const int pPbLumi2016 = 437; // ub
const int ppLumi2015  = 26;  // pb

//===================================
//         COMMON FUNCTIONS
//===================================
bool CT::AnalysisTools::EpsilonEqual( double a, double b )
{ return fabs( a - b ) < constants::EPSILON; }

// returns dphi in range 0<dphi<2pi
double CT::AnalysisTools::DPhiFC
( const TLorentzVector& jet1, const TLorentzVector& jet2 ){

  double phi1 = jet1.Phi(); double phi2 = jet2.Phi();
  double deltaPhi = phi1 - phi2;
  while( deltaPhi < 0 ) deltaPhi += 2*constants::PI;
  return deltaPhi;
}

double CT::AnalysisTools::DeltaPhi
( const TLorentzVector& jet1, const TLorentzVector& jet2 ){

  double phi1 = jet1.Phi(); double phi2 = jet2.Phi();
  double deltaPhi = TMath::Abs(phi1 - phi2);
  if( deltaPhi > constants::PI ){ deltaPhi = 2*constants::PI - deltaPhi; };
  return deltaPhi;
}

double CT::AnalysisTools::DeltaR
( const TLorentzVector& jet1, const TLorentzVector& jet2 ){
  
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

std::pair< double, double > CT::AnalysisTools::GetRMS
( TH1* h, double xLow, double xUp, double ref ){

  // rms is standard rms by definition
  // rmsError is yi * xi ^ 2 * dyi 
  double yixi2         = 0;
  double yixi2Error    = 0;
  
  for ( int xBin = 1; xBin <= h->GetNbinsX(); xBin++ ){

    double xBinCenter = h->GetBinCenter( xBin );
        if ( xBinCenter < xLow ||
	 xBinCenter > xUp ) continue;
    
    double yi  = h->GetBinContent( xBin );
    double eYi = h->GetBinError  ( xBin );
    
    // avoid negative values
    if ( yi < 0 ){ yi = 0; }

    yixi2      += yi * std::pow
      ( ( xBinCenter - ref ), 2 );

    yixi2Error += yi * std::pow
      ( ( xBinCenter - ref ), 2 ) * eYi;
  }

  int xBinLow = h->FindBin( xLow );
  int xBinUp  = h->FindBin( xUp  );
  
  double integral      = 0;
  double integralError = 0;

  integral = h->IntegralAndError
    ( xBinLow, xBinUp, integralError );

  double rms      = 0;
  double rmsError = 0;

  if ( integral > 0 && yixi2 > 0 ){
    rms      = std::sqrt( yixi2 / integral);
    rmsError = std::sqrt( std::pow( rmsError / integral, 2 ) +
			  std::pow( rms * integralError / integral , 2 ) );
      
    std::cout << " ---- " << integral << " " << integralError
	      << " " << rms << " " << rmsError << std::endl;
  } else {
    rms      = 0;
    rmsError = 0;
  }

  return std::make_pair( rms, rmsError );
}

bool CT::AnalysisTools::AverageOver( TH2* h2,
				     const std::string& rc ){

  TAxis* a1 = NULL;
  TAxis* a2 = NULL;

  // i, j are the loop indexes.
  // i->a1, j->a2. i.e. for each bin in
  // a1 we sum over all other bins in a2 (j)
  // x or y can point to either i or j
  int  i, j;;
  int* x = NULL; int* y = NULL;

  // average over rows. a1,i are y,
  // sum over a2,j which are x
  if( !rc.compare( "row" ) ){
    a1 = h2->GetYaxis(); y = &i;
    a2 = h2->GetXaxis(); x = &j;
  } else if( !rc.compare( "column" ) ){
    a1 = h2->GetXaxis(); x = &i;
    a2 = h2->GetYaxis(); y = &j;
  } else{
    return false;
  }

  for( i = 1; i <= a1->GetNbins(); i++ ){
    double iSum = 0;
    for( j = 1; j <= a2->GetNbins(); j++ )
      { iSum += h2->GetBinContent( *x, *y ); }

    if( !iSum ){ continue; }
    
    for( j = 1 ; j <= a2->GetNbins(); j++ ){
      double fraction = h2->GetBinContent( *x, *y ) / iSum;
      h2->SetBinContent( *x, *y, fraction );
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

std::vector<int> CT::AnalysisTools::vectoriseI
(TString str, TString sep) {

  std::vector<int> result;
   std::vector<std::string> vecS = vectorise(str,sep);
   for (uint i=0;i<vecS.size();++i)
     {result.push_back(std::stoi(vecS[i]));}
   return result;
}   

std::vector<double> CT::AnalysisTools::vectoriseD
(TString str, TString sep) {

  std::vector<double> result;
   std::vector<std::string> vecS = vectorise(str,sep);
   for (uint i=0;i<vecS.size();++i)
     {result.push_back(std::stod(vecS[i]));}
   return result;
}   

bool CT::AnalysisTools::SubtractCombinatoric( TH1* hProj, double xLow, double xHigh ){

  auto combFit = [&]( double* x, double* par){ return par[0]; };

  // to estimate some combinatoric
  // contribution on range [0,1], where
  // distribution has very few entries
  // only do if there are bins with entries
  // on that range.
  if( !hProj->Integral( 1, hProj->FindBin( 1. ) ) ){ return false; }

  TF1 cFit( "cFit", combFit, xLow, xHigh, 1 );
  hProj->Fit( "cFit", "NQ", "", xLow, xHigh ); 
    
  double comb      = cFit.GetParameter(0);
  double combError = cFit.GetParError (0);
    
  for( int xBin = 1; xBin <= hProj->GetNbinsX(); xBin++ ){
    double val      = hProj->GetBinContent( xBin );
    double binError = hProj->GetBinError  ( xBin );
    double newVal   = val - comb > 0 ? val - comb : 0;
    double newError = 
      TMath::Sqrt( binError*binError + combError*combError );	
    hProj->SetBinContent( xBin, newVal   );
    hProj->SetBinError  ( xBin, newError );
  }

  return true;
}

void CT::AnalysisTools::UndoWidthScaling( TH1* h ){

  // undo h->Scale( 1., "width" );
  for( int xBin = 1; xBin <= h->GetNbinsX(); xBin++ ){
    double binWidth   = h->GetBinWidth ( xBin );
    h->SetBinContent( xBin, h->GetBinContent( xBin ) * binWidth );
    h->SetBinError  ( xBin, h->GetBinError  ( xBin ) * binWidth );
  }
}

TF1* CT::AnalysisTools::FitDphi( TH1* histo, double xLow, double xHigh ){

  auto EMG = [&]( double* x, double* par){
    return par[0]*TMath::Exp(par[2]*par[2]/(2*par[1]*par[1])) *
    ( TMath::Exp((x[0]-constants::PI)/par[1]) * 0.5 *
      TMath::Erfc( (1/1.41) * (x[0]-constants::PI)/par[2] + par[2]/par[1]) +
      TMath::Exp((constants::PI-x[0])/par[1]) *
      ( 1 - 0.5 * TMath::Erfc( (1/1.41) * (x[0]-constants::PI)/par[2] - par[2]/par[1])));
  };
  
  // set range to be in range of fit
  histo->GetXaxis()->SetRange( 1, -1 );
  
  TF1* dPhiFit = new TF1( Form("f_%s", histo->GetName()), EMG, xLow, xHigh, 3 );

  if( !histo->GetEntries() )
    { return dPhiFit; }

  double amp = histo->GetBinContent( histo->GetMaximumBin () );
  
  dPhiFit->SetParameters( amp, 0.20, 0.20 );
  dPhiFit->SetParLimits ( 1, 1E-3, 1 );
  dPhiFit->SetParLimits ( 2, 1E-2, 1 );
  
  int status = histo->Fit( dPhiFit->GetName(), "NQI", "", xLow, xHigh );

  if( status ){ std::cout << " +++++++++++++ " << status
			  << " " << dPhiFit->GetName() << std::endl; }
  
  return dPhiFit;
}

TF1* CT::AnalysisTools::FitDphiWC( TH1* histo, double xLow, double xHigh ){

  auto EMG = [&]( double* x, double* par){
    return par[0]*TMath::Exp(par[2]*par[2]/(2*par[1]*par[1])) *
    ( TMath::Exp((x[0]-constants::PI)/par[1]) * 0.5 *
      TMath::Erfc( (1/1.41) * (x[0]-constants::PI)/par[2] + par[2]/par[1]) +
      TMath::Exp((constants::PI-x[0])/par[1]) *
      ( 1 - 0.5 * TMath::Erfc( (1/1.41) * (x[0]-constants::PI)/par[2] - par[2]/par[1])))
    + par[3];
  };
  
  // set range to be in range of fit
  histo->GetXaxis()->SetRange( 1, -1 );
  
  TF1* dPhiFit = new TF1( Form("f_%s", histo->GetName()), EMG, xLow, xHigh, 4 );

  if( !histo->GetEntries() )
    { return dPhiFit; }

  double amp = histo->GetBinContent( histo->GetMaximumBin () );
  
  dPhiFit->SetParameters( amp, 0.20, 0.1, 0);
  dPhiFit->SetParLimits ( 1, 1E-3, 1 );
  dPhiFit->SetParLimits ( 2, 1E-2, 1 );
  dPhiFit->SetParLimits ( 3, 0,    1 );
  
  int status = histo->Fit( dPhiFit->GetName(), "NQI", "", xLow, xHigh );

  if( status ){ std::cout << " +++++++++++++ " << status
			  << " " << dPhiFit->GetName() << std::endl; }
  
  return dPhiFit;
}

TF1* CT::AnalysisTools::FitDphi( TGraph* graph, double xLow, double xHigh ){

  auto EMG = [&]( double* x, double* par){
    return par[0]*TMath::Exp(par[2]*par[2]/(2*par[1]*par[1])) *
    ( TMath::Exp((x[0]-constants::PI)/par[1]) * 0.5 *
      TMath::Erfc( (1/1.41) * (x[0]-constants::PI)/par[2] + par[2]/par[1]) +
      TMath::Exp((constants::PI-x[0])/par[1]) *
      ( 1 - 0.5 * TMath::Erfc( (1/1.41) * (x[0]-constants::PI)/par[2] - par[2]/par[1])));
  };

  TF1* dPhiFit = new TF1( Form("f_%s", graph->GetName()), EMG, xLow, xHigh, 3 );

  if( !graph->GetN() )
    { return dPhiFit; }

  double amp = CT::AnalysisTools::GetGraphMax( graph, xLow, xHigh );
  
  dPhiFit->SetParameters( amp, 0.20, 0.20 );
  dPhiFit->SetParLimits ( 2, 1E-2, 1 );
  
  int status = graph->Fit( dPhiFit->GetName(), "NQEX0", "", xLow, xHigh );

  if( status ){ std::cout << " +++++++++++++ " << status
			  << " " << dPhiFit->GetName() << std::endl; }
  
  return dPhiFit;
}

TF1* CT::AnalysisTools::FitGaussian( TH1* histo, double xLow, double xHigh){

  TF1* fit  = new TF1( Form("f_%s", histo->GetName()), "gaus(0)" );

  double hXmin  = histo->GetXaxis()->GetXmin();
  double hXmax  = histo->GetXaxis()->GetXmax();
  
  // fit once 
  double mean   = histo->GetMean();
  double rms    = histo->GetRMS();

  double fitMin = mean - 2.0 * rms;
  double fitMax = mean + 2.0 * rms;

  if( histo->GetEntries() < 5 || fitMin < hXmin || fitMax > hXmax )
    { return fit; }
  
  histo->Fit( fit->GetName(), "NQ", "", fitMin, fitMax );

  // fit second time with better parameters
  mean   = fit->GetParameter(1);
  rms    = fit->GetParameter(2);

  fitMin = mean - 2.0 * rms;
  fitMax = mean + 2.0 * rms;

  if( histo->GetEntries() < 5 || fitMin < hXmin || fitMax > hXmax )
    { return fit; }
  
  fit->SetRange( fitMin, fitMax );
  
  histo->Fit( fit->GetName(), "NQ", "", fitMin, fitMax );

  return fit;
}

TF1* CT::AnalysisTools::FitPol1( TH1* histo, double xLow, double xHigh){
  TF1* fit  = new TF1( Form("f_%s", histo->GetName()), "pol1(0)", xLow, xHigh );

  fit->SetParameters( 1, 1 );

  if( histo->GetEntries() > 5 ){
    histo->Fit( fit->GetName(), "NQ", "", xLow, xHigh );
  }
  
  return fit;
}

TF1* CT::AnalysisTools::FitPol2( TH1* histo, double xLow, double xHigh){
  TF1* fit  = new TF1( Form("f_%s", histo->GetName()), "pol2(0)", xLow, xHigh );

  fit->SetParameters( 1, 1, 1 );

  if( histo->GetEntries() > 5 ){
    histo->Fit( fit->GetName(), "NQ", "", xLow, xHigh );
  }
  
  return fit;
}

TF1* CT::AnalysisTools::FitLogPol2( TH1* histo, double xLow, double xHigh){

  auto logPol2 = [&]( double* x, double* par){
    return par[0] * TMath::Log( x[0] ) +
    par[1] * TMath::Power( TMath::Log( x[1] ), 2 );
  };
  
  TF1* fit  = new TF1( Form("f_%s", histo->GetName()), logPol2, xLow, xHigh, 2 );

  fit->SetParameters( 1, 1 );

  histo->Fit( fit->GetName(), "NQ", "", xLow, xHigh );
  
  return fit;
}

double CT::AnalysisTools::GetGraphMax( TGraph* g, double xMin, double xMax ){

  double yMax = -1;

  for( int i = 0; i < g->GetN(); i++ ){
    double x = g->GetX()[i];
    double y = g->GetY()[i];
    if( x > xMax || x < xMin ){ continue; }
    if( y > yMax ){ yMax = y; }
  }

  return yMax;
}

TGraph* CT::AnalysisTools::Barycenters( TH1* hEnt, TH1* hSum ){

  TGraphAsymmErrors* g = new TGraphAsymmErrors( hEnt );

  for( int i = 0; i < g->GetN(); i++ ){
    double gX = g->GetX()[ i ];
    double gY = g->GetY()[ i ];

    int  xBin = hEnt->FindBin( gX );

    double sum         = hSum->GetBinContent( xBin );
    double gXnew       = sum / gY;
    double deltaCenter = gXnew - hEnt->GetBinCenter( xBin );
    double binWidth    = hEnt->GetBinWidth( xBin );
    
    g->SetPoint( i, gXnew, gY );
    g->SetPointEXlow ( i, 0.5 * binWidth + deltaCenter );
    g->SetPointEXhigh( i, 0.5 * binWidth - deltaCenter );
  }

  return g;
}

void CT::AnalysisTools::MatchGraphGraphX( TGraph* g1, TGraph * g2 ){

  TGraphAsymmErrors* gg1 = dynamic_cast< TGraphAsymmErrors* >( g1 );
  TGraphAsymmErrors* gg2 = dynamic_cast< TGraphAsymmErrors* >( g2 );

  if( !gg1 || !gg2 ){ return; }
  
  for( int i  = 0; i < gg1->GetN(); i++ ){
    gg1->SetPoint( i, gg2->GetX()[i], gg1->GetY()[i] );
    gg1->SetPointEYlow ( i, gg2->GetErrorXlow ( i ) );
    gg1->SetPointEYhigh( i, gg2->GetErrorXhigh( i ) );
  }
}

void CT::AnalysisTools::MatchGraphHistoY( TGraph* g, TH1* h ){

  TGraphAsymmErrors* gg = dynamic_cast< TGraphAsymmErrors* >( g );

  if( !gg ){ return; }
  
  for( int i  = 0; i < gg->GetN(); i++ ){
    double x  = gg->GetX()[i];
    int xBin  = h->FindBin( x );
    double y  = h->GetBinContent( xBin );
    double ey = h->GetBinError( xBin );
    if( y < 0 ){ gg->RemovePoint( i ); }
    gg->SetPoint( i, x, y );
    gg->SetPointEYlow ( i, ey );
    gg->SetPointEYhigh( i, ey );
  }
}

void CT::AnalysisTools::GetBinRange( TAxis* a,
				     int b1, int b2,
				     double& x1, double& x2 ){
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
					      std::string label){
  std::stringstream ss;
  if( is_pPb )
    { ss << boost::format("%3.2g<%s<%3.2g") % ystarMin % label % ystarMax;}
  else {
    if( std::abs(ystarMin) > std::abs(ystarMax) )
      { ss << boost::format("%3.2g<|%s|<%3.2g")
	  % std::abs(ystarMax) % label % std::abs(ystarMin); }
    else
      { ss << boost::format("%3.2g<|%s|<%3.2g")
	  % std::abs(ystarMin) % label % std::abs(ystarMax); }
  }
  return ss.str();
}

std::string CT::AnalysisTools::GetLabel
( double vMin,double vMax,
  const std::string& var, const std::string& unitIn ){

  // for now, want [GeV] by all pT labels
  std::string unit = unitIn;
  if( var.find("{p}_{T}") != std::string::npos )
    { unit = "[GeV]"; }
  
  std::stringstream ss;

  if( vMin != vMax ){
    ss << boost::format("%3.2g<%s<%3.2g")
      % vMin % var % vMax;
  } else {
    ss << boost::format("%s>%3.2g")
      % var % vMin ;  
  }
  if( !unit.empty() ){ ss << boost::format(" %s") % unit; }
  
  return ss.str();
}

double CT::AnalysisTools::GetLogMaximum( double max ){
  double power = std::log10(max);
  power = std::ceil(power);
  return std::pow( 10, power );
}

double CT::AnalysisTools::GetLogMinimum( double min ){
  double power = std::log10(min);
  power = std::floor(power);
  return std::pow( 10, power );
}
void CT::AnalysisTools::ResetAxisRanges( TH1* h ){

  h->GetXaxis()->SetRange( 0, 0 );
  h->GetYaxis()->SetRange( 0, 0 );
  h->GetZaxis()->SetRange( 0, 0 );
}

void CT::AnalysisTools::ResetAxisRanges( THnSparse* h ){

  for( int i = 0; i < h->GetNdimensions(); i++ ){
    h->GetAxis( i )->SetRange( 0, 0 );
  }
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

std::vector<std::string> CT::AnalysisTools::ListFiles
(const char* dirname, const char* ext){

  std::vector<std::string> vFiles;
   TSystemDirectory dir(dirname, dirname);
   TList *files = dir.GetListOfFiles();
   if (files) {
      TSystemFile *file;
      TString fname;
      TIter next(files);
      while ((file=(TSystemFile*)next())) {
         fname = file->GetName();
         if (!file->IsDirectory() && fname.EndsWith(ext)) {
	   vFiles.push_back( fname.Data() );
         }
      }
   }
   return vFiles;
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
    his->SetMarkerSize(1.4);
  } else if(iflag == 1 ){
    his->SetLineColor(kRed);
    his->SetMarkerColor(kRed);
    his->SetMarkerStyle(21);
    his->SetMarkerSize(1.3);
  } else if(iflag == 2 ){
    his->SetLineColor(kAzure-3);
    his->SetMarkerColor(kAzure-3);
    his->SetMarkerStyle(33);
    his->SetMarkerSize(2.0);
  } else if(iflag == 3 ){
    his->SetLineColor(kSpring-6);
    his->SetMarkerColor(kSpring-6);
    his->SetMarkerStyle(34);
    his->SetMarkerSize(1.7);
  } else if(iflag == 4 ){
    his->SetLineColor(kOrange+1);
    his->SetLineWidth(2);
    his->SetMarkerColor(kOrange+1);
    his->SetMarkerStyle(29);
    his->SetMarkerSize(1.8);
  } else if(iflag == 5 ){
    his->SetLineColor(kBlack);
    his->SetMarkerColor(kBlack);
    his->SetMarkerStyle(24);
    his->SetMarkerSize(1.5);
  } else if(iflag == 6 ){
    his->SetLineColor(kRed);
    his->SetMarkerColor(kRed);
    his->SetMarkerStyle(25);
    his->SetMarkerSize(1.4);
  } else if(iflag == 7 ){
    his->SetLineColor(kAzure-3);
    his->SetMarkerColor(kAzure-3);
    his->SetMarkerStyle(27);
    his->SetMarkerSize(2.1);
  } else if(iflag == 8 ){
    his->SetLineColor(kSpring-6);
    his->SetMarkerColor(kSpring-6);
    his->SetMarkerStyle(28);
    his->SetMarkerSize(1.8);
  } else if(iflag == 9 ){
    his->SetLineColor(kOrange+1);
    his->SetLineWidth(2);
    his->SetMarkerColor(kOrange+1);
    his->SetMarkerStyle(30);
    his->SetMarkerSize(1.9);
  }
}

void CT::StyleTools::SetCustomMarkerStyle( TGraph* graph , int iflag ){
  //Set Color
  graph->SetLineWidth(2);
  if( iflag == 0 ){
    graph->SetLineColor(kBlack);
    graph->SetMarkerColor(kBlack);
    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(1.5);
    graph->SetFillColor( kGray + 1 );
    graph->SetFillStyle( 3001 );
  } 
  else if(iflag == 1 ){
    graph->SetLineColor(kRed);
    graph->SetMarkerColor(kRed);
    graph->SetMarkerStyle(21);
    graph->SetMarkerSize(1.4);
    graph->SetFillColor( kRed - 4  );
    graph->SetFillStyle( 3001 );
  }
  else if(iflag == 2 ){
    graph->SetLineColor( kAzure - 3 );
    graph->SetMarkerColor( kAzure - 3 );
    graph->SetMarkerStyle(33);
    graph->SetMarkerSize(2.1);
    graph->SetFillColor( kAzure - 4 );
    graph->SetFillStyle( 3001 );
  }
  else if(iflag == 3 ){
    graph->SetLineColor( kSpring - 6 );
    graph->SetMarkerColor( kSpring - 6 );
    graph->SetMarkerStyle(34);
    graph->SetMarkerSize(1.8);
    graph->SetFillColor( kSpring - 5 );
    graph->SetFillStyle( 3001 );
  }
  else if(iflag == 4 ){
    graph->SetLineColor(kBlack);
    graph->SetMarkerColor(kBlack);
    graph->SetMarkerStyle(24);
    graph->SetMarkerSize(1.4);
    graph->SetFillColor( kGray + 1 );
    graph->SetFillStyle( 3002 );
  }
  else if(iflag == 5 ){
    graph->SetLineColor(kRed);
    graph->SetMarkerColor(kRed);
    graph->SetMarkerStyle(25);
    graph->SetMarkerSize(1.7);
    graph->SetFillColor( kRed - 4 );
    graph->SetFillStyle( 3002 );
   }  
  else if(iflag == 6 ){
    graph->SetLineColor( kAzure - 3 );
    graph->SetMarkerColor( kAzure - 3 );
    graph->SetMarkerStyle(27);
    graph->SetMarkerSize(2.1);
    graph->SetFillColor( kAzure -4 );
    graph->SetFillStyle( 3002 );
  }
  else if(iflag == 7 ){
    graph->SetLineColor( kSpring - 6 );
    graph->SetMarkerColor( kSpring - 6 );
    graph->SetMarkerStyle(28);
    graph->SetMarkerSize(1.8);
    graph->SetFillColor( kSpring - 5 );
    graph->SetFillStyle( 3002 );
  }
}

void CT::StyleTools::SetHStyle( TH1* his, int iflag, double scale)
{
  
  his->GetXaxis()->SetNdivisions( 504 );
  his->GetYaxis()->SetNdivisions( 504 );
  
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

  graph->GetXaxis()->SetNdivisions( 504 );
  graph->GetYaxis()->SetNdivisions( 504 );
  
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

void CT::StyleTools::SetHStyleRatio( TH1* his, int iflag, double scale ){
  SetHStyle( his, iflag );
  his->SetTitle("");
  his->SetMaximum( 1.5 );
  his->SetMinimum( 0.5 );
  his->GetYaxis()->SetNdivisions(503);
}

TH1F* CT::StyleTools::SetCStyleGraph ( TPad& c,
				       double x0, double y0,
				       double x1, double y1,
				       const std::string& title ){
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

void CT::StyleTools::HideAxis( TH1* h, const std::string& xyz ){
  if( !xyz.compare("x") ){
    h->GetXaxis()->SetLabelOffset(999);
    h->GetXaxis()->SetLabelSize(0);
  } else if( !xyz.compare("y") ){
    h->GetYaxis()->SetLabelOffset(999);
    h->GetYaxis()->SetLabelSize(0);
  } // continue like this
}

const double CT::StyleTools::lSS = 0.95;

const double CT::StyleTools::hSS = 1.00;

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
    (0.88 , 0.96, "#bf{#font[72]{ATLAS}} Internal", scale, 1 );
}

void CT::DrawTools::DrawAtlasInternalDataRight
( double x0, double y0, bool is_pPb, double scale ){
  DrawRightLatex
    (0.90 , 0.96,"#bf{#font[72]{ATLAS}} Internal", scale, 1 );
  if( is_pPb ){
    DrawRightLatex
      (0.87 + x0, 0.87 + y0, Form("#it{p}+Pb 2016, %i #mub^{-1}",
			     pPbLumi2016), scale, 1 );
  } else {
    DrawRightLatex
      (0.87 + x0, 0.87 + y0,Form("#it{pp} 2015, %i pb^{-1}",
			ppLumi2015), scale, 1 );
  }
  DrawRightLatex(0.87 + x0, 0.79 + y0, 
		 "#sqrt{s_{NN}}=5.02 TeV", scale, 1 );
}

void CT::DrawTools::DrawAtlasInternalDataLeft
( double x0, double y0, bool is_pPb, double scale ){
  DrawRightLatex(0.90, 0.96, 
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
  DrawLeftLatex(0.18 + x0, 0.79 + y0, 
		"#sqrt{s_{NN}}=5.02 TeV", scale, 1 ); 
}

// ============ MC ================

void CT::DrawTools::DrawAtlasInternalMCRight
( double x0, double y0, const std::string& mcType, int mode, double scale ){ 

  std::string system = "";
  
  // !is_pPb => mode = 0
  if( mode == 0 ){
    system = "#it{pp}";
  } else if( mode == 1 ){
    system = "#it{p}+Pb";
  }
  
  DrawRightLatex(0.88, 0.96, 
		 "#bf{#font[72]{ATLAS}} Simulation Internal", scale, 1 );
  DrawRightLatex( 0.87 + x0, 0.87,
		  Form( "%s %s", system.c_str(), mcType.c_str() ), scale, 1 );
 }


void CT::DrawTools::DrawAtlasInternalMCLeft
( double x0, double y0, const std::string& mcType, int mode, double scale ){ 

  std::string system = "";
  
  // !is_pPb => mode = 0
  if( mode == 0 ){
    system = "#it{pp}";
  } else if( mode == 1 ){
    system = "#it{p}+Pb";
  }

  DrawRightLatex(0.88, 0.96, 
		 "#bf{#font[72]{ATLAS}} Simulation Internal", scale, 1 );

  DrawLeftLatex( 0.18 + x0, 0.87,
		 Form( "%s %s", system.c_str(), mcType.c_str() ), scale, 1 );
}
