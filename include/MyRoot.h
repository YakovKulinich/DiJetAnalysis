#ifndef MYROOT_H
#define MYROOT_H

#include <TFile.h>
#include <TKey.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TH2D.h>
#include <TH3D.h>
#include <THStack.h>
#include <THnSparse.h>
#include <TF1.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TMarker.h>
#include <TLine.h>
#include <TText.h>
#include <TLatex.h>
#include <TLorentzVector.h>
#include <TLegend.h>
#include <TMath.h>
#include <TTree.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>

#include <iostream>
#include <string>
#include <sstream>


//===================================
//         CONSTANTS
//===================================

namespace constants{
  // FCal eta ranges
  const double FETAMIN = 3.2;
  const double FETAMAX = 4.9;
  // barrel range
  const double CETAMAX = 3.2;

  // betaz shift in pPb
  const double BETAZ = 0.5;
  
  // Forward ranges in ystar
  const double FYSTARMIN = FETAMIN - BETAZ;
  const double FYSTARMAX = FETAMAX - BETAZ;
  // central range in ystar
  const double CYSTARMAX = CETAMAX - BETAZ;
  
  // atlas maximum eta
  const double ETAMIN = -5.0;
  const double ETAMAX = 5.0;

  const double PI = TMath::Pi();

  const double EPSILON = 1e-6;
}

namespace CT{
  //===================================
  //         COMMON FUNCTIONS
  //===================================
  class AnalysisTools{
  public:
    bool EpsilonEqual( double a, double b );

    // returns dphi in range 0<dphi<2pi
    double DPhiFC( double phi1, double phi2 );

    double DeltaPhi( double phi1, double phi2 );

    double DeltaR(  const TLorentzVector& jet1,
		    const TLorentzVector& jet2 );

    static bool sortByDecendingPt
      ( const TLorentzVector& jet1,
	const TLorentzVector& jet2 );

    bool TruncateHistoBins( THnSparse* , THnSparse* );
    
    bool TruncateHistoBins( TH3* h3 );

    bool DoPrint( int ev );

    std::vector<std::string> vectorise(TString str, TString sep);

    std::vector<double> vectoriseD(TString str, TString sep);

    TF1* FitGaussian( TH1* hProj, double = 0, double = 0 );

    TF1* FitDphi( TH1* hProj, double = 0, double = 0 );

    void GetBinRange( TAxis*, int, int, double&, double& );

    std::string GetName( double, double, const std::string& );

    std::string GetEtaLabel( double, double, bool = false );

    std::string GetYstarLabel( double, double, bool = false, const std::string = "#it{y}*" );

    std::string GetLabel( double, double,
			const std::string& = "",
			const std::string& = "" );

    // this shuold be in miscallaneous
    void CheckWriteDir( const char* ); 
  };

  //===================================
  //          STYLE STUFF
  //===================================
  class StyleTools{
  public:
    static const double lSS; 
    static const double hSS;
  
    // Histogram/Graph/Function styles
    void SetCustomMarkerStyle( TH1* his    , int iflag );
    void SetCustomMarkerStyle( TGraph* his , int iflag );
    
    void SetHStyle( TH1*    his, int iflag, double scale = hSS);
    void SetHStyle( TGraph* his, int iflag, double scale = hSS);
    void SetHStyle( TF1*    his, int iflag, double scale = hSS);

    // For Canvas where you draw efficiency
    TH1F* SetCStyleEff( TCanvas&, double, double, double, double,
			const std::string& );
  
    // Legend Style
    void SetLegendStyle(TLegend * legend, double scale = lSS);
  };
  
  //===================================
  //          DRAWING STUFF
  //===================================
  class DrawTools{
  public:
    void DrawRightLatex ( double x, double y, const std::string&,
			  double scale = StyleTools::lSS,
			  int color = 1 );
    void DrawLeftLatex  ( double x, double y, const std::string&,
			  double scale = StyleTools::lSS,
			  int color = 1 );
    void DrawCenterLatex( double x, double y, const std::string&,
			  double scale = StyleTools::lSS,
			  int color = 1 );
    void DrawAtlasInternal( double scale = StyleTools::lSS );

    void DrawAtlasInternalDataRight( double x0, double y0, bool is_pPb,
				     double scale = StyleTools::lSS );
    void DrawAtlasInternalDataLeft ( double x0, double y0, bool is_pPb,
				     double scale = StyleTools::lSS );
    void DrawAtlasInternalMCRight  ( double x0, double y0,
				     const std::string& mcType,
				     double scale = StyleTools::lSS );
    void DrawAtlasInternalMCLeft   ( double x0, double y0,
				     const std::string& mcType,
				     double scale = StyleTools::lSS );
  };
}

#endif
