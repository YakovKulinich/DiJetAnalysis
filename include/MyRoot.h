#include <TFile.h>
#include <TKey.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TH2D.h>
#include <TH3D.h>
#include <THStack.h>
#include <TF1.h>
#include <TPad.h>
#include <TCanvas.h>
#include <iostream>
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
#include <sstream>

//===================================
//         CONSTANTS
//===================================

namespace constants{
  // barrel range
  const double CETAMAX = 2.8;

  // FCal eta ranges
  const double FETAMIN = 3.2;
  const double FETAMAX = 4.4;

  // trigger eta range
  const double TETAMIN = 3.2;
  const double TETAMAX = 4.8;
  
  // atlas maximum eta
  const double ETAMIN = -5.0;
  const double ETAMAX = 5.0;

  const double DELTAPHIMIN = 2.5;

  const double DELTA = 0.001;

  const double PI = TMath::Pi();

  const double EPSILON = 1e-6;

  const double JETPTCUT = 10;
}

//===================================
//         COMMON FUNCTIONS
//===================================
namespace AnalysisTools{
  bool isForward( const double& eta );
  bool isCentral( const double& eta );
  bool isInTriggerEtaRange( const double& eta );
  bool EpsilonEqual( double a, double b );
  // returns dphi in range 0<dphi<2pi
  double DPhiFC( double phi1, double phi2 );
  double DeltaPhi( double phi1, double phi2 );
  double DeltaR(  const TLorentzVector& jet1, const TLorentzVector& jet2 );
  bool sortByDecendingPt( const TLorentzVector& jet1, const TLorentzVector& jet2 );
  bool TruncateHistoBins( TH3* h3 );
  bool DoPrint( int ev );
  std::vector<std::string> vectorise(TString str, TString sep);
  void FitGaussian( TH1* hProj, TF1* fit );
}

//===================================
//          DRAWING STUFF
//===================================
namespace DrawTools{
  void DrawRightLatex ( double x, double y ,
			const char* s, float scale, int color);
  void DrawLeftLatex  ( double x, double y ,
			const char* s, float scale, int color);
  void DrawCenterLatex( double x, double y ,
			const char* s, float scale, int color);
  void DrawAtlasInternalDataRight( double x0, double y0,
				   double scale, bool is_pPb);
  void DrawAtlasInternalDataLeft ( double x0, double y0,
				   double scale, bool is_pPb );
  void DrawAtlasInternalMCRight  ( double x0, double y0,
				   double scale, bool isReco );
  void DrawAtlasInternalMCLeft   ( double x0, double y0,
				   double scale, bool isReco );
}

namespace StyleTools{
  // Histogram/Graph/Function styles
  void SetCustomMarkerStyle( TH1* his , int iflag );
  void SetCustomMarkerStyle( TGraph* his , int iflag );
  void SetHStyle( TH1*    his, int iflag, float scale);
  void SetHStyle( TGraph* his, int iflag, float scale);
  void SetHStyle( TF1*    his, int iflag, float scale);

  // Legend Style
  void SetLegendStyle(TLegend * legend, float scale);
}
