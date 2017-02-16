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
  const double CETAMAX = 2.8;
  const double FETAMAX = 4.4;
  const double FETAMIN = 3.2;

  // trigger bounds
  const double TETAMIN = 3.2;
  const double TETAMAX = 4.9;

  const double DELTAPHIMIN = 2.5;

  const double DELTA = 0.001;

  const double PI = TMath::Pi();

  const double EPSILON = 1e-6;

  const double JETPTCUT = 10;
}

//===================================
//         COMMON FUNCTIONS
//===================================

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

//===================================
//          DRAWING STUFF
//===================================
void DrawRightLatex( double x, double y , const char* s, float scale, int color);
void DrawLeftLatex( double x, double y , const char* s, float scale, int color);
void DrawCenterLatex( double x, double y , const char* s, float scale, int color);
void DrawAtlasInternalDataRight_pPb( double x0, double y0, double scale );
void DrawAtlasInternalDataLeft_pPb( double x0, double y0, double scale );
void DrawAtlasInternalDataRight_pp( double x0, double y0, double scale );
void DrawAtlasInternalDataLeft_pp( double x0, double y0, double scale );
void DrawAtlasInternalMC( bool isReco, double scale );
void SetCustomMarkerStyle( TH1* his , int iflag );
void SetHStyle( TH1* his, int iflag, float scale);
void SetHStyle( TF1* his, int iflag);
void SetHStyle_open( TH1* his, int iflag);
void SetHStyle_graph( TGraphErrors* his, int iflag);
void SetHStyle_TF1( TF1* his, int iflag);
void SetHStyle_graph( TGraphAsymmErrors* his, int iflag);
void SetHStyle_graph_open( TGraphErrors* his, int iflag);
void SetLegendStyle(TLegend * legend, float scale);
