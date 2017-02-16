//
// ATLAS Style, based on a style file from BaBar
//

#include <iostream>

#include "AtlasStyle.h"

#include "TROOT.h"

#include "TStyle.h"
#include "TColor.h"

void SetAtlasStyle ()
{
  static TStyle* atlasStyle = 0;
  std::cout << "\nApplying ATLAS style settings...\n" << std::endl ;
  if ( atlasStyle==0 ) atlasStyle = AtlasStyle();
  gROOT->SetStyle("ATLAS");
  gROOT->ForceStyle();
}


void Set_Palette_Style(TStyle* style)
{
    const Int_t NRGBs = 5;
    const Int_t NCont = 255;

    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    style->SetNumberContours(NCont);
}


TStyle* AtlasStyle() 
{
   TStyle *atlasStyle = new TStyle("ATLAS","Atlas style");
   
  // use plain black on white colors
  Int_t icol=0; // WHITE
  atlasStyle->SetFrameBorderMode(icol);
  atlasStyle->SetFrameFillColor(icol);
  atlasStyle->SetCanvasBorderMode(icol);
  atlasStyle->SetCanvasColor(icol);
  atlasStyle->SetPadBorderMode(icol);
  atlasStyle->SetPadColor(icol);
  atlasStyle->SetStatColor(icol);
  // atlasStyle->SetFillColor(icol); // don't use: white fill color for *all* objects
  
  // set the paper & margin sizes
  // atlasStyle->SetPaperSize(20,26);

  // set margin sizes
  atlasStyle->SetPadTopMargin(0.075);
  atlasStyle->SetPadRightMargin(0.1);
  atlasStyle->SetPadBottomMargin(0.2);
  atlasStyle->SetPadLeftMargin(0.16);

  // set title offsets (for axis label)
  atlasStyle->SetTitleOffset(1., "xy");

  // use large fonts
  //Int_t font=72; // Helvetica italics
  Int_t font=42; // Helvetica
  Double_t tsize=0.07;
  atlasStyle->SetTextFont(font);
  atlasStyle->SetTextSize(tsize);

  atlasStyle->SetLabelFont(font,"xyz");
  atlasStyle->SetLabelSize(tsize,"xyz");

  atlasStyle->SetTitleFont(font,"xyz");
  atlasStyle->SetTitleSize(tsize,"xyz");

  // use bold lines and markers
  atlasStyle->SetMarkerStyle(20);
  atlasStyle->SetMarkerSize(1);
  atlasStyle->SetHistLineWidth(1);
  atlasStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes

  // get rid of X error bars 
  // atlasStyle->SetErrorX(0.001);
  // get rid of error bar caps
  atlasStyle->SetEndErrorSize(0.);

  // title things
  atlasStyle->SetTitleStyle(0);
  atlasStyle->SetTitleBorderSize(0);
  atlasStyle->SetLabelFont(font,"t");
  atlasStyle->SetTitleFont(font,"t");
  atlasStyle->SetLabelSize(tsize, "t");
  atlasStyle->SetTitleSize(tsize, "t");
  atlasStyle->SetTitleX(0.1f);
  atlasStyle->SetTitleW(0.8f);

  // label things
  atlasStyle->SetLegendBorderSize(0);

  // do not display any of the standard histogram decorations
  atlasStyle->SetOptTitle(0);
  atlasStyle->SetOptStat(0);
  atlasStyle->SetOptFit(0);

  // put tick marks on top and RHS of plots
  atlasStyle->SetPadTickX(1);
  atlasStyle->SetPadTickY(1);

  Set_Palette_Style(atlasStyle);

  return atlasStyle;

}
