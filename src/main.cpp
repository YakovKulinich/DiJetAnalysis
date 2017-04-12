#include <TROOT.h>
#include <TApplication.h>

#include <iostream>

#include "AtlasStyle.h"
#include "DiJetAnalysisData.h"
#include "DiJetAnalysisMC.h"

int main(int argc, char *argv[])
{ 
  int mode             = 1;
  int nEvents          = -1;
  int startEvent       = 0;
  bool isData          = true;
  bool is_pPb          = true;
  int  mcType          = 0;    // pythia8+powheg

  if(argc >= 2){
    mode       = atoi( argv[1] );
  } if(argc >= 3){
    nEvents    = atoi( argv[2] );
  } if(argc >= 4){
    startEvent = atoi( argv[3] );
  } if(argc >= 5){
    isData     = atoi( argv[4] );
  } if(argc >= 6){
    is_pPb     = atoi( argv[5] );
  } if(argc >= 7){
    mcType     = atoi( argv[6] );
  }

  std::cerr << "argc: " << argc
	    << "    mode: " << mode
	    << "    nEvents: " << nEvents 
	    << "    startEvent: " << startEvent
	    << "    isData: " << isData
	    << "    is_pPb: " << is_pPb
    	    << "    mcType: " << mcType
	    << std::endl;

  TApplication* rootapp = NULL;
  
  DiJetAnalysis* analysis = isData ?
    static_cast< DiJetAnalysis* >
    ( new DiJetAnalysisData( isData, is_pPb ) ) :
    static_cast< DiJetAnalysis* >
    ( new DiJetAnalysisMC  ( isData, is_pPb, mcType ) );

  analysis->Initialize();

  if( mode == 1){
    analysis->RunOverTreeFillHistos( nEvents, startEvent ); 
  } else if( mode == 0 ) {
    rootapp = new TApplication("JetAnalysis",&argc, argv);
    gROOT->SetBatch(kTRUE);
    analysis->ProcessPlotHistos();
    rootapp->Run();
  } else if( mode == 2 ){
    rootapp = new TApplication("JetAnalysis",&argc, argv);
    gROOT->SetBatch(kTRUE);
    analysis->PlotDataTogether();
    rootapp->Run();
  }
  
  delete analysis;
  return 0;
}
