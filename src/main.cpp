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

  if(argc == 2){
    mode       = atoi( argv[1] );
  } else if(argc == 3){
    mode       = atoi( argv[1] );
    nEvents    = atoi( argv[2] );
  } else if(argc == 4){
    mode       = atoi( argv[1] );
    nEvents    = atoi( argv[2] );
    startEvent = atoi( argv[3] );
  } else if(argc == 5){
    mode       = atoi( argv[1] );
    nEvents    = atoi( argv[2] );
    startEvent = atoi( argv[3] );
    isData     = atoi( argv[4] );
  } else if(argc == 6){
    mode       = atoi( argv[1] );
    nEvents    = atoi( argv[2] );
    startEvent = atoi( argv[3] );
    isData     = atoi( argv[4] );
    is_pPb     = atoi( argv[5] );
  }

  std::cerr << "argc: " << argc
	    << "    mode: " << mode
	    << "    nEvents: " << nEvents 
	    << "    startEvent: " << startEvent
	    << "    isData: " << isData
	    << "    is_pPb: " << is_pPb
	    << std::endl;

  TApplication* rootapp = NULL;
  
  DiJetAnalysis* analysis = isData ?
    static_cast< DiJetAnalysis* >
    ( new DiJetAnalysisData( isData, is_pPb ) ) :
    static_cast< DiJetAnalysis* >
    ( new DiJetAnalysisMC  ( isData, is_pPb ) );

  analysis->Initialize();

  if( mode ){
    analysis->RunOverTreeFillHistos( nEvents, startEvent ); 
  } else if( !mode ) {
    rootapp = new TApplication("JetAnalysis",&argc, argv);
    gROOT->SetBatch(kTRUE);
    analysis->ProcessPlotHistos();
    rootapp->Run();
  }
  
  delete analysis;
  return 0;
}
