#include <TROOT.h>
#include <TApplication.h>

#include <iostream>

#include "AtlasStyle.h"
#include "DiJetAnalysisData.h"
#include "DiJetAnalysisMC.h"
#include "DiJetAnalysisBoth.h"

int main(int argc, char *argv[])
{ 
  int mode             = 1;
  int nEvents          = -1;
  int startEvent       = 0;
  bool isData          = false;
  bool is_pPb          = false;
  int  mcType          = 0;    
  bool isReco          = false;
  
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

  if( mode == 3 ){
    if(argc >= 3){
      is_pPb    = atoi( argv[2] );
    } if(argc >= 4){
      isReco    = atoi( argv[3] );
    } 
  }
  
  std::cerr << "argc: " << argc
	    << "   mode: " << mode
	    << "   nEvents: " << nEvents 
	    << "   startEvent: " << startEvent
	    << "   isData: " << isData
	    << "   is_pPb: " << is_pPb
    	    << "   mcType: " << mcType
    	    << "   mcMode: " << isReco
	    << std::endl;

  TApplication*   rootapp = NULL;  
  DiJetAnalysis* analysis = NULL;

  // clean up
  if( mode < 3 ){
    analysis = isData ?
      static_cast< DiJetAnalysis* >
      ( new DiJetAnalysisData( is_pPb ) ) :
      static_cast< DiJetAnalysis* >
      ( new DiJetAnalysisMC  ( is_pPb, mcType ) );
  } else if( mode == 3 ){
    analysis = static_cast< DiJetAnalysis* >
      ( new DiJetAnalysisBoth( is_pPb, isReco ) );
  }
  
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
    analysis->MakeDphiTogether();
    rootapp->Run();
  } else if( mode == 3 ){
    rootapp = new TApplication("JetAnalysis",&argc, argv);
    gROOT->SetBatch(kTRUE);
    analysis->MakeDphiTogether();
    rootapp->Run();
  }
  
  delete analysis;
  return 0;
}
