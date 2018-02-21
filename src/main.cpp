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
  bool is_pPb          = false;
  bool isData          = false;
  int  mcType          = 0;    
  int  uncertComp      = 0;
  
  bool isReco          = false;
  
  if(argc >= 2){
    mode       = atoi( argv[1] );
  } if(argc >= 3){
    nEvents    = atoi( argv[2] );
  } if(argc >= 4){
    startEvent = atoi( argv[3] );
  } if(argc >= 5){
    is_pPb     = atoi( argv[4] );
  } if(argc >= 6){
    isData     = atoi( argv[5] );
  } if(argc >= 7){
    mcType     = atoi( argv[6] );
  }  if(argc >= 8){
    uncertComp = atoi( argv[7] );
  }

  // reserve mode 8 for plotting MC and Data together
  if( mode == 8 ){
    if(argc >= 3){
      is_pPb    = atoi( argv[2] );
    } if(argc >= 4 ){
      isReco    = atoi( argv[3] );
    } 
  }
  
  std::cerr << "   argc: " << argc
	    << "   mode: " << mode
	    << "   nEvents: " << nEvents 
	    << "   startEvent: " << startEvent
	    << "   is_pPb: " << is_pPb
    	    << "   isData: " << isData
	    << "   mcType: " << mcType
	    << "   uncertComp: " << uncertComp
	    << std::endl;

  DiJetAnalysis* analysis = NULL;

  SetAtlasStyle();
  
  // reserve mode 8 for plotting MC and Data together
  if( mode < 8 ){
    analysis = isData ?
      static_cast< DiJetAnalysis* >
      ( new DiJetAnalysisData( is_pPb, mcType, uncertComp ) ) :
      static_cast< DiJetAnalysis* >
      ( new DiJetAnalysisMC  ( is_pPb, mcType, uncertComp ) );
  } else if( mode == 8 ){
    analysis = static_cast< DiJetAnalysis* >
      ( new DiJetAnalysisBoth( is_pPb, isReco ) );
  }

  gROOT->SetBatch(kTRUE);
  
  analysis->Initialize();

  switch( mode ){
  case -1:
    analysis->RunOverTreeMakeWeights( nEvents, startEvent ); 
    break;
  case -2:
    analysis->ProcessWeights();
    break;
  case 0:
    analysis->RunOverTreeFillHistos( nEvents, startEvent ); 
    break;
  case 1:
    analysis->ProcessPerformance();
    break;
  case 2:
    analysis->UnfoldPerformance();
    break;
  case 3:
    analysis->ProcessPhysics();
    break;
  case 4:
    analysis->UnfoldPhysics();
    break;
  case 5:
    analysis->MakeResultsTogether();
    break;
  case 6:
    analysis->ProcessSystematics();
    break;
  case 8:
    analysis->MakeResultsTogether();
  default:
    std::cout << "Wrong Mode. Exiting." << std::endl;
    return 1;
  }
  
  delete analysis;
  return 0;
}
