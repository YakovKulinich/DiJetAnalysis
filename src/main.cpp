#include <TROOT.h>
#include <TApplication.h>

#include <iostream>

#include "AtlasStyle.h"
#include "DiJetAnalysis.h"

int main(int argc, char *argv[])
{ 
  int mode             = 1;
  int nEvents          = -1;
  int startEvent       = 0;
  std::string labelOut = "";
  std::string fNameOut = "myOut.root";

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
    labelOut   = argv[4];
  }

  std::cerr << "argc: " << argc
	    << "    mode: " << mode
	    << "    nEvents: " << nEvents 
	    << "    startEvent: " << startEvent
	    << "    labelOut: " << labelOut 
	    << std::endl;

  TApplication* rootapp = new TApplication("JetAnalysis",&argc, argv);
  SetAtlasStyle();
 
  std::string fNameIn = "/home/yakov/Projects/atlas/data/pPb.root";
  
  DiJetAnalysis* analysis = new DiJetAnalysis();
  if( mode ){ 
      analysis->RunOverTreeFillHistos( 1, nEvents, startEvent, fNameIn, fNameOut ); 
      return 0;
  } else if( !mode ) { 
    analysis->PlotExistingHistos( 1, labelOut, fNameOut ); 
  }
  
  rootapp->Run();

  delete analysis;
  return 0;
}
