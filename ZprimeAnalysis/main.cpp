#include "analysis.h"

int main(int argc, char* argv[])
{
  // set atlas style
  AtlasROOTStyle ars;
  ars.SetStyle();

  // make root app
  // TApplication* RootApp = new TApplication("RootApp",&argc,argv);
  
  // Select Channel : el = electron, mu = muon
  
  // TString channel("el"); 
  TString outputFileName("histos.root");
  
  AnalysisZprime analysis(outputFileName);

  // run root app
  // RootApp->Run( kTRUE );
    
  return 0;  
}
