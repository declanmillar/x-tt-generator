#include "analysis.h"

int main(int argc, char* argv[])
{
  // set atlas style
  AtlasROOTStyle ars;
  ars.SetStyle();

  // make root app
  // TApplication* RootApp = new TApplication("RootApp",&argc,argv);
  
  // Select Channel : 2to2 = pp->tt, 2to6 = pp->tt->bbffff
  
  TString channel("2to2"); 
  TString model("PS");
  TString outputFileName(channel + "_" + model + "_R_histos.root");
  
  AnalysisZprime analysis(channel, model, outputFileName);

  // run root app
  // RootApp->Run( kTRUE );
    
  return 0;  
}
