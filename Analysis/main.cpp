#include "analysis.h"
#include "superpose.h"

int main(int argc, char* argv[])
{
  AtlasROOTStyle atlasStyle;
  atlasStyle.SetStyle();

  TString channel("2to2"); 
  TString model("PS");
  TString options("_");
  TString energy("13");
  TString points("1x10000000");

  TString base(channel + "_" + model + "_" + energy + options + points);
  TString inputFileName(base + ".root");
  TString outputFileName(base + "_histos.root");
  
  AnalysisZprime* analysis = new AnalysisZprime(channel, model, inputFileName, outputFileName);

  TString channel2("2to2"); 
  TString model2("PS");
  TString options2("_R_");
  TString points2("1x10000000");

  TString base2(channel2 + "_" + model2 + "_" + energy + options2 + points2);
  TString inputFileName2(base2 + ".root");
  TString outputFileName2(base2 + "_histos.root");

  AnalysisZprime* analysis2 = new AnalysisZprime(channel2, model2, inputFileName2, outputFileName2);

  TApplication* RootApp = new TApplication("RootApp",&argc,argv);

  TCanvas* c_HT = superpose("Mtt", outputFileName, "Mtt", outputFileName2);
  RootApp->Run(kTRUE);

  RootApp->Run(kTRUE);
    
  return 0;  
}