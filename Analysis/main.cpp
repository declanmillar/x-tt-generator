#include "analysis.h"
#include "superpose.h"

int main(int argc, char* argv[])
{
  AtlasROOTStyle atlasStyle;
  atlasStyle.SetStyle();

  TString channel("2to6"); 
  TString model("SM");
  TString options("_");
  TString energy("13");
  TString points("1x10000000");

  TString base(channel + "_" + model + "_" + energy + options + points);
  TString inputFileName(base + ".root");
  TString outputFileName(base + "_histos.root");
  
  AnalysisZprime* analysis = new AnalysisZprime(channel, model, inputFileName, outputFileName);

  TString channel2("2to6"); 
  TString model2("SSM-m2500-w500_Zp");
  TString options2("_");
  TString points2("1x10000000");

  TString base2(channel2 + "_" + model2 + "_" + energy + options2 + points2);
  TString inputFileName2(base2 + ".root");
  TString outputFileName2(base2 + "_histos.root");

  AnalysisZprime* analysis2 = new AnalysisZprime(channel2, model2, inputFileName2, outputFileName2);

  TApplication* RootApp = new TApplication("RootApp",&argc,argv);

  TCanvas* c_HT2 = superpose("HT", outputFileName, "HT", outputFileName2);
  c_HT2->Draw();
  RootApp->Run(kTRUE);
  TCanvas* c_KTbbll2 = superpose("KTbbll", outputFileName, "KTbbll", outputFileName2);
  c_KTbbll2->Draw();

  RootApp->Run(kTRUE);
    
  return 0;  
}
