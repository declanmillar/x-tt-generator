#include "analysis.h"

int main(int argc, char* argv[])
{
  AtlasROOTStyle atlasStyle;
  atlasStyle.SetStyle();

  // TApplication* RootApp = new TApplication("RootApp",&argc,argv);
  
  TString channel("2to2"); 
  TString model("SM");
  TString options("_");
  TString energy("13");
  TString points("1x500000");

  TString base(channel + "_" + model + "_" + energy + options + points);
  TString inputFileName(base + ".root");
  TString outputFileName(base + "_histos.root");
  
  AnalysisZprime analysis(channel, model, inputFileName, outputFileName);

  // AnalysisTransverse SManalysis(TinputFileName1, ToutputFileName1);
  // AnalysisTransverse SSMm2000w200analysis(TinputFileName2, ToutputFileName2);

  // TCanvas* canvas = superpose("mtt", outputFileName1, "mtt", outputFileName2);

  // TFile* superposition = new TFile("superposition.root", "RECREATE");
  // canvas->Write();
  // run root app
  // RootApp->Run( kTRUE );
  // superposition->Close();

  // RootApp->Run( kTRUE );
    
  return 0;  
}
