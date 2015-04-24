#include "analysis_transverse.h"
#include "superpose.h"

int main(int argc, char* argv[])
{
  // set atlas style
  AtlasROOTStyle ars;
  ars.SetStyle();

  // make root app
  // TApplication* RootApp = new TApplication("RootApp",&argc,argv);
  
  std::string inputFileName1("SM_13_2to6_DL_1x5000000.root");
  TString TinputFileName1(inputFileName1);
  std::string outputFileName1("SM.root");
  TString ToutputFileName1(outputFileName1);

  std::string inputFileName2("SSM-m2000-w200_Zp_13_2to6_DL_1x5000000.root");
  TString TinputFileName2(inputFileName2);
  std::string outputFileName2("SSM-m2000-w200.root");
  TString ToutputFileName2(outputFileName2);

  
  AnalysisTransverse SManalysis(TinputFileName1, ToutputFileName1);
  AnalysisTransverse SSMm2000w200analysis(TinputFileName2, ToutputFileName2);

  TCanvas* canvas = superpose("mtt", outputFileName1, "mtt", outputFileName2);

  TFile* superposition = new TFile("superposition.root", "RECREATE");
  canvas->Write();
  // run root app
  // RootApp->Run( kTRUE );
  superposition->Close();
    
  return 0;  
}
