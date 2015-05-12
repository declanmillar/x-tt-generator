#include "analysis.h"
#include "superpose.h"
#include "superpose2d.h"

int main(int argc, char* argv[])
{
  AtlasROOTStyle atlasStyle;
  atlasStyle.SetStyle();

  TString inDir("/afs/cern.ch/work/d/demillar/Ntuples_Zprime/");
  // TString inDir("/Users/declan/Data/Ntuples_Zprime/");
  TString outDir("Histos/");
  TString channel("2to6"); 
  TString model("SM");
  TString options("_");
  TString energy("13");
  TString points("5x1000000");

  TString base(channel + "_" + model + "_" + energy + options + points);
  TString inputFileName(inDir + base + ".root");
  TString outputFileName(base + "_histos.root");
  
  AnalysisZprime analysis(channel, model, inputFileName, outputFileName);

  TString channel2("2to6"); 
  TString model2("SSM-m2500-w500_Zp");
  TString options2("_");
  TString points2("5x1500000");

  TString base2(channel2 + "_" + model2 + "_" + energy + options2 + points2);
  TString inputFileName2(inDir + base2 + ".root");
  TString outputFileName2(base2 + "_histos.root");

  AnalysisZprime analysis2(channel2, model2, inputFileName2, outputFileName2);

  TApplication* RootApp = new TApplication("RootApp",&argc,argv);

  if (channel == "2to2") {
  TCanvas* c_Mtt = superpose("Mtt", "|M|^{2}=|M_{Z'}|^{2}+2R(M_{#gamma}^{*}M_{Z'})+2R(M^{*}_{Z}M_{Z'})", outputFileName, 
                            "Mtt", "|M|^{2} = |M_{Z'}|^{2}", outputFileName2);
    c_Mtt->Draw();
    RootApp->Run(kTRUE);
  }

  if (channel == "2to6") {
    // TCanvas* c_KTbbll = superpose("KTbbll", "|M_{SM}|^{2}", outputFileName, 
    //                               "KTbbll", "|M_{Z'}|^{2} (M = 2.5TeV, #Gamma = 0.5 TeV)", outputFileName2);
    // c_KTbbll->Draw();

    // RootApp->Run(kTRUE);

    // TCanvas* c_HT = superpose("HT", "|M_{SM}|^{2}", outputFileName, 
    //                           "HT", "|M_{Z'}|^{2} (M = 2.5TeV, #Gamma = 0.5 TeV)", outputFileName2);
    // c_HT->Draw();

    // RootApp->Run(kTRUE);

    // TCanvas* c_dphi_HT = superpose2d("dphi_HT", "|M_{SM}|^{2}", outputFileName, 
    //                           "dphi_HT", "|M_{Z'}|^{2} (M = 2.5TeV, #Gamma = 0.5 TeV)", outputFileName2);
    // c_dphi_HT->Draw();

    // RootApp->Run(kTRUE);

    TCanvas* c_dphi_KTbbll = superpose2d("dphi_KTbbll", "|M_{SM}|^{2}", outputFileName, 
                              "dphi_KTbbll", "|M_{Z'}|^{2} (M = 2.5TeV, #Gamma = 0.5 TeV)", outputFileName2);
    c_dphi_KTbbll->Draw();

    RootApp->Run(kTRUE);
  }

  return 0;  
}
