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
  TString channel("2to2"); 
  TString model("SM");
  TString options("_xc_");
  TString energy("13");
  TString points("5x500000");

  TString base(channel + "_" + model + "_" + energy + options + points);
  TString inputFileName(inDir + base + ".root");
  TString outputFileName(base + "_histos.root");
  
  // AnalysisZprime analysis(channel, model, inputFileName, outputFileName);

  TString channel2("2to2"); 
  TString model2("GSM-SM");
  TString options2("_xc_");
  TString points2("5x500000");

  TString base2(channel2 + "_" + model2 + "_" + energy + options2 + points2);
  TString inputFileName2(inDir + base2 + ".root");
  TString outputFileName2(base2 + "_histos.root");

  // AnalysisZprime analysis2(channel2, model2, inputFileName2, outputFileName2);

  TApplication* RootApp = new TApplication("RootApp",&argc,argv);

  TCanvas* c_Mtt = superpose(false, "Mtt", "|M|^{2}=|M_{Z'}|^{2}+2R(M_{#gamma}^{*}M_{Z'})+2R(M^{*}_{Z}M_{Z'})", outputFileName, 
                            "Mtt", "|M|^{2} = |M_{SM}|^{2}", outputFileName2);
  c_Mtt->Draw();
  RootApp->Run(kTRUE);

  if (channel == "2to2") {
    TCanvas* c_ALL = superpose(false,"ALL", "SM", outputFileName, 
                               "ALL", "GSM-SM", outputFileName2);
    c_ALL->Draw();
    RootApp->Run(kTRUE);

    TCanvas* c_AL = superpose(false,"AL", "SM", outputFileName, 
                               "AL", "GSM-SM", outputFileName2);
    c_AL->Draw();
    RootApp->Run(kTRUE);

    TCanvas* c_AFBstar = superpose(false,"AFBstar", "SM", outputFileName, 
                               "AFBstar", "GSM-SM", outputFileName2);
    c_AFBstar->Draw();
    RootApp->Run(kTRUE);

  }


  if (channel == "2to6") {
    TCanvas* c_KTbbll = superpose(true, "KTbbll", "SM", outputFileName, 
                                  "KTbbll", "(Z',SM) (M = 2.5TeV, #Gamma = 0.5 TeV)", outputFileName2);
    c_KTbbll->Draw();
    RootApp->Run(kTRUE);

    TCanvas* c_HT = superpose(true, "HT", "SM", outputFileName, 
                              "HT", "(Z',SM) (M = 2.5TeV, #Gamma = 0.5 TeV)", outputFileName2);
    c_HT->Draw();
    RootApp->Run(kTRUE);

    TCanvas* c_dphi_HT = superpose2d(true, "dphi_HT", "|M_{SM}|^{2}", outputFileName, 
                              "dphi_HT", "|M_{Z'}|^{2} (M = 2.5TeV, #Gamma = 0.5 TeV)", outputFileName2);
    c_dphi_HT->Draw();

    RootApp->Run(kTRUE);

    TCanvas* c_dphi_KTbbll = superpose2d(true, "dphi_KTbbll", "|M_{SM}|^{2}", outputFileName, 
                              "dphi_KTbbll", "|M_{Z'}|^{2} (M = 2.5TeV, #Gamma = 0.5 TeV)", outputFileName2);
    c_dphi_KTbbll->Draw();

    RootApp->Run(kTRUE);
  }

  return 0;  
}
