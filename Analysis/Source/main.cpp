#include "analysis.h"
#include "overlay.h"
#include "superpose2d.h"

int main(int argc, char* argv[])
{
  AtlasROOTStyle atlasStyle;
  atlasStyle.SetStyle();

  // TString inDir("/afs/cern.ch/work/d/demillar/Ntuples_Zprime/");
  TString inDir("/Users/declan/Data/Ntuples_Zprime/");
  TString outDir("/Users/declan/Data/Histos_Zprime/");

  TString channel("2to6"); 
  TString model("SM");
  TString options("_xc_");
  TString energy("13");
  TString points("5x1000000");
  TString base(channel + "_" + model + "_" + energy + options + points);
  TString inputFileName(inDir + "/" + channel + "/" + base + ".root");
  TString outputFileName(outDir + base + "_histos.root");  

  TString channel2("2to6"); 
  TString model2("GSM-SM");
  TString options2("_xc_");
  TString points2("5x1000000");
  TString base2(channel2 + "_" + model2 + "_" + energy + options2 + points2);
  TString inputFileName2(inDir + "/" + channel + "/" + base2 + ".root");
  TString outputFileName2(outDir + base2 + "_histos.root");

  TString channel3("2to6"); 
  TString model3("GLR-R");
  TString options3("_xc_");
  TString points3("5x1000000");
  TString base3(channel3 + "_" + model3 + "_" + energy + options3 + points3);
  TString inputFileName3(inDir + "/" + channel + "/" + base3 + ".root");
  TString outputFileName3(outDir + base3 + "_histos.root");

  TString channel4("2to6"); 
  TString model4("E6-eta");
  TString options4("_xc_");
  TString points4("5x1000000");
  TString base4(channel4 + "_" + model4 + "_" + energy + options4 + points4);
  TString inputFileName4(inDir + "/" + channel + "/" + base4 + ".root");
  TString outputFileName4(outDir + base4 + "_histos.root");

  // AnalysisZprime analysis(channel, model, inputFileName, outputFileName);
  // AnalysisZprime analysis2(channel2, model2, inputFileName2, outputFileName2);
  // AnalysisZprime analysis3(channel3, model3, inputFileName3, outputFileName3);
  // AnalysisZprime analysis4(channel4, model4, inputFileName4, outputFileName4);

  bool overlayTransverse = false;

  TApplication* RootApp = new TApplication("RootApp",&argc,argv);

    // TCanvas* c_AttC =overlay(false, false, "AttC", model, outputFileName, 
    //                                    "AttC", model2, outputFileName2,
    //                                    "AFBstar", model3, outputFileName3,
    //                                    "AFBstar", model4, outputFileName4);
    // c_AttC->Draw();
    // RootApp->Run(kTRUE);

    TCanvas* c_AFBstar =overlay(false, true, "AFBstar", model, outputFileName, 
                                          // "AFBstar", model2, outputFileName2,)
                                          "AFBstar", model3, outputFileName3);
                                          // "AFBstar", model4, outputFileName4);
    c_AFBstar->Draw();
    RootApp->Run(kTRUE);


  if (channel == "2to2") {

    TCanvas* c_Mtt =overlay(false, false, "Mtt", model, outputFileName, 
                                      "Mtt", model2, outputFileName2,
                                      "Mtt", model3, outputFileName3,
                                      "Mtt", model4, outputFileName4
                                      );
    c_Mtt->Draw();
    RootApp->Run(kTRUE);

    TCanvas* c_AL =overlay(false, false,  "AL", model, outputFileName, 
                                      "AL", model2, outputFileName2,
                                      "AL", model3, outputFileName3,
                                      "AL", model4, outputFileName4);
    c_AL->Draw();
    RootApp->Run(kTRUE);


    TCanvas* c_ALL =overlay(false, false, "ALL", model, outputFileName, 
                                      "ALL", model2, outputFileName2,
                                      "ALL", model3, outputFileName3,
                                      "ALL", model4, outputFileName4);
    c_ALL->Draw();
    RootApp->Run(kTRUE);
  }

  if (channel == "2to6") {
    if (overlayTransverse == true) {
      TCanvas* c_KTbbll =overlay(true, "KTbbll", "SM", outputFileName, 
                                    "KTbbll", "(Z',SM) (M = 2.5TeV, #Gamma = 0.5 TeV)", outputFileName2);
      c_KTbbll->Draw();
      RootApp->Run(kTRUE);

      TCanvas* c_HT =overlay(true, "HT", "SM", outputFileName, 
                                "HT", "(Z',SM) (M = 2.5TeV, #Gamma = 0.5 TeV)", outputFileName2);
      c_HT->Draw();
      RootApp->Run(kTRUE);

      TCanvas* c_dphi_HT =superpose2d(true, "dphi_HT", "|M_{SM}|^{2}", outputFileName, 
                                "dphi_HT", "|M_{Z'}|^{2} (M = 2.5TeV, #Gamma = 0.5 TeV)", outputFileName2);
      c_dphi_HT->Draw();

      RootApp->Run(kTRUE);

      TCanvas* c_dphi_KTbbll =superpose2d(true, "dphi_KTbbll", "|M_{SM}|^{2}", outputFileName, 
                                "dphi_KTbbll", "|M_{Z'}|^{2} (M = 2.5TeV, #Gamma = 0.5 TeV)", outputFileName2);
      c_dphi_KTbbll->Draw();

      RootApp->Run(kTRUE);
    }

    TCanvas* c_AFBstarNuReco =overlay(false, true, "AFBstarNuReco", model, outputFileName, 
                                                // "AFBstarNuReco", model2, outputFileName2,
                                                "AFBstarNuReco", model3, outputFileName3);
                                                // "AFBstarNuReco", model4, outputFileName4);
    c_AFBstarNuReco->Draw();
    RootApp->Run(kTRUE);

    TCanvas* c_PzNu =overlay(false, false, "PzNu", "truth", outputFileName, 
                          "PzNuReco", "reconstructed", outputFileName);
    c_PzNu->Draw();
    RootApp->Run(kTRUE);

    TCanvas* c_MttReco =overlay(false, false, "Mtt", "truth", outputFileName, 
                              "MttReco", "reconstructed", outputFileName);
    c_MttReco->Draw();
    RootApp->Run(kTRUE);

    TCanvas* c_yttReco =overlay(false, false, "ytt", "truth", outputFileName, 
                              "yttReco", "reconstructed", outputFileName);
    c_yttReco->Draw();
    RootApp->Run(kTRUE);

    TCanvas* c_CosThetaReco =overlay(false, false, "CosTheta", "truth", outputFileName, 
                              "CosThetaReco", "reconstructed", outputFileName);
    c_CosThetaReco->Draw();
    RootApp->Run(kTRUE);

    TCanvas* c_CosThetaStarReco =overlay(false, false, "CosThetaStar", "truth", outputFileName, 
                              "CosThetaStarReco", "reconstructed", outputFileName);
    c_CosThetaStarReco->Draw();
    RootApp->Run(kTRUE);
  }

  return 0;  
}
