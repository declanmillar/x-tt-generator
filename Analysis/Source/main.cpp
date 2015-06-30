#include "analysis.h"
#include "overlay.h"
#include "superpose2d.h"

int main(int argc, char* argv[])
{
  AtlasROOTStyle atlasStyle;
  atlasStyle.SetStyle();

  TString dataDirectory;
  #if __linux || __linux__
    dataDirectory = "/afs/cern.ch/work/d/demillar/Zp-tt_pheno";
  #elif __APPLE__ || __MACH__
    dataDirectory = "/Users/declan/Data/Zp-tt_pheno";
  #endif

  TString ntupleDirectory = dataDirectory + "/NTuples";
  TString weightsDirectory = dataDirectory + "/Weights";
  TString histogramDirectory = dataDirectory + "/Histograms";

  TString channel("tt"); 
  TString model("SM");
  TString options("_xc_");
  TString energy("13");
  TString points("5x50000");
  TString filename(channel + "_" + model + "_" + energy + options + points);
  TString inputFilename(ntupleDirectory + "/" + channel + "/" + filename + ".root");
  TString weightFilename(weightsDirectory + "/" + channel + "/" + filename + ".txt");
  TString outputFilename(histogramDirectory + "/" + channel + "/" + filename + "_hist.root");

  TString channel2("tt"); 
  TString model2("NUSU2-40");
  TString options2("_xc_");
  TString points2("5x50000");
  TString filename2(channel2 + "_" + model2 + "_" + energy + options2 + points2);
  TString inputFilename2(ntupleDirectory + "/" + channel2 + "/" + filename2 + ".root");
  TString weightFilename2(weightsDirectory + "/" + channel + "/" + filename2 + ".txt");
  TString outputFilename2(histogramDirectory + "/" + channel + "/" + filename2 + "_hist.root");

  // TString channel3("2to6"); 
  // TString model3("GLR-R");
  // TString options3("_xc_");
  // TString points3("5x1000000");
  // TString filename3(channel3 + "_" + model3 + "_" + energy + options3 + points3);
  // TString inputFilename3(ntupleDirectory + "/" + channel3 + "/" + filename3);
  // TString outputFilename3(histogramDirectory + filename3 + "_hist.root");

  // TString channel4("2to6"); 
  // TString model4("E6-eta");
  // TString options4("_xc_");
  // TString points4("5x5000000");
  // TString filename4(channel4 + "_" + model4 + "_" + energy + options4 + points4);
  // TString inputFilename4(ntupleDirectory + "/" + channel4 + "/" + filename4);
  // TString outputFilename4(histogramDirectory + "/"  filename4 + "_hist.root");

  // AnalysisZprime analysis(channel, model, inputFilename, weightFilename, outputFilename);
  // AnalysisZprime analysis2(channel2, model2, inputFilename2, weightFilename2, outputFilename2);
  // AnalysisZprime analysis3(channel3, model3, inputFilename3, outputFilename3);
  // AnalysisZprime analysis4(channel4, model4, inputFilename4, outputFilename4);

  // bool overlayTransverse = false;

  TApplication* RootApp = new TApplication("RootApp",&argc,argv);

  if (channel == "ll" or channel == "tt") {
    TCanvas* c_Mtt =overlay(false, false, "Mtt", model, outputFilename,
                                      "Mtt", model2, outputFilename2
                                      // "Mtt", model3, outputFilename3,
                                      // "Mtt", model4, outputFilename4
                                      );
    c_Mtt->Draw();
    RootApp->Run(kTRUE);
  }

  //   // TCanvas* c_AttC =overlay(false, false, "AttC", model, outputFilename, 
  //   //                                    "AttC", model2, outputFilename2,
  //   //                                    "AFBstar", model3, outputFilename3,
  //   //                                    "AFBstar", model4, outputFilename4);
  //   // c_AttC->Draw();
  //   // RootApp->Run(kTRUE);

  if (channel == "tt") {


    TCanvas* c_AFBstar =overlay(false, false, "AFBstar", model, outputFilename, 
                                              "AFBstar", model2, outputFilename2);
                                          // "AFBstar", model3, outputFilename3);
                                          // "AFBstar", model4, outputFilename4);
    c_AFBstar->Draw();
    RootApp->Run(kTRUE);
  }


  // if (channel == "2to2") {

  //   TCanvas* c_Mtt =overlay(false, false, "Mtt", model, outputFilename, 
  //                                     "Mtt", model2, outputFilename2,
  //                                     "Mtt", model3, outputFilename3,
  //                                     "Mtt", model4, outputFilename4
  //                                     );
  //   c_Mtt->Draw();
  //   RootApp->Run(kTRUE);

  //   TCanvas* c_AL =overlay(false, false,  "AL", model, outputFilename, 
  //                                     "AL", model2, outputFilename2,
  //                                     "AL", model3, outputFilename3,
  //                                     "AL", model4, outputFilename4);
  //   c_AL->Draw();
  //   RootApp->Run(kTRUE);


  //   TCanvas* c_ALL =overlay(false, false, "ALL", model, outputFilename, 
  //                                     "ALL", model2, outputFilename2,
  //                                     "ALL", model3, outputFilename3,
  //                                     "ALL", model4, outputFilename4);
  //   c_ALL->Draw();
  //   RootApp->Run(kTRUE);
  // }

  // if (channel == "2to6") {
  //   if (overlayTransverse == true) {
  //     TCanvas* c_KTbbll =overlay(true, "KTbbll", "SM", outputFilename, 
  //                                   "KTbbll", "(Z',SM) (M = 2.5TeV, #Gamma = 0.5 TeV)", outputFilename2);
  //     c_KTbbll->Draw();
  //     RootApp->Run(kTRUE);

  //     TCanvas* c_HT =overlay(true, "HT", "SM", outputFilename, 
  //                               "HT", "(Z',SM) (M = 2.5TeV, #Gamma = 0.5 TeV)", outputFilename2);
  //     c_HT->Draw();
  //     RootApp->Run(kTRUE);

  //     TCanvas* c_dphi_HT =superpose2d(true, "dphi_HT", "|M_{SM}|^{2}", outputFilename, 
  //                               "dphi_HT", "|M_{Z'}|^{2} (M = 2.5TeV, #Gamma = 0.5 TeV)", outputFilename2);
  //     c_dphi_HT->Draw();

  //     RootApp->Run(kTRUE);

  //     TCanvas* c_dphi_KTbbll =superpose2d(true, "dphi_KTbbll", "|M_{SM}|^{2}", outputFilename, 
  //                               "dphi_KTbbll", "|M_{Z'}|^{2} (M = 2.5TeV, #Gamma = 0.5 TeV)", outputFilename2);
  //     c_dphi_KTbbll->Draw();

  //     RootApp->Run(kTRUE);
  //   }

  //   TCanvas* c_AFBstarNuReco =overlay(false, false, "AFBstarNuReco", model, outputFilename, 
  //                                               "AFBstarNuReco", model2, outputFilename2);
  //                                               // "AFBstarNuReco", model3, outputFilename3);
  //                                               // "AFBstarNuReco", model4, outputFilename4);
  //   c_AFBstarNuReco->Draw();
  //   RootApp->Run(kTRUE);

  //   TCanvas* c_PzNu =overlay(false, false, "PzNu", "truth", outputFilename, 
  //                         "PzNuReco", "reconstructed", outputFilename);
  //   c_PzNu->Draw();
  //   RootApp->Run(kTRUE);

  //   TCanvas* c_MttReco =overlay(false, false, "Mtt", "truth", outputFilename, 
  //                             "MttReco", "reconstructed", outputFilename);
  //   c_MttReco->Draw();
  //   RootApp->Run(kTRUE);

  //   TCanvas* c_yttReco =overlay(false, false, "ytt", "truth", outputFilename, 
  //                             "yttReco", "reconstructed", outputFilename);
  //   c_yttReco->Draw();
  //   RootApp->Run(kTRUE);

  //   TCanvas* c_CosThetaReco =overlay(false, false, "CosTheta", "truth", outputFilename, 
  //                             "CosThetaReco", "reconstructed", outputFilename);
  //   c_CosThetaReco->Draw();
  //   RootApp->Run(kTRUE);

  //   TCanvas* c_CosThetaStarReco =overlay(false, false, "CosThetaStar", "truth", outputFilename, 
  //                             "CosThetaStarReco", "reconstructed", outputFilename);
  //   c_CosThetaStarReco->Draw();
  //   RootApp->Run(kTRUE);
  // }

  return 0;  
}
