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
  TString points("5x1000000");
  TString base(channel + "_" + model + "_" + energy + options + points);
  TString inputFileName(inDir + base + ".root");
  TString outputFileName(outDir + base + "_histos_test.root");  

  TString channel2("2to2"); 
  TString model2("GSM-SM");
  TString options2("_xc_");
  TString points2("5x1000000");
  TString base2(channel2 + "_" + model2 + "_" + energy + options2 + points2);
  TString inputFileName2(inDir + base2 + ".root");
  TString outputFileName2(outDir + base2 + "_histos.root");

  TString channel3("2to2"); 
  TString model3("GSM-Q");
  TString options3("_xc_");
  TString points3("5x1000000");
  TString base3(channel3 + "_" + model3 + "_" + energy + options3 + points3);
  TString inputFileName3(inDir + base3 + ".root");
  TString outputFileName3(outDir + base3 + "_histos.root");

  TString channel4("2to2"); 
  TString model4("GSM-T3L");
  TString options4("_xc_");
  TString points4("5x1000000");
  TString base4(channel4 + "_" + model4 + "_" + energy + options4 + points4);
  TString inputFileName4(inDir + base4 + ".root");
  TString outputFileName4(outDir + base4 + "_histos.root");

  // AnalysisZprime analysis(channel, model, inputFileName, outputFileName);
  // AnalysisZprime analysis2(channel2, model2, inputFileName2, outputFileName2);
  // AnalysisZprime analysis3(channel3, model3, inputFileName3, outputFileName3);
  AnalysisZprime analysis4(channel4, model4, inputFileName4, outputFileName4);

  TApplication* RootApp = new TApplication("RootApp",&argc,argv);


  if (channel == "2to2") {
    TCanvas* c_Mtt = superpose(false, "Mtt", model, outputFileName, 
                                      "Mtt", model2, outputFileName2,
                                      "Mtt", model3, outputFileName3,
                                      "Mtt", model4, outputFileName4);
    c_Mtt->Draw();
    RootApp->Run(kTRUE);


    TCanvas* c_ALL = superpose(false, "ALL", model, outputFileName, 
                                      "ALL", model2, outputFileName2,
                                      "ALL", model3, outputFileName3,
                                      "ALL", model4, outputFileName4);
    c_ALL->Draw();
    RootApp->Run(kTRUE);

    TCanvas* c_AL = superpose(false, "AL", model, outputFileName, 
                                     "AL", model2, outputFileName2,
                                     "AL", model3, outputFileName3,
                                     "AL", model4, outputFileName4);
    c_AL->Draw();
    RootApp->Run(kTRUE);

    TCanvas* c_AFBstar = superpose(false, "AFBstar", model, outputFileName, 
                                          "AFBstar", model2, outputFileName2,
                                          "AFBstar", model3, outputFileName3,
                                          "AFBstar", model4, outputFileName4);
    c_AFBstar->Draw();
    RootApp->Run(kTRUE);

    TCanvas* c_ARFB = superpose(false, "ARFB", model, outputFileName, 
                                       "ARFB", model2, outputFileName2,
                                       "ARFB", model3, outputFileName3,
                                       "ARFB", model4, outputFileName4);
    c_ARFB->Draw();
    RootApp->Run(kTRUE);

  }


  if (channel == "2to6") {
    // TCanvas* c_KTbbll = superpose(true, "KTbbll", "SM", outputFileName, 
    //                               "KTbbll", "(Z',SM) (M = 2.5TeV, #Gamma = 0.5 TeV)", outputFileName2);
    // c_KTbbll->Draw();
    // RootApp->Run(kTRUE);

    // TCanvas* c_HT = superpose(true, "HT", "SM", outputFileName, 
    //                           "HT", "(Z',SM) (M = 2.5TeV, #Gamma = 0.5 TeV)", outputFileName2);
    // c_HT->Draw();
    // RootApp->Run(kTRUE);

    // TCanvas* c_dphi_HT = superpose2d(true, "dphi_HT", "|M_{SM}|^{2}", outputFileName, 
    //                           "dphi_HT", "|M_{Z'}|^{2} (M = 2.5TeV, #Gamma = 0.5 TeV)", outputFileName2);
    // c_dphi_HT->Draw();

    // RootApp->Run(kTRUE);

    // TCanvas* c_dphi_KTbbll = superpose2d(true, "dphi_KTbbll", "|M_{SM}|^{2}", outputFileName, 
    //                           "dphi_KTbbll", "|M_{Z'}|^{2} (M = 2.5TeV, #Gamma = 0.5 TeV)", outputFileName2);
    // c_dphi_KTbbll->Draw();

    // RootApp->Run(kTRUE);

    TCanvas* c_PzNu = superpose2d(true, "PzNu", "NULL", outputFileName, 
                          "PzNuReco", "NULL", outputFileName2);
		c_PzNu->Draw();
		RootApp->Run(kTRUE);

		TCanvas* c_dphi_KTbbll = superpose2d(true, "MttPzNuReco", "NULL", outputFileName, 
		                          "MttPzNuReco", "NULL", outputFileName2);
		c_dphi_KTbbll->Draw();
		RootApp->Run(kTRUE);
  }

  return 0;  
}
