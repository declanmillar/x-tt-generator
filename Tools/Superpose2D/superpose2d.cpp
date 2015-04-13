// superpose
// Declan Millar

#include <iostream>
#include <fstream>
#include <string>
#include <boost/program_options.hpp>

#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TApplication.h>
#include <TAxis.h>
#include <TH2.h>
#include <TFile.h>

#include "term_colours.cpp"

namespace po = boost::program_options;

int main(int argc, char *argv[])
{
  bool epsOutput;
  bool logY;
  bool normalize;

  po::options_description desc("Options for my program");
    desc.add_options()
        ("epsOutput,e", po::value<bool>(& epsOutput)->default_value(false),
            "If true save canvas directly to an eps file.")

        ("logY,l", po::value<bool>(& logY)->default_value(false),
            "If true make y-axis logarithmic scale.")

        ("normalize,n", po::value<bool>(& normalize)->default_value(false),
            "If true normalize histos")
        ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

  // copy arguments (RootApp modifies them)
  const int args = argc;
  char *argsv[args];
  for (int i = 0; i < args; i++) argsv[i] = argv[i];

  std::string histName1;
  std::string fileName1; 
  std::string histName2;   
  std::string fileName2;
  std::string histName3;   
  std::string fileName3;

  TFile *f1, *f2, *f3;
  TH2D *h1, *h2, *h3;
  std::string name1, name2, name3;
  std::string epsFileName;

  int nFiles = (args - 1)/2;   
  printf("Superimposing %i 2D histograms...\n", nFiles);
  
  TApplication* rootApp = new TApplication("rootApp", &argc, argv);

  gStyle->SetOptStat(0);
  gStyle->SetLegendFont(132);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetOptTitle(0);
  
  TCanvas *canvas = new TCanvas("Superposition", "Superposition", 1280, 751);
  canvas->cd();

  histName1 = argsv[1];
  fileName1 = argsv[2];
  f1 = new TFile(fileName1.c_str(), "READ");
  if (!f1->IsOpen()) printf("Failed to open %s\n", fileName1.c_str());
  h1 = (TH2D*) f1->Get(histName1.c_str());
  name1 = histName1 + '@' + fileName1;
  h1->SetTitle(name1.c_str());
  h1->Draw("LEGO");
  h1->SetLineColor(kBlack);
  h1->GetYaxis()->SetTitleOffset(2.0);
  h1->GetXaxis()->SetTitleOffset(2.0);
  h1->GetYaxis()->SetRangeUser(0,4000);


  if (nFiles > 1) {
    histName2 = argsv[3];
    fileName2 = argsv[4];
    f2 = new TFile(fileName2.c_str(), "READ");
    if (!f2->IsOpen()) printf("Failed to open %s\n", fileName2.c_str());
    name2 = histName2 + '@' + fileName2;
    h2 = (TH2D*) f2->Get(histName2.c_str());
    h2->SetTitle(name2.c_str());
    h2->Draw("LEGOSAME");
    h2->SetLineColor(kRed);
    h2->GetYaxis()->SetRangeUser(0,4000);
  }

  if (nFiles > 2) {
    histName3 = argsv[5];
    fileName3 = argsv[6];
    f3 = new TFile(fileName3.c_str(), "READ");
    if (!f3->IsOpen()) printf("Failed to open %s\n", fileName3.c_str());
    h3 = (TH2D*) f3->Get(histName3.c_str());
    name3 = histName3 + '@' + fileName3;
    h3->SetTitle(name3.c_str());
    h3->Draw("LEGOSAME");
    h3->SetLineColor(kBlue);
    h3->GetYaxis()->SetRangeUser(0,4000);
  }

  // normalize histograms
  if (normalize == true) { 
    std::string yTitle;
    yTitle = h1->GetYaxis()->GetTitle();
    yTitle = "1/#sigma #times " + yTitle;
    h1->GetYaxis()->SetTitle(yTitle.c_str());
    if (nFiles > 0) h1->Scale(1.0/h1->Integral());
    if (nFiles > 1) h2->Scale(1.0/h2->Integral());
    if (nFiles > 2) h3->Scale(1.0/h3->Integral());
  }

  // find overlapping area (histograms must be have the same user ranges and same number of bins)
  double overlap = 0;
  for (int i = 0; i < h1->GetNbinsX(); i++) {
    for (int j = 0; j < h1->GetNbinsY(); j++) {
      int bin1 = h1->GetBin(i, j);
      int bin2 = h2->GetBin(i, j);
      double one = h1->GetBinContent(bin1);
      double two = h2->GetBinContent(bin2);
      if (one <= two) overlap += one;
      else if (one > two) overlap += two;
    }
  }

  printf("Overlapping volume: %f\n", overlap);

  if (logY == true) canvas->SetLogy();

  canvas->BuildLegend(0.7, 0.9, 1.0, 1.0, "");
  
  if (epsOutput == false) {
    printf("Running ROOT app...\n");
    canvas->SetEditable(true);
    canvas->SetBit(2048, 0);
    canvas->Show();
    rootApp->Run();
  }
  else if (epsOutput == true) {
    epsFileName = histName1 + "_2to6_SMMvsSM_EW.eps";
    printf("Saving to %s\n",epsFileName.c_str());
    canvas->SaveAs(epsFileName.c_str());
  }

  printf("Superposition complete. Have a nice day.\n");
  return 0; 
}