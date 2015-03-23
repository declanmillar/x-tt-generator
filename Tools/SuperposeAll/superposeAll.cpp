// superpose
// author: Declan Millar

#include <iostream>
#include <fstream>
#include <string>
#include <boost/program_options.hpp>

#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TApplication.h>
#include <TAxis.h>
#include <TH1.h>
#include <TFile.h>
#include <TObject.h>

namespace po = boost::program_options;

int main(int argc, char *argv[])
{
  bool epsOutput;
  bool logY;
  bool normalize;

  po::options_description desc("Options for my program");
    desc.add_options()
        ("epsOutput, e", po::value<bool>(& epsOutput)->default_value(false),
            "If true save canvas directly to an eps file.")

        ("logY, l", po::value<bool>(& logY)->default_value(false),
            "If true make y-axis logarithmic scale.")

        ("normalize, n", po::value<bool>(& normalize)->default_value(false),
            "If true normalize histos")
        ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

  std::vector<TCanvas*> vcanvas;

  std::string histName;
  std::string fileName1;    
  std::string fileName2; 
  std::string fileName3;

  TFile *f1, *f2, *f3;
  TH1D *h1, *h2, *h3;
  std::string outName;

  // copy arguments (RootApp modifies them)
  const int args = argc;
  char *argsv[args];
  for (int i = 0; i < args; i++) argsv[i] = argv[i];

  int nFiles = (args - 2);  

  outName = argv[1];
  outName += ".root";

  fileName1 = argsv[2];
  if (nFiles > 1) fileName2 = argsv[3];
  else fileName2 = "";
  if (nFiles > 2) fileName3 = argsv[4];
  else fileName3 = "";

  printf("Superimposing all histograms in %i files...\n", nFiles);

  f1 = new TFile(fileName1.c_str(), "READ");
  TList *list = f1->GetListOfKeys();
  int size = list->GetSize();
  printf("here74 %i\n",size);
  for (int i = 1; i < size; i++) {

    printf("here77\n");

    const char* histName = list->At(i)->GetName();
    TCanvas *canvas = new TCanvas(histName, "", 1280, 751);
    canvas->cd(); 
    h1 = (TH1D*) f1->Get(list->At(i)->GetName());
    h1->SetTitle(fileName1.c_str());
    h1->Draw();
    h1->SetStats(false);
    h1->SetLineColor(kBlack);
    h1->GetYaxis()->SetTitleOffset(1.3);
    h1->GetXaxis()->SetTitleOffset(1.2);

    printf("here92\n");

    if (nFiles > 1) { 
      f2 = new TFile(fileName2.c_str(), "READ");    
      h2 = (TH1D*) f2->Get(list->At(i)->GetName());
      f2->Close();
      h2->SetTitle(fileName2.c_str());
      h2->Draw("SAME");
      h2->SetLineColor(kRed);
    }

    if (nFiles > 2) {   
      f3 = new TFile(fileName3.c_str(), "READ");
      h3 = (TH1D*) f3->Get(list->At(i)->GetName());
      f3->Close();
      h3->SetTitle(fileName3.c_str());
      h3->Draw("SAME");
      h3->SetLineColor(kBlue);
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

    if (logY == true) canvas->SetLogy();

    canvas->BuildLegend(0.50, 0.60, 0.88, 0.88, "");

    printf("here125\n");

    vcanvas.push_back(canvas);

    printf("here129\n");    

    // h1->Delete();
    // h2->Delete();
    // h3->Delete();
    // canvas->Delete();
    printf("here135\n");   
  }

  f1->Close();

  printf("here140%i\n",size); 

  TFile *outFile = new TFile(outName.c_str(), "RECREATE");
  for (int i = 1; i < vcanvas.size(); i++) {
    printf("test%i\n", i);
    vcanvas.at(i)->Write();
  }
  outFile->Close();

  printf("Superposition complete. Have a nice day.\n");
  return 0; 
}