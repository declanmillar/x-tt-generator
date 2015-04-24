#include "superpose.h"

TCanvas* superpose(int ndimensions,
                   std::string histName1, std::string fileName1, 
                   std::string histName2, std::string fileName2, 
                   std::string histName3, std::string fileName3) {

  TFile *f1, *f2, *f3;
  TH1D *h1, *h2, *h3;
  TH2D *h1, *h2, *h3;
  std::string name1, name2, name3;
  std::string epsFileName;



  // gStyle->SetOptStat(0);
  // gStyle->SetLegendFont(132);
  // gStyle->SetLegendBorderSize(0);
  // gStyle->SetOptTitle(0);
  
  TCanvas *canvas = new TCanvas("test", "test", 1280, 751);
  canvas->cd();

  f1 = new TFile(fileName1.c_str(), "READ");
  if (!f1->IsOpen()) printf("Failed to open %s\n", fileName1.c_str());
  h1 = (TH1D*) f1->Get(histName1.c_str());
  name1 = histName1 + '@' + fileName1;
  h1->SetTitle(name1.c_str());
  h1->SetMarkerStyle(20);
  h1->Draw();//("E1x0P");
  h1->SetMarkerColor(kBlack);
  h1->SetLineColor(kBlack);
  h1->GetYaxis()->SetTitleOffset(1.3);
  h1->GetXaxis()->SetTitleOffset(1.2);

  if (histName2 != "NULL") {
    f2 = new TFile(fileName2.c_str(), "READ");
    if (!f2->IsOpen()) printf("Failed to open %s\n", fileName2.c_str());
    name2 = histName2 + '@' + fileName2;
    h2 = (TH1D*) f2->Get(histName2.c_str());
    h2->SetTitle(name2.c_str());
    h2->SetMarkerStyle(20);
    h2->Draw("SAME");//("E1x0PSAME");
    h2->SetLineColor(kRed);
    h2->SetMarkerColor(kRed);
  }

  if (histName3 != "NULL") {
    f3 = new TFile(fileName3.c_str(), "READ");
    if (!f3->IsOpen()) printf("Failed to open %s\n", fileName3.c_str());
    h3 = (TH1D*) f3->Get(histName3.c_str());
    name3 = histName3 + '@' + fileName3;
    h3->SetTitle(name3.c_str());
    h3->SetMarkerStyle(20);
    h3->Draw("SAME");//("E1x0PSAME");
    h3->SetLineColor(kBlue);
    h3->SetMarkerColor(kBlue);
  }

  // normalize histograms
  bool normalize = true;
  if (normalize == true) { 
    std::string yTitle;
    yTitle = h1->GetYaxis()->GetTitle();
    yTitle = "1/#sigma #times " + yTitle;
    h1->GetYaxis()->SetTitle(yTitle.c_str());
    h1->Scale(1.0/h1->Integral("width"));
    if (histName2 != "NULL") h2->Scale(1.0/h2->Integral("width"));
    if (histName3 != "NULL") h3->Scale(1.0/h3->Integral("width"));
  }

  // set range user
  // if(rangeMax == rangeMax && rangeMax == rangeMax) {
  //   h1->GetXaxis()->SetRangeUser(rangeMin, rangeMax); 
  //   if (nFiles > 1) h2->GetXaxis()->SetRangeUser(rangeMin, rangeMax);
  //   if (nFiles > 2) h3->GetXaxis()->SetRangeUser(rangeMin, rangeMax);
  // }

  // if (logY == true) canvas->SetLogy();

  canvas->BuildLegend(0.50, 0.60, 0.88, 0.88, "");

  return canvas;
}