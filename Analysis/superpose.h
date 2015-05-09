#ifndef _SUPERPOSE_
#define _SUPERPOSE_

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

TCanvas* superpose(TString histName1, TString fileName1, 
                   TString histName2 = "NULL", TString fileName2 = "NULL", 
                   TString histName3 = "NULL", TString fileName3 = "NULL");

#endif