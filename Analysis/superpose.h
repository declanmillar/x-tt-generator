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

TCanvas* superpose(std::string histName1, std::string fileName1, 
                   std::string histName2 = "NULL", std::string fileName2 = "NULL", 
                   std::string histName3 = "NULL", std::string fileName3 = "NULL");

#endif