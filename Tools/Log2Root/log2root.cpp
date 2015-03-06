// ================================== header ===================================
// program log2root
// -------------------
// Author: Declan Millar
// ================================== start ====================================
// includes
#include <iostream>
#include <fstream>
#include <string>
#include <boost/program_options.hpp>

#include <stdlib.h>
#include <TVector.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TApplication.h>
#include <TAxis.h>
#include <TH1.h>
#include <TFile.h>

#include "plotDistribution.cpp"
#include "plotAsymmetry.cpp"
#include "AtlasROOTStyle.cpp"

namespace po = boost::program_options;
// ================================== main =====================================
int main(int argc, char *argv[])
{
  // Atlas style
  // AtlasROOTStyle ars;
  // ars.SetStyle();

  // check number of arguments is correct.
  if ( argc < 2 )
  {
    // argc should be 2 for correct execution.
    // print argv[0] assuming it is the program name.
    std::cout << "usage: " << argv[0] << " <filename>\n";
    return 0;
  }

  // Declarations
  std::string fileName; // File name without extension
  std::string logFileName; // adds extension to file name  
  std::string outputFileName; // output .root file name
  std::string dist,asym; // starts reading a histogram
  std::string end; // stops reading the text file
  std::string mainString; // Information strings preceding distributions
  int nHistos; // number of histograms read from text file
  int lastIndex; // location of . in string
  double luminosity,efficiency;

  // assume argv[1] is a filename to open.
  logFileName = argv[1];
  lastIndex = logFileName.find_last_of("."); 
  fileName = logFileName.substr(0, lastIndex); 
  outputFileName = fileName + ".root";

  po::options_description desc("Options for my program");
    desc.add_options()
      ("help,h", "Display this help message and quit.")
      ("luminosity,L", po::value<double>(& luminosity)->default_value(0),
      "Multiply the cross section by the luminosity in pb-1 for the number of events.")
      ("efficiency,e", po::value<double>(& efficiency)->default_value(0),
      "Multiply the Number of events by the reconstruction efficiency.")
        // ("logY,l", po::value<bool>(& logY)->default_value(false),
        //     "If true make y-axis logarithmic scale.")

        // ("normalize,n", po::value<bool>(& normalize)->default_value(false),
        //     "If true normalize histos")

        // both rangeMin and rangeMax must be set simultaneously for any effect.
        ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

  // Help message
  if (vm.count("help")) {
    std::cout << desc << "\n";
    return 1;
  }

  // Triggers
  dist = "DISTRIBUTION";
  asym = "ASYMMETRY";
  end  = "CLOSE";

  // convert luminosity in fb to pb
  luminosity = luminosity*1000;

  // input stream.
  ifstream logStream(logFileName);  
  if (!logStream.is_open()){
   printf("Failed to open: %s\n",logFileName.c_str());
   return 0;
  }    

  // Output ROOT file
  TFile rootFile(outputFileName.c_str(), "RECREATE");
  printf("Output file name: %s\n",outputFileName.c_str());
  
  // Loop over all Histograms
  nHistos=0;
  while (1)
  {
    logStream >> mainString;
    // printf("Main loop string: %s\n",mainString.c_str());
    if (mainString==end || !logStream.good())
    {
      // printf("Main loop string: %s\n",mainString.c_str());
      break;
    }
    if (mainString==dist)
    {
      nHistos+=1;
      printf("Histogram %i\n",nHistos);
      TH1D* Distribution = plotDistribution(efficiency,luminosity,&logStream);
      Distribution -> Write();
      delete Distribution;
    }
    else if (mainString==asym)
    {
      nHistos+=1;
      printf("Histogram %i\n",nHistos);
      TH1D* Asymmetry = plotAsymmetry(efficiency,luminosity,&logStream);
      Asymmetry -> Write();
      delete Asymmetry;
    }
  }

  // Wrapping up
  logStream.close();
  rootFile.Close();
  printf("Program complete.\n");
  return 0;
}
// ================================ end main ===================================