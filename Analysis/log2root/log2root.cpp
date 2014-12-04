// ================================== header ===================================
// program log2root
// -------------------
// Author: Declan Millar
// ================================== start ====================================
// includes
#include <iostream>
#include <fstream>
#include <string>

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
#include "AtlasROOTStyle.cpp"
// ================================== main =====================================
int main(int argc, char *argv[])
{
  // Atlas style
  // AtlasROOTStyle ars;
  // ars.SetStyle();
    
  // check number of arguments is correct.
  if ( argc != 2 )
  {
    // argc should be 2 for correct execution.
    // print argv[0] assuming it is the program name.
    std::cout << "usage: " << argv[0] << " <filename>\n";
    return 0;
  }
  else
  {    
    // Declarations
    std::string fileName; // File name without extension
    std::string logFileName; // adds extension to file name  
    std::string outputFileName; // output .root file name
    std::string start; // starts reading a histogram
    std::string end; // stops reading the text file
    std::string mainString; // Information strings preceding distributions
    int nHistos; // number of histograms read from text file
    int lastIndex; // location of . in string

    // assume argv[1] is a filename to open.
    logFileName = argv[1];
    lastIndex = logFileName.find_last_of("."); 
    fileName = logFileName.substr(0, lastIndex); 
    outputFileName = fileName + ".root";

    // Triggers
    start = "HISTOGRAM";
    end = "CLOSE";

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
      if (mainString==start)
      {
        nHistos+=1;
        printf("Histogram %i\n",nHistos);
        TH1D* Distribution = plotDistribution(&logStream);
        Distribution -> Write();
        delete Distribution;
      }
    }

    // Wrapping up
    logStream.close();
    rootFile.Close();
    printf("Program complete.\n");
    return 0;
  }
}
// ================================ end main ===================================