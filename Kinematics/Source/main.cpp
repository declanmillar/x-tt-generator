#include "kinematics.h"

int main(int argc, char* argv[])
{
  AtlasROOTStyle atlasStyle;
  atlasStyle.SetStyle();

  TString inDir("/afs/cern.ch/work/d/demillar/Ntuples_Zprime/");
  TString outDir("/afs/cern.ch/work/d/demillar/ProcessedNtuples_Zprime/");
  // TString inDir("/Users/declan/Data/Ntuples_Zprime/");
  // TString outDir("/Users/declan/Data/Histos_Zprime/");
  if (argc != 2) printf("Usage: ./go <inputFileName>.root\n");

  TString inputFileName(inDir + argv[1]);
  TString outputFileName(outDir + argv[1] + ".p");  

  KinematicsZprime("2to6", inputFileName, outputFileName);

  return 0;  
}
