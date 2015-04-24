#ifndef _ANALYSIS_ZPRIME_H_
#define _ANALYSIS_ZPRIME_H_

#include "RootTuple.h"
#include "atlas_style.h"
#include "TCanvas.h"
#include "TApplication.h"
#include <cmath> 
#include <TString.h>
#include <TH2.h>

class AnalysisTransverse{
public:
  AnalysisTransverse(const TString& outputFileName);
  virtual ~AnalysisTransverse();
  
protected:
  // Internal for Event Looping
  Long64_t TotalEvents();
  Long64_t IncrementEvent(Long64_t i);
  void SetupTreesForNewFile(const TString& s);
  void CleanUp();
  
  void SetupInputFiles();
  
  void PreLoop();
  void Loop();
  void EachEvent();
  void PostLoop();
  void MakeGraphs();
  
  bool PassCuts() const;
  // bool PassCuts_NJets() const;
  // bool PassCuts_MET() const;
  // bool PassCuts_MWT() const;
    
  // Delta R
  // float deltaR(const float& eta1,const float& eta2,const float& phi1,const float& phi2) const;
  // float deltaPhi(const float& phi1,const float& phi2) const;    
  
  
private:
  AnalysisTransverse();
  AnalysisTransverse(const AnalysisTransverse& rhs);  
  void operator=(const AnalysisTransverse& rhs);  
  
  float m_pi;
  float m_GeV;
  // TString m_channel;  
  TString m_outputFileName;

  // Input Data
  vector<TString>* m_inputFiles;
  RootTuple* m_ntup;
  TChain* m_chainNtup; 
  
  // OutputFile
  TFile* m_outputFile;
  
  // Define Histograms
  TH1D* h_positron_pz;
  TH1D* h_mtt;
  TH1D* h_mttLL;
  TH1D* h_mttLR;
  TH1D* h_mttRL;
  TH1D* h_mttRR;
  TH1D* h_mttALLnum;
  TH1D* h_mttALLden;
  TH1D* h_mttALL;
  typedef vector<TString>::const_iterator Itr_s;
};
#endif
