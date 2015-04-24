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
  AnalysisTransverse(const TString& inputFileName, const TString& outputFileName);
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
  TString m_inputFileName;
  TString m_outputFileName;

  // Input Data
  vector<TString>* m_inputFiles;
  RootTuple* m_ntup;
  TChain* m_chainNtup; 
  
  // OutputFile
  TFile* m_outputFile;
  
  // Define Histograms
  TH1D* h_mtt;
  TH1D* h_ht;
  TH1D* h_mttvis;
  TH1D* h_mt1;
  TH1D* h_mt2;
  TH1D* h_mt3;
  TH1D* h_mct1;
  TH1D* h_mct2;
  TH1D* h_mct3;
  TH1D* h_mlt;
  TH1D* h_mlct;
  TH2D* h_dphi_mtt;
  TH2D* h_dphi_ht;
  TH2D* h_dphi_mttvis;
  TH2D* h_dphi_mt1;
  TH2D* h_dphi_mt2;
  TH2D* h_dphi_mt3;
  TH2D* h_dphi_mct1;
  TH2D* h_dphi_mct2;
  TH2D* h_dphi_mct3;
  TH2D* h_dphi_mlt;
  TH2D* h_dphi_mlct;
  typedef vector<TString>::const_iterator Itr_s;
};
#endif
