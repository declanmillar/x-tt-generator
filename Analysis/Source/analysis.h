#ifndef _ANALYSIS_ZPRIME_H_
#define _ANALYSIS_ZPRIME_H_

#include "RootTuple.h"
#include "atlas_style.h"
#include "TCanvas.h"
#include "TApplication.h"
#include <cmath> 
#include <TString.h>
#include <TH2.h>
#include "TLorentzVector.h"
#include "TVector2.h"
#include "TBrowser.h"
#include <fstream>
#include "TLegend.h"
#include <complex>

class AnalysisZprime{
public:
  AnalysisZprime(const TString channel, const TString model, const TString& inputFileName, const TString& weightsFileName, const TString& outputFileName);
  virtual ~AnalysisZprime();
  
protected:
  // Internal for Event Looping
  Long64_t TotalEvents();
  Long64_t IncrementEvent(Long64_t i);
  void SetupTreesForNewFile(const TString& s);
  void CleanUp();
  
  void SetupInputFiles();
  void SetupOutputFiles();
  
  void PreLoop();
  void Loop();
  void PostLoop();
  void EachEvent();
  void CreateHistograms();
  void MakeGraphs();
  void WriteHistograms();
  double TotalAsymmetry(TH1D* h_A, TH1D* h_B);
  void TotalSpinAsymmetries();
  TH1D* Asymmetry(TString name, TString title, TH1D* h_A, TH1D* h_B);
  void AsymmetryUncertainty(TH1D* h_Asymmetry, TH1D* h_A, TH1D* h_B);
  void ALL2to6();
  TH1D* MttALL();
  TH1D* MttAL();
  double deltaPhi(const double& phi1,const double& phi2) const;
  vector<std::complex<double> > SolveQuadratic(double a, double b, double c);
  double ResolveNeutrinoPz(TLorentzVector p_l, TVector2 pT_nu);
  
  bool PassCuts();
  bool PassCuts_MET();
  bool PassCuts_Mtt();
  bool PassCutsFiducial();
  bool PassCutsYtt();

  void InitialiseCutflow();
  void PrintCutflow();
  const void UpdateCutflow(int cut, bool passed);
     
private:
  AnalysisZprime();
  AnalysisZprime(const AnalysisZprime& rhs);  
  void operator=(const AnalysisZprime& rhs);  

  m_bAssignAttempts;
  m_bAssignSuccesses;
  
  float m_pi;
  float m_GeV;
  double m_intLumi;
  double m_Wmass;
  TString m_channel;  
  TString m_model;
  TString m_inputFileName;
  TString m_weightsFileName;
  TString m_outputFileName;
  double m_sigma;
  vector<double> m_weights;
  enum m_cutlist{
    c_Event,
    c_Mtt,
    c_MET,
    c_Ytt,
    c_Fiducial,
    m_cuts // keep as last entry
  };


  // Input Data
  vector<TString>* m_inputFiles;
  RootTuple* m_ntup;
  TChain* m_chainNtup; 
  
  // OutputFile
  TFile* m_outputFile;

  // cutflow
  TH1D* h_cutflow;
  std::vector<int> m_cutflow;
  std::vector<TString> m_cutNames;

  // regualar
  TH1D* h_Mff;
  TH1D* h_ytt;
  TH1D* h_Pz_nu;
  TH1D* h_CosTheta;
  TH1D* h_CosThetaStar;
  
  // reconstruction
  TH1D* h_Pz_nu_r;
  TH1D* h_Mtt_r;
  TH1D* h_CosTheta_r;
  TH1D* h_CosThetaStar_r;
  TH1D* h_ytt_r;

  // polarisation weighted
  TH1D* h_Mff_LL;
  TH1D* h_Mff_LR;
  TH1D* h_Mff_RL;
  TH1D* h_Mff_RR;

  // Spin Asymmetries
  TH1D* h_ALL;
  TH1D* h_AL;

  // FB-type histograms
  TH1D* h_AFBstar;
  TH1D* h_AFBstarF;
  TH1D* h_AFBstarB;

  TH1D* h_AFBstar_r;
  TH1D* h_AFBstar_rF;
  TH1D* h_AFBstar_rB;

  TH1D* h_AttC;
  TH1D* h_AttCF;
  TH1D* h_AttCB;

  TH1D* h_AllC;  
  TH1D* h_AllCF;
  TH1D* h_AllCB;

  TH1D* h_AlL;
  TH1D* h_AlLF;
  TH1D* h_AlLB;

  // kinematics
  vector<TLorentzVector> p;
  vector<TLorentzVector> pcm;
  vector<TLorentzVector> p_r1;
  vector<TLorentzVector> p_r2;
  vector<TLorentzVector> pcm_r1;
  vector<TLorentzVector> pcm_r2;

  vector<TLorentzVector> P;
  vector<TLorentzVector> Pcm;
  vector<TLorentzVector> P_r1;
  vector<TLorentzVector> P_r2;
  vector<TLorentzVector> Pcm_r1;
  vector<TLorentzVector> Pcm_r2;

  vector<TLorentzVector> p_t;
  vector<TLorentzVector> p_tb;
  vector<TLorentzVector> p_t_r1;
  vector<TLorentzVector> p_tb_r1;
  vector<TLorentzVector> p_t_r2;
  vector<TLorentzVector> p_tb_r2;

  typedef vector<TString>::const_iterator Itr_s;
};
#endif
