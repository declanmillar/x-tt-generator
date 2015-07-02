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
  vector<double> m_cnorm;
  // int m_cuts = 2;
  enum m_cuts{
    cutEvent,
    cutMtt,
    cutMET,
    cutYtt,
    cutFiducial,
    nCuts // keep as last entry
  };


  // Input Data
  vector<TString>* m_inputFiles;
  RootTuple* m_ntup;
  TChain* m_chainNtup; 
  
  // OutputFile
  TFile* m_outputFile;

  // cutflow
  TH1D* h_cutflow;
  std::vector<int> cutflow;
  std::vector<TString> cutNames;
  
  // Angular histograms
  TH1D* h_dphi;
  TH1D* h_costheta5;
  TH1D* h_ct7ct5; 
  TH1D* h_CosTheta;
  TH1D* h_CosThetaReco;
  TH1D* h_CosThetaStar;

  // Rapidity
  TH1D* h_ytt;
  TH1D* h_yttReco;

  // asymmetry histograms
  TH1D* h_AFBstar;
  TH1D* h_AFBstar1;
  TH1D* h_AFBstar2;

  TH1D* h_AFBstarReco;
  TH1D* h_AFBstarReco1;
  TH1D* h_AFBstarReco2;

  TH1D* h_AttC;
  TH1D* h_AttC1;
  TH1D* h_AttC2;

  TH1D* h_AllC;  
  TH1D* h_AllCF;
  TH1D* h_AllCB;

  TH1D* h_AL;
  TH1D* h_AlLF;
  TH1D* h_AlLB;

  TH1D* h_ALL;
  TH1D* h_ALLF;
  TH1D* h_ALLB;

  // neutrino resolution
  TH1D* h_PzNu;
  TH1D* h_PzNuReco;
  TH1D* h_Mtt;
  TH1D* h_MttReco;
  TH1D* h_CosThetaStarReco;

  // polarisation weighted
  TH1D* h_MttLL;
  TH1D* h_MttLR;
  TH1D* h_MttRL;
  TH1D* h_MttRR;

  // transverse variables 
  TH1D* h_MET;
  TH1D* h_HT;
  TH1D* h_Mbbll;
  TH1D* h_mll;
  TH1D* h_ETbbll;
  TH1D* h_KTbbll;
  TH1D* h_MTll;
  TH1D* h_MCTll;
  TH1D* h_MTblbl;
  TH1D* h_MCTblbl;

  // 2d Histograms
  TH2D* h_dphi_HT;
  TH2D* h_dphi_Mbbll;
  TH2D* h_dphi_mll;
  TH2D* h_dphi_ETbbll;
  TH2D* h_dphi_KTbbll;
  TH2D* h_dphi_MTll;
  TH2D* h_dphi_MCTll;
  TH2D* h_dphi_MTblbl;
  TH2D* h_dphi_MCTblbl;

  // kinematics
  vector<TLorentzVector> pcol;
  vector<TLorentzVector> p;
  vector<TLorentzVector> pcolRecot;
  vector<TLorentzVector> pRecot;
  vector<TLorentzVector> pcolRecotb;
  vector<TLorentzVector> pRecotb;
  double Mff;
  vector<double> ycol;
  vector<double> etacol;
  vector<TVector2> pTcol;

  typedef vector<TString>::const_iterator Itr_s;
};
#endif
