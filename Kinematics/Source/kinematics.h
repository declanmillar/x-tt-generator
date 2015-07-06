#ifndef _KINEMATICS_ZPRIME_H_
#define _KINEMATICS_ZPRIME_H_

#include "RootTuple.h"
#include "atlas_style.h"
#include <cmath> 
#include <TString.h>
#include "TLorentzVector.h"
#include "TVector2.h"
#include <fstream>
#include <complex>

class KinematicsZprime{
public:
  KinematicsZprime(const TString& channel, const TString& inputFileName, const TString& outputFileName);
  virtual ~KinematicsZprime();
  
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

  double deltaPhi(const double& phi1,const double& phi2) const;
  vector<std::complex<double> > solveQuadratic(double a, double b, double c);
  double resolveNeutrinoPz(TLorentzVector p_l, TVector2 pT_nu);
     
private:
  KinematicsZprime();
  KinematicsZprime(const KinematicsZprime& rhs);  
  void operator=(const KinematicsZprime& rhs);  
  
  float m_pi;
  float m_GeV;
  double m_intLumi;
  double m_Wmass;
  TString m_channel;  
  TString m_model;
  TString m_inputFileName;
  TString m_outputFileName;
  double m_sigma;
  vector<double> m_cnorm;

  // Input Data
  vector<TString>* m_inputFiles;
  RootTuple* m_ntup;
  TChain* m_chainNtup; 
  
  // OutputFile
  TFile* m_outputFile;

  // Branches
  tree->Branch("weight")

  tree->Branch("e_chargeQ",m_e_Q);
  tree->Branch("e_E",m_e_E);
  tree->Branch("e_Px",m_e_Px);
  tree->Branch("e_Px",m_e_Py);
  tree->Branch("e_Px",m_e_Pz);

  tree->Branch("nu_E",m_nu_E);
  tree->Branch("nu_Px",m_nu_Px);
  tree->Branch("nu_Px",m_nu_Py);
  tree->Branch("nu_Px",m_nu_Pz);

  tree->Branch("b_Q",m_b_Q);
  tree->Branch("b_E",m_b_E);
  tree->Branch("b_Pz",m_b_Py);
  tree->Branch("b_Py",m_b_Pz);
  tree->Branch("b_Px",m_b_Px);

  tree->Branch("q_Q",m_q_Q);
  tree->Branch("q_E",m_q_E);
  tree->Branch("q_Px",m_q_Px);
  tree->Branch("q_Pz",m_q_Py);
  tree->Branch("q_Py",m_q_Pz);

  tree->Branch("MET_E",m_MET_E);
  tree->Branch("MET_Px",m_MET_Px);
  tree->Branch("MET_Py",m_MET_Py);

  typedef vector<TString>::const_iterator Itr_s;
};
#endif
