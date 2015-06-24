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

  typedef vector<TString>::const_iterator Itr_s;
};
#endif
