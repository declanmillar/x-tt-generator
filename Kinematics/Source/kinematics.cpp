#include "kinematics.h"

KinematicsZprime::KinematicsZprime(const TString& channel, const TString& inputFileName, const TString& outputFileName) :
  m_pi(3.14159265),
  m_GeV(1000.0),
  m_intLumi(0),//(300000.0),
  m_Wmass(80.23),
  m_channel(channel),
  m_inputFileName(inputFileName),
  m_outputFileName(outputFileName),
  m_inputFiles(NULL),
  m_ntup(NULL),
  m_chainNtup(NULL),
  m_outputFile(NULL)
{  
  this->PreLoop();
  this->Loop();
  this->PostLoop();
}

void KinematicsZprime::EachEvent()
{
  // 0 = b, 1 = bbar, 2 = l+, 3 = nu, 4 = l-, 5 = nubar

  vector<TLorentzVector> pcol(6);
  vector<TLorentzVector> p(6);
  vector<TLorentzVector> pcolReco(6);
  vector<TLorentzVector> pReco(6);
  TLorentzVector pcoltot, ptot;
  vector<double> mass(6);
  vector<double> ETcol(6);
  vector<TVector2> pTcol(6);
  vector<double> ycol(6);
  vector<double> etacol(6);
  vector<double> phicol(6);
  vector<double> ET(6);
  vector<TVector2> pT(6);
  vector<double> y(6);
  vector<double> eta(6);
  vector<double> phi(6);
  vector<TVector2> pTReco(6);
  vector<double> yNuReco(6);
  vector<double> etaNuReco(6);
  vector<double> phiNuReco(6);
  vector<double> ETNuReco(6);
  TVector2 pTtotNuReco;
  TVector2 pTtot;
  TLorentzVector ptcol;
  TLorentzVector ptbcol;
  TLorentzVector pt;
  TLorentzVector ptb;
  TLorentzVector ptcolReco;
  TLorentzVector ptReco;

  // final particle collider variables
  TVector2 pTcoltot;
  for (int i = 0; i < (int) m_ntup->E()->size(); i++) {
    pcol[i].SetPxPyPzE(m_ntup->Px()->at(i), m_ntup->Py()->at(i), m_ntup->Pz()->at(i), m_ntup->E()->at(i));
    p[i] = pcol[i];
    pReco[i] = pcol[i];
    pcolReco[i] = pcol[i];
    pTcol[i].Set(pcol[i].Px(), pcol[i].Py());
    ycol[i] = pcol[i].Rapidity();
    etacol[i] = pcol[i].PseudoRapidity();
    phicol[i] = pcol[i].Phi();   
    mass[i] = pcol[i].M();
    ETcol[i] = std::sqrt(mass[i]*mass[i] + pTcol[i].Mod2());
    pcoltot += pcol[i];
    pTcoltot += pTcol[i];
  }
  double ytt = pcoltot.Rapidity();

  // negative velocity of full system in collider frame
  TVector3 vcoltot = -1*pcoltot.BoostVector();

  // true final particle parton CoM variables
  for (int i = 0; i < (int) m_ntup->E()->size(); i++) {
    p[i].Boost(vcoltot);
    pT[i].Set(p[i].Px(), p[i].Py());
    y[i] = p[i].Rapidity();
    eta[i] = p[i].PseudoRapidity();
    phi[i] = p[i].Phi();   
    ET[i] = std::sqrt(mass[i]*mass[i] + pT[i].Mod2());
    ptot += p[i];
    pTtot += pT[i];
  }  

  double PzNuRecoCol = -999;
  double MttReco =-999;
  double yttReco =-999;
  if (m_channel =="2to6") {
    // resolve longitudinal neutrino momentum in the semihadronic case
    PzNuRecoCol = resolveNeutrinoPz(pcol[2], pTcol[3]);
    pcolReco[3].SetPxPyPzE(pcol[3].Px(), pcol[3].Py(), PzNuRecoCol, std::sqrt(pcol[3].Px()*pcol[3].Px()+pcol[3].Py()*pcol[3].Py()+PzNuRecoCol*PzNuRecoCol));
    pReco[3].SetPxPyPzE(pcol[3].Px(), pcol[3].Py(), PzNuRecoCol, std::sqrt(pcol[3].Px()*pcol[3].Px()+pcol[3].Py()*pcol[3].Py()+PzNuRecoCol*PzNuRecoCol));
    TLorentzVector pcoltotReco = pcoltot - pcol[3] + pcolReco[3];
    MttReco = pcoltotReco.M();
    yttReco = pcoltotReco.Rapidity();
    // reconstructed final particle parton CoM variables
    TVector3 vcoltotReco = -1*pcoltotReco.BoostVector();
    for (int i = 0; i < (int) m_ntup->E()->size(); i++) {
      pReco[i].Boost(vcoltotReco);
      // pReco[i].Print();
      pTReco[i].Set(pReco[i].Px(), pReco[i].Py());
      yNuReco[i] = pReco[i].Rapidity();
      // etaNuReco[i] = pReco[i].PseudoRapidity();
      phiNuReco[i] = pReco[i].Phi();   
      ETNuReco[i] = std::sqrt(mass[i]*mass[i] + pTReco[i].Mod2());
      pTtotNuReco += pTReco[i];
    }  
  }

  // top and antitop
  if (m_channel == "2to2") {
    ptcol = pcol[0];
    ptbcol = pcol[1];
    pt = p[0];
    ptb = p[1];
  }
  else if (m_channel == "2to6") {
    ptcol = pcol[0] + pcol[2] + pcol[3];
    ptbcol = pcol[1] + pcol[4] + pcol[5];
    pt = p[0] + p[2] + p[3];
    ptb = p[1] + p[4] + p[5];
    ptcolReco = pcolReco[0] + pcolReco[2] + pcolReco[3];
    ptReco = pReco[0] + pReco[2] + pReco[3];
  }

  double ytcol = ptcol.Rapidity();
  double ytbcol = ptbcol.Rapidity();
  // double etatcol = ptcol.PseudoRapidity();
  // double etatbcol = ptbcol.PseudoRapidity();
  // double phitcol = ptcol.Phi();
  // double phitbcol = ptbcol.Phi();
  // double pTtcol = ptcol.Pt();
  // double PTtbcol = ptbcol.Pt();
  // double yt = pt.Rapidity();
  // double ytb = ptb.Rapidity();
  // double etat = pt.PseudoRapidity();
  // double etatb = ptb.PseudoRapidity();
  // double phit = pt.Phi();
  // double phitb = ptb.Phi();
  // double pTt = pt.Pt();
  // double PTtb = ptb.Pt();
  double deltaYcol = std::abs(ytcol) - std::abs(ytbcol);
  // double CosThetaCol = ptcol.CosTheta();
  double CosTheta = pt.CosTheta();
  double CosThetaReco = ptReco.CosTheta();
  double CosThetaStar = int(ytt/std::abs(ytt))*CosTheta;
  double CosThetaStarReco = int(yttReco/std::abs(yttReco))*CosThetaReco;

  // negative velocity of t/tb in collider frame
  TVector3 vtcol = -1*ptcol.BoostVector();
  TVector3 vtbcol = -1*ptbcol.BoostVector();

  TLorentzVector plepptop;
  TLorentzVector plepmtop;

  int it = m_ntup->iteration();
  double Mtt = std::abs(pcoltot.M());
  double costhetalpcol = -999;
  double costhetalmcol = -999;
  double costhetalptop = -999;
  double costhetalmtop = -999;
  double clpclmcol = -999;
  double clpclmtop = -999;
  double dphi = -999;
 //  double MET = -999;
 //  double mll = -999;
	// double Mbbll = -999;
	// double HT = -999;
	// double ETbbll = -999;
	// double KTbbll = -999;
	// double ET5 = -999;
	// double ET7 = -999;
	// double MTll = -999;
	// double MCTll = -999;
	// double m35 = -999;
	// double m47 = -999;
	// double ET35 = -999;
	// double ET47 = -999;
	// double MTblbl = -999;
	// double MCTblbl = -999;
  double deltaEtal = -999;

  if (m_channel == "2to6") {
    TLorentzVector plepptop = pcol[2];
    TLorentzVector plepmtop = pcol[4];
    plepptop.Boost(vtcol);
    plepmtop.Boost(vtbcol);

    costhetalpcol = pcol[2].CosTheta();
    costhetalmcol = pcol[4].CosTheta();
    costhetalptop = plepptop.CosTheta();
    costhetalmtop = plepmtop.CosTheta();
    clpclmcol = costhetalpcol*costhetalmcol;
    clpclmtop = costhetalptop*costhetalmtop;
    dphi = deltaPhi(phicol[2],phicol[4]);

    // Delta |eta| = |eta_l+| - |eta_l-|
    deltaEtal = std::abs(eta[2]) - std::abs(eta[4]);

    // TRANSVERSE VARIABLES

    // TVector2 ETmiss = -1*pTcol[0] - pTcol[1] - pTcol[2] - pTcol[4];
    // MET = ETmiss.Mod();

    // mll = (pcol[2] + pcol[4]).M();

    // // calculate invariant mass of visible decay products
    // Mbbll = (pcol[0] + pcol[1] + pcol[2] + pcol[4]).M();

    // // calculate total scalar sum of transverse energy
    // HT = ETcol[0] + ETcol[1] + ETcol[2] + ETcol[4] + MET;

    // // ET of visible decay products
    // TVector2 pTbbll = pTcol[0] + pTcol[1] + pTcol[2] + pTcol[4];
    // ETbbll = std::sqrt(Mbbll*Mbbll + pTbbll.Mod2());

    // // scalar sum of visible decay products and MET
    // KTbbll = ETbbll + MET;

    // ET5 = std::sqrt(mass[2]*mass[2] + pTcol[2].Mod2());
    // ET7 = std::sqrt(mass[4]*mass[4] + pTcol[4].Mod2());

    // MTll = std::sqrt((ET5 + ET7)*(ET5 + ET7) - (pTcol[2] + pTcol[4]).Mod2());
    // MCTll = std::sqrt((ET5 + ET7)*(ET5 + ET7) - (pTcol[2] - pTcol[4]).Mod2());

    // m35 = (pcol[0] + pcol[2]).M();
    // m47 = (pcol[1] + pcol[3]).M();
    // TVector2 pT35 = pTcol[0] + pTcol[2];
    // TVector2 pT47 = pTcol[1] + pTcol[4];
    // ET35 = std::sqrt(m35*m35 - (pTcol[0] + pTcol[2]).Mod2());
    // ET47 = std::sqrt(m47*m47 - (pTcol[1] + pTcol[4]).Mod2());

    // MTblbl = std::sqrt((ET35 + ET47)*(ET35 + ET47) - (pTcol[0] + pTcol[2] + pTcol[1] + pTcol[4]).Mod2());
    // MCTblbl = std::sqrt((ET35 + ET47)*(ET35 + ET47) - (pTcol[0] + pTcol[2] - pTcol[1] - pTcol[4]).Mod2());

    // magnitiude of l+ in collider frame
    // double plepmag = std::sqrt(pcol[2].Px()*pcol[2].Px() + pcol[2].Py()*pcol[2].Py() + pcol[2].Pz()*pcol[2].Pz());

    // // Create basis vectors for the GRRS based on the top direction
    // TVector3 x2, y2, z2;
    // TLorentzVector plepGRRS = pcol[2];
    // z2[1] = 0;
    // z2[2] = 0;
    // z2[3] = 1;
    // double p5xp, p5yp, p5zp;

    // // top direction is x'
    // for (int i = 1; i < 3; i++) {
    //   x2[i] = ptcol[i]/sqrt(ptcol[1]*ptcol[1] + ptcol[2]*ptcol[2]);
    // }
    // x2[3] = 0;

    // // y' obtained using cross product
    // y2[1] = z2[2]*x2[3] - z2[3]*x2[2];
    // y2[2] = z2[3]*x2[1] - z2[1]*x2[3];
    // y2[3] = z2[1]*x2[2] - z2[2]*x2[1];

    // double plepGRRSx, plepGRRSy, plepGRRSz;
    // for (int i = 1; i < 4; i++) {
    //   p5xp = p5xp + pcol[3][i]*x2[i];
    //   p5yp = p5yp + pcol[3][i]*y2[i];
    //   p5zp = p5zp + pcol[3][i]*z2[i];
    // }

    // double plepmagGRRS = std::sqrt(p5xp*p5xp + p5yp*p5yp + p5zp*p5zp);

    // if (std::abs(plepmag - plepmagGRRS) >= 1e-11) printf("Error: coordinate transform mismatch.\n");

    // phi_l = atan2(p5yp,p5xp)

    // if (phi_l < 0)phi_l = phi_l + 2*pi
    // call rootadddouble(phi_l, "fl")

    // cosfl = cos(phi_l)
    // call rootadddouble(cosfl, "cosphil")
    //         end if
  }
}

void KinematicsZprime::PostLoop()
{

}

double KinematicsZprime::deltaPhi(const double& phi1,const double& phi2) const
{
  double dPhi = std::fabs(phi1 - phi2);
  // if (dPhi > m_pi) dPhi = 2.*m_pi - dPhi;
  return dPhi;
}

void KinematicsZprime::PreLoop()
{
  this->SetupInputFiles();
  this->SetupOutputFiles();
}

void KinematicsZprime::Loop()
{
  // Loop over all files
  for(Itr_s i = m_inputFiles->begin(); i != m_inputFiles->end(); ++i)
  {
    cout<<"Processing File = "<<(*i)<<endl;
    
    this->SetupTreesForNewFile((*i));
 
    // The Event Loop
    Long64_t nEvents = this->TotalEvents();
    printf("--- Start of event loop. ---\n");
    for(Long64_t jentry=0; jentry<nEvents;++jentry) 
    {
      // printf("Processing entry %lli\n", jentry);
      Long64_t ientry = this->IncrementEvent(jentry);
      if (ientry < 0) break;
      this->EachEvent();
    }
    printf("---  End of event loop.  ---\n");
    this->CleanUp();
  }
}

KinematicsZprime::~KinematicsZprime()
{
  delete m_inputFiles;
}

void KinematicsZprime::SetupOutputFiles()
{
  m_outputFile = new TFile(m_outputFileName,"RECREATE");
  printf("Processed TTree will be saved to %s.\n", m_outputFileName.Data());
  TTree *tree = new TTree("T","Processed Z->tt data");
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

}

void KinematicsZprime::SetupInputFiles()
{
  m_inputFiles = new vector<TString>;  
  m_inputFiles->push_back(m_inputFileName);
}

Long64_t KinematicsZprime::TotalEvents()
{
  // Internal for Event Looping
  if (m_ntup != 0){return m_ntup->totalEvents();}
  return -999;  
}

Long64_t KinematicsZprime::IncrementEvent(Long64_t i)
{
  Long64_t ev(-1);
  if (m_ntup != 0){ev = m_ntup->LoadTree(i);}
  return ev;  
}

void KinematicsZprime::SetupTreesForNewFile(const TString& s)
{
  TString treeToUse = "RootTuple";
  
  m_chainNtup = new TChain(treeToUse,"");
  TString TStringNtuple = s + "/" + treeToUse;
  m_chainNtup->Add(TStringNtuple);
  m_ntup = new RootTuple(m_chainNtup);  
}

void KinematicsZprime::CleanUp()
{
  delete m_chainNtup;
  delete m_ntup;
}