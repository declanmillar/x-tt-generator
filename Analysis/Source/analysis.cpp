#include "analysis.h"

AnalysisZprime::AnalysisZprime(const TString channel, const TString model, const TString& inputFileName, const TString& outputFileName) :
  m_pi(3.14159265),
  m_GeV(1000.0),
  m_Wmass(80.23),
  m_channel(channel),
  m_model(model),
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

void AnalysisZprime::EachEvent()
{
  // 0 = b, 1 = bbar, 2 = l+, 3 = nu, 4 = l-, 5 = nubar

  vector<TLorentzVector> pcol(6);
  vector<TLorentzVector> p(6);
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
  TLorentzVector ptcol;
  TLorentzVector ptbcol;
  TLorentzVector pt;
  TLorentzVector ptb;

  // final particle collider variables
  TVector2 pTcoltot;
  for (int i = 0; i < (int) m_ntup->E()->size(); i++) {
    pcol[i].SetPxPyPzE(m_ntup->Px()->at(i), m_ntup->Py()->at(i), m_ntup->Pz()->at(i), m_ntup->E()->at(i));
    p[i].SetPxPyPzE(m_ntup->Px()->at(i), m_ntup->Py()->at(i), m_ntup->Pz()->at(i), m_ntup->E()->at(i));
    pTcol[i].Set(pcol[i].Px(), pcol[i].Py());
    ycol[i] = pcol[i].Rapidity();
    etacol[i] = pcol[i].PseudoRapidity();
    phicol[i] = pcol[i].Phi();   
    mass[i] = pcol[i].M();
    ETcol[i] = sqrt(mass[i]*mass[i] + pTcol[i].Mod2());
    pcoltot += pcol[i];
    pTcoltot += pTcol[i];
  }
  // negative velocity of full system in collider frame
  TVector3 vcoltot = -1*pcoltot.BoostVector();
  double ytt = pcoltot.Rapidity();
  double PzNuReco = -9999;
  double MttPzNuReco =-9999;
  if (m_channel =="2to6") {
    PzNuReco = resolveNeutrinoPz(pcol[2], pT[3]);

    // printf("PzNu = %f\n", pcol[3].Pz());
    // printf("PzNuReco = %f\n", PzNuReco);


    TLorentzVector pNuReco;
    pNuReco.SetPxPyPzE(pcol[3].Px(), pcol[3].Py(), PzNuReco, sqrt(pcol[3].Px()*pcol[3].Px()+pcol[3].Py()*pcol[3].Px()+PzNuReco*PzNuReco));
    TLorentzVector pcoltotPzNuReco = pcoltot - pcol[3] + pNuReco;
    MttPzNuReco = pcoltotPzNuReco.M();
  }

   // final particle parton CoM variables
  TVector2 pTtot;
  for (int i = 0; i < (int) m_ntup->E()->size(); i++) {
    p[i].Boost(vcoltot);
    pT[i].Set(p[i].Px(), p[i].Py());
    y[i] = p[i].Rapidity();
    eta[i] = p[i].PseudoRapidity();
    phi[i] = p[i].Phi();   
    ET[i] = sqrt(mass[i]*mass[i] + pT[i].Mod2());
    ptot += p[i];
    pTtot += pT[i];
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
  }

  double ytcol = ptcol.Rapidity();
  double ytbcol = ptbcol.Rapidity();
  double etatcol = ptcol.PseudoRapidity();
  double etatbcol = ptbcol.PseudoRapidity();
  double phitcol = ptcol.Phi();
  double phitbcol = ptbcol.Phi();
  double pTtcol = ptcol.Pt();
  double PTtbcol = ptbcol.Pt();
  double yt = pt.Rapidity();
  double ytb = ptb.Rapidity();
  double etat = pt.PseudoRapidity();
  double etatb = ptb.PseudoRapidity();
  double phit = pt.Phi();
  double phitb = ptb.Phi();
  double pTt = pt.Pt();
  double PTtb = ptb.Pt();
  double deltaabsycol = std::abs(ytcol) - std::abs(ytbcol);
  double costhetatcol = ptcol.CosTheta();
  double costhetat = pt.CosTheta();
  double costhetastar = int(ytt/std::abs(ytt))*costhetat;
  // printf("costhetastar = %f\n", costhetastar);
  // printf("\nytt             = %f\n", ytt);
  // printf("yt              = %f\n", yt);
  // printf("costhetat       = %f\n", costhetat);
  // printf("ytb              = %f\n", ytb);
  // printf("ytcol-ytt       = %f\n", ytcol-ytt);
  // printf("deltaabsycol    = %f\n", deltaabsycol);
  // printf("ytcol           = %f\n", ytcol);
  // printf("ytbcol           = %f\n", ytbcol);
  // printf("2yt             = %f\n", 2*yt);

  // negative velocity of t/tb in collider frame
  TVector3 vtcol = -1*ptcol.BoostVector();
  TVector3 vtbcol = -1*ptbcol.BoostVector();

  TLorentzVector plepptop;
  TLorentzVector plepmtop;

  int it = m_ntup->iteration();

  // printf("it = %i\n", it);

  double Mtt = std::abs(pcoltot.M());
  // printf("Mtt = %f\n", Mtt);
  double costhetalpcol = -9999;
  double costhetalmcol = -9999;
  double costhetalptop = -9999;
  double costhetalmtop = -9999;
  double clpclmcol = -9999;
  double clpclmtop = -9999;
  double dphi = -9999;
  double MET = -9999;
  double mll = -9999;
	double Mbbll = -9999;
	double HT = -9999;
	double ETbbll = -9999;
	double KTbbll = -9999;
	double ET5 = -9999;
	double ET7 = -9999;
	double MTll = -9999;
	double MCTll = -9999;
	double m35 = -9999;
	double m47 = -9999;
	double ET35 = -9999;
	double ET47 = -9999;
	double MTblbl = -9999;
	double MCTblbl = -9999;

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

    // TRANSVERSE VARIABLES

    TVector2 ETmiss = -1*pTcol[0] - pTcol[1] - pTcol[2] - pTcol[4];
    MET = ETmiss.Mod();

    mll = (pcol[2] + pcol[4]).M();

    // calculate invariant mass of visible decay products
    Mbbll = (pcol[0] + pcol[1] + pcol[2] + pcol[4]).M();

    // calculate total scalar sum of transverse energy
    HT = ETcol[0] + ETcol[1] + ETcol[2] + ETcol[4] + MET;

    // ET of visible decay products
    TVector2 pTbbll = pTcol[0] + pTcol[1] + pTcol[2] + pTcol[4];
    ETbbll = sqrt(Mbbll*Mbbll + pTbbll.Mod2());

    // scalar sum of visible decay products and MET
    KTbbll = ETbbll + MET;

    ET5 = sqrt(mass[2]*mass[2] + pTcol[2].Mod2());
    ET7 = sqrt(mass[4]*mass[4] + pTcol[4].Mod2());

    MTll = sqrt((ET5 + ET7)*(ET5 + ET7) - (pTcol[2] + pTcol[4]).Mod2());
    MCTll = sqrt((ET5 + ET7)*(ET5 + ET7) - (pTcol[2] - pTcol[4]).Mod2());

    m35 = (pcol[0] + pcol[2]).M();
    m47 = (pcol[1] + pcol[3]).M();
    TVector2 pT35 = pTcol[0] + pTcol[2];
    TVector2 pT47 = pTcol[1] + pTcol[4];
    ET35 = sqrt(m35*m35 - (pTcol[0] + pTcol[2]).Mod2());
    ET47 = sqrt(m47*m47 - (pTcol[1] + pTcol[4]).Mod2());

    MTblbl = sqrt((ET35 + ET47)*(ET35 + ET47) - (pTcol[0] + pTcol[2] + pTcol[1] + pTcol[4]).Mod2());
    MCTblbl = sqrt((ET35 + ET47)*(ET35 + ET47) - (pTcol[0] + pTcol[2] - pTcol[1] - pTcol[4]).Mod2());
  }


  if (this->PassCuts())
  {    
    double weight = m_ntup->weight();
    double weight_ee = m_ntup->weight_ee();
    double weight_eq = m_ntup->weight_eq();

    // re-weight for different iterations
    weight = weight*m_sigma/m_cnorm[it-1];
    weight_ee = weight_ee*m_sigma/m_cnorm[it-1];
    weight_eq = weight_eq*m_sigma/m_cnorm[it-1];

    // Fill Histograms (assumes fixed bin width!)
    h_Mtt->Fill(Mtt, weight/h_Mtt->GetXaxis()->GetBinWidth(1));

    if (costhetastar > 0) {
      h_AFstar->Fill(Mtt, weight/h_AFstar->GetXaxis()->GetBinWidth(1));
    }

    if (costhetastar < 0) {
      h_ABstar->Fill(Mtt, weight/h_ABstar->GetXaxis()->GetBinWidth(1));
    }

    if (deltaabsycol > 0) {
      h_RF->Fill(Mtt, weight/h_RF->GetXaxis()->GetBinWidth(1));
    }

    if (deltaabsycol < 0) {
      h_RB->Fill(Mtt, weight/h_RB->GetXaxis()->GetBinWidth(1));
    }

    if (m_channel == "2to6") {
      h_costheta5_eq->Fill(costhetalptop, weight);
      h_costheta5_ee->Fill(costhetalptop, weight);
      h_ct7ct5->Fill(clpclmtop, weight);
      h_dphi->Fill(dphi, weight/h_dphi->GetXaxis()->GetBinWidth(1));
      h_PzNu->Fill(pcol[3].Pz(), weight_eq);
      h_PzNuReco->Fill(PzNuReco, weight_eq);
      h_MttPzNuReco->Fill(MttPzNuReco, weight_eq);
      h_MET->Fill(MET, weight_ee);
			h_HT->Fill(HT, weight_ee);
			h_Mbbll->Fill(Mbbll, weight_ee);
			h_mll->Fill(mll, weight_ee);
			h_ETbbll->Fill(ETbbll, weight_ee);
			h_KTbbll->Fill(KTbbll, weight_ee);
			h_MTll->Fill(MTll, weight_ee);
			h_MCTll->Fill(MCTll, weight_ee);
			h_MTblbl->Fill(MTblbl, weight_ee);
			h_MCTblbl->Fill(MCTblbl, weight_ee);
      h_dphi_HT->Fill(dphi, HT, weight_ee);
      h_dphi_Mbbll->Fill(dphi, Mbbll, weight_ee);
      h_dphi_mll->Fill(dphi, mll, weight_ee);
      h_dphi_ETbbll->Fill(dphi, ETbbll, weight_ee);
      h_dphi_KTbbll->Fill(dphi, KTbbll, weight_ee);
      h_dphi_MTll->Fill(dphi, MTll, weight_ee);
      h_dphi_MCTll->Fill(dphi, MCTll, weight_ee);
      h_dphi_MTblbl->Fill(dphi, MTblbl, weight_ee);
      h_dphi_MCTblbl->Fill(dphi, MCTblbl, weight_ee);

      if (costhetalptop > 0) {
        h_AlLF->Fill(Mtt, weight/h_AlLF->GetXaxis()->GetBinWidth(1));
      }

      if (costhetalptop < 0) {
        h_AlLB->Fill(Mtt, weight/h_AlLB->GetXaxis()->GetBinWidth(1));
      }

      if (clpclmtop > 0) {
        h_ALLF->Fill(Mtt, weight/h_ALLF->GetXaxis()->GetBinWidth(1));
      }

      if (clpclmtop < 0) {
        h_ALLB->Fill(Mtt, weight/h_ALLB->GetXaxis()->GetBinWidth(1));
      }
    }

    if (m_channel == "2to2") {
      h_MttLL->Fill(Mtt, m_ntup->weightLL()/h_MttLL->GetXaxis()->GetBinWidth(1));
      h_MttLR->Fill(Mtt, m_ntup->weightLR()/h_MttLR->GetXaxis()->GetBinWidth(1));
      h_MttRL->Fill(Mtt, m_ntup->weightRL()/h_MttRL->GetXaxis()->GetBinWidth(1));
      h_MttRR->Fill(Mtt, m_ntup->weightRR()/h_MttRR->GetXaxis()->GetBinWidth(1));
    }    
  }
}

void AnalysisZprime::PostLoop()
{
  double sigma = h_Mtt->Integral("width");
  if (sigma != m_sigma) {
    printf("Cross section from generation and analysis stage do not match!\n");
    printf("sigma_generation = %f\n", m_sigma);
    printf("sigma_analysis   = %f\n", sigma);
  }
  else printf("Total cross section = %f\n", sigma);
  if (m_channel == "2to2") {
    h_MttALL = this->MttALL();
    h_MttAL = this->MttAL();
    this->TotalSpinAsymmetries();
  }

  h_AFBstar = this->Asymmetry("AFBstar", "A^{*}_{FB}", h_AFstar, h_ABstar);
  h_ARFB = this->Asymmetry("ARFB", "A_{RFB}", h_RF, h_RB);

  double AFBstar = this->TotalAsymmetry(h_AFstar,h_ABstar);
  double ARFB = this->TotalAsymmetry(h_RF,h_RB);
  double ALL = this->TotalAsymmetry(h_ALLF,h_ALLB);
  double AL = this->TotalAsymmetry(h_AlLF,h_AlLF);
  printf("AFBstar = %f\n", AFBstar);
  printf("ARFB = %f\n", ARFB);
  printf("ALL = %f\n", ALL);
  printf("AL = %f\n", AL);

  this->MakeGraphs();

  if (m_channel == "2to6") {
    this->ALL2to6();
    h_AL = this->Asymmetry("AL", "A_{L}", h_AlLF, h_AlLB);
    h_ALL = this->Asymmetry("ALL", "A_{LL}", h_ALLF, h_ALLB);
  }
  
  this->WriteHistograms();
}

double AnalysisZprime::deltaPhi(const double& phi1,const double& phi2) const
{
  double dPhi = std::fabs(phi1 - phi2);
  // if (dPhi > m_pi) dPhi = 2.*m_pi - dPhi;
  return dPhi;
}

TH1D* AnalysisZprime::MttALL()
{
  TH1D* h_A = (TH1D*) h_MttLL->Clone();
  TH1D* h_B = (TH1D*) h_MttLR->Clone();

  h_A->Add(h_MttRR);
  h_B->Add(h_MttRL);

  double ALL = this->TotalAsymmetry(h_A,h_B);
  printf("ALL' = %f\n", ALL);

  TH1D* h_ALL = this->Asymmetry("ALL", "A_{LL}", h_A, h_B);
  delete h_A;
  delete h_B;
  return h_ALL;
}

TH1D* AnalysisZprime::MttAL()
{
  TH1D* h_A = (TH1D*) h_MttLL->Clone();
  TH1D* h_B = (TH1D*) h_MttRR->Clone();

  h_A->Add(h_MttLR);
  h_B->Add(h_MttRL);

  TH1D* h_AL = this->Asymmetry("AL", "A_{L}", h_A, h_B);
  delete h_A;
  delete h_B;
  return h_AL;
}

void AnalysisZprime::TotalSpinAsymmetries()
{
  double sigma = h_Mtt->Integral("width");
  double sigmaLL = h_MttLL->Integral("width");
  double sigmaLR = h_MttLR->Integral("width");
  double sigmaRL = h_MttRL->Integral("width");
  double sigmaRR = h_MttRR->Integral("width");

  double ALL = (sigmaLL + sigmaRR - sigmaRL - sigmaLR)/
               (sigmaLL + sigmaRR + sigmaRL + sigmaLR);

  double AL = (sigmaRR + sigmaRL - sigmaLR - sigmaLL)/
              (sigmaLL + sigmaRR + sigmaRL + sigmaLR);

  printf("sigma = %f\n", sigma);
  printf("ALL = %f\n", ALL);
  printf("AL = %f\n", AL);  
}

double AnalysisZprime::TotalAsymmetry(TH1D* h_A, TH1D* h_B)
{
  double A = h_A->Integral("width");
  double B = h_B->Integral("width");
  double Atot = (A - B)/(A + B);
  return Atot;
}

void AnalysisZprime::ALL2to6() 
{
  double mean = h_ct7ct5->GetMean();

  double ALL = -9*mean;

  printf("ALL = %f\n", ALL);
}

TH1D* AnalysisZprime::Asymmetry(TString name, TString title, TH1D* h_A, TH1D* h_B)
{
  TH1D* h_numerator = (TH1D*) h_A->Clone(name);
  TH1D* h_denominator = (TH1D*) h_A->Clone();
  h_numerator->SetTitle(title);
  h_numerator->Add(h_B, -1);
  h_denominator->Add(h_B, 1);
  h_numerator->Divide(h_denominator);
  delete h_denominator;
  return h_numerator;
}

void AnalysisZprime::CreateHistograms() 
{

  h_Mtt = new TH1D("Mtt", "M_{tt}", 100, 0.0, 13000.0);
  h_AFstar = new TH1D("AFstar", "AFstar", 50, 0.0, 13000.0);
  h_ABstar = new TH1D("ABstar", "ABstar", 50, 0.0, 13000.0);
  h_RF = new TH1D("RF", "RF", 50, 0.0, 13000.0);
  h_RB = new TH1D("RB", "RB", 50, 0.0, 13000.0);


  if (m_channel == "2to6") {
    h_PzNuReco = new TH1D("PzNuReco", "p_{z}^{#nu} (reconstructed)", 100, -1000.0, 1000.0);
    h_MttPzNuReco = new TH1D("MttPzNuReco", "M_{tt} (p_{z}^{#nu} reconstructed)", 100, 0.0, 13000.0);
    h_PzNu = new TH1D("PzNu", "p_{z}^{#nu}", 100,-1000.0, 1000.0);
    h_costheta5_eq = new TH1D("costheta5_eq", "cos#theta_{l^{+}} (eq)", 50, -1.0, 1.0);
    h_costheta5_ee = new TH1D("costheta5_ee", "cos#theta_{l^{+}} (ee)", 50, -1.0, 1.0);
    h_ct7ct5 = new TH1D("ct7ct5", "cos#theta_{l^{+}}cos#theta_{l^{-}}", 50, -1.0, 1.0);
    h_dphi = new TH1D("dphi", "dphi", 100, 0, 2*m_pi);
    h_MET = new TH1D("MET", "MET", 20, 0, 4000);
    h_HT = new TH1D("HT", "H_{T}", 20, 0, 4000);
    h_Mbbll = new TH1D("Mbbll", "Mbbll", 20, 0, 4000);
    h_mll = new TH1D("mll", "mll", 20, 0, 4000);
    h_ETbbll = new TH1D("ETbbll", "ETbbll", 20, 0, 4000);
    h_KTbbll = new TH1D("KTbbll", "K_{T}^{bbll}", 20, 0, 4000);
    h_MTll = new TH1D("MTll", "MTll", 20, 0, 4000);
    h_MCTll = new TH1D("MCTll", "MCTll", 20, 0, 4000);
    h_MTblbl = new TH1D("MTblbl", "MTblbl", 20, 0, 4000);
    h_MCTblbl = new TH1D("MCTblbl", "MCTblbl", 20, 0, 4000);
    h_dphi_HT = new TH2D("dphi_HT", "dphi_HT", 20, 0, 2*m_pi, 20, 0, 4000);
    h_dphi_Mbbll = new TH2D("dphi_Mbbll", "dphi_Mbbll", 20, 0, 2*m_pi, 20, 0, 4000);
    h_dphi_mll = new TH2D("dphi_mll", "dphi_mll", 20, 0, 2*m_pi, 20, 0, 4000);
    h_dphi_ETbbll = new TH2D("dphi_ETbbll", "dphi_ETbbll", 20, 0, 2*m_pi, 20, 0, 4000);
    h_dphi_KTbbll = new TH2D("dphi_KTbbll", "dphi_KTbbll", 20, 0, 2*m_pi, 20, 0, 4000);
    h_dphi_MTll = new TH2D("dphi_MTll", "dphi_MTll", 20, 0, 2*m_pi, 20, 0, 4000);
    h_dphi_MCTll = new TH2D("dphi_MCTll", "dphi_MCTll", 20, 0, 2*m_pi, 20, 0, 4000);
    h_dphi_MTblbl = new TH2D("dphi_MTblbl", "dphi_MTblbl", 20, 0, 2*m_pi, 20, 0, 4000);
    h_dphi_MCTblbl = new TH2D("dphi_MCTblbl", "dphi_MCTblbl", 20, 0, 2*m_pi, 20, 0, 4000);
    h_AlLF = new TH1D("AlLF", "AlLF", 50, 0.0, 13000.0);
    h_AlLB = new TH1D("AlLB", "AlLB", 50, 0.0, 13000.0);
    h_ALLF = new TH1D("ALLF", "ALLF", 50, 0.0, 13000.0);
    h_ALLB = new TH1D("ALLB", "ALLB", 50, 0.0, 13000.0);
  }

  if (m_channel == "2to2") {
    h_MttLL = new TH1D("MttLL", "MttLL", 130, 0.0, 13000.0);
    h_MttLR = new TH1D("MttLR", "MttLR", 130, 0.0, 13000.0);
    h_MttRL = new TH1D("MttRL", "MttRL", 130, 0.0, 13000.0);
    h_MttRR = new TH1D("MttRR", "MttRR", 130, 0.0, 13000.0);
  }
}

void AnalysisZprime::MakeGraphs()
{
  printf("Making Graphs...\n");
  TString numBase;
  if (m_channel == "2to2") numBase = "d#sigma(pp->t#bar{t}) / d"; 
  if (m_channel == "2to6") numBase = "d#sigma / d"; //pp->t#bar{t}->b#bar{b}l^{+}l^{-}#nu#bar{#nu}
  TString units = "pb";

  // TCanvas *c_Mtt   = new TCanvas(h_Mtt->GetName(), h_Mtt->GetTitle());
  // c_Mtt->cd(); 
  // h_Mtt->Draw("hist"); 
  h_Mtt->GetXaxis()->SetTitle(h_Mtt->GetTitle());
  h_Mtt->GetYaxis()->SetTitle(numBase + h_Mtt->GetTitle() + " [" + units +"/GeV]");

  // TCanvas *c_AFBstar   = new TCanvas(h_AFBstar->GetName(), h_AFBstar->GetTitle());
  // c_AFBstar->cd(); 
  // h_AFBstar->Draw("hist"); 
  // h_AFBstar->GetYaxis()->SetTitle("");

  if (m_channel == "2to6") {
    // TCanvas *c_pz5 = new TCanvas(h_pz5->GetName(), h_pz5->GetTitle());
    // c_pz5->cd();
    // h_pz5->Draw("hist"); 
    // h_pz5->GetYaxis()->SetTitle(numBase + h_pz5->GetTitle() + " [" + units +"/GeV]");

    // TCanvas *c_costheta5_eq = new TCanvas(h_costheta5_eq->GetName(), h_costheta5_eq->GetTitle());
    // c_costheta5_eq->cd(); 
    // h_costheta5_eq->Draw("hist"); 
    // h_costheta5_eq->GetYaxis()->SetTitle(numBase + h_costheta5_eq->GetTitle() + " [" + units +"]");

    // TCanvas *c_costheta5_ee = new TCanvas(h_costheta5_ee->GetName(), h_costheta5_ee->GetTitle());
    // c_costheta5_ee->cd(); 
    // h_costheta5_ee->Draw("hist"); 
    // h_costheta5_ee->GetYaxis()->SetTitle(numBase + h_costheta5_ee->GetTitle() + " [" + units +"]");

    // TCanvas *c_ct7ct5 = new TCanvas(h_ct7ct5->GetName(), h_ct7ct5->GetTitle());
    // c_ct7ct5->cd(); 
    // h_ct7ct5->Draw("hist"); 
    // h_ct7ct5->GetYaxis()->SetTitle(numBase + h_ct7ct5->GetTitle() + " [" + units +"]");

    // TCanvas *c_KTbbll = new TCanvas(h_KTbbll->GetName(), h_KTbbll->GetTitle());
    // c_KTbbll->cd(); 
    // h_KTbbll->Draw("hist"); 
    h_KTbbll->GetYaxis()->SetTitle(numBase + h_KTbbll->GetTitle() + " [" + units +"/GeV]");
    h_KTbbll->GetXaxis()->SetTitle(h_KTbbll->GetTitle());

    h_HT->GetYaxis()->SetTitle(numBase + h_HT->GetTitle() + " [" + units +"/GeV]");
    h_HT->GetXaxis()->SetTitle(h_HT->GetTitle());

    h_dphi_KTbbll->GetYaxis()->SetTitle("K_{T}^{bbll}");
    h_dphi_KTbbll->GetXaxis()->SetTitle("#Delta#phi");

    h_dphi_HT->GetYaxis()->SetTitle("H_{T}");
    h_dphi_HT->GetXaxis()->SetTitle("#Delta#phi");

    h_PzNu->GetYaxis()->SetTitle(numBase + h_PzNu->GetTitle() + " [" + units +"/GeV]");
    h_PzNu->GetXaxis()->SetTitle(h_PzNu->GetTitle());

    h_PzNuReco->GetYaxis()->SetTitle(numBase + h_PzNuReco->GetTitle() + " [" + units +"/GeV]");
    h_PzNuReco->GetXaxis()->SetTitle(h_PzNuReco->GetTitle());

    h_MttPzNuReco->GetYaxis()->SetTitle(numBase + h_MttPzNuReco->GetTitle() + " [" + units +"/GeV]");
    h_MttPzNuReco->GetXaxis()->SetTitle(h_MttPzNuReco->GetTitle());

  }

  if (m_channel == "2to2") {
    h_MttLL->GetYaxis()->SetTitle(numBase + h_MttLL->GetTitle() + " [" + units +"/GeV]");
    h_MttLR->GetYaxis()->SetTitle(numBase + h_MttLR->GetTitle() + " [" + units +"/GeV]");
    h_MttRL->GetYaxis()->SetTitle(numBase + h_MttRL->GetTitle() + " [" + units +"/GeV]");
    h_MttRR->GetYaxis()->SetTitle(numBase + h_MttRR->GetTitle() + " [" + units +"/GeV]");
    h_MttALL->GetYaxis()->SetTitle(h_MttALL->GetTitle());
    h_MttAL->GetYaxis()->SetTitle(h_MttAL->GetTitle());
    h_AFBstar->GetYaxis()->SetTitle(h_AFBstar->GetTitle());
    h_ARFB->GetYaxis()->SetTitle(h_ARFB->GetTitle());
    h_MttLL->GetXaxis()->SetTitle("M_{tt} [GeV]");
    h_MttLR->GetXaxis()->SetTitle("M_{tt} [GeV]");
    h_MttRL->GetXaxis()->SetTitle("M_{tt} [GeV]");
    h_MttRR->GetXaxis()->SetTitle("M_{tt} [GeV]");
    h_MttALL->GetXaxis()->SetTitle("M_{tt} [GeV]");
    h_MttAL->GetXaxis()->SetTitle("M_{tt} [GeV]");
    h_AFBstar->GetXaxis()->SetTitle("M_{tt} [GeV]");
    h_ARFB->GetXaxis()->SetTitle("M_{tt} [GeV]");
  }
  printf("...complete.\n");
}

void AnalysisZprime::WriteHistograms() 
{
  m_outputFile->cd();
  m_outputFile->cd("/");

  printf("Writing histograms...\n");

  // // Save histograms
  h_Mtt->Write();
  h_AFBstar->Write();
  h_ARFB->Write();

  if (m_channel == "2to2") {
    h_MttLL->Write();
    h_MttLR->Write();
    h_MttRL->Write();
    h_MttRR->Write();
    h_MttALL->Write();
    h_MttAL->Write();
  }

  if (m_channel == "2to6") {
    h_dphi->Write();
    h_AL->Write();
    h_ALL->Write();
    h_PzNu->Write();
    h_PzNuReco->Write();
    h_MttPzNuReco->Write();
    h_costheta5_eq->Write();
    h_costheta5_ee->Write();
    h_ct7ct5->Write();
    h_MET->Write();
    h_HT->Write();
    h_Mbbll->Write();
    h_mll->Write();
    h_ETbbll->Write();
    h_KTbbll->Write();
    h_MTll->Write();
    h_MCTll->Write();
    h_MTblbl->Write();
    h_MCTblbl->Write();
    h_dphi_HT->Write();
    h_dphi_Mbbll->Write();
    h_dphi_mll->Write();
    h_dphi_ETbbll->Write();
    h_dphi_KTbbll->Write();
    h_dphi_MTll->Write();
    h_dphi_MCTll->Write();
    h_dphi_MTblbl->Write();
    h_dphi_MCTblbl->Write();
  }

  printf("...complete\n");
  m_outputFile->Close();
  delete m_outputFile;
}

bool AnalysisZprime::PassCuts() const
{
  if (this->PassCuts_Mtt())
  {
    if (this->PassCuts_MET())
      {
        return true;
      }
  }
  return false;
}

bool AnalysisZprime::PassCuts_MET() const
{
  if (m_channel == "2to2")
  {
    return true;
  }
  if (m_channel == "2to6")
  {
    // if (m_ntup->Etmiss() > 0.0 * m_GeV)
    // {
      return true;
    // }
  }  
  return false;
}

bool AnalysisZprime::PassCuts_Mtt() const
{
  // if (m_ntup->Mtt() > 1200.0)
  //   {
  //     if (m_ntup->Mtt() < 2200.0)
  //       {
          return true;
  //       }
  //   }
  // return false;
}

void AnalysisZprime::PreLoop()
{
  TString cnormsName(m_inputFileName + ".txt");
  ifstream cnorms(cnormsName.Data());
  if (!cnorms.is_open()) printf("Failed to open %s\n", cnormsName.Data());
  cnorms >> m_sigma;
  string line;
  double cnorm;
  while ( cnorms >> cnorm ) m_cnorm.push_back(cnorm);
  cnorms.close();

  this->SetupInputFiles();
  this->SetupOutputFiles();
  this->CreateHistograms();
}

void AnalysisZprime::Loop()
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

AnalysisZprime::~AnalysisZprime()
{

  delete m_inputFiles;
}

void AnalysisZprime::SetupOutputFiles()
{
  m_outputFile = new TFile(m_outputFileName,"RECREATE");
  printf("Histograms will be saved to %s.\n", m_outputFileName.Data());
}

void AnalysisZprime::SetupInputFiles()
{
  m_inputFiles = new vector<TString>;  
  m_inputFiles->push_back(m_inputFileName);
}

Long64_t AnalysisZprime::TotalEvents()
{
  // Internal for Event Looping
  if (m_ntup != 0){return m_ntup->totalEvents();}
  return -9999;  
}

Long64_t AnalysisZprime::IncrementEvent(Long64_t i)
{
  Long64_t ev(-1);
  if (m_ntup != 0){ev = m_ntup->LoadTree(i);}
  return ev;  
}

void AnalysisZprime::SetupTreesForNewFile(const TString& s)
{
  TString treeToUse = "RootTuple";
  
  m_chainNtup = new TChain(treeToUse,"");
  TString TStringNtuple = s + "/" + treeToUse;
  m_chainNtup->Add(TStringNtuple);
  m_ntup = new RootTuple(m_chainNtup);  
}

void AnalysisZprime::CleanUp()
{
  delete m_chainNtup;
  delete m_ntup;
}

double AnalysisZprime::resolveNeutrinoPz(TLorentzVector p_l, TVector2 pT_nu) {

  // finds the longitudinal neutrino momentum for semi-hadronic decay
  // assuming all particles are massless

  double PzNu;
  std::vector<std::complex<double> > root, root2;
  double a = -9999, b = -9999, c = -9999, k = -9999;
  // double a2 = -9999, b2 = -9999, c2 = -9999;

  // recalculate lepton energy in zero mass approximation
  double p_l0 = sqrt(p_l.Px()*p_l.Px() + p_l.Py()*p_l.Py() + p_l.Pz()*p_l.Pz());

  // check this matches the 4-vector energy
  // printf("p_l.E() = %f\n", p_l.E());
  // printf("p_l0 = %f\n", p_l0);
  if ( std::abs(p_l0 - p_l.E() > 0.00001)) printf("p_l0 doesn't match\n");

  // // Alternative calculation from Ruth
  // double lx = p_l.Px();
  // double ly = p_l.Py();
  // double lz = p_l.Pz();
  // double nx = pT_nu.Px();
  // double ny = pT_nu.Py();
  // double nz = 0.;

  // a2 = (4.*lx*lx
  //    +4.*ly*ly );
  // b2 = (-4.*m_Wmass*m_Wmass*lz
  //     -8.*lx*nx*lz
  //     -8.*ly*ny*lz );
  // c2  = ( 4.*ly*ly*nx*nx
  //      + 4.*lz*lz*nx*nx
  //      + 4.*lx*lx*ny*ny
  //      + 4.*lz*lz*ny*ny
  //      - 8.*lx*nx*ly*ny
  //      - 4.*lx*nx*m_Wmass*m_Wmass
  //      - 4.*ly*ny*m_Wmass*m_Wmass
  //      - m_Wmass*m_Wmass*m_Wmass*m_Wmass );

  k = m_Wmass*m_Wmass/2 + p_l.Px()*pT_nu.Px() + p_l.Py()*pT_nu.Py();

  a = p_l.Px()*p_l.Px() + p_l.Py()*p_l.Py();

  b = -2*k*(p_l.Pz());

  c = (pT_nu.Px()*pT_nu.Px() + pT_nu.Py()*pT_nu.Py())*p_l.E()*p_l.E() - k*k;

  // printf("a = %f, %f\n", a, a2);
  // printf("b = %f, %f\n", b, b2);
  // printf("c = %f, %f\n", c, c2);

  root = this->solveQuadratic(a, b, c);
  // root2 = this->solveQuadratic(a2, b2, c2);

  // for (int i = 0; i < 2; i++) cout << "root " << root[i] << ", " << root2[i] << endl;

  // select single solution
  if (root[0].imag() == 0 and root[1].imag() == 0) {
    // two real solutions - pick smallest one
    if (std::abs(root[0].real()) < std::abs(root[1].real())) {
      // solution 1 < than solution 2
      PzNu = root[0].real();
    }
    else if (std::abs(root[0].real()) > std::abs(root[1].real())) { 
      // solution 1 > than solution 2
      PzNu = root[1].real();
    }
    else {
      // solutions are equal pick 1
      PzNu = root[0].real();
    }
  }
  else {
    // no real solutions - take the real part of 1
    PzNu = root[0].real();
  }

  return PzNu;

}

std::vector<std::complex<double> > AnalysisZprime::solveQuadratic(double a, double b, double c) {
    // solves quadratic for both roots 
    // returns both as complex values in a complex vector x(2)

    std::vector<std::complex<double> > roots;
    std::complex<double> term1;
    std::complex<double> term2;
    std::complex<double> discriminator;

    term1 = -b/(2*a);
    discriminator = b*b - 4*a*c;
    term2 = sqrt(discriminator)/(2*a);

    roots.push_back(term1 + term2);
    roots.push_back(term1 - term2);
    
    return roots;
    
}