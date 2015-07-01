#include "analysis.h"

AnalysisZprime::AnalysisZprime(const TString channel, const TString model, const TString& inputFileName, const TString& weightsFileName, const TString& outputFileName) :
  m_pi(3.14159265),
  m_GeV(1000.0),
  m_intLumi(0),//(300000.0),
  m_Wmass(80.23),
  m_channel(channel),
  m_model(model),
  m_inputFileName(inputFileName),
  m_weightsFileName(weightsFileName),
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
  UpdateCutflow(cutEvent,true);
  pcol = vector<TLorentzVector>(6);
  p = vector<TLorentzVector>(6);
  pcolReco = vector<TLorentzVector>(6);
  pReco = vector<TLorentzVector>(6);

  TLorentzVector pcoltot, ptot;
  vector<double> mass(6);
  vector<double> ETcol(6);
  pTcol = vector<TVector2>(6);
  ycol = vector<double>(6);
  etacol = vector<double>(6);
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
  if (m_channel =="bbllnn") {
    // resolve longitudinal neutrino momentum in the semihadronic case
    PzNuRecoCol = ResolveNeutrinoPz(pcol[2], pTcol[3]);
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
  if (m_channel == "tt" or "ll") {
    ptcol = pcol[0];
    ptbcol = pcol[1];
    pt = p[0];
    ptb = p[1];
  }
  else if (m_channel == "bbllnn") {
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
  Mff = std::abs(pcoltot.M());
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

  if (m_channel == "bbllnn") {
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


  if (this->PassCuts())
  {    
    double weight = m_ntup->weight();

    // re-weight for different iterations
    weight = weight*m_sigma/m_cnorm[it-1];

    Mff = Mff/1000; // convert to TeV

    // Fill Histograms (assumes fixed bin width!)
    h_Mtt->Fill(Mff, weight/h_Mtt->GetXaxis()->GetBinWidth(1));
    h_ytt->Fill(ytt, weight/h_ytt->GetXaxis()->GetBinWidth(1));
    h_CosTheta->Fill(CosTheta, weight/h_CosTheta->GetXaxis()->GetBinWidth(1));
    h_CosThetaStar->Fill(CosThetaStar, weight/h_CosThetaStar->GetXaxis()->GetBinWidth(1));

    // asymmetries
    if (CosThetaStar > 0) {
      h_AFBstar1->Fill(Mff, weight/h_AFBstar1->GetXaxis()->GetBinWidth(1));
    }

    if (CosThetaStar < 0) {
      h_AFBstar2->Fill(Mff, weight/h_AFBstar2->GetXaxis()->GetBinWidth(1));
    }

    if (deltaYcol > 0) {
      h_AttC1->Fill(Mff, weight/h_AttC1->GetXaxis()->GetBinWidth(1));
    }

    if (deltaYcol < 0) {
      h_AttC2->Fill(Mff, weight/h_AttC2->GetXaxis()->GetBinWidth(1));
    }

    if (m_channel == "tt") {
      h_MttLL->Fill(Mff, m_ntup->weightLL()/h_MttLL->GetXaxis()->GetBinWidth(1));
      h_MttLR->Fill(Mff, m_ntup->weightLR()/h_MttLR->GetXaxis()->GetBinWidth(1));
      h_MttRL->Fill(Mff, m_ntup->weightRL()/h_MttRL->GetXaxis()->GetBinWidth(1));
      h_MttRR->Fill(Mff, m_ntup->weightRR()/h_MttRR->GetXaxis()->GetBinWidth(1));
    }    
    else if (m_channel == "bbllnn") {
      h_yttReco->Fill(yttReco, weight/h_yttReco->GetXaxis()->GetBinWidth(1));
      h_CosThetaReco->Fill(CosThetaReco, weight/h_CosThetaReco->GetXaxis()->GetBinWidth(1));
      h_CosThetaStarReco->Fill(CosThetaStarReco, weight/h_CosThetaStarReco->GetXaxis()->GetBinWidth(1));
      h_ct7ct5->Fill(clpclmtop, weight/h_ct7ct5->GetXaxis()->GetBinWidth(1));
      h_dphi->Fill(dphi, weight/h_dphi->GetXaxis()->GetBinWidth(1)/h_dphi->GetXaxis()->GetBinWidth(1));
      h_PzNu->Fill(pcol[3].Pz(), weight/h_PzNu->GetXaxis()->GetBinWidth(1));
      h_PzNuReco->Fill(PzNuRecoCol, weight/h_PzNuReco->GetXaxis()->GetBinWidth(1));
      h_MttReco->Fill(MttReco, weight/h_MttReco->GetXaxis()->GetBinWidth(1));
   //    h_MET->Fill(MET, weight/h_MET->GetXaxis()->GetBinWidth(1));
			// h_HT->Fill(HT, weight/h_HT->GetXaxis()->GetBinWidth(1));
			// h_Mbbll->Fill(Mbbll, weight/h_Mbbll->GetXaxis()->GetBinWidth(1));
			// h_mll->Fill(mll, weight/h_mll->GetXaxis()->GetBinWidth(1));
			// h_ETbbll->Fill(ETbbll, weight/h_ETbbll->GetXaxis()->GetBinWidth(1));
			// h_KTbbll->Fill(KTbbll, weight/h_KTbbll->GetXaxis()->GetBinWidth(1));
			// h_MTll->Fill(MTll, weight/h_MTll->GetXaxis()->GetBinWidth(1));
			// h_MCTll->Fill(MCTll, weight/h_MCTll->GetXaxis()->GetBinWidth(1));
			// h_MTblbl->Fill(MTblbl, weight/h_MTblbl->GetXaxis()->GetBinWidth(1));
			// h_MCTblbl->Fill(MCTblbl, weight/h_MCTblbl->GetXaxis()->GetBinWidth(1));
      // h_dphi_HT->Fill(dphi, HT, weight/h_dphi_HT->GetXaxis()->GetBinWidth(1));
      // h_dphi_Mbbll->Fill(dphi, Mbbll, weight/h_dphi_Mbbll->GetXaxis()->GetBinWidth(1));
      // h_dphi_mll->Fill(dphi, mll, weight/h_dphi_mll->GetXaxis()->GetBinWidth(1));
      // h_dphi_ETbbll->Fill(dphi, ETbbll, weight/h_dphi_ETbbll->GetXaxis()->GetBinWidth(1));
      // h_dphi_KTbbll->Fill(dphi, KTbbll, weight/h_dphi_KTbbll->GetXaxis()->GetBinWidth(1));
      // h_dphi_MTll->Fill(dphi, MTll, weight/h_dphi_MTll->GetXaxis()->GetBinWidth(1));
      // h_dphi_MCTll->Fill(dphi, MCTll, weight/h_dphi_MCTll->GetXaxis()->GetBinWidth(1));
      // h_dphi_MTblbl->Fill(dphi, MTblbl, weight/h_dphi_MTblbl->GetXaxis()->GetBinWidth(1));
      // h_dphi_MCTblbl->Fill(dphi, MCTblbl, weight/h_dphi_MCTblbl->GetXaxis()->GetBinWidth(1));


      // asymmetries
      if (costhetalptop > 0) {
        h_AlLF->Fill(Mff, weight/h_AlLF->GetXaxis()->GetBinWidth(1));
      }

      if (costhetalptop < 0) {
        h_AlLB->Fill(Mff, weight/h_AlLB->GetXaxis()->GetBinWidth(1));
      }

      if (clpclmtop > 0) {
        h_ALLF->Fill(Mff, weight/h_ALLF->GetXaxis()->GetBinWidth(1));
      }

      if (clpclmtop < 0) {
        h_ALLB->Fill(Mff, weight/h_ALLB->GetXaxis()->GetBinWidth(1));
      }

      if (deltaEtal > 0) {
        h_AllCF->Fill(Mff, weight/h_AllCF->GetXaxis()->GetBinWidth(1));
      }

      if (deltaEtal < 0) {
        h_AllCB->Fill(Mff, weight/h_AllCB->GetXaxis()->GetBinWidth(1));
      }

      if (CosThetaStarReco > 0) {
        h_AFBstarReco1->Fill(MttReco, weight/h_AFBstarReco1->GetXaxis()->GetBinWidth(1));
      }

      if (CosThetaStarReco < 0) {
        h_AFBstarReco2->Fill(MttReco, weight/h_AFBstarReco2->GetXaxis()->GetBinWidth(1));
      }
    }
  }
}

void AnalysisZprime::PostLoop()
{
  printf("------\n");
  double sigma = h_Mtt->Integral("width");
  if (std::abs(sigma - m_sigma) > 10e-11) {
    printf("Cross section from generation and analysis stage do not match!\n");
    printf("sigma_generation = %f\n", m_sigma);
    printf("sigma_analysis   = %f\n", sigma);
  }
  else printf("sigma = %f [pb]\n", sigma);
  if (m_channel == "tt") {
    h_ALL = this->MttALL();
    h_AL = this->MttAL();
    this->TotalSpinAsymmetries();
  }
  printf("------\n");

  h_AFBstar = this->Asymmetry("AFBstar", "A^{*}_{FB}", h_AFBstar1, h_AFBstar2);
  h_AttC = this->Asymmetry("AttC", "A_{C}", h_AttC1, h_AttC2);

  // double AFBstar = this->TotalAsymmetry(h_AFBstar1,h_AFBstar2);
  // double AttC = this->TotalAsymmetry(h_AttC1,h_AttC2);
  // double ALL = this->TotalAsymmetry(h_ALLF,h_ALLB);
  // double AL = this->TotalAsymmetry(h_AlLF,h_AlLF);
  // printf("AFBstar = %f\n", AFBstar);
  // printf("AttC = %f\n", AttC);
  // printf("ALL = %f\n", ALL);
  // printf("AL = %f\n", AL);

  if (m_channel == "bbllnn") {
    this->ALL2to6();
    h_AL = this->Asymmetry("AL", "A_{L}", h_AlLF, h_AlLB);
    h_ALL = this->Asymmetry("ALL", "A_{LL}", h_ALLF, h_ALLB);
    h_AllC = this->Asymmetry("AllC", "A^{ll}_{C}", h_AllCF, h_AllCB);
    h_AFBstarReco = this->Asymmetry("AFBstarNuReco", "A_{FB}^* (#nu reconstructed)", h_AFBstarReco1, h_AFBstarReco2);
  }
  this->MakeGraphs();
  this->PrintCutflow();
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
  // note that root has a GetAsymmetry method also!
  TH1D* h_numerator = (TH1D*) h_A->Clone(name);
  TH1D* h_denominator = (TH1D*) h_A->Clone();
  h_numerator->SetTitle(title);
  h_numerator->Add(h_B, -1);
  h_denominator->Add(h_B, 1);
  h_numerator->Divide(h_denominator);
  delete h_denominator;
  if (m_intLumi > 0) this->AsymmetryUncertainty(h_numerator, h_A, h_B);
  return h_numerator;
}

void AnalysisZprime::AsymmetryUncertainty(TH1D* h_Asymmetry, TH1D* h_A, TH1D* h_B)
{
  // Find errors
  double efficiency = 1.0;
  double A;
  double sigmaA, sigmaB, sigma;
  double deltaA;
  double N;
  for (int i = 1; i < h_Asymmetry->GetNbinsX(); i++)
  {
    A = h_Asymmetry->GetBinContent(i);
    sigmaA = h_A->GetBinContent(i);
    sigmaB = h_B->GetBinContent(i);
    sigma = (sigmaA + sigmaB)*(h_Asymmetry->GetBinWidth(i));
    N = m_intLumi*efficiency*sigma;
    if (N > 0) deltaA = std::sqrt((1.0 - A*A)/N);
    else deltaA = 0;
    printf("A = %f, dA= %f, N= %f\n", A, deltaA, N);
    h_Asymmetry->SetBinError(i, deltaA);
  }
}

void AnalysisZprime::CreateHistograms() 
{
  h_Mtt = new TH1D("Mff", "M_{tt}", 1300, 0.0, 13.0);
  h_AFBstar1 = new TH1D("AFstar", "AFstar", 130, 0.0, 13.0);
  h_AFBstar2 = new TH1D("ABstar", "ABstar", 130, 0.0, 13.0);
  h_AttC1 = new TH1D("AttC1", "AttC1", 130, 0.0, 13.0);
  h_AttC2 = new TH1D("AttC2", "AttC2", 130, 0.0, 13.0);
  h_ytt = new TH1D("ytt", "ytt", 50, -2.5, 2.5);
  h_CosTheta = new TH1D("CosTheta", "cos#theta", 50, -1.0, 1.0);
  h_CosThetaStar = new TH1D("CosThetaStar", "cos#theta^*", 50, -1.0, 1.0);


  if (m_channel == "tt") {
    h_MttLL = new TH1D("MttLL", "MttLL", 130, 0.0, 13.0);
    h_MttLR = new TH1D("MttLR", "MttLR", 130, 0.0, 13.0);
    h_MttRL = new TH1D("MttRL", "MttRL", 130, 0.0, 13.0);
    h_MttRR = new TH1D("MttRR", "MttRR", 130, 0.0, 13.0);
  }

  if (m_channel == "bbllnn") {
    h_yttReco = new TH1D("yttReco", "yttReco", 100, -2.5, 2.5);
    h_PzNuReco = new TH1D("PzNuReco", "p_{z}^{#nu} (reconstructed)", 100, -1000.0, 1000.0);
    h_MttReco = new TH1D("MttReco", "M^{reco}_{tt}", 100, 0.0, 13.0);
    h_PzNu = new TH1D("PzNu", "p_{z}^{#nu}", 100,-1000.0, 1000.0);
    h_costheta5 = new TH1D("costheta5", "cos#theta_{l^{+}} (eq)", 50, -1.0, 1.0);
    h_CosThetaReco = new TH1D("CosThetaReco", "cos#theta_{t}^{reco}", 50, -1.0, 1.0);
    h_CosThetaStarReco = new TH1D("CosThetaStarReco", "cos#theta_{t}^{*reco}", 50, -1.0, 1.0);
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
    h_AlLF = new TH1D("AlLF", "AlLF", 20, 0.0, 13.0);
    h_AlLB = new TH1D("AlLB", "AlLB", 20, 0.0, 13.0);
    h_ALLF = new TH1D("ALLF", "ALLF", 50, 0.0, 13.0);
    h_ALLB = new TH1D("ALLB", "ALLB", 50, 0.0, 13.0);
    h_AllCF = new TH1D("AllCF", "AllCF", 50, 0.0, 13.0);
    h_AllCB = new TH1D("AllCB", "AllCB", 50, 0.0, 13.0);
    h_AFBstarReco1 = new TH1D("AFBstarNuReco1", "AFBstarNuReco1", 65, 0.0, 13.0);
    h_AFBstarReco2 = new TH1D("AFBstarNuReco2", "AFBstarNuReco2", 65, 0.0, 13.0);
  }
}

void AnalysisZprime::MakeGraphs()
{
  printf("Making Graphs...\n");
  TString numBase;
  if (m_channel == "tt") numBase = "d#sigma(pp->t#bar{t}) / d"; 
  if (m_channel == "bbllnn") numBase = "d#sigma / d"; //pp->t#bar{t}->b#bar{b}l^{+}l^{-}#nu#bar{#nu}
  TString units = "pb";

  // TCanvas *c_Mtt   = new TCanvas(h_Mtt->GetName(), h_Mtt->GetTitle());
  // c_Mtt->cd(); 
  // h_Mtt->Draw("hist"); 
  h_Mtt->GetXaxis()->SetTitle("M_{tt} [TeV]");
  h_Mtt->GetYaxis()->SetTitle(numBase + "M_{tt}" + " [" + units +"/TeV]");

  // TCanvas *c_AFBstar   = new TCanvas(h_AFBstar->GetName(), h_AFBstar->GetTitle());
  // c_AFBstar->cd(); 
  // h_AFBstar->Draw("hist"); 
  // h_AFBstar->GetYaxis()->SetTitle("");

  if (m_channel == "bbllnn") {
    // TCanvas *c_pz5 = new TCanvas(h_pz5->GetName(), h_pz5->GetTitle());
    // c_pz5->cd();
    // h_pz5->Draw("hist"); 
    // h_pz5->GetYaxis()->SetTitle(numBase + h_pz5->GetTitle() + " [" + units +"/GeV]");

    // TCanvas *c_costheta5 = new TCanvas(h_costheta5->GetName(), h_costheta5->GetTitle());
    // c_costheta5->cd(); 
    // h_costheta5->Draw("hist"); 
    // h_costheta5->GetYaxis()->SetTitle(numBase + h_costheta5->GetTitle() + " [" + units +"]");

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

    h_MttReco->GetYaxis()->SetTitle(numBase + h_MttReco->GetTitle() + " [" + units +"/GeV]");
    h_MttReco->GetXaxis()->SetTitle(h_MttReco->GetTitle());

    h_AllC->GetYaxis()->SetTitle(h_AllC->GetTitle());
    h_AllC->GetXaxis()->SetTitle("M_{tt} [GeV]");
  }

  if (m_channel == "tt") {
    h_MttLL->GetYaxis()->SetTitle(numBase + h_MttLL->GetTitle() + " [" + units +"/GeV]");
    h_MttLR->GetYaxis()->SetTitle(numBase + h_MttLR->GetTitle() + " [" + units +"/GeV]");
    h_MttRL->GetYaxis()->SetTitle(numBase + h_MttRL->GetTitle() + " [" + units +"/GeV]");
    h_MttRR->GetYaxis()->SetTitle(numBase + h_MttRR->GetTitle() + " [" + units +"/GeV]");
    h_ALL->GetYaxis()->SetTitle(h_ALL->GetTitle());
    h_AL->GetYaxis()->SetTitle(h_AL->GetTitle());
    h_AFBstar->GetYaxis()->SetTitle(h_AFBstar->GetTitle());
    h_AttC->GetYaxis()->SetTitle(h_AttC->GetTitle());
    h_MttLL->GetXaxis()->SetTitle("M_{tt} [GeV]");
    h_MttLR->GetXaxis()->SetTitle("M_{tt} [GeV]");
    h_MttRL->GetXaxis()->SetTitle("M_{tt} [GeV]");
    h_MttRR->GetXaxis()->SetTitle("M_{tt} [GeV]");
    h_ALL->GetXaxis()->SetTitle("M_{tt} [GeV]");
    h_AL->GetXaxis()->SetTitle("M_{tt} [GeV]");
    h_AFBstar->GetXaxis()->SetTitle("M_{tt} [GeV]");
    h_AttC->GetXaxis()->SetTitle("M_{tt} [GeV]");
  }
  printf("...complete.\n");
}

void AnalysisZprime::WriteHistograms() 
{
  m_outputFile->cd();
  m_outputFile->cd("/");

  printf("Writing histograms...\n");

  // Save histograms

  if (m_channel == "tt" or "bbllnn" or "ll" ) {
    h_Mtt->Write();
  }

  if (m_channel == "tt" || m_channel == "bbllnn") {
    h_AFBstar->Write();
    h_AttC->Write();
    h_ytt->Write();
    h_CosTheta->Write();
    h_CosThetaStar->Write();
  }
  
  if (m_channel == "tt") {
    h_MttLL->Write();
    h_MttLR->Write();
    h_MttRL->Write();
    h_MttRR->Write();
    h_ALL->Write();
    h_AL->Write();
  }

  if (m_channel == "bbllnn") {
    h_dphi->Write();
    h_AL->Write();
    h_ALL->Write();
    h_AllC->Write();
    h_AFBstarReco->Write();
    h_PzNu->Write();
    h_PzNuReco->Write();
    h_MttReco->Write();
    h_costheta5->Write();
    h_ct7ct5->Write();
    h_yttReco->Write();
    h_CosThetaReco->Write();
    h_CosThetaStarReco->Write();
    // h_MET->Write();
    // h_HT->Write();
    // h_Mbbll->Write();
    // h_mll->Write();
    // h_ETbbll->Write();
    // h_KTbbll->Write();
    // h_MTll->Write();
    // h_MCTll->Write();
    // h_MTblbl->Write();
    // h_MCTblbl->Write();
    // h_dphi_HT->Write();
    // h_dphi_Mbbll->Write();
    // h_dphi_mll->Write();
    // h_dphi_ETbbll->Write();
    // h_dphi_KTbbll->Write();
    // h_dphi_MTll->Write();
    // h_dphi_MCTll->Write();
    // h_dphi_MTblbl->Write();
    // h_dphi_MCTblbl->Write();
  }
  h_cutflow->Write();
  printf("...complete\n");
  m_outputFile->Close();
  delete m_outputFile;
}

bool AnalysisZprime::PassCuts() 
{
  // if (this->PassCuts_Mtt())
  // {
  //   if (this->PassCuts_MET())
  //     {
  //       if (this->PassCutsFiducial())
  //       {
        return true;
  //       }
  //     }
  // }
  // return false;
}

bool AnalysisZprime::PassCuts_MET()
{
  bool pass;
  if (m_channel == "ll") pass = true;
  else if (m_channel == "tt") pass = true;
  else if (m_channel == "bbllnn") pass = true;
  else pass = true;

  this->UpdateCutflow(cutMET, pass);
  return pass;
}

bool AnalysisZprime::PassCuts_Mtt()
{
  // if (Mff > 1200.0)
  //   {
  //     if (Mff < 2200.0)
  //       {
          UpdateCutflow(cutMtt, true);
          return true;
  //       }
  //   }
  // UpdateCutflow(cutMtt, false);
  // return false;
}

bool AnalysisZprime::PassCutsFiducial()
{ 
  for (int i = 0; i < 6; i++) {
    bool outsideCrack = etacol[i] <= 1.37 || etacol[i] >= 1.52;
    bool central      = etacol[i] <= 2.47;
    bool passesFiducialCuts = outsideCrack && central;
    if (passesFiducialCuts == false){
      UpdateCutflow(cutFiducial, false);
      return false;
    }
    else continue;
  }
  UpdateCutflow(cutFiducial, true);
  return true;
}

bool AnalysisZprime::PassCutsYtt()
{ 
  if (Mff > 1200.0)
  {
      UpdateCutflow(cutYtt, true);
      return true;
  }
  UpdateCutflow(cutMtt, false);
  return false;
}

void AnalysisZprime::PreLoop()
{
  TString cnormsName(m_weightsFileName);
  ifstream cnorms(cnormsName.Data());
  if (!cnorms.is_open()) printf("Failed to open %s\n", cnormsName.Data());
  cnorms >> m_sigma;
  string line;
  double cnorm;
  while ( cnorms >> cnorm ) m_cnorm.push_back(cnorm);
  cnorms.close();

  this->SetupInputFiles();
  this->InitialiseCutflow();
  this->SetupOutputFiles();
  this->CreateHistograms();
}

void AnalysisZprime::Loop()
{
  // Loop over all files
  for (Itr_s i = m_inputFiles->begin(); i != m_inputFiles->end(); ++i)
  {
    cout << "Processing File '" << (*i) << "'." << endl;
    
    this->SetupTreesForNewFile((*i));
 
    // The Event Loop
    Long64_t nEvents = this->TotalEvents();
    printf("Looping over events...\n");
    for (Long64_t jentry=0; jentry<nEvents;++jentry) 
    {
      // printf("Processing entry %lli\n", jentry);
      Long64_t ientry = this->IncrementEvent(jentry);
      if (ientry < 0) break;
      this->EachEvent();
    }
    printf("...complete.\n");
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
  printf("Histograms will be written to '%s'.\n", m_outputFileName.Data());
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
  return -999;  
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

double AnalysisZprime::ResolveNeutrinoPz(TLorentzVector p_l, TVector2 pT_nu) 
{

  // finds the longitudinal neutrino momentum for semi-hadronic decay
  // assuming all particles are massless

  double PzNu;
  std::vector<std::complex<double> > root, root2;
  double a = -999, b = -999, c = -999, k = -999;
  // double a2 = -999, b2 = -999, c2 = -999;

  // recalculate lepton energy in zero mass approximation
  double p_l0 = std::sqrt(p_l.Px()*p_l.Px() + p_l.Py()*p_l.Py() + p_l.Pz()*p_l.Pz());

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

  root = this->SolveQuadratic(a, b, c);
  // root2 = this->SolveQuadratic(a2, b2, c2);

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

std::vector<std::complex<double> > AnalysisZprime::SolveQuadratic(double a, double b, double c) 
{
    // solves quadratic for both roots 
    // returns both as complex values in a complex vector x(2)

    std::vector<std::complex<double> > roots;
    std::complex<double> term1;
    std::complex<double> term2;
    std::complex<double> discriminator;

    term1 = -b/(2*a);
    discriminator = b*b - 4*a*c;
    term2 = std::sqrt(discriminator)/(2*a);

    roots.push_back(term1 + term2);
    roots.push_back(term1 - term2);
    
    return roots;
}

const void AnalysisZprime::UpdateCutflow(int cut, bool passed)
{
  if (cutflow[cut] == -999) cutflow[cut] = 0;
  if (passed) cutflow[cut] +=1;
}

void AnalysisZprime::InitialiseCutflow() 
{
  cutflow = std::vector<int>(nCuts, -999);
  cutNames = std::vector<TString>(nCuts, "no name");
  cutNames[cutEvent] = "Event";
  cutNames[cutMET] = "MET";
  cutNames[cutMtt] = "Mff";
  cutNames[cutFiducial] = "Fiducial";


  h_cutflow = new TH1D("Cutflow", "Cutflow", nCuts, 0, nCuts);
}

void AnalysisZprime::PrintCutflow() 
{
  printf("------\n");
  for (int cut = 0; cut < nCuts; cut++) {
    if (cutflow[cut] == -999) continue;

    h_cutflow->SetBinContent(cut+1, cutflow[cut]);
    h_cutflow->GetXaxis()->SetBinLabel(cut+1, cutNames[cut]);

    printf("%s cut: %i pass\n", cutNames[cut].Data(), cutflow[cut]);
  }
  printf("------\n");
}
