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

void AnalysisZprime::EachEvent(TString l_Q)
{
  UpdateCutflow(cutEvent, true);

  p = vector<TLorentzVector>(6);
  for (int i = 0; i < (int) m_ntup->E()->size(); i++) {
    p[i].SetPxPyPzE(m_ntup->Px()->at(i), m_ntup->Py()->at(i), m_ntup->Pz()->at(i), m_ntup->E()->at(i));
  }

  for (int i = 0; i < p.size(); i++) {
    Pcm += pcm[i];
    pcm[i] = p[i];
    p_r1[i] = p[i];
    p_r2[i] = p[i];
  }
   
  TVector3 V = -1*P.BoostVector();

  pcm = vector<TLorentzVector>(6);
  for (int i = 0; i < p.size(); i++) {
    pcm[i].Boost(V);
    Pcm += pcm[i];
  }  

  if (m_channel == "bbllnn") {
    // resolve longitudinal neutrino momentum in the semihadronic case
    this->Resolvebbnu(p,1);
    this->Resolvebbnu(p,1);

    
    // reconstructed final particle parton CoM variables
    TVector3 V_r1 = -1*P_r1.BoostVector();
    TVector3 V_r2 = -1*P_r2.BoostVector();
    for (int i = 0; i < p.size(); i++) {
      pcm_r1[i].Boost(V_r1);
      pcm_r2[i].Boost(V_r2);
    }  
  }

  // top and antitop
  if (m_channel == "tt" or m_channel == "ll") {
    p_t = pcm[0];
    p_tb = pcm[1];
  }
  else if (m_channel == "bbllnn") {
    p_t = pcm[0] + pcm[2] + pcm[3];
    p_tb = pcm[1] + pcm[4] + pcm[5];
    p_t_r1 = pcm_r1[0] + pcm_r1[2] + pcm_r1[3];
    p_tb_r1 = pcm_r1[1] + pcm_r1[4] + pcm_r1[5];
    p_t_r2 = pcm_r2[0] + pcm_r2[2] + pcm_r2[3];
    p_tb_r2 = pcm_r2[1] + pcm_r2[4] + pcm_r2[5];
  }

  double y_t = p_t.Rapidity();
  double y_tb = p_tb.Rapidity();
  double dy = std::abs(y_t) - std::abs(y_tb);
  double ytt = P.Rapidity();
  double CosTheta = p_t.CosTheta();
  double CosThetaStar = int(ytt/std::abs(ytt))*CosTheta;
  double ytt_r1 = -999;
  double ytt_r2 = -999;
  double CosTheta_r1 = -999;
  double CosTheta_r2 = -999;
  double CosThetaStar_r1 = -999;
  double CosThetaStar_r2 = -999;
  if (m_channel == "bbllnn"){
    ytt_r1 = P_r1.Rapidity();
    ytt_r2 = P_r2.Rapidity();
    CosTheta_r1 = p_t_r1.CosTheta();
    CosTheta_r2 = p_t_r2.CosTheta();
    CosThetaStar_r1 = int(ytt_r1/std::abs(ytt_r1))*CosTheta_r1;
    CosThetaStar_r2 = int(ytt_r2/std::abs(ytt_r2))*CosTheta_r2;
  }
    
  if (this->PassCuts())
  {    
    // re-weight for different iterations
    double it = m_ntup->iteration();
    double weight = m_ntup->weight();
    weight = weight*m_sigma/m_weights[it-1];

    // convert to TeV
    Mff = P.M()/1000;
    Mtt_r1 = P_r1.M()/1000;
    Mtt_r2 = P_r2,M()/1000;

    // fill histograms (assumes fixed bin width!)
    h_Mff->Fill(Mff, weight/h_Mtt->GetXaxis()->GetBinWidth(1));
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
      h_Pz_nu->Fill(pcol[3].Pz(), weight/h_PzNu->GetXaxis()->GetBinWidth(1));

      h_ytt_r->Fill(ytt_r1, weight/2/h_ytt_r->GetXaxis()->GetBinWidth(1));
      h_ytt_r->Fill(ytt_r2, weight/2/h_ytt_r->GetXaxis()->GetBinWidth(1));

      h_CosTheta_r->Fill(CosTheta_r, weight/2/h_CosTheta_r->GetXaxis()->GetBinWidth(1));
      h_CosTheta_r->Fill(CosTheta_r, weight/2/h_CosTheta_r->GetXaxis()->GetBinWidth(1));

      h_CosThetaStar_r->Fill(CosThetaStar_r1, weight/2/h_CosThetaStar_r->GetXaxis()->GetBinWidth(1));
      h_CosThetaStar_r->Fill(CosThetaStar_r2, weight/2/h_CosThetaStar_r->GetXaxis()->GetBinWidth(1));

      h_Pz_nu_r->Fill(p_r1[3] weight/2/h_PzNu_r->GetXaxis()->GetBinWidth(1));
      h_Pz_nu_r->Fill(p_r2[5], weight/2/h_PzNu_r->GetXaxis()->GetBinWidth(1));

      h_Mtt_r->Fill(Mtt_r1, weight/2/h_Mtt_r->GetXaxis()->GetBinWidth(1));
      h_Mtt_r->Fill(Mtt_r2, weight/2/h_Mtt_r->GetXaxis()->GetBinWidth(1));

      if (CosThetaStar_r1 > 0) {
        h_AFBstarReco1->Fill(Mtt_r1, weight/2/h_AFBstarReco1->GetXaxis()->GetBinWidth(1));
      }
      if (CosThetaStar_r2 > 0) {
        h_AFBstarReco1->Fill(Mtt_r2, weight/2/h_AFBstarReco1->GetXaxis()->GetBinWidth(1));
      }

      if (CosThetaStar_r1 < 0) {
        h_AFBstarReco2->Fill(Mtt_r1, weight/2/h_AFBstarReco2->GetXaxis()->GetBinWidth(1));
      }
      if (CosThetaStar_r2 < 0) {
        h_AFBstarReco2->Fill(Mtt_r2, weight/2/h_AFBstarReco2->GetXaxis()->GetBinWidth(1));
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

  printf("b-assigniment success rate: %f\n", m_bAssignSuccesses/m_bAssignAttempts);

  if (m_channel == "bbllnn") {
    this->ALL2to6();
    h_AlL = this->Asymmetry("AL", "A_{L}", h_AlLF, h_AlLB);
    h_AllC = this->Asymmetry("AllC", "A^{ll}_{C}", h_AllCF, h_AllCB);
    h_AFBstarReco = this->Asymmetry("AFBstarNuReco", "A_{FB}^* (#nu reconstructed)", h_AFBstarReco1, h_AFBstarReco2);
  }
  this->MakeGraphs();
  this->PrintCutflow();
  this->WriteHistograms();
}

TH1D* AnalysisZprime::PlotALL()
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

TH1D* AnalysisZprime::PlotAL()
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
  h_Mff = new TH1D("Mff", "M_{tt}", 1300, 0.0, 13.0);
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
    h_Pz_nu = new TH1D("PzNu", "p_{z}^{#nu}", 100,-1000.0, 1000.0);

    h_ytt_r = new TH1D("yttReco", "y_{tt}^{reco}", 100, -2.5, 2.5);
    h_Pz_nu_r = new TH1D("PzNuReco", "p_{z}^{#nu} (reconstructed)", 100, -1000.0, 1000.0);
    h_Mtt_r = new TH1D("MttReco", "M^{reco}_{tt}", 100, 0.0, 13.0);

    h_AlLF = new TH1D("AlLF", "AlLF", 20, 0.0, 13.0);
    h_AlLB = new TH1D("AlLB", "AlLB", 20, 0.0, 13.0);
    h_AllCF = new TH1D("AllCF", "AllCF", 50, 0.0, 13.0);
    h_AllCB = new TH1D("AllCB", "AllCB", 50, 0.0, 13.0);
    h_AFBstar_rF = new TH1D("AFBstarNuReco1", "AFBstarNuReco1", 65, 0.0, 13.0);
    h_AFBstar_rB = new TH1D("AFBstarNuReco2", "AFBstarNuReco2", 65, 0.0, 13.0);
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

bool AnalysisZprime::PassCutsMET()
{
  bool pass;
  if (m_channel == "ll") pass = true;
  else if (m_channel == "tt") pass = true;
  else if (m_channel == "bbllnn") pass = true;
  else pass = true;

  this->UpdateCutflow(c_MET, pass);
  return pass;
}

bool AnalysisZprime::PassCutsMtt()
{
  // if (Mff > 1200.0)
  //   {
  //     if (Mff < 2200.0)
  //       {
          UpdateCutflow(c_Mtt, true);
          return true;
  //       }
  //   }
  // UpdateCutflow(c_Mtt, false);
  // return false;
}

bool AnalysisZprime::PassCutsFiducial()
{ 
  for (int i = 0; i < p.size[i]; i++) {
    bool outsideCrack = p[i].PseudoRapidity() <= 1.37 || p[i].PseudoRapidity() >= 1.52;
    bool central      = p[i].PseudoRapidity() <= 2.47;
    bool passesFiducialCuts = outsideCrack && central;
    if (passesFiducialCuts == false){
      UpdateCutflow(c_Fiducial, false);
      return false;
    }
    else continue;
  }
  UpdateCutflow(c_Fiducial, true);
  return true;
}

bool AnalysisZprime::PassCutsYtt()
{ 
  if (std::abs(P.Rapidity()) > 0)
  {
    UpdateCutflow(c_Ytt, true);
    return true;
  }
  UpdateCutflow(c_Mtt, false);
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
  while ( cnorms >> cnorm ) m_weights.push_back(cnorm);
  cnorms.close();
  m_bAssignAttempts = 0;
  m_bAssignSuccesses = 0;

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

// std::vector<TLorentzVector> AnalysisZprime::ResolveNeutrinoPz(std::vector<TLorentzVector>,,3) 
// {

//   // finds the longitudinal neutrino momentum for semi-hadronic decay
//   // assuming all particles are massless

//   double PzNu;
//   std::vector<std::complex<double> > root, root2;
//   double a = -999, b = -999, c = -999, k = -999;
//   // double a2 = -999, b2 = -999, c2 = -999;

//   // recalculate lepton energy in zero mass approximation
//   double p_l0 = std::sqrt(p_l.Px()*p_l.Px() + p_l.Py()*p_l.Py() + p_l.Pz()*p_l.Pz());

//   if ( std::abs(p_l0 - p_l.E() > 0.00001)) printf("p_l0 doesn't match\n");


//   k = m_Wmass*m_Wmass/2 + p_l.Px()*pT_nu.Px() + p_l.Py()*pT_nu.Py();

//   a = p_l.Px()*p_l.Px() + p_l.Py()*p_l.Py();

//   b = -2*k*(p_l.Pz());

//   c = (pT_nu.Px()*pT_nu.Px() + pT_nu.Py()*pT_nu.Py())*p_l.E()*p_l.E() - k*k;

//   root = this->SolveQuadratic(a, b, c);

//   // select single solution
//   if (root[0].imag() == 0 and root[1].imag() == 0) {
//     // two real solutions - pick smallest one
//     if (std::abs(root[0].real()) < std::abs(root[1].real())) {
//       // solution 1 < than solution 2
//       PzNu = root[0].real();
//     }
//     else if (std::abs(root[0].real()) > std::abs(root[1].real())) { 
//       // solution 1 > than solution 2
//       PzNu = root[1].real();
//     }
//     else {
//       // solutions are equal pick 1
//       PzNu = root[0].real();
//     }
//   }
//   else {
//     // no real solutions - take the real part of 1
//     PzNu = root[0].real();
//   }
// }

std::vector<TLorentzVector> KinematicsZprime::resolvebbnu(std::vector<TLorentzVector> p, int l_Q) 
{

  // finds the longitudinal neutrino momentum and matches b-quarks to the leptonially or hadronically decayin top quark
  //  for semi-leptonic decay

  TLorentzVector p_l;
  TLorentzVector p_nu;
  std::vector<TLorentzVector> p_b(2);
  std::vector<TLorentzVector> p_q(2);

  p_b[0] = p[0];
  p_b[1] = p[1]
  if (l_Q = 1) {
    p_l = p[2];
    p_nu = p[3];
    p_q[0] = p[4];
    p_q[1]= p[5];
  }
  else if )l_Q = 1) {
    p_l = p[4];
    p_nu = p[5];
    p_q[0] = p[2];
    p_q[1]= p[3];
  }
  else {
    printf("Invalid charge.\n");
    return 0;
  }

  double PzNu;
  double X2, X2min;

  double px_nu = p_nu.Px();
  double py_nu = p_nu.Py();
  double pz_nu = -999
  std::vector<std::complex<double>> root;
  double a = -999, b = -999, c = -999, k = -999;

  // recalculate lepton energy in zero mass approximation
  double p_l0 = std::sqrt(p_l.Px()*p_l.Px() + p_l.Py()*p_l.Py() + p_l.Pz()*p_l.Pz());

  if ( std::abs(p_l0 - p_l.E() > 0.00001)) printf("p_l0 doesn't match\n");

  k = m_Wmass*m_Wmass/2 + p_l.Px()*pT_nu.Px() + p_l.Py()*pT_nu.Py();

  a = p_l.Px()*p_l.Px() + p_l.Py()*p_l.Py();

  b = -2*k*(p_l.Pz());

  c = (pT_nu.Px()*pT_nu.Px() + pT_nu.Py()*pT_nu.Py())*p_l.E()*p_l.E() - k*k;

  root = this->solveQuadratic(a, b, c);

  // select single solution

  if (root[0].imag() == 0 and root[1].imag() == 0) {
    for (int j = 0; j < 2; i++) {
      for (int j = 0; j < 2; j++) {
        E_nu = sqrt(px_nu*px_nu + py_nu*py_nu + root[i]*root[i])
        p_nu.SetPxPyPzE(px_nu, py_nu, root[i], E_nu);
        mblv = (p_b[std::abs(j-1)] + p_l + p_nu).M();
        mjjb = (p_b[std::abs(j)] + p_q[0] + p_q[1]).M();
        dh = mjjb - tmass;
        dl = mblv - tmass;
        X2 = dh*dh + dl*dl;
        if (it == 0){
          X2min = X2;
          imin = i;
          jmin = j;
        if (X2 < X2min){
          X2min = X2;
          imin = i;
          jmin = j;
        } 
      }
    }
    
  else {
    // no real solutions - take the real part of 1
    printf("NO REAL SOLUTIONS!\n");;
  }
  pz_nu = root[imin];
  pz_nu_truth = p_nu.Pz();
  bool RecoMatchesTruth(jmin = 0);
  if (RecoMatchesTruth) printf("b-assignment: correct. \n");
  else printf("b-assignment: incorrect. \n");

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
  cutflow = std::vector<int>(m_cuts, -999);
  cutNames = std::vector<TString>(m_cuts, "no name");
  cutNames[c_Event] = "Event";
  cutNames[c_MET] = "MET";
  cutNames[c_Mtt] = "Mff";
  cutNames[c_Fiducial] = "Fiducial";


  h_cutflow = new TH1D("Cutflow", "Cutflow", m_cuts, 0, m_cuts);
}

void AnalysisZprime::PrintCutflow() 
{
  printf("------\n");
  for (int cut = 0; cut < m_cuts; cut++) {
    if (cutflow[cut] == -999) continue;

    h_cutflow->SetBinContent(cut+1, cutflow[cut]);
    h_cutflow->GetXaxis()->SetBinLabel(cut+1, cutNames[cut]);

    printf("%s cut: %i pass\n", cutNames[cut].Data(), cutflow[cut]);
  }
  printf("------\n");
}
