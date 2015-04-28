#include "analysis.h"

AnalysisZprime::AnalysisZprime(const TString channel, const TString& outputFileName) :
  m_pi(3.14159265),
  m_GeV(1000.0),
  m_channel(channel),
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
  if (this->PassCuts())
  {    
    double weight = m_ntup->weight();
    double weight_ee = m_ntup->weight_ee();
    double weight_eq = m_ntup->weight_eq();

    // Fill Histograms (assumes fixed bin width!)

    h_Mtt->Fill(m_ntup->Mtt(), weight/h_Mtt->GetXaxis()->GetBinWidth(1));

    if (m_channel == "2to6") {
      if (m_ntup->barcode()->at(2) == -11)  h_pz5->Fill(m_ntup->Pz()->at(2), m_ntup->weight());
      // if (m_ntup->Pz()->size() > 0) {
      //   for(unsigned int i = 0; i < m_ntup->Pz()->size(); ++i) {
      //     h_pz5->Fill(m_ntup->Pz()->at(i) / m_GeV, eventWeight);
      //   }
      // }
      h_costheta5_eq->Fill(m_ntup->costheta5(), weight_eq);
      h_costheta5_ee->Fill(m_ntup->costheta5(), weight_ee);
      h_ct7ct5->Fill(m_ntup->ct7ct5(), weight_ee);
    }

    if (m_channel == "2to2") {
      h_MttLL->Fill(m_ntup->Mtt(), m_ntup->weightLL()/h_MttLL->GetXaxis()->GetBinWidth(1));
      h_MttLR->Fill(m_ntup->Mtt(), m_ntup->weightLR()/h_MttLR->GetXaxis()->GetBinWidth(1));
      h_MttRL->Fill(m_ntup->Mtt(), m_ntup->weightRL()/h_MttRL->GetXaxis()->GetBinWidth(1));
      h_MttRR->Fill(m_ntup->Mtt(), m_ntup->weightRR()/h_MttRR->GetXaxis()->GetBinWidth(1));
      // printf("%e %e\n", m_ntup->weightLL() + m_ntup->weightLR() + m_ntup->weightRL() + m_ntup->weightRR(), eventWeight);
    }    
  }
}

void AnalysisZprime::PreLoop()
{
  this->SetupInputFiles();
  
  // Define output file 
  m_outputFile = new TFile(m_outputFileName,"RECREATE");
  
  // Define Histograms    
  h_Mtt = new TH1D("Mtt", "M_{tt}", 100, 0.0, 14000.0);

  if (m_channel == "2to6") {  
    h_pz5 = new TH1D("pz5", "p_{z}^{5}", 50,0.0, 10000.0);
    // h_pz5->Sumw2();

    h_costheta5_eq = new TH1D("costheta5_eq", "cos#theta_{l^{+}} (eq)", 50, -1.0, 1.0);
    // h_costheta5_eq->Sumw2();

    h_costheta5_ee = new TH1D("costheta5_ee", "cos#theta_{l^{+}} (ee)", 50, -1.0, 1.0);
    // h_costheta5_ee->Sumw2();

    h_ct7ct5 = new TH1D("ct7ct5", "cos#theta_{l^{+}}cos#theta_{l^{-}}", 50, -1.0, 1.0);
    // h_ct7ct5->Sumw2();
  }

  if (m_channel == "2to2") {
    h_MttLL = new TH1D("MttLL", "MttLL", 100, 0.0, 14000.0);
    h_MttLR = new TH1D("MttLR", "MttLR", 100, 0.0, 14000.0);
    h_MttRL = new TH1D("MttRL", "MttRL", 100, 0.0, 14000.0);
    h_MttRR = new TH1D("MttRR", "MttRR", 100, 0.0, 14000.0);
    h_MttALLnum = new TH1D("MttALLnum", "MttALLnum", 100, 0.0, 14000.0);
    h_MttALLden = new TH1D("MttALLden", "MttALLden", 100, 0.0, 14000.0);
    h_MttALL = new TH1D("MttALL", "MttALL", 100, 0.0, 14000.0);
    h_Mtt->Sumw2();
    h_MttLL->Sumw2();
    h_MttLR->Sumw2();
    h_MttRL->Sumw2();
    h_MttRR->Sumw2();
    h_MttALLnum->Sumw2();
    h_MttALLden->Sumw2();
    h_MttALL->Sumw2();
  }
}

void AnalysisZprime::PostLoop()
{
  if (m_channel == "2to2") this->TotalAsymmetries();

  this->MakeGraphs();

  this->ALL2to6();
  
  m_outputFile->cd();
  m_outputFile->cd("/");

  // Save histograms
  h_Mtt->Write();

  if (m_channel == "2to2") {
    h_MttLL->Write();
    h_MttLR->Write();
    h_MttRL->Write();
    h_MttRR->Write();
    h_MttALL->Write();
  }

  if (m_channel == "2to6") {
    h_pz5->Write();
    h_costheta5_eq->Write();
    h_costheta5_ee->Write();
    h_ct7ct5->Write();
  }

  m_outputFile->Close();
  
}

void AnalysisZprime::TotalAsymmetries()
{
  h_MttALL->Add(h_MttRR, 1);
  h_MttALL->Add(h_MttLL, 1);
  h_MttALL->Add(h_MttRL,-1);
  h_MttALL->Add(h_MttLR,-1);

  h_MttALLden->Add(h_MttRR, 1);
  h_MttALLden->Add(h_MttLL, 1);
  h_MttALLden->Add(h_MttRL, 1);
  h_MttALLden->Add(h_MttLR, 1);

  h_MttALL->Divide(h_MttALLden);

  double sigma = h_Mtt->Integral("width");
  double sigmaLL = h_MttLL->Integral("width");
  double sigmaLR = h_MttLR->Integral("width");
  double sigmaRL = h_MttRL->Integral("width");
  double sigmaRR = h_MttRR->Integral("width");

  double ALL = (sigmaLL + sigmaRR - sigmaRL - sigmaLR)/
               (sigmaLL + sigmaRR + sigmaRL + sigmaLR);

  printf("sigma = %f\n", sigma);
  printf("ALL = %f\n", ALL);  
}

void AnalysisZprime::ALL2to6() {
  double mean = h_ct7ct5->GetMean();

  double ALL = -9*mean;

  printf("ALL = %f\n", ALL);
}

void AnalysisZprime::MakeGraphs()
{
  printf("Making Graphs...\n");

  TString numBase = "d#sigma / d"; 
  TString units = "";

  TCanvas *c_Mtt   = new TCanvas(h_Mtt->GetName(), h_Mtt->GetTitle());
  c_Mtt->cd(); 
  h_Mtt->Draw("hist"); 
  h_Mtt->GetYaxis()->SetTitle("");

  if (m_channel == "2to6") {
    TCanvas *c_pz5 = new TCanvas(h_pz5->GetName(), h_pz5->GetTitle());
    c_pz5->cd();
    h_pz5->Draw("hist"); 
    h_pz5->GetYaxis()->SetTitle(numBase + h_pz5->GetTitle() + " [" + units +"/GeV]");

    TCanvas *c_costheta5_eq = new TCanvas(h_costheta5_eq->GetName(), h_costheta5_eq->GetTitle());
    c_costheta5_eq->cd(); 
    h_costheta5_eq->Draw("hist"); 
    h_costheta5_eq->GetYaxis()->SetTitle(numBase + h_costheta5_eq->GetTitle() + " [" + units +"]");

    TCanvas *c_costheta5_ee = new TCanvas(h_costheta5_ee->GetName(), h_costheta5_ee->GetTitle());
    c_costheta5_ee->cd(); 
    h_costheta5_ee->Draw("hist"); 
    h_costheta5_ee->GetYaxis()->SetTitle(numBase + h_costheta5_ee->GetTitle() + " [" + units +"]");

    TCanvas *c_ct7ct5 = new TCanvas(h_ct7ct5->GetName(), h_ct7ct5->GetTitle());
    c_ct7ct5->cd(); 
    h_ct7ct5->Draw("hist"); 
    h_ct7ct5->GetYaxis()->SetTitle(numBase + h_ct7ct5->GetTitle() + " [" + units +"]");
  }

  if (m_channel == "2to2") {
    TCanvas *c_MttLL   = new TCanvas("MttLL " ,"MttLL " );
    c_MttLL->cd(); 
    h_MttLL->Draw("hist"); 
    h_MttLL->GetYaxis()->SetTitle("Events");

    TCanvas *c_MttLR   = new TCanvas("MttLR " ,"MttLR " );
    c_MttLR->cd(); 
    h_MttLR->Draw("hist"); 
    h_MttLR->GetYaxis()->SetTitle("Events");

    TCanvas *c_MttRL   = new TCanvas("MttRL " ,"MttRL " );
    c_MttRL->cd(); 
    h_MttRL->Draw("hist"); 
    h_MttRL->GetYaxis()->SetTitle("Events");

    TCanvas *c_MttRR   = new TCanvas("MttRR " ,"MttRR " );
    c_MttRR->cd(); 
    h_MttRR->Draw("hist"); 
    h_MttRR->GetYaxis()->SetTitle("Events");

    TCanvas *c_MttALL   = new TCanvas("MttALL " ,"MttALL " );
    c_MttALL->cd(); 
    h_MttALL->Draw("hist"); 
    h_MttALL->GetYaxis()->SetTitle("Events");
  }
  printf("...complete.\n");
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
    if (m_ntup->Etmiss() > 0.0 * m_GeV)
    {
      return true;
    }
  }  
  return false;
}

bool AnalysisZprime::PassCuts_Mtt() const
{
  if (m_ntup->Mtt() > 1200.0)
    {
      if (m_ntup->Mtt() < 2200.0)
        {
          return true;
        }
    }
  return false;
}

void AnalysisZprime::Loop()
{
  // Loop over all files
  for(Itr_s i = m_inputFiles->begin(); i != m_inputFiles->end(); ++i)
  {
    cout<<"  Processing File = "<<(*i)<<endl;
    
    this->SetupTreesForNewFile((*i));
 
    // The Event Loop
    Long64_t nEvents = this->TotalEvents();
    for(Long64_t jentry=0; jentry<nEvents;++jentry) 
    {
      printf("Processing entry %lli\n", jentry);
      Long64_t ientry = this->IncrementEvent(jentry);
      if (ientry < 0) break;
      this->EachEvent();
    }
    
    this->CleanUp();
  }
}  

AnalysisZprime::~AnalysisZprime()
{       
  delete m_inputFiles;
}

void AnalysisZprime::SetupInputFiles()
{   
  m_inputFiles = new vector<TString>;
  TString base("/afs/cern.ch/work/d/demillar/Ntuples_Zprime/");
  // TString base("/Users/declan/Data/Ntuples_Zprime/");
  
  m_inputFiles->push_back(base + "SM_13_2to6_xc_1x1000000.root");
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