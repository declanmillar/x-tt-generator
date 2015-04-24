#include "analysis_transverse.h"

AnalysisTransverse::AnalysisTransverse(const TString& inputFileName, const TString& outputFileName) :
  m_pi(3.14159265),
  m_GeV(1000.0),
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

void AnalysisTransverse::EachEvent()
{
  if( this->PassCuts() )
  {
    // cout << "Event = " << m_ntup->eventNumber()<<endl;
    
    double eventWeight = m_ntup->weight();

    h_mtt->Fill(m_ntup->mtt(), eventWeight/h_mtt->GetXaxis()->GetBinWidth(1));
    h_ht->Fill(m_ntup->ht(), eventWeight/h_ht->GetXaxis()->GetBinWidth(1));
    h_mttvis->Fill(m_ntup->mttvis(), eventWeight/h_mttvis->GetXaxis()->GetBinWidth(1));
    h_mt1->Fill(m_ntup->mt1(), eventWeight/h_mt1->GetXaxis()->GetBinWidth(1));
    h_mt2->Fill(m_ntup->mt2(), eventWeight/h_mt2->GetXaxis()->GetBinWidth(1));
    h_mt3->Fill(m_ntup->mt3(), eventWeight/h_mt3->GetXaxis()->GetBinWidth(1));
    h_mct1->Fill(m_ntup->mct1(), eventWeight/h_mct1->GetXaxis()->GetBinWidth(1));
    h_mct2->Fill(m_ntup->mct2(), eventWeight/h_mct2->GetXaxis()->GetBinWidth(1));
    h_mct3->Fill(m_ntup->mct3(), eventWeight/h_mct3->GetXaxis()->GetBinWidth(1));
    h_mlt->Fill(m_ntup->mlt(), eventWeight/h_mlt->GetXaxis()->GetBinWidth(1));
    h_mlct->Fill(m_ntup->mlct(), eventWeight/h_mlct->GetXaxis()->GetBinWidth(1));

    h_dphi_mtt->Fill(m_ntup->dphi(), m_ntup->mtt(), eventWeight/h_dphi_mtt->GetXaxis()->GetBinWidth(1)/h_dphi_mtt->GetYaxis()->GetBinWidth(1));
    h_dphi_ht->Fill(m_ntup->dphi(), m_ntup->ht(), eventWeight/h_dphi_ht->GetXaxis()->GetBinWidth(1)/h_dphi_ht->GetYaxis()->GetBinWidth(1));
    h_dphi_mttvis->Fill(m_ntup->dphi(), m_ntup->mttvis(), eventWeight/h_dphi_mttvis->GetXaxis()->GetBinWidth(1)/h_dphi_mttvis->GetYaxis()->GetBinWidth(1));
    h_dphi_mt1->Fill(m_ntup->dphi(), m_ntup->mt1(), eventWeight/h_dphi_mt1->GetXaxis()->GetBinWidth(1)/h_dphi_mt1->GetYaxis()->GetBinWidth(1));
    h_dphi_mt2->Fill(m_ntup->dphi(), m_ntup->mt2(), eventWeight/h_dphi_mt2->GetXaxis()->GetBinWidth(1)/h_dphi_mt2->GetYaxis()->GetBinWidth(1));
    h_dphi_mt3->Fill(m_ntup->dphi(), m_ntup->mt3(), eventWeight/h_dphi_mt3->GetXaxis()->GetBinWidth(1)/h_dphi_mt3->GetYaxis()->GetBinWidth(1));
    h_dphi_mct1->Fill(m_ntup->dphi(), m_ntup->mct1(), eventWeight/h_dphi_mct1->GetXaxis()->GetBinWidth(1)/h_dphi_mct1->GetYaxis()->GetBinWidth(1));
    h_dphi_mct2->Fill(m_ntup->dphi(), m_ntup->mct2(), eventWeight/h_dphi_mct2->GetXaxis()->GetBinWidth(1)/h_dphi_mct2->GetYaxis()->GetBinWidth(1));
    h_dphi_mct3->Fill(m_ntup->dphi(), m_ntup->mct3(), eventWeight/h_dphi_mct3->GetXaxis()->GetBinWidth(1)/h_dphi_mct3->GetYaxis()->GetBinWidth(1));
    h_dphi_mlt->Fill(m_ntup->dphi(), m_ntup->mlt(), eventWeight/h_dphi_mlt->GetXaxis()->GetBinWidth(1)/h_dphi_mlt->GetYaxis()->GetBinWidth(1));
    h_dphi_mlct->Fill(m_ntup->dphi(), m_ntup->mlct(), eventWeight/h_dphi_mlct->GetXaxis()->GetBinWidth(1)/h_dphi_mlct->GetYaxis()->GetBinWidth(1));
  }
}

void AnalysisTransverse::PreLoop()
{
  this->SetupInputFiles();
  
  // Define output file 
  m_outputFile = new TFile(m_outputFileName,"RECREATE");
  
  // Define Histograms  
  h_mtt = new TH1D("mtt", "mtt", 100, 0.0, 14000.0);
  h_ht = new TH1D("ht", "ht", 100, 0.0, 14000.0);
  h_mttvis = new TH1D("mvis", "mvis", 100, 0.0, 14000.0);
  h_mt1 = new TH1D("mt1", "mt", 100, 0.0, 14000.0);
  h_mt2 = new TH1D("mt2", "mt", 100, 0.0, 14000.0);
  h_mt3 = new TH1D("mt3", "mt", 100, 0.0, 14000.0);
  h_mct1 = new TH1D("mct1", "mct", 100, 0.0, 14000.0);
  h_mct2 = new TH1D("mct2", "mct", 100, 0.0, 14000.0);
  h_mct3 = new TH1D("mct3", "mct", 100, 0.0, 14000.0);
  h_mlt = new TH1D("mlt", "ml", 100, 0.0, 14000.0);
  h_mlct = new TH1D("mlct", "mlc", 100, 0.0, 14000.0);

  h_dphi_mtt = new TH2D("dphi_mtt", "dphi_mtt", 10, 0, 2*m_pi, 100, 0.0, 14000.0);
  h_dphi_ht = new TH2D("dphi_ht", "dphi_ht", 10, 0, 2*m_pi, 100, 0.0, 14000.0);
  h_dphi_mttvis = new TH2D("dphi_mttvis", "dphi_mttvis", 10, 0, 2*m_pi, 100, 0.0, 14000.0);
  h_dphi_mt1 = new TH2D("dphi_mt1", "dphi_mt1", 10, 0, 2*m_pi, 100, 0.0, 14000.0);
  h_dphi_mt2 = new TH2D("dphi_mt2", "dphi_mt2", 10, 0, 2*m_pi, 100, 0.0, 14000.0);
  h_dphi_mt3 = new TH2D("dphi_mt3", "dphi_mt3", 10, 0, 2*m_pi, 100, 0.0, 14000.0);
  h_dphi_mct1 = new TH2D("dphi_mct1", "dphi_mct1", 10, 0, 2*m_pi, 100, 0.0, 14000.0);
  h_dphi_mct2 = new TH2D("dphi_mct2", "dphi_mct2", 10, 0, 2*m_pi, 100, 0.0, 14000.0);
  h_dphi_mct3 = new TH2D("dphi_mct3", "dphi_mct3", 10, 0, 2*m_pi, 100, 0.0, 14000.0);
  h_dphi_mlt = new TH2D("dphi_mlt", "dphi_mlt", 10, 0, 2*m_pi, 100, 0.0, 14000.0);
  h_dphi_mlct = new TH2D("dphi_mlct", "dphi_mlct", 10, 0, 2*m_pi, 100, 0.0, 14000.0);
}

void AnalysisTransverse::PostLoop()
{
  this->MakeGraphs();
  
  m_outputFile->cd();
  m_outputFile->cd("/");

  // Save histograms
  h_mtt->Write();
  h_ht->Write();
  h_mttvis->Write();
  h_mt1->Write();
  h_mt2->Write();
  h_mt3->Write();
  h_mct1->Write();
  h_mct2->Write();
  h_mct3->Write();
  h_mlt->Write();
  h_mlct->Write();
  h_dphi_mtt->Write();
  h_dphi_ht->Write();
  h_dphi_mttvis->Write();
  h_dphi_mt1->Write();
  h_dphi_mt2->Write();
  h_dphi_mt3->Write();
  h_dphi_mct1->Write();
  h_dphi_mct2->Write();
  h_dphi_mct3->Write();
  h_dphi_mlt->Write();
  h_dphi_mlct->Write();

  m_outputFile->Close();
  
}

void AnalysisTransverse::MakeGraphs()
{
  printf("Making Graphs...\n");

  // TCanvas *c_mtt   = new TCanvas( "mtt " ,"mtt "  );
  // c_mtt->cd(); 
  // h_mtt->Draw(); 
  // h_mtt->GetYaxis()->SetTitle( "Events" );

  // TCanvas *c_ht   = new TCanvas( "ht " ,"ht "  );
  // c_ht->cd(); 
  // h_ht->Draw(); 
  // h_ht->GetYaxis()->SetTitle( "Events" );

  // TCanvas *c_mttvis   = new TCanvas( " " ," "  );
  // c_mttvis->cd(); 
  // h_mttvis->Draw(); 
  // h_mttvis->GetYaxis()->SetTitle( "Events" );

  // TCanvas *c_mt1   = new TCanvas( "mt1 " ,"mt1 "  );
  // c_mt1->cd(); 
  // h_mt1->Draw(); 
  // h_mt1->GetYaxis()->SetTitle( "Events" );

  // TCanvas *c_mt2   = new TCanvas( "mt2 " ,"mt2 "  );
  // c_mt2->cd(); 
  // h_mt2->Draw(); 
  // h_mt2->GetYaxis()->SetTitle( "Events" );

  // TCanvas *c_mt3   = new TCanvas( "mt3 " ,"mt3 "  );
  // c_mt3->cd(); 
  // h_mt3->Draw(); 
  // h_mt3->GetYaxis()->SetTitle( "Events" );

  // TCanvas *c_mct1   = new TCanvas( "mct1 " ,"mct1 "  );
  // c_mct1->cd(); 
  // h_mct1->Draw(); 
  // h_mct1->GetYaxis()->SetTitle( "Events" );

  // TCanvas *c_mct2   = new TCanvas( "mct2 " ,"mct2 "  );
  // c_mct2->cd(); 
  // h_mct2->Draw(); 
  // h_mct2->GetYaxis()->SetTitle( "Events" );

  // TCanvas *c_mct3   = new TCanvas( "mct3 " ,"mct3 "  );
  // c_mct3->cd(); 
  // h_mct3->Draw(); 
  // h_mct3->GetYaxis()->SetTitle( "Events" );

  // TCanvas *c_mlt   = new TCanvas( "mlt " ,"mlt "  );
  // c_mlt->cd(); 
  // h_mlt->Draw(); 
  // h_mlt->GetYaxis()->SetTitle( "Events" );

  // TCanvas *c_mlct   = new TCanvas( "mlct " ,"mlct "  );
  // c_mlct->cd(); 
  // h_mlct->Draw(); 
  // h_mlct->GetYaxis()->SetTitle( "Events" );

}

bool AnalysisTransverse::PassCuts() const
{
  // if( this->PassCuts_NJets() )
  // {
  //   if( this->PassCuts_MET() )
  //   {
  //     if( this->PassCuts_MWT() )
  //     {
        return true;
  //     }
  //   }
  // }
  // return false;
}

// bool AnalysisTransverse::PassCuts_NJets() const
// {
//   // if( m_ntup->jet_n() >= 4 )
//   // {
//     return true;
//   // }
//   // return false;
// }

// bool AnalysisTransverse::PassCuts_MET() const
// {
//   if( m_channel == "el" )
//   {
//     if( m_ntup->met_met() > 30.0 * m_GeV )
//     {
//       return true;
//     }
//   }
//   if( m_channel == "mu" )
//   {
//     if( m_ntup->met_met() > 20.0 * m_GeV )
//     {
//       return true;
//     }
//   }  
//   return false;
// }

// bool AnalysisTransverse::PassCuts_MWT() const
// {
//   if( m_channel == "el" )
//   {
//     if( m_ntup->met_met() > 30.0 * m_GeV ) // changed from MWT!!!!!!!!!!!!!
//     {
//       return true;
//     }
//   }
//   if( m_channel == "mu" ){
//     if( ( m_ntup->met_met() + m_ntup->met_met() ) > 60.0 * m_GeV ) // changed from MWT!!!!!!!!!
//     {             ///met_mwt          mwt
//       return true;
//     }
//   }  
//   return false;  
// }

void AnalysisTransverse::Loop()
{
  // Loop over all files
  for(Itr_s i = m_inputFiles->begin(); i != m_inputFiles->end(); ++i)
  {
    cout<<"  Processing File = "<<(*i)<<endl;
    
    this->SetupTreesForNewFile( (*i) );
 
    // The Event Loop
    Long64_t nEvents = this->TotalEvents();
    for(Long64_t jentry=0; jentry<nEvents;++jentry) 
    {
      Long64_t ientry = this->IncrementEvent(jentry);
      if (ientry < 0) break;
      this->EachEvent();
    }
    
    this->CleanUp();
  }
}  

AnalysisTransverse::~AnalysisTransverse()
{       
  delete m_inputFiles;
}

void AnalysisTransverse::SetupInputFiles()
{   
  m_inputFiles = new vector<TString>;
  TString base("/afs/cern.ch/work/d/demillar/Ntuples_Zprime/");
  // TString base("/Users/declan/Data/Ntuples_Zprime/");
  
  // m_inputFiles->push_back(base + "SSM-m2000-w200_Zp_13_2to6_DL_1x5000000.root");
  m_inputFiles->push_back(base + m_inputFileName);
}

Long64_t AnalysisTransverse::TotalEvents()
{
  // Internal for Event Looping
  if(m_ntup != 0){return m_ntup->totalEvents();}
  return -9999;  
}

Long64_t AnalysisTransverse::IncrementEvent(Long64_t i)
{
  Long64_t ev(-1);
  if(m_ntup != 0){ev = m_ntup->LoadTree(i);}
  return ev;  
}

void AnalysisTransverse::SetupTreesForNewFile(const TString& s)
{
  TString treeToUse = "RootTuple";
  
  m_chainNtup = new TChain(treeToUse,"");
  TString TStringNtuple = s + "/" + treeToUse;
  m_chainNtup->Add(TStringNtuple);
  m_ntup = new RootTuple(m_chainNtup);  
}

void AnalysisTransverse::CleanUp()
{
  delete m_chainNtup;
  delete m_ntup;  
}
