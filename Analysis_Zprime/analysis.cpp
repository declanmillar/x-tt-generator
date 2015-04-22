#include "analysis.h"

AnalysisZprime::AnalysisZprime(const TString& outputFileName) :
  m_pi(3.14159265),
  m_GeV(1000.0),
  // m_channel(channel),  
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
  if( this->PassCuts() )
  {
    // cout << "Event = " << m_ntup->eventNumber()<<endl;
    
    float eventWeight = m_ntup->weight();
    
    // Fill Histograms
    if(m_ntup->barcode()->at(2) == -11)  h_positron_pz->Fill(m_ntup->Pz()->at(2), m_ntup->weight());

    h_mtt->Fill(m_ntup->Mtt(), eventWeight/h_mtt->GetXaxis()->GetBinWidth(1));

    // if(m_ntup->Pz()->size() > 0)
    // {
    //   for(unsigned int i = 0; i < m_ntup->Pz()->size(); ++i)
    //   {
    //     h_positron_pz->Fill( m_ntup->Pz()->at(i) / m_GeV, EventWeight );
    //   }
    // }
  }
}

void AnalysisZprime::PreLoop()
{
  this->SetupInputFiles();
  
  // Define output file 
  m_outputFile = new TFile(m_outputFileName,"RECREATE");
  
  // Define Histograms  

  h_positron_pz = new TH1F( "pz", "pz" ,50 ,0.0 , 10000.0);
  h_mtt = new TH1F( "mtt", "mtt" ,100 ,0.0 , 14000.0); 

  // h_positron_pz->Sumw2(); 
}

void AnalysisZprime::PostLoop()
{
  this->MakeGraphs();
  
  m_outputFile->cd();
  m_outputFile->cd("/");

  // Save histograms
  h_positron_pz->Write(); 
  h_mtt->Write();

  m_outputFile->Close();
  
}

void AnalysisZprime::MakeGraphs()
{
  printf("Making Graphs...\n");

  TCanvas *c_positron_pz   = new TCanvas( "positron_pz " ,"positron_pz "  );

  c_positron_pz->cd(); 
  h_positron_pz->Draw(); 

  h_positron_pz->GetYaxis()->SetTitle( "Events" );

  TCanvas *c_mtt   = new TCanvas( "mtt " ,"mtt "  );

  c_mtt->cd(); 
  h_mtt->Draw(); 

  h_mtt->GetYaxis()->SetTitle( "Events" );

}

bool AnalysisZprime::PassCuts() const
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

// bool AnalysisZprime::PassCuts_NJets() const
// {
//   // if( m_ntup->jet_n() >= 4 )
//   // {
//     return true;
//   // }
//   // return false;
// }

// bool AnalysisZprime::PassCuts_MET() const
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

// bool AnalysisZprime::PassCuts_MWT() const
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

void AnalysisZprime::Loop()
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

AnalysisZprime::~AnalysisZprime()
{       
  delete m_inputFiles;
}

void AnalysisZprime::SetupInputFiles()
{   
  m_inputFiles = new vector<TString>;
  TString base("../Ntuples_Zprime/");
  
  m_inputFiles->push_back(base + "SM_13_2to6_DL_1x500000.root");
}

Long64_t AnalysisZprime::TotalEvents()
{
  // Internal for Event Looping
  if(m_ntup != 0){return m_ntup->totalEvents();}
  return -9999;  
}

Long64_t AnalysisZprime::IncrementEvent(Long64_t i)
{
  Long64_t ev(-1);
  if(m_ntup != 0){ev = m_ntup->LoadTree(i);}
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