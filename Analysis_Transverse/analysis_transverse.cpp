#include "analysis_transverse.h"

AnalysisTransverse::AnalysisTransverse(const TString& outputFileName) :
  m_pi(3.14159265),
  m_GeV(1000.0),
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

  }
}

void AnalysisTransverse::PreLoop()
{
  this->SetupInputFiles();
  
  // Define output file 
  m_outputFile = new TFile(m_outputFileName,"RECREATE");
  
  // Define Histograms  
  h_mtt = new TH1D( "mtt", "mtt", 100, 0.0, 14000.0);

}

void AnalysisTransverse::PostLoop()
{
  this->MakeGraphs();
  
  m_outputFile->cd();
  m_outputFile->cd("/");

  // Save histograms
  h_mtt->Write();

  m_outputFile->Close();
  
}

void AnalysisTransverse::MakeGraphs()
{
  printf("Making Graphs...\n");

  TCanvas *c_mtt   = new TCanvas( "mtt " ,"mtt "  );
  c_mtt->cd(); 
  h_mtt->Draw(); 
  h_mtt->GetYaxis()->SetTitle( "Events" );
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
  // TString base("/afs/cern.ch/work/d/demillar/Ntuples_Zprime/");
  TString base("/Users/declan/Data/Ntuples_Zprime/");
  
  m_inputFiles->push_back(base + "SM_13_2to6_DL_1x500000.root");
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
