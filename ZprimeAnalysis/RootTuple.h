//  *************************************************************************** 
//  *                                                                         * 
//  *   This program is free software; you can redistribute it and/or modify  * 
//  *   it under the terms of the GNU General Public License as published by  * 
//  *   the Free Software Foundation; either version 2 of the License, or     * 
//  *   (at your option) any later version.                                   * 
//  *                                                                         * 
//  *   Author: John Morris (john.morris@cern.ch)                             * 
//  *           Queen Mary University of London                               * 
//  *   Editor: Declan Millar (declan.millar@cern.ch)                         * 
//  *   File Generated on Thu Apr 23 17:05:45 2015                            * 
//  *                                                                         * 
//  ***************************************************************************/ 

// This class is for accessing the RootTuple Ntuples  
// You should not have to edit this file 

#ifndef _NTUPLE_ROOTTUPLE_H_ 
#define _NTUPLE_ROOTTUPLE_H_ 

#include <TROOT.h> 
#include <TChain.h> 
#include <TFile.h> 
#include <vector> 
using std::vector; 
#include <string> 
using std::string; 
#include <map> 
using std::map; 
#include <iostream> 
using std::cout; 
using std::endl; 

class RootTuple{ 
  public: 
    explicit RootTuple(TTree* tree); 
    RootTuple(TTree* tree,const bool& isMC,const bool& isAFII); 
    virtual ~RootTuple(); 
    Long64_t totalEvents(); 
    Long64_t LoadTree(Long64_t entry); 

    // public inline member functions -- Use these to get access to the TTree variables 
    inline vector<double>*  E() const {b_E->GetEntry(m_currentEvent);return m_E;} 
    inline vector<double>*  Px() const {b_Px->GetEntry(m_currentEvent);return m_Px;} 
    inline vector<double>*  Py() const {b_Py->GetEntry(m_currentEvent);return m_Py;} 
    inline vector<double>*  Pz() const {b_Pz->GetEntry(m_currentEvent);return m_Pz;} 
    inline vector<int>*  barcode() const {b_barcode->GetEntry(m_currentEvent);return m_barcode;} 
    inline Double_t  beta() const {b_beta->GetEntry(m_currentEvent);return m_beta;} 
    inline Double_t  cost() const {b_cost->GetEntry(m_currentEvent);return m_cost;} 
    inline Double_t  delta_y() const {b_delta_y->GetEntry(m_currentEvent);return m_delta_y;} 
    inline Double_t  et() const {b_et->GetEntry(m_currentEvent);return m_et;} 
    inline Double_t  etat() const {b_etat->GetEntry(m_currentEvent);return m_etat;} 
    inline Double_t  etatb() const {b_etatb->GetEntry(m_currentEvent);return m_etatb;} 
    inline Double_t  mtt() const {b_mtt->GetEntry(m_currentEvent);return m_mtt;} 
    inline Double_t  phit() const {b_phit->GetEntry(m_currentEvent);return m_phit;} 
    inline Double_t  phitb() const {b_phitb->GetEntry(m_currentEvent);return m_phitb;} 
    inline Double_t  ptt() const {b_ptt->GetEntry(m_currentEvent);return m_ptt;} 
    inline Double_t  pttb() const {b_pttb->GetEntry(m_currentEvent);return m_pttb;} 
    inline Double_t  weight() const {b_weight->GetEntry(m_currentEvent);return m_weight;} 
    inline Double_t  weightLL() const {b_weightLL->GetEntry(m_currentEvent);return m_weightLL;} 
    inline Double_t  weightLR() const {b_weightLR->GetEntry(m_currentEvent);return m_weightLR;} 
    inline Double_t  weightRL() const {b_weightRL->GetEntry(m_currentEvent);return m_weightRL;} 
    inline Double_t  weightRR() const {b_weightRR->GetEntry(m_currentEvent);return m_weightRR;} 
    inline Double_t  ycolt() const {b_ycolt->GetEntry(m_currentEvent);return m_ycolt;} 
    inline Double_t  ycoltb() const {b_ycoltb->GetEntry(m_currentEvent);return m_ycoltb;} 

    inline Long64_t currentEvent() const {return m_currentEvent;} 

  protected: 
    Int_t    GetEntry(Long64_t entry); 
    void     Init(TTree *tree); 

  private: 
    RootTuple(); 
    RootTuple(const RootTuple& rhs);  
    void operator=(const RootTuple& rhs); 

    bool m_isMC; 
    bool m_isAFII; 

    TTree          *fChain; 
    int             fCurrent; 

    Long64_t m_currentEvent; 

    vector<double>*  m_E; 
    vector<double>*  m_Px; 
    vector<double>*  m_Py; 
    vector<double>*  m_Pz; 
    vector<int>*  m_barcode; 
    Double_t  m_beta; 
    Double_t  m_cost; 
    Double_t  m_delta_y; 
    Double_t  m_et; 
    Double_t  m_etat; 
    Double_t  m_etatb; 
    Double_t  m_mtt; 
    Double_t  m_phit; 
    Double_t  m_phitb; 
    Double_t  m_ptt; 
    Double_t  m_pttb; 
    Double_t  m_weight; 
    Double_t  m_weightLL; 
    Double_t  m_weightLR; 
    Double_t  m_weightRL; 
    Double_t  m_weightRR; 
    Double_t  m_ycolt; 
    Double_t  m_ycoltb; 

    TBranch*  b_E; 
    TBranch*  b_Px; 
    TBranch*  b_Py; 
    TBranch*  b_Pz; 
    TBranch*  b_barcode; 
    TBranch*  b_beta; 
    TBranch*  b_cost; 
    TBranch*  b_delta_y; 
    TBranch*  b_et; 
    TBranch*  b_etat; 
    TBranch*  b_etatb; 
    TBranch*  b_mtt; 
    TBranch*  b_phit; 
    TBranch*  b_phitb; 
    TBranch*  b_ptt; 
    TBranch*  b_pttb; 
    TBranch*  b_weight; 
    TBranch*  b_weightLL; 
    TBranch*  b_weightLR; 
    TBranch*  b_weightRL; 
    TBranch*  b_weightRR; 
    TBranch*  b_ycolt; 
    TBranch*  b_ycoltb; 
}; 
#endif 

