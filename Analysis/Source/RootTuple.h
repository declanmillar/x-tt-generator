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
//  *   File Generated on Mon May 11 12:24:29 2015                            * 
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
    inline Double_t  Delta_y() const {b_Delta_y->GetEntry(m_currentEvent);return m_Delta_y;} 
    inline vector<double>*  E() const {b_E->GetEntry(m_currentEvent);return m_E;} 
    inline Double_t  Mtt() const {b_Mtt->GetEntry(m_currentEvent);return m_Mtt;} 
    inline vector<double>*  Px() const {b_Px->GetEntry(m_currentEvent);return m_Px;} 
    inline vector<double>*  Py() const {b_Py->GetEntry(m_currentEvent);return m_Py;} 
    inline vector<double>*  Pz() const {b_Pz->GetEntry(m_currentEvent);return m_Pz;} 
    inline vector<int>*  barcode() const {b_barcode->GetEntry(m_currentEvent);return m_barcode;} 
    inline Double_t  beta() const {b_beta->GetEntry(m_currentEvent);return m_beta;} 
    inline Double_t  costhetastar() const {b_costhetastar->GetEntry(m_currentEvent);return m_costhetastar;} 
    inline Double_t  costhetat() const {b_costhetat->GetEntry(m_currentEvent);return m_costhetat;} 
    inline Double_t  costhetatcol() const {b_costhetatcol->GetEntry(m_currentEvent);return m_costhetatcol;} 
    inline Int_t  iteration() const {b_iteration->GetEntry(m_currentEvent);return m_iteration;} 
    inline Double_t  weight() const {b_weight->GetEntry(m_currentEvent);return m_weight;} 
    inline Double_t  weightLL() const {b_weightLL->GetEntry(m_currentEvent);return m_weightLL;} 
    inline Double_t  weightLR() const {b_weightLR->GetEntry(m_currentEvent);return m_weightLR;} 
    inline Double_t  weightRL() const {b_weightRL->GetEntry(m_currentEvent);return m_weightRL;} 
    inline Double_t  weightRR() const {b_weightRR->GetEntry(m_currentEvent);return m_weightRR;} 
    inline Double_t  weight_ee() const {b_weight_ee->GetEntry(m_currentEvent);return m_weight_ee;} 
    inline Double_t  weight_emu() const {b_weight_emu->GetEntry(m_currentEvent);return m_weight_emu;} 
    inline Double_t  weight_eq() const {b_weight_eq->GetEntry(m_currentEvent);return m_weight_eq;} 
    inline Double_t  weight_qq() const {b_weight_qq->GetEntry(m_currentEvent);return m_weight_qq;} 
    inline Double_t  yt() const {b_yt->GetEntry(m_currentEvent);return m_yt;} 
    inline Double_t  ytbar() const {b_ytbar->GetEntry(m_currentEvent);return m_ytbar;} 
    inline Double_t  ytt() const {b_ytt->GetEntry(m_currentEvent);return m_ytt;} 

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

    Double_t  m_Delta_y; 
    vector<double>*  m_E; 
    Double_t  m_Mtt; 
    vector<double>*  m_Px; 
    vector<double>*  m_Py; 
    vector<double>*  m_Pz; 
    vector<int>*  m_barcode; 
    Double_t  m_beta; 
    Double_t  m_costhetastar; 
    Double_t  m_costhetat; 
    Double_t  m_costhetatcol; 
    Int_t  m_iteration; 
    Double_t  m_weight; 
    Double_t  m_weightLL; 
    Double_t  m_weightLR; 
    Double_t  m_weightRL; 
    Double_t  m_weightRR; 
    Double_t  m_weight_ee; 
    Double_t  m_weight_emu; 
    Double_t  m_weight_eq; 
    Double_t  m_weight_qq; 
    Double_t  m_yt; 
    Double_t  m_ytbar; 
    Double_t  m_ytt; 

    TBranch*  b_Delta_y; 
    TBranch*  b_E; 
    TBranch*  b_Mtt; 
    TBranch*  b_Px; 
    TBranch*  b_Py; 
    TBranch*  b_Pz; 
    TBranch*  b_barcode; 
    TBranch*  b_beta; 
    TBranch*  b_costhetastar; 
    TBranch*  b_costhetat; 
    TBranch*  b_costhetatcol; 
    TBranch*  b_iteration; 
    TBranch*  b_weight; 
    TBranch*  b_weightLL; 
    TBranch*  b_weightLR; 
    TBranch*  b_weightRL; 
    TBranch*  b_weightRR; 
    TBranch*  b_weight_ee; 
    TBranch*  b_weight_emu; 
    TBranch*  b_weight_eq; 
    TBranch*  b_weight_qq; 
    TBranch*  b_yt; 
    TBranch*  b_ytbar; 
    TBranch*  b_ytt; 
}; 
#endif 

