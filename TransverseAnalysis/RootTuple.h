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
//  *   File Generated on Fri Apr 24 17:09:51 2015                            * 
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
    inline Double_t  cosfl() const {b_cosfl->GetEntry(m_currentEvent);return m_cosfl;} 
    inline Double_t  cost() const {b_cost->GetEntry(m_currentEvent);return m_cost;} 
    inline Double_t  cost5() const {b_cost5->GetEntry(m_currentEvent);return m_cost5;} 
    inline Double_t  cost7() const {b_cost7->GetEntry(m_currentEvent);return m_cost7;} 
    inline Double_t  ct7ct5() const {b_ct7ct5->GetEntry(m_currentEvent);return m_ct7ct5;} 
    inline Double_t  delta_y() const {b_delta_y->GetEntry(m_currentEvent);return m_delta_y;} 
    inline Double_t  dphi() const {b_dphi->GetEntry(m_currentEvent);return m_dphi;} 
    inline Double_t  et() const {b_et->GetEntry(m_currentEvent);return m_et;} 
    inline Double_t  etab() const {b_etab->GetEntry(m_currentEvent);return m_etab;} 
    inline Double_t  etabb() const {b_etabb->GetEntry(m_currentEvent);return m_etabb;} 
    inline Double_t  etalm() const {b_etalm->GetEntry(m_currentEvent);return m_etalm;} 
    inline Double_t  etalp() const {b_etalp->GetEntry(m_currentEvent);return m_etalp;} 
    inline Double_t  etanu() const {b_etanu->GetEntry(m_currentEvent);return m_etanu;} 
    inline Double_t  etanub() const {b_etanub->GetEntry(m_currentEvent);return m_etanub;} 
    inline Double_t  etat() const {b_etat->GetEntry(m_currentEvent);return m_etat;} 
    inline Double_t  etatb() const {b_etatb->GetEntry(m_currentEvent);return m_etatb;} 
    inline Double_t  etmiss() const {b_etmiss->GetEntry(m_currentEvent);return m_etmiss;} 
    inline Double_t  fl() const {b_fl->GetEntry(m_currentEvent);return m_fl;} 
    inline Double_t  ht() const {b_ht->GetEntry(m_currentEvent);return m_ht;} 
    inline Double_t  mct1() const {b_mct1->GetEntry(m_currentEvent);return m_mct1;} 
    inline Double_t  mct2() const {b_mct2->GetEntry(m_currentEvent);return m_mct2;} 
    inline Double_t  mct3() const {b_mct3->GetEntry(m_currentEvent);return m_mct3;} 
    inline Double_t  mlct() const {b_mlct->GetEntry(m_currentEvent);return m_mlct;} 
    inline Double_t  mll() const {b_mll->GetEntry(m_currentEvent);return m_mll;} 
    inline Double_t  mlt() const {b_mlt->GetEntry(m_currentEvent);return m_mlt;} 
    inline Double_t  mt1() const {b_mt1->GetEntry(m_currentEvent);return m_mt1;} 
    inline Double_t  mt2() const {b_mt2->GetEntry(m_currentEvent);return m_mt2;} 
    inline Double_t  mt3() const {b_mt3->GetEntry(m_currentEvent);return m_mt3;} 
    inline Double_t  mtt() const {b_mtt->GetEntry(m_currentEvent);return m_mtt;} 
    inline Double_t  mttvis() const {b_mttvis->GetEntry(m_currentEvent);return m_mttvis;} 
    inline Double_t  phib() const {b_phib->GetEntry(m_currentEvent);return m_phib;} 
    inline Double_t  phibb() const {b_phibb->GetEntry(m_currentEvent);return m_phibb;} 
    inline Double_t  philm() const {b_philm->GetEntry(m_currentEvent);return m_philm;} 
    inline Double_t  philp() const {b_philp->GetEntry(m_currentEvent);return m_philp;} 
    inline Double_t  phinu() const {b_phinu->GetEntry(m_currentEvent);return m_phinu;} 
    inline Double_t  phinub() const {b_phinub->GetEntry(m_currentEvent);return m_phinub;} 
    inline Double_t  phit() const {b_phit->GetEntry(m_currentEvent);return m_phit;} 
    inline Double_t  phitb() const {b_phitb->GetEntry(m_currentEvent);return m_phitb;} 
    inline Double_t  ptb() const {b_ptb->GetEntry(m_currentEvent);return m_ptb;} 
    inline Double_t  ptbb() const {b_ptbb->GetEntry(m_currentEvent);return m_ptbb;} 
    inline Double_t  ptlm() const {b_ptlm->GetEntry(m_currentEvent);return m_ptlm;} 
    inline Double_t  ptlp() const {b_ptlp->GetEntry(m_currentEvent);return m_ptlp;} 
    inline Double_t  ptnu() const {b_ptnu->GetEntry(m_currentEvent);return m_ptnu;} 
    inline Double_t  ptnub() const {b_ptnub->GetEntry(m_currentEvent);return m_ptnub;} 
    inline Double_t  ptt() const {b_ptt->GetEntry(m_currentEvent);return m_ptt;} 
    inline Double_t  pttb() const {b_pttb->GetEntry(m_currentEvent);return m_pttb;} 
    inline Double_t  weight() const {b_weight->GetEntry(m_currentEvent);return m_weight;} 
    inline Double_t  weightLL() const {b_weightLL->GetEntry(m_currentEvent);return m_weightLL;} 
    inline Double_t  weightLR() const {b_weightLR->GetEntry(m_currentEvent);return m_weightLR;} 
    inline Double_t  weightRL() const {b_weightRL->GetEntry(m_currentEvent);return m_weightRL;} 
    inline Double_t  weightRR() const {b_weightRR->GetEntry(m_currentEvent);return m_weightRR;} 
    inline Double_t  ycolb() const {b_ycolb->GetEntry(m_currentEvent);return m_ycolb;} 
    inline Double_t  ycolbb() const {b_ycolbb->GetEntry(m_currentEvent);return m_ycolbb;} 
    inline Double_t  ycollm() const {b_ycollm->GetEntry(m_currentEvent);return m_ycollm;} 
    inline Double_t  ycollp() const {b_ycollp->GetEntry(m_currentEvent);return m_ycollp;} 
    inline Double_t  ycolnu() const {b_ycolnu->GetEntry(m_currentEvent);return m_ycolnu;} 
    inline Double_t  ycolnub() const {b_ycolnub->GetEntry(m_currentEvent);return m_ycolnub;} 
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
    Double_t  m_cosfl; 
    Double_t  m_cost; 
    Double_t  m_cost5; 
    Double_t  m_cost7; 
    Double_t  m_ct7ct5; 
    Double_t  m_delta_y; 
    Double_t  m_dphi; 
    Double_t  m_et; 
    Double_t  m_etab; 
    Double_t  m_etabb; 
    Double_t  m_etalm; 
    Double_t  m_etalp; 
    Double_t  m_etanu; 
    Double_t  m_etanub; 
    Double_t  m_etat; 
    Double_t  m_etatb; 
    Double_t  m_etmiss; 
    Double_t  m_fl; 
    Double_t  m_ht; 
    Double_t  m_mct1; 
    Double_t  m_mct2; 
    Double_t  m_mct3; 
    Double_t  m_mlct; 
    Double_t  m_mll; 
    Double_t  m_mlt; 
    Double_t  m_mt1; 
    Double_t  m_mt2; 
    Double_t  m_mt3; 
    Double_t  m_mtt; 
    Double_t  m_mttvis; 
    Double_t  m_phib; 
    Double_t  m_phibb; 
    Double_t  m_philm; 
    Double_t  m_philp; 
    Double_t  m_phinu; 
    Double_t  m_phinub; 
    Double_t  m_phit; 
    Double_t  m_phitb; 
    Double_t  m_ptb; 
    Double_t  m_ptbb; 
    Double_t  m_ptlm; 
    Double_t  m_ptlp; 
    Double_t  m_ptnu; 
    Double_t  m_ptnub; 
    Double_t  m_ptt; 
    Double_t  m_pttb; 
    Double_t  m_weight; 
    Double_t  m_weightLL; 
    Double_t  m_weightLR; 
    Double_t  m_weightRL; 
    Double_t  m_weightRR; 
    Double_t  m_ycolb; 
    Double_t  m_ycolbb; 
    Double_t  m_ycollm; 
    Double_t  m_ycollp; 
    Double_t  m_ycolnu; 
    Double_t  m_ycolnub; 
    Double_t  m_ycolt; 
    Double_t  m_ycoltb; 

    TBranch*  b_E; 
    TBranch*  b_Px; 
    TBranch*  b_Py; 
    TBranch*  b_Pz; 
    TBranch*  b_barcode; 
    TBranch*  b_beta; 
    TBranch*  b_cosfl; 
    TBranch*  b_cost; 
    TBranch*  b_cost5; 
    TBranch*  b_cost7; 
    TBranch*  b_ct7ct5; 
    TBranch*  b_delta_y; 
    TBranch*  b_dphi; 
    TBranch*  b_et; 
    TBranch*  b_etab; 
    TBranch*  b_etabb; 
    TBranch*  b_etalm; 
    TBranch*  b_etalp; 
    TBranch*  b_etanu; 
    TBranch*  b_etanub; 
    TBranch*  b_etat; 
    TBranch*  b_etatb; 
    TBranch*  b_etmiss; 
    TBranch*  b_fl; 
    TBranch*  b_ht; 
    TBranch*  b_mct1; 
    TBranch*  b_mct2; 
    TBranch*  b_mct3; 
    TBranch*  b_mlct; 
    TBranch*  b_mll; 
    TBranch*  b_mlt; 
    TBranch*  b_mt1; 
    TBranch*  b_mt2; 
    TBranch*  b_mt3; 
    TBranch*  b_mtt; 
    TBranch*  b_mttvis; 
    TBranch*  b_phib; 
    TBranch*  b_phibb; 
    TBranch*  b_philm; 
    TBranch*  b_philp; 
    TBranch*  b_phinu; 
    TBranch*  b_phinub; 
    TBranch*  b_phit; 
    TBranch*  b_phitb; 
    TBranch*  b_ptb; 
    TBranch*  b_ptbb; 
    TBranch*  b_ptlm; 
    TBranch*  b_ptlp; 
    TBranch*  b_ptnu; 
    TBranch*  b_ptnub; 
    TBranch*  b_ptt; 
    TBranch*  b_pttb; 
    TBranch*  b_weight; 
    TBranch*  b_weightLL; 
    TBranch*  b_weightLR; 
    TBranch*  b_weightRL; 
    TBranch*  b_weightRR; 
    TBranch*  b_ycolb; 
    TBranch*  b_ycolbb; 
    TBranch*  b_ycollm; 
    TBranch*  b_ycollp; 
    TBranch*  b_ycolnu; 
    TBranch*  b_ycolnub; 
    TBranch*  b_ycolt; 
    TBranch*  b_ycoltb; 
}; 
#endif 

