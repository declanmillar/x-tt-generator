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

#include "RootTuple.h" 

RootTuple::RootTuple(TTree* tree) : 
  m_isMC(false), 
  m_isAFII(false) 
{ 
  m_currentEvent = 0; 
  this->Init(tree); 
} 

RootTuple::RootTuple(TTree* tree,const bool& isMC,const bool& isAFII) : 
  m_isMC(isMC), 
  m_isAFII(isAFII) 
{ 
  m_currentEvent = 0; 
  this->Init(tree); 
} 

RootTuple::~RootTuple(){ 
  if (!fChain) return; 
  delete m_E; 
  delete m_Px; 
  delete m_Py; 
  delete m_Pz; 
  delete m_barcode; 
} 

int RootTuple::GetEntry(Long64_t entry){ 
  if (!fChain) return 0; 
  return fChain->GetEntry(entry); 
} 

Long64_t RootTuple::totalEvents(){ 
  return fChain->GetEntriesFast(); 
} 

Long64_t RootTuple::LoadTree(Long64_t entry){ 
  m_currentEvent = entry; 
  if (!fChain) return -5; 
  Long64_t centry = fChain->LoadTree(entry); 
  if (centry < 0) return centry; 
  if (!fChain->InheritsFrom(TChain::Class()))  return centry; 
  TChain *chain = (TChain*)fChain; 
  if (chain->GetTreeNumber() != fCurrent) { 
    fCurrent = chain->GetTreeNumber(); 
  } 
  return centry; 
} 

void RootTuple::Init(TTree* tree){ 
  m_E = 0; 
  m_Px = 0; 
  m_Py = 0; 
  m_Pz = 0; 
  m_barcode = 0; 
  m_beta = 0; 
  m_cosfl = 0; 
  m_cost = 0; 
  m_cost5 = 0; 
  m_cost7 = 0; 
  m_ct7ct5 = 0; 
  m_delta_y = 0; 
  m_dphi = 0; 
  m_et = 0; 
  m_etab = 0; 
  m_etabb = 0; 
  m_etalm = 0; 
  m_etalp = 0; 
  m_etanu = 0; 
  m_etanub = 0; 
  m_etat = 0; 
  m_etatb = 0; 
  m_etmiss = 0; 
  m_fl = 0; 
  m_ht = 0; 
  m_mct1 = 0; 
  m_mct2 = 0; 
  m_mct3 = 0; 
  m_mlct = 0; 
  m_mll = 0; 
  m_mlt = 0; 
  m_mt1 = 0; 
  m_mt2 = 0; 
  m_mt3 = 0; 
  m_mtt = 0; 
  m_mttvis = 0; 
  m_phib = 0; 
  m_phibb = 0; 
  m_philm = 0; 
  m_philp = 0; 
  m_phinu = 0; 
  m_phinub = 0; 
  m_phit = 0; 
  m_phitb = 0; 
  m_ptb = 0; 
  m_ptbb = 0; 
  m_ptlm = 0; 
  m_ptlp = 0; 
  m_ptnu = 0; 
  m_ptnub = 0; 
  m_ptt = 0; 
  m_pttb = 0; 
  m_weight = 0; 
  m_weightLL = 0; 
  m_weightLR = 0; 
  m_weightRL = 0; 
  m_weightRR = 0; 
  m_ycolb = 0; 
  m_ycolbb = 0; 
  m_ycollm = 0; 
  m_ycollp = 0; 
  m_ycolnu = 0; 
  m_ycolnub = 0; 
  m_ycolt = 0; 
  m_ycoltb = 0; 

  if (!tree) return; 
  fChain = tree; 
  fCurrent = -1; 
  fChain->SetMakeClass(1); 

  fChain->SetBranchAddress("E", &m_E, &b_E); 
  fChain->SetBranchAddress("Px", &m_Px, &b_Px); 
  fChain->SetBranchAddress("Py", &m_Py, &b_Py); 
  fChain->SetBranchAddress("Pz", &m_Pz, &b_Pz); 
  fChain->SetBranchAddress("barcode", &m_barcode, &b_barcode); 
  fChain->SetBranchAddress("beta", &m_beta, &b_beta); 
  fChain->SetBranchAddress("cosfl", &m_cosfl, &b_cosfl); 
  fChain->SetBranchAddress("cost", &m_cost, &b_cost); 
  fChain->SetBranchAddress("cost5", &m_cost5, &b_cost5); 
  fChain->SetBranchAddress("cost7", &m_cost7, &b_cost7); 
  fChain->SetBranchAddress("ct7ct5", &m_ct7ct5, &b_ct7ct5); 
  fChain->SetBranchAddress("delta_y", &m_delta_y, &b_delta_y); 
  fChain->SetBranchAddress("dphi", &m_dphi, &b_dphi); 
  fChain->SetBranchAddress("et", &m_et, &b_et); 
  fChain->SetBranchAddress("etab", &m_etab, &b_etab); 
  fChain->SetBranchAddress("etabb", &m_etabb, &b_etabb); 
  fChain->SetBranchAddress("etalm", &m_etalm, &b_etalm); 
  fChain->SetBranchAddress("etalp", &m_etalp, &b_etalp); 
  fChain->SetBranchAddress("etanu", &m_etanu, &b_etanu); 
  fChain->SetBranchAddress("etanub", &m_etanub, &b_etanub); 
  fChain->SetBranchAddress("etat", &m_etat, &b_etat); 
  fChain->SetBranchAddress("etatb", &m_etatb, &b_etatb); 
  fChain->SetBranchAddress("etmiss", &m_etmiss, &b_etmiss); 
  fChain->SetBranchAddress("fl", &m_fl, &b_fl); 
  fChain->SetBranchAddress("ht", &m_ht, &b_ht); 
  fChain->SetBranchAddress("mct1", &m_mct1, &b_mct1); 
  fChain->SetBranchAddress("mct2", &m_mct2, &b_mct2); 
  fChain->SetBranchAddress("mct3", &m_mct3, &b_mct3); 
  fChain->SetBranchAddress("mlct", &m_mlct, &b_mlct); 
  fChain->SetBranchAddress("mll", &m_mll, &b_mll); 
  fChain->SetBranchAddress("mlt", &m_mlt, &b_mlt); 
  fChain->SetBranchAddress("mt1", &m_mt1, &b_mt1); 
  fChain->SetBranchAddress("mt2", &m_mt2, &b_mt2); 
  fChain->SetBranchAddress("mt3", &m_mt3, &b_mt3); 
  fChain->SetBranchAddress("mtt", &m_mtt, &b_mtt); 
  fChain->SetBranchAddress("mttvis", &m_mttvis, &b_mttvis); 
  fChain->SetBranchAddress("phib", &m_phib, &b_phib); 
  fChain->SetBranchAddress("phibb", &m_phibb, &b_phibb); 
  fChain->SetBranchAddress("philm", &m_philm, &b_philm); 
  fChain->SetBranchAddress("philp", &m_philp, &b_philp); 
  fChain->SetBranchAddress("phinu", &m_phinu, &b_phinu); 
  fChain->SetBranchAddress("phinub", &m_phinub, &b_phinub); 
  fChain->SetBranchAddress("phit", &m_phit, &b_phit); 
  fChain->SetBranchAddress("phitb", &m_phitb, &b_phitb); 
  fChain->SetBranchAddress("ptb", &m_ptb, &b_ptb); 
  fChain->SetBranchAddress("ptbb", &m_ptbb, &b_ptbb); 
  fChain->SetBranchAddress("ptlm", &m_ptlm, &b_ptlm); 
  fChain->SetBranchAddress("ptlp", &m_ptlp, &b_ptlp); 
  fChain->SetBranchAddress("ptnu", &m_ptnu, &b_ptnu); 
  fChain->SetBranchAddress("ptnub", &m_ptnub, &b_ptnub); 
  fChain->SetBranchAddress("ptt", &m_ptt, &b_ptt); 
  fChain->SetBranchAddress("pttb", &m_pttb, &b_pttb); 
  fChain->SetBranchAddress("weight", &m_weight, &b_weight); 
  fChain->SetBranchAddress("weightLL", &m_weightLL, &b_weightLL); 
  fChain->SetBranchAddress("weightLR", &m_weightLR, &b_weightLR); 
  fChain->SetBranchAddress("weightRL", &m_weightRL, &b_weightRL); 
  fChain->SetBranchAddress("weightRR", &m_weightRR, &b_weightRR); 
  fChain->SetBranchAddress("ycolb", &m_ycolb, &b_ycolb); 
  fChain->SetBranchAddress("ycolbb", &m_ycolbb, &b_ycolbb); 
  fChain->SetBranchAddress("ycollm", &m_ycollm, &b_ycollm); 
  fChain->SetBranchAddress("ycollp", &m_ycollp, &b_ycollp); 
  fChain->SetBranchAddress("ycolnu", &m_ycolnu, &b_ycolnu); 
  fChain->SetBranchAddress("ycolnub", &m_ycolnub, &b_ycolnub); 
  fChain->SetBranchAddress("ycolt", &m_ycolt, &b_ycolt); 
  fChain->SetBranchAddress("ycoltb", &m_ycoltb, &b_ycoltb); 
} 

