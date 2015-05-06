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
//  *   File Generated on Wed May  6 13:50:24 2015                            * 
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
  m_Delta_phi_l = 0; 
  m_Delta_y = 0; 
  m_Delta_y_reco = 0; 
  m_E = 0; 
  m_E6_reco = 0; 
  m_Et = 0; 
  m_MET = 0; 
  m_Mtt = 0; 
  m_Mtt_reco = 0; 
  m_Px = 0; 
  m_Py = 0; 
  m_Pz = 0; 
  m_Pz6_reco = 0; 
  m_barcode = 0; 
  m_beta = 0; 
  m_cosphil = 0; 
  m_costheta5 = 0; 
  m_costheta5cm = 0; 
  m_costheta5col = 0; 
  m_costheta7 = 0; 
  m_costheta7cm = 0; 
  m_costheta7col = 0; 
  m_costhetastar = 0; 
  m_costhetastar_reco = 0; 
  m_costhetat = 0; 
  m_costhetat_reco = 0; 
  m_costhetatcol = 0; 
  m_ct7ct5 = 0; 
  m_ct7ct5cm = 0; 
  m_ct7ct5col = 0; 
  m_eta3 = 0; 
  m_eta356 = 0; 
  m_eta4 = 0; 
  m_eta478 = 0; 
  m_eta5 = 0; 
  m_eta6 = 0; 
  m_eta7 = 0; 
  m_eta8 = 0; 
  m_fl = 0; 
  m_m478 = 0; 
  m_mt_reco = 0; 
  m_pT3 = 0; 
  m_pT356 = 0; 
  m_pT4 = 0; 
  m_pT478 = 0; 
  m_pT5 = 0; 
  m_pT6 = 0; 
  m_pT7 = 0; 
  m_pT8 = 0; 
  m_phi3 = 0; 
  m_phi356 = 0; 
  m_phi4 = 0; 
  m_phi478 = 0; 
  m_phi5 = 0; 
  m_phi6 = 0; 
  m_phi7 = 0; 
  m_phi8 = 0; 
  m_weight = 0; 
  m_weightLL = 0; 
  m_weightLR = 0; 
  m_weightRL = 0; 
  m_weightRR = 0; 
  m_weight_ee = 0; 
  m_weight_emu = 0; 
  m_weight_eq = 0; 
  m_weight_qq = 0; 
  m_xMET = 0; 
  m_yMET = 0; 
  m_ycol3 = 0; 
  m_ycol356 = 0; 
  m_ycol4 = 0; 
  m_ycol478 = 0; 
  m_ycol5 = 0; 
  m_ycol6 = 0; 
  m_ycol7 = 0; 
  m_ycol8 = 0; 
  m_yt = 0; 
  m_yt_reco = 0; 
  m_ytbar = 0; 
  m_ytt = 0; 

  if (!tree) return; 
  fChain = tree; 
  fCurrent = -1; 
  fChain->SetMakeClass(1); 

  fChain->SetBranchAddress("Delta_phi_l", &m_Delta_phi_l, &b_Delta_phi_l); 
  fChain->SetBranchAddress("Delta_y", &m_Delta_y, &b_Delta_y); 
  fChain->SetBranchAddress("Delta_y_reco", &m_Delta_y_reco, &b_Delta_y_reco); 
  fChain->SetBranchAddress("E", &m_E, &b_E); 
  fChain->SetBranchAddress("E6_reco", &m_E6_reco, &b_E6_reco); 
  fChain->SetBranchAddress("Et", &m_Et, &b_Et); 
  fChain->SetBranchAddress("MET", &m_MET, &b_MET); 
  fChain->SetBranchAddress("Mtt", &m_Mtt, &b_Mtt); 
  fChain->SetBranchAddress("Mtt_reco", &m_Mtt_reco, &b_Mtt_reco); 
  fChain->SetBranchAddress("Px", &m_Px, &b_Px); 
  fChain->SetBranchAddress("Py", &m_Py, &b_Py); 
  fChain->SetBranchAddress("Pz", &m_Pz, &b_Pz); 
  fChain->SetBranchAddress("Pz6_reco", &m_Pz6_reco, &b_Pz6_reco); 
  fChain->SetBranchAddress("barcode", &m_barcode, &b_barcode); 
  fChain->SetBranchAddress("beta", &m_beta, &b_beta); 
  fChain->SetBranchAddress("cosphil", &m_cosphil, &b_cosphil); 
  fChain->SetBranchAddress("costheta5", &m_costheta5, &b_costheta5); 
  fChain->SetBranchAddress("costheta5cm", &m_costheta5cm, &b_costheta5cm); 
  fChain->SetBranchAddress("costheta5col", &m_costheta5col, &b_costheta5col); 
  fChain->SetBranchAddress("costheta7", &m_costheta7, &b_costheta7); 
  fChain->SetBranchAddress("costheta7cm", &m_costheta7cm, &b_costheta7cm); 
  fChain->SetBranchAddress("costheta7col", &m_costheta7col, &b_costheta7col); 
  fChain->SetBranchAddress("costhetastar", &m_costhetastar, &b_costhetastar); 
  fChain->SetBranchAddress("costhetastar_reco", &m_costhetastar_reco, &b_costhetastar_reco); 
  fChain->SetBranchAddress("costhetat", &m_costhetat, &b_costhetat); 
  fChain->SetBranchAddress("costhetat_reco", &m_costhetat_reco, &b_costhetat_reco); 
  fChain->SetBranchAddress("costhetatcol", &m_costhetatcol, &b_costhetatcol); 
  fChain->SetBranchAddress("ct7ct5", &m_ct7ct5, &b_ct7ct5); 
  fChain->SetBranchAddress("ct7ct5cm", &m_ct7ct5cm, &b_ct7ct5cm); 
  fChain->SetBranchAddress("ct7ct5col", &m_ct7ct5col, &b_ct7ct5col); 
  fChain->SetBranchAddress("eta3", &m_eta3, &b_eta3); 
  fChain->SetBranchAddress("eta356", &m_eta356, &b_eta356); 
  fChain->SetBranchAddress("eta4", &m_eta4, &b_eta4); 
  fChain->SetBranchAddress("eta478", &m_eta478, &b_eta478); 
  fChain->SetBranchAddress("eta5", &m_eta5, &b_eta5); 
  fChain->SetBranchAddress("eta6", &m_eta6, &b_eta6); 
  fChain->SetBranchAddress("eta7", &m_eta7, &b_eta7); 
  fChain->SetBranchAddress("eta8", &m_eta8, &b_eta8); 
  fChain->SetBranchAddress("fl", &m_fl, &b_fl); 
  fChain->SetBranchAddress("m478", &m_m478, &b_m478); 
  fChain->SetBranchAddress("mt_reco", &m_mt_reco, &b_mt_reco); 
  fChain->SetBranchAddress("pT3", &m_pT3, &b_pT3); 
  fChain->SetBranchAddress("pT356", &m_pT356, &b_pT356); 
  fChain->SetBranchAddress("pT4", &m_pT4, &b_pT4); 
  fChain->SetBranchAddress("pT478", &m_pT478, &b_pT478); 
  fChain->SetBranchAddress("pT5", &m_pT5, &b_pT5); 
  fChain->SetBranchAddress("pT6", &m_pT6, &b_pT6); 
  fChain->SetBranchAddress("pT7", &m_pT7, &b_pT7); 
  fChain->SetBranchAddress("pT8", &m_pT8, &b_pT8); 
  fChain->SetBranchAddress("phi3", &m_phi3, &b_phi3); 
  fChain->SetBranchAddress("phi356", &m_phi356, &b_phi356); 
  fChain->SetBranchAddress("phi4", &m_phi4, &b_phi4); 
  fChain->SetBranchAddress("phi478", &m_phi478, &b_phi478); 
  fChain->SetBranchAddress("phi5", &m_phi5, &b_phi5); 
  fChain->SetBranchAddress("phi6", &m_phi6, &b_phi6); 
  fChain->SetBranchAddress("phi7", &m_phi7, &b_phi7); 
  fChain->SetBranchAddress("phi8", &m_phi8, &b_phi8); 
  fChain->SetBranchAddress("weight", &m_weight, &b_weight); 
  fChain->SetBranchAddress("weightLL", &m_weightLL, &b_weightLL); 
  fChain->SetBranchAddress("weightLR", &m_weightLR, &b_weightLR); 
  fChain->SetBranchAddress("weightRL", &m_weightRL, &b_weightRL); 
  fChain->SetBranchAddress("weightRR", &m_weightRR, &b_weightRR); 
  fChain->SetBranchAddress("weight_ee", &m_weight_ee, &b_weight_ee); 
  fChain->SetBranchAddress("weight_emu", &m_weight_emu, &b_weight_emu); 
  fChain->SetBranchAddress("weight_eq", &m_weight_eq, &b_weight_eq); 
  fChain->SetBranchAddress("weight_qq", &m_weight_qq, &b_weight_qq); 
  fChain->SetBranchAddress("xMET", &m_xMET, &b_xMET); 
  fChain->SetBranchAddress("yMET", &m_yMET, &b_yMET); 
  fChain->SetBranchAddress("ycol3", &m_ycol3, &b_ycol3); 
  fChain->SetBranchAddress("ycol356", &m_ycol356, &b_ycol356); 
  fChain->SetBranchAddress("ycol4", &m_ycol4, &b_ycol4); 
  fChain->SetBranchAddress("ycol478", &m_ycol478, &b_ycol478); 
  fChain->SetBranchAddress("ycol5", &m_ycol5, &b_ycol5); 
  fChain->SetBranchAddress("ycol6", &m_ycol6, &b_ycol6); 
  fChain->SetBranchAddress("ycol7", &m_ycol7, &b_ycol7); 
  fChain->SetBranchAddress("ycol8", &m_ycol8, &b_ycol8); 
  fChain->SetBranchAddress("yt", &m_yt, &b_yt); 
  fChain->SetBranchAddress("yt_reco", &m_yt_reco, &b_yt_reco); 
  fChain->SetBranchAddress("ytbar", &m_ytbar, &b_ytbar); 
  fChain->SetBranchAddress("ytt", &m_ytt, &b_ytt); 
} 

