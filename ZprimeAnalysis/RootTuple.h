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
//  *   File Generated on Tue Apr 28 11:35:59 2015                            * 
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
    inline Double_t  Delta_phi_l() const {b_Delta_phi_l->GetEntry(m_currentEvent);return m_Delta_phi_l;} 
    inline Double_t  Delta_y() const {b_Delta_y->GetEntry(m_currentEvent);return m_Delta_y;} 
    inline Double_t  Delta_y_reco() const {b_Delta_y_reco->GetEntry(m_currentEvent);return m_Delta_y_reco;} 
    inline vector<double>*  E() const {b_E->GetEntry(m_currentEvent);return m_E;} 
    inline Double_t  E6_reco() const {b_E6_reco->GetEntry(m_currentEvent);return m_E6_reco;} 
    inline Double_t  Et() const {b_Et->GetEntry(m_currentEvent);return m_Et;} 
    inline Double_t  Etmiss() const {b_Etmiss->GetEntry(m_currentEvent);return m_Etmiss;} 
    inline Double_t  Etmissx() const {b_Etmissx->GetEntry(m_currentEvent);return m_Etmissx;} 
    inline Double_t  Etmissy() const {b_Etmissy->GetEntry(m_currentEvent);return m_Etmissy;} 
    inline Double_t  HT() const {b_HT->GetEntry(m_currentEvent);return m_HT;} 
    inline Double_t  MCT1() const {b_MCT1->GetEntry(m_currentEvent);return m_MCT1;} 
    inline Double_t  MCT2() const {b_MCT2->GetEntry(m_currentEvent);return m_MCT2;} 
    inline Double_t  MCT3() const {b_MCT3->GetEntry(m_currentEvent);return m_MCT3;} 
    inline Double_t  MT1() const {b_MT1->GetEntry(m_currentEvent);return m_MT1;} 
    inline Double_t  MT2() const {b_MT2->GetEntry(m_currentEvent);return m_MT2;} 
    inline Double_t  MT3() const {b_MT3->GetEntry(m_currentEvent);return m_MT3;} 
    inline Double_t  MlCT() const {b_MlCT->GetEntry(m_currentEvent);return m_MlCT;} 
    inline Double_t  MlT() const {b_MlT->GetEntry(m_currentEvent);return m_MlT;} 
    inline Double_t  Mll() const {b_Mll->GetEntry(m_currentEvent);return m_Mll;} 
    inline Double_t  Mtt() const {b_Mtt->GetEntry(m_currentEvent);return m_Mtt;} 
    inline Double_t  Mtt_reco() const {b_Mtt_reco->GetEntry(m_currentEvent);return m_Mtt_reco;} 
    inline Double_t  Mttvis() const {b_Mttvis->GetEntry(m_currentEvent);return m_Mttvis;} 
    inline vector<double>*  Px() const {b_Px->GetEntry(m_currentEvent);return m_Px;} 
    inline vector<double>*  Py() const {b_Py->GetEntry(m_currentEvent);return m_Py;} 
    inline vector<double>*  Pz() const {b_Pz->GetEntry(m_currentEvent);return m_Pz;} 
    inline Double_t  Pz6_reco() const {b_Pz6_reco->GetEntry(m_currentEvent);return m_Pz6_reco;} 
    inline vector<int>*  barcode() const {b_barcode->GetEntry(m_currentEvent);return m_barcode;} 
    inline Double_t  beta() const {b_beta->GetEntry(m_currentEvent);return m_beta;} 
    inline Double_t  cosphil() const {b_cosphil->GetEntry(m_currentEvent);return m_cosphil;} 
    inline Double_t  costheta5() const {b_costheta5->GetEntry(m_currentEvent);return m_costheta5;} 
    inline Double_t  costheta5cm() const {b_costheta5cm->GetEntry(m_currentEvent);return m_costheta5cm;} 
    inline Double_t  costheta5col() const {b_costheta5col->GetEntry(m_currentEvent);return m_costheta5col;} 
    inline Double_t  costheta7() const {b_costheta7->GetEntry(m_currentEvent);return m_costheta7;} 
    inline Double_t  costheta7cm() const {b_costheta7cm->GetEntry(m_currentEvent);return m_costheta7cm;} 
    inline Double_t  costheta7col() const {b_costheta7col->GetEntry(m_currentEvent);return m_costheta7col;} 
    inline Double_t  costhetastar() const {b_costhetastar->GetEntry(m_currentEvent);return m_costhetastar;} 
    inline Double_t  costhetastar_reco() const {b_costhetastar_reco->GetEntry(m_currentEvent);return m_costhetastar_reco;} 
    inline Double_t  costhetat() const {b_costhetat->GetEntry(m_currentEvent);return m_costhetat;} 
    inline Double_t  costhetat_reco() const {b_costhetat_reco->GetEntry(m_currentEvent);return m_costhetat_reco;} 
    inline Double_t  costhetatcol() const {b_costhetatcol->GetEntry(m_currentEvent);return m_costhetatcol;} 
    inline Double_t  ct7ct5() const {b_ct7ct5->GetEntry(m_currentEvent);return m_ct7ct5;} 
    inline Double_t  ct7ct5cm() const {b_ct7ct5cm->GetEntry(m_currentEvent);return m_ct7ct5cm;} 
    inline Double_t  ct7ct5col() const {b_ct7ct5col->GetEntry(m_currentEvent);return m_ct7ct5col;} 
    inline Double_t  eta3() const {b_eta3->GetEntry(m_currentEvent);return m_eta3;} 
    inline Double_t  eta356() const {b_eta356->GetEntry(m_currentEvent);return m_eta356;} 
    inline Double_t  eta4() const {b_eta4->GetEntry(m_currentEvent);return m_eta4;} 
    inline Double_t  eta478() const {b_eta478->GetEntry(m_currentEvent);return m_eta478;} 
    inline Double_t  eta5() const {b_eta5->GetEntry(m_currentEvent);return m_eta5;} 
    inline Double_t  eta6() const {b_eta6->GetEntry(m_currentEvent);return m_eta6;} 
    inline Double_t  eta7() const {b_eta7->GetEntry(m_currentEvent);return m_eta7;} 
    inline Double_t  eta8() const {b_eta8->GetEntry(m_currentEvent);return m_eta8;} 
    inline Double_t  fl() const {b_fl->GetEntry(m_currentEvent);return m_fl;} 
    inline Double_t  m478() const {b_m478->GetEntry(m_currentEvent);return m_m478;} 
    inline Double_t  mt_reco() const {b_mt_reco->GetEntry(m_currentEvent);return m_mt_reco;} 
    inline Double_t  pT3() const {b_pT3->GetEntry(m_currentEvent);return m_pT3;} 
    inline Double_t  pT356() const {b_pT356->GetEntry(m_currentEvent);return m_pT356;} 
    inline Double_t  pT4() const {b_pT4->GetEntry(m_currentEvent);return m_pT4;} 
    inline Double_t  pT478() const {b_pT478->GetEntry(m_currentEvent);return m_pT478;} 
    inline Double_t  pT5() const {b_pT5->GetEntry(m_currentEvent);return m_pT5;} 
    inline Double_t  pT6() const {b_pT6->GetEntry(m_currentEvent);return m_pT6;} 
    inline Double_t  pT7() const {b_pT7->GetEntry(m_currentEvent);return m_pT7;} 
    inline Double_t  pT8() const {b_pT8->GetEntry(m_currentEvent);return m_pT8;} 
    inline Double_t  phi3() const {b_phi3->GetEntry(m_currentEvent);return m_phi3;} 
    inline Double_t  phi356() const {b_phi356->GetEntry(m_currentEvent);return m_phi356;} 
    inline Double_t  phi4() const {b_phi4->GetEntry(m_currentEvent);return m_phi4;} 
    inline Double_t  phi478() const {b_phi478->GetEntry(m_currentEvent);return m_phi478;} 
    inline Double_t  phi5() const {b_phi5->GetEntry(m_currentEvent);return m_phi5;} 
    inline Double_t  phi6() const {b_phi6->GetEntry(m_currentEvent);return m_phi6;} 
    inline Double_t  phi7() const {b_phi7->GetEntry(m_currentEvent);return m_phi7;} 
    inline Double_t  phi8() const {b_phi8->GetEntry(m_currentEvent);return m_phi8;} 
    inline Double_t  weight() const {b_weight->GetEntry(m_currentEvent);return m_weight;} 
    inline Double_t  weightLL() const {b_weightLL->GetEntry(m_currentEvent);return m_weightLL;} 
    inline Double_t  weightLR() const {b_weightLR->GetEntry(m_currentEvent);return m_weightLR;} 
    inline Double_t  weightRL() const {b_weightRL->GetEntry(m_currentEvent);return m_weightRL;} 
    inline Double_t  weightRR() const {b_weightRR->GetEntry(m_currentEvent);return m_weightRR;} 
    inline Double_t  weight_ee() const {b_weight_ee->GetEntry(m_currentEvent);return m_weight_ee;} 
    inline Double_t  weight_emu() const {b_weight_emu->GetEntry(m_currentEvent);return m_weight_emu;} 
    inline Double_t  weight_eq() const {b_weight_eq->GetEntry(m_currentEvent);return m_weight_eq;} 
    inline Double_t  weight_qq() const {b_weight_qq->GetEntry(m_currentEvent);return m_weight_qq;} 
    inline Double_t  ycol3() const {b_ycol3->GetEntry(m_currentEvent);return m_ycol3;} 
    inline Double_t  ycol356() const {b_ycol356->GetEntry(m_currentEvent);return m_ycol356;} 
    inline Double_t  ycol4() const {b_ycol4->GetEntry(m_currentEvent);return m_ycol4;} 
    inline Double_t  ycol478() const {b_ycol478->GetEntry(m_currentEvent);return m_ycol478;} 
    inline Double_t  ycol5() const {b_ycol5->GetEntry(m_currentEvent);return m_ycol5;} 
    inline Double_t  ycol6() const {b_ycol6->GetEntry(m_currentEvent);return m_ycol6;} 
    inline Double_t  ycol7() const {b_ycol7->GetEntry(m_currentEvent);return m_ycol7;} 
    inline Double_t  ycol8() const {b_ycol8->GetEntry(m_currentEvent);return m_ycol8;} 
    inline Double_t  yt() const {b_yt->GetEntry(m_currentEvent);return m_yt;} 
    inline Double_t  yt_reco() const {b_yt_reco->GetEntry(m_currentEvent);return m_yt_reco;} 
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

    Double_t  m_Delta_phi_l; 
    Double_t  m_Delta_y; 
    Double_t  m_Delta_y_reco; 
    vector<double>*  m_E; 
    Double_t  m_E6_reco; 
    Double_t  m_Et; 
    Double_t  m_Etmiss; 
    Double_t  m_Etmissx; 
    Double_t  m_Etmissy; 
    Double_t  m_HT; 
    Double_t  m_MCT1; 
    Double_t  m_MCT2; 
    Double_t  m_MCT3; 
    Double_t  m_MT1; 
    Double_t  m_MT2; 
    Double_t  m_MT3; 
    Double_t  m_MlCT; 
    Double_t  m_MlT; 
    Double_t  m_Mll; 
    Double_t  m_Mtt; 
    Double_t  m_Mtt_reco; 
    Double_t  m_Mttvis; 
    vector<double>*  m_Px; 
    vector<double>*  m_Py; 
    vector<double>*  m_Pz; 
    Double_t  m_Pz6_reco; 
    vector<int>*  m_barcode; 
    Double_t  m_beta; 
    Double_t  m_cosphil; 
    Double_t  m_costheta5; 
    Double_t  m_costheta5cm; 
    Double_t  m_costheta5col; 
    Double_t  m_costheta7; 
    Double_t  m_costheta7cm; 
    Double_t  m_costheta7col; 
    Double_t  m_costhetastar; 
    Double_t  m_costhetastar_reco; 
    Double_t  m_costhetat; 
    Double_t  m_costhetat_reco; 
    Double_t  m_costhetatcol; 
    Double_t  m_ct7ct5; 
    Double_t  m_ct7ct5cm; 
    Double_t  m_ct7ct5col; 
    Double_t  m_eta3; 
    Double_t  m_eta356; 
    Double_t  m_eta4; 
    Double_t  m_eta478; 
    Double_t  m_eta5; 
    Double_t  m_eta6; 
    Double_t  m_eta7; 
    Double_t  m_eta8; 
    Double_t  m_fl; 
    Double_t  m_m478; 
    Double_t  m_mt_reco; 
    Double_t  m_pT3; 
    Double_t  m_pT356; 
    Double_t  m_pT4; 
    Double_t  m_pT478; 
    Double_t  m_pT5; 
    Double_t  m_pT6; 
    Double_t  m_pT7; 
    Double_t  m_pT8; 
    Double_t  m_phi3; 
    Double_t  m_phi356; 
    Double_t  m_phi4; 
    Double_t  m_phi478; 
    Double_t  m_phi5; 
    Double_t  m_phi6; 
    Double_t  m_phi7; 
    Double_t  m_phi8; 
    Double_t  m_weight; 
    Double_t  m_weightLL; 
    Double_t  m_weightLR; 
    Double_t  m_weightRL; 
    Double_t  m_weightRR; 
    Double_t  m_weight_ee; 
    Double_t  m_weight_emu; 
    Double_t  m_weight_eq; 
    Double_t  m_weight_qq; 
    Double_t  m_ycol3; 
    Double_t  m_ycol356; 
    Double_t  m_ycol4; 
    Double_t  m_ycol478; 
    Double_t  m_ycol5; 
    Double_t  m_ycol6; 
    Double_t  m_ycol7; 
    Double_t  m_ycol8; 
    Double_t  m_yt; 
    Double_t  m_yt_reco; 
    Double_t  m_ytbar; 
    Double_t  m_ytt; 

    TBranch*  b_Delta_phi_l; 
    TBranch*  b_Delta_y; 
    TBranch*  b_Delta_y_reco; 
    TBranch*  b_E; 
    TBranch*  b_E6_reco; 
    TBranch*  b_Et; 
    TBranch*  b_Etmiss; 
    TBranch*  b_Etmissx; 
    TBranch*  b_Etmissy; 
    TBranch*  b_HT; 
    TBranch*  b_MCT1; 
    TBranch*  b_MCT2; 
    TBranch*  b_MCT3; 
    TBranch*  b_MT1; 
    TBranch*  b_MT2; 
    TBranch*  b_MT3; 
    TBranch*  b_MlCT; 
    TBranch*  b_MlT; 
    TBranch*  b_Mll; 
    TBranch*  b_Mtt; 
    TBranch*  b_Mtt_reco; 
    TBranch*  b_Mttvis; 
    TBranch*  b_Px; 
    TBranch*  b_Py; 
    TBranch*  b_Pz; 
    TBranch*  b_Pz6_reco; 
    TBranch*  b_barcode; 
    TBranch*  b_beta; 
    TBranch*  b_cosphil; 
    TBranch*  b_costheta5; 
    TBranch*  b_costheta5cm; 
    TBranch*  b_costheta5col; 
    TBranch*  b_costheta7; 
    TBranch*  b_costheta7cm; 
    TBranch*  b_costheta7col; 
    TBranch*  b_costhetastar; 
    TBranch*  b_costhetastar_reco; 
    TBranch*  b_costhetat; 
    TBranch*  b_costhetat_reco; 
    TBranch*  b_costhetatcol; 
    TBranch*  b_ct7ct5; 
    TBranch*  b_ct7ct5cm; 
    TBranch*  b_ct7ct5col; 
    TBranch*  b_eta3; 
    TBranch*  b_eta356; 
    TBranch*  b_eta4; 
    TBranch*  b_eta478; 
    TBranch*  b_eta5; 
    TBranch*  b_eta6; 
    TBranch*  b_eta7; 
    TBranch*  b_eta8; 
    TBranch*  b_fl; 
    TBranch*  b_m478; 
    TBranch*  b_mt_reco; 
    TBranch*  b_pT3; 
    TBranch*  b_pT356; 
    TBranch*  b_pT4; 
    TBranch*  b_pT478; 
    TBranch*  b_pT5; 
    TBranch*  b_pT6; 
    TBranch*  b_pT7; 
    TBranch*  b_pT8; 
    TBranch*  b_phi3; 
    TBranch*  b_phi356; 
    TBranch*  b_phi4; 
    TBranch*  b_phi478; 
    TBranch*  b_phi5; 
    TBranch*  b_phi6; 
    TBranch*  b_phi7; 
    TBranch*  b_phi8; 
    TBranch*  b_weight; 
    TBranch*  b_weightLL; 
    TBranch*  b_weightLR; 
    TBranch*  b_weightRL; 
    TBranch*  b_weightRR; 
    TBranch*  b_weight_ee; 
    TBranch*  b_weight_emu; 
    TBranch*  b_weight_eq; 
    TBranch*  b_weight_qq; 
    TBranch*  b_ycol3; 
    TBranch*  b_ycol356; 
    TBranch*  b_ycol4; 
    TBranch*  b_ycol478; 
    TBranch*  b_ycol5; 
    TBranch*  b_ycol6; 
    TBranch*  b_ycol7; 
    TBranch*  b_ycol8; 
    TBranch*  b_yt; 
    TBranch*  b_yt_reco; 
    TBranch*  b_ytbar; 
    TBranch*  b_ytt; 
}; 
#endif 

