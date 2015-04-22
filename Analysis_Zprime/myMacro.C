#include "TChain.h"
#include <iostream>
#include "Riostream.h"
#include <vector>
#include <TROOT.h>

void myMacro() {
	TChain* t_tuple = new TChain("RootTuple");
	t_tuple->Add("output.root");
	Double_t *mtt;
	t_tuple->SetBranchAddress("Mtt", mtt);
	for(int event = 0; event < t_tuple->GetEntries(); event++){
		t_tuple->GetEntry(event);
		cout << "Mtt of " << event << ": " << mtt << endl;
	}
}

// void myMacro() {
// 	TChain* t_tuple = new TChain("RootTuple");
// 	t_tuple->Add("output.root");
// 	vector<int>* pz = 0;
// 	t_tuple->SetBranchAddress("Pz", &pz);
// 	// t_tuple->GetEntry(100);
// 	// cout << "Pz: " << pz;
// 	for(int event = 0; event < t_tuple->GetEntries(); event++){
// 		t_tuple->GetEntry(event);
// 		cout << "Pz of " << event << ": " << pz->at(0) << endl;
// 	}

// }