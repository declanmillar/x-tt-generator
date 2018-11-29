#include "RootTuple.h"


RootTuple::RootTuple() :
	m_filename("output.root"),
	m_treename("tree")
{
}


RootTuple::RootTuple(std::string filename, std::string treename, std::string treename2) :
	m_filename(filename),
	m_treename(treename),
	m_treename2(treename2)
{
	if (treename.empty()){
		std::cout << "RootTuple:: Warning: Using default tree name" << std::endl;
		m_treename = "tree";
	}
	std::string fileExtension = ".root";
	if (filename.length() <= fileExtension.length()){
		std::cout << "RootTuple:: Warning: Using default file name" << std::endl;
		m_filename = "output.root";
	}
	if (filename.compare(filename.length() - fileExtension.length(), fileExtension.length(), fileExtension) != 0){
		m_filename = filename + ".root";
	}
}


RootTuple::~RootTuple()
{
}


void RootTuple::Initialise()
{
	m_file = new TFile(m_filename.c_str(), "RECREATE");
	if (!m_file){
		std::cout << "RootTuple:: Error: Cannot create ROOT file" << std::endl;
		return;
	}

	m_tree2 = new TTree(m_treename2.c_str(), m_treename2.c_str());
	if (!m_tree2) std::cout << "RootTuple:: Error: Cannot create process tree" << std::endl;

	m_tree = new TTree(m_treename.c_str(), m_treename.c_str());
	if (!m_tree) std::cout << "RootTuple:: Error: Cannot create events tree" << std::endl;

	DeclareBranches();
}


void RootTuple::AddEvent()
{
	if ((int)m_Px.size() != (int)m_barcode.size() ||
		(int)m_Py.size() != (int)m_barcode.size() ||
		(int)m_Pz.size() != (int)m_barcode.size() ||
		(int)m_E.size()  != (int)m_barcode.size())
		std::cout << "RootTuple:: Warning: Inconsistent vector sizes" << std::endl;

	FillBranches();
	ClearVectors();
}


void RootTuple::Write()
{
	FillProcessBranches();
	if (m_file) m_file->Write();
	else std::cout << "RootTuple:: Error: No ROOT file was opened" << std::endl;
}


void RootTuple::Close()
{
	FillProcessBranches();
	if (m_file){
		m_file->Write();
		m_file->Close();
		delete m_file;
	}
	else std::cout << "RootTuple:: Error: No ROOT file was opened" << std::endl;
}


void RootTuple::AddSingleDouble(std::string branchname, double *ptr)
{
	if (m_tree2->GetBranch(branchname.c_str())) m_tree2->GetBranch(branchname.c_str())->SetAddress(ptr);
	else m_tree2->Branch(branchname.c_str(), ptr);
}


void RootTuple::AddSingleFloat(std::string branchname, float *ptr)
{
	if (m_tree2->GetBranch(branchname.c_str())) m_tree2->GetBranch(branchname.c_str())->SetAddress(ptr);
	else m_tree2->Branch(branchname.c_str(), ptr);
}


void RootTuple::AddSingleInt(std::string branchname, int *ptr)
{
	if (m_tree2->GetBranch(branchname.c_str())) m_tree2->GetBranch(branchname.c_str())->SetAddress(ptr);
	else m_tree2->Branch(branchname.c_str(), ptr);
}


void RootTuple::AddSingleBool(std::string branchname, bool *ptr)
{
	if (m_tree2->GetBranch(branchname.c_str())) m_tree2->GetBranch(branchname.c_str())->SetAddress(ptr);
	else m_tree2->Branch(branchname.c_str(), ptr);
}


void RootTuple::AddParticle(int barcode, double px, double py, double pz, double e)
{
	m_barcode.push_back(barcode);
	m_Px.push_back(px);
	m_Py.push_back(py);
	m_Pz.push_back(pz);
	m_E.push_back(e);
}


void RootTuple::SetWeight(double weight)
{
	m_weight = weight;
}


void RootTuple::SetDoubleBranch(std::string branchname, double *ptr)
{
	if (m_tree->GetBranch(branchname.c_str())) m_tree->GetBranch(branchname.c_str())->SetAddress(ptr);
	else m_tree->Branch(branchname.c_str(), ptr);
}


void RootTuple::SetFloatBranch(std::string branchname, float *ptr)
{
	if (m_tree->GetBranch(branchname.c_str())) m_tree->GetBranch(branchname.c_str())->SetAddress(ptr);
	else m_tree->Branch(branchname.c_str(), ptr);
}


void RootTuple::SetIntBranch(std::string branchname, int *ptr)
{
	if (m_tree->GetBranch(branchname.c_str())) m_tree->GetBranch(branchname.c_str())->SetAddress(ptr);
	else m_tree->Branch(branchname.c_str(), ptr);
}


void RootTuple::SetBoolBranch(std::string branchname, bool *ptr)
{
	if (m_tree->GetBranch(branchname.c_str())) m_tree->GetBranch(branchname.c_str())->SetAddress(ptr);
	else m_tree->Branch(branchname.c_str(), ptr);
}


void RootTuple::DeclareBranches()
{
	m_tree->Branch("weight", &m_weight);
	m_tree->Branch("barcode", &m_barcode);
	m_tree->Branch("Px", &m_Px);
	m_tree->Branch("Py", &m_Py);
	m_tree->Branch("Pz", &m_Pz);
	m_tree->Branch("E", &m_E);
}


void RootTuple::FillBranches()
{
	if (m_tree) m_tree->Fill();
}


void RootTuple::FillProcessBranches()
{
	if (m_tree2) m_tree2->Fill();
}


void RootTuple::ClearVectors()
{
	m_barcode.clear();
	m_Px.clear();
	m_Py.clear();
	m_Pz.clear();
	m_E.clear();
}
