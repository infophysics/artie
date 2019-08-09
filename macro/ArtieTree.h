//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Aug  7 15:53:25 2019 by ROOT version 6.08/00
// from TTree artie/
// found on file: argon.root
//////////////////////////////////////////////////////////

#ifndef ArtieTree_h
#define ArtieTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class ArtieTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Double_t        gen_energy;
   Double_t        arrival_time;
   Double_t        arrival_e;
   Int_t           num_elastic;
   Int_t           num_inelastic;
   Int_t           num_ncapture;
   Int_t           num_fission;
   Int_t           num_scatter;
   Int_t           num_scatout;
   Int_t           z_scatter;
   Int_t           first_gas;
   Double_t        max_dphi;
   Double_t        max_dp;
   Double_t        max_de;

   // List of branches
   TBranch        *b_gen_energy;   //!
   TBranch        *b_arrival_time;   //!
   TBranch        *b_arrival_e;   //!
   TBranch        *b_num_elastic;   //!
   TBranch        *b_num_inelastic;   //!
   TBranch        *b_num_ncapture;   //!
   TBranch        *b_num_fission;   //!
   TBranch        *b_num_scatter;   //!
   TBranch        *b_num_scatout;   //!
   TBranch        *b_z_scatter;   //!
   TBranch        *b_first_gas;   //!
   TBranch        *b_max_dphi;   //!
   TBranch        *b_max_dp;   //!
   TBranch        *b_max_de;   //!

   //ArtieTree(TTree *tree=0);
   ArtieTree(TString file);
   virtual ~ArtieTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

//#ifdef ArtieTree_cxx
ArtieTree::ArtieTree(TString file) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   TTree* tree=0;
   TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(file);
   if (!f || !f->IsOpen()) {
     f = new TFile(file);
     f->GetObject("artie",tree);
   }
   Init(tree);
}

ArtieTree::~ArtieTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t ArtieTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t ArtieTree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void ArtieTree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("gen_energy", &gen_energy, &b_gen_energy);
   fChain->SetBranchAddress("arrival_time", &arrival_time, &b_arrival_time);
   fChain->SetBranchAddress("arrival_e", &arrival_e, &b_arrival_e);
   fChain->SetBranchAddress("num_elastic", &num_elastic, &b_num_elastic);
   fChain->SetBranchAddress("num_inelastic", &num_inelastic, &b_num_elastic);
   fChain->SetBranchAddress("num_ncapture", &num_ncapture, &b_num_elastic);
   fChain->SetBranchAddress("num_fission", &num_fission, &b_num_elastic);
   fChain->SetBranchAddress("num_scatter", &num_scatter, &b_num_scatter);
   fChain->SetBranchAddress("num_scatout", &num_scatout, &b_num_scatout);
   fChain->SetBranchAddress("z_scatter", &z_scatter, &b_z_scatter);
   fChain->SetBranchAddress("first_gas", &first_gas, &b_first_gas);
   fChain->SetBranchAddress("max_dphi", &max_dphi, &b_max_dphi);
   fChain->SetBranchAddress("max_dp", &max_dp, &b_max_dp);
   fChain->SetBranchAddress("max_de", &max_de, &b_max_de);
   Notify();
}

Bool_t ArtieTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void ArtieTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t ArtieTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
//#endif // #ifdef ArtieTree_cxx
