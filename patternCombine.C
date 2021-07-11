#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <fstream>
#include "AnalyzeMps.h"
//root -l patternCombine.C'(4043,240)'
//HelFreq = Helicity frequency = 2*LED flashing frequency
void patternCombine(Int_t run,Int_t HelFreq,float npat=2.0){

   TString root_file = Form("~/scratch/rootfiles/r%d.root",run);
   TChain* tree = new TChain("evt");
   tree->Add(root_file);
   if(root_file==NULL){
   cout<<"Rootfile doesn't exit in "<<root_file<<endl;
   }
   if(tree==NULL){
   cout<<"Tree "<<tree<<" doesn't exit"<<endl;
   }

   tree->SetMakeClass(1);
  
   tree->SetBranchAddress("CodaEventNumber", &CodaEventNumber, &b_CodaEventNumber);
   tree->SetBranchAddress("CodaEventType", &CodaEventType, &b_CodaEventType);
   tree->SetBranchAddress("Coda_CleanData", &Coda_CleanData, &b_Coda_CleanData);
   tree->SetBranchAddress("Coda_ScanData1", &Coda_ScanData1, &b_Coda_ScanData1);
   tree->SetBranchAddress("Coda_ScanData2", &Coda_ScanData2, &b_Coda_ScanData2);
   tree->SetBranchAddress("qwk_mod0ch0", &qwk_mod0ch0_hw_sum, &b_qwk_mod0ch0);
   tree->SetBranchAddress("qwk_mod0ch1", &qwk_mod0ch1_hw_sum, &b_qwk_mod0ch1);
   tree->SetBranchAddress("qwk_mod0ch2", &qwk_mod0ch2_hw_sum, &b_qwk_mod0ch2);
   tree->SetBranchAddress("qwk_mod0ch3", &qwk_mod0ch3_hw_sum, &b_qwk_mod0ch3);
   tree->SetBranchAddress("qwk_mod0ch4", &qwk_mod0ch4_hw_sum, &b_qwk_mod0ch4);
   tree->SetBranchAddress("qwk_mod0ch5", &qwk_mod0ch5_hw_sum, &b_qwk_mod0ch5);
   tree->SetBranchAddress("qwk_mod0ch6", &qwk_mod0ch6_hw_sum, &b_qwk_mod0ch6);
   tree->SetBranchAddress("qwk_mod0ch7", &qwk_mod0ch7_hw_sum, &b_qwk_mod0ch7);
   tree->SetBranchAddress("ErrorFlag", &ErrorFlag, &b_ErrorFlag);
   Long64_t nentries = tree->GetEntries();
   cout<<"Total Entries: \t"<<nentries<<endl;
   gStyle->SetOptStat(1111);

   TFile f(Form("./rootfiles/patternCombine_%d_run%d.root",int(npat),run),"recreate");
   TTree patternSum("patternSum","patternSum");
   TBranch* NumVal;
   TBranch* SumVal;
   Double_t octet[int(npat)];
   vector<Double_t> patsum[int(npat)];
   Int_t count = 0;;
   vector<Int_t> octetNo[int(npat)];

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     tree->GetEntry(jentry);
    for(int ipat=0;ipat<int(npat);ipat++){
     if((CodaEventNumber-ipat)/npat - int((CodaEventNumber-ipat)/int(npat)) == 0) {
       octet[ipat] = qwk_mod0ch4_hw_sum_raw;
     if(octet[ipat]!=0){
       count++;
       octetNo[ipat].push_back(count);
       patsum[ipat].push_back(octet[ipat]);
     }
     }
     }
     }
     Double_t sum0=0;
     Double_t sum1=0;
     Int_t pattern_number;
     SumVal = patternSum.Branch("pattern_number",&pattern_number);
     SumVal = patternSum.Branch("pattern_sum_0",&sum0);
     SumVal = patternSum.Branch("pattern_sum_1",&sum1);

    for(int i=0;i<nentries/int(npat);i++){
        pattern_number = i;
      for(int ipat=0;ipat<npat;ipat++){
       if(ipat/2.0-int(ipat/2)==0)
	sum0 = +patsum[ipat][i]/npat*2;
       else
	sum1 = +patsum[ipat][i]/npat*2;
      }
        cout<<pattern_number<<endl;
        patternSum.Fill();
    }
     f.Write();
     f.Close();
   }
