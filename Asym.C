#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <fstream>
#include "AnalyzeMps.h"
//root -l Asym.C'(3899,5407,10.0,745,0.15,0.671,120,0.005,0.055,1550000,1555000,0.005,0.055,0.010,0.045,0.010,0.045,-0.002,0.002,0.020,0.040)'
void Asym(Int_t run,Int_t pmt,Double_t preAmp,Int_t HV,Double_t LL,Double_t Iout,Int_t FF,Double_t asym_X1,Double_t asym_X2,Int_t chsum_min,Int_t chsum,Double_t asym1_X1,Double_t asym1_X2,Double_t asym8_X1,Double_t asym8_X2,Double_t asymF_X1,Double_t asymF_X2,Double_t asymP_X1,Double_t asymP_X2,Double_t asym_Y1,Double_t asym_Y2){

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
     
     double max_adc_axis = 0;
     if(FF==120)
        max_adc_axis = 262000000;
     else if(FF==240)
        max_adc_axis = 131000000;
     else if(FF==480)
        max_adc_axis = 65500000;
     else if(FF==960)
        max_adc_axis = 32750000;
     TH2F *h1_raw = new TH2F("h1_raw", "", nentries/100, 0, nentries, 1000, 0,max_adc_axis);
     TH2F *h1_raw_even = new TH2F("h1_raw_even", "", nentries/100, 0, nentries, 1000, 0, max_adc_axis);
     TH2F *h1_raw_odd = new TH2F("h1_raw_odd", "", nentries/100, 0, nentries, 1000, 0, max_adc_axis);
     TH1F *h1_ped = new TH1F("h1_ped","",100,chsum_min,chsum);;
     TH1F *h1_ped_even = new TH1F("h1_ped_even","",100,chsum_min,chsum);;
     TH1F *h1_ped_odd = new TH1F("h1_ped_odd","",100,chsum_min,chsum);;
     TH2F *h1_rawCor = new TH2F("h1_rawCor", "", nentries/10, 0, nentries, 1000, 0, max_adc_axis);
     TH2F *h1_rawCor_even = new TH2F("h1_rawCor_even", "", nentries/10, 0, nentries, 1000, 0, max_adc_axis);
     TH2F *h1_rawCor_odd = new TH2F("h1_rawCor_odd", "", nentries/10, 0, nentries, 1000, 0, max_adc_axis);
     TH2F *h1_asym_my = new TH2F("h1_asym_my", "", nentries/20, 0, nentries/2, 800, asym_X1, asym_X2);


   for (Long64_t jentry=10; jentry<nentries;jentry++) {
	tree->GetEntry(jentry);
	h1_raw->Fill(CodaEventNumber,qwk_mod0ch4_hw_sum_raw);

	if(CodaEventNumber/2.0 - int(CodaEventNumber/2)==0)
		h1_raw_even->Fill(CodaEventNumber,qwk_mod0ch4_hw_sum_raw);
	else
		h1_raw_odd->Fill(CodaEventNumber,qwk_mod0ch4_hw_sum_raw);

	if(qwk_mod0ch4_hw_sum_raw<chsum){
		h1_ped->Fill(qwk_mod0ch4_hw_sum_raw);
		if(CodaEventNumber/2.0 - int(CodaEventNumber/2)==0)
			h1_ped_even->Fill(qwk_mod0ch4_hw_sum_raw);
		else
			h1_ped_odd->Fill(qwk_mod0ch4_hw_sum_raw);
	}
     }
   gSystem->Exec(Form("rm -rf asymmetry.csv"));
   ofstream outfile("asymmetry.txt");
   TFile f("./rootfiles/asymmetry.root","recreate");
   TTree asym("asym","asym");
   TBranch* pairNum = asym.Branch("pairNum",&pairNumber);
   TBranch* FilterEvtNum = asym.Branch("Filter_Event_Num",&Filter_Event_Num);
   TBranch* AsymVal = asym.Branch("Filter_Asymmetry",&Filter_Asymmetry);
   TBranch* SumVal = asym.Branch("pairSum",&pairSum);

   Double_t pedAve = h1_ped->GetMean();
   Double_t pair[2];
   pair[0] = 0;
   pair[1] = 0;
   for (Long64_t jentry=10; jentry<nentries;jentry++) {
     tree->GetEntry(jentry);
     if(CodaEventNumber/2.0 - int(CodaEventNumber/2) == 0) {
       h1_rawCor_even->Fill(CodaEventNumber, qwk_mod0ch4_hw_sum_raw-pedAve);
       pair[1] = qwk_mod0ch4_hw_sum_raw-pedAve;
     } else {
       h1_rawCor_odd->Fill(CodaEventNumber, qwk_mod0ch4_hw_sum_raw-pedAve);
       pair[0] = qwk_mod0ch4_hw_sum_raw-pedAve;
     }
       h1_rawCor->Fill(CodaEventNumber, qwk_mod0ch4_hw_sum_raw-pedAve);
     if(pair[0] != 0&&pair[1] != 0) {
       pairNumber++;
	  sum = pair[1] + pair[0];
          pairSum = sum;
	  diff = pair[1] - pair[0];
	  Asymmetry = fabs(diff/sum);
	  Asymmetry_Ped = (pair[1] - pair[0])/(pair[1] + pair[0] +2*pedAve);
     if(sum > 1.0E6){
      h1_asym_my->Fill(pairNumber, Asymmetry);
	  Filter_Event_Num = pairNumber;
	  Filter_Asymmetry = Asymmetry;
     }else{
	  Filter_Event_Num = 0.0;
	  Filter_Asymmetry = 0.0;
     }
     outfile << pairNumber <<"\t"<< Filter_Event_Num << "\t"<< Filter_Asymmetry << "\t"<<pairSum <<endl;
//     pairNum->Fill();
//     FilterEvtNum->Fill();
//     AsymVal->Fill();
//     SumVal->Fill();
      asym.Fill();
      pair[0] = 0;
      pair[1] = 0;
     }
     }
     outfile.close();
     asym.SetEntries();
     f.Write();
     f.Close();

     h1_raw->SetTitle(Form("Filter wheel scan: pmt#%4d @ %5d V, I_{cath} = %2.2f nA;Event No.;ADC (raw ch sum)", pmt, HV, LL));
     h1_ped->SetTitle("Shutter closed (pedestal data);ADC (raw ch sum);Events");
     h1_rawCor->SetTitle("Pedestal-corrected;Event No.;ADC (raw, ped-corr ch sum)");
     h1_asym_my->SetTitle("Pair-wise Asymmetry vs. pair No.;Pair No.;A_{LED}");
     h1_raw_even->SetMarkerColor(kBlue);
     h1_raw_odd->SetMarkerColor(kRed);
     h1_ped->SetLineColor(kBlack);
     h1_ped_even->SetLineColor(kBlue);
     h1_ped_odd->SetLineColor(kRed);
     h1_rawCor_even->SetMarkerColor(kBlue);
     h1_rawCor_odd->SetMarkerColor(kRed);
 
    TCanvas *c1 = new TCanvas("c1","Mps hwsum Tree Distributions",1000,650);
    c1->Clear();
    c1->Divide(2,2); 
    c1->cd(1);
    gPad->SetGridx(1);
    gPad->SetGridy(1);
    h1_raw->Draw();
    h1_raw_odd->Draw("same");
    h1_raw_even->Draw("same");
    TPaveLabel *pt_rl = new TPaveLabel(0.25,0.85,0.58,0.89,Form("run %4d, PreAmp @ %3.2f M#Omega LED flash @ %4d Hz:", run, preAmp,FF),"NDC");
    pt_rl->SetBorderSize(0);
    pt_rl->SetTextColor(kBlack);
    pt_rl->SetTextSize(0.80);
    pt_rl->SetFillColor(0);
    pt_rl->Draw();
    c1->cd(2);
    gPad->SetGridx(1);
    gPad->SetGridy(1);
    h1_ped->Draw();
    h1_ped_odd->Draw("same");
    h1_ped_even->Draw("same");
    c1->cd(3);
    gPad->SetGridx(1);
    gPad->SetGridy(1);
    h1_rawCor->Draw();
    h1_rawCor_odd->Draw("same");
    h1_rawCor_even->Draw("same");
    c1->cd(4);
    gPad->SetGridx(1);
    gPad->SetGridy(1);
    h1_asym_my->Draw();
    c1->SaveAs(Form("plots/r%dplotsMps.png",run));
   }
