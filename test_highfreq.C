#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <fstream>
void test_highfreq(Int_t run,Int_t pmt,Double_t preAmp,Int_t HV,Double_t LL,Double_t 
Iout,Int_t FF,Double_t asym_X1,Double_t asym_X2,Int_t chsum,Double_t asym1_X1,Double_t asym1_X2,Double_t asym8_X1,Double_t asym8_X2,Double_t asymF_X1,Double_t asymF_X2,Double_t asymP_X1,Double_t asymP_X2,Double_t asym_Y1,Double_t asym_Y2){

     TString root_file = Form("/home/daq/scratch/rootfiles/r%d.root",run);
     TChain* tree = new TChain("evt");
     tree->Add(root_file);
     if(root_file==NULL){
     cout<<"Rootfile doesn't exit in "<<root_file<<endl;
     }
     if(tree==NULL){
     cout<<"Tree "<<tree<<" doesn't exit"<<endl;
     }
     Long64_t nentries = tree->GetEntries();
     cout<<"Total Entries: \t"<<nentries<<endl;
     gStyle->SetOptStat(1111);

     double max_adc_axis = 131000000;//for 240Hz flipping
     TH2F *h1_raw = new TH2F("h1_raw", "", nentries/100, 0, nentries, 1000, 0, max_adc_axis);
     TH2F *h1_raw_even = new TH2F("h1_raw_even", "", nentries/100, 0, nentries, 1000, 0, max_adc_axis);
     TH2F *h1_raw_odd = new TH2F("h1_raw_odd", "", nentries/100, 0, nentries, 1000, 0, max_adc_axis);
     TH1F *h1_ped;
     TH1F *h1_ped_even;
     TH1F *h1_ped_odd;
     TH2F *h1_rawCor = new TH2F("h1_rawCor", "", nentries/10, 0, nentries, 1000, 0, max_adc_axis);
     TH2F *h1_rawCor_even = new TH2F("h1_rawCor_even", "", nentries/10, 0, nentries, 1000, 0, max_adc_axis);
     TH2F *h1_rawCor_odd = new TH2F("h1_rawCor_odd", "", nentries/10, 0, nentries, 1000, 0, max_adc_axis);
     TH2F *h1_asym_my = new TH2F("h1_asym_my", "", nentries/20, 0, nentries/2, 800, asym_X1, asym_X2);

     const int positions = 160;
     const int filters = 8;
     const int wheelCycles = 20;
     TH1S *h[positions];
     TH1S *hPed;
     char name[20];
     sprintf(name,"hPed_asym");
     hPed = new TH1S(name,"",100, asymP_X1, asymP_X2);
     TH1S *havg[positions];
     Int_t f;

    for(f=0;f<positions;f++){
        sprintf(name,"h_asym%d",f);
     if((f-2)/8.0 - int((f-2)/8) == 0){
        h[f] = new TH1S(name,"",100, asym1_X1, asym1_X2);
     }else if((f-1)/8.0 - int((f-1)/8) == 0){
        h[f] = new TH1S(name,"",100, asym8_X1, asym8_X2);
     }else{
        h[f] = new TH1S(name,"",100, asymF_X1, asymF_X2);
     }
        h[f]->SetYTitle("Events");
        h[f]->SetXTitle("A_{LED}");
     }
    for(f=0;f<filters;f++){
        sprintf(name,"havg_asym%d",f);
        if(f==0)
        havg[f] = new TH1S(name,"",100, asym1_X1, asym1_X2);
        else
        havg[f] = new TH1S(name,"",100, asymF_X1, asymF_X2);
        havg[f]->SetYTitle("Events");
        havg[f]->SetXTitle("A_{LED}");
     }

     Int_t samples = 200;
     Int_t blocks = 4;
     Double_t Vperch = 76.29E-6; //Volts per channel
     const Int_t arraysize = nentries;
     const Int_t pairsize = arraysize/2;
//     TBranch *b_qwk_mod0ch4;
     TBranch *b_CodaEventNumber;
     Double_t CodaEventNumber;
//     Double_t hw_sum_raw;
     tree->SetMakeClass(1);
//     tree->SetBranchAddress("qwk_mod0ch4",&hw_sum_raw,&b_qwk_mod0ch4);
     tree->SetBranchAddress("CodaEventNumber",&CodaEventNumber,&b_CodaEventNumber);
 
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
	tree->GetEntry(jentry);
//	cout<<jentry<<"\t"<<CodaEventNumber<<"\t"<<hw_sum_raw<<endl;
	cout<<jentry<<"\t"<<CodaEventNumber<<endl;
     }
/*
     TCut ped_cut = Form("CodaEventNumber>10&&qwk_mod0ch4.hw_sum_raw<%d",chsum);
     TCut even_cut = Form("qwk_mod0ch4.hw_sum_raw<%d&&Entry$>10&&(Entry$/2.0-int(Entry$/2))==0",chsum); 
     TCut odd_cut = Form("qwk_mod0ch4.hw_sum_raw<%d&&Entry$>10&&(Entry$/2.0-int(Entry$/2))!=0",chsum); 
     tree->Draw("qwk_mod0ch4.hw_sum_raw>>hped",ped_cut,"goff");
     h1_ped = (TH1F*)gDirectory->FindObject("hped");
     h1_ped->SetDirectory(gROOT);
     tree->Draw("qwk_mod0ch4.hw_sum_raw>>hpedE",even_cut,"gmff");
     h1_ped_even = (TH1F*)gDirectory->FindObject("hpedE");
     h1_ped_even->SetDirectory(gROOT);
     tree->Draw("qwk_mod0ch4.hw_sum_raw>>hpedO",odd_cut,"goff");
     h1_ped_odd = (TH1F*)gDirectory->FindObject("hpedO");
     h1_ped_odd->SetDirectory(gROOT);

     Double_t pedAve = h1_ped->GetMean();
     Long64_t nbytes2 = 0, nb2 = 0;
     Double_t pair[2];
     Long64_t pairNumber = 0;
     Double_t Navg[positions] = {0};
     Double_t NavgPed = 0;
     Double_t Navgmean[filters] = {0};
     Double_t Navgmean_Cal[filters] = {0};
     Double_t dummyError[filters] = {0};
     Double_t Aled[positions] = {0};
     Double_t Aledmean[filters] = {0};
     Double_t AledPed = 0;
     Double_t AledError[filters] = {0};
     Double_t AledErrormean[filters] = {0};
     Double_t sum, diff, Asymmetry, Asymmetry_Ped;
     Double_t pairSum[pairsize];
     Double_t Filter_Event_Num[pairsize];
     Double_t Ped_Event_Num[pairsize];
     int group = 0;
     int group_data_count = 0;
     pair[0] = 0;
     pair[1] = 0;
   for (Long64_t jentry=10; jentry<nentries;jentry++) {
     if(EventNumArray[jentry]/2.0 - int(EventNumArray[jentry]/2) == 0) {
       h1_rawCor_even->Fill(EventNumArray[jentry], SumArray[jentry]-pedAve);
       pair[1] = SumArray[jentry]-pedAve;
     } else {
       h1_rawCor_odd->Fill(EventNumArray[jentry], SumArray[jentry]-pedAve);
       pair[0] = SumArray[jentry]-pedAve;
     }
       h1_rawCor->Fill(EventNumArray[jentry], SumArray[jentry]-pedAve);
     if(pair[0] != 0&&pair[1] != 0) {
       pairNumber++;
	  sum = pair[1] + pair[0];
          pairSum[pairNumber] = sum;
	  diff = pair[1] - pair[0];
	  Asymmetry = fabs(diff/sum);
	  Asymmetry_Ped = (pair[1] - pair[0])/(pair[1] + pair[0] +2*pedAve);
	  outfile << pairNumber << "	" << sum << endl;
     if(sum > 1.0E6){
          h1_asym_my->Fill(pairNumber, Asymmetry);
	  Filter_Event_Num[pairNumber] = pairNumber;
     }else{
	  Ped_Event_Num[pairNumber] = pairNumber; 
     }
     if(Filter_Event_Num[pairNumber]!=Filter_Event_Num[pairNumber-1]+1&&Filter_Event_Num[pairNumber]!=0){
	  group++;
	  cout << group << "	" << pairNumber << "	" << group_data_count << "	" << Asymmetry << "	" << sum << endl;
	  group_data_count = 0;
     }else if(Filter_Event_Num[pairNumber]!=0){
	  group_data_count++;
    
     }
         pair[0] = 0;
         pair[1] = 0;
     }
     }
	  outfile.close();
     h1_raw->SetTitle(Form("Filter wheel scan: pmt#%4d @ %5d V, I_{cath} = %2.2f nA;Event No.;ADC (raw ch sum)", pmt, HV, LL));
     h1_ped->SetTitle("Shutter closed (pedestal data);ADC (raw ch sum);Events");
     h1_rawCor->SetTitle("Pedestal-corrected;Event No.;ADC (raw, ped-corr ch sum)");
     h1_asym_my->SetTitle("Pair-wise Asymmetry vs. pair No.;Pair No.;A_{LED}");
     h1_ped->SetLineColor(kBlack);
     h1_raw_even->SetMarkerColor(kBlue);
     h1_raw_odd->SetMarkerColor(kRed);
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
  */  
   }
