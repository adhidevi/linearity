#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <fstream>

void AnalyzeMps_octet(Int_t run,Int_t pmt,Double_t preAmp,Int_t HV,Double_t LL,Double_t Iout,Int_t FF,Double_t asym_X1,Double_t asym_X2,Double_t chsum,Double_t asym1_X1,Double_t asym1_X2,Double_t asym8_X1,Double_t asym8_X2,Double_t asymF_X1,Double_t asymF_X2,Double_t asymP_X1,Double_t asymP_X2,Double_t asym_Y1,Double_t asym_Y2)
{

     TString root_file = Form("/home/daq/scratch/rootfiles/r%d.root",run);
     TChain* tree = new TChain("evt");
     tree->Add(root_file);
     if(root_file==NULL){
     cout <<"Rootfile "<<root_file<<" doesn't exist"<<endl;
     exit(0);
     }
     if(tree==NULL){
     cout<<"Tree "<<tree<<" doesn't exit"<<endl;
     }
     Long64_t nentries = tree->GetEntries(); 
     cout<<"Total Entries: \t"<<nentries<<endl;
     gStyle->SetOptStat(1111);

     TH2F *h1_raw = new TH2F("h1_raw", "", nentries/100, 0, nentries, 1000, 0, 262000000);
     TH2F *h1_raw_even = new TH2F("h1_raw_even", "", nentries/100, 0, nentries, 1000, 0, 262000000);
     TH2F *h1_raw_odd = new TH2F("h1_raw_odd", "", nentries/100, 0, nentries, 1000, 0, 262000000);
     TH1F *h1_ped;
     TH1F *h1_ped_even;
     TH1F *h1_ped_odd;
     TH2F *h1_rawCor = new TH2F("h1_rawCor", "", nentries/10, 0, nentries, 1000, 0, 262000000);
     TH2F *h1_rawCor_even = new TH2F("h1_rawCor_even", "", nentries/10, 0, nentries, 1000, 0, 262000000);
     TH2F *h1_rawCor_odd = new TH2F("h1_rawCor_odd", "", nentries/10, 0, nentries, 1000, 0, 262000000);
     TH2F *h1_asym_my = new TH2F("h1_asym_my", "", nentries/20, 0, nentries/8, 800, asym_X1, asym_X2);

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
        h[f] = new TH1S(name,"",50, asym1_X1, asym1_X2);
     }else if((f-1)/8.0 - int((f-1)/8) == 0){
        h[f] = new TH1S(name,"",50, asym8_X1, asym8_X2);
     }else{
        h[f] = new TH1S(name,"",50, asymF_X1, asymF_X2);
     }
//        h[f]->SetMinimum(0.0);
//        h[f]->SetMaximum(120);
        h[f]->SetYTitle("Events");
        h[f]->SetXTitle("A_{LED}");
     }
    for(f=0;f<filters;f++){
        sprintf(name,"havg_asym%d",f);
        if(f==0)
        havg[f] = new TH1S(name,"",100, asym1_X1, asym1_X2);
        else
        havg[f] = new TH1S(name,"",100, asymF_X1, asymF_X2);

//        havg[f]->SetMinimum(0.0);
//        havg[f]->SetMaximum(8);
        havg[f]->SetYTitle("Events");
        havg[f]->SetXTitle("A_{LED}");
     }

     Int_t samples = 400;
     Int_t blocks = 4;
     Double_t Vperch = 76.29E-6; //Volts per channel
     const Int_t arraysize = nentries;
     const Int_t pairsize = arraysize/2;
     const Int_t octsize = pairsize;
     Double_t* hw_sum_raw = new Double_t[nentries];
     Double_t* EventNum = new Double_t[nentries];
     ofstream outfile("diff.txt");
     tree->Draw("qwk_mod0ch0.hw_sum_raw:CodaEventNumber","","goff");

   for (Long64_t jentry=10; jentry<nentries;jentry++) {
          EventNum[jentry] = *(tree->GetV2()+jentry);
          hw_sum_raw[jentry] = *(tree->GetV1()+jentry);
	  h1_raw->Fill(EventNum[jentry], hw_sum_raw[jentry]);
	if(EventNum[jentry]/2.0-int(EventNum[jentry]/2) ==0){
	  h1_raw_even->Fill(EventNum[jentry], hw_sum_raw[jentry]);
	}else{
	  h1_raw_odd->Fill(EventNum[jentry], hw_sum_raw[jentry]);
	}
	}

	TCut ped_cut = Form("CodaEventNumber>10&&qwk_mod0ch0.hw_sum_raw<%f",chsum);
	TCut even_cut = Form("qwk_mod0ch0.hw_sum_raw<%f&&Entry$>10&&(Entry$/2.0-int(Entry$/2))==0",chsum);	  
	TCut odd_cut = Form("qwk_mod0ch0.hw_sum_raw<%f&&Entry$>10&&(Entry$/2.0-int(Entry$/2))!=0",chsum);	  
	tree->Draw("qwk_mod0ch0.hw_sum_raw>>hped",ped_cut,"goff");
	h1_ped = (TH1F*)gDirectory->FindObject("hped");
	h1_ped->SetDirectory(gROOT);
	tree->Draw("qwk_mod0ch0.hw_sum_raw>>hpedE",even_cut,"goff");
	h1_ped_even = (TH1F*)gDirectory->FindObject("hpedE");
	h1_ped_even->SetDirectory(gROOT);
	tree->Draw("qwk_mod0ch0.hw_sum_raw>>hpedO",odd_cut,"goff");
	h1_ped_odd = (TH1F*)gDirectory->FindObject("hpedO");
	h1_ped_odd->SetDirectory(gROOT);

	Double_t pedAve = h1_ped->GetMean();
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
	Double_t Filter_Event_Num[pairsize];
	Double_t Ped_Event_Num[pairsize];
	int group = 0;
	int group_data_count = 0;

	Double_t ChSum_even = 0;
	Double_t ChSum_odd = 0;
	Double_t* EvenSum = new Double_t[octsize];
	Double_t* OddSum = new Double_t[octsize];
	Double_t* raw_sum0 = new Double_t[octsize];
	Double_t* raw_sum1 = new Double_t[octsize];
	Double_t* raw_sum2 = new Double_t[octsize];
	Double_t* raw_sum3 = new Double_t[octsize];
	Double_t* raw_sum4 = new Double_t[octsize];
	Double_t* raw_sum5 = new Double_t[octsize];
	Double_t* raw_sum6 = new Double_t[octsize];
	Double_t* raw_sum7 = new Double_t[octsize];
	Int_t Q0 = 0;
	Int_t Q1 = 0;
	Int_t Q2 = 0;
	Int_t Q3 = 0;
	Int_t Q4 = 0;
	Int_t Q5 = 0;
	Int_t Q6 = 0;
	Int_t Q7 = 0;

	for (Long64_t jentry=0; jentry<arraysize;jentry++) {
	  h1_rawCor->Fill(EventNum[jentry], hw_sum_raw[jentry]-pedAve);
	if(EventNum[jentry]/2.0-int(EventNum[jentry]/2) ==0){
	  h1_rawCor_even->Fill(EventNum[jentry], hw_sum_raw[jentry]-pedAve);
	}else{
	  h1_rawCor_odd->Fill(EventNum[jentry], hw_sum_raw[jentry]-pedAve);
	}
	  if(EventNum[jentry]/8.0 - int(EventNum[jentry]/8) ==0){
	  Q0++;
	  raw_sum0[Q0] = hw_sum_raw[jentry]-pedAve;
	  }else if((EventNum[jentry]-1)/8.0 - int((EventNum[jentry]-1)/8) ==0){
	  Q1++;
	  raw_sum1[Q1] = hw_sum_raw[jentry]-pedAve;
	  }else if((EventNum[jentry]-2)/8.0 - int((EventNum[jentry]-2)/8) ==0){
	  Q2++;
	  raw_sum2[Q2] = hw_sum_raw[jentry]-pedAve;
	  }else if((EventNum[jentry]-3)/8.0 - int((EventNum[jentry]-3)/8) ==0){
	  Q3++;
	  raw_sum3[Q3] = hw_sum_raw[jentry]-pedAve;
	  }else if((EventNum[jentry]-4)/8.0 - int((EventNum[jentry]-4)/8) ==0){
	  Q4++;
	  raw_sum4[Q4] = hw_sum_raw[jentry]-pedAve;
	  }else if((EventNum[jentry]-5)/8.0 - int((EventNum[jentry]-5)/8) ==0){
	  Q5++;
	  raw_sum5[Q5] = hw_sum_raw[jentry]-pedAve;
	  }else if((EventNum[jentry]-6)/8.0 - int((EventNum[jentry]-6)/8) ==0){
	  Q6++;
	  raw_sum6[Q6] = hw_sum_raw[jentry]-pedAve;
	  }else if((EventNum[jentry]-7)/8.0 - int((EventNum[jentry]-7)/8) ==0){
	  Q7++;
	  raw_sum7[Q7] = hw_sum_raw[jentry]-pedAve;
	  }
	}
	
    for (int i=1; i<octsize/4; i++){
    	  EvenSum[i] = raw_sum0[i]+raw_sum2[i]+raw_sum4[i]+raw_sum6[i];
    	  ChSum_even = EvenSum[i];
//	  outfile << i << "	" << EvenSum[i] << endl;
    	  OddSum[i] = raw_sum1[i]+raw_sum3[i]+raw_sum5[i]+raw_sum7[i];
    	  ChSum_odd = OddSum[i];
	  sum = EvenSum[i] + OddSum[i];
	  diff = EvenSum[i] - OddSum[i];
	  Asymmetry = fabs(diff/sum);
	  Asymmetry_Ped = (EvenSum[i] - OddSum[i])/(EvenSum[i] + OddSum[i]+8*pedAve);
          h1_asym_my->Fill(i, Asymmetry);
//          outfile2 << i << "	" << sum << endl;
       if(sum > 2.0E6){
	  Filter_Event_Num[i] = i;
       }else{
	  Ped_Event_Num[i] = i; 
       }
       if(Filter_Event_Num[i]!=Filter_Event_Num[i-1]+1 && Filter_Event_Num[i]!=0){
	  group++;
	  cout << group << "	" << i << "	" << group_data_count << "	" << Asymmetry << "	" << sum << endl;
	  group_data_count = 0;
       }else if(Filter_Event_Num[i]!=0){
	  group_data_count++;
       if(group_data_count>10 && group_data_count<=290){
//	outfile << pairNumber << "	" << group << "	" << group_data_count << "	" << Asymmetry << endl;
	  h[group-1]->Fill(Asymmetry);
	  Navg[group-1] += sum/8/280/samples/blocks;
       }
       }
       if(Ped_Event_Num[i]!=0){
	  hPed->Fill(Asymmetry_Ped);
       }
       }
      
	outfile.close();
//	  outfile2.close();

     h1_raw->SetTitle(Form("Filter wheel scan: pmt#ZK%4d @ %5d V, I_{cath} = %2.2f nA;Event No.;ADC (raw ch sum)", pmt, HV, LL));
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
   for(int igp = 0;igp<8;igp++)
     cout<<Navg[igp]<<endl;

    TCanvas *c = new TCanvas("c","Mps hwsum Tree Distributions",1000,650);
    c->Clear();
    c->Divide(2,2); 
    c->cd(1);
    gPad->SetGridx(1);
    gPad->SetGridy(1);
    h1_raw->Draw();
    h1_raw_odd->Draw("same");
    h1_raw_even->Draw("same");
    TPaveLabel *pt_rl = new TPaveLabel(0.25,0.85,0.58,0.89,Form("run %4d, PreAmp @ %3.2f M#Ome    ga LED flash @ %4d Hz:", run, preAmp,FF),"NDC");
    pt_rl->SetBorderSize(0);
    pt_rl->SetTextColor(kBlack);
    pt_rl->SetTextSize(0.80);
    pt_rl->SetFillColor(0);
    pt_rl->Draw();
    c->cd(2);
    gPad->SetGridx(1);
    gPad->SetGridy(1);
    h1_ped->Draw();
    h1_ped_odd->Draw("same");
    h1_ped_even->Draw("same");
    c->cd(3);
    gPad->SetGridx(1);
    gPad->SetGridy(1);
    h1_rawCor->Draw();
    h1_rawCor_odd->Draw("same");
    h1_rawCor_even->Draw("same");
    c->cd(4);
    gPad->SetGridx(1);
    gPad->SetGridy(1);
    h1_asym_my->Draw();
    c->SaveAs(Form("plots/r%dplotsMps.png",run));
 
    TCanvas *c21 = new TCanvas("c21","Mps Tree Distributions_auto_20cyc",1000,650);
    c21->Clear();
    c21->Divide(3,3); 
    gStyle->SetStatW(0.30);
    gStyle->SetStatH(0.30);
    for(int i=0; i<wheelCycles; i++){
    c21->cd(1);
    gPad->SetGridx(1);
    gPad->SetGridy(1);
    h[8*i+2]->Draw();
    char label[filters][256], labelF[filters][256], labelFPed[256], labelRun[256];;
    sprintf(label[0],"N_avg = %4.3e (%4.1f %%)" ,Navg[8*i+2], Navg[8*i+2]/Navg[8*i+0]*100);
    TPaveLabel *pt1 = new TPaveLabel(0.15,0.84, 0.50,0.89,label[0],"NDC");
    pt1->SetBorderSize(0);
    pt1->SetTextColor(kBlack);
    pt1->SetTextSize(1.0);
    pt1->SetFillColor(0);
    pt1->Draw();
    TPaveLabel *ptF1 = new TPaveLabel(0.15,0.75, 0.40,0.81,"1 #pm 0.1 %% ND filter","NDC");
    ptF1->SetBorderSize(0);
    ptF1->SetTextColor(kMagenta);
    ptF1->SetTextSize(1.0);
    ptF1->SetFillColor(0);
    ptF1->Draw();
    TPaveLabel *ptR1 = new TPaveLabel(0.12,0.91, 0.68,0.99,Form("run %4d, %5d V, I_pe #approx %2.2f nA, %3d Hz", run, HV, LL, FF),"NDC");
    ptR1->SetBorderSize(0);
    ptR1->SetTextColor(kRed);
    ptR1->SetTextSize(0.80);
    ptR1->SetFillColor(0);
    ptR1->Draw();

    c21->cd(2);
    gPad->SetGridx(1);
    gPad->SetGridy(1);
    h[8*i+3]->Draw();
    TPaveLabel *pt2 = new TPaveLabel(0.15,0.84, 0.55,0.89,Form("N_avg = %4.3e (%4.1f %%)" ,Navg[8*i+3], Navg[8*i+3]/Navg[8*i+0]*100),"NDC");
    pt2->SetBorderSize(0);
    pt2->SetTextColor(kBlack);
    pt2->SetTextSize(1.0);
    pt2->SetFillColor(0);
    pt2->Draw();
    TPaveLabel *ptF2 = new TPaveLabel(0.15,0.75, 0.40,0.81,"50 #pm 5%% ND filter","NDC");
    ptF2->SetBorderSize(0);
    ptF2->SetTextColor(kMagenta);
    ptF2->SetTextSize(1.0);
    ptF2->SetFillColor(0);
    ptF2->Draw();

    c21->cd(3);
    gPad->SetGridx(1);
    gPad->SetGridy(1);
    h[8*i+4]->Draw();
    TPaveLabel *pt3 = new TPaveLabel(0.15,0.84, 0.55,0.89,Form("N_avg = %4.3e (%4.1f %%)" ,Navg[8*i+4], Navg[8*i+4]/Navg[8*i+0]*100),"NDC");
    pt3->SetBorderSize(0);
    pt3->SetTextColor(kBlack);
    pt3->SetTextSize(1.0);
    pt3->SetFillColor(0);
    pt3->Draw();
    TPaveLabel *ptF3 = new TPaveLabel(0.15,0.75, 0.40,0.81,"63 #pm 5%% ND filter","NDC");
    ptF3->SetBorderSize(0);
    ptF3->SetTextColor(kMagenta);
    ptF3->SetTextSize(1.0);
    ptF3->SetFillColor(0);
    ptF3->Draw();

    c21->cd(4);
    gPad->SetGridx(1);
    gPad->SetGridy(1);
    h[8*i+5]->Draw();
    TPaveLabel *pt4 = new TPaveLabel(0.15,0.84, 0.55,0.89,Form("N_avg = %4.3e (%4.1f %%)" ,Navg[8*i+5], Navg[8*i+5]/Navg[8*i+0]*100),"NDC");
    pt4->SetBorderSize(0);
    pt4->SetTextColor(kBlack);
    pt4->SetTextSize(1.0);
    pt4->SetFillColor(0);
    pt4->Draw();
    TPaveLabel *ptF4 = new TPaveLabel(0.15,0.75, 0.40,0.81,"25 #pm 2.5%% ND filter","NDC");
    ptF4->SetBorderSize(0);
    ptF4->SetTextColor(kMagenta);
    ptF4->SetTextSize(1.0);
    ptF4->SetFillColor(0);
    ptF4->Draw();

    c21->cd(5);
    gPad->SetGridx(1);
    gPad->SetGridy(1);
    h[8*i+6]->Draw();
    TPaveLabel *pt5 = new TPaveLabel(0.15,0.84, 0.55,0.89,Form("N_avg = %4.3e (%4.1f %%)" ,Navg[8*i+6], Navg[8*i+6]/Navg[8*i+0]*100),"NDC");
    pt5->SetBorderSize(0);
    pt5->SetTextColor(kBlack);
    pt5->SetTextSize(1.0);
    pt5->SetFillColor(0);
    pt5->Draw();
    TPaveLabel *ptF5 = new TPaveLabel(0.15,0.75, 0.40,0.81,"79 #pm 5%% ND filter","NDC");
    ptF5->SetBorderSize(0);
    ptF5->SetTextColor(kMagenta);
    ptF5->SetTextSize(1.0);
    ptF5->SetFillColor(0);
    ptF5->Draw();

    c21->cd(6);
    gPad->SetGridx(1);
    gPad->SetGridy(1);
    h[8*i+7]->Draw();
    TPaveLabel *pt6 = new TPaveLabel(0.15,0.84, 0.55,0.89,Form("N_avg = %4.3e (%4.1f %%)" ,Navg[8*i+7], Navg[8*i+7]/Navg[8*i+0]*100),"NDC");
    pt6->SetBorderSize(0);
    pt6->SetTextColor(kBlack);
    pt6->SetTextSize(1.0);
    pt6->SetFillColor(0);
    pt6->Draw();
    TPaveLabel *ptF6 = new TPaveLabel(0.15,0.75, 0.40,0.81,"40 #pm 4%% ND filter","NDC");
    ptF6->SetBorderSize(0);
    ptF6->SetTextColor(kMagenta);
    ptF6->SetTextSize(1.0);
    ptF6->SetFillColor(0);
    ptF6->Draw();

    c21->cd(7);
    gPad->SetGridx(1);
    gPad->SetGridy(1);
    h[8*i+0]->Draw();
    TPaveLabel *pt7 = new TPaveLabel(0.15,0.84, 0.55,0.89,Form("N_avg = %4.3e (%4.1f %%)" ,Navg[8*i+0], Navg[8*i+0]/Navg[8*i+0]*100),"NDC");
    pt7->SetBorderSize(0);
    pt7->SetTextColor(kBlack);
    pt7->SetTextSize(1.0);
    pt7->SetFillColor(0);
    pt7->Draw();
    TPaveLabel *ptF7 = new TPaveLabel(0.15,0.75, 0.40,0.81,"100%% (No filter)","NDC");
    ptF7->SetBorderSize(0);
    ptF7->SetTextColor(kMagenta);
    ptF7->SetTextSize(1.0);
    ptF7->SetFillColor(0);
    ptF7->Draw();

    c21->cd(8);
    gPad->SetGridx(1);
    gPad->SetGridy(1);
    h[8*i+1]->Draw();
    TPaveLabel *pt8 = new TPaveLabel(0.15,0.84, 0.55,0.89,Form("N_avg = %4.3e (%4.1f %%)" ,Navg[8*i+1], Navg[8*i+1]/Navg[8*i+0]*100),"NDC");
    pt8->SetBorderSize(0);
    pt8->SetTextColor(kBlack);
    pt8->SetTextSize(1.0);
    pt8->SetFillColor(0);
    pt8->Draw();
    TPaveLabel *ptF8 = new TPaveLabel(0.15,0.75, 0.40,0.81,"10 #pm 1%% ND filter","NDC");
    ptF8->SetBorderSize(0);
    ptF8->SetTextColor(kMagenta);
    ptF8->SetTextSize(1.0);
    ptF8->SetFillColor(0);
    ptF8->Draw();

    c21->cd(9);
    gPad->SetGridx(1);
    gPad->SetGridy(1);
    hPed->Draw();

    TPaveLabel *ptF9 = new TPaveLabel(0.15,0.75, 0.40,0.81,"Shutter Closerd","NDC");
    ptF9->SetBorderSize(0);
    ptF9->SetTextColor(kMagenta);
    ptF9->SetTextSize(1.0);
    ptF9->SetFillColor(0);
    ptF9->Draw();
    c21->SaveAs(Form("plots/r%d_asyms_auto_20cyc_c%d.png",run,i));
	}
	
   for(f=0;f<positions;f++){
   	 Aled[f] = h[f]->GetMean();
     if((f-2)/8.0-int((f-2)/8)==0){
	 havg[0]->Fill(Aled[f]);
	 Navgmean[0] += Navg[f]/wheelCycles;
     } else if((f-3)/8.0-int((f-3)/8)==0){
	 havg[1]->Fill(Aled[f]);
	 Navgmean[1] += Navg[f]/wheelCycles;
     } else if((f-4)/8.0-int((f-4)/8)==0){
	 havg[2]->Fill(Aled[f]);
	 Navgmean[2] += Navg[f]/wheelCycles;
     } else if((f-5)/8.0-int((f-5)/8)==0){
	 havg[3]->Fill(Aled[f]);
	 Navgmean[3] += Navg[f]/wheelCycles;
     } else if((f-6)/8.0-int((f-6)/8)==0){
	 havg[4]->Fill(Aled[f]);
	 Navgmean[4] += Navg[f]/wheelCycles;
     } else if((f-7)/8.0-int((f-7)/8)==0){
	 havg[5]->Fill(Aled[f]);
	 Navgmean[5] += Navg[f]/wheelCycles;
     } else if((f-0)/8.0-int((f-0)/8)==0){
	 havg[6]->Fill(Aled[f]);
	 Navgmean[6] += Navg[f]/wheelCycles;
     } else if((f-1)/8.0-int((f-1)/8)==0){
	 havg[7]->Fill(Aled[f]);
	 Navgmean[7] += Navg[f]/wheelCycles;
     }
     }

   for(int i=0;i<filters;i++){
      Navgmean_Cal[i] = Navgmean[i]*Vperch;
      }

    TCanvas *c_avg = new TCanvas("c_avg","Mps Tree Distributions_auto_20cyc_c2",1000,600);
    c_avg->Clear();
    c_avg->Divide(4,2); 
    gStyle->SetStatW(0.30);
    gStyle->SetStatH(0.30);
    c_avg->cd(1);
    gPad->SetGridx(1);
    gPad->SetGridy(1);
    havg[0]->Draw();
    Aledmean[0] = havg[0]->GetMean();
    AledErrormean[0] = sqrt(pow(havg[0]->GetMeanError(),2));
    dummyError[0] = 0.0;
    TPaveLabel *pl1 = new TPaveLabel(0.25,0.85, 0.55,0.89,Form("N_avg = %4.3e (%4.1f %%)" ,Navgmean[0], Navgmean[0]/Navgmean[6]*100),"NDC");
    pl1->SetBorderSize(0);
    pl1->SetTextColor(kBlack);
    pl1->SetTextSize(1.0);
    pl1->SetFillColor(0);
    pl1->Draw();
    TPaveLabel *plF1 = new TPaveLabel(0.25,0.80, 0.38,0.84,"1 #pm 0.1% ND filter","NDC");
    plF1->SetBorderSize(0);
    plF1->SetTextColor(kMagenta);
    plF1->SetTextSize(1.0);
    plF1->SetFillColor(0);
    plF1->Draw();
    TPaveLabel *plR1 = new TPaveLabel(0.32,0.94, 0.65,0.99,Form("run %4d, %5d V, I_pe #approx %2.2f nA, %3d Hz", run, HV, LL, FF),"NDC");
    plR1->SetBorderSize(0);
    plR1->SetTextColor(kRed);
    plR1->SetTextSize(0.80);
    plR1->SetFillColor(0);
    plR1->Draw();

    c_avg->cd(2);
    gPad->SetGridx(1);
    gPad->SetGridy(1);
    havg[1]->Draw();
    Aledmean[1] = havg[1]->GetMean();
    AledErrormean[1] = sqrt(pow(havg[1]->GetMeanError(),2));
    dummyError[1] = 0.0;
    TPaveLabel *pl2 = new TPaveLabel(0.25,0.85, 0.55,0.89,Form("N_avg = %4.3e (%4.1f %%)" ,Navgmean[1], Navgmean[1]/Navgmean[6]*100),"NDC");
    pl2->SetBorderSize(0);
    pl2->SetTextColor(kBlack);
    pl2->SetTextSize(1.0);
    pl2->SetFillColor(0);
    pl2->Draw();
    TPaveLabel *plF2 = new TPaveLabel(0.25,0.80, 0.38,0.84,"50 #pm 5%% ND filter","NDC");
    plF2->SetBorderSize(0);
    plF2->SetTextColor(kMagenta);
    plF2->SetTextSize(1.0);
    plF2->SetFillColor(0);
    plF2->Draw();

    c_avg->cd(3);
    gPad->SetGridx(1);
    gPad->SetGridy(1);
    havg[2]->Draw();
    Aledmean[2] = havg[2]->GetMean();
    AledErrormean[2] = sqrt(pow(havg[2]->GetMeanError(),2));
    dummyError[2] = 0.0;
    TPaveLabel *pl3 = new TPaveLabel(0.25,0.85, 0.55,0.89,Form("N_avg = %4.3e (%4.1f %%)" ,Navgmean[2], Navgmean[2]/Navgmean[6]*100),"NDC");
    pl3->SetBorderSize(0);
    pl3->SetTextColor(kBlack);
    pl3->SetTextSize(1.0);
    pl3->SetFillColor(0);
    pl3->Draw();
    TPaveLabel *plF3 = new TPaveLabel(0.25,0.80, 0.38,0.84,"63 #pm 5%% ND filter","NDC");
    plF3->SetBorderSize(0);
    plF3->SetTextColor(kMagenta);
    plF3->SetTextSize(1.0);
    plF3->SetFillColor(0);
    plF3->Draw();

    c_avg->cd(4);
    gPad->SetGridx(1);
    gPad->SetGridy(1);
    havg[3]->Draw();
    Aledmean[3] = havg[3]->GetMean();
    AledErrormean[3] = sqrt(pow(havg[3]->GetMeanError(),2));
    dummyError[3] = 0.0;
    TPaveLabel *pl4 = new TPaveLabel(0.25,0.85, 0.55,0.89,Form("N_avg = %4.3e (%4.1f %%)" ,Navgmean[3], Navgmean[3]/Navgmean[6]*100),"NDC");
    pl4->SetBorderSize(0);
    pl4->SetTextColor(kBlack);
    pl4->SetTextSize(1.0);
    pl4->SetFillColor(0);
    pl4->Draw();
    TPaveLabel *plF4 = new TPaveLabel(0.25,0.80, 0.38,0.84,"25 #pm 2.5%% ND filter","NDC");
    plF4->SetBorderSize(0);
    plF4->SetTextColor(kMagenta);
    plF4->SetTextSize(1.0);
    plF4->SetFillColor(0);
    plF4->Draw();

    c_avg->cd(5);
    gPad->SetGridx(1);
    gPad->SetGridy(1);
    havg[4]->Draw();
    Aledmean[4] = havg[4]->GetMean();
    AledErrormean[4] = sqrt(pow(havg[4]->GetMeanError(),2));
    dummyError[4] = 0.0;
    TPaveLabel *pl5 = new TPaveLabel(0.25,0.85, 0.55,0.89,Form("N_avg = %4.3e (%4.1f %%)" ,Navgmean[4], Navgmean[4]/Navgmean[6]*100),"NDC");
    pl5->SetBorderSize(0);
    pl5->SetTextColor(kBlack);
    pl5->SetTextSize(1.0);
    pl5->SetFillColor(0);
    pl5->Draw();
    TPaveLabel *plF5 = new TPaveLabel(0.25,0.80, 0.38,0.84,"79 #pm 5%% ND filter","NDC");
    plF5->SetBorderSize(0);
    plF5->SetTextColor(kMagenta);
    plF5->SetTextSize(1.0);
    plF5->SetFillColor(0);
    plF5->Draw();

    c_avg->cd(6);
    gPad->SetGridx(1);
    gPad->SetGridy(1);
    havg[5]->Draw();
    Aledmean[5] = havg[5]->GetMean();
    AledErrormean[5] = sqrt(pow(havg[5]->GetMeanError(),2));
    dummyError[5] = 0.0;
    TPaveLabel *pl6 = new TPaveLabel(0.25,0.85, 0.55,0.89,Form("N_avg = %4.3e (%4.1f %%)" ,Navgmean[5], Navgmean[5]/Navgmean[6]*100),"NDC");
    pl6->SetBorderSize(0);
    pl6->SetTextColor(kBlack);
    pl6->SetTextSize(1.0);
    pl6->SetFillColor(0);
    pl6->Draw();
    TPaveLabel *plF6 = new TPaveLabel(0.25,0.80, 0.38,0.84,"40 #pm 4%% ND filter","NDC");
    plF6->SetBorderSize(0);
    plF6->SetTextColor(kMagenta);
    plF6->SetTextSize(1.0);
    plF6->SetFillColor(0);
    plF6->Draw();

    c_avg->cd(7);
    gPad->SetGridx(1);
    gPad->SetGridy(1);
    havg[6]->Draw();
    Aledmean[6] = havg[6]->GetMean();
    AledErrormean[6] = sqrt(pow(havg[6]->GetMeanError(),2));
    dummyError[6] = 0.0;
    TPaveLabel *pl7 = new TPaveLabel(0.25,0.85, 0.55,0.89,Form("N_avg = %4.3e (%4.1f %%)" ,Navgmean[6], Navgmean[6]/Navgmean[6]*100),"NDC");
    pl7->SetBorderSize(0);
    pl7->SetTextColor(kBlack);
    pl7->SetTextSize(1.0);
    pl7->SetFillColor(0);
    pl7->Draw();
    TPaveLabel *plF7 = new TPaveLabel(0.25,0.80, 0.38,0.84,"100%% (No filter)","NDC");
    plF7->SetBorderSize(0);
    plF7->SetTextColor(kMagenta);
    plF7->SetTextSize(1.0);
    plF7->SetFillColor(0);
    plF7->Draw();

    c_avg->cd(8);
    gPad->SetGridx(1);
    gPad->SetGridy(1);
    havg[7]->Draw();
    Aledmean[7] = havg[7]->GetMean();
    AledErrormean[7] = sqrt(pow(havg[7]->GetMeanError(),2));
    dummyError[7] = 0.0;
    TPaveLabel *pl8 = new TPaveLabel(0.25,0.85, 0.55,0.89,Form("N_avg = %4.3e (%4.1f %%)" ,Navgmean[7], Navgmean[7]/Navgmean[6]*100),"NDC");
    pl8->SetBorderSize(0);
    pl8->SetTextColor(kBlack);
    pl8->SetTextSize(1.0);
    pl8->SetFillColor(0);
    pl8->Draw();
    TPaveLabel *plF8 = new TPaveLabel(0.25,0.80, 0.38,0.84,"10 #pm 1 %% ND filter","NDC");
    plF8->SetBorderSize(0);
    plF8->SetTextColor(kMagenta);
    plF8->SetTextSize(1.0);
    plF8->SetFillColor(0);
    plF8->Draw();
    c_avg->SaveAs(Form("plots/r%d_asyms_auto_20cyc_avg.png",run));

   for(int i=0;i<filters;i++){
      cout<<Navgmean[i]<<"	"<<Navgmean_Cal[i]<<"	"<< Aledmean[i] <<"	"<<AledErrormean[i]<<endl;
      }
      
   TCanvas *c3 = new TCanvas("c3","Non-linearity fit_auto_20cyc_avg",0,0,900,600);
   c3->SetGridx(1);
   c3->SetGridy(1);
   gStyle->SetOptStat(kFALSE);
   TH2F *hr1;
   hr1 = new TH2F("hr1","",100,0,131072,100,asym_Y1,asym_Y2);
   hr1->Draw();
   TGraphErrors *gr1 = new TGraphErrors(8, Navgmean, Aledmean, dummyError, AledErrormean);
   gr1->SetMarkerColor(kBlue);
   gr1->SetMarkerStyle(20);   
   gr1->SetMarkerSize(0.75);
   gStyle->SetOptFit(111);
   gStyle->SetStatW(0.2);
   gStyle->SetStatH(0.15);
   gr1->Draw("P");
   TF1 *f1 = new TF1("f1","pol1",0,131000);
   f1->SetLineWidth(2);
   f1->SetLineColor(kMagenta);
   f1->SetParName(0,"A_{true}");
   f1->SetParName(1,"A_{true}*#beta");
   Double_t f1par[2];
   gr1->Fit("f1","R"); 
   f1->GetParameters(&f1par[0]);
   Double_t Error00 = f1->GetParError(0);
   Double_t Error01 = f1->GetParError(1);
   Double_t scale = sqrt(f1->GetChisquare()/f1->GetNDF());
   cout << "Scale from Chi2 value is: " << scale << endl;
   char fitLabel[256], LabelCal[256];
   sprintf(fitLabel,"run %4d, %3d Hz, %3.2f M#Omega, I_{anode} = %4.2f #mu A", run, FF, preAmp, Iout);
   TPaveLabel *pt_rlfit = new TPaveLabel(0.20,0.84,0.50,0.89,fitLabel,"NDC");
   pt_rlfit->SetBorderSize(0);
   pt_rlfit->SetTextColor(kBlack);
   pt_rlfit->SetTextSize(0.75);
   pt_rlfit->SetFillColor(0);
   pt_rlfit->Draw();
   sprintf(LabelCal,"%% non-lin = %5.3f V^{-1}", f1par[1]/f1par[0]*Navgmean[6]*100/Navgmean_Cal[6]);
   TPaveLabel *pt_cal = new TPaveLabel(0.68,0.70,0.85,0.74,LabelCal,"NDC");
   pt_cal->SetBorderSize(0);
   pt_cal->SetTextColor(kMagenta);
   pt_cal->SetTextSize(0.75);
   pt_cal->SetFillColor(0);
   pt_cal->Draw();

   char labelB[256], labelChi[256], labelP0[256];
   sprintf(labelB,"#beta*N_{avg} @ %4.3e = %5.3f #pm %4.3f %% non-lin",Navgmean[6], f1par[1]/f1par[0]*Navgmean[6]*100, fabs(sqrt(pow(Error00/f1par[0],2) + pow(Error01/f1par[1],2))*f1par[1]/f1par[0]*Navgmean[6]*100)*scale);
   TPaveLabel *ptB = new TPaveLabel(0.11,0.78, 0.53,0.82,labelB,"NDC");
   ptB->SetBorderSize(0);
   ptB->SetTextColor(kMagenta);
   ptB->SetTextSize(0.80);
   ptB->SetFillColor(0);
   ptB->Draw();
   sprintf(labelChi,"#chi^{2}                                                %4.4f",f1->GetChisquare()/f1->GetNDF());
   TPaveLabel *ptChi = new TPaveLabel(0.63,0.888, 0.97,0.928,labelChi,"NDC");
   ptChi->SetBorderSize(0);
   ptChi->SetTextColor(kBlack);
   ptChi->SetTextSize(0.80);
   ptChi->SetFillColor(0);
   ptChi->Draw();
   hr1->SetTitle(Form("Filter wheel scan: pmt#ZK%4d @ %5d V, I_{cath} = %2.2f nA;N_{avg} (adc channels);A_{LED}", pmt, HV, LL));   
   TLegend *leg = new TLegend(0.55,0.10,0.89,0.30);
   leg->SetBorderSize(0);
   leg->SetFillColor(0);
   leg->SetFillStyle(0);
   leg->AddEntry(gr1,"data","pe");
   leg->AddEntry(f1,"fit: A_{LED} = A_{true}(1+#beta*N_{avg})","l");
   leg->Draw();
   c3->SaveAs(Form("plots/r%d_asymTrends_auto_20cyc.png",run));   
 
   TCanvas *c3Cal = new TCanvas("c3Cal","Non-linearity fit_auto_20cyc_avg_Cal",0,0,900,600);
   c3Cal->SetGridx(1);
   c3Cal->SetGridy(1);
   gStyle->SetOptStat(kFALSE);
   TH2F *hr2;
   hr2 = new TH2F("hr2","",100,0,10,100,asym_Y1,asym_Y2);
   hr2->Draw();
   TGraphErrors *gr3 = new TGraphErrors(8, Navgmean_Cal, Aledmean, dummyError, AledErrormean);
   gr3->SetMarkerColor(kBlue);
   gr3->SetMarkerStyle(20);   
   gr3->SetMarkerSize(0.75);
   gStyle->SetOptFit(111);
   gStyle->SetStatW(0.2);
   gStyle->SetStatH(0.15);
   gr3->Draw("P");
   TF1 *f3 = new TF1("f3","pol1",0,10);
   f3->SetLineWidth(2);
   f3->SetLineColor(kMagenta);
   f3->SetParName(0,"A_{true}");
   f3->SetParName(1,"A_{true}*#beta");
   gr3->Fit("f3","R"); 
   Double_t f3par[2];
   f3->GetParameters(&f3par[0]);
   Double_t Error20 = f3->GetParError(0);
   Double_t Error21 = f3->GetParError(1);
   char fitLabelA[256], LabelCalA[256];
   sprintf(LabelCalA,"%% non-lin = %5.3f V^{-1}", f1par[1]/f1par[0]*Navgmean[6]*100/Navgmean_Cal[6]);
   TPaveLabel *pt_calA = new TPaveLabel(0.68,0.70,0.85,0.74,LabelCalA,"NDC");
   pt_calA->SetBorderSize(0);
   pt_calA->SetTextColor(kMagenta);
   pt_calA->SetTextSize(0.75);
   pt_calA->SetFillColor(0);
   pt_calA->Draw();
   sprintf(fitLabelA,"run %4d, %3d Hz, %3.2f M#Omega, I_{anode} = %4.2f #mu A", run, FF, preAmp, Iout);
   TPaveLabel *pt_rlcal = new TPaveLabel(0.20,0.84,0.50,0.89,fitLabelA,"NDC");
   pt_rlcal->SetBorderSize(0);
   pt_rlcal->SetTextColor(kBlack);
   pt_rlcal->SetTextSize(0.75);
   pt_rlcal->SetFillColor(0);
   pt_rlcal->Draw();
   char labelC[256], labelCC[256];
   sprintf(labelC,"#beta*N_{avg} @ %4.3f = %5.3f #pm %4.3f %% non-lin",Navgmean_Cal[6], f3par[1]/f3par[0]*Navgmean_Cal[6]*100, fabs(sqrt(pow(Error20/f3par[0],2) + pow(Error21/f3par[1],2))*f3par[1]/f3par[0]*Navgmean_Cal[6]*100)*scale);
   TPaveLabel *plB = new TPaveLabel(0.11,0.78, 0.53,0.82,labelC,"NDC");
   plB->SetBorderSize(0);
   plB->SetTextColor(kMagenta);
   plB->SetTextSize(0.80);
   plB->SetFillColor(0);
   plB->Draw();
   sprintf(labelCC,"#chi^{2}                                                %4.4f",f3->GetChisquare()/f3->GetNDF());
   TPaveLabel *ptBB = new TPaveLabel(0.63,0.888, 0.97,0.928,labelCC,"NDC");
   ptBB->SetBorderSize(0);
   ptBB->SetTextColor(kBlack);
   ptBB->SetTextSize(0.80);
   ptBB->SetFillColor(0);
   ptBB->Draw();
   hr2->SetTitle(Form("Filter wheel scan: pmt#ZK%4d @ %5d V, I_{cath} = %2.2f nA;ADC Voltage (V);A_{LED}", pmt, HV, LL));
   TLegend *legA = new TLegend(0.55,0.10,0.89,0.30);
   legA->SetBorderSize(0);
   legA->SetFillColor(0);
   legA->SetFillStyle(0);
   legA->AddEntry(gr1,"data","pe");
   legA->AddEntry(f3,"fit: A_{LED} = A_{true}(1+#beta*N_{avg})","l");
   legA->Draw();
   c3Cal->SaveAs(Form("plots/r%d_asymTrends_auto_20cyc_Cal.png",run));

   }
