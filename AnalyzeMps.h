#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

   // Declaration of leaf types
   Double_t        CodaEventNumber;
   Double_t        CodaEventType;
   Double_t        Coda_CleanData;
   Double_t        Coda_ScanData1;
   Double_t        Coda_ScanData2;
   Double_t        qwk_mod0ch0_hw_sum;
   Double_t        qwk_mod0ch0_block0;
   Double_t        qwk_mod0ch0_block1;
   Double_t        qwk_mod0ch0_block2;
   Double_t        qwk_mod0ch0_block3;
   Double_t        qwk_mod0ch0_num_samples;
   Double_t        qwk_mod0ch0_Device_Error_Code;
   Double_t        qwk_mod0ch0_hw_sum_raw;
   Double_t        qwk_mod0ch0_block0_raw;
   Double_t        qwk_mod0ch0_block1_raw;
   Double_t        qwk_mod0ch0_block2_raw;
   Double_t        qwk_mod0ch0_block3_raw;
   Double_t        qwk_mod0ch0_sequence_number;
   Double_t        qwk_mod0ch1_hw_sum;
   Double_t        qwk_mod0ch1_block0;
   Double_t        qwk_mod0ch1_block1;
   Double_t        qwk_mod0ch1_block2;
   Double_t        qwk_mod0ch1_block3;
   Double_t        qwk_mod0ch1_num_samples;
   Double_t        qwk_mod0ch1_Device_Error_Code;
   Double_t        qwk_mod0ch1_hw_sum_raw;
   Double_t        qwk_mod0ch1_block0_raw;
   Double_t        qwk_mod0ch1_block1_raw;
   Double_t        qwk_mod0ch1_block2_raw;
   Double_t        qwk_mod0ch1_block3_raw;
   Double_t        qwk_mod0ch1_sequence_number;
   Double_t        qwk_mod0ch2_hw_sum;
   Double_t        qwk_mod0ch2_block0;
   Double_t        qwk_mod0ch2_block1;
   Double_t        qwk_mod0ch2_block2;
   Double_t        qwk_mod0ch2_block3;
   Double_t        qwk_mod0ch2_num_samples;
   Double_t        qwk_mod0ch2_Device_Error_Code;
   Double_t        qwk_mod0ch2_hw_sum_raw;
   Double_t        qwk_mod0ch2_block0_raw;
   Double_t        qwk_mod0ch2_block1_raw;
   Double_t        qwk_mod0ch2_block2_raw;
   Double_t        qwk_mod0ch2_block3_raw;
   Double_t        qwk_mod0ch2_sequence_number;
   Double_t        qwk_mod0ch3_hw_sum;
   Double_t        qwk_mod0ch3_block0;
   Double_t        qwk_mod0ch3_block1;
   Double_t        qwk_mod0ch3_block2;
   Double_t        qwk_mod0ch3_block3;
   Double_t        qwk_mod0ch3_num_samples;
   Double_t        qwk_mod0ch3_Device_Error_Code;
   Double_t        qwk_mod0ch3_hw_sum_raw;
   Double_t        qwk_mod0ch3_block0_raw;
   Double_t        qwk_mod0ch3_block1_raw;
   Double_t        qwk_mod0ch3_block2_raw;
   Double_t        qwk_mod0ch3_block3_raw;
   Double_t        qwk_mod0ch3_sequence_number;
   Double_t        qwk_mod0ch4_hw_sum;
   Double_t        qwk_mod0ch4_block0;
   Double_t        qwk_mod0ch4_block1;
   Double_t        qwk_mod0ch4_block2;
   Double_t        qwk_mod0ch4_block3;
   Double_t        qwk_mod0ch4_num_samples;
   Double_t        qwk_mod0ch4_Device_Error_Code;
   Double_t        qwk_mod0ch4_hw_sum_raw;
   Double_t        qwk_mod0ch4_block0_raw;
   Double_t        qwk_mod0ch4_block1_raw;
   Double_t        qwk_mod0ch4_block2_raw;
   Double_t        qwk_mod0ch4_block3_raw;
   Double_t        qwk_mod0ch4_sequence_number;
   Double_t        qwk_mod0ch5_hw_sum;
   Double_t        qwk_mod0ch5_block0;
   Double_t        qwk_mod0ch5_block1;
   Double_t        qwk_mod0ch5_block2;
   Double_t        qwk_mod0ch5_block3;
   Double_t        qwk_mod0ch5_num_samples;
   Double_t        qwk_mod0ch5_Device_Error_Code;
   Double_t        qwk_mod0ch5_hw_sum_raw;
   Double_t        qwk_mod0ch5_block0_raw;
   Double_t        qwk_mod0ch5_block1_raw;
   Double_t        qwk_mod0ch5_block2_raw;
   Double_t        qwk_mod0ch5_block3_raw;
   Double_t        qwk_mod0ch5_sequence_number;
   Double_t        qwk_mod0ch6_hw_sum;
   Double_t        qwk_mod0ch6_block0;
   Double_t        qwk_mod0ch6_block1;
   Double_t        qwk_mod0ch6_block2;
   Double_t        qwk_mod0ch6_block3;
   Double_t        qwk_mod0ch6_num_samples;
   Double_t        qwk_mod0ch6_Device_Error_Code;
   Double_t        qwk_mod0ch6_hw_sum_raw;
   Double_t        qwk_mod0ch6_block0_raw;
   Double_t        qwk_mod0ch6_block1_raw;
   Double_t        qwk_mod0ch6_block2_raw;
   Double_t        qwk_mod0ch6_block3_raw;
   Double_t        qwk_mod0ch6_sequence_number;
   Double_t        qwk_mod0ch7_hw_sum;
   Double_t        qwk_mod0ch7_block0;
   Double_t        qwk_mod0ch7_block1;
   Double_t        qwk_mod0ch7_block2;
   Double_t        qwk_mod0ch7_block3;
   Double_t        qwk_mod0ch7_num_samples;
   Double_t        qwk_mod0ch7_Device_Error_Code;
   Double_t        qwk_mod0ch7_hw_sum_raw;
   Double_t        qwk_mod0ch7_block0_raw;
   Double_t        qwk_mod0ch7_block1_raw;
   Double_t        qwk_mod0ch7_block2_raw;
   Double_t        qwk_mod0ch7_block3_raw;
   Double_t        qwk_mod0ch7_sequence_number;
   Double_t        ErrorFlag;

   // List of branches
   TBranch        *b_CodaEventNumber;
   TBranch        *b_CodaEventType; 
   TBranch        *b_Coda_CleanData;
   TBranch        *b_Coda_ScanData1;
   TBranch        *b_Coda_ScanData2;
   TBranch        *b_qwk_mod0ch0;
   TBranch        *b_qwk_mod0ch1;
   TBranch        *b_qwk_mod0ch2;
   TBranch        *b_qwk_mod0ch3;
   TBranch        *b_qwk_mod0ch4;
   TBranch        *b_qwk_mod0ch5;
   TBranch        *b_qwk_mod0ch6;
   TBranch        *b_qwk_mod0ch7;
   TBranch        *b_ErrorFlag;
   
   //User defined variable

   Int_t samples = 200;
   Int_t blocks = 4;
   Double_t Vperch = 76.29E-6; //Volts per channel
   const int positions = 160;
   const int filters = 8;
   const int wheelCycles = 20;
   vector<Double_t> EventNumber;
   vector<Double_t> sum_raw;

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
//   vector<Double_t> pairSum;
//   vector<Double_t> Filter_Event_Num;
//   vector<Double_t> Filter_Asymmetry;
//   vector<Double_t> Ped_Event_Num;
   Float_t patternNumber, pattern_hw_sum;
   Float_t Filter_Event_Num, Filter_Asymmetry,Ped_Event_Num,pairSum;
   int group = 0;
   int group_data_count = 0;
   int filter_count = 0;

