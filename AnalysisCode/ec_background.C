#include <iostream>
#include <cmath>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TLine.h>
#include <TGraph.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TChain.h>
#include "TFile.h"
#include "TApplication.h"
#include <TROOT.h>

using namespace std; 
TApplication gui("GUI",0,NULL);

int main() {

   cout << "Setting Style" << endl;
   
   gStyle->SetCanvasBorderMode(0);
   gStyle->SetCanvasColor(10);

   gStyle->SetPadBorderMode(0);
   gStyle->SetPadLeftMargin(0.1);
   gStyle->SetPadRightMargin(0.12);
   gStyle->SetPadTopMargin(0.1);
   gStyle->SetPadBottomMargin(0.13);
   gStyle->SetPadColor(10);


   gStyle->SetTitleFont(72,"X");
   gStyle->SetTitleFont(72,"Y");
   gStyle->SetTitleOffset(0.7,"X");
   gStyle->SetTitleOffset(0.5,"Y");
   gStyle->SetTitleSize(0.05,"X");
   gStyle->SetTitleSize(0.05,"Y");

   gStyle->SetLabelFont(72,"X");
   gStyle->SetLabelFont(72,"Y");
   gStyle->SetLabelFont(72,"Z");
   //   gStyle->SetLabelSize(0.04,"X");
   //   gStyle->SetLabelSize(0.04,"Y");
   //   gStyle->SetLabelSize(0.04,"Z");
   gStyle->SetPalette(1);
   gStyle->SetOptFit(111);
   gStyle->SetOptStat("nemriou");
   gStyle->SetOptStat("");

   int   evn;
   int   ngen;

   const int NTOFMAX = 100000;
// pcal
//   int npcalhit;
//   vector<double> *pcal_sector=new vector<double>;
//   vector<double> *pcal_layer=new vector<double>;
//   vector<double> *pcal_view=new vector<double>;
//   vector<double> *pcal_strip=new vector<double>;
//   vector<double> *pcal_ADC=new vector<double>;
//   vector<double> *pcal_TDC=new vector<double>;
//   vector<double> *pcal_Edep=new vector<double>;
//   vector<double> *pcal_E=new vector<double>;
//   vector<double> *pcal_x=new vector<double>;
//   vector<double> *pcal_y=new vector<double>;
//   vector<double> *pcal_z=new vector<double>;
//   vector<double> *pcal_lx=new vector<double>;
//   vector<double> *pcal_ly=new vector<double>;
//   vector<double> *pcal_lz=new vector<double>;
//   vector<double> *pcal_t=new vector<double>;
//   vector<double> *pcal_pid=new vector<double>;
//   vector<double> *pcal_vx=new vector<double>;
//   vector<double> *pcal_vy=new vector<double>;
//   vector<double> *pcal_vz=new vector<double>;
// ec
   int nechit;
   vector<double> *ec_sector=new vector<double>;
   vector<double> *ec_layer=new vector<double>;
   vector<double> *ec_view=new vector<double>;
   vector<double> *ec_strip=new vector<double>;
   vector<double> *ec_ADC=new vector<double>;
   vector<double> *ec_TDC=new vector<double>;
   vector<double> *ec_Edep=new vector<double>;
   vector<double> *ec_E=new vector<double>;
   vector<double> *ec_x=new vector<double>;
   vector<double> *ec_y=new vector<double>;
   vector<double> *ec_z=new vector<double>;
   vector<double> *ec_lx=new vector<double>;
   vector<double> *ec_ly=new vector<double>;
   vector<double> *ec_lz=new vector<double>;
   vector<double> *ec_t=new vector<double>;
   vector<double> *ec_pid=new vector<double>;
   vector<double> *ec_vx=new vector<double>;
   vector<double> *ec_vy=new vector<double>;
   vector<double> *ec_vz=new vector<double>;


   cout << "Creating Tree chains" << endl;
//   TChain *pcal= new TChain("pcal");
//   pcal->Add("out*.root");
   TChain *ec= new TChain("ec");
   ec->Add("out*.root");


// PCAL
//   pcal->SetBranchAddress("sector" ,&pcal_sector);
//   pcal->SetBranchAddress("module" ,&pcal_layer);
//   pcal->SetBranchAddress("view"   ,&pcal_view);
//   pcal->SetBranchAddress("strip"  ,&pcal_strip);
//   pcal->SetBranchAddress("ADC"    ,&pcal_ADC);
//   pcal->SetBranchAddress("TDC"    ,&pcal_TDC);
//   pcal->SetBranchAddress("trackE" ,&pcal_E);
//   pcal->SetBranchAddress("totEdep",&pcal_Edep);
//   pcal->SetBranchAddress("avg_t"  ,&pcal_t);
//   pcal->SetBranchAddress("pid"    ,&pcal_pid);
//   pcal->SetBranchAddress("avg_x"  ,&pcal_x);
//   pcal->SetBranchAddress("avg_y"  ,&pcal_y);
//   pcal->SetBranchAddress("avg_z"  ,&pcal_z);
//   pcal->SetBranchAddress("avg_lx" ,&pcal_lx);
//   pcal->SetBranchAddress("avg_ly" ,&pcal_ly);
//   pcal->SetBranchAddress("avg_lz" ,&pcal_lz); 
//   pcal->SetBranchAddress("vx"     ,&pcal_vx);
//   pcal->SetBranchAddress("vy"     ,&pcal_vy);
//   pcal->SetBranchAddress("vz"     ,&pcal_vz);

// EC
   ec->SetBranchAddress("sector" ,&ec_sector);
   ec->SetBranchAddress("stack"  ,&ec_layer);
   ec->SetBranchAddress("view"   ,&ec_view);
   ec->SetBranchAddress("strip"  ,&ec_strip);
   ec->SetBranchAddress("ADC"    ,&ec_ADC);
   ec->SetBranchAddress("TDC"    ,&ec_TDC);
   ec->SetBranchAddress("trackE" ,&ec_E);
   ec->SetBranchAddress("totEdep",&ec_Edep);
   ec->SetBranchAddress("avg_t"  ,&ec_t);
   ec->SetBranchAddress("pid"    ,&ec_pid);
   ec->SetBranchAddress("avg_x"  ,&ec_x);
   ec->SetBranchAddress("avg_y"  ,&ec_y);
   ec->SetBranchAddress("avg_z"  ,&ec_z);
   ec->SetBranchAddress("avg_lx" ,&ec_lx);
   ec->SetBranchAddress("avg_ly" ,&ec_ly);
   ec->SetBranchAddress("avg_lz" ,&ec_lz); 
   ec->SetBranchAddress("vx"     ,&ec_vx);
   ec->SetBranchAddress("vy"     ,&ec_vy);
   ec->SetBranchAddress("vz"     ,&ec_vz);



   Long64_t nentries = ec->GetEntries();
   cout << "N. entries:" << nentries << " " << ec->GetEntries() << endl;



// Create histos
   TH2F *hi_ecal_occ1[6];
   TH2F *hi_ecal_occ_cut1[6]; 
   TH1F *hi_ecal_vz_all1[6];
   TH1F *hi_ecal_vz_e1[6];
   TH1F *hi_ecal_vz_g1[6];
   TH1F *hi_ecal_vz_h1[6];
   TH1F *hi_ecal_vz_n1[6];
   TH2F *hi_ecal_origin_all1[6]; 
   TH1F *hi_ecal_edep1[6];
   TH2F *hi_ecal_occ2[6];
   TH2F *hi_ecal_occ_cut2[6];
   TH1F *hi_ecal_vz_all2[6];
   TH1F *hi_ecal_vz_e2[6];
   TH1F *hi_ecal_vz_g2[6];
   TH1F *hi_ecal_vz_h2[6];
   TH1F *hi_ecal_vz_n2[6];
   TH2F *hi_ecal_origin_all2[6];
   TH1F *hi_ecal_edep2[6];
   TH2F *hi_ecal_occ3[6];
   TH2F *hi_ecal_occ_cut3[6];
   TH1F *hi_ecal_vz_all3[6];
   TH1F *hi_ecal_vz_e3[6];
   TH1F *hi_ecal_vz_g3[6];
   TH1F *hi_ecal_vz_h3[6];
   TH1F *hi_ecal_vz_n3[6];
   TH2F *hi_ecal_origin_all3[6];
   TH1F *hi_ecal_edep3[6]; 

   for(int i=0; i<6; i++) {
       int nstrip=36;
       hi_ecal_occ1[i] = new TH2F(Form("hi_ecal_occ_layer1 %d",i+1), Form("hi_ecal_occ_layer1 %d",i+1),nstrip, 1.,nstrip*1.+1.,3,1.,4.);
       hi_ecal_occ1[i]->GetXaxis()->SetTitle("Strip");
       hi_ecal_occ1[i]->GetYaxis()->SetTitle("View");
       hi_ecal_occ_cut1[i] = new TH2F(Form("hi_ecal_occ_cut_layer1 %d",i+1), Form("hi_ecal_occ_cut_layer1 %d",i+1),nstrip, 1.,nstrip*1.+1.,3,1.,4.);
       hi_ecal_occ_cut1[i]->GetXaxis()->SetTitle("Strip");
       hi_ecal_occ_cut1[i]->GetYaxis()->SetTitle("View");
       hi_ecal_vz_all1[i] = new TH1F(Form("hi_ecal_vz_all_layer1 %d",i+1), Form("hi_ecal_vz_all_layer1 %d",i+1),300, 0., 900.);
       hi_ecal_vz_all1[i]->GetXaxis()->SetTitle("Z_{vertex} (cm)");
       hi_ecal_vz_all1[i]->GetYaxis()->SetTitle("Rate (MHz/3cm)");
       hi_ecal_vz_e1[i] = new TH1F(Form("hi_ecal_vz_e_layer1 %d",i+1), Form("hi_ecal_vz_e_layer1 %d",i+1),300, 0., 900.);
       hi_ecal_vz_e1[i]->GetXaxis()->SetTitle("Z_{vertex} (cm)");
       hi_ecal_vz_e1[i]->GetYaxis()->SetTitle("Rate (MHz/3cm)");
       hi_ecal_vz_g1[i] = new TH1F(Form("hi_ecal_vz_g_layer1 %d",i+1), Form("hi_ecal_vz_g_layer1 %d",i+1),300, 0., 900.);
       hi_ecal_vz_g1[i]->GetXaxis()->SetTitle("Z_{vertex} (cm)");
       hi_ecal_vz_g1[i]->GetYaxis()->SetTitle("Rate (MHz/3cm)");
       hi_ecal_vz_h1[i] = new TH1F(Form("hi_ecal_vz_h_layer1 %d",i+1), Form("hi_ecal_vz_h_layer1 %d",i+1),300, 0., 900.);
       hi_ecal_vz_h1[i]->GetXaxis()->SetTitle("Z_{vertex} (cm)");
       hi_ecal_vz_h1[i]->GetYaxis()->SetTitle("Rate (MHz/3cm)");
       hi_ecal_vz_n1[i] = new TH1F(Form("hi_ecal_vz_n_layer1 %d",i+1), Form("hi_ecal_vz_n_layer1 %d",i+1),300, 0., 900.);
       hi_ecal_vz_n1[i]->GetXaxis()->SetTitle("Z_{vertex} (cm)");
       hi_ecal_vz_n1[i]->GetYaxis()->SetTitle("Rate (MHz/3cm)");
       hi_ecal_origin_all1[i] = new TH2F(Form("hi_ecal_origin_all_layer1 %d",i+1), Form("hi_ecal_origin_all_layer1 %d",i+1),300, 0., 900.,200, 0., 600.);
       hi_ecal_origin_all1[i]->GetXaxis()->SetTitle("Z_{vertex} (cm)");
       hi_ecal_origin_all1[i]->GetYaxis()->SetTitle("R_{vertex} (cm)");
       hi_ecal_edep1[i] = new TH1F(Form("hi_ecal_edep_layer1 %d",i+1), Form("hi_ecal_edep_layer1 %d",i+1),200, 0.,10.);
       hi_ecal_edep1[i]->GetXaxis()->SetTitle("E(MeV)");
       hi_ecal_edep1[i]->GetYaxis()->SetTitle("Rate (MHz/50keV)");
       hi_ecal_occ2[i] = new TH2F(Form("hi_ecal_occ_layer2 %d",i+1), Form("hi_ecal_occ_layer2 %d",i+1),nstrip, 1.,nstrip*1.+1.,3,1.,4.);
       hi_ecal_occ2[i]->GetXaxis()->SetTitle("Strip");
       hi_ecal_occ2[i]->GetYaxis()->SetTitle("View");
       hi_ecal_occ_cut2[i] = new TH2F(Form("hi_ecal_occ_cut_layer2 %d",i+1), Form("hi_ecal_occ_cut_layer2 %d",i+1),nstrip, 1.,nstrip*1.+1.,3,1.,4.);
       hi_ecal_occ_cut2[i]->GetXaxis()->SetTitle("Strip");
       hi_ecal_occ_cut2[i]->GetYaxis()->SetTitle("View");
       hi_ecal_vz_all2[i] = new TH1F(Form("hi_ecal_vz_all_layer2 %d",i+1), Form("hi_ecal_vz_all_layer2 %d",i+1),300, 0., 900.);
       hi_ecal_vz_all2[i]->GetXaxis()->SetTitle("Z_{vertex} (cm)");
       hi_ecal_vz_all2[i]->GetYaxis()->SetTitle("Rate (MHz/3cm)");
       hi_ecal_vz_e2[i] = new TH1F(Form("hi_ecal_vz_e_layer2 %d",i+1), Form("hi_ecal_vz_e_layer2 %d",i+1), 300, 0., 900.);
       hi_ecal_vz_e2[i]->GetXaxis()->SetTitle("Z_{vertex} (cm)");
       hi_ecal_vz_e2[i]->GetYaxis()->SetTitle("Rate (MHz/3cm)");
       hi_ecal_vz_g2[i] = new TH1F(Form("hi_ecal_vz_g_layer2 %d",i+1), Form("hi_ecal_vz_g_layer2 %d",i+1),300, 0., 900.);
       hi_ecal_vz_g2[i]->GetXaxis()->SetTitle("Z_{vertex} (cm)");
       hi_ecal_vz_g2[i]->GetYaxis()->SetTitle("Rate (MHz/3cm)");
       hi_ecal_vz_h2[i] = new TH1F(Form("hi_ecal_vz_h_layer2 %d",i+1), Form("hi_ecal_vz_h_layer2 %d",i+1),300, 0., 900.);
       hi_ecal_vz_h2[i]->GetXaxis()->SetTitle("Z_{vertex} (cm)");
       hi_ecal_vz_h2[i]->GetYaxis()->SetTitle("Rate (MHz/3cm)");
       hi_ecal_vz_n2[i] = new TH1F(Form("hi_ecal_vz_n_layer2 %d",i+1), Form("hi_ecal_vz_n_layer2 %d",i+1),300, 0., 900.);
       hi_ecal_vz_n2[i]->GetXaxis()->SetTitle("Z_{vertex} (cm)");
       hi_ecal_vz_n2[i]->GetYaxis()->SetTitle("Rate (MHz/3cm)");
       hi_ecal_origin_all2[i] = new TH2F(Form("hi_ecal_origin_all_layer2 %d",i+1), Form("hi_ecal_origin_all_layer2 %d",i+1),300, 0., 900.,200, 0., 600.);
       hi_ecal_origin_all2[i]->GetXaxis()->SetTitle("Z_{vertex} (cm)");
       hi_ecal_origin_all2[i]->GetYaxis()->SetTitle("R_{vertex} (cm)");
       hi_ecal_edep2[i] = new TH1F(Form("hi_ecal_edep_layer2 %d",i+1), Form("hi_ecal_edep_layer2 %d",i+1),200, 0.,10.);
       hi_ecal_edep2[i]->GetXaxis()->SetTitle("E(MeV)");
       hi_ecal_edep2[i]->GetYaxis()->SetTitle("Rate (MHz/50keV)");
       hi_ecal_occ3[i] = new TH2F(Form("hi_ecal_occ_layer3 %d",i+1), Form("hi_ecal_occ_layer3 %d",i+1),nstrip, 1.,nstrip*1.+1.,3,1.,4.);
       hi_ecal_occ3[i]->GetXaxis()->SetTitle("Strip");
       hi_ecal_occ3[i]->GetYaxis()->SetTitle("View");
       hi_ecal_occ_cut3[i] = new TH2F(Form("hi_ecal_occ_cut_layer3 %d",i+1), Form("hi_ecal_occ_cut_layer3 %d",i+1),nstrip, 1.,nstrip*1.+1.,3,1.,4.);
       hi_ecal_occ_cut3[i]->GetXaxis()->SetTitle("Strip");
       hi_ecal_occ_cut3[i]->GetYaxis()->SetTitle("View");
       hi_ecal_vz_all3[i] = new TH1F(Form("hi_ecal_vz_all_layer3 %d",i+1), Form("hi_ecal_vz_all_layer3 %d",i+1),300, 0., 900.);
       hi_ecal_vz_all3[i]->GetXaxis()->SetTitle("Z_{vertex} (cm)");
       hi_ecal_vz_all3[i]->GetYaxis()->SetTitle("Rate (MHz/3cm)");
       hi_ecal_vz_e3[i] = new TH1F(Form("hi_ecal_vz_e_layer3 %d",i+1), Form("hi_ecal_vz_e_layer3 %d",i+1), 300, 0., 900.);
       hi_ecal_vz_e3[i]->GetXaxis()->SetTitle("Z_{vertex} (cm)");
       hi_ecal_vz_e3[i]->GetYaxis()->SetTitle("Rate (MHz/3cm)");
       hi_ecal_vz_g3[i] = new TH1F(Form("hi_ecal_vz_g_layer3 %d",i+1), Form("hi_ecal_vz_g_layer3 %d",i+1),300, 0., 900.);
       hi_ecal_vz_g3[i]->GetXaxis()->SetTitle("Z_{vertex} (cm)");
       hi_ecal_vz_g3[i]->GetYaxis()->SetTitle("Rate (MHz/3cm)");
       hi_ecal_vz_h3[i] = new TH1F(Form("hi_ecal_vz_h_layer3 %d",i+1), Form("hi_ecal_vz_h_layer3 %d",i+1),300, 0., 900.);
       hi_ecal_vz_h3[i]->GetXaxis()->SetTitle("Z_{vertex} (cm)");
       hi_ecal_vz_h3[i]->GetYaxis()->SetTitle("Rate (MHz/3cm)");
       hi_ecal_vz_n3[i] = new TH1F(Form("hi_ecal_vz_n_layer3 %d",i+1), Form("hi_ecal_vz_n_layer3 %d",i+1),300, 0., 900.);
       hi_ecal_vz_n3[i]->GetXaxis()->SetTitle("Z_{vertex} (cm)");
       hi_ecal_vz_n3[i]->GetYaxis()->SetTitle("Rate (MHz/3cm)");
       hi_ecal_origin_all3[i] = new TH2F(Form("hi_ecal_origin_all_layer3 %d",i+1), Form("hi_ecal_origin_all_layer3 %d",i+1),300, 0., 900.,200, 0., 600.);
       hi_ecal_origin_all3[i]->GetXaxis()->SetTitle("Z_{vertex} (cm)");
       hi_ecal_origin_all3[i]->GetYaxis()->SetTitle("R_{vertex} (cm)");
       hi_ecal_edep3[i] = new TH1F(Form("hi_ecal_edep_layer3 %d",i+1), Form("hi_ecal_edep_layer3 %d",i+1),200, 0.,10.);
       hi_ecal_edep3[i]->GetXaxis()->SetTitle("E(MeV)");
       hi_ecal_edep3[i]->GetYaxis()->SetTitle("Rate (MHz/50keV)");
   }
   TH2F *hi_ecal_y_vs_x [2];
   for(int i=0; i<2; i++){
   	hi_ecal_y_vs_x[i] = new TH2F(Form("ecal_y_vs_x%i",i+1),Form("ecal_y_vs_x%1",i+1),200, -5000.,5000., 200, -5000.,5000.);
   	hi_ecal_y_vs_x[i]->GetXaxis()->SetTitle("x(mm)");
   	hi_ecal_y_vs_x[i]->GetYaxis()->SetTitle("y(mm)");
   }
// in standard run mode (1 electron at a time)
//   float lumi=nentries*5.*0.07*0.602/1000; // (mbarn^-1 or 10^27 cm^-2)
// in luminosity mode (50*2300 or 118600 electrons in 250ns window)
   float levents=124000;
   int   nsum=124000./levents;
   float norm=124000/levents/(nentries);
   float lumi=nentries*250/10*levents/124000; // (mbarn^-1 or 10^27 cm^-2)
   float time=250/norm/1000;

   cout << "number of events to integrate = " << nsum  << endl;
   cout << "normalization factor          = " << norm  << endl;
   cout << "Run time                      = " << time  << " us" << endl;


   float Ethr=0.4;            // energy thresholds in MeV
   int nint=0;
   int ngoodentries=0;
   for (Long64_t jentry=0; jentry < nentries; jentry++) {
  
//       pcal_sector->clear();
//       pcal_layer->clear();
//       pcal_view->clear();
//       pcal_strip->clear();
//       pcal_ADC->clear();
//       pcal_TDC->clear();
//       pcal_E->clear();
//       pcal_Edep->clear();
//       pcal_t->clear();
//       pcal_pid->clear();
//       pcal_x->clear();
//       pcal_y->clear();
//       pcal_z->clear();
//       pcal_lx->clear();
//       pcal_ly->clear();
//       pcal_lz->clear(); 
//       pcal_vx->clear();
//       pcal_vy->clear();
//       pcal_vz->clear();
     
       ec_sector->clear();
       ec_layer->clear();
       ec_view->clear();
       ec_strip->clear();
       ec_ADC->clear();
       ec_TDC->clear();
       ec_E->clear();
       ec_Edep->clear();
       ec_t->clear();
       ec_pid->clear();
       ec_x->clear();
       ec_y->clear();
       ec_z->clear();
       ec_lx->clear();
       ec_ly->clear();
       ec_lz->clear(); 
       ec_vx->clear();
       ec_vy->clear();
       ec_vz->clear();
     
       int nb = ec->GetEntry(jentry); 
       ec->GetEntry(jentry); 
       ngoodentries++;
       int pi = 3.14159265359;
       double phi;
       if(int(jentry/1000)*1000==jentry) cout << "Analyzed " << jentry << " events of " << nentries << endl;
	
//       npcalhit=pcal_pid->size();
       nechit=ec_pid->size();

//       for(int i=0; i<npcalhit; i++) {
//           phi = atan2((*ec_vy)[i],(*ec_vx)[i])*180/pi;
//           phi += (phi<0) ? 360:0;
//           int nsect = floor((phi-30)/60) + 2;
//           if(nsect == 7) nsect = 1;
//  	   hi_ecal_occ[0]->Fill((*pcal_strip)[i]*1.,(*pcal_view)[i]*1.);
//	   hi_ecal_vz_all[0]->Fill((*pcal_vz)[i]/10.);
//	   if(abs((*pcal_pid)[i])==11) hi_ecal_vz_e[0]->Fill((*pcal_vz)[i]/10.);
//	   else if((*pcal_pid)[i]==22) hi_ecal_vz_g[0]->Fill((*pcal_vz)[i]/10.);
//	   else                       hi_ecal_vz_h[0]->Fill((*pcal_vz)[i]/10.);
//	   if((*pcal_pid)[i]==2112)    hi_ecal_vz_n[0]->Fill((*pcal_vz)[i]/10.);
//	   hi_ecal_origin_all[0]->Fill((*pcal_vz)[i]/10.,sqrt((*pcal_vx)[i]*(*pcal_vx)[i]/100.+(*pcal_vy)[i]*(*pcal_vy)[i]/100.));
//	   hi_ecal_edep[0]->Fill((*pcal_Edep)[i]);
//	   if((*pcal_Edep)[i]>Ethr) {
//	     hi_ecal_occ_cut[0]->Fill((*pcal_strip)[i]*1.,(*pcal_view)[i]*1.);
//	   }
//       }
       for(int i=0; i<nechit; i++) {
	   int istack=(*ec_layer)[i];
           phi = atan2((*ec_vy)[i],(*ec_vx)[i])*180/pi;
	   phi += (phi<0) ? 360:0;
	   int nsect = floor((phi-30)/60) + 2;
	   if(nsect == 7) nsect = 1;
           hi_ecal_y_vs_x[istack-1]->Fill((*ec_vx)[i],(*ec_vy)[i]);
           if(istack == 1){
  	      hi_ecal_occ1[nsect-1]->Fill((*ec_strip)[i]*1.,(*ec_view)[i]*1.);
	      hi_ecal_vz_all1[nsect-1]->Fill((*ec_vz)[i]/10.);
	      if(abs((*ec_pid)[i])==11) hi_ecal_vz_e1[nsect-1]->Fill((*ec_vz)[i]/10.);
	      else if((*ec_pid)[i]==22) hi_ecal_vz_g1[nsect-1]->Fill((*ec_vz)[i]/10.);
	      else                       hi_ecal_vz_h1[nsect-1]->Fill((*ec_vz)[i]/10.);
	      if((*ec_pid)[i]==2112)    hi_ecal_vz_n1[nsect-1]->Fill((*ec_vz)[i]/10.);
	      hi_ecal_origin_all1[nsect-1]->Fill((*ec_vz)[i]/10.,sqrt((*ec_vx)[i]*(*ec_vx)[i]/100.+(*ec_vy)[i]*(*ec_vy)[i]/100.));
//	      hi_ecal_edep1[nsect-1]->Fill((*pcal_Edep)[i]);
	      if((*ec_Edep)[i]>Ethr) {
	      hi_ecal_occ_cut1[nsect-1]->Fill((*ec_strip)[i]*1.,(*ec_view)[i]*1.);
               }
           }
           if(istack == 2){
              hi_ecal_occ2[nsect-1]->Fill((*ec_strip)[i]*1.,(*ec_view)[i]*1.);
              hi_ecal_vz_all2[nsect-1]->Fill((*ec_vz)[i]/10.);
              if(abs((*ec_pid)[i])==11) hi_ecal_vz_e2[nsect-1]->Fill((*ec_vz)[i]/10.);
              else if((*ec_pid)[i]==22) hi_ecal_vz_g2[nsect-1]->Fill((*ec_vz)[i]/10.);
              else                       hi_ecal_vz_h2[nsect-1]->Fill((*ec_vz)[i]/10.);
              if((*ec_pid)[i]==2112)    hi_ecal_vz_n2[nsect-1]->Fill((*ec_vz)[i]/10.);
              hi_ecal_origin_all2[nsect-1]->Fill((*ec_vz)[i]/10.,sqrt((*ec_vx)[i]*(*ec_vx)[i]/100.+(*ec_vy)[i]*(*ec_vy)[i]/100.));
//            hi_ecal_edep2[nsect-1]->Fill((*pcal_Edep)[i]);
              if((*ec_Edep)[i]>Ethr) {
              hi_ecal_occ_cut2[nsect-1]->Fill((*ec_strip)[i]*1.,(*ec_view)[i]*1.);
               }
            }
           if(istack == 3){
              hi_ecal_occ3[nsect-1]->Fill((*ec_strip)[i]*1.,(*ec_view)[i]*1.);
              hi_ecal_vz_all3[nsect-1]->Fill((*ec_vz)[i]/10.);
              if(abs((*ec_pid)[i])==11) hi_ecal_vz_e3[nsect-1]->Fill((*ec_vz)[i]/10.);
              else if((*ec_pid)[i]==22) hi_ecal_vz_g3[nsect-1]->Fill((*ec_vz)[i]/10.);
              else                       hi_ecal_vz_h3[nsect-1]->Fill((*ec_vz)[i]/10.);
              if((*ec_pid)[i]==2112)    hi_ecal_vz_n3[nsect-1]->Fill((*ec_vz)[i]/10.);
              hi_ecal_origin_all3[nsect-1]->Fill((*ec_vz)[i]/10.,sqrt((*ec_vx)[i]*(*ec_vx)[i]/100.+(*ec_vy)[i]*(*ec_vy)[i]/100.));
//            hi_ecal_edep3[nsect-1]->Fill((*pcal_Edep)[i]);
              if((*ec_Edep)[i]>Ethr) {
              hi_ecal_occ_cut3[nsect-1]->Fill((*ec_strip)[i]*1.,(*ec_view)[i]*1.);

	     }
          }
       }
   }
   
// in standard run mode (1 electron at a time)
//   float lumi=nentries*5.*0.07*0.602/1000; // (mbarn^-1 or 10^27 cm^-2)
// in luminosity mode (50*2300 or 118600 electrons in 250ns window)
   norm=124000/levents/ngoodentries;
   lumi=ngoodentries*250/10*levents/124000; // (mbarn^-1 or 10^27 cm^-2)
   time=250/norm/1000;

   cout << "number of events to integrate = " << nsum  << endl;
   cout << "normalization factor          = " << norm  << endl;
   cout << "Run time                      = " << time  << " us" << endl;

   // normalizing rate histogram to 10^35 luminosity
   for(int i=0; i<6; i++) {
       hi_ecal_occ1[i]->Scale(1000/time/6);
       hi_ecal_occ_cut1[i]->Scale(1000/time/6);
       hi_ecal_vz_all1[i]->Scale(1/time);
       hi_ecal_vz_e1[i]->Scale(1/time);
       hi_ecal_vz_g1[i]->Scale(1/time);
       hi_ecal_vz_h1[i]->Scale(1/time);
       hi_ecal_vz_n1[i]->Scale(1/time);
       hi_ecal_origin_all1[i]->Scale(1/time);
       hi_ecal_edep1[i]->Scale(1/time);
       hi_ecal_occ2[i]->Scale(1000/time/6);
       hi_ecal_occ_cut2[i]->Scale(1000/time/6);
       hi_ecal_vz_all2[i]->Scale(1/time);
       hi_ecal_vz_e2[i]->Scale(1/time);
       hi_ecal_vz_g2[i]->Scale(1/time);
       hi_ecal_vz_h2[i]->Scale(1/time);
       hi_ecal_vz_n2[i]->Scale(1/time);
       hi_ecal_origin_all2[i]->Scale(1/time);
       hi_ecal_edep2[i]->Scale(1/time);
       hi_ecal_occ3[i]->Scale(1000/time/6);
       hi_ecal_occ_cut3[i]->Scale(1000/time/6);
       hi_ecal_vz_all3[i]->Scale(1/time);
       hi_ecal_vz_e3[i]->Scale(1/time);
       hi_ecal_vz_g3[i]->Scale(1/time);
       hi_ecal_vz_h3[i]->Scale(1/time);
       hi_ecal_vz_n3[i]->Scale(1/time);
       hi_ecal_origin_all3[i]->Scale(1/time);
       hi_ecal_edep3[i]->Scale(1/time);
   }
   

   TCanvas *c_occ1=new TCanvas("c_occ1","Occupancy",750,1000);
   c_occ1->Divide(3,2);
   FILE *fp = fopen("ecal_occupancy.txt","w");
   for(int i=0; i<6; i++) {
       c_occ1->cd(i+1);
       hi_ecal_occ1[i]->Draw("COLZ");
       for(int iv=0; iv<hi_ecal_occ1[i]->GetNbinsY(); iv++) {
  	   for(int is=0; is<hi_ecal_occ1[i]->GetNbinsX(); is++) {
  	     int layer = i*3+iv+1;
  	     int strip = is +1;
  	     float value = hi_ecal_occ1[i]->GetBinContent(is+1,iv+1);
  	     fprintf(fp,"%d\t%d\t%8.1f\n",layer,strip,value);
  	   }
         }
   }
   fclose(fp);
   c_occ1->Print("ecal_occupancy.pdf(");

   TCanvas *c_occ2=new TCanvas("c_occ2","Occupancy",750,1000);
   c_occ2->Divide(3,2);
   fp = fopen("ecal_occupancy.txt","w");
   for(int i=0; i<6; i++) {
       c_occ2->cd(i+1);
       hi_ecal_occ2[i]->Draw("COLZ");
       for(int iv=0; iv<hi_ecal_occ2[i]->GetNbinsY(); iv++) {
           for(int is=0; is<hi_ecal_occ2[i]->GetNbinsX(); is++) {
             int layer = i*3+iv+1;
             int strip = is +1;
             float value = hi_ecal_occ2[i]->GetBinContent(is+1,iv+1);
             fprintf(fp,"%d\t%d\t%8.1f\n",layer,strip,value);
           }
       }
   }
   fclose(fp);
   c_occ2->Print("ecal_occupancy.pdf(");

 //  TCanvas *c_occ3=new TCanvas("c_occ3","Occupancy",750,1000);
 //  c_occ3->Divide(3,2);
//   fp = fopen("ecal_occupancy.txt","w");
//   for(int i=0; i<6; i++) {
//       c_occ3->cd(i+1);
//       hi_ecal_occ3[i]->Draw("COLZ");
//       for(int iv=0; iv<hi_ecal_occ3[i]->GetNbinsY(); iv++) {
//           for(int is=0; is<hi_ecal_occ3[i]->GetNbinsX(); is++) {
//             int layer = i*3+iv+1;
//             int strip = is +1;
//             float value = hi_ecal_occ3[i]->GetBinContent(is+1,iv+1);
//             fprintf(fp,"%d\t%d\t%8.1f\n",layer,strip,value);
//           }
//       }
//   }
//   fclose(fp);
//   c_occ3->Print("ecal_occupancy.pdf");


   TCanvas *c_occ_cut1=new TCanvas("c_occ_cut1","Occupancy_Cuts",750,1000);
   c_occ_cut1->Divide(3,2);
   fp = fopen("ecal_occupancy_cut.txt","w");
   for(int i=0; i<6; i++) {
       c_occ_cut1->cd(i+1);
       hi_ecal_occ_cut1[i]->Draw("COLZ");
       for(int iv=0; iv<hi_ecal_occ_cut1[i]->GetNbinsY(); iv++) {
	   for(int is=0; is<hi_ecal_occ_cut1[i]->GetNbinsX(); is++) {
	     int layer = i*3+iv+1;
	     int strip = is +1;
	     float value = hi_ecal_occ_cut1[i]->GetBinContent(is+1,iv+1);
	     fprintf(fp,"%d\t%d\t%8.1f\n",layer,strip,value);
	   }
       }
   }
   fclose(fp);
   c_occ_cut1->Print("ecal_occupancy.pdf");

   TCanvas *c_occ_cut2=new TCanvas("c_occ_cut2","Occupancy_Cuts",750,1000);
   c_occ_cut2->Divide(3,2);
   fp = fopen("ecal_occupancy_cut.txt","w");
   for(int i=0; i<6; i++) {
       c_occ_cut2->cd(i+1);
       hi_ecal_occ_cut2[i]->Draw("COLZ");
       for(int iv=0; iv<hi_ecal_occ_cut2[i]->GetNbinsY(); iv++) {
           for(int is=0; is<hi_ecal_occ_cut2[i]->GetNbinsX(); is++) {
             int layer = i*3+iv+1;
             int strip = is +1;
             float value = hi_ecal_occ_cut2[i]->GetBinContent(is+1,iv+1);
             fprintf(fp,"%d\t%d\t%8.1f\n",layer,strip,value);
           }
       }
   }
   fclose(fp);
   c_occ_cut2->Print("ecal_occupancy.pdf");

//   TCanvas *c_occ_cut3=new TCanvas("c_occ_cut3","Occupancy_Cuts",750,1000);
//   c_occ_cut3->Divide(3,2);
//   fp = fopen("ecal_occupancy_cut.txt","w");
//   for(int i=0; i<6; i++) {
//       c_occ_cut3->cd(i);
//       hi_ecal_occ_cut3[i]->Draw("COLZ");
//       for(int iv=0; iv<hi_ecal_occ_cut3[i]->GetNbinsY(); iv++) {
//           for(int is=0; is<hi_ecal_occ_cut3[i]->GetNbinsX(); is++) {
//             int layer = i*3+iv+1;
//             int strip = is +1;
//             float value = hi_ecal_occ_cut3[i]->GetBinContent(is+1,iv+1);
//             fprintf(fp,"%d\t%d\t%8.1f\n",layer,strip,value);
//           }
//       }
//   }
//   fclose(fp);
//   c_occ_cut3->Print("ecal_occupancy.pdf");

   TCanvas *c_origin1=new TCanvas("c_origin1","Origin",750,1000);
   c_origin1->Divide(3,2);
   for(int i=0; i<6; i++) {
       c_origin1->cd(i+1);
       gPad->SetLogz();
       hi_ecal_origin_all1[i]->Draw("COLZ");
   }
   c_origin1->Print("ecal_occupancy.pdf");

   TCanvas *c_origin2=new TCanvas("c_origin2","Origin",750,1000);
   c_origin2->Divide(3,2);
   for(int i=0; i<6; i++) {
       c_origin2->cd(i+1);
       gPad->SetLogz();
       hi_ecal_origin_all2[i]->Draw("COLZ");
   }
   c_origin2->Print("ecal_occupancy.pdf");

//   TCanvas *c_origin3=new TCanvas("c_origin3","Origin",750,1000);
//   c_origin3->Divide(3,2);
//   for(int i=0; i<6; i++) {
//       c_origin3->cd(i+1);
//       gPad->SetLogz();
//       hi_ecal_origin_all1[i]->Draw("COLZ");
//   }
//   c_origin3->Print("ecal_occupancy.pdf");

   TCanvas *c_vz1=new TCanvas("c_vz1","VZ",750,1000);
   c_vz1->Divide(3,2);
   c_vz1->cd(1);
   hi_ecal_vz_all1[0]->SetMinimum(0.001);
   hi_ecal_vz_all1[0]->Draw("H");
   hi_ecal_vz_e1[0]->SetLineColor(2);
   hi_ecal_vz_e1[0]->Draw("SAME");
   hi_ecal_vz_g1[0]->SetLineColor(4);
   hi_ecal_vz_g1[0]->Draw("SAME");
   hi_ecal_vz_h1[0]->SetLineColor(kGreen);
   hi_ecal_vz_h1[0]->Draw("SAME");
   hi_ecal_vz_n1[0]->SetLineColor(kGreen+2);
   hi_ecal_vz_n1[0]->Draw("SAME");

    TLegend *leg1 = new TLegend(0.7,0.75,0.96,0.96);
     leg1->SetTextSize(.04);
     leg1->AddEntry(hi_ecal_vz_all1[0],"All","l");
     leg1->AddEntry(hi_ecal_vz_e1[0],"electrons","l");
     leg1->AddEntry(hi_ecal_vz_g1[0],"photons","l");
     leg1->AddEntry(hi_ecal_vz_n1[0],"neutron","l");
     leg1->AddEntry(hi_ecal_vz_h1[0],"hadron","l");
     leg1->Draw();

   c_vz1->cd(2);
   hi_ecal_vz_all1[1]->SetMinimum(0.001);
   hi_ecal_vz_all1[1]->Draw("H");
   hi_ecal_vz_e1[1]->SetLineColor(2);
   hi_ecal_vz_e1[1]->Draw("SAME");
   hi_ecal_vz_g1[1]->SetLineColor(4);
   hi_ecal_vz_g1[1]->Draw("SAME");
   hi_ecal_vz_h1[1]->SetLineColor(kGreen);
   hi_ecal_vz_h1[1]->Draw("SAME");
   hi_ecal_vz_n1[1]->SetLineColor(kGreen+2);
   hi_ecal_vz_n1[1]->Draw("SAME");

    TLegend *leg2 = new TLegend(0.7,0.75,0.96,0.96);
     leg2->SetTextSize(.04);
     leg2->AddEntry(hi_ecal_vz_all1[1],"All","l");
     leg2->AddEntry(hi_ecal_vz_e1[1],"electrons","l");
     leg2->AddEntry(hi_ecal_vz_g1[1],"photons","l");
     leg2->AddEntry(hi_ecal_vz_n1[1],"neutron","l");
     leg2->AddEntry(hi_ecal_vz_h1[1],"hadron","l");
     leg2->Draw();

   c_vz1->cd(3);
   hi_ecal_vz_all1[2]->SetMinimum(0.001);
   hi_ecal_vz_all1[2]->Draw("H");
   hi_ecal_vz_e1[2]->SetLineColor(2);
   hi_ecal_vz_e1[2]->Draw("SAME");
   hi_ecal_vz_g1[2]->SetLineColor(4);
   hi_ecal_vz_g1[2]->Draw("SAME");
   hi_ecal_vz_h1[2]->SetLineColor(kGreen);
   hi_ecal_vz_h1[2]->Draw("SAME");
   hi_ecal_vz_n1[2]->SetLineColor(kGreen+2);
   hi_ecal_vz_n1[2]->Draw("SAME");

    TLegend *leg3 = new TLegend(0.7,0.75,0.96,0.96);
     leg3->SetTextSize(.04);
     leg3->AddEntry(hi_ecal_vz_all1[2],"All","l");
     leg3->AddEntry(hi_ecal_vz_e1[2],"electrons","l");
     leg3->AddEntry(hi_ecal_vz_g1[2],"photons","l");
     leg3->AddEntry(hi_ecal_vz_n1[2],"neutron","l");
     leg3->AddEntry(hi_ecal_vz_h1[2],"hadron","l");
     leg3->Draw();

   c_vz1->cd(4);
   hi_ecal_vz_all1[3]->SetMinimum(0.001);
   hi_ecal_vz_all1[3]->Draw("H");
   hi_ecal_vz_e1[3]->SetLineColor(2);
   hi_ecal_vz_e1[3]->Draw("SAME");
   hi_ecal_vz_g1[3]->SetLineColor(4);
   hi_ecal_vz_g1[3]->Draw("SAME");
   hi_ecal_vz_h1[3]->SetLineColor(kGreen);
   hi_ecal_vz_h1[3]->Draw("SAME");
   hi_ecal_vz_n1[3]->SetLineColor(kGreen+2);
   hi_ecal_vz_n1[3]->Draw("SAME");

    TLegend *leg4 = new TLegend(0.7,0.75,0.96,0.96);
     leg4->SetTextSize(.04);
     leg4->AddEntry(hi_ecal_vz_all1[3],"All","l");
     leg4->AddEntry(hi_ecal_vz_e1[3],"electrons","l");
     leg4->AddEntry(hi_ecal_vz_g1[3],"photons","l");
     leg4->AddEntry(hi_ecal_vz_n1[3],"neutron","l");
     leg4->AddEntry(hi_ecal_vz_h1[3],"hadron","l");
     leg4->Draw();

   c_vz1->cd(5);
   hi_ecal_vz_all1[4]->SetMinimum(0.001);
   hi_ecal_vz_all1[4]->Draw("H");
   hi_ecal_vz_e1[4]->SetLineColor(2);
   hi_ecal_vz_e1[4]->Draw("SAME");
   hi_ecal_vz_g1[4]->SetLineColor(4);
   hi_ecal_vz_g1[4]->Draw("SAME");
   hi_ecal_vz_h1[4]->SetLineColor(kGreen);
   hi_ecal_vz_h1[4]->Draw("SAME");
   hi_ecal_vz_n1[4]->SetLineColor(kGreen+2);
   hi_ecal_vz_n1[4]->Draw("SAME");

    TLegend *leg5 = new TLegend(0.7,0.75,0.96,0.96);
     leg5->SetTextSize(.04);
     leg5->AddEntry(hi_ecal_vz_all1[4],"All","l");
     leg5->AddEntry(hi_ecal_vz_e1[4],"electrons","l");
     leg5->AddEntry(hi_ecal_vz_g1[4],"photons","l");
     leg5->AddEntry(hi_ecal_vz_n1[4],"neutron","l");
     leg5->AddEntry(hi_ecal_vz_h1[4],"hadron","l");
     leg5->Draw();

   c_vz1->cd(6);
   hi_ecal_vz_all1[5]->SetMinimum(0.001);
   hi_ecal_vz_all1[5]->Draw("H");
   hi_ecal_vz_e1[5]->SetLineColor(2);
   hi_ecal_vz_e1[5]->Draw("SAME");
   hi_ecal_vz_g1[5]->SetLineColor(4);
   hi_ecal_vz_g1[5]->Draw("SAME");
   hi_ecal_vz_h1[5]->SetLineColor(kGreen);
   hi_ecal_vz_h1[5]->Draw("SAME");
   hi_ecal_vz_n1[5]->SetLineColor(kGreen+2);
   hi_ecal_vz_n1[5]->Draw("SAME");

    TLegend *leg6 = new TLegend(0.7,0.75,0.96,0.96);
     leg6->SetTextSize(.04);
     leg6->AddEntry(hi_ecal_vz_all1[5],"All","l");
     leg6->AddEntry(hi_ecal_vz_e1[5],"electrons","l");
     leg6->AddEntry(hi_ecal_vz_g1[5],"photons","l");
     leg6->AddEntry(hi_ecal_vz_n1[5],"neutron","l");
     leg6->AddEntry(hi_ecal_vz_h1[5],"hadron","l");
     leg6->Draw();
   c_vz1->Print("ecal_occupancy.pdf");

   TCanvas *c_vz2=new TCanvas("c_vz2","VZ",750,1000);
   c_vz2->Divide(3,2);
   c_vz2->cd(1);
   hi_ecal_vz_all2[0]->SetMinimum(0.001);
   hi_ecal_vz_all2[0]->Draw("H");
   hi_ecal_vz_e2[0]->SetLineColor(2);
   hi_ecal_vz_e2[0]->Draw("SAME");
   hi_ecal_vz_g2[0]->SetLineColor(4);
   hi_ecal_vz_g2[0]->Draw("SAME");
   hi_ecal_vz_h2[0]->SetLineColor(kGreen);
   hi_ecal_vz_h2[0]->Draw("SAME");
   hi_ecal_vz_n2[0]->SetLineColor(kGreen+2);
   hi_ecal_vz_n2[0]->Draw("SAME");

    TLegend *leg21 = new TLegend(0.7,0.75,0.96,0.96);
     leg21->SetTextSize(.04);
     leg21->AddEntry(hi_ecal_vz_all2[0],"All","l");
     leg21->AddEntry(hi_ecal_vz_e2[0],"electrons","l");
     leg21->AddEntry(hi_ecal_vz_g2[0],"photons","l");
     leg21->AddEntry(hi_ecal_vz_n2[0],"neutron","l");
     leg21->AddEntry(hi_ecal_vz_h2[0],"hadron","l");
     leg21->Draw();

   c_vz2->cd(2);
   hi_ecal_vz_all2[1]->SetMinimum(0.001);
   hi_ecal_vz_all2[1]->Draw("H");
   hi_ecal_vz_e2[1]->SetLineColor(2);
   hi_ecal_vz_e2[1]->Draw("SAME");
   hi_ecal_vz_g2[1]->SetLineColor(4);
   hi_ecal_vz_g2[1]->Draw("SAME");
   hi_ecal_vz_h2[1]->SetLineColor(kGreen);
   hi_ecal_vz_h2[1]->Draw("SAME");
   hi_ecal_vz_n2[1]->SetLineColor(kGreen+2);
   hi_ecal_vz_n2[1]->Draw("SAME");

    TLegend *leg22 = new TLegend(0.7,0.75,0.96,0.96);
     leg22->SetTextSize(.04);
     leg22->AddEntry(hi_ecal_vz_all2[1],"All","l");
     leg22->AddEntry(hi_ecal_vz_e2[1],"electrons","l");
     leg22->AddEntry(hi_ecal_vz_g2[1],"photons","l");
     leg22->AddEntry(hi_ecal_vz_n2[1],"neutron","l");
     leg22->AddEntry(hi_ecal_vz_h2[1],"hadron","l");
     leg22->Draw();

   c_vz2->cd(3);
   hi_ecal_vz_all2[2]->SetMinimum(0.001);
   hi_ecal_vz_all2[2]->Draw("H");
   hi_ecal_vz_e2[2]->SetLineColor(2);
   hi_ecal_vz_e2[2]->Draw("SAME");
   hi_ecal_vz_g2[2]->SetLineColor(4);
   hi_ecal_vz_g2[2]->Draw("SAME");
   hi_ecal_vz_h2[2]->SetLineColor(kGreen);
   hi_ecal_vz_h2[2]->Draw("SAME");
   hi_ecal_vz_n2[2]->SetLineColor(kGreen+2);
   hi_ecal_vz_n2[2]->Draw("SAME");

    TLegend *leg23 = new TLegend(0.7,0.75,0.96,0.96);
     leg23->SetTextSize(.04);
     leg23->AddEntry(hi_ecal_vz_all2[2],"All","l");
     leg23->AddEntry(hi_ecal_vz_e2[2],"electrons","l");
     leg23->AddEntry(hi_ecal_vz_g2[2],"photons","l");
     leg23->AddEntry(hi_ecal_vz_n2[2],"neutron","l");
     leg23->AddEntry(hi_ecal_vz_h2[2],"hadron","l");
     leg23->Draw();

   c_vz2->cd(4);
   hi_ecal_vz_all2[3]->SetMinimum(0.001);
   hi_ecal_vz_all2[3]->Draw("H");
   hi_ecal_vz_e2[3]->SetLineColor(2);
   hi_ecal_vz_e2[3]->Draw("SAME");
   hi_ecal_vz_g2[3]->SetLineColor(4);
   hi_ecal_vz_g2[3]->Draw("SAME");
   hi_ecal_vz_h2[3]->SetLineColor(kGreen);
   hi_ecal_vz_h2[3]->Draw("SAME");
   hi_ecal_vz_n2[3]->SetLineColor(kGreen+2);
   hi_ecal_vz_n2[3]->Draw("SAME");

    TLegend *leg24 = new TLegend(0.7,0.75,0.96,0.96);
     leg24->SetTextSize(.04);
     leg24->AddEntry(hi_ecal_vz_all2[3],"All","l");
     leg24->AddEntry(hi_ecal_vz_e2[3],"electrons","l");
     leg24->AddEntry(hi_ecal_vz_g2[3],"photons","l");
     leg24->AddEntry(hi_ecal_vz_n2[3],"neutron","l");
     leg24->AddEntry(hi_ecal_vz_h2[3],"hadron","l");
     leg24->Draw();

   c_vz2->cd(5);
   hi_ecal_vz_all2[4]->SetMinimum(0.001);
   hi_ecal_vz_all2[4]->Draw("H");
   hi_ecal_vz_e2[4]->SetLineColor(2);
   hi_ecal_vz_e2[4]->Draw("SAME");
   hi_ecal_vz_g2[4]->SetLineColor(4);
   hi_ecal_vz_g2[4]->Draw("SAME");
   hi_ecal_vz_h2[4]->SetLineColor(kGreen);
   hi_ecal_vz_h2[4]->Draw("SAME");
   hi_ecal_vz_n2[4]->SetLineColor(kGreen+2);
   hi_ecal_vz_n2[4]->Draw("SAME");

    TLegend *leg25 = new TLegend(0.7,0.75,0.96,0.96);
     leg25->SetTextSize(.04);
     leg25->AddEntry(hi_ecal_vz_all2[4],"All","l");
     leg25->AddEntry(hi_ecal_vz_e2[4],"electrons","l");
     leg25->AddEntry(hi_ecal_vz_g2[4],"photons","l");
     leg25->AddEntry(hi_ecal_vz_n2[4],"neutron","l");
     leg25->AddEntry(hi_ecal_vz_h2[4],"hadron","l");
     leg25->Draw();

   c_vz2->cd(6);
   hi_ecal_vz_all2[5]->SetMinimum(0.001);
   hi_ecal_vz_all2[5]->Draw("H");
   hi_ecal_vz_e2[5]->SetLineColor(2);
   hi_ecal_vz_e2[5]->Draw("SAME");
   hi_ecal_vz_g2[5]->SetLineColor(4);
   hi_ecal_vz_g2[5]->Draw("SAME");
   hi_ecal_vz_h2[5]->SetLineColor(kGreen);
   hi_ecal_vz_h2[5]->Draw("SAME");
   hi_ecal_vz_n2[5]->SetLineColor(kGreen+2);
   hi_ecal_vz_n2[5]->Draw("SAME");

    TLegend *leg26 = new TLegend(0.7,0.75,0.96,0.96);
     leg26->SetTextSize(.04);
     leg26->AddEntry(hi_ecal_vz_all2[5],"All","l");
     leg26->AddEntry(hi_ecal_vz_e2[5],"electrons","l");
     leg26->AddEntry(hi_ecal_vz_g2[5],"photons","l");
     leg26->AddEntry(hi_ecal_vz_n2[5],"neutron","l");
     leg26->AddEntry(hi_ecal_vz_h2[5],"hadron","l");
     leg26->Draw();
   c_vz2->Print("ecal_occupancy.pdf");
   
   TCanvas *c_vx_vy=new TCanvas("c_vx_vy","VX_VY",750,1000);
   c_vx_vy->Divide(1,2);
   c_vx_vy->cd(1);
   gPad->SetLogz();
   hi_ecal_y_vs_x[0]->Draw("COLZ");
   c_vx_vy->cd(2);
   gPad->SetLogz();
   hi_ecal_y_vs_x[1]->Draw("COLZ");
   c_vx_vy->Print("ecal_occupancy.pdf)");
 //  TCanvas *c_vz3=new TCanvas("c_vz3","VZ",750,1000);
 //  c_vz3->Divide(3,2);
 //  c_vz3->cd(1);
 //  hi_ecal_vz_all3[0]->SetMinimum(0.001);
 //  hi_ecal_vz_all3[0]->Draw("H");
 //  hi_ecal_vz_e3[0]->SetLineColor(2);
 //  hi_ecal_vz_e3[0]->Draw("SAME");
 //  hi_ecal_vz_g3[0]->SetLineColor(4);
 //  hi_ecal_vz_g3[0]->Draw("SAME");
 //  hi_ecal_vz_h3[0]->SetLineColor(kGreen);
 //  hi_ecal_vz_h3[0]->Draw("SAME");
 //  hi_ecal_vz_n3[0]->SetLineColor(kGreen+2);
 //  hi_ecal_vz_n3[0]->Draw("SAME");

 //   TLegend *leg31 = new TLegend(0.7,0.75,0.96,0.96);
 //    leg31->SetTextSize(.04);
 //    leg31->AddEntry(hi_ecal_vz_all3[0],"All","l");
 //    leg31->AddEntry(hi_ecal_vz_e3[0],"electrons","l");
 //    leg31->AddEntry(hi_ecal_vz_g3[0],"photons","l");
 //    leg31->AddEntry(hi_ecal_vz_n3[0],"neutron","l");
 //    leg31->AddEntry(hi_ecal_vz_h3[0],"hadron","l");
 //    leg31->Draw();

 //  c_vz3->cd(2);
 //  hi_ecal_vz_all3[1]->SetMinimum(0.001);
 //  hi_ecal_vz_all3[1]->Draw("H");
 //  hi_ecal_vz_e3[1]->SetLineColor(2);
 //  hi_ecal_vz_e3[1]->Draw("SAME");
 //  hi_ecal_vz_g3[1]->SetLineColor(4);
 //  hi_ecal_vz_g3[1]->Draw("SAME");
 //  hi_ecal_vz_h3[1]->SetLineColor(kGreen);
 //  hi_ecal_vz_h3[1]->Draw("SAME");
 //  hi_ecal_vz_n3[1]->SetLineColor(kGreen+2);
 //  hi_ecal_vz_n3[1]->Draw("SAME");

 //   TLegend *leg32 = new TLegend(0.7,0.75,0.96,0.96);
 //    leg32->SetTextSize(.04);
 //    leg32->AddEntry(hi_ecal_vz_all3[1],"All","l");
 //    leg32->AddEntry(hi_ecal_vz_e3[1],"electrons","l");
 //    leg32->AddEntry(hi_ecal_vz_g3[1],"photons","l");
 //    leg32->AddEntry(hi_ecal_vz_n3[1],"neutron","l");
 //    leg32->AddEntry(hi_ecal_vz_h3[1],"hadron","l");
 //   leg32->Draw();

 //  c_vz3->cd(3);
 //  hi_ecal_vz_all3[2]->SetMinimum(0.001);
 //  hi_ecal_vz_all3[2]->Draw("H");
 //  hi_ecal_vz_e3[2]->SetLineColor(2);
 //  hi_ecal_vz_e3[2]->Draw("SAME");
 //  hi_ecal_vz_g3[2]->SetLineColor(4);
 //  hi_ecal_vz_g3[2]->Draw("SAME");
 //  hi_ecal_vz_h3[2]->SetLineColor(kGreen);
 //  hi_ecal_vz_h3[2]->Draw("SAME");
 //  hi_ecal_vz_n3[2]->SetLineColor(kGreen+2);
 //  hi_ecal_vz_n3[2]->Draw("SAME");

 //   TLegend *leg33 = new TLegend(0.7,0.75,0.96,0.96);
 //    leg33->SetTextSize(.04);
 //    leg33->AddEntry(hi_ecal_vz_all3[2],"All","l");
 //    leg33->AddEntry(hi_ecal_vz_e3[2],"electrons","l");
 //    leg33->AddEntry(hi_ecal_vz_g3[2],"photons","l");
 //    leg33->AddEntry(hi_ecal_vz_n3[2],"neutron","l");
 //    leg33->AddEntry(hi_ecal_vz_h3[2],"hadron","l");
 //    leg33->Draw();

 //  c_vz3->cd(4);
 //  hi_ecal_vz_all3[3]->SetMinimum(0.001);
 //  hi_ecal_vz_all3[3]->Draw("H");
 //  hi_ecal_vz_e3[3]->SetLineColor(2);
 //  hi_ecal_vz_e3[3]->Draw("SAME");
 //  hi_ecal_vz_g3[3]->SetLineColor(4);
 //  hi_ecal_vz_g3[3]->Draw("SAME");
 //  hi_ecal_vz_h3[3]->SetLineColor(kGreen);
 //  hi_ecal_vz_h3[3]->Draw("SAME");
 //  hi_ecal_vz_n3[3]->SetLineColor(kGreen+2);
 //  hi_ecal_vz_n3[3]->Draw("SAME");

 //   TLegend *leg34 = new TLegend(0.7,0.75,0.96,0.96);
 //    leg34->SetTextSize(.04);
 //    leg34->AddEntry(hi_ecal_vz_all3[3],"All","l");
 //    leg34->AddEntry(hi_ecal_vz_e3[3],"electrons","l");
 //    leg34->AddEntry(hi_ecal_vz_g3[3],"photons","l");
 //    leg34->AddEntry(hi_ecal_vz_n3[3],"neutron","l");
 //    leg34->AddEntry(hi_ecal_vz_h3[3],"hadron","l");
 //    leg34->Draw();

 //  c_vz3->cd(5);
 //  hi_ecal_vz_all3[4]->SetMinimum(0.001);
 //  hi_ecal_vz_all3[4]->Draw("H");
 //  hi_ecal_vz_e3[4]->SetLineColor(2);
 //  hi_ecal_vz_e3[4]->Draw("SAME");
 //  hi_ecal_vz_g3[4]->SetLineColor(4);
 //  hi_ecal_vz_g3[4]->Draw("SAME");
 //  hi_ecal_vz_h3[4]->SetLineColor(kGreen);
 //  hi_ecal_vz_h3[4]->Draw("SAME");
 //  hi_ecal_vz_n3[4]->SetLineColor(kGreen+2);
 //  hi_ecal_vz_n3[4]->Draw("SAME");

 //   TLegend *leg35 = new TLegend(0.7,0.75,0.96,0.96);
 //    leg35->SetTextSize(.04);
 //    leg35->AddEntry(hi_ecal_vz_all3[4],"All","l");
 //    leg35->AddEntry(hi_ecal_vz_e3[4],"electrons","l");
 //    leg35->AddEntry(hi_ecal_vz_g3[4],"photons","l");
 //    leg35->AddEntry(hi_ecal_vz_n3[4],"neutron","l");
 //    leg35->AddEntry(hi_ecal_vz_h3[4],"hadron","l");
 //    leg35->Draw();

 //  c_vz3->cd(6);
 //  hi_ecal_vz_all3[5]->SetMinimum(0.001);
 //  hi_ecal_vz_all3[5]->Draw("H");
 //  hi_ecal_vz_e3[5]->SetLineColor(2);
 //  hi_ecal_vz_e3[5]->Draw("SAME");
 //  hi_ecal_vz_g3[5]->SetLineColor(4);
 //  hi_ecal_vz_g3[5]->Draw("SAME");
 //  hi_ecal_vz_h3[5]->SetLineColor(kGreen);
 //  hi_ecal_vz_h3[5]->Draw("SAME");
 //  hi_ecal_vz_n3[5]->SetLineColor(kGreen+2);
 //  hi_ecal_vz_n3[5]->Draw("SAME");

 //   TLegend *leg36 = new TLegend(0.7,0.75,0.96,0.96);
 //    leg36->SetTextSize(.04);
 //    leg36->AddEntry(hi_ecal_vz_all3[5],"All","l");
 //    leg36->AddEntry(hi_ecal_vz_e3[5],"electrons","l");
 //    leg36->AddEntry(hi_ecal_vz_g3[5],"photons","l");
 //    leg36->AddEntry(hi_ecal_vz_n3[5],"neutron","l");
 //    leg36->AddEntry(hi_ecal_vz_h3[5],"hadron","l");
 //    leg36->Draw();
 //  c_vz3->Print("ecal_occupancy.pdf)");

   gui.Run(1);

}  
