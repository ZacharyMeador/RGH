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
   int npcalhit;
   vector<double> *pcal_sector=new vector<double>;
   vector<double> *pcal_layer=new vector<double>;
   vector<double> *pcal_view=new vector<double>;
   vector<double> *pcal_strip=new vector<double>;
   vector<double> *pcal_ADC=new vector<double>;
   vector<double> *pcal_TDC=new vector<double>;
   vector<double> *pcal_Edep=new vector<double>;
   vector<double> *pcal_E=new vector<double>;
   vector<double> *pcal_x=new vector<double>;
   vector<double> *pcal_y=new vector<double>;
   vector<double> *pcal_z=new vector<double>;
   vector<double> *pcal_lx=new vector<double>;
   vector<double> *pcal_ly=new vector<double>;
   vector<double> *pcal_lz=new vector<double>;
   vector<double> *pcal_t=new vector<double>;
   vector<double> *pcal_pid=new vector<double>;
   vector<double> *pcal_vx=new vector<double>;
   vector<double> *pcal_vy=new vector<double>;
   vector<double> *pcal_vz=new vector<double>;
// ec
//   int nechit;
//   vector<double> *ec_sector=new vector<double>;
//   vector<double> *ec_layer=new vector<double>;
//   vector<double> *ec_view=new vector<double>;
//   vector<double> *ec_strip=new vector<double>;
//   vector<double> *ec_ADC=new vector<double>;
//   vector<double> *ec_TDC=new vector<double>;
//   vector<double> *ec_Edep=new vector<double>;
//   vector<double> *ec_E=new vector<double>;
//   vector<double> *ec_x=new vector<double>;
//   vector<double> *ec_y=new vector<double>;
//   vector<double> *ec_z=new vector<double>;
//   vector<double> *ec_lx=new vector<double>;
//   vector<double> *ec_ly=new vector<double>;
//   vector<double> *ec_lz=new vector<double>;
//   vector<double> *ec_t=new vector<double>;
//   vector<double> *ec_pid=new vector<double>;
//   vector<double> *ec_vx=new vector<double>;
//  vector<double> *ec_vy=new vector<double>;
//   vector<double> *ec_vz=new vector<double>;


   cout << "Creating Tree chains" << endl;
   TChain *pcal= new TChain("pcal");
   pcal->Add("out*.root");
//   TChain *ec= new TChain("ec");
//   ec->Add("out*.root");


// PCAL
   pcal->SetBranchAddress("sector" ,&pcal_sector);
   pcal->SetBranchAddress("module" ,&pcal_layer);
   pcal->SetBranchAddress("view"   ,&pcal_view);
   pcal->SetBranchAddress("strip"  ,&pcal_strip);
   pcal->SetBranchAddress("ADC"    ,&pcal_ADC);
   pcal->SetBranchAddress("TDC"    ,&pcal_TDC);
   pcal->SetBranchAddress("trackE" ,&pcal_E);
   pcal->SetBranchAddress("totEdep",&pcal_Edep);
   pcal->SetBranchAddress("avg_t"  ,&pcal_t);
   pcal->SetBranchAddress("pid"    ,&pcal_pid);
   pcal->SetBranchAddress("avg_x"  ,&pcal_x);
   pcal->SetBranchAddress("avg_y"  ,&pcal_y);
   pcal->SetBranchAddress("avg_z"  ,&pcal_z);
   pcal->SetBranchAddress("avg_lx" ,&pcal_lx);
   pcal->SetBranchAddress("avg_ly" ,&pcal_ly);
   pcal->SetBranchAddress("avg_lz" ,&pcal_lz); 
   pcal->SetBranchAddress("vx"     ,&pcal_vx);
   pcal->SetBranchAddress("vy"     ,&pcal_vy);
   pcal->SetBranchAddress("vz"     ,&pcal_vz);

// EC
//   ec->SetBranchAddress("sector" ,&ec_sector);
//   ec->SetBranchAddress("stack"  ,&ec_layer);
//   ec->SetBranchAddress("view"   ,&ec_view);
//   ec->SetBranchAddress("strip"  ,&ec_strip);
//   ec->SetBranchAddress("ADC"    ,&ec_ADC);
//   ec->SetBranchAddress("TDC"    ,&ec_TDC);
//   ec->SetBranchAddress("trackE" ,&ec_E);
//   ec->SetBranchAddress("totEdep",&ec_Edep);
//   ec->SetBranchAddress("avg_t"  ,&ec_t);
//   ec->SetBranchAddress("pid"    ,&ec_pid);
//   ec->SetBranchAddress("avg_x"  ,&ec_x);
//   ec->SetBranchAddress("avg_y"  ,&ec_y);
//   ec->SetBranchAddress("avg_z"  ,&ec_z);
//   ec->SetBranchAddress("avg_lx" ,&ec_lx);
//   ec->SetBranchAddress("avg_ly" ,&ec_ly);
//   ec->SetBranchAddress("avg_lz" ,&ec_lz); 
//   ec->SetBranchAddress("vx"     ,&ec_vx);
//   ec->SetBranchAddress("vy"     ,&ec_vy);
//   ec->SetBranchAddress("vz"     ,&ec_vz);



   Long64_t nentries = pcal->GetEntries();
   cout << "N. entries:" << nentries << " " << pcal->GetEntries() << endl;



// Create histos
   TH2F *hi_pcal_occ[3];
   TH2F *hi_pcal_occ_cut[3]; 
   TH1F *hi_pcal_vz_all[3];
   TH1F *hi_pcal_vz_e[3];
   TH1F *hi_pcal_vz_g[3];
   TH1F *hi_pcal_vz_h[3];
   TH1F *hi_pcal_vz_n[3];
   TH2F *hi_pcal_origin_all[3]; 
   TH1F *hi_pcal_edep[3]; 

   for(int i=0; i<3; i++) {
       int nstrip=68;
       if(i>0) nstrip=36;
       hi_pcal_occ[i] = new TH2F(Form("hi_pcal_occ %d",i+1), "",nstrip, 1.,nstrip*1.+1.,3,1.,4.);
       hi_pcal_occ[i]->GetXaxis()->SetTitle("Strip");
       hi_pcal_occ[i]->GetYaxis()->SetTitle("View");
       hi_pcal_occ_cut[i] = new TH2F(Form("hi_pcal_occ_cut %d",i+1), "",nstrip, 1.,nstrip*1.+1.,3,1.,4.);
       hi_pcal_occ_cut[i]->GetXaxis()->SetTitle("Strip");
       hi_pcal_occ_cut[i]->GetYaxis()->SetTitle("View");
       hi_pcal_vz_all[i] = new TH1F(Form("hi_pcal_vz_all %d",i+1), "",300, 0., 900.);
       hi_pcal_vz_all[i]->GetXaxis()->SetTitle("Z_{vertex} (cm)");
       hi_pcal_vz_all[i]->GetYaxis()->SetTitle("Rate (MHz/3cm)");
       hi_pcal_vz_e[i] = new TH1F(Form("hi_pcal_vz_e %d",i+1), "",300, 0., 900.);
       hi_pcal_vz_e[i]->GetXaxis()->SetTitle("Z_{vertex} (cm)");
       hi_pcal_vz_e[i]->GetYaxis()->SetTitle("Rate (MHz/3cm)");
       hi_pcal_vz_g[i] = new TH1F(Form("hi_pcal_vz_g %d",i+1), "",300, 0., 900.);
       hi_pcal_vz_g[i]->GetXaxis()->SetTitle("Z_{vertex} (cm)");
       hi_pcal_vz_g[i]->GetYaxis()->SetTitle("Rate (MHz/3cm)");
       hi_pcal_vz_h[i] = new TH1F(Form("hi_pcal_vz_h %d",i+1), "",300, 0., 900.);
       hi_pcal_vz_h[i]->GetXaxis()->SetTitle("Z_{vertex} (cm)");
       hi_pcal_vz_h[i]->GetYaxis()->SetTitle("Rate (MHz/3cm)");
       hi_pcal_vz_n[i] = new TH1F(Form("hi_pcal_vz_n %d",i+1), "",300, 0., 900.);
       hi_pcal_vz_n[i]->GetXaxis()->SetTitle("Z_{vertex} (cm)");
       hi_pcal_vz_n[i]->GetYaxis()->SetTitle("Rate (MHz/3cm)");
       hi_pcal_origin_all[i] = new TH2F(Form("hi_pcal_origin_all %d",i+1), "",300, 0., 900.,200, 0., 600.);
       hi_pcal_origin_all[i]->GetXaxis()->SetTitle("Z_{vertex} (cm)");
       hi_pcal_origin_all[i]->GetYaxis()->SetTitle("R_{vertex} (cm)");
       hi_pcal_edep[i] = new TH1F(Form("hi_pcal_edep %d",i+1), "",200, 0.,10.);
       hi_pcal_edep[i]->GetXaxis()->SetTitle("E(MeV)");
       hi_pcal_edep[i]->GetYaxis()->SetTitle("Rate (MHz/50keV)");
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
  
       pcal_sector->clear();
       pcal_layer->clear();
       pcal_view->clear();
       pcal_strip->clear();
       pcal_ADC->clear();
       pcal_TDC->clear();
       pcal_E->clear();
       pcal_Edep->clear();
       pcal_t->clear();
       pcal_pid->clear();
       pcal_x->clear();
       pcal_y->clear();
       pcal_z->clear();
       pcal_lx->clear();
       pcal_ly->clear();
       pcal_lz->clear(); 
       pcal_vx->clear();
       pcal_vy->clear();
       pcal_vz->clear();
     
//       ec_sector->clear();
//       ec_layer->clear();
//       ec_view->clear();
//       ec_strip->clear();
//       ec_ADC->clear();
//       ec_TDC->clear();
//       ec_E->clear();
//       ec_Edep->clear();
//       ec_t->clear();
//       ec_pid->clear();
//       ec_x->clear();
//       ec_y->clear();
//       ec_z->clear();
//       ec_lx->clear();
//       ec_ly->clear();
//       ec_lz->clear(); 
//       ec_vx->clear();
//       ec_vy->clear();
//       ec_vz->clear();
     
       int nb = pcal->GetEntry(jentry); 
       pcal->GetEntry(jentry); 
       ngoodentries++;
       if(int(jentry/1000)*1000==jentry) cout << "Analyzed " << jentry << " events of " << nentries << endl;
	
       npcalhit=pcal_pid->size();
//       nechit=ec_pid->size();

       for(int i=0; i<npcalhit; i++) {
  	    hi_pcal_occ[0]->Fill((*pcal_strip)[i]*1.,(*pcal_view)[i]*1.);
	      hi_pcal_vz_all[0]->Fill((*pcal_vz)[i]/10.);
        if(abs((*pcal_pid)[i])==11) hi_pcal_vz_e[0]->Fill((*pcal_vz)[i]/10.);
        else if((*pcal_pid)[i]==22) hi_pcal_vz_g[0]->Fill((*pcal_vz)[i]/10.);
   	    else                       hi_pcal_vz_h[0]->Fill((*pcal_vz)[i]/10.);
        if((*pcal_pid)[i]==2112)    hi_pcal_vz_n[0]->Fill((*pcal_vz)[i]/10.);
	      hi_pcal_origin_all[0]->Fill((*pcal_vz)[i]/10.,sqrt((*pcal_vx)[i]*(*pcal_vx)[i]/100.+(*pcal_vy)[i]*(*pcal_vy)[i]/100.));
	      hi_pcal_edep[0]->Fill((*pcal_Edep)[i]);
	      if((*pcal_Edep)[i]>Ethr) {
	        hi_pcal_occ_cut[0]->Fill((*pcal_strip)[i]*1.,(*pcal_view)[i]*1.);
	      }
       }
//       for(int i=0; i<nechit; i++) {
//	   int istack=(*pcal_layer)[i];
//  	   hi_ecal_occ[istack]->Fill((*ec_strip)[i]*1.,(*ec_view)[i]*1.);
//	   hi_ecal_vz_all[istack]->Fill((*ec_vz)[i]/10.);
//	   if(abs((*ec_pid)[i])==11) hi_ecal_vz_e[istack]->Fill((*ec_vz)[i]/10.);
//	   else if((*ec_pid)[i]==22) hi_ecal_vz_g[istack]->Fill((*ec_vz)[i]/10.);
//	   else                       hi_ecal_vz_h[istack]->Fill((*ec_vz)[i]/10.);
//	   if((*ec_pid)[i]==2112)    hi_ecal_vz_n[istack]->Fill((*ec_vz)[i]/10.);
//	   hi_ecal_origin_all[istack]->Fill((*ec_vz)[i]/10.,sqrt((*ec_vx)[i]*(*ec_vx)[i]/100.+(*ec_vy)[i]*(*ec_vy)[i]/100.));
//	   hi_pcal_edep[istack]->Fill((*pcal_Edep)[i]);
//	   if((*ec_Edep)[i]>Ethr) {
//	     hi_ecal_occ_cut[istack]->Fill((*ec_strip)[i]*1.,(*ec_view)[i]*1.);
//	   }
//       }
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
   for(int i=0; i<3; i++) {
       hi_pcal_occ[i]->Scale(1000/time/6);
       hi_pcal_occ_cut[i]->Scale(1000/time/6);
       hi_pcal_vz_all[i]->Scale(1/time);
       hi_pcal_vz_e[i]->Scale(1/time);
       hi_pcal_vz_g[i]->Scale(1/time);
       hi_pcal_vz_h[i]->Scale(1/time);
       hi_pcal_vz_n[i]->Scale(1/time);
       hi_pcal_origin_all[i]->Scale(1/time);
       hi_pcal_edep[i]->Scale(1/time);
   }
   

   TCanvas *c_occ=new TCanvas("c_occ","Occupancy",750,1000);
   c_occ->Divide(1,3);
   FILE *fp = fopen("pcal_occupancy.txt","w");
   for(int i=0; i<3; i++) {
       c_occ->cd(i+1);
       hi_pcal_occ[i]->Draw("COLZ");
       for(int iv=0; iv<hi_pcal_occ[i]->GetNbinsY(); iv++) {
	   for(int is=0; is<hi_pcal_occ[i]->GetNbinsX(); is++) {
	     int layer = i*3+iv+1;
	     int strip = is +1;
	     float value = hi_pcal_occ[i]->GetBinContent(is+1,iv+1);
	     fprintf(fp,"%d\t%d\t%8.1f\n",layer,strip,value);
	   }
       }
   }
   fclose(fp);
   c_occ->Print("pcal_occupancy.pdf(");

   TCanvas *c_occ_cut=new TCanvas("c_occ_cut","Occupancy_Cuts",750,1000);
   c_occ_cut->Divide(1,3);
   fp = fopen("pcal_occupancy_cut.txt","w");
   for(int i=0; i<3; i++) {
       c_occ_cut->cd(i+1);
       hi_pcal_occ_cut[i]->Draw("COLZ");
       for(int iv=0; iv<hi_pcal_occ_cut[i]->GetNbinsY(); iv++) {
	   for(int is=0; is<hi_pcal_occ_cut[i]->GetNbinsX(); is++) {
	     int layer = i*3+iv+1;
	     int strip = is +1;
	     float value = hi_pcal_occ_cut[i]->GetBinContent(is+1,iv+1);
	     fprintf(fp,"%d\t%d\t%8.1f\n",layer,strip,value);
	   }
       }
   }
   fclose(fp);
   c_occ_cut->Print("pcal_occupancy.pdf");


   TCanvas *c_origin=new TCanvas("c_origin","Origin",750,1000);
   c_origin->Divide(1,3);
   for(int i=0; i<3; i++) {
       c_origin->cd(i+1);
       gPad->SetLogz();
       hi_pcal_origin_all[i]->Draw("COLZ");
   }
   c_origin->Print("pcal_occupancy.pdf");

   TCanvas *c_vz=new TCanvas("c_vz","VZ",750,1000);
   c_vz->Divide(1,3);
   for(int i=0; i<3; i++) {
       c_vz->cd(i+1);
       hi_pcal_vz_all[i]->SetMinimum(0.001);
       hi_pcal_vz_all[i]->Draw("H");
       hi_pcal_vz_e[i]->SetLineColor(2);
       hi_pcal_vz_e[i]->Draw("SAME");
       hi_pcal_vz_g[i]->SetLineColor(4);
       hi_pcal_vz_g[i]->Draw("SAME");
       hi_pcal_vz_h[i]->SetLineColor(kGreen);
       hi_pcal_vz_h[i]->Draw("SAME");
       hi_pcal_vz_n[i]->SetLineColor(kGreen+2);
       hi_pcal_vz_n[i]->Draw("SAME");
   }
   c_vz->Print("pcal_occupancy.pdf");

   TCanvas *c_edep=new TCanvas("c_edep","Deposited Energy",750,1000);
   c_edep->Divide(1,3);
   for(int i=0; i<3; i++) {
       c_edep->cd(i+1);
       gPad->SetLogy();
       hi_pcal_edep[i]->Draw();
   }
   c_edep->Print("pcal_occupancy.pdf)");


   gui.Run(1);

}
