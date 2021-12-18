#include <cstdlib>
#include <iostream>
#include <cmath>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TF1.h>
#include <TLine.h>
#include <TGraph.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TColor.h>
#include "TFile.h"
#include <TLegend.h>
#include "TApplication.h"
#include <TROOT.h>

using namespace std;
TApplication gui("GUI",0,NULL);

int main() {

    gStyle->SetCanvasBorderMode(0);
    gStyle->SetCanvasColor(10);

    gStyle->SetPadBorderMode(0);
    gStyle->SetPadLeftMargin(0.12);
    gStyle->SetPadRightMargin(0.12);
    gStyle->SetPadColor(10);


    gStyle->SetTitleFont(72,"X");
    gStyle->SetTitleFont(72,"Y");
    gStyle->SetTitleOffset(0.9,"X");
    gStyle->SetTitleOffset(1.2,"Y");
    gStyle->SetTitleSize(0.05,"X");
    gStyle->SetTitleSize(0.05,"Y");
    
    gStyle->SetLabelFont(72,"X");
    gStyle->SetLabelFont(72,"Y");
    gStyle->SetLabelFont(72,"Z");
//    gStyle->SetOptFit(111);
    gStyle->SetOptStat("nemriou");
    gStyle->SetOptStat("");

    gStyle->SetPalette(1);

    // generated
    int   ngen;
    vector<double>  *pid_gen = new vector<double>;
    // dc
    const int NDCHITMAX = 100000;
    int   ndchit;
    vector<double> *sector = new vector<double>;
    vector<double>  *layer = new vector<double>;
    vector<double>   *wire = new vector<double>;
    vector<double>   *hitn = new vector<double>;
    vector<double>   *pid  = new vector<double>;
    vector<double>   *mpid = new vector<double>;
    vector<double>   *Edep = new vector<double>;
    vector<double>      *E = new vector<double>;
    vector<double>     *vx = new vector<double>;
    vector<double>     *vy = new vector<double>;
    vector<double>     *vz = new vector<double>;
    vector<double>    *mvx = new vector<double>;
    vector<double>    *mvy = new vector<double>;
    vector<double>    *mvz = new vector<double>;
    vector<double> *avg_x = new vector<double>;
    vector<double> *avg_y = new vector<double>;
	

    
    // GENERATED
    TChain *g= new TChain("generated");
    g->Add("out*.root");
    g->SetBranchAddress("pid",&pid_gen);
    
    
    // DC
    TChain *d= new TChain("dc");
    d->Add("out*.root");
    d->SetBranchAddress("wire",&wire);
    d->SetBranchAddress("layer",&layer);
    d->SetBranchAddress("sector",&sector);
    d->SetBranchAddress("hitn",&hitn);
    d->SetBranchAddress("trackE",&E);
    d->SetBranchAddress("totEdep",&Edep);
    d->SetBranchAddress("pid",&pid);
    d->SetBranchAddress("mpid",&mpid);
    d->SetBranchAddress("mvx",&mvx);
    d->SetBranchAddress("mvy",&mvy);
    d->SetBranchAddress("mvz",&mvz);
    d->SetBranchAddress("vx",&vx);
    d->SetBranchAddress("vy",&vy);
    d->SetBranchAddress("vz",&vz);
    d->SetBranchAddress("avg_x",&avg_x);
    d->SetBranchAddress("avg_y",&avg_y);
   
    
    Long64_t nentries = d->GetEntries();
    cout << nentries << endl;

    char *histname = new char[50];
    TH2F *hi_dcocc[6];

    // Lines 109-351 create histos
    // creating histos for wire vs layer
    for(int i=0; i<6; i++){
        hi_dcocc[i]= new TH2F(Form("DC Occ. E > 50 eV Sector%i",i+1),Form("DC Occ. E > 50 eV Sector%i",i+1),112, 1., 113., 36, 1., 37.);
        hi_dcocc[i]->GetXaxis()->SetTitle("wire");
        hi_dcocc[i]->GetYaxis()->SetTitle("layer");
    }
    TH2F *hi_dcocc_tgt   = new TH2F("DC Occ. Target", "DC Occ. Target",112, 1.,113., 36, 1.,37.);
    hi_dcocc_tgt->GetXaxis()->SetTitle("wire");
    hi_dcocc_tgt->GetYaxis()->SetTitle("layer");

    //creating histo for average occupancy sector vs %
    TH1F *hi_dcocc_region[3];
    for(int i=0; i<3; i++) {
      //  sprintf(histname,"hi_dcocc_region%i",i);
	hi_dcocc_region[i]= new TH1F(Form("hi_dcocc_region%i",i+1),"",6,0.5,6.5);
	hi_dcocc_region[i]->GetXaxis()->SetTitle("Sector");
	hi_dcocc_region[i]->GetYaxis()->SetTitle("Occupancy (%)");
    }

    TH2F *hi_bg_origin   = new TH2F("Origin of Bg", "",160, -500.,3500., 160, 0.,800.);
    hi_bg_origin->SetTitle("Origin of Bg");
    hi_bg_origin->GetXaxis()->SetTitle("z(mm)");                                        
    hi_bg_origin->GetYaxis()->SetTitle("r(mm)");                                        

    // ceating histo for hit detection for z vs r: lines 133 - 201
    TH2F *hi_bg_r_vs_z_reg1[6];
    TH2F *hi_bg_r_vs_z_vs_ene_reg1[6];
    TH2F *hi_bg_r_vs_z_vs_ene_reg_temp1[6];
    
    TH2F *hi_bg_r_vs_z_reg2[6];
    TH2F *hi_bg_r_vs_z_vs_ene_reg2[6];
    TH2F *hi_bg_r_vs_z_vs_ene_reg_temp2[6];
    
    TH2F *hi_bg_r_vs_z_reg3[6];
    TH2F *hi_bg_r_vs_z_vs_ene_reg3[6];
    TH2F *hi_bg_r_vs_z_vs_ene_reg_temp3[6];
    
    TH2F *hi_bg_y_vs_x_reg[3];
    
    for(int i=0; i<6; i++) {
       // sprintf(histname,"hi_bg_r_vs_z_region%i",i);
        hi_bg_r_vs_z_reg1[i]= new TH2F(Form("hi_bg_r_vs_z_reg1_sec%i",i+1),Form("hi_bg_r_vs_z_reg1_sec%i",i+1),200, -200., 6500., 200, 0.,2000.);
        hi_bg_r_vs_z_reg1[i]->GetXaxis()->SetTitle("z(mm)");
        hi_bg_r_vs_z_reg1[i]->GetYaxis()->SetTitle("r(mm)");
      //  sprintf(histname,"hi_bg_y_vs_x_region%i",i);
        
        
        hi_bg_r_vs_z_vs_ene_reg1[i]= new TH2F(Form("hi_bg_r_vs_z_vs_ene1.%i",i+1),Form("hi_bg_r_vs_z_vs_ene1.%i",i+1),200, -200., 6500., 200, 0.,2000.);
        hi_bg_r_vs_z_vs_ene_reg1[i]->GetXaxis()->SetTitle("z(mm)");
        hi_bg_r_vs_z_vs_ene_reg1[i]->GetYaxis()->SetTitle("r(mm)");
        hi_bg_r_vs_z_vs_ene_reg1[i]->GetZaxis()->SetTitle("Energy(MeV)");
        

        hi_bg_r_vs_z_vs_ene_reg_temp1[i]= new TH2F(Form("hi_bg_r_vs_z_vs_ene_temp1.%i",i+1),Form("hi_bg_r_vs_z_vs_ene_temp1.%i",i+1),200, -200., 6500., 200, 0.,2000.);
        hi_bg_r_vs_z_vs_ene_reg_temp1[i]->GetXaxis()->SetTitle("z(mm)");
        hi_bg_r_vs_z_vs_ene_reg_temp1[i]->GetYaxis()->SetTitle("r(mm)");
        hi_bg_r_vs_z_vs_ene_reg_temp1[i]->GetZaxis()->SetTitle("Energy(MeV)");
	
	// sprintf(histname,"hi_bg_r_vs_z_region%i",i);
        hi_bg_r_vs_z_reg2[i]= new TH2F(Form("hi_bg_r_vs_z_reg2_sec%i",i+1),Form("hi_bg_r_vs_z_reg2_sec%i",i+1),200, -200., 6500., 200, 0.,2000.);
        hi_bg_r_vs_z_reg2[i]->GetXaxis()->SetTitle("z(mm)");
        hi_bg_r_vs_z_reg2[i]->GetYaxis()->SetTitle("r(mm)");
      //  sprintf(histname,"hi_bg_y_vs_x_region%i",i);
        
        
        hi_bg_r_vs_z_vs_ene_reg2[i]= new TH2F(Form("hi_bg_r_vs_z_vs_ene2.%i",i+1),Form("hi_bg_r_vs_z_vs_ene2.%i",i+1),200, -200., 6500., 200, 0.,2000.);
        hi_bg_r_vs_z_vs_ene_reg2[i]->GetXaxis()->SetTitle("z(mm)");
        hi_bg_r_vs_z_vs_ene_reg2[i]->GetYaxis()->SetTitle("r(mm)");
        hi_bg_r_vs_z_vs_ene_reg2[i]->GetZaxis()->SetTitle("Energy(MeV)");
        

        hi_bg_r_vs_z_vs_ene_reg_temp2[i]= new TH2F(Form("hi_bg_r_vs_z_vs_ene_temp2.%i",i+1),Form("hi_bg_r_vs_z_vs_ene_temp2.%i",i+1),200, -200., 6500., 200, 0.,2000.);
        hi_bg_r_vs_z_vs_ene_reg_temp2[i]->GetXaxis()->SetTitle("z(mm)");
        hi_bg_r_vs_z_vs_ene_reg_temp2[i]->GetYaxis()->SetTitle("r(mm)");
        hi_bg_r_vs_z_vs_ene_reg_temp2[i]->GetZaxis()->SetTitle("Energy(MeV)");
        
	// sprintf(histname,"hi_bg_r_vs_z_region%i",i);
        hi_bg_r_vs_z_reg3[i]= new TH2F(Form("hi_bg_r_vs_z_reg3_sec%i",i+1),Form("hi_bg_r_vs_z_reg2_sec%i",i+1),200, -200., 6500., 200, 0.,2000.);
        hi_bg_r_vs_z_reg3[i]->GetXaxis()->SetTitle("z(mm)");
        hi_bg_r_vs_z_reg3[i]->GetYaxis()->SetTitle("r(mm)");
      //  sprintf(histname,"hi_bg_y_vs_x_region%i",i);
        
        
        hi_bg_r_vs_z_vs_ene_reg3[i]= new TH2F(Form("hi_bg_r_vs_z_vs_ene3.%i",i+1),Form("hi_bg_r_vs_z_vs_ene3.%i",i+1),200, -200., 6500., 200, 0.,2000.);
        hi_bg_r_vs_z_vs_ene_reg3[i]->GetXaxis()->SetTitle("z(mm)");
        hi_bg_r_vs_z_vs_ene_reg3[i]->GetYaxis()->SetTitle("r(mm)");
        hi_bg_r_vs_z_vs_ene_reg3[i]->GetZaxis()->SetTitle("Energy(MeV)");
        

        hi_bg_r_vs_z_vs_ene_reg_temp3[i]= new TH2F(Form("hi_bg_r_vs_z_vs_ene_temp3.%i",i+1),Form("hi_bg_r_vs_z_vs_ene_temp3.%i",i+1),200, -200., 6500., 200, 0.,2000.);
        hi_bg_r_vs_z_vs_ene_reg_temp3[i]->GetXaxis()->SetTitle("z(mm)");
        hi_bg_r_vs_z_vs_ene_reg_temp3[i]->GetYaxis()->SetTitle("r(mm)");
        hi_bg_r_vs_z_vs_ene_reg_temp3[i]->GetZaxis()->SetTitle("Energy(MeV)");
    }
    // ceating histos for hit distribution x(mm) vs y(mm)
    for(int i=0; i<3; i++){
        hi_bg_y_vs_x_reg[i]= new TH2F(Form("hi_bg_y_vs_x_region%i",i+1), Form("hi_bg_y_vs_x_region%i",i+1),200, -1000.,1000., 200, -1000.,1000.);
            hi_bg_y_vs_x_reg[i]->GetXaxis()->SetTitle("x(mm)");
            hi_bg_y_vs_x_reg[i]->GetYaxis()->SetTitle("y(mm)");
    }
    TH2F *hi_bg_energy_tmp   = new TH2F("Energy of Bg tmp", "",160, -500.,3500., 160, 0.,800.);
    hi_bg_energy_tmp->GetXaxis()->SetTitle("z(mm)");
    hi_bg_energy_tmp->GetYaxis()->SetTitle("r(mm)");

    TH2F *hi_bg_energy   = new TH2F("Energy of Bg", "Energy of Bg",160, -500.,3500., 160, 0.,800.);
    hi_bg_energy->GetXaxis()->SetTitle("z(mm)");
    hi_bg_energy->GetYaxis()->SetTitle("r(mm)");
    hi_bg_energy->GetZaxis()->SetTitle("Energy(MeV)");
    
    TH1F *hi_bg_z   = new TH1F("Vz of Bg", "Vz of Bg",400, -500.,3500.);
    //hi_bg_z->SetTitle("Vz of Bg");
    hi_bg_z->GetXaxis()->SetTitle("z(mm)");
    hi_bg_z->GetYaxis()->SetTitle("Rate (MHz)");
    
    // creating histos for origin of particle z(mm) vs Rates(MHz) 
    TH1F *hi_bg_z_reg1[6];
    TH1F *hi_bg_z_e_reg1[6];
    TH1F *hi_bg_z_g_reg1[6];
    TH1F *hi_bg_z_p_reg1[6];
    TH1F *hi_bg_z_pi_reg1[6];
    TH1F *hi_bg_z_n_reg1[6];
    TH1F *hi_bg_z_o_reg1[6];
    TH1F *hi_bg_z_reg2[6];
    TH1F *hi_bg_z_e_reg2[6];
    TH1F *hi_bg_z_g_reg2[6];
    TH1F *hi_bg_z_p_reg2[6];
    TH1F *hi_bg_z_pi_reg2[6];
    TH1F *hi_bg_z_n_reg2[6];
    TH1F *hi_bg_z_o_reg2[6];
    TH1F *hi_bg_z_reg3[6];
    TH1F *hi_bg_z_e_reg3[6];
    TH1F *hi_bg_z_g_reg3[6];
    TH1F *hi_bg_z_p_reg3[6];
    TH1F *hi_bg_z_pi_reg3[6];
    TH1F *hi_bg_z_n_reg3[6];
    TH1F *hi_bg_z_o_reg3[6];
    
    for(int i=0; i<6; i++) {
    //    sprintf(histname,"hi_bg_z_region%i",i);
    hi_bg_z_reg1[i]= new TH1F(Form("hi_bg_z_reg1_sec%i",i+1), Form("hi_bg_z_reg1_sec%i",i+1),200, -200.,6500.);
	hi_bg_z_reg1[i]->GetXaxis()->SetTitle("z(mm)");
	hi_bg_z_reg1[i]->GetYaxis()->SetTitle("Rate (MHz)");

        
    hi_bg_z_e_reg1[i]= new TH1F(Form("hi_bg_z_e1.%i",i+1), Form("hi_bg_z_e1.%i",i+1),200, -200.,6500.);
    hi_bg_z_e_reg1[i]->GetXaxis()->SetTitle("z(mm)");
    hi_bg_z_e_reg1[i]->GetYaxis()->SetTitle("Rate (MHz)");
        
    hi_bg_z_g_reg1[i]= new TH1F(Form("hi_bg_z_g%i",i+1), Form("hi_bg_z_g%i",i+1),200, -200.,6500.);
    hi_bg_z_g_reg1[i]->GetXaxis()->SetTitle("z(mm)");
    hi_bg_z_g_reg1[i]->GetYaxis()->SetTitle("Rate (MHz)");
    
    hi_bg_z_n_reg1[i]= new TH1F(Form("hi_bg_z_n1.%i",i+1), Form("hi_bg_z_n1.%i",i+1),200, -200.,6500.);
    hi_bg_z_n_reg1[i]->GetXaxis()->SetTitle("z(mm)");
    hi_bg_z_n_reg1[i]->GetYaxis()->SetTitle("Rate (MHz)");
        
    hi_bg_z_p_reg1[i]= new TH1F(Form("hi_bg_z_p1.%i",i+1), Form("hi_bg_z_p1.%i",i+1),200, -200.,6500.);
    hi_bg_z_p_reg1[i]->GetXaxis()->SetTitle("z(mm)");
    hi_bg_z_p_reg1[i]->GetYaxis()->SetTitle("Rate (MHz)");

    hi_bg_z_pi_reg1[i]= new TH1F(Form("hi_bg_z_pi1.%i",i+1), Form("hi_bg_z_pi1.%i",i+1),200, -200.,6500.);
    hi_bg_z_pi_reg1[i]->GetXaxis()->SetTitle("z(mm)");
    hi_bg_z_pi_reg1[i]->GetYaxis()->SetTitle("Rate (MHz)");
        
    hi_bg_z_o_reg1[i]= new TH1F(Form("hi_bg_z_o1.%i",i+1), Form("hi_bg_z_o1.%i",i+1),200, -200.,6500.);
    hi_bg_z_o_reg1[i]->GetXaxis()->SetTitle("z(mm)");
    hi_bg_z_o_reg1[i]->GetYaxis()->SetTitle("Rate (MHz)");
    
    hi_bg_z_reg2[i]= new TH1F(Form("hi_bg_z_reg2_sec%i",i+1), Form("hi_bg_z_reg2_sec%i",i+1),200, -200.,6500.);
	hi_bg_z_reg2[i]->GetXaxis()->SetTitle("z(mm)");
	hi_bg_z_reg2[i]->GetYaxis()->SetTitle("Rate (MHz)");

        
    hi_bg_z_e_reg2[i]= new TH1F(Form("hi_bg_z_e2.%i",i+1), Form("hi_bg_z_e2.%i",i+1),200, -200.,6500.);
    hi_bg_z_e_reg2[i]->GetXaxis()->SetTitle("z(mm)");
    hi_bg_z_e_reg2[i]->GetYaxis()->SetTitle("Rate (MHz)");
        
    hi_bg_z_g_reg2[i]= new TH1F(Form("hi_bg_z_g2.%i",i+1), Form("hi_bg_z_g2.%i",i+1),200, -200.,6500.);
    hi_bg_z_g_reg2[i]->GetXaxis()->SetTitle("z(mm)");
    hi_bg_z_g_reg2[i]->GetYaxis()->SetTitle("Rate (MHz)");
    
    hi_bg_z_n_reg2[i]= new TH1F(Form("hi_bg_z_n2.%i",i+1), Form("hi_bg_z_n2.%i",i+1),200, -200.,6500.);
    hi_bg_z_n_reg2[i]->GetXaxis()->SetTitle("z(mm)");
    hi_bg_z_n_reg2[i]->GetYaxis()->SetTitle("Rate (MHz)");
        
    hi_bg_z_p_reg2[i]= new TH1F(Form("hi_bg_z_p2.%i",i+1), Form("hi_bg_z_p2.%i",i+1),200, -200.,6500.);
    hi_bg_z_p_reg2[i]->GetXaxis()->SetTitle("z(mm)");
    hi_bg_z_p_reg2[i]->GetYaxis()->SetTitle("Rate (MHz)");

    hi_bg_z_pi_reg2[i]= new TH1F(Form("hi_bg_z_pi2.%i",i+1), Form("hi_bg_z_pi2.%i",i+1),200, -200.,6500.);
    hi_bg_z_pi_reg2[i]->GetXaxis()->SetTitle("z(mm)");
    hi_bg_z_pi_reg2[i]->GetYaxis()->SetTitle("Rate (MHz)");
        
    hi_bg_z_o_reg2[i]= new TH1F(Form("hi_bg_z_o2.%i",i+1), Form("hi_bg_z_o2.%i",i+1),200, -200.,6500.);
    hi_bg_z_o_reg2[i]->GetXaxis()->SetTitle("z(mm)");
    hi_bg_z_o_reg2[i]->GetYaxis()->SetTitle("Rate (MHz)");
        
    hi_bg_z_reg3[i]= new TH1F(Form("hi_bg_z_reg3_sec%i",i+1), Form("hi_bg_z_reg3_sec%i",i+1),200, -200.,6500.);
	hi_bg_z_reg3[i]->GetXaxis()->SetTitle("z(mm)");
	hi_bg_z_reg3[i]->GetYaxis()->SetTitle("Rate (MHz)");

        
    hi_bg_z_e_reg3[i]= new TH1F(Form("hi_bg_z_e3.%i",i+1), Form("hi_bg_z_e3.%i",i+1),200, -200.,6500.);
    hi_bg_z_e_reg3[i]->GetXaxis()->SetTitle("z(mm)");
    hi_bg_z_e_reg3[i]->GetYaxis()->SetTitle("Rate (MHz)");
        
    hi_bg_z_g_reg3[i]= new TH1F(Form("hi_bg_z_g3.%i",i+1), Form("hi_bg_zg3.%i",i+1),200, -200.,6500.);
    hi_bg_z_g_reg3[i]->GetXaxis()->SetTitle("z(mm)");
    hi_bg_z_g_reg3[i]->GetYaxis()->SetTitle("Rate (MHz)");
    
    hi_bg_z_n_reg3[i]= new TH1F(Form("hi_bg_z_n3.%i",i+1), Form("hi_bg_z_n3.%i",i+1),200, -200.,6500.);
    hi_bg_z_n_reg3[i]->GetXaxis()->SetTitle("z(mm)");
    hi_bg_z_n_reg3[i]->GetYaxis()->SetTitle("Rate (MHz)");
        
    hi_bg_z_p_reg3[i]= new TH1F(Form("hi_bg_z_p3.%i",i+1), Form("hi_bg_z_p3.%i",i+1),200, -200.,6500.);
    hi_bg_z_p_reg3[i]->GetXaxis()->SetTitle("z(mm)");
    hi_bg_z_p_reg3[i]->GetYaxis()->SetTitle("Rate (MHz)");

    hi_bg_z_pi_reg3[i]= new TH1F(Form("hi_bg_z_pi3.%i",i+1), Form("hi_bg_z_pi3.%i",i+1),200, -200.,6500.);
    hi_bg_z_pi_reg3[i]->GetXaxis()->SetTitle("z(mm)");
    hi_bg_z_pi_reg3[i]->GetYaxis()->SetTitle("Rate (MHz)");
        
    hi_bg_z_o_reg3[i]= new TH1F(Form("hi_bg_z_o3.%i",i+1), Form("hi_bg_z_o3.%i",i+1),200, -200.,6500.);
    hi_bg_z_o_reg3[i]->GetXaxis()->SetTitle("z(mm)");
    hi_bg_z_o_reg3[i]->GetYaxis()->SetTitle("Rate (MHz)");
    }


    TH1F *hi_bg_z_e   = new TH1F("Vz of Bg electrons", "",400, -500.,3500.);
    hi_bg_z_e->GetXaxis()->SetTitle("z(mm)");
    hi_bg_z_e->GetYaxis()->SetTitle("Rate (MHz)");

    TH1F *hi_bg_z_g   = new TH1F("Vz of Bg Photons", "",400, -500.,3500.);
    hi_bg_z_g->GetXaxis()->SetTitle("z(mm)");
    hi_bg_z_g->GetYaxis()->SetTitle("Rate (MHz)");
    
    TH1F *hi_bg_z_o   = new TH1F("Vz of Bg Other", "",400, -500.,3500.);
    hi_bg_z_o->GetXaxis()->SetTitle("z(mm)");
    hi_bg_z_o->GetYaxis()->SetTitle("Rate (MHz)");

    TH1F *hi_bg_E  = new TH1F("momentum of Bg", "momentum of Bg",200,0.,10.);
   // hi_bg_E->SetTitle("momentum of Bg");
    hi_bg_E->GetXaxis()->SetTitle("p(MeV)");
    hi_bg_E->GetYaxis()->SetTitle("Rate (MHz)");


// in standard run mode (1 electron at a time)
// float lumi=nentries*5.*0.07*0.602/1000; // (mbarn^-1 or 10^27 cm^-2)
// in luminosity mode (124000 electrons in 250ns window)
    float levents=124000;
    int   nsum=124000./levents;
    int   nfull=0;
 
    cout << "number of events to integrate = " << nsum  << endl;

    int nint=0;
    int ngoodentries=0;
    double mass=0;
    double dc_weight;
    const double pi = 3.14159265358979323846;
    double phi;
    
    for(Long64_t jentry=0; jentry < nentries; jentry++) {
//      for(Long64_t jentry=0; jentry < 1000; jentry++) {

        
        // clear vectors
        pid_gen->clear();
        sector->clear();
        layer->clear();
        wire->clear();
        pid->clear();
        mpid->clear();
        E->clear();
        Edep->clear();
        vx->clear();
        vy->clear();
        vz->clear();
        mvx->clear();
        mvy->clear();
        mvz->clear();
        avg_x->clear();
        avg_y->clear();
        
        d->GetEntry(jentry);
        g->GetEntry(jentry);
        ngen=pid_gen->size();
        ndchit=sector->size();
        if(ngen>0) ngoodentries++;
        if(int(jentry/1000)*1000==jentry) cout << "Analyzed " << jentry << " events of " << nentries << endl;
            
	for(int i=0; i<ndchit; i++) {
	    int it=(*hitn)[i]-1;
         // extracting phi using the hit index where (*vx) and (*vy) are the x- and y- vertix components (respectively)
         phi = atan2((*vy)[it],(*vx)[it])*180/pi;	// atan2() returns [0, pi] if (*vy) >= 0 and (-pi, 0) if (*vy) < 0; 180/pi converts radian to degree	 
         phi += (phi<0) ? 360:0;	// using ternary operator to find actual angel if (*vy) < 0
	 // sorting hit into sectors
         int nsect = floor((phi-30)/60) + 2; // subtracting 30 from phi to account for sector 1 range being [0, 30) and [330, 0]. Each sector is 60deg. Adding 2 to get sector# 
         if(nsect == 7) nsect = 1; // sector 1 will return either 1 or 7
	 // defining mass of particles pid is particle id in gemc
	    if((*pid)[it]==2112 || (*pid)[it]==2212) {
                mass=938;
            }
            else{
                mass=0;
            }
	  // defining layers: region 1 = layers 1-12, region 2 = 13-24, and region 25-36
            int dc_reg=int(((*layer)[i]-1)/12)+1;
            if(dc_reg==1) dc_weight=1;
            else          dc_weight=2;
            
	    if((*Edep)[it]>0.00005) {
		// Filling arrays
                hi_dcocc[nsect-1]->Fill((*wire)[i],(*layer)[i],dc_weight);
	        if((*vz)[it]<100)hi_dcocc_tgt->Fill((*wire)[i],(*layer)[i],dc_weight);
		hi_dcocc_region[dc_reg-1]->Fill((*sector)[i],dc_weight);
	        if(dc_reg==1) {
                    hi_bg_origin->Fill((*vz)[it],sqrt((*vx)[it]*(*vx)[it]+(*vy)[it]*(*vy)[it]));
                    hi_bg_z->Fill((*vz)[it]);
                    if((*pid)[it]==11) {
                        hi_bg_z_e->Fill((*vz)[it]);
                    }
                    else if((*pid)[it]==22) {
                        hi_bg_z_g->Fill((*vz)[it]);
                    }
                    else {
                        hi_bg_z_o->Fill((*vz)[it]);
                    }
                    if((*E)[it]>mass) hi_bg_energy_tmp->Fill((*vz)[it],sqrt((*vx)[it]*(*vx)[it]+(*vy)[it]*(*vy)[it]),sqrt((*E)[it]*(*E)[it]-mass*mass));
                  
                    hi_bg_E->Fill(sqrt((*E)[it]*(*E)[it]-mass*mass));
                }
	    // Lines 441-481 Filling arrays for histos by region
            if(dc_reg==1){
                    hi_bg_r_vs_z_vs_ene_reg_temp1[nsect-1]->Fill((*vz)[it],sqrt((*vx)[it]*(*vx)[it]+(*vy)[it]*(*vy)[it]),sqrt((*E)[it]*(*E)[it]-mass*mass));
                    hi_bg_r_vs_z_reg1[nsect-1]->Fill((*vz)[it],sqrt((*vx)[it]*(*vx)[it]+(*vy)[it]*(*vy)[it]));

                    if((*pid)[it]==11) hi_bg_z_e_reg1[nsect-1]->Fill((*vz)[it]);
                    if((*pid)[it]==22) hi_bg_z_g_reg1[nsect-1]->Fill((*vz)[it]);
                    if((*pid)[it]==2112) hi_bg_z_n_reg1[nsect-1]->Fill((*vz)[it]);
                    if((*pid)[it]==2212) hi_bg_z_p_reg1[nsect-1]->Fill((*vz)[it]);
                    if((*pid)[it]==211 || (*pid)[it]==-211 ) hi_bg_z_pi_reg1[nsect-1]->Fill((*vz)[it]);
                    if((*pid)[it]!=211 && (*pid)[it]!=-211 && (*pid)[it]!=11 && (*pid)[it]!=22 && (*pid)[it]!=2212) hi_bg_z_o_reg1[nsect-1]->Fill((*vz)[it]);
                    
			    hi_bg_z_reg1[nsect-1]->Fill((*vz)[it]);
                }
            if(dc_reg==2){
                    hi_bg_r_vs_z_vs_ene_reg_temp2[nsect-1]->Fill((*vz)[it],sqrt((*vx)[it]*(*vx)[it]+(*vy)[it]*(*vy)[it]),sqrt((*E)[it]*(*E)[it]-mass*mass));
                    hi_bg_r_vs_z_reg2[nsect-1]->Fill((*vz)[it],sqrt((*vx)[it]*(*vx)[it]+(*vy)[it]*(*vy)[it]));

                    if((*pid)[it]==11) hi_bg_z_e_reg2[nsect-1]->Fill((*vz)[it]);
                    if((*pid)[it]==22) hi_bg_z_g_reg2[nsect-1]->Fill((*vz)[it]);
                    if((*pid)[it]==2112) hi_bg_z_n_reg2[nsect-1]->Fill((*vz)[it]);
                    if((*pid)[it]==2212) hi_bg_z_p_reg2[nsect-1]->Fill((*vz)[it]);
                    if((*pid)[it]==211 || (*pid)[it]==-211 ) hi_bg_z_pi_reg2[nsect-1]->Fill((*vz)[it]);
                    if((*pid)[it]!=211 && (*pid)[it]!=-211 && (*pid)[it]!=11 && (*pid)[it]!=22 && (*pid)[it]!=2212) hi_bg_z_o_reg2[nsect-1]->Fill((*vz)[it]);

                    hi_bg_z_reg2[nsect-1]->Fill((*vz)[it]);
                }

            if(dc_reg==3){
                    hi_bg_r_vs_z_vs_ene_reg_temp3[nsect-1]->Fill((*vz)[it],sqrt((*vx)[it]*(*vx)[it]+(*vy)[it]*(*vy)[it]),sqrt((*E)[it]*(*E)[it]-mass*mass));
                    hi_bg_r_vs_z_reg3[nsect-1]->Fill((*vz)[it],sqrt((*vx)[it]*(*vx)[it]+(*vy)[it]*(*vy)[it]));

                    if((*pid)[it]==11) hi_bg_z_e_reg3[nsect-1]->Fill((*vz)[it]);
                    if((*pid)[it]==22) hi_bg_z_g_reg3[nsect-1]->Fill((*vz)[it]);
                    if((*pid)[it]==2112) hi_bg_z_n_reg3[nsect-1]->Fill((*vz)[it]);
                    if((*pid)[it]==2212) hi_bg_z_p_reg3[nsect-1]->Fill((*vz)[it]);
                    if((*pid)[it]==211 || (*pid)[it]==-211 ) hi_bg_z_pi_reg3[nsect-1]->Fill((*vz)[it]);
                    if((*pid)[it]!=211 && (*pid)[it]!=-211 && (*pid)[it]!=11 && (*pid)[it]!=22 && (*pid)[it]!=2212) hi_bg_z_o_reg3[nsect-1]->Fill((*vz)[it]);

                    hi_bg_z_reg3[nsect-1]->Fill((*vz)[it]);
                }       
            if((*vz)[it]>1000*(dc_reg+1)) hi_bg_y_vs_x_reg[dc_reg-1]->Fill((*vx)[it],(*vy)[it]);
            
           }
        }
    }

    // calculating normalization factors based on number of good events
    float norm=124000/levents/ngoodentries;
    float lumi=ngoodentries*250/10*levents/124000; // (mbarn^-1 or 10^27 cm^-2)
    float time=250/norm;
    cout << norm/1. << " "<< ngoodentries << endl;

    cout << "normalization factor          = " << norm  << endl;
    cout << "Run time                      = " << time  << " ns" << endl;
    cout << nfull << " " << 1/norm << endl;

    // normalizing rate histogram
    // If histo is integrated over all 6 sectors divide by 6
    // If histo is integrate by sector divide by 1
    for(int i=0; i<6; i++){
        hi_dcocc[i]->Scale(norm/1.);
    }
        hi_dcocc_tgt->Scale(norm/6.);
    
    hi_bg_z->Scale(1000./time);
    hi_bg_z_e->Scale(1000./time);
    hi_bg_z_g->Scale(1000./time);
    hi_bg_z_o->Scale(1000./time);
    hi_bg_E->Scale(1000./time);
    for(int i=0; i<6; i++) {
      hi_bg_z_reg1[i]->Scale(1000./time);
      hi_bg_z_e_reg1[i]->Scale(1000./time);
      hi_bg_z_g_reg1[i]->Scale(1000./time);
      hi_bg_z_p_reg1[i]->Scale(1000./time);
      hi_bg_z_pi_reg1[i]->Scale(1000./time);
      hi_bg_z_n_reg1[i]->Scale(1000./time);
      hi_bg_z_o_reg1[i]->Scale(1000./time);
    
      hi_bg_z_reg2[i]->Scale(1000./time);
      hi_bg_z_e_reg2[i]->Scale(1000./time);
      hi_bg_z_g_reg2[i]->Scale(1000./time);
      hi_bg_z_p_reg2[i]->Scale(1000./time);
      hi_bg_z_pi_reg2[i]->Scale(1000./time);
      hi_bg_z_n_reg2[i]->Scale(1000./time);
      hi_bg_z_o_reg2[i]->Scale(1000./time);
        
      hi_bg_z_reg3[i]->Scale(1000./time);
      hi_bg_z_e_reg3[i]->Scale(1000./time);
      hi_bg_z_g_reg3[i]->Scale(1000./time);
      hi_bg_z_p_reg3[i]->Scale(1000./time);
      hi_bg_z_pi_reg3[i]->Scale(1000./time);
      hi_bg_z_n_reg3[i]->Scale(1000./time);
      hi_bg_z_o_reg3[i]->Scale(1000./time);
      
      
      // Lines 537-539 break code with DifferentDimension exception raised
//       hi_bg_r_vs_z_vs_ene_reg1[i]->Divide(hi_bg_r_vs_z_vs_ene_reg_temp1[i],hi_bg_r_vs_z_reg1[i]);
//       hi_bg_r_vs_z_vs_ene_reg2[i]->Divide(hi_bg_r_vs_z_vs_ene_reg_temp2[i],hi_bg_r_vs_z_reg2[i]);
//       hi_bg_r_vs_z_vs_ene_reg3[i]->Divide(hi_bg_r_vs_z_vs_ene_reg_temp3[i],hi_bg_r_vs_z_reg3[i]);
    }
    for(int i=0; i<3; i++){
      hi_dcocc_region[i]->Sumw2();
      hi_dcocc_region[i]->Scale(100*norm/112/12);
    }
    
    hi_bg_energy->Divide(hi_bg_energy_tmp,hi_bg_origin);
// following lines print histos to pdf
// Lines 548-572 print wire vs layer histos
    TCanvas *c1=new TCanvas("c1","Occupancy",750,1000);
    c1->Divide(3,3);
    c1->cd(1);
//     gPad->SetLogz();
    hi_dcocc[0]->Draw("COLZ");
    c1->cd(2);
//     gPad->SetLogz();
    hi_dcocc[1]->Draw("COLZ");
    c1->cd(3);
//     gPad->SetLogz();
    hi_dcocc[2]->Draw("COLZ");
    c1->cd(4);
//     gPad->SetLogz();
    hi_dcocc[3]->Draw("COLZ");
    c1->cd(5);
//     gPad->SetLogz();
    hi_dcocc[4]->Draw("COLZ");
    c1->cd(6);
//     gPad->SetLogz();
    hi_dcocc[5]->Draw("COLZ");
    c1->cd(7);
//     gPad->SetLogz();
    hi_dcocc_tgt->Draw("COLZ");
    c1->Print("dc_occ.pdf(");
// lines 574-594 print region 1 z(mm) vs r(mm)
    TCanvas *m1=new TCanvas("m1","Background Origin Region 1.1",750,1000);
    m1->Divide(2,3);
    m1->cd(1);
    gPad->SetLogz();
    hi_bg_r_vs_z_reg1[0]->Draw("COLZ");
    m1->cd(2);
    gPad->SetLogz();
    hi_bg_r_vs_z_reg1[1]->Draw("COLZ");
    m1->cd(3);
    gPad->SetLogz();
    hi_bg_r_vs_z_reg1[2]->Draw("COLZ");
    m1->cd(4);
    gPad->SetLogz();
    hi_bg_r_vs_z_reg1[3]->Draw("COLZ");
    m1->cd(5);
    gPad->SetLogz();
    hi_bg_r_vs_z_reg1[4]->Draw("COLZ");
    m1->cd(6);
    gPad->SetLogz();
    hi_bg_r_vs_z_reg1[5]->Draw("COLZ");
    m1->Print("dc_occ.pdf");
// //     TCanvas *m2=new TCanvas("m3","Background Origin Region 1.2",750,1000);
// //     m2->Divide(2,3);
// //     for(int i=0; i<6; i++){
// // 	    m2->cd(i);
// // 	    hi_bg_r_vs_z_vs_ene_reg1[i]->SetMinimum(0);
// // 	    hi_bg_r_vs_z_vs_ene_reg1[i]->SetMaximum(100);
// // 	    hi_bg_r_vs_z_vs_ene_reg1[i]->Draw("COLZ");
// // 	    }
//     m2->Print("dc_occ.pdf");
// lines 605-805 print z(mm) vs Rate(MHz) region 1 
    TCanvas *m3=new TCanvas("m3","Background Origin Region 1.3",750,1000);
    m3->Divide(2,3);
    m3->cd(1);
    gPad->SetLogy();
    hi_bg_z_reg1[0]->Draw("H");

    hi_bg_z_e_reg1[0]->SetLineColor(2);
    hi_bg_z_e_reg1[0]->Draw("SAMEH");

    hi_bg_z_g_reg1[0]->SetLineColor(4);
    hi_bg_z_g_reg1[0]->Draw("SAMEH");

    hi_bg_z_p_reg1[0]->SetLineColor(91);
    hi_bg_z_p_reg1[0]->Draw("SAMEH");

    hi_bg_z_n_reg1[0]->SetLineColor(8);
    hi_bg_z_n_reg1[0]->Draw("SAMEH");

    hi_bg_z_pi_reg1[0]->SetLineColor(3);
    hi_bg_z_pi_reg1[0]->Draw("SAMEH");

    hi_bg_z_o_reg1[0]->SetLineColor(6);
    hi_bg_z_o_reg1[0]->Draw("SAMEH");


    TLegend *leg1= new TLegend(0.7,0.75,0.96,0.96);
    leg1->SetTextSize(.04);
    leg1->AddEntry(hi_bg_z_reg1[0],"All","l");
    leg1->AddEntry(hi_bg_z_e_reg1[0],"electrons","l");
    leg1->AddEntry(hi_bg_z_g_reg1[0],"photons","l");
    leg1->AddEntry(hi_bg_z_p_reg1[0],"proton","l");
    leg1->AddEntry(hi_bg_z_pi_reg1[0],"pion","l");
    leg1->AddEntry(hi_bg_z_n_reg1[0],"neutron","l");
    leg1->AddEntry(hi_bg_z_o_reg1[0],"other","l");
    leg1->Draw();
    m3->cd(2);
    gPad->SetLogy();
    hi_bg_z_reg1[1]->Draw("H");

    hi_bg_z_e_reg1[1]->SetLineColor(2);
    hi_bg_z_e_reg1[1]->Draw("SAMEH");

    hi_bg_z_g_reg1[1]->SetLineColor(4);
    hi_bg_z_g_reg1[1]->Draw("SAMEH");

    hi_bg_z_p_reg1[1]->SetLineColor(91);
    hi_bg_z_p_reg1[1]->Draw("SAMEH");

    hi_bg_z_n_reg1[1]->SetLineColor(8);
    hi_bg_z_n_reg1[1]->Draw("SAMEH");

    hi_bg_z_pi_reg1[1]->SetLineColor(3);
    hi_bg_z_pi_reg1[1]->Draw("SAMEH");

    hi_bg_z_o_reg1[1]->SetLineColor(6);
    hi_bg_z_o_reg1[1]->Draw("SAMEH");


    TLegend *leg2= new TLegend(0.7,0.75,0.96,0.96);
    leg2->SetTextSize(.04);
    leg2->AddEntry(hi_bg_z_reg1[1],"All","l");
    leg2->AddEntry(hi_bg_z_e_reg1[1],"electrons","l");
    leg2->AddEntry(hi_bg_z_g_reg1[1],"photons","l");
    leg2->AddEntry(hi_bg_z_p_reg1[1],"proton","l");
    leg2->AddEntry(hi_bg_z_pi_reg1[1],"pion","l");
    leg2->AddEntry(hi_bg_z_n_reg1[1],"neutron","l");
    leg2->AddEntry(hi_bg_z_o_reg1[1],"other","l");
    leg2->Draw();
    m3->cd(3);
    gPad->SetLogy();
    hi_bg_z_reg1[2]->Draw("H");

    hi_bg_z_e_reg1[2]->SetLineColor(2);
    hi_bg_z_e_reg1[2]->Draw("SAMEH");

    hi_bg_z_g_reg1[2]->SetLineColor(4);
    hi_bg_z_g_reg1[2]->Draw("SAMEH");

    hi_bg_z_p_reg1[2]->SetLineColor(91);
    hi_bg_z_p_reg1[2]->Draw("SAMEH");

    hi_bg_z_n_reg1[2]->SetLineColor(8);
    hi_bg_z_n_reg1[2]->Draw("SAMEH");

    hi_bg_z_pi_reg1[2]->SetLineColor(3);
    hi_bg_z_pi_reg1[2]->Draw("SAMEH");

    hi_bg_z_o_reg1[2]->SetLineColor(6);
    hi_bg_z_o_reg1[2]->Draw("SAMEH");


    TLegend *leg3= new TLegend(0.7,0.75,0.96,0.96);
    leg3->SetTextSize(.04);
    leg3->AddEntry(hi_bg_z_reg1[2],"All","l");
    leg3->AddEntry(hi_bg_z_e_reg1[2],"electrons","l");
    leg3->AddEntry(hi_bg_z_g_reg1[2],"photons","l");
    leg3->AddEntry(hi_bg_z_p_reg1[2],"proton","l");
    leg3->AddEntry(hi_bg_z_pi_reg1[2],"pion","l");
    leg3->AddEntry(hi_bg_z_n_reg1[2],"neutron","l");
    leg3->AddEntry(hi_bg_z_o_reg1[2],"other","l");
    leg3->Draw();
    m3->cd(4);
    gPad->SetLogy();
    hi_bg_z_reg1[3]->Draw("H");

    hi_bg_z_e_reg1[3]->SetLineColor(2);
    hi_bg_z_e_reg1[3]->Draw("SAMEH");

    hi_bg_z_g_reg1[3]->SetLineColor(4);
    hi_bg_z_g_reg1[3]->Draw("SAMEH");

    hi_bg_z_p_reg1[3]->SetLineColor(91);
    hi_bg_z_p_reg1[3]->Draw("SAMEH");

    hi_bg_z_n_reg1[3]->SetLineColor(8);
    hi_bg_z_n_reg1[3]->Draw("SAMEH");

    hi_bg_z_pi_reg1[3]->SetLineColor(3);
    hi_bg_z_pi_reg1[3]->Draw("SAMEH");

    hi_bg_z_o_reg1[3]->SetLineColor(6);
    hi_bg_z_o_reg1[3]->Draw("SAMEH");


    TLegend *leg4= new TLegend(0.7,0.75,0.96,0.96);
    leg4->SetTextSize(.04);
    leg4->AddEntry(hi_bg_z_reg1[3],"All","l");
    leg4->AddEntry(hi_bg_z_e_reg1[3],"electrons","l");
    leg4->AddEntry(hi_bg_z_g_reg1[3],"photons","l");
    leg4->AddEntry(hi_bg_z_p_reg1[3],"proton","l");
    leg4->AddEntry(hi_bg_z_pi_reg1[3],"pion","l");
    leg4->AddEntry(hi_bg_z_n_reg1[3],"neutron","l");
    leg4->AddEntry(hi_bg_z_o_reg1[3],"other","l");
    leg4->Draw();
    m3->cd(5);
    gPad->SetLogy();
    hi_bg_z_reg1[0]->Draw("H");

    hi_bg_z_e_reg1[4]->SetLineColor(2);
    hi_bg_z_e_reg1[4]->Draw("SAMEH");

    hi_bg_z_g_reg1[4]->SetLineColor(4);
    hi_bg_z_g_reg1[4]->Draw("SAMEH");

    hi_bg_z_p_reg1[4]->SetLineColor(91);
    hi_bg_z_p_reg1[4]->Draw("SAMEH");

    hi_bg_z_n_reg1[4]->SetLineColor(8);
    hi_bg_z_n_reg1[4]->Draw("SAMEH");

    hi_bg_z_pi_reg1[4]->SetLineColor(3);
    hi_bg_z_pi_reg1[4]->Draw("SAMEH");

    hi_bg_z_o_reg1[4]->SetLineColor(6);
    hi_bg_z_o_reg1[4]->Draw("SAMEH");


    TLegend *leg5= new TLegend(0.7,0.75,0.96,0.96);
    leg5->SetTextSize(.04);
    leg5->AddEntry(hi_bg_z_reg1[4],"All","l");
    leg5->AddEntry(hi_bg_z_e_reg1[4],"electrons","l");
    leg5->AddEntry(hi_bg_z_g_reg1[4],"photons","l");
    leg5->AddEntry(hi_bg_z_p_reg1[4],"proton","l");
    leg5->AddEntry(hi_bg_z_pi_reg1[4],"pion","l");
    leg5->AddEntry(hi_bg_z_n_reg1[4],"neutron","l");
    leg5->AddEntry(hi_bg_z_o_reg1[4],"other","l");
    leg5->Draw();
    m3->cd(6);
    gPad->SetLogy();
    hi_bg_z_reg1[5]->Draw("H");

    hi_bg_z_e_reg1[5]->SetLineColor(2);
    hi_bg_z_e_reg1[5]->Draw("SAMEH");

    hi_bg_z_g_reg1[5]->SetLineColor(4);
    hi_bg_z_g_reg1[5]->Draw("SAMEH");

    hi_bg_z_p_reg1[5]->SetLineColor(91);
    hi_bg_z_p_reg1[5]->Draw("SAMEH");

    hi_bg_z_n_reg1[5]->SetLineColor(8);
    hi_bg_z_n_reg1[5]->Draw("SAMEH");

    hi_bg_z_pi_reg1[5]->SetLineColor(3);
    hi_bg_z_pi_reg1[5]->Draw("SAMEH");

    hi_bg_z_o_reg1[5]->SetLineColor(6);
    hi_bg_z_o_reg1[5]->Draw("SAMEH");


    TLegend *leg6= new TLegend(0.7,0.75,0.96,0.96);
    leg6->SetTextSize(.04);
    leg6->AddEntry(hi_bg_z_reg1[5],"All","l");
    leg6->AddEntry(hi_bg_z_e_reg1[5],"electrons","l");
    leg6->AddEntry(hi_bg_z_g_reg1[5],"photons","l");
    leg6->AddEntry(hi_bg_z_p_reg1[5],"proton","l");
    leg6->AddEntry(hi_bg_z_pi_reg1[5],"pion","l");
    leg6->AddEntry(hi_bg_z_n_reg1[5],"neutron","l");
    leg6->AddEntry(hi_bg_z_o_reg1[5],"other","l");
    leg6->Draw();
    m3->Print("dc_occ.pdf");
    
// Lines 808-828 print z(mm) vs r(mm) region 2	
    TCanvas *m4=new TCanvas("m4","Background Origin Region 2.1",750,1000);
    m4->Divide(2,3);
    m4->cd(1);
    gPad->SetLogz();
    hi_bg_r_vs_z_reg2[0]->Draw("COLZ");
    m4->cd(2);
    gPad->SetLogz();
    hi_bg_r_vs_z_reg2[1]->Draw("COLZ");
    m4->cd(3);
    gPad->SetLogz();
    hi_bg_r_vs_z_reg2[2]->Draw("COLZ");
    m4->cd(4);
    gPad->SetLogz();
    hi_bg_r_vs_z_reg2[3]->Draw("COLZ");
    m4->cd(5);
    gPad->SetLogz();
    hi_bg_r_vs_z_reg2[4]->Draw("COLZ");
    m4->cd(6);
    gPad->SetLogz();
    hi_bg_r_vs_z_reg2[5]->Draw("COLZ");
    m4->Print("dc_occ.pdf");
// //     TCanvas *m5=new TCanvas("m5","Background Origin Region 2.2",750,1000);
// //     m5->Divide(2,3);
// //     for(int i=6; i<12; i++){
// // 	    m5->cd(i);
// // 	    hi_bg_r_vs_z_vs_ene_reg2[i]->SetMinimum(0);
// // 	    hi_bg_r_vs_z_vs_ene_reg2[i]->SetMaximum(100);
// // 	    hi_bg_r_vs_z_vs_ene_reg2[i]->Draw("COLZ");
// // 	    }
//     m5->Print("dc_occ.pdf");
// Lines 839-1039 print z(mm) vs Rate(MHz) region 2
    TCanvas *m6=new TCanvas("m6","Background Origin Region 2.3",750,1000);
    m6->Divide(2,3);
    m6->cd(1);
    gPad->SetLogy();
    hi_bg_z_reg2[0]->Draw("H");

    hi_bg_z_e_reg2[0]->SetLineColor(2);
    hi_bg_z_e_reg2[0]->Draw("SAMEH");

    hi_bg_z_g_reg2[0]->SetLineColor(4);
    hi_bg_z_g_reg2[0]->Draw("SAMEH");

    hi_bg_z_p_reg2[0]->SetLineColor(91);
    hi_bg_z_p_reg2[0]->Draw("SAMEH");

    hi_bg_z_n_reg2[0]->SetLineColor(8);
    hi_bg_z_n_reg2[0]->Draw("SAMEH");

    hi_bg_z_pi_reg2[0]->SetLineColor(3);
    hi_bg_z_pi_reg2[0]->Draw("SAMEH");

    hi_bg_z_o_reg2[0]->SetLineColor(6);
    hi_bg_z_o_reg2[0]->Draw("SAMEH");


    TLegend *leg12= new TLegend(0.7,0.75,0.96,0.96);
    leg12->SetTextSize(.04);
    leg12->AddEntry(hi_bg_z_reg2[0],"All","l");
    leg12->AddEntry(hi_bg_z_e_reg2[0],"electrons","l");
    leg12->AddEntry(hi_bg_z_g_reg2[0],"photons","l");
    leg12->AddEntry(hi_bg_z_p_reg2[0],"proton","l");
    leg12->AddEntry(hi_bg_z_pi_reg2[0],"pion","l");
    leg12->AddEntry(hi_bg_z_n_reg2[0],"neutron","l");
    leg12->AddEntry(hi_bg_z_o_reg2[0],"other","l");
    leg12->Draw();
    m6->cd(2);
    gPad->SetLogy();
    hi_bg_z_reg2[1]->Draw("H");

    hi_bg_z_e_reg2[1]->SetLineColor(2);
    hi_bg_z_e_reg2[1]->Draw("SAMEH");

    hi_bg_z_g_reg2[1]->SetLineColor(4);
    hi_bg_z_g_reg2[1]->Draw("SAMEH");

    hi_bg_z_p_reg2[1]->SetLineColor(91);
    hi_bg_z_p_reg2[1]->Draw("SAMEH");

    hi_bg_z_n_reg2[1]->SetLineColor(8);
    hi_bg_z_n_reg2[1]->Draw("SAMEH");

    hi_bg_z_pi_reg2[1]->SetLineColor(3);
    hi_bg_z_pi_reg2[1]->Draw("SAMEH");

    hi_bg_z_o_reg2[1]->SetLineColor(6);
    hi_bg_z_o_reg2[1]->Draw("SAMEH");


    TLegend *leg22= new TLegend(0.7,0.75,0.96,0.96);
    leg22->SetTextSize(.04);
    leg22->AddEntry(hi_bg_z_reg1[1],"All","l");
    leg22->AddEntry(hi_bg_z_e_reg1[1],"electrons","l");
    leg22->AddEntry(hi_bg_z_g_reg1[1],"photons","l");
    leg22->AddEntry(hi_bg_z_p_reg1[1],"proton","l");
    leg22->AddEntry(hi_bg_z_pi_reg1[1],"pion","l");
    leg22->AddEntry(hi_bg_z_n_reg1[1],"neutron","l");
    leg22->AddEntry(hi_bg_z_o_reg1[1],"other","l");
    leg22->Draw();
    m6->cd(3);
    gPad->SetLogy();
    hi_bg_z_reg2[2]->Draw("H");

    hi_bg_z_e_reg2[2]->SetLineColor(2);
    hi_bg_z_e_reg2[2]->Draw("SAMEH");

    hi_bg_z_g_reg2[2]->SetLineColor(4);
    hi_bg_z_g_reg2[2]->Draw("SAMEH");

    hi_bg_z_p_reg2[2]->SetLineColor(91);
    hi_bg_z_p_reg2[2]->Draw("SAMEH");

    hi_bg_z_n_reg2[2]->SetLineColor(8);
    hi_bg_z_n_reg2[2]->Draw("SAMEH");

    hi_bg_z_pi_reg2[2]->SetLineColor(3);
    hi_bg_z_pi_reg2[2]->Draw("SAMEH");

    hi_bg_z_o_reg2[2]->SetLineColor(6);
    hi_bg_z_o_reg2[2]->Draw("SAMEH");


    TLegend *leg32= new TLegend(0.7,0.75,0.96,0.96);
    leg32->SetTextSize(.04);
    leg32->AddEntry(hi_bg_z_reg2[2],"All","l");
    leg32->AddEntry(hi_bg_z_e_reg2[2],"electrons","l");
    leg32->AddEntry(hi_bg_z_g_reg2[2],"photons","l");
    leg32->AddEntry(hi_bg_z_p_reg2[2],"proton","l");
    leg32->AddEntry(hi_bg_z_pi_reg2[2],"pion","l");
    leg32->AddEntry(hi_bg_z_n_reg2[2],"neutron","l");
    leg32->AddEntry(hi_bg_z_o_reg2[2],"other","l");
    leg32->Draw();
    m6->cd(4);
    gPad->SetLogy();
    hi_bg_z_reg2[3]->Draw("H");

    hi_bg_z_e_reg2[3]->SetLineColor(2);
    hi_bg_z_e_reg2[3]->Draw("SAMEH");

    hi_bg_z_g_reg2[3]->SetLineColor(4);
    hi_bg_z_g_reg2[3]->Draw("SAMEH");

    hi_bg_z_p_reg2[3]->SetLineColor(91);
    hi_bg_z_p_reg2[3]->Draw("SAMEH");

    hi_bg_z_n_reg2[3]->SetLineColor(8);
    hi_bg_z_n_reg2[3]->Draw("SAMEH");

    hi_bg_z_pi_reg2[3]->SetLineColor(3);
    hi_bg_z_pi_reg2[3]->Draw("SAMEH");

    hi_bg_z_o_reg2[3]->SetLineColor(6);
    hi_bg_z_o_reg2[3]->Draw("SAMEH");


    TLegend *leg42= new TLegend(0.7,0.75,0.96,0.96);
    leg42->SetTextSize(.04);
    leg42->AddEntry(hi_bg_z_reg2[3],"All","l");
    leg42->AddEntry(hi_bg_z_e_reg2[3],"electrons","l");
    leg42->AddEntry(hi_bg_z_g_reg2[3],"photons","l");
    leg42->AddEntry(hi_bg_z_p_reg2[3],"proton","l");
    leg42->AddEntry(hi_bg_z_pi_reg2[3],"pion","l");
    leg42->AddEntry(hi_bg_z_n_reg2[3],"neutron","l");
    leg42->AddEntry(hi_bg_z_o_reg2[3],"other","l");
    leg42->Draw();
    m6->cd(5);
    gPad->SetLogy();
    hi_bg_z_reg2[0]->Draw("H");

    hi_bg_z_e_reg2[4]->SetLineColor(2);
    hi_bg_z_e_reg2[4]->Draw("SAMEH");

    hi_bg_z_g_reg2[4]->SetLineColor(4);
    hi_bg_z_g_reg2[4]->Draw("SAMEH");

    hi_bg_z_p_reg2[4]->SetLineColor(91);
    hi_bg_z_p_reg2[4]->Draw("SAMEH");

    hi_bg_z_n_reg2[4]->SetLineColor(8);
    hi_bg_z_n_reg2[4]->Draw("SAMEH");

    hi_bg_z_pi_reg2[4]->SetLineColor(3);
    hi_bg_z_pi_reg2[4]->Draw("SAMEH");

    hi_bg_z_o_reg2[4]->SetLineColor(6);
    hi_bg_z_o_reg2[4]->Draw("SAMEH");


    TLegend *leg52= new TLegend(0.7,0.75,0.96,0.96);
    leg52->SetTextSize(.04);
    leg52->AddEntry(hi_bg_z_reg2[4],"All","l");
    leg52->AddEntry(hi_bg_z_e_reg2[4],"electrons","l");
    leg52->AddEntry(hi_bg_z_g_reg2[4],"photons","l");
    leg52->AddEntry(hi_bg_z_p_reg2[4],"proton","l");
    leg52->AddEntry(hi_bg_z_pi_reg2[4],"pion","l");
    leg52->AddEntry(hi_bg_z_n_reg2[4],"neutron","l");
    leg52->AddEntry(hi_bg_z_o_reg2[4],"other","l");
    leg52->Draw();
    m6->cd(6);
    gPad->SetLogy();
    hi_bg_z_reg2[5]->Draw("H");

    hi_bg_z_e_reg2[5]->SetLineColor(2);
    hi_bg_z_e_reg2[5]->Draw("SAMEH");

    hi_bg_z_g_reg2[5]->SetLineColor(4);
    hi_bg_z_g_reg2[5]->Draw("SAMEH");

    hi_bg_z_p_reg2[5]->SetLineColor(91);
    hi_bg_z_p_reg2[5]->Draw("SAMEH");

    hi_bg_z_n_reg2[5]->SetLineColor(8);
    hi_bg_z_n_reg2[5]->Draw("SAMEH");

    hi_bg_z_pi_reg2[5]->SetLineColor(3);
    hi_bg_z_pi_reg2[5]->Draw("SAMEH");

    hi_bg_z_o_reg2[5]->SetLineColor(6);
    hi_bg_z_o_reg2[5]->Draw("SAMEH");


    TLegend *leg62= new TLegend(0.7,0.75,0.96,0.96);
    leg62->SetTextSize(.04);
    leg62->AddEntry(hi_bg_z_reg2[5],"All","l");
    leg62->AddEntry(hi_bg_z_e_reg2[5],"electrons","l");
    leg62->AddEntry(hi_bg_z_g_reg2[5],"photons","l");
    leg62->AddEntry(hi_bg_z_p_reg2[5],"proton","l");
    leg62->AddEntry(hi_bg_z_pi_reg2[5],"pion","l");
    leg62->AddEntry(hi_bg_z_n_reg2[5],"neutron","l");
    leg62->AddEntry(hi_bg_z_o_reg2[5],"other","l");
    leg62->Draw();
    m6->Print("dc_occ.pdf");

// lines 1042-1062 print z(mm) vs r(mm)	
    TCanvas *m7=new TCanvas("m7","Background Origin Region 3.1",750,1000);
    m7->Divide(2,3);
    m7->cd(1);
    gPad->SetLogz();
    hi_bg_r_vs_z_reg3[0]->Draw("COLZ");
    m7->cd(2);
    gPad->SetLogz();
    hi_bg_r_vs_z_reg3[1]->Draw("COLZ");
    m7->cd(3);
    gPad->SetLogz();
    hi_bg_r_vs_z_reg3[2]->Draw("COLZ");
    m7->cd(4);
    gPad->SetLogz();
    hi_bg_r_vs_z_reg3[3]->Draw("COLZ");
    m7->cd(5);
    gPad->SetLogz();
    hi_bg_r_vs_z_reg3[4]->Draw("COLZ");
    m7->cd(6);
    gPad->SetLogz();
    hi_bg_r_vs_z_reg3[5]->Draw("COLZ");
    m7->Print("dc_occ.pdf");
// //     TCanvas *m8=new TCanvas("m8","Background Origin Region 3.2",750,1000);
// //     m8->Divide(2,3);
// //     for(int i=12; i<18; i++){
// // 	    m8->cd(i);
// // 	    hi_bg_r_vs_z_vs_ene_reg3[i]->SetMinimum(0);
// // 	    hi_bg_r_vs_z_vs_ene_reg3[i]->SetMaximum(100);
// // 	    hi_bg_r_vs_z_vs_ene_reg3[i]->Draw("COLZ");
// // 	    }
// //     m8->Print("dc_occ.pdf");
//lines 1073-1273 print z(mm) vs Rate(MHz) region 3
    TCanvas *m9=new TCanvas("m9","Background Origin Region 3.3",750,1000);
    m9->Divide(2,3);
    m9->cd(1);
    gPad->SetLogy();
    hi_bg_z_reg3[0]->Draw("H");

    hi_bg_z_e_reg3[0]->SetLineColor(2);
    hi_bg_z_e_reg3[0]->Draw("SAMEH");

    hi_bg_z_g_reg3[0]->SetLineColor(4);
    hi_bg_z_g_reg3[0]->Draw("SAMEH");

    hi_bg_z_p_reg3[0]->SetLineColor(91);
    hi_bg_z_p_reg3[0]->Draw("SAMEH");

    hi_bg_z_n_reg3[0]->SetLineColor(8);
    hi_bg_z_n_reg3[0]->Draw("SAMEH");

    hi_bg_z_pi_reg3[0]->SetLineColor(3);
    hi_bg_z_pi_reg3[0]->Draw("SAMEH");

    hi_bg_z_o_reg3[0]->SetLineColor(6);
    hi_bg_z_o_reg3[0]->Draw("SAMEH");


    TLegend *leg13= new TLegend(0.7,0.75,0.96,0.96);
    leg13->SetTextSize(.04);
    leg13->AddEntry(hi_bg_z_reg3[0],"All","l");
    leg13->AddEntry(hi_bg_z_e_reg3[0],"electrons","l");
    leg13->AddEntry(hi_bg_z_g_reg3[0],"photons","l");
    leg13->AddEntry(hi_bg_z_p_reg3[0],"proton","l");
    leg13->AddEntry(hi_bg_z_pi_reg3[0],"pion","l");
    leg13->AddEntry(hi_bg_z_n_reg3[0],"neutron","l");
    leg13->AddEntry(hi_bg_z_o_reg3[0],"other","l");
    leg13->Draw();
    m9->cd(2);
    gPad->SetLogy();
    hi_bg_z_reg3[1]->Draw("H");

    hi_bg_z_e_reg3[1]->SetLineColor(2);
    hi_bg_z_e_reg3[1]->Draw("SAMEH");

    hi_bg_z_g_reg3[1]->SetLineColor(4);
    hi_bg_z_g_reg3[1]->Draw("SAMEH");

    hi_bg_z_p_reg3[1]->SetLineColor(91);
    hi_bg_z_p_reg3[1]->Draw("SAMEH");

    hi_bg_z_n_reg3[1]->SetLineColor(8);
    hi_bg_z_n_reg3[1]->Draw("SAMEH");

    hi_bg_z_pi_reg3[1]->SetLineColor(3);
    hi_bg_z_pi_reg3[1]->Draw("SAMEH");

    hi_bg_z_o_reg3[1]->SetLineColor(6);
    hi_bg_z_o_reg3[1]->Draw("SAMEH");


    TLegend *leg23= new TLegend(0.7,0.75,0.96,0.96);
    leg23->SetTextSize(.04);
    leg23->AddEntry(hi_bg_z_reg3[1],"All","l");
    leg23->AddEntry(hi_bg_z_e_reg3[1],"electrons","l");
    leg23->AddEntry(hi_bg_z_g_reg3[1],"photons","l");
    leg23->AddEntry(hi_bg_z_p_reg3[1],"proton","l");
    leg23->AddEntry(hi_bg_z_pi_reg3[1],"pion","l");
    leg23->AddEntry(hi_bg_z_n_reg3[1],"neutron","l");
    leg23->AddEntry(hi_bg_z_o_reg3[1],"other","l");
    leg23->Draw();
    m9->cd(3);
    gPad->SetLogy();
    hi_bg_z_reg3[2]->Draw("H");

    hi_bg_z_e_reg3[2]->SetLineColor(2);
    hi_bg_z_e_reg3[2]->Draw("SAMEH");

    hi_bg_z_g_reg3[2]->SetLineColor(4);
    hi_bg_z_g_reg3[2]->Draw("SAMEH");

    hi_bg_z_p_reg3[2]->SetLineColor(91);
    hi_bg_z_p_reg3[2]->Draw("SAMEH");

    hi_bg_z_n_reg3[2]->SetLineColor(8);
    hi_bg_z_n_reg3[2]->Draw("SAMEH");

    hi_bg_z_pi_reg3[2]->SetLineColor(3);
    hi_bg_z_pi_reg3[2]->Draw("SAMEH");

    hi_bg_z_o_reg3[2]->SetLineColor(6);
    hi_bg_z_o_reg3[2]->Draw("SAMEH");


    TLegend *leg33= new TLegend(0.7,0.75,0.96,0.96);
    leg33->SetTextSize(.04);
    leg33->AddEntry(hi_bg_z_reg3[2],"All","l");
    leg33->AddEntry(hi_bg_z_e_reg3[2],"electrons","l");
    leg33->AddEntry(hi_bg_z_g_reg3[2],"photons","l");
    leg33->AddEntry(hi_bg_z_p_reg3[2],"proton","l");
    leg33->AddEntry(hi_bg_z_pi_reg3[2],"pion","l");
    leg33->AddEntry(hi_bg_z_n_reg3[2],"neutron","l");
    leg33->AddEntry(hi_bg_z_o_reg3[2],"other","l");
    leg33->Draw();
    m9->cd(4);
    gPad->SetLogy();
    hi_bg_z_reg3[3]->Draw("H");

    hi_bg_z_e_reg3[3]->SetLineColor(2);
    hi_bg_z_e_reg3[3]->Draw("SAMEH");

    hi_bg_z_g_reg3[3]->SetLineColor(4);
    hi_bg_z_g_reg3[3]->Draw("SAMEH");

    hi_bg_z_p_reg3[3]->SetLineColor(91);
    hi_bg_z_p_reg3[3]->Draw("SAMEH");

    hi_bg_z_n_reg3[3]->SetLineColor(8);
    hi_bg_z_n_reg3[3]->Draw("SAMEH");

    hi_bg_z_pi_reg3[3]->SetLineColor(3);
    hi_bg_z_pi_reg3[3]->Draw("SAMEH");

    hi_bg_z_o_reg3[3]->SetLineColor(6);
    hi_bg_z_o_reg3[3]->Draw("SAMEH");


    TLegend *leg43= new TLegend(0.7,0.75,0.96,0.96);
    leg43->SetTextSize(.04);
    leg43->AddEntry(hi_bg_z_reg3[3],"All","l");
    leg43->AddEntry(hi_bg_z_e_reg3[3],"electrons","l");
    leg43->AddEntry(hi_bg_z_g_reg3[3],"photons","l");
    leg43->AddEntry(hi_bg_z_p_reg3[3],"proton","l");
    leg43->AddEntry(hi_bg_z_pi_reg3[3],"pion","l");
    leg43->AddEntry(hi_bg_z_n_reg3[3],"neutron","l");
    leg43->AddEntry(hi_bg_z_o_reg3[3],"other","l");
    leg43->Draw();
    m9->cd(5);
    gPad->SetLogy();
    hi_bg_z_reg3[0]->Draw("H");

    hi_bg_z_e_reg3[4]->SetLineColor(2);
    hi_bg_z_e_reg3[4]->Draw("SAMEH");

    hi_bg_z_g_reg3[4]->SetLineColor(4);
    hi_bg_z_g_reg3[4]->Draw("SAMEH");

    hi_bg_z_p_reg3[4]->SetLineColor(91);
    hi_bg_z_p_reg3[4]->Draw("SAMEH");

    hi_bg_z_n_reg3[4]->SetLineColor(8);
    hi_bg_z_n_reg3[4]->Draw("SAMEH");

    hi_bg_z_pi_reg3[4]->SetLineColor(3);
    hi_bg_z_pi_reg3[4]->Draw("SAMEH");

    hi_bg_z_o_reg3[4]->SetLineColor(6);
    hi_bg_z_o_reg3[4]->Draw("SAMEH");


    TLegend *leg53= new TLegend(0.7,0.75,0.96,0.96);
    leg53->SetTextSize(.04);
    leg53->AddEntry(hi_bg_z_reg3[4],"All","l");
    leg53->AddEntry(hi_bg_z_e_reg3[4],"electrons","l");
    leg53->AddEntry(hi_bg_z_g_reg3[4],"photons","l");
    leg53->AddEntry(hi_bg_z_p_reg3[4],"proton","l");
    leg53->AddEntry(hi_bg_z_pi_reg3[4],"pion","l");
    leg53->AddEntry(hi_bg_z_n_reg3[4],"neutron","l");
    leg53->AddEntry(hi_bg_z_o_reg3[4],"other","l");
    leg53->Draw();
    m9->cd(6);
    gPad->SetLogy();
    hi_bg_z_reg3[5]->Draw("H");

    hi_bg_z_e_reg3[5]->SetLineColor(2);
    hi_bg_z_e_reg3[5]->Draw("SAMEH");

    hi_bg_z_g_reg3[5]->SetLineColor(4);
    hi_bg_z_g_reg3[5]->Draw("SAMEH");

    hi_bg_z_p_reg3[5]->SetLineColor(91);
    hi_bg_z_p_reg3[5]->Draw("SAMEH");

    hi_bg_z_n_reg3[5]->SetLineColor(8);
    hi_bg_z_n_reg3[5]->Draw("SAMEH");

    hi_bg_z_pi_reg3[5]->SetLineColor(3);
    hi_bg_z_pi_reg3[5]->Draw("SAMEH");

    hi_bg_z_o_reg3[5]->SetLineColor(6);
    hi_bg_z_o_reg3[5]->Draw("SAMEH");


    TLegend *leg63= new TLegend(0.7,0.75,0.96,0.96);
    leg63->SetTextSize(.04);
    leg63->AddEntry(hi_bg_z_reg3[5],"All","l");
    leg63->AddEntry(hi_bg_z_e_reg3[5],"electrons","l");
    leg63->AddEntry(hi_bg_z_g_reg3[5],"photons","l");
    leg63->AddEntry(hi_bg_z_p_reg3[5],"proton","l");
    leg63->AddEntry(hi_bg_z_pi_reg3[5],"pion","l");
    leg63->AddEntry(hi_bg_z_n_reg3[5],"neutron","l");
    leg63->AddEntry(hi_bg_z_o_reg3[5],"other","l");
    leg63->Draw();
    m9->Print("dc_occ.pdf");
    
    TCanvas *c3=new TCanvas("c3","Background Origin",750,1000);
    c3->Divide(1,2);
    c3->cd(1);
    gPad->SetLogz();
    hi_bg_origin->Draw("COLZ");
    c3->cd(2);
    gPad->SetLogy();
    hi_bg_z->SetLineColor(1);
    //    hi_bg_z->SetMaximum(15);
    hi_bg_z->Draw("H");
    hi_bg_z_e->SetLineColor(2);
    hi_bg_z_e->Draw("HSAME");
    hi_bg_z_g->SetLineColor(4);
    hi_bg_z_g->Draw("HSAME");
    hi_bg_z_o->SetLineColor(3);
    hi_bg_z_o->Draw("HSAME");
    
    
    TLegend *leg = new TLegend(0.7,0.75,0.96,0.96);

    // leg->SetFillStyle(0);
   //  leg ->SetBorderSize(0);
     leg->SetTextSize(.04);
     leg->AddEntry(hi_bg_z,"Vz of Bg","l");
     leg->AddEntry(hi_bg_z_e,"Vz of Bg electrons","l");
     leg->AddEntry(hi_bg_z_g,"Vz of Bg photons","l");
     leg->AddEntry(hi_bg_z_o,"Vz of Bg other","l");
     leg->Draw();
    c3->Print("dc_occ.pdf");
    
    TCanvas *c31=new TCanvas("c31","Background Origin",750,500);
    gPad->SetLogz();
    hi_bg_origin->Draw("COLZ");
    c3->Print("dc_occ_origin.pdf");

    TCanvas *c4=new TCanvas("c4","Background Energy",750,1000);
    c4->Divide(1,2);
    c4->cd(1);
    gPad->SetLogz();
    hi_bg_energy->Draw("COLZ");
    c4->cd(2);
    hi_bg_E->Draw("");
    c4->Print("dc_occ.pdf");
    
    TCanvas *c5=new TCanvas("c5","Region 1 Background Origin",750,1000);
    c5->Divide(1,3);
    c5->cd(1);
    gPad->SetLogz();
    hi_bg_y_vs_x_reg[0]->Draw("COLZ");
    c5->cd(2);
    gPad->SetLogz();
    hi_bg_y_vs_x_reg[1]->Draw("COLZ");
    c5->cd(3);
    gPad->SetLogz();
    hi_bg_y_vs_x_reg[2]->Draw("COLZ");
    c5->Print("dc_occ.pdf");
    
/*
    TCanvas *c8=new TCanvas("c8","Occupancy",750,500);
    hi_dcocc_ecut->Draw("COLZ");
    c8->Print("dc_occ_map.pdf");
    c8->Print("dc_occ.pdf");
*/
    
    TCanvas *c9=new TCanvas("c9","Sector Occupancy",500,500);
    double region_ave_occ[3]={0};
    double region_max_occ[3]={0};
    for(int i=0; i<6; i++){
        for(int iy=0; iy<hi_dcocc[i]->GetYaxis()->GetNbins(); iy++) {
            int dc_reg=int(iy/12);
            for(int ix=0; ix<hi_dcocc[i]->GetXaxis()->GetNbins(); ix++) {
                region_ave_occ[dc_reg]+=hi_dcocc[i]->GetBinContent(ix+1,iy+1);
                if(hi_dcocc[i]->GetBinContent(ix+1,iy+1)>region_max_occ[dc_reg]) region_max_occ[dc_reg]=hi_dcocc[i]->GetBinContent(ix+1,iy+1);
            }
        }    
    }
    TF1 *mypol0[3];
    for(int i=0; i<3; i++) {
        hi_dcocc_region[i]->SetMarkerStyle(20+i);
        hi_dcocc_region[i]->SetMinimum(0.);
        hi_dcocc_region[i]->SetMaximum(6);
        mypol0[i]= new TF1(Form("mypol%i",i),"pol0",-0.5,6.5);
        mypol0[i]->SetLineWidth(1);
        mypol0[i]->SetLineColor(1);
    }
    hi_dcocc_region[0]->SetMarkerColor(kBlue+2);
    hi_dcocc_region[0]->SetLineColor(kBlue+2);
    hi_dcocc_region[0]->Draw("PE");
    mypol0[0]->SetLineColor(kBlue+2);
    hi_dcocc_region[0]->Fit("mypol0","SAME");
    hi_dcocc_region[1]->SetMarkerColor(kRed+2);
    hi_dcocc_region[1]->SetLineColor(kRed+2);
    hi_dcocc_region[1]->Draw("SAMEPE");
    mypol0[1]->SetLineColor(kRed+2);
    hi_dcocc_region[1]->Fit("mypol1","SAME");
    hi_dcocc_region[2]->SetMarkerColor(kGreen+2);
    hi_dcocc_region[2]->SetLineColor(kGreen+2);
    hi_dcocc_region[2]->Draw("SAMEPE");
    mypol0[2]->SetLineColor(kGreen+2);
    hi_dcocc_region[2]->Fit("mypol2","SAME");
    TLegend *l_occ=new TLegend(0.5,0.75,0.96,0.96);
    for(int i=0; i<3; i++) {
      l_occ->AddEntry(hi_dcocc_region[i],Form("region %i occupancy: %.2f %%",i+1,mypol0[i]->GetParameter(0)),"p");
    }
    l_occ->SetTextFont(72);
    l_occ->Draw();
    c9->Print("dc_region_occ.pdf");
    c9->Print("dc_occ.pdf)");
    FILE *fp = fopen("dc_occ.txt","w");
    for(int i=0; i<3;i++) {
        region_ave_occ[i]=region_ave_occ[i]/12.0/hi_dcocc[i]->GetXaxis()->GetNbins();
	fprintf(fp,"%d   %5.3f  %5.3f  %5.3f  \n",i+1,100*region_ave_occ[i],100*region_max_occ[i],mypol0[i]->GetParameter(0));
    }
    fclose(fp);    

    // saving histograms to file
	TIter next(gDirectory->GetList());
	TObject* obj;
	TH1 *my_hi[1000];
	int ihi=0;
	while(obj= ( TObject* ) next()){
		if(obj->InheritsFrom(TH1::Class())){
			my_hi[ihi] = (TH1*)obj;
			ihi++;
		}
	}
	TFile myfile("dc.root", "recreate");
	for(int i=0; i<ihi; i++) my_hi[i]->Write();
	myfile.Close();  
	
    gui.Run(1);

}   
