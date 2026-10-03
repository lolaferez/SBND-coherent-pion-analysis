// Definition of the variables that will be used to access the data in the TTree.
#include "tree_utils.cpp"
// Necessary includes to make this work
#include "Includes.h"
#include <cmath>

double deltaPhiTransverso(double px1, double py1,
                          double px2, double py2) {

    double phi1 = atan2(py1, px1);
    double phi2 = atan2(py2, px2);

    double dphi = phi1 - phi2;

    while (dphi <= -TMath::Pi()) {
        dphi += 2.0*TMath::Pi();
    }

    while (dphi > TMath::Pi()) {
        dphi -= 2.0*TMath::Pi();
    }

    return fabs(dphi);
}

void angulogenvistas2(){
    TFile *input_file;               // ROOT file
    TDirectoryFile *tree_dir;     // Directory of the data inside the ROOT file
    TTree *subrun_tree;
    TTree *event_tree;
    std::vector<string> filenames;
    filenames.push_back("./Data/SiFSI/COHpi_SiFSI.root");
    filenames.push_back("./Data/NoFSI/COHpi_NoFSI.root");

    double pxp, pyp, pzp, pxn, pyn, pzn, pxpp, pypp, pzpp;

    double p_pion, p_nu, p_muon, tetha1, tetha2, tetha, thetaori,aper;

   TH1F *tetha_Si = new TH1F("Efpi_SiFSI", "", 100, 0,180.3);
   TH1F *tetha_No = new TH1F("Efpi_NoFSI", "", 100, 0, 180.3);
   
   for(string file : filenames){
        // Open the file
        // Declare varibales as null pointers to avoid memory leaks
        input_file = nullptr;
        tree_dir = nullptr;
        subrun_tree = nullptr;
        event_tree = nullptr;

        input_file = new TFile(file.c_str());
        // Get the directory of the data inside the ROOT file. In this case, it is called "ana"
        tree_dir = (TDirectoryFile*) input_file->Get("ana");

        // Retrieve the TTrees from the ROOT file in its directory "ana"
        event_tree = (TTree*)tree_dir->Get("tree");
        subrun_tree = (TTree*)tree_dir->Get("subrun_tree");

        // Use the code from tree_utils.cpp to define all variables and set the branches of the TTree to access the data
        set_branch(event_tree);
        set_branch_subtree(subrun_tree);

        //////////////// DATA HAS BEEN LOADED ALREADY ////////////////





    
        int num_events = event_tree->GetEntries();
        for(int i=0; i<num_events; i++){
            event_tree->GetEntry(i);
            int num_gen_particles = gen_part_trackID->size();
            
            if(nu_CC_NC == 0){
                if(num_gen_particles != 1){

                    for(int j=0; j<num_gen_particles; j++){

                            //busco los piones cargados estado final
                            if(gen_part_statusCode->at(j)==1 && (abs(gen_part_PDGcode->at(j))==211)){
                                pxp = gen_part_P0_X->at(j);
                                pyp = gen_part_P0_Y->at(j);
                                pzp = gen_part_P0_Z->at(j);
                            }
                            //localizo los neutrinos inicial
                            if(gen_part_statusCode->at(j)==0 && (abs(gen_part_PDGcode->at(j))==14)){
                                pxn = gen_part_P0_X->at(j);
                                pyn = gen_part_P0_Y->at(j);
                                pzn = gen_part_P0_Z->at(j);

                            }
                            //guardo los datos del vector del lepton producto
                            if(gen_part_statusCode->at(j)==1 && (abs(gen_part_PDGcode->at(j))==13)){
                                pxpp = gen_part_P0_X->at(j);
                                pypp = gen_part_P0_Y->at(j);
                                pzpp = gen_part_P0_Z->at(j);
                            }
                    }

                    //calculo las normas
                    p_pion=sqrt(pxp*pxp+pyp*pyp+pzp*pzp);
                    p_nu=sqrt(pxn*pxn+pyn*pyn+pzn*pzn);
                    p_muon = sqrt(pxpp*pxpp+pzpp*pzpp+pypp*pypp);

                    double aperT = deltaPhiTransverso(pxp, pzp, pxpp, pzpp);

                    aperT = (aperT) * 360 / (2*(TMath::Pi())); //paso a grados 
                    
                    if(file.find("SiFSI") != string::npos){ // If the filename contains "SiFSI", we fill the histogram with FSI
                        tetha_Si->Fill(aperT);
                    }
                    else if(file.find("NoFSI") != string::npos){ // If the filename contains "NoFSI", we fill the histogram without FSI
                        tetha_No->Fill(aperT);
                    }

                }

            }
        }
    }

TCanvas *c1 = new TCanvas("c1", "", 900, 700);
// Escala logaritmica en el eje Y
c1->SetLogy();

c1->SetLeftMargin(0.15);
c1->SetBottomMargin(0.15);
c1->SetRightMargin(0.05);
c1->SetTopMargin(0.08);

tetha_No->SetStats(0);
tetha_No->GetXaxis()->SetTitleSize(0.043);
tetha_No->GetYaxis()->SetTitleSize(0.043);
tetha_No->GetXaxis()->SetLabelSize(0.04);
tetha_No->GetYaxis()->SetLabelSize(0.04);
tetha_No->GetXaxis()->SetTitle("#theta transversal [grados]");
tetha_No->GetYaxis()->SetTitle("Total eventos");
tetha_No->GetYaxis()->SetTitleOffset(1.47);

gStyle->SetTitleFontSize(0.07);

tetha_No->SetLineStyle(1);
tetha_No->SetLineWidth(2);
tetha_No->SetLineColor(kRed);
tetha_No->SetFillStyle(1001);
tetha_No->SetMarkerColor(kRed);
tetha_No->SetMarkerStyle(21);
tetha_No->SetMarkerSize(2.);

tetha_Si->SetStats(0);
tetha_Si->GetXaxis()->SetTitleSize(0.04);
tetha_Si->GetYaxis()->SetTitleSize(0.04);
tetha_Si->GetXaxis()->SetLabelSize(0.03);
tetha_Si->GetYaxis()->SetLabelSize(0.03);
tetha_Si->GetXaxis()->SetTitle("Momento del #nu_{$mu} incidente");
tetha_Si->GetYaxis()->SetTitle("Total Eventos");
tetha_Si->GetYaxis()->SetTitleOffset(1.2);

tetha_Si->SetLineStyle(1);
tetha_Si->SetLineWidth(2);
tetha_Si->SetLineColor(kBlue);
tetha_Si->SetFillStyle(1001);
tetha_Si->SetMarkerColor(kBlue);
tetha_Si->SetMarkerStyle(21);
tetha_Si->SetMarkerSize(2.);

// Draw the histograms
tetha_No->Draw("HIST");
tetha_Si->Draw("SAME");

// Estadísticos
double mean_No = tetha_No->GetMean();
double sigma_No = tetha_No->GetStdDev();

double mean_Si = tetha_Si->GetMean();
double sigma_Si = tetha_Si->GetStdDev();

// Leyenda
TLegend *legend = new TLegend(0.6, 0.72, 0.9, 0.88);

legend->AddEntry(tetha_No, "NO FSI", "l");
legend->AddEntry((TObject*)0,
                 Form("#mu = %.2f, #sigma = %.2f", mean_No, sigma_No),
                 "");

legend->AddEntry(tetha_Si, "SI FSI", "l");
legend->AddEntry((TObject*)0,
                 Form("#mu = %.2f, #sigma = %.2f", mean_Si, sigma_Si),
                 "");

// Aspecto de la caja
legend->SetBorderSize(1);   // borde visible
legend->SetFillStyle(1001); // fondo blanco
legend->SetTextSize(0.035);

legend->Draw();

}    

            