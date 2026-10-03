#include "tree_utils.cpp"
#include "Includes.h"
#include <cmath>


void momentogen(){
    TFile *input_file;            // ROOT file
    TDirectoryFile *tree_dir;     // Directory of the data inside the ROOT file
    TTree *subrun_tree;
    TTree *event_tree;
    std::vector<string> filenames;
    filenames.push_back("./Data/SiFSI/COHpi_SiFSI.root");
    filenames.push_back("./Data/NoFSI/COHpi_NoFSI.root");
    double px0, py0, pz0, px1, py1, pz1, px2, py2, pz2;
    double px3, py3, pz3, px4, py4, pz4;
    double Enu, Emu, Epi, qA;

    bool corriente = true; //para seleccionar si quiero un angulo en cc (true) o cn (false)
    double pt, pt0, pt1, pt2, ptz=0, pty=0, ptx=0;
    int n=0, m=0;

   TH1F *tetha_Si = new TH1F("Efpi_SiFSI", "", 80, 0, 2);
   TH1F *tetha_No = new TH1F("Efpi_NoFSI", "", 80, 0, 2);
   
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

            if(num_gen_particles != 1){ //para quitarme simulaciones mal generadas
                //selecciono corriente cargada
                if(nu_CC_NC == 0){
                     
                    for(int j=0; j<num_gen_particles; j++){

                        if( abs(gen_part_PDGcode->at(j)) == 14 && gen_part_statusCode->at(j)==0){

                            px0= gen_part_P0_X->at(j);
                            py0= gen_part_P0_Y->at(j);
                            pz0= gen_part_P0_Z->at(j);
                            Enu=gen_part_E0->at(j);

                        }
                        if(gen_part_statusCode->at(j)!=0 && gen_part_statusCode->at(j)!=(-1) && abs(gen_part_PDGcode->at(j)) == 211 ){

                            px1= gen_part_P0_X->at(j);
                            py1= gen_part_P0_Y->at(j);
                            pz1= gen_part_P0_Z->at(j);
                            Epi=gen_part_E0->at(j);

                            ptz = ptz + pz1;
                            ptx = ptx + px1;
                            pty = pty + py1;

                        }

                         if(gen_part_statusCode->at(j)!=0 && gen_part_statusCode->at(j)!=(-1) && abs(gen_part_PDGcode->at(j)) == 13 ){

                            px2= gen_part_P0_X->at(j);
                            py2= gen_part_P0_Y->at(j);
                            pz2= gen_part_P0_Z->at(j);
                            Emu=gen_part_E0->at(j);

                        }
                        

                    }

                    pt0 = sqrt(px0*px0+py0*py0+pz0*pz0);//debe ser aprox nulo
                    pt2= sqrt((px1+px2)*(px1+px2)+(py1+py2)*(py1+py2));

                    pt1 = sqrt(px1*px1+py1*py1+pz1*pz1);
                    pt=sqrt((ptx)*(ptx)+(pty*pty));


                    qA=(Enu-Epi-Emu)*(Enu-Epi-Emu)-((px0-px1-px2)*(px0-px1-px2)+(py0-py1-py2)*(py0-py1-py2)+(pz0-pz1-pz2)*(pz0-pz1-pz2));
                    

                    if(file.find("SiFSI") != string::npos){ // If the filename contains "SiFSI", we fill the histogram with FSI
                        tetha_Si->Fill(abs(qA));
                    }
                    else if(file.find("NoFSI") != string::npos){ // If the filename contains "NoFSI", we fill the histogram without FSI
                        tetha_No->Fill(abs(qA));
                    }
                    ptz= pty = ptx= 0;

                }
               
            }
        
        }

    }

TCanvas *c1 = new TCanvas("c1", "", 900, 700);
// Escala logaritmica en el eje Y
//c1->SetLogy();

c1->SetLeftMargin(0.15);
c1->SetBottomMargin(0.15);
c1->SetRightMargin(0.05);
c1->SetTopMargin(0.08);

tetha_No->SetStats(0);
tetha_No->GetXaxis()->SetTitleSize(0.043);
tetha_No->GetYaxis()->SetTitleSize(0.043);
tetha_No->GetXaxis()->SetLabelSize(0.04);
tetha_No->GetYaxis()->SetLabelSize(0.04);
tetha_No->GetXaxis()->SetTitle("| t | [GeV^{2}]");
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
tetha_Si->GetXaxis()->SetTitle("| t | [GeV^{2}]");
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

// Moda
int binModa_No = tetha_No->GetMaximumBin();
double moda_No = tetha_No->GetXaxis()->GetBinCenter(binModa_No);

int binModa_Si = tetha_Si->GetMaximumBin();
double moda_Si = tetha_Si->GetXaxis()->GetBinCenter(binModa_Si);

// Leyenda
TLegend *legend = new TLegend(0.55, 0.72, 0.93, 0.88);

legend->AddEntry(tetha_No, "NO FSI", "l");
legend->AddEntry((TObject*)0,
                 Form("#mu = %.2f, #sigma = %.2f, M = %.2f",
                      mean_No, sigma_No, moda_No),
                 "");

legend->AddEntry(tetha_Si, "SI FSI", "l");
legend->AddEntry((TObject*)0,
                 Form("#mu = %.2f, #sigma = %.2f, M = %.2f",
                      mean_Si, sigma_Si, moda_Si),
                 "");

// Aspecto de la caja
legend->SetBorderSize(1);   // borde visible
legend->SetFillStyle(1001); // fondo blanco
legend->SetTextSize(0.03);

legend->Draw();

}    

            