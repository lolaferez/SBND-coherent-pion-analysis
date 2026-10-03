#include "tree_utils.cpp"
#include "Includes.h"
#include <cmath>

// CODIGO PARA VER LA CANTIDAD DE PARTICULAS Y SUS PDGS EN GENIE
// SOLO EVENTOS DE CORRIENTE CARGADA

void energiasgenie(){
    TFile *input_file;            // ROOT file
    TDirectoryFile *tree_dir;     // Directory of the data inside the ROOT file
    TTree *subrun_tree;
    TTree *event_tree;
    std::vector<string> filenames;
    filenames.push_back("./Data/SiFSI/COHpi_SiFSI.root");
    filenames.push_back("./Data/NoFSI/COHpi_NoFSI.root");
    double pdg, pdgmom;

    bool pioni = false; 
    double p_pion, p_nu, tetha;
    int tot = 0, n=0, m=0, elec=0, pion=0, muoni=0, prot=0, neu=0, rara=0, nu1=0, nu2=0, nu3=0, nu4=0;

   TH1F *tetha_Si = new TH1F("Efpi_SiFSI", "", 70, 0, 1500);
   TH1F *tetha_No = new TH1F("Efpi_NoFSI", "", 70, 0, 1500);
   
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
                    
                    for(int k=0; k<num_gen_particles; k++){
                         

                        if(gen_part_statusCode->at(k)==1 && abs(gen_part_PDGcode->at(k))==211){

                            pioni = true;
                        }
                           

                    }

                    if(pioni == true){

                        for(int k=0; k<num_gen_particles; k++){

                            if(gen_part_statusCode->at(k)==1 && abs(gen_part_PDGcode->at(k))==2112){
                                double Tn=gen_part_E0->at(k)-gen_part_mass->at(k);
                                Tn=Tn*1000;

                                tetha_Si->Fill(Tn);
                                
                            }

                            if(gen_part_statusCode->at(k)==1 && abs(gen_part_PDGcode->at(k))==2212){
                                double Tp=gen_part_E0->at(k)-gen_part_mass->at(k);
                                Tp=Tp*1000;

                                tetha_No->Fill(Tp);
                            }

                        }
                    }

                    pioni = false;

                    

                      

                }

            }
               
        }

    }



if (tetha_No->Integral() > 0) {
    tetha_No->Scale(1.0 / tetha_No->Integral());
}

if (tetha_Si->Integral() > 0) {
    tetha_Si->Scale(1.0 / tetha_Si->Integral());
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
tetha_No->GetXaxis()->SetTitle("T [MeV] con #pi no absorbido");
//tetha_No->GetYaxis()->SetTitle("Total eventos");
tetha_No->GetYaxis()->SetTitleOffset(1.47);

gStyle->SetTitleFontSize(0.07);

tetha_No->SetLineStyle(1);
tetha_No->SetLineWidth(2);
tetha_No->SetLineColor(kOrange+7);
tetha_No->SetFillStyle(1001);
tetha_No->SetMarkerColor(kOrange+7);
tetha_No->SetMarkerStyle(21);
tetha_No->SetMarkerSize(2.);

tetha_Si->SetStats(0);
tetha_Si->GetXaxis()->SetTitleSize(0.04);
tetha_Si->GetYaxis()->SetTitleSize(0.04);
tetha_Si->GetXaxis()->SetLabelSize(0.03);
tetha_Si->GetYaxis()->SetLabelSize(0.03);
tetha_Si->GetXaxis()->SetTitle("T [MeV] con #pi no absorbido");
//tetha_Si->GetYaxis()->SetTitle("Total Eventos");
tetha_Si->GetYaxis()->SetTitleOffset(1.2);

tetha_Si->SetLineStyle(1);
tetha_Si->SetLineWidth(2);
tetha_Si->SetLineColor(kGreen+2);
tetha_Si->SetFillStyle(1001);
tetha_Si->SetMarkerColor(kGreen+2);
tetha_Si->SetMarkerStyle(21);
tetha_Si->SetMarkerSize(2.);

tetha_No->GetYaxis()->SetTitle("Eventos normalizados");
tetha_Si->GetYaxis()->SetTitle("Eventos normalizados");



// Para que el eje Y cubra ambos histogramas

double max_No = tetha_No->GetMaximum();

double max_Si = tetha_Si->GetMaximum();

double max_y = std::max(max_No, max_Si);

tetha_No->SetMaximum(1.2 * max_y);
// Draw the histograms
tetha_Si->Draw("HIST");
tetha_No->Draw("SAME");

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
TLegend *legend = new TLegend(0.5, 0.72, 0.93, 0.88);

legend->AddEntry(tetha_No, "SI FSI: protones", "l");
legend->AddEntry((TObject*)0,
                 Form("#mu = %.2f, #sigma = %.2f, M = %.2f",
                      mean_No, sigma_No, moda_No),
                 "");

legend->AddEntry(tetha_Si, "SI FSI: neutrones", "l");
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

            