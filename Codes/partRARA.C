#include "tree_utils.cpp"
#include "Includes.h"
#include <cmath>

// CODIGO PARA VER LA CANTIDAD DE PARTICULAS Y SUS PDGS EN GENIE
// SOLO EVENTOS DE CORRIENTE CARGADA

void partRARA(){
    TFile *input_file;            // ROOT file
    TDirectoryFile *tree_dir;     // Directory of the data inside the ROOT file
    TTree *subrun_tree;
    TTree *event_tree;
    std::vector<string> filenames;
    filenames.push_back("./Data/SiFSI/COHpi_SiFSI.root");
    filenames.push_back("./Data/NoFSI/COHpi_NoFSI.root");
    double pdg, pdgmom;

    bool corriente = true, pineutro =false, picargado = false, muoni= false, proton = false, neutron = false; //para seleccionar si quiero un angulo en cc (true) o cn (false)
    double px, py , pz;
    int n=0, m=0, elec=0, pion=0, muon=0, rara=0;

   TH1F *tetha_Si = new TH1F("Efpi_SiFSI", "", 100, 0, 16.3);
   TH1F *tetha_No = new TH1F("Efpi_NoFSI", "", 100, 0, 16.3);
   
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

                        if(abs(gen_part_PDGcode->at(j))==211){
                            picargado = true;
                        }

                        if(abs(gen_part_PDGcode->at(j))==2112){
                            neutron = true;
                        }

                        if(abs(gen_part_PDGcode->at(j))==2212){
                            proton = true;
                        }

                        if(abs(gen_part_PDGcode->at(j))==13){
                            muoni = true;
                        }

                        if(gen_part_PDGcode->at(j)==111){
                            pineutro = true;

                            for(int i=0; i<num_gen_particles; i++){

                                if(gen_part_statusCode->at(i)!=0 && gen_part_statusCode->at(i)!=-1 && gen_part_PDGcode->at(i) < 10000){

                                    if(file.find("SiFSI") != string::npos){ // If the filename contains "SiFSI", we fill the histogram with FSI
                                        //cout << gen_part_PDGcode->at(i) << endl;
                                    }

                                }
                            }
                            //cout << endl;
                    
                        }
                        

                    }
                    if(pineutro == true){m++;}

                    if(pineutro == true && (proton == true || neutron == true)){
                        n++;
                    }

                    pineutro = false;
                    picargado = false;
                    muoni = false;
                    proton = false;
                    neutron = false;

                }
               
            }
        
        }

        cout << " evento con piones neutros + pi cargado + muon " << n << endl;
        cout << " evento con piones cargados " << m << endl;
        n = m = 0;
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
tetha_No->GetXaxis()->SetTitle("Cantidad de p + n en SiFSI");
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
tetha_Si->GetXaxis()->SetTitle("Cantidad de p + n en SiFSI con #pi absorbido");
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
TLegend *legend = new TLegend(0.67, 0.75, 0.93, 0.88);

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

            