// Definition of the variables that will be used to access the data in the TTree.
#include "tree_utils.cpp"
// Necessary includes to make this work
#include "Includes.h"
#include <cmath>

void tethagenie(){
    TFile *input_file;               // ROOT file
    TDirectoryFile *tree_dir;     // Directory of the data inside the ROOT file
    TTree *subrun_tree;
    TTree *event_tree;
    std::vector<string> filenames;
    filenames.push_back("./Data/SiFSI/COHpi_SiFSI.root");
    filenames.push_back("./Data/NoFSI/COHpi_NoFSI.root");
    double motherE0GEN;
    double motherE0G4;
    double pxp, pyp, pzp, pxn, pyn, pzn = 0;

    double p_pion, p_nu, tetha;

   TH1F *tetha_Si = new TH1F("Efpi_SiFSI", "", 100, -1.3,100.3);
   TH1F *tetha_No = new TH1F("Efpi_NoFSI", "", 100, -1.3, 100.3);
   
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

                    //localizo los piones del estado final y guardo su momento en el instante posterior a la colision
                        if(gen_part_statusCode->at(j)==1 && (abs(gen_part_PDGcode->at(j))==13)){
                            pxp = gen_part_P0_X->at(j);
                            pyp = gen_part_P0_Y->at(j);
                            pzp = gen_part_P0_Z->at(j);
                        }
                
                        if(gen_part_statusCode->at(j)==0 && abs(gen_part_PDGcode->at(j))==14){
                            pxn = gen_part_P0_X->at(j);
                            pyn = gen_part_P0_Y->at(j);
                            pzn = gen_part_P0_Z->at(j);
                        }

                        //calculo las normas
                        p_pion=sqrt(pxp*pxp+pyp*pyp+pzp*pzp);
                        p_nu=sqrt(pxn*pxn+pyn*pyn+pzn*pzn);

                    }

                    //muestro las normas en pantalla
                    // cout << "NORMA P PION: " << p_pion << endl;
                    // cout << "NORMA P NU: " << p_nu << endl;

                    //calculo del angulo que forman
                    tetha = acos((pxp*pxn + pyp*pyn+ pzp*pzn)/(p_pion * p_nu));
                    tetha = tetha * 360 / (2*(TMath::Pi())); //paso a grados 
                    
                    if(file.find("SiFSI") != string::npos){ // If the filename contains "SiFSI", we fill the histogram with FSI
                        tetha_Si->Fill(tetha);
                    }
                    else if(file.find("NoFSI") != string::npos){ // If the filename contains "NoFSI", we fill the histogram without FSI
                        tetha_No->Fill(tetha);
                    }

                }

            }
        }
    }

// Mean and standard deviation of both histograms
std::ostringstream statsText;
statsText << std::fixed << std::setprecision(4);

double mean_No = tetha_No->GetMean();
double std_No  = tetha_No->GetStdDev();
double mean_Si = tetha_Si->GetMean();
double std_Si  = tetha_Si->GetStdDev();

TPaveText *statsBox = new TPaveText(0.55, 0.23, 0.8, 0.33, "NDC");
statsBox->SetTextSize(0.025);
statsBox->SetFillStyle(0);
statsBox->SetBorderSize(0);
statsBox->SetTextAlign(12);
statsBox->SetTextSize(0.032);

statsText.str("");
statsText << "NOFSI:  #mu = " << mean_No << ",  #sigma = " << std_No;
statsBox->AddText(statsText.str().c_str());

statsText.str("");
statsText << "SIFSI:  #mu = " << mean_Si << ",  #sigma = " << std_Si;
statsBox->AddText(statsText.str().c_str());

TCanvas *c1 = new TCanvas("c1", "", 800, 600);  

tetha_No->SetStats(0);
tetha_No->GetXaxis()->SetTitleSize(0.04);
tetha_No->GetYaxis()->SetTitleSize(0.04);
tetha_No->GetXaxis()->SetLabelSize(0.03);
tetha_No->GetYaxis()->SetLabelSize(0.03);
tetha_No->GetXaxis()->SetTitle("#theta [grados]");
tetha_No->GetYaxis()->SetTitle("# Events");
tetha_No->GetYaxis()->SetTitleOffset(1.2);

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
tetha_Si->GetXaxis()->SetTitle("#theta [grados]");
tetha_Si->GetYaxis()->SetTitle("# Events");
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

// Add a legend
TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
legend->AddEntry(tetha_No, "NOFSI DATA", "l");
legend->AddEntry(tetha_Si, "SIFSI DATA", "l");
legend->SetBorderSize(0);
legend->SetFillStyle(0);
legend->Draw("SAME");
statsBox->Draw("SAME");

c1->Update();
}    

            