#include "tree_utils.cpp"
#include "Includes.h"
#include <cmath>

//CODIGO QUE ME DA EL ANGULO EN RAD QUE FORMAN LOS MOMENTOS INICIALES/FINALES DEL PION Y LEPTON 
//PRODUCTO DE LA INTERACCION 
//CON DATOS DE GEANT4

void angulopimuG4(){
    TFile *input_file;            // ROOT file
    TDirectoryFile *tree_dir;     // Directory of the data inside the ROOT file
    TTree *subrun_tree;
    TTree *event_tree;
    std::vector<string> filenames;
    filenames.push_back("./Data/SiFSI/COHpi_SiFSI.root");
    filenames.push_back("./Data/NoFSI/COHpi_NoFSI.root");
    double motherE0GEN;
    double motherE0G4;
    double pxp, pyp, pzp, pxn, pyn, pzn = 0;

    bool corriente = true; //para seleccionar si quiero un angulo en cc (true) o cn (false)
    double p_pion, p_nu, tetha;

   TH1F *tetha_Si = new TH1F("Efpi_SiFSI", "", 100, -6.3,6.3);
   TH1F *tetha_No = new TH1F("Efpi_NoFSI", "", 100, -6.3, 6.3);
   
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
            int num_g4_particles = g4_part_trackID->size();
            
            if(num_g4_particles != 1){
                //Separo los suceso en CC y CN
                if(corriente){  //CC

                    if(nu_CC_NC == 0){
                        for(int j=0; j<num_g4_particles; j++){
                            //busco los piones cargados del estado final en el proceso primario
                            if(g4_part_process->at(j)== "primary" && (g4_part_PDGcode->at(j)==211)){
                                pxp = g4_part_Pf_X->at(j);
                                pyp = g4_part_Pf_Y->at(j);
                                pzp = g4_part_Pf_Z->at(j);
                            }
                            //localizo los neutrinos del muon producto en el estado final
                            if(g4_part_process->at(j)== "primary" && (g4_part_PDGcode->at(j)==abs(13))){
                                pxn = g4_part_Pf_X->at(j);
                                pyn = g4_part_Pf_Y->at(j);
                                pzn = g4_part_Pf_Z->at(j);
                            }

                            //calculo las normas
                            p_pion=sqrt(pxp*pxp+pyp*pyp+pzp*pzp);
                            p_nu=sqrt(pxn*pxn+pyn*pyn+pzn*pzn);

                        }

                    //muestro las normas en pantalla
                    cout << "NORMA P PION: " << p_pion << endl;
                    cout << "NORMA P NU: " << p_nu << endl;

                        //calculo del angulo que forman
                        tetha = acos((pxp*pxn + pyp*pyn+ pzp*pzn)/(p_pion * p_nu));

                        if(file.find("SiFSI") != string::npos){ // If the filename contains "SiFSI", we fill the histogram with FSI
                            tetha_Si->Fill(tetha);
                        }
                        else if(file.find("NoFSI") != string::npos){ // If the filename contains "NoFSI", we fill the histogram without FSI
                            tetha_No->Fill(tetha);
                        }
                    }

                }else{ //CN
                    if(nu_CC_NC == 1){
                        for(int j=0; j<num_g4_particles; j++){
                            //busco los piones neutros dle estado final
                            if(g4_part_process->at(j)== "primary" && (g4_part_PDGcode->at(j)==111)){
                                pxp = g4_part_Pf_X->at(j);
                                pyp = g4_part_Pf_Y->at(j);
                                pzp = g4_part_Pf_Z->at(j);
                            }
                            //localizo los neutrinos del muon producto en el estado final
                            if(g4_part_process->at(j)== "primary" && (g4_part_PDGcode->at(j)==abs(14))){
                                pxn = g4_part_Pf_X->at(j);
                                pyn = g4_part_Pf_Y->at(j);
                                pzn = g4_part_Pf_Z->at(j);
                            }

                            //calculo las normas
                            p_pion=sqrt(pxp*pxp+pyp*pyp+pzp*pzp);
                            p_nu=sqrt(pxn*pxn+pyn*pyn+pzn*pzn);

                        }

                    //muestro las normas en pantalla
                    //cout << "NORMA P PION: " << p_pion << endl;
                    //cout << "NORMA P NU: " << p_nu << endl;

                        //calculo del angulo que forman
                        tetha = acos((pxp*pxn + pyp*pyn+ pzp*pzn)/(p_pion * p_nu));

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
    }

//Calculo media y desviacion de mis histogramas
double media_No = tetha_No->GetMean();
double media_Si = tetha_Si->GetMean();
double desv_No = tetha_No->GetStdDev();
double desv_Si = tetha_Si->GetStdDev();

cout << "Media NOFSI = " << media_No << " rad, en grados:  " << (media_No*360/(2*3.141592)) << endl;
cout << "Media SIFSI = " << media_Si <<  " rad, en grados: " << (media_Si*360/(2*3.141592)) << endl;
cout << "Desv NOFSI = " << desv_No <<  " rad, en grados:  " << (desv_No*360/(2*3.141592)) << endl;
cout << "Desv SIFSI = " << desv_Si <<  " rad, en grados:  " << (desv_Si*360/(2*3.141592)) << endl;



    TCanvas *c1 = new TCanvas("c1", "", 800, 600);  

tetha_No->SetStats(0);
tetha_No->GetXaxis()->SetTitleSize(0.04);
tetha_No->GetYaxis()->SetTitleSize(0.04);
tetha_No->GetXaxis()->SetLabelSize(0.03);
tetha_No->GetYaxis()->SetLabelSize(0.03);
tetha_No->GetXaxis()->SetTitle("#theta [rad]");
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
tetha_Si->GetXaxis()->SetTitle("#theta [rad]");
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
}    

            