#include "tree_utils.cpp"
#include "Includes.h"
#include <cmath>

struct Vec3 {

    double x;

    double y;

    double z;

};

double dot(const Vec3& a, const Vec3& b) {

    return a.x*b.x + a.y*b.y + a.z*b.z;

}

Vec3 cross(const Vec3& a, const Vec3& b) {

    Vec3 c;

    c.x = a.y*b.z - a.z*b.y;

    c.y = a.z*b.x - a.x*b.z;

    c.z = a.x*b.y - a.y*b.x;

    return c;

}

double norm(const Vec3& a) {

    return sqrt(dot(a,a));

}

Vec3 normalize(const Vec3& a) {

    double n = norm(a);

    Vec3 u;

    u.x = a.x/n;

    u.y = a.y/n;

    u.z = a.z/n;

    return u;

}

void construirBaseConZ(Vec3 v, Vec3& ex, Vec3& ey, Vec3& ez) {

    // Nuevo eje Z: direccion del vector v
    ez = normalize(v);

    // Vector auxiliar no paralelo a ez
    Vec3 aux;

    if (fabs(ez.z) < 0.9) {
        aux = {0.0, 0.0, 1.0};
    } else {
        aux = {1.0, 0.0, 0.0};
    }

    // Nuevo eje X perpendicular a ez
    ex = cross(aux, ez);
    ex = normalize(ex);

    // Nuevo eje Y para completar base ortonormal derecha
    ey = cross(ez, ex);
    ey = normalize(ey);
}

Vec3 cambiarCoordenadas(Vec3 p, Vec3 ex, Vec3 ey, Vec3 ez) {

    Vec3 p_new;

    p_new.x = dot(p, ex);
    p_new.y = dot(p, ey);
    p_new.z = dot(p, ez);

    return p_new;
}

// Angulo de apertura entre dos vectores proyectados en el plano XZ

double anguloAperturaXZ(double px1, double pz1,

                        double px2, double pz2) {

    // Producto escalar en el plano XZ

    double dot = px1*px2 + pz1*pz2;

    // Modulos de las proyecciones en XZ

    double mod1 = sqrt(px1*px1 + pz1*pz1);

    double mod2 = sqrt(px2*px2 + pz2*pz2);

    // Evitar division por cero

    if (mod1 == 0 || mod2 == 0) {

        return -999.0;

    }

    double cosang = dot / (mod1 * mod2);

    // Proteger por errores numericos

    if (cosang >  1.0) cosang =  1.0;

    if (cosang < -1.0) cosang = -1.0;

    // Angulo entre 0 y pi

    return acos(cosang);

}

void angulocambiosist(){
    TFile *input_file;            // ROOT file
    TDirectoryFile *tree_dir;     // Directory of the data inside the ROOT file
    TTree *subrun_tree;
    TTree *event_tree;
    std::vector<string> filenames;
    filenames.push_back("./Data/SiFSI/COHpi_SiFSI.root");
    filenames.push_back("./Data/NoFSI/COHpi_NoFSI.root");
    double px0, py0, pz0, px1, py1, pz1, px2, py2, pz2;
    double px3, py3, pz3, px4, py4, pz4;

    bool corriente = true; //para seleccionar si quiero un angulo en cc (true) o cn (false)
    double pt, pt0, pt1, pt2, ptz=0, pty=0, ptx=0;
    int n=0, m=0;

   TH1F *tetha_Si = new TH1F("Efpi_SiFSI", "", 65, 0, 180);
   TH1F *tetha_No = new TH1F("Efpi_NoFSI", "", 65, 0, 180);
   
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
            px0 = py0 = pz0 = 0.0;

            px1 = py1 = pz1 = 0.0;

            px2 = py2 = pz2 = 0.0;

            int num_gen_particles = gen_part_trackID->size();

            if(num_gen_particles != 1){ //para quitarme simulaciones mal generadas
                //selecciono corriente cargada
                if(nu_CC_NC == 0){
                     
                    for(int j=0; j<num_gen_particles; j++){

                        if( abs(gen_part_PDGcode->at(j)) == 14 && gen_part_statusCode->at(j)==0){

                            px0= gen_part_P0_X->at(j);
                            py0= gen_part_P0_Y->at(j);
                            pz0= gen_part_P0_Z->at(j);

                        }
                        if(gen_part_statusCode->at(j)==1 && abs(gen_part_PDGcode->at(j)) == 211 ){

                            px1= gen_part_P0_X->at(j);
                            py1= gen_part_P0_Y->at(j);
                            pz1= gen_part_P0_Z->at(j);

                        }

                         if(gen_part_statusCode->at(j)!=0 && gen_part_statusCode->at(j)!=(-1) && abs(gen_part_PDGcode->at(j)) == 13 ){

                            px2= gen_part_P0_X->at(j);
                            py2= gen_part_P0_Y->at(j);
                            pz2= gen_part_P0_Z->at(j);

                        }
                        

                    }

                    Vec3 pnu = {px0, py0, pz0};
                    Vec3 ex, ey, ez;
                    construirBaseConZ(pnu, ex, ey, ez);

                    Vec3 ppi = {px1, py1, pz1};
                    Vec3 ppi_new = cambiarCoordenadas(ppi, ex, ey, ez);
                    double pxp_new = ppi_new.x;
                    double pyp_new = ppi_new.y;
                    double pzp_new = ppi_new.z;

                    Vec3 pnu_new = cambiarCoordenadas(pnu, ex, ey, ez);
                    double pnux_new = pnu_new.x;
                    double pnuy_new = pnu_new.y;
                    double pnuz_new = pnu_new.z;
                    
                    Vec3 pmu = {px2, py2, pz2};
                    Vec3 pmu_new = cambiarCoordenadas(pmu, ex, ey, ez);
                    double pmux_new = pmu_new.x;
                    double pmuy_new = pmu_new.y;
                    double pmuz_new = pmu_new.z;

                    double apertura_XZ = anguloAperturaXZ(pxp_new, pyp_new,
                                      pmux_new, pmuy_new);

                    if (apertura_XZ != -999.0) {

                        double apertura_XZ_deg = apertura_XZ * 180.0 / TMath::Pi();

                        if(file.find("SiFSI") != string::npos){
                            tetha_Si->Fill(apertura_XZ_deg);
                        }
                        else if(file.find("NoFSI") != string::npos){
                            tetha_No->Fill(apertura_XZ_deg);
                        }
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
tetha_No->GetXaxis()->SetTitle("#theta_{xy}^{#mu#pi} [grados]");
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
tetha_Si->GetXaxis()->SetTitle("#theta_{z}^{#mu#pi} [grados]");
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
TLegend *legend = new TLegend(0.25, 0.72, 0.67, 0.88);

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

            