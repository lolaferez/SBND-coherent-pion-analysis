#include "tree_utils.cpp"
#include "Includes.h"
#include <cmath>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TProfile.h>
#include <iostream>

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

void cambiosist_scatterGRAD1(){
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

    // Cortes sobre el pion, igual que en el codigo TP:
    // - elimina piones con momento nulo
    // - elimina piones con energia mayor que 2 GeV
    const double masaPionGeV = 0.13957039;
    const double energiaMaxPionGeV = 2.0;
    const double momentoMinPionGeV = 1.0e-12;

    // Corte angular para cuantificar eventos SiFSI que quedan fuera
    const double corteAnguloDeg = 50.0;
    int nSiTotales = 0;
    int nSiFueraCorte60 = 0;

   TGraph *gr_Si = new TGraph();
   TGraph *gr_No = new TGraph();
   gr_Si->SetName("pttot_vs_pnu_SiFSI");
   gr_No->SetName("pttot_vs_pnu_NoFSI");

   int nSi = 0;
   int nNo = 0;
   
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

            ptx = pty = ptz = 0.0;

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

                    //cambio de base las coordenadas de los momentos lineales de las particulas
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

                    // Momento total del pion
                    pt0 = sqrt(pxp_new*pxp_new + pyp_new*pyp_new + pzp_new*pzp_new);

                    // Energia del pion calculada a partir de su momento total.
                    double energiaPion = sqrt(pt0*pt0 + masaPionGeV*masaPionGeV);

                    // Cortes: descarta piones sin momento reconstruido/generado y piones
                    // con energia por encima de 2 GeV.
                    if(pt0 <= momentoMinPionGeV || energiaPion > energiaMaxPionGeV){
                        continue;
                    }

                    //Momento transverso total (pion + muon)
                    pt1= sqrt((pxp_new+pmux_new)*(pxp_new+pmux_new)+(pyp_new+pmuy_new)*(pyp_new+pmuy_new));
                    pt2= pzp_new + pmuz_new;
                    
                    // Angulo entre el pion y el neutrino incidente usando TVector3::Angle()

                    TVector3 pnu_vec(pnux_new, pnuy_new, pnuz_new);
                    TVector3 ppi_vec(pxp_new, pyp_new, pzp_new);

                    //Angulo pion muon 
                    double aperT = pnu_vec.Angle(ppi_vec);

                    aperT = (aperT) * 360 / (2*(TMath::Pi())); //paso a grados 
                    

                    
                    double theta_nupi = pnu_vec.Angle(ppi_vec);  // radianes
                    double theta_nupi_deg = theta_nupi * 180.0 / TMath::Pi();

                    // Scatter: eje X = p_{#pi} total, eje Y = angulo transverso pion-muon
                    if (file.find("SiFSI") != string::npos) {
                            gr_Si->SetPoint(nSi, pt0, aperT);
                            nSi++;

                            nSiTotales++;
                            if(aperT > corteAnguloDeg){
                                nSiFueraCorte60++;
                            }
                    } else if (file.find("NoFSI") != string::npos) {
                            gr_No->SetPoint(nNo, pt0, aperT);
                            nNo++;
                    }

                    
                    
                    
    
                    
                  
                }
                    
            }
               
        }
        

    }
    


// ============================================================
// DIBUJO Y GUARDADO DE LAS TRES GRAFICAS
// ============================================================

gStyle->SetOptTitle(0);
gStyle->SetOptStat(0);

// Aspecto de los puntos, mismo formato que en TP
// NoFSI
gr_No->SetMarkerColor(kRed-7);
gr_No->SetMarkerStyle(20);
gr_No->SetMarkerSize(0.35);

// SiFSI
gr_Si->SetMarkerColor(kBlue-7);
gr_Si->SetMarkerStyle(21);
gr_Si->SetMarkerSize(0.35);

// ------------------------------------------------------------
// Rango comun de los ejes
// ------------------------------------------------------------
auto actualizarRangoY = [](TGraph *gr, double &yMin, double &yMax){
    double x, y;
    for(int i = 0; i < gr->GetN(); ++i){
        gr->GetPoint(i, x, y);
        if(y < yMin) yMin = y;
        if(y > yMax) yMax = y;
    }
};

auto mediaY = [](TGraph *gr){
    double x, y;
    double suma = 0.0;
    for(int i = 0; i < gr->GetN(); ++i){
        gr->GetPoint(i, x, y);
        suma += y;
    }
    return (gr->GetN() > 0) ? suma/gr->GetN() : 0.0;
};

auto crearPerfil = [](const char *nombre, TGraph *gr, int nBins, double xMin, double xMax){
    TProfile *prof = new TProfile(nombre, "", nBins, xMin, xMax);
    double x, y;
    for(int i = 0; i < gr->GetN(); ++i){
        gr->GetPoint(i, x, y);
        if(x >= xMin && x <= xMax){
            prof->Fill(x, y);
        }
    }
    return prof;
};

// Eje horizontal fijado hasta 2 GeV, igual que en TP.
const double xMinComun = 0.0;
const double xMaxComun = 2.0;

// Angulos en grados: se fija un rango comun natural.
double yMinComun = 0.0;
double yMaxComun = 180.0;

// Perfil promedio de cada nube de puntos.
const int nBinsPerfil = 75;
TProfile *prof_Si = crearPerfil("prof_grad_vs_ppi_SiFSI", gr_Si, nBinsPerfil, xMinComun, xMaxComun);
TProfile *prof_No = crearPerfil("prof_grad_vs_ppi_NoFSI", gr_No, nBinsPerfil, xMinComun, xMaxComun);

// Formato de los TProfile, igual que en TP.
prof_Si->SetStats(0);
prof_Si->SetLineColor(kBlue+2);
prof_Si->SetMarkerColor(kBlue+2);
prof_Si->SetLineWidth(5);
prof_Si->SetLineStyle(1);
prof_Si->SetMarkerStyle(20);
prof_Si->SetMarkerSize(1.2);

prof_No->SetStats(0);
prof_No->SetLineColor(kRed+2);
prof_No->SetMarkerColor(kRed+2);
prof_No->SetLineWidth(5);
prof_No->SetLineStyle(1);
prof_No->SetMarkerStyle(22);
prof_No->SetMarkerSize(1.2);

// -------------------------
// 1) Grafica SiFSI sola
// -------------------------
TCanvas *cSi = new TCanvas("cSi", "SiFSI", 900, 700);
cSi->SetLeftMargin(0.15);
cSi->SetBottomMargin(0.15);
cSi->SetRightMargin(0.05);
cSi->SetTopMargin(0.05);

gr_Si->SetTitle(";p_{#pi} [GeV];#theta_{#nu#pi} [#circ]");
gr_Si->SetMinimum(yMinComun);
gr_Si->SetMaximum(yMaxComun);
gr_Si->Draw("AP");
gr_Si->GetXaxis()->SetLimits(xMinComun, xMaxComun);

gr_Si->GetXaxis()->SetTitleSize(0.055);
gr_Si->GetYaxis()->SetTitleSize(0.055);
gr_Si->GetXaxis()->SetLabelSize(0.045);
gr_Si->GetYaxis()->SetLabelSize(0.045);
gr_Si->GetXaxis()->SetTitleOffset(1.10);
gr_Si->GetYaxis()->SetTitleOffset(1.25);

prof_Si->Draw("SAME HIST P");

TLegend *legSi = new TLegend(0.65, 0.7, 0.93, 0.90);
legSi->SetBorderSize(1);
legSi->SetLineColor(kBlack);
legSi->SetLineWidth(1);
legSi->SetFillColor(kWhite);
legSi->SetFillStyle(1001);
legSi->SetTextSize(0.04);
legSi->AddEntry(gr_Si, "SiFSI", "p");
legSi->AddEntry(prof_Si, "TProfile SiFSI", "lp");
legSi->Draw();

cSi->Update();
cSi->SaveAs("scatterGRADpiT_SiFSI.png");
cSi->SaveAs("scatterGRADpiT_SiFSI.pdf");

// -------------------------
// 2) Grafica NoFSI sola
// -------------------------
TCanvas *cNo = new TCanvas("cNo", "NoFSI", 900, 700);
cNo->SetLeftMargin(0.15);
cNo->SetBottomMargin(0.15);
cNo->SetRightMargin(0.05);
cNo->SetTopMargin(0.05);

gr_No->SetTitle(";p_{#pi} [GeV];#theta_{#nu#pi} [#circ]");
gr_No->SetMinimum(yMinComun);
gr_No->SetMaximum(yMaxComun);
gr_No->Draw("AP");
gr_No->GetXaxis()->SetLimits(xMinComun, xMaxComun);

gr_No->GetXaxis()->SetTitleSize(0.055);
gr_No->GetYaxis()->SetTitleSize(0.055);
gr_No->GetXaxis()->SetLabelSize(0.045);
gr_No->GetYaxis()->SetLabelSize(0.045);
gr_No->GetXaxis()->SetTitleOffset(1.10);
gr_No->GetYaxis()->SetTitleOffset(1.25);

prof_No->Draw("SAME HIST P");

TLegend *legNo = new TLegend(0.65, 0.7, 0.93, 0.90);
legNo->SetBorderSize(1);
legNo->SetLineColor(kBlack);
legNo->SetLineWidth(1);
legNo->SetFillColor(kWhite);
legNo->SetFillStyle(1001);
legNo->SetTextSize(0.04);
legNo->AddEntry(gr_No, "NoFSI", "p");
legNo->AddEntry(prof_No, "TProfile NoFSI", "lp");
legNo->Draw();

cNo->Update();
cNo->SaveAs("scatterGRADpiT_NoFSI.png");
cNo->SaveAs("scatterGRADpiT_NoFSI.pdf");

// -------------------------
// 3) Grafica comparativa
// -------------------------
TCanvas *cComp = new TCanvas("cComp", "Comparacion", 900, 700);
cComp->SetLeftMargin(0.15);
cComp->SetBottomMargin(0.15);
cComp->SetRightMargin(0.05);
cComp->SetTopMargin(0.05);

TMultiGraph *mg = new TMultiGraph();
mg->SetTitle(";p_{#pi} [GeV];#theta_{#nu#pi} [#circ]");

// En la comparativa, ROOT dibuja encima el ultimo TGraph agregado.
double yMediaSi = mediaY(gr_Si);
double yMediaNo = mediaY(gr_No);
if(yMediaSi < yMediaNo){
    mg->Add(gr_No, "P");
    mg->Add(gr_Si, "P");
} else {
    mg->Add(gr_Si, "P");
    mg->Add(gr_No, "P");
}

mg->SetMinimum(yMinComun);
mg->SetMaximum(yMaxComun);
mg->Draw("A");
mg->GetXaxis()->SetLimits(xMinComun, xMaxComun);

mg->GetXaxis()->SetTitleSize(0.055);
mg->GetYaxis()->SetTitleSize(0.055);
mg->GetXaxis()->SetLabelSize(0.045);
mg->GetYaxis()->SetLabelSize(0.045);
mg->GetXaxis()->SetTitleOffset(1.10);
mg->GetYaxis()->SetTitleOffset(1.25);

// Dibujamos tambien los perfiles.
if(yMediaSi < yMediaNo){
    prof_No->Draw("SAME HIST P");
    prof_Si->Draw("SAME HIST P");
} else {
    prof_Si->Draw("SAME HIST P");
    prof_No->Draw("SAME HIST P");
}

TLegend *legComp = new TLegend(0.65, 0.7, 0.93, 0.90);
legComp->SetBorderSize(1);
legComp->SetLineColor(kBlack);
legComp->SetLineWidth(1);
legComp->SetFillStyle(1001);
legComp->SetFillColor(kWhite);
legComp->SetTextSize(0.04);
legComp->AddEntry(gr_Si, "SiFSI", "p");
legComp->AddEntry(gr_No, "NoFSI", "p");
legComp->AddEntry(prof_Si, "TProfile SiFSI", "lp");
legComp->AddEntry(prof_No, "TProfile NoFSI", "lp");
legComp->Draw();

cComp->Update();
cComp->SaveAs("scatterGRADpiT_comparacion.png");
cComp->SaveAs("scatterGRADpiT_comparacion.pdf");

// ============================================================
// PORCENTAJE DE EVENTOS SiFSI QUE SALEN DEL CORTE A 60 GRADOS
// ============================================================
double porcentajeSiFuera60 = 0.0;
if(nSiTotales > 0){
    porcentajeSiFuera60 = 100.0 * static_cast<double>(nSiFueraCorte60) / static_cast<double>(nSiTotales);
}

std::cout << "============================================" << std::endl;
std::cout << "Eventos SiFSI usados en GRAD: " << nSiTotales << std::endl;
std::cout << "Eventos SiFSI con theta > " << corteAnguloDeg << " grados: "
          << nSiFueraCorte60 << std::endl;
std::cout << "Porcentaje SiFSI fuera del corte de " << corteAnguloDeg << " grados: "
          << porcentajeSiFuera60 << " %" << std::endl;
std::cout << "============================================" << std::endl;

}    
            
