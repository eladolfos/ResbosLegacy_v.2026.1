#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include <cmath>

double calc_cos_theta(TLorentzVector mum, TLorentzVector mup) {
    TLorentzVector Z = mum+mup;
    double m = Z.M();
    double mT = sqrt(m*m + Z.Pt()*Z.Pt());
    double sign = Z.Pz()/std::abs(Z.Pz());
    double eplus = (mum.E() + mum.Pz())/sqrt(2.0);
    double eminus = (mum.E() - mum.Pz())/sqrt(2.0);
    double pplus = (mup.E() + mup.Pz())/sqrt(2.0);
    double pminus = (mup.E() - mup.Pz())/sqrt(2.0);

    double costheta = sign*2.0*(eplus*pminus - eminus*pplus);
    return costheta/(m*mT);
}

double calc_phi(TLorentzVector mum, TLorentzVector mup) {
    TLorentzVector boostedMu = mum;
    TLorentzVector Z = (mum+mup);
    TVector3 CSAxis, xAxis, yAxis;
    boostedMu.Boost(-Z.BoostVector());
    double ProtonMass = 0.938272;
    double BeamEnergy = 4000;
    double sign = fabs(Z.Pz())/Z.Pz();
    TLorentzVector p1, p2;
    p1.SetPxPyPzE(0, 0, sign*BeamEnergy, BeamEnergy);
    p2.SetPxPyPzE(0, 0, -sign*BeamEnergy, BeamEnergy);

    p1.Boost(-Z.BoostVector());
    p2.Boost(-Z.BoostVector());
    CSAxis = (p1.Vect().Unit()-p2.Vect().Unit()).Unit();
    yAxis = ((p1.Vect().Unit()).Cross(p2.Vect().Unit())).Unit();
    xAxis = -(yAxis.Cross(CSAxis)).Unit();

    return atan2(boostedMu.Vect().Dot(yAxis), boostedMu.Vect().Dot(xAxis));
}

void angular() {
    TFile *f = TFile::Open("resbos.root", "READ");
    TTree *tree = (TTree*)f->Get("h10");
    float Px_d1, Py_d1, Pz_d1, E_d1;
    float Px_d2, Py_d2, Pz_d2, E_d2;
    float Px_V, Py_V, Pz_V, E_V;
    float WT00;

    // Set branchess Addresses 
    tree->SetBranchAddress("Px_d1", &Px_d1);
    tree->SetBranchAddress("Py_d1", &Py_d1);
    tree->SetBranchAddress("Pz_d1", &Pz_d1);
    tree->SetBranchAddress( "E_d1",  &E_d1);
    tree->SetBranchAddress("Px_d2", &Px_d2);
    tree->SetBranchAddress("Py_d2", &Py_d2);
    tree->SetBranchAddress("Pz_d2", &Pz_d2);
    tree->SetBranchAddress( "E_d2",  &E_d2);
    tree->SetBranchAddress("Px_V", &Px_V);
    tree->SetBranchAddress("Py_V", &Py_V);
    tree->SetBranchAddress("Pz_V", &Pz_V);
    tree->SetBranchAddress( "E_V",  &E_V);
    tree->SetBranchAddress("WT00", &WT00);

    // Angular histograms
    TFile *fOut = TFile::Open("angular_kf.root", "RECREATE");
    std::vector<double> bins = {1E-1, 2.5, 5, 8, 11.4, 14.9, 18.5, 22, 25.5, 29, 32.6, 36.4, 40.4,
                                44.9, 50.2, 56.4, 63.9, 73.4, 85.4, 105, 132, 173, 253, 600};
    TH1F *h_pt = new TH1F("h_pt", "pt", bins.size()-1, &bins[0]);
    TH1F *h_A0 = new TH1F("h_A0", "A0", bins.size()-1, &bins[0]);
    TH1F *h_A02 = new TH1F("h_A02", "A0", bins.size()-1, &bins[0]);
    TH1F *h_A1 = new TH1F("h_A1", "A1", bins.size()-1, &bins[0]);
    TH1F *h_A2 = new TH1F("h_A2", "A2", bins.size()-1, &bins[0]);
    TH1F *h_A3 = new TH1F("h_A3", "A3", bins.size()-1, &bins[0]);
    TH1F *h_A4 = new TH1F("h_A4", "A4", bins.size()-1, &bins[0]);
    TH1F *h_A5 = new TH1F("h_A5", "A5", bins.size()-1, &bins[0]);
    

    size_t nentries = (size_t)tree->GetEntries();
    std::cout << nentries << std::endl;
    // Fill the histograms
    TLorentzVector electron, positron, zboson;
    for (size_t i = 0; i < nentries; i++){
        tree->GetEntry(i);
        if(i % 1000000 == 0) {
            std::cout << "\r" << i << std::flush;
        }

        electron.SetPxPyPzE(Px_d2, Py_d2, Pz_d2, E_d2);
        positron.SetPxPyPzE(Px_d1, Py_d1, Pz_d1, E_d1);
        zboson.SetPxPyPzE(Px_V, Py_V, Pz_V, E_V);

        if((zboson.M() > 80) && (zboson.M() < 100)) {
            double ptZ = zboson.Pt();
            
            double costheta = calc_cos_theta(electron, positron);
            if(costheta > 1)
                costheta = 1;
            else if(costheta < -1)
                costheta = -1;
            double sintheta = sqrt(1-costheta*costheta);
            double phi = calc_phi(electron, positron);
            if(phi < 0) phi += 2*M_PI;

            const double A0 = (10.0/3.0*(1.0-3.0*costheta*costheta)+2.0/3.0)*WT00;
            const double A02 = 10*(-(1.0+costheta*costheta)+7./5.)*WT00; //(10.0/3.0*(1.0-3.0*costheta*costheta)+2.0/3.0)*WT00;
            const double A1 = 10*sintheta*costheta*cos(phi)*WT00;
            const double A2 = 10*sintheta*sintheta*cos(2*phi)*WT00;
            const double A3 = 4*sintheta*cos(phi)*WT00;
            const double A4 = 4*costheta*WT00;

            h_pt -> Fill(ptZ, WT00);
            h_A0 -> Fill(ptZ, A0);
            h_A02 -> Fill(ptZ, A02);
            h_A1 -> Fill(ptZ, A1);
            h_A2 -> Fill(ptZ, A2);
            h_A3 -> Fill(ptZ, A3);
            h_A4 -> Fill(ptZ, A4);
        }
    }

    //for(size_t i = 1; i <= h_pt -> GetNbinsX(); ++i) {
    //    h_pt -> SetBinContent(i, h_pt -> GetBinContent(i)/h_pt -> GetBinWidth(i));
    //    h_A0 -> SetBinContent(i, h_A0 -> GetBinContent(i)/h_A0 -> GetBinWidth(i));
    //    h_A1 -> SetBinContent(i, h_A1 -> GetBinContent(i)/h_A1 -> GetBinWidth(i));
    //    h_A2 -> SetBinContent(i, h_A2 -> GetBinContent(i)/h_A2 -> GetBinWidth(i));
    //    h_A3 -> SetBinContent(i, h_A3 -> GetBinContent(i)/h_A3 -> GetBinWidth(i));
    //    h_A4 -> SetBinContent(i, h_A4 -> GetBinContent(i)/h_A4 -> GetBinWidth(i));
    //}

    h_A0 -> Divide(h_pt);
    h_A02 -> Divide(h_pt);
    h_A1 -> Divide(h_pt);
    h_A2 -> Divide(h_pt);
    h_A3 -> Divide(h_pt);
    h_A4 -> Divide(h_pt);

    std::ofstream outfile("A0.txt");
    outfile << "pt, pt low, pt high, A0, stat+, stat-" << std::endl;
    for(size_t i = 1; i <= h_A0 -> GetNbinsX(); ++i) {
        outfile << h_A0 -> GetBinCenter(i) << "," << h_A0 -> GetBinLowEdge(i) << ",";
        outfile << h_A0 -> GetBinLowEdge(i) + h_A0 -> GetBinWidth(i) << ",";
        outfile << h_A0 -> GetBinContent(i) << ",";
        outfile << -h_A0 -> GetBinError(i) << "," << h_A0 -> GetBinError(i) << std::endl;
    }
    outfile.close();

    outfile.open("A1.txt");
    outfile << "pt, pt low, pt high, A1, stat+, stat-" << std::endl;
    for(size_t i = 1; i <= h_A1 -> GetNbinsX(); ++i) {
        outfile << h_A1 -> GetBinCenter(i) << "," << h_A1 -> GetBinLowEdge(i) << ",";
        outfile << h_A1 -> GetBinLowEdge(i) + h_A1 -> GetBinWidth(i) << ",";
        outfile << h_A1 -> GetBinContent(i) << ",";
        outfile << -h_A1 -> GetBinError(i) << "," << h_A1 -> GetBinError(i) << std::endl;
    }
    outfile.close();

    outfile.open("A2.txt");
    outfile << "pt, pt low, pt high, A2, stat+, stat-" << std::endl;
    for(size_t i = 1; i <= h_A2 -> GetNbinsX(); ++i) {
        outfile << h_A2 -> GetBinCenter(i) << "," << h_A2 -> GetBinLowEdge(i) << ",";
        outfile << h_A2 -> GetBinLowEdge(i) + h_A2 -> GetBinWidth(i) << ",";
        outfile << h_A2 -> GetBinContent(i) << ",";
        outfile << -h_A2 -> GetBinError(i) << "," << h_A2 -> GetBinError(i) << std::endl;
    }
    outfile.close();

    outfile.open("A3.txt");
    outfile << "pt, pt low, pt high, A3, stat+, stat-" << std::endl;
    for(size_t i = 1; i <= h_A3 -> GetNbinsX(); ++i) {
        outfile << h_A3 -> GetBinCenter(i) << "," << h_A3 -> GetBinLowEdge(i) << ",";
        outfile << h_A3 -> GetBinLowEdge(i) + h_A3 -> GetBinWidth(i) << ",";
        outfile << h_A3 -> GetBinContent(i) << ",";
        outfile << -h_A3 -> GetBinError(i) << "," << h_A3 -> GetBinError(i) << std::endl;
    }
    outfile.close();

    outfile.open("A4.txt");
    outfile << "pt, pt low, pt high, A4, stat+, stat-" << std::endl;
    for(size_t i = 1; i <= h_A4 -> GetNbinsX(); ++i) {
        outfile << h_A4 -> GetBinCenter(i) << "," << h_A4 -> GetBinLowEdge(i) << ",";
        outfile << h_A4 -> GetBinLowEdge(i) + h_A4 -> GetBinWidth(i) << ",";
        outfile << h_A4 -> GetBinContent(i) << ",";
        outfile << -h_A4 -> GetBinError(i) << "," << h_A4 -> GetBinError(i) << std::endl;
    }
    outfile.close();

    fOut -> Write();
}
