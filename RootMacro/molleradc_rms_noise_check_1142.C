#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TH1.h>
#include <TStyle.h>
#include <TMath.h>

#include <iostream>
#include <fstream>
#include <vector>

void molleradc_rms_noise_check_1142()
{
    TString mollerFile = "../RootFiles/isu_sample_1142.001.root";

    TFile *fMoller = TFile::Open(mollerFile);

    if (!fMoller || fMoller->IsZombie()) {
        std::cout << "ERROR: Could not open file: " << mollerFile << std::endl;
        return;
    }

    TTree *evt = (TTree*) fMoller->Get("evt");

    if (!evt) {
        std::cout << "ERROR: Could not find evt tree" << std::endl;
        return;
    }

    TString ch = "tq01_r1";
    const int nBlocks = 4;

    if (!evt->GetBranch(ch)) {
        std::cout << "ERROR: Could not find channel branch " << ch << std::endl;
        return;
    }

    gStyle->SetOptStat(1110);

    TString pdfName = "molleradc_rms_noise_check_1142.pdf";
    TString csvName = "molleradc_rms_noise_check_1142_summary.csv";

    std::ofstream out(csvName.Data());

    out << "Quantity,Block,Mean,RMS,Entries,"
        << "N_gt_5sigma,Fraction_gt_5sigma,"
        << "N_gt_10sigma,Fraction_gt_10sigma\n";

    std::cout << "\nRun 1142 RMS/Noise Diagnostic\n";
    std::cout << "Channel: " << ch << "\n";
    std::cout << "Input file: " << mollerFile << "\n";
    std::cout << "Blocks: " << nBlocks << "\n\n";

    TString fullRMS = ch + ".hw_sum_rms";

    TCanvas *c = new TCanvas("c", "Run 1142 RMS/Noise Checks", 1800, 1100);

    // ------------------------------------------------------------
    // Page 1: RMS Distribution for Every Block
    // ------------------------------------------------------------
    c->Clear();
    c->Divide(2, 2);

    for (int i = 0; i < nBlocks; i++) {
        c->cd(i + 1);

        TString expr = Form("%s.block%d_rms", ch.Data(), i);
        evt->Draw(expr);

        TH1 *h = (TH1*) gPad->GetPrimitive("htemp");
        if (h) {
            h->SetTitle(Form("Block %d RMS Distribution;%s;Entries", i, expr.Data()));

            out << "BlockRMS," << i << ","
                << h->GetMean() << ","
                << h->GetRMS() << ","
                << h->GetEntries() << ","
                << "NA,NA,NA,NA\n";
        }
    }

    c->Print(pdfName + "(");

    // ------------------------------------------------------------
    // Page 2: Mean RMS vs Block Number
    // ------------------------------------------------------------
    TGraph *gMeanRMS = new TGraph();

    for (int i = 0; i < nBlocks; i++) {
        TString expr = Form("%s.block%d_rms", ch.Data(), i);
        evt->Draw(expr, "", "goff");

        int n = evt->GetSelectedRows();
        double meanRMS = 0.0;

        if (n > 0)
            meanRMS = TMath::Mean(n, evt->GetV1());

        gMeanRMS->SetPoint(i, i, meanRMS);

        std::cout << "Block " << i
                  << " mean RMS = " << meanRMS
                  << std::endl;
    }

    c->Clear();
    gMeanRMS->SetTitle(Form("%s Mean RMS vs Block Number;Block Number;Mean RMS", ch.Data()));
    gMeanRMS->SetMarkerStyle(20);
    gMeanRMS->SetMarkerSize(1.3);
    gMeanRMS->Draw("APL");

    c->Print(pdfName);

    // ------------------------------------------------------------
    // Page 3: Full Window RMS Distribution
    // ------------------------------------------------------------
    c->Clear();

    evt->Draw(fullRMS);

    TH1 *hFull = (TH1*) gPad->GetPrimitive("htemp");
    if (hFull) {
        hFull->SetTitle(Form("%s Full Window RMS Distribution;%s;Entries", ch.Data(), fullRMS.Data()));

        out << "FullWindowRMS,-1,"
            << hFull->GetMean() << ","
            << hFull->GetRMS() << ","
            << hFull->GetEntries() << ","
            << "NA,NA,NA,NA\n";
    }

    c->Print(pdfName);

    // ------------------------------------------------------------
    // Page 4: Width Distribution for Every Block
    // width = abs(RawMax_i - RawMin_i)
    // ------------------------------------------------------------
    c->Clear();
    c->Divide(2, 2);

    for (int i = 0; i < nBlocks; i++) {
        c->cd(i + 1);

        TString expr = Form("abs(%s.RawMax_%d-%s.RawMin_%d)",
                            ch.Data(), i,
                            ch.Data(), i);

        evt->Draw(expr);

        TH1 *h = (TH1*) gPad->GetPrimitive("htemp");
        if (h) {
            h->SetTitle(Form("Block %d Width = |RawMax - RawMin|;Width;Entries", i));

            out << "Width," << i << ","
                << h->GetMean() << ","
                << h->GetRMS() << ","
                << h->GetEntries() << ","
                << "NA,NA,NA,NA\n";
        }
    }

    c->Print(pdfName);

    // ------------------------------------------------------------
    // Page 5: Width vs RMS Correlation
    // ------------------------------------------------------------
    c->Clear();
    c->Divide(2, 2);

    for (int i = 0; i < nBlocks; i++) {
        c->cd(i + 1);

        TString expr = Form("abs(%s.RawMax_%d-%s.RawMin_%d):%s.block%d_rms",
                            ch.Data(), i,
                            ch.Data(), i,
                            ch.Data(), i);

        evt->Draw(expr, "", "*");

        TH1 *h = (TH1*) gPad->GetPrimitive("htemp");
        if (h) {
            h->SetTitle(Form("Block %d Width vs RMS;RMS;Width", i));
        }
    }

    c->Print(pdfName);

    // ------------------------------------------------------------
    // Page 6: Min/Max Distance from Mean
    // ------------------------------------------------------------
    c->Clear();
    c->Divide(2, 2);

    for (int i = 0; i < nBlocks; i++) {
        c->cd(i + 1);

        TString expr = Form("abs(%s.RawMax_%d-%s.block%d):abs(%s.block%d-%s.RawMin_%d)",
                            ch.Data(), i,
                            ch.Data(), i,
                            ch.Data(), i,
                            ch.Data(), i);

        evt->Draw(expr, "", "*");

        TH1 *h = (TH1*) gPad->GetPrimitive("htemp");
        if (h) {
            h->SetTitle(Form("Block %d: |RawMax-Mean| vs |Mean-RawMin|;|Mean-RawMin|;|RawMax-Mean|", i));
        }
    }

    c->Print(pdfName);

    // ------------------------------------------------------------
    // Page 7: Min/Max Distance from Mean in RMS Units
    // ------------------------------------------------------------
    c->Clear();
    c->Divide(2, 2);

    for (int i = 0; i < nBlocks; i++) {
        c->cd(i + 1);

        TString expr = Form("(abs(%s.RawMax_%d-%s.block%d)/%s.block%d_rms):(abs(%s.block%d-%s.RawMin_%d)/%s.block%d_rms)",
                            ch.Data(), i,
                            ch.Data(), i,
                            ch.Data(), i,
                            ch.Data(), i,
                            ch.Data(), i,
                            ch.Data(), i);

        evt->Draw(expr, "", "*");

        TH1 *h = (TH1*) gPad->GetPrimitive("htemp");
        if (h) {
            h->SetTitle(Form("Block %d: Min/Max Distance in RMS Units;|Mean-RawMin|/RMS;|RawMax-Mean|/RMS", i));
        }
    }

    c->Print(pdfName);

    // ------------------------------------------------------------
    // Page 8: Largest Min/Max Excursion in RMS Units
    // ------------------------------------------------------------
    c->Clear();
    c->Divide(2, 2);

    for (int i = 0; i < nBlocks; i++) {
        c->cd(i + 1);

        TString expr = Form("max(abs(%s.RawMax_%d-%s.block%d),abs(%s.block%d-%s.RawMin_%d))/%s.block%d_rms",
                            ch.Data(), i,
                            ch.Data(), i,
                            ch.Data(), i,
                            ch.Data(), i,
                            ch.Data(), i);

        evt->Draw(expr);

        TH1 *h = (TH1*) gPad->GetPrimitive("htemp");
        if (h) {
            h->SetTitle(Form("Block %d: Largest Min/Max Distance;Largest distance / RMS;Entries", i));

            double nTotal = h->GetEntries();

            double nAbove5 = h->Integral(h->FindBin(5.0), h->GetNbinsX());
            double nAbove10 = h->Integral(h->FindBin(10.0), h->GetNbinsX());

            double fAbove5 = 0.0;
            double fAbove10 = 0.0;

            if (nTotal > 0) {
                fAbove5 = nAbove5 / nTotal;
                fAbove10 = nAbove10 / nTotal;
            }

            std::cout << "Block " << i
                      << " largest excursion/RMS: mean = " << h->GetMean()
                      << "  RMS = " << h->GetRMS()
                      << "  >5sigma = " << nAbove5
                      << " (" << fAbove5 << ")"
                      << "  >10sigma = " << nAbove10
                      << " (" << fAbove10 << ")"
                      << std::endl;

            out << "LargestDistanceOverRMS," << i << ","
                << h->GetMean() << ","
                << h->GetRMS() << ","
                << h->GetEntries() << ","
                << nAbove5 << ","
                << fAbove5 << ","
                << nAbove10 << ","
                << fAbove10 << "\n";
        }
    }

    c->Print(pdfName);

    // ------------------------------------------------------------
    // Page 9: Positive-side excursion, |RawMax - Mean| / RMS
    // ------------------------------------------------------------
    c->Clear();
    c->Divide(2, 2);

    for (int i = 0; i < nBlocks; i++) {
        c->cd(i + 1);

        TString expr = Form("abs(%s.RawMax_%d-%s.block%d)/%s.block%d_rms",
                            ch.Data(), i,
                            ch.Data(), i,
                            ch.Data(), i);

        evt->Draw(expr);

        TH1 *h = (TH1*) gPad->GetPrimitive("htemp");
        if (h) {
            h->SetTitle(Form("Block %d: |RawMax-Mean| / RMS;|RawMax-Mean| / RMS;Entries", i));

            double nTotal = h->GetEntries();
            double nAbove5 = h->Integral(h->FindBin(5.0), h->GetNbinsX());
            double nAbove10 = h->Integral(h->FindBin(10.0), h->GetNbinsX());

            double fAbove5 = 0.0;
            double fAbove10 = 0.0;

            if (nTotal > 0) {
                fAbove5 = nAbove5 / nTotal;
                fAbove10 = nAbove10 / nTotal;
            }

            out << "RawMaxDistanceOverRMS," << i << ","
                << h->GetMean() << ","
                << h->GetRMS() << ","
                << h->GetEntries() << ","
                << nAbove5 << ","
                << fAbove5 << ","
                << nAbove10 << ","
                << fAbove10 << "\n";
        }
    }

    c->Print(pdfName);

    // ------------------------------------------------------------
    // Page 10: Negative-side excursion, |Mean - RawMin| / RMS
    // ------------------------------------------------------------
    c->Clear();
    c->Divide(2, 2);

    for (int i = 0; i < nBlocks; i++) {
        c->cd(i + 1);

        TString expr = Form("abs(%s.block%d-%s.RawMin_%d)/%s.block%d_rms",
                            ch.Data(), i,
                            ch.Data(), i,
                            ch.Data(), i);

        evt->Draw(expr);

        TH1 *h = (TH1*) gPad->GetPrimitive("htemp");
        if (h) {
            h->SetTitle(Form("Block %d: |Mean-RawMin| / RMS;|Mean-RawMin| / RMS;Entries", i));

            double nTotal = h->GetEntries();
            double nAbove5 = h->Integral(h->FindBin(5.0), h->GetNbinsX());
            double nAbove10 = h->Integral(h->FindBin(10.0), h->GetNbinsX());

            double fAbove5 = 0.0;
            double fAbove10 = 0.0;

            if (nTotal > 0) {
                fAbove5 = nAbove5 / nTotal;
                fAbove10 = nAbove10 / nTotal;
            }

            out << "RawMinDistanceOverRMS," << i << ","
                << h->GetMean() << ","
                << h->GetRMS() << ","
                << h->GetEntries() << ","
                << nAbove5 << ","
                << fAbove5 << ","
                << nAbove10 << ","
                << fAbove10 << "\n";
        }
    }

    c->Print(pdfName + ")");

    out.close();

    std::cout << "\nSaved PDF: " << pdfName << std::endl;
    std::cout << "Saved CSV: " << csvName << std::endl;
    std::cout << "Done.\n";
}
