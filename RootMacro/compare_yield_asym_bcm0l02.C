#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <iostream>

void compare_yield_asym_bcm0l02()
{
    gStyle->SetOptStat(1110);

    TFile *fMoller = TFile::Open("../RootFiles/isu_sample_1267.000.root");
    TFile *fVQWK   = TFile::Open("/chafs2/work1/parity/japanOutput/prexinj_20941.000.root");

    if (!fMoller || fMoller->IsZombie()) {
        std::cout << "Could not open MOLLER file" << std::endl;
        return;
    }

    if (!fVQWK || fVQWK->IsZombie()) {
        std::cout << "Could not open VQWK file" << std::endl;
        return;
    }

    TTree *tMoller = (TTree*) fMoller->Get("mul");
    TTree *tVQWK   = (TTree*) fVQWK->Get("mul");

    if (!tMoller || !tVQWK) {
        std::cout << "Could not find mul tree" << std::endl;
        return;
    }

    // ==============================
    // Calibration constants
    // ==============================
    double yield_gain     = -0.00098;
    double yield_pedestal = 1.35;

    double asym_gain       = 1.13;
    double asym_offset_ppm = 0.0;

    // ==========================================================
    // Yield comparison
    // ==========================================================
    TH1D *hYieldMoller = new TH1D(
        "hYieldMoller",
        "yield_bcm0l02 comparison;yield_bcm0l02;Entries",
        200, 26.5, 30.5
    );

    TH1D *hYieldVQWK = new TH1D(
        "hYieldVQWK",
        "yield_bcm0l02 comparison;yield_bcm0l02;Entries",
        200, 26.5, 30.5
    );

    tMoller->Draw(
        Form("(%f*yield_bcm0l02.hw_sum + %f)>>hYieldMoller",
             yield_gain, yield_pedestal),
        Form("(%f*yield_bcm0l02.hw_sum + %f)>26.5 && (%f*yield_bcm0l02.hw_sum + %f)<30.5",
             yield_gain, yield_pedestal,
             yield_gain, yield_pedestal),
        "goff"
    );

    tVQWK->Draw(
        "yield_bcm0l02>>hYieldVQWK",
        "yield_bcm0l02>26.5 && yield_bcm0l02<30.5",
        "goff"
    );

    hYieldMoller->SetLineColor(kRed);
    hYieldMoller->SetLineWidth(2);

    hYieldVQWK->SetLineColor(kBlue);
    hYieldVQWK->SetLineWidth(2);

    TCanvas *cYield = new TCanvas("cYield", "Yield Comparison", 1000, 700);

    hYieldMoller->Draw("hist");
    hYieldVQWK->Draw("hist same");

    TLegend *legYield = new TLegend(0.65, 0.75, 0.88, 0.88);
    legYield->AddEntry(hYieldMoller, "MOLLERADC", "l");
    legYield->AddEntry(hYieldVQWK, "VQWK", "l");
    legYield->Draw();

    cYield->SaveAs("yield_bcm0l02_moller_vqwk_compare.pdf");

    // ==========================================================
    // Asymmetry comparison in ppm
    // ==========================================================
    TH1D *hAsymMoller = new TH1D(
        "hAsymMoller",
        "asym_bcm0l02 comparison;asym_bcm0l02 ppm;Entries",
        200, -6000, 10000
    );

    TH1D *hAsymVQWK = new TH1D(
        "hAsymVQWK",
        "asym_bcm0l02 comparison;asym_bcm0l02 ppm;Entries",
        200, -6000, 10000
    );

    tMoller->Draw(
        Form("(%f*asym_bcm0l02.hw_sum*1000000 + %f)>>hAsymMoller",
             asym_gain, asym_offset_ppm),
        Form("(%f*asym_bcm0l02.hw_sum*1000000 + %f)>-6000 && (%f*asym_bcm0l02.hw_sum*1000000 + %f)<10000",
             asym_gain, asym_offset_ppm,
             asym_gain, asym_offset_ppm),
        "goff"
    );

    tVQWK->Draw(
        "(asym_bcm0l02*1000000)>>hAsymVQWK",
        "asym_bcm0l02*1000000>-6000 && asym_bcm0l02*1000000<10000",
        "goff"
    );

    hAsymMoller->SetLineColor(kRed);
    hAsymMoller->SetLineWidth(2);

    hAsymVQWK->SetLineColor(kBlue);
    hAsymVQWK->SetLineWidth(2);

    TCanvas *cAsym = new TCanvas("cAsym", "Asymmetry Comparison", 1000, 700);

    hAsymMoller->Draw("hist");
    hAsymVQWK->Draw("hist same");

    TLegend *legAsym = new TLegend(0.65, 0.75, 0.88, 0.88);
    legAsym->AddEntry(hAsymMoller, "MOLLERADC", "l");
    legAsym->AddEntry(hAsymVQWK, "VQWK", "l");
    legAsym->Draw();

    cAsym->SaveAs("asym_bcm0l02_moller_vqwk_compare.pdf");

    std::cout << "Saved output PDFs." << std::endl;
}
