#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TGraph.h>
#include <TF1.h>
#include <TCanvas.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <iomanip>
#include <cmath>

struct Stats {
    int n = 0;
    double sum = 0.0;

    void add(double x) {
        n++;
        sum += x;
    }

    double mean() const {
        if (n <= 0) return 0.0;
        return sum / n;
    }
};

void compare_all_molleradc_vqwk_channels()
{
    TString mollerFile = "../RootFiles/isu_sample_1269.000.root";
    TString vqwkFile   = "/adaqfs/home/apar/PREX/japan_feedback-SBS/japanOutput/prexinj_20942.000.root";

    TFile *fMoller = TFile::Open(mollerFile);
    TFile *fVQWK   = TFile::Open(vqwkFile);

    if (!fMoller || fMoller->IsZombie()) {
        std::cout << "ERROR: Could not open MOLLER file\n";
        return;
    }

    if (!fVQWK || fVQWK->IsZombie()) {
        std::cout << "ERROR: Could not open VQWK file\n";
        return;
    }

    TTree *mEvt = (TTree*) fMoller->Get("evt");
    TTree *vEvt = (TTree*) fVQWK->Get("evt");

    if (!mEvt || !vEvt) {
        std::cout << "ERROR: Could not find evt tree\n";
        return;
    }

    std::vector<TString> channels = {
        "bpm2i00XP","bpm2i00XM","bpm2i00YP","bpm2i00YM",
        "bpm2i00aXP","bpm2i00aXM","bpm2i00aYP","bpm2i00aYM",
        "bpm2i01aXP","bpm2i01aXM","bpm2i01aYP","bpm2i01aYM",
        "bpm1i02XP","bpm1i02XM","bpm1i02YP","bpm1i02YM",

        "bpm1i04XP","bpm1i04XM","bpm1i04YP","bpm1i04YM",
        "bpm1i05XP","bpm1i05XM","bpm1i05YP","bpm1i05YM",
        "bpm1i07XP","bpm1i07XM","bpm1i07YP","bpm1i07YM",
        "bpm0i07XP","bpm0i07XM","bpm0i07YP","bpm0i07YM",

        "bpm0l01XP","bpm0l01XM","bpm0l01YP","bpm0l01YM",
        "bpm0l02XP","bpm0l02XM","bpm0l02YP","bpm0l02YM",
        "bpm0l03XP","bpm0l03XM","bpm0l03YP","bpm0l03YM",
        "bpm0l04XP","bpm0l04XM","bpm0l04YP","bpm0l04YM",

        "bpm0l05XP","bpm0l05XM","bpm0l05YP","bpm0l05YM",
        "bpm0l06XP","bpm0l06XM","bpm0l06YP","bpm0l06YM",
        "bpm0r05XP","bpm0r05XM","bpm0r05YP","bpm0r05YM",

        "bpm0l07XP","bpm0l07XM","bpm0l07YP","bpm0l07YM",
        "bpm0l08XP","bpm0l08XM","bpm0l08YP","bpm0l08YM",
        "bpm0l09XP","bpm0l09XM","bpm0l09YP","bpm0l09YM",
        "bpm0l10XP","bpm0l10XM","bpm0l10YP","bpm0l10YM",

        "bpm0r03XP","bpm0r03XM","bpm0r03YP","bpm0r03YM",
        "bpm0i01XP","bpm0i01XM","bpm0i01YP","bpm0i01YM",
        "bpm0i01bXP","bpm0i01bXM","bpm0i01bYP","bpm0i01bYM",
        "bpm0i05XP","bpm0i05XM","bpm0i05YP","bpm0i05YM"
    };

    TString mPattern = "hd_pattern_number";
    TString vPattern = "pattern_number";

    const int patternOffset = 2777;

    // Use only stable beam states, excluding transition.
    const int highStart = 4584;
    const int highEnd   = 4638;

    const int lowStart  = 4655;
    const int lowEnd    = 4697;

    const double minCorrGood = 0.95;
    const double minCorrCheck = 0.80;

    TString pdfName = "all_molleradc_vqwk_channel_comparison.pdf";
    std::ofstream out("all_molleradc_vqwk_channel_comparison.csv");

    out << "Channel,Slope,Intercept,Corr,Npoints,Status\n";

    TCanvas *c = new TCanvas("c", "MOLLERADC vs VQWK comparison", 1600, 1000);
    bool firstPage = true;

    std::cout << std::left
              << std::setw(18) << "Channel"
              << std::setw(14) << "Slope"
              << std::setw(14) << "Intercept"
              << std::setw(12) << "Corr"
              << std::setw(10) << "N"
              << std::setw(14) << "Status"
              << std::endl;

    for (auto ch : channels)
    {
        TString status = "GOOD";

        if (!mEvt->GetBranch(ch)) {
            status = "MISSING_MOLLER";
            out << ch << ",NA,NA,NA,0," << status << "\n";
            continue;
        }

        if (!vEvt->GetBranch(ch)) {
            status = "MISSING_VQWK";
            out << ch << ",NA,NA,NA,0," << status << "\n";
            continue;
        }

        std::map<int, Stats> mMap;
        std::map<int, Stats> vMap;

        TString mExpr = Form("%s:%s", ch.Data(), mPattern.Data());
        TString vExpr = Form("%s:%s", ch.Data(), vPattern.Data());

        TString mCut = Form("((%s>=%d && %s<=%d) || (%s>=%d && %s<=%d))",
                            mPattern.Data(), highStart + patternOffset,
                            mPattern.Data(), highEnd + patternOffset,
                            mPattern.Data(), lowStart + patternOffset,
                            mPattern.Data(), lowEnd + patternOffset);

        TString vCut = Form("((%s>=%d && %s<=%d) || (%s>=%d && %s<=%d))",
                            vPattern.Data(), highStart,
                            vPattern.Data(), highEnd,
                            vPattern.Data(), lowStart,
                            vPattern.Data(), lowEnd);

        Long64_t nM = mEvt->Draw(mExpr, mCut, "goff");

        for (Long64_t i = 0; i < nM; i++) {
            double val = mEvt->GetV1()[i];
            int pat    = int(std::round(mEvt->GetV2()[i]));
            mMap[pat].add(val);
        }

        Long64_t nV = vEvt->Draw(vExpr, vCut, "goff");

        for (Long64_t i = 0; i < nV; i++) {
            double val = vEvt->GetV1()[i];
            int pat    = int(std::round(vEvt->GetV2()[i]));
            vMap[pat].add(val);
        }

        TGraph *gMpat = new TGraph();
        TGraph *gVpat = new TGraph();
        TGraph *gMV   = new TGraph();

        int iM = 0;
        int iV = 0;
        int iG = 0;

        for (auto const &vp : vMap) {
            int vKey = vp.first;
            int mKey = vKey + patternOffset;

            if (mMap.find(mKey) == mMap.end()) continue;

            double vMean = vp.second.mean();
            double mMean = mMap[mKey].mean();

            gVpat->SetPoint(iV, vKey, vMean);
            iV++;

            gMpat->SetPoint(iM, mKey, mMean);
            iM++;

            gMV->SetPoint(iG, vMean, mMean);
            iG++;
        }

        if (gMV->GetN() < 5) {
            status = "LOW_STATS";
            out << ch << ",NA,NA,NA," << gMV->GetN() << "," << status << "\n";

            delete gMpat;
            delete gVpat;
            delete gMV;
            continue;
        }

        TF1 *fit = new TF1(Form("fit_%s", ch.Data()), "pol1");

        gMV->Fit(fit, "Q");

        double intercept = fit->GetParameter(0);
        double slope     = fit->GetParameter(1);
        double corr      = gMV->GetCorrelationFactor();

        if (std::abs(corr) < minCorrCheck)
            status = "BAD";
        else if (std::abs(corr) < minCorrGood)
            status = "CHECK";
        else
            status = "GOOD";

        TGraph *gRes = new TGraph();

        for (int i = 0; i < gMV->GetN(); i++) {
            double x, y;
            gMV->GetPoint(i, x, y);
            double residual = y - fit->Eval(x);
            gRes->SetPoint(i, x, residual);
        }

        std::cout << std::left
                  << std::setw(18) << ch
                  << std::setw(14) << slope
                  << std::setw(14) << intercept
                  << std::setw(12) << corr
                  << std::setw(10) << gMV->GetN()
                  << std::setw(14) << status
                  << std::endl;

        out << ch << ","
            << slope << ","
            << intercept << ","
            << corr << ","
            << gMV->GetN() << ","
            << status << "\n";

        c->Clear();
        c->Divide(2, 2);

        c->cd(1);
        gVpat->SetTitle(Form("%s VQWK vs pattern;VQWK pattern;VQWK mean", ch.Data()));
        gVpat->SetMarkerStyle(20);
        gVpat->SetMarkerSize(0.7);
        gVpat->Draw("AP");

        c->cd(2);
        gMpat->SetTitle(Form("%s MOLLERADC vs pattern;MOLLER pattern;MOLLER mean", ch.Data()));
        gMpat->SetMarkerStyle(20);
        gMpat->SetMarkerSize(0.7);
        gMpat->Draw("AP");

        c->cd(3);
        gMV->SetTitle(Form("%s MOLLERADC vs VQWK | %s | corr=%.3f;VQWK mean;MOLLERADC mean",
                           ch.Data(), status.Data(), corr));
        gMV->SetMarkerStyle(20);
        gMV->SetMarkerSize(0.7);
        gMV->Draw("AP");
        fit->Draw("same");

        c->cd(4);
        gRes->SetTitle(Form("%s residual: M - fit(V);VQWK mean;Residual", ch.Data()));
        gRes->SetMarkerStyle(20);
        gRes->SetMarkerSize(0.7);
        gRes->Draw("AP");

        c->Update();

        if (firstPage) {
            c->Print(pdfName + "(");
            firstPage = false;
        } else {
            c->Print(pdfName);
        }

        delete gMpat;
        delete gVpat;
        delete gMV;
        delete gRes;
        delete fit;
    }

    c->Print(pdfName + ")");
    c->Close();

    out.close();

    std::cout << "\nSaved PDF: " << pdfName << std::endl;
    std::cout << "Saved CSV: all_molleradc_vqwk_channel_comparison.csv\n";
    std::cout << "Done.\n";
}
