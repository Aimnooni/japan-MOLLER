#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TString.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <iomanip>
#include <cmath>

struct AvgData {
    double sum = 0.0;
    int count = 0;

    void add(double x) {
        sum += x;
        count++;
    }

    double mean() const {
        return sum / count;
    }
};

void compare_moller_vqwk_bpm()
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

    std::map<TString, double> vPed = {
        {"bpm2i00XP", -894.067}, {"bpm2i00XM", -826.56},
        {"bpm2i00YP", -400.037}, {"bpm2i00YM", -439.255},

        {"bpm2i00aXP", -4748.41}, {"bpm2i00aXM", -5067.03},
        {"bpm2i00aYP", -3537.54}, {"bpm2i00aYM", -3321.07},

        {"bpm2i01aXP", -3530.8}, {"bpm2i01aXM", -3203.5},
        {"bpm2i01aYP", -2295.62}, {"bpm2i01aYM", -2042.09},

        {"bpm1i02XP", -4992.25}, {"bpm1i02XM", -4417.44},
        {"bpm1i02YP", -3804.62}, {"bpm1i02YM", -3298.87},

        {"bpm1i04XP", -3452.05}, {"bpm1i04XM", -3289.27},
        {"bpm1i04YP", -2813.34}, {"bpm1i04YM", -3551.81},

        {"bpm1i05XP", -5120.36}, {"bpm1i05XM", -4660.85},
        {"bpm1i05YP", -4138.57}, {"bpm1i05YM", -4034.45},

        {"bpm1i07XP", -4810}, {"bpm1i07XM", -4861},
        {"bpm1i07YP", -3567.44}, {"bpm1i07YM", -3238.15},

        {"bpm0i07XP", -215.624}, {"bpm0i07XM", -1224.9},
        {"bpm0i07YP", -826.275}, {"bpm0i07YM", 320.886},

        {"bpm0l01XP", 956.9}, {"bpm0l01XM", 695.673},
        {"bpm0l01YP", 1396.42}, {"bpm0l01YM", 2192.06},

        {"bpm0l02XP", -1600.1}, {"bpm0l02XM", 58.7412},
        {"bpm0l02YP", -2422.98}, {"bpm0l02YM", -208.547},

        {"bpm0l03XP", -349.076}, {"bpm0l03XM", -690.771},
        {"bpm0l03YP", 156.417}, {"bpm0l03YM", -740.477},

        {"bpm0l04XP", 10847.8}, {"bpm0l04XM", 10838.5},
        {"bpm0l04YP", 10582.2}, {"bpm0l04YM", 10657.9},

        {"bpm0l05XP", 9365.08}, {"bpm0l05XM", 9908.91},
        {"bpm0l05YP", 12726.7}, {"bpm0l05YM", 11577.4},

        {"bpm0l06XP", 5380.28}, {"bpm0l06XM", 5589.29},
        {"bpm0l06YP", 5355.09}, {"bpm0l06YM", 5293.64},

        {"bpm0r05XP", 3353.69}, {"bpm0r05XM", 3785.37},
        {"bpm0r05YP", 3184.72}, {"bpm0r05YM", 3321.5},

        {"bpm0l07XP", 14674.1}, {"bpm0l07XM", 14947.1},
        {"bpm0l07YP", 9902.09}, {"bpm0l07YM", 9866.24},

        {"bpm0l08XP", 7715.82}, {"bpm0l08XM", 7887.57},
        {"bpm0l08YP", 3878.59}, {"bpm0l08YM", 4278.31},

        {"bpm0l09XP", 6925.28}, {"bpm0l09XM", 7430.81},
        {"bpm0l09YP", 6899.51}, {"bpm0l09YM", 6731.63},

        {"bpm0l10XP", 5979.61}, {"bpm0l10XM", 5924.27},
        {"bpm0l10YP", 6878.67}, {"bpm0l10YM", 6794.42},

        {"bpm0r03XP", 7549.79}, {"bpm0r03XM", 7831.56},
        {"bpm0r03YP", 6828.55}, {"bpm0r03YM", 6872.59},

        {"bpm0i01XP", -7863.07}, {"bpm0i01XM", -4674.84},
        {"bpm0i01YP", -5402.72}, {"bpm0i01YM", 5218.1},

        {"bpm0i01bXP", -5909.88}, {"bpm0i01bXM", -5395.24},
        {"bpm0i01bYP", -2291.45}, {"bpm0i01bYM", -1876.6},

        {"bpm0i05XP", -6829.98}, {"bpm0i05XM", -4522.86},
        {"bpm0i05YP", -4185.67}, {"bpm0i05YM", -4586}
    };

    TString mPatternBranch = "hd_pattern_number";
    TString vPatternBranch = "pattern_number";

    double mPatternMin = 7360;
    double mPatternMax = 7475;

    double vPatternMin = 4583;
    double vPatternMax = 4698;

    int patternOffset = int(mPatternMin - vPatternMin);

    double ratioMin = 1.5;
    double ratioMax = 3.5;

    double goodCorr  = 0.95;
    double checkCorr = 0.90;

    TString pdfName = "all_bpm_pattern_averaged_cleaned_fits.pdf";

    TCanvas *cPage = new TCanvas("cPage", "Cleaned Pattern-Averaged BPM Fits", 1600, 1000);
    cPage->Divide(2, 2);

    int plotCounter = 0;
    bool firstPage = true;
    bool madePlot = false;

    std::ofstream out("all_bpm_pattern_averaged_cleaned_results.csv");
    out << "Channel,Vped,Slope,Intercept,Mped,Corr,Ncandidate,Nused,Nskipped,Status\n";

    std::cout << std::left
              << std::setw(18) << "Channel"
              << std::setw(12) << "Slope"
              << std::setw(14) << "Intercept"
              << std::setw(12) << "Mped"
              << std::setw(10) << "Corr"
              << std::setw(8)  << "Ncand"
              << std::setw(8)  << "Nused"
              << std::setw(8)  << "Nskip"
              << std::setw(14) << "Status"
              << std::endl;

    for (auto ch : channels)
    {
        TString status = "GOOD";

        if (!mEvt->GetBranch(ch) || !vEvt->GetBranch(ch)) {
            status = "MISSING";
            out << ch << ",NA,NA,NA,NA,NA,0,0,0," << status << "\n";
            std::cout << std::left << std::setw(18) << ch << status << std::endl;
            continue;
        }

        if (!mEvt->GetBranch(mPatternBranch) || !vEvt->GetBranch(vPatternBranch)) {
            std::cout << "ERROR: Missing pattern branch\n";
            return;
        }

        if (vPed.find(ch) == vPed.end()) {
            status = "NO_VPED";
            out << ch << ",NA,NA,NA,NA,NA,0,0,0," << status << "\n";
            std::cout << std::left << std::setw(18) << ch << status << std::endl;
            continue;
        }

        TString mExpr = Form("%s:%s", ch.Data(), mPatternBranch.Data());
        TString vExpr = Form("%s:%s", ch.Data(), vPatternBranch.Data());

        TString mCut = Form("%s>%f && %s<%f",
                            mPatternBranch.Data(), mPatternMin,
                            mPatternBranch.Data(), mPatternMax);

        TString vCut = Form("%s>%f && %s<%f",
                            vPatternBranch.Data(), vPatternMin,
                            vPatternBranch.Data(), vPatternMax);

        std::map<int, AvgData> mMap;
        std::map<int, AvgData> vMap;

        Long64_t nM = mEvt->Draw(mExpr, mCut, "goff");

        for (Long64_t i = 0; i < nM; i++) {
            double value = mEvt->GetV1()[i];
            int pattern  = int(std::round(mEvt->GetV2()[i]));
            mMap[pattern].add(value);
        }

        Long64_t nV = vEvt->Draw(vExpr, vCut, "goff");

        for (Long64_t i = 0; i < nV; i++) {
            double value = vEvt->GetV1()[i];
            int pattern  = int(std::round(vEvt->GetV2()[i]));
            vMap[pattern].add(value);
        }

        TGraph *gr = new TGraph();

        int nCandidate = 0;
        int nUsed = 0;
        int nSkipped = 0;

        for (auto const &vp : vMap) {
            int vKey = vp.first;
            int mKey = vKey + patternOffset;

            if (mMap.find(mKey) == mMap.end()) continue;

            double vMean = vp.second.mean();
            double mMean = mMap[mKey].mean();

            if (std::abs(vMean) < 1e-9) continue;

            nCandidate++;

            double ratio = mMean / vMean;

            if (ratio < ratioMin || ratio > ratioMax) {
                nSkipped++;
                continue;
            }

            gr->SetPoint(nUsed, vMean, mMean);
            nUsed++;
        }

        if (nUsed < 5) {
            status = "LOW_STATS";

            out << ch << "," << vPed[ch]
                << ",NA,NA,NA,NA,"
                << nCandidate << ","
                << nUsed << ","
                << nSkipped << ","
                << status << "\n";

            std::cout << std::left
                      << std::setw(18) << ch
                      << std::setw(12) << "NA"
                      << std::setw(14) << "NA"
                      << std::setw(12) << "NA"
                      << std::setw(10) << "NA"
                      << std::setw(8)  << nCandidate
                      << std::setw(8)  << nUsed
                      << std::setw(8)  << nSkipped
                      << std::setw(14) << status
                      << std::endl;

            delete gr;
            continue;
        }

        TF1 *fit = new TF1(Form("fit_%s", ch.Data()), "pol1");
        gr->Fit(fit, "Q");

        double intercept = fit->GetParameter(0);
        double slope     = fit->GetParameter(1);
        double corr      = gr->GetCorrelationFactor();
        double mped      = slope * vPed[ch] + intercept;

        if (std::abs(corr) < checkCorr)
            status = "BAD";
        else if (std::abs(corr) < goodCorr)
            status = "CHECK";
        else
            status = "GOOD";

        std::cout << std::left
                  << std::setw(18) << ch
                  << std::setw(12) << slope
                  << std::setw(14) << intercept
                  << std::setw(12) << mped
                  << std::setw(10) << corr
                  << std::setw(8)  << nCandidate
                  << std::setw(8)  << nUsed
                  << std::setw(8)  << nSkipped
                  << std::setw(14) << status
                  << std::endl;

        out << ch << ","
            << vPed[ch] << ","
            << slope << ","
            << intercept << ","
            << mped << ","
            << corr << ","
            << nCandidate << ","
            << nUsed << ","
            << nSkipped << ","
            << status << "\n";

        int pad = (plotCounter % 4) + 1;
        cPage->cd(pad);

        gr->SetTitle(Form("%s | %s | Corr = %.3f | skipped = %d;VQWK pattern avg;MOLLERADC pattern avg",
                          ch.Data(), status.Data(), corr, nSkipped));

        gr->SetMarkerStyle(20);
        gr->SetMarkerSize(0.6);
        gr->Draw("AP");
        fit->Draw("same");

        plotCounter++;
        madePlot = true;

        if (plotCounter % 4 == 0) {
            cPage->Update();

            if (firstPage) {
                cPage->Print(pdfName + "(");
                firstPage = false;
            } else {
                cPage->Print(pdfName);
            }

            cPage->Clear();
            cPage->Divide(2, 2);
        }

        delete gr;
        delete fit;
    }

    if (madePlot) {
        if (plotCounter % 4 != 0) {
            cPage->Update();

            if (firstPage) {
                cPage->Print(pdfName + "(");
                firstPage = false;
            } else {
                cPage->Print(pdfName);
            }
        }

        cPage->Print(pdfName + ")");
        cPage->Close();

        std::cout << "\nSaved plots to: " << pdfName << std::endl;
    }

    out.close();

    std::cout << "Saved results to: all_bpm_pattern_averaged_cleaned_results.csv\n";
    std::cout << "Done.\n";
}
