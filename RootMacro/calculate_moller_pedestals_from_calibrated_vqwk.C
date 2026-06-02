#include <TFile.h>
#include <TTree.h>
#include <TString.h>

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

void calculate_moller_pedestals_from_calibrated_vqwk()
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

    const int highStart = 4584;
    const int highEnd   = 4638;

    const int lowStart  = 4655;
    const int lowEnd    = 4697;

    const double minVSeparation = 100.0;
    const double gainMin = 1.0;
    const double gainMax = 5.0;

    std::ofstream out("moller_pedestal_from_calibrated_vqwk.csv");

    out << "Channel,Von_cal,Voff_cal,Mon,Moff,Gain,Mped,Status\n";

    std::cout << std::left
              << std::setw(18) << "Channel"
              << std::setw(14) << "Gain"
              << std::setw(14) << "Mped"
              << std::setw(16) << "Status"
              << std::endl;

    for (auto ch : channels)
    {
        TString status = "GOOD";

        if (!mEvt->GetBranch(ch)) {
            status = "MISSING_MOLLER";
            out << ch << ",NA,NA,NA,NA,NA,NA," << status << "\n";
            continue;
        }

        if (!vEvt->GetBranch(ch)) {
            status = "MISSING_VQWK";
            out << ch << ",NA,NA,NA,NA,NA,NA," << status << "\n";
            continue;
        }

        if (!mEvt->GetBranch(mPattern) || !vEvt->GetBranch(vPattern)) {
            std::cout << "ERROR: Missing pattern branch\n";
            return;
        }

        std::map<int, Stats> mMap;
        std::map<int, Stats> vMap;

        TString mExpr = Form("%s:%s", ch.Data(), mPattern.Data());
        TString vExpr = Form("%s:%s", ch.Data(), vPattern.Data());

        TString mCut = Form("%s>=%d && %s<=%d",
                            mPattern.Data(), highStart + patternOffset,
                            mPattern.Data(), lowEnd + patternOffset);

        TString vCut = Form("%s>=%d && %s<=%d",
                            vPattern.Data(), highStart,
                            vPattern.Data(), lowEnd);

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

        Stats highV, highM, lowV, lowM;

        for (auto const &vp : vMap) {
            int vKey = vp.first;
            int mKey = vKey + patternOffset;

            if (mMap.find(mKey) == mMap.end()) continue;

            double vMean = vp.second.mean();
            double mMean = mMap[mKey].mean();

            if (vKey >= highStart && vKey <= highEnd) {
                highV.add(vMean);
                highM.add(mMean);
            }
            else if (vKey >= lowStart && vKey <= lowEnd) {
                lowV.add(vMean);
                lowM.add(mMean);
            }
        }

        double Von  = highV.mean();
        double Voff = lowV.mean();
        double Mon  = highM.mean();
        double Moff = lowM.mean();

        double gain = 0.0;
        double Mped = 0.0;

        if (highV.n < 5 || lowV.n < 5) {
            status = "LOW_STATS";
        }
        else if (std::abs(Von - Voff) < minVSeparation) {
            status = "LOW_V_SEPARATION";
        }
        else {
            gain = (Mon - Moff) / (Von - Voff);

            // Calibrated VQWK values already have pedestal removed.
            // Therefore Vped = 0.
            Mped = Mon - gain * Von;

            if (gain < gainMin || gain > gainMax)
                status = "GAIN_OUTLIER";
            else
                status = "GOOD";
        }

        std::cout << std::left
                  << std::setw(18) << ch
                  << std::setw(14) << gain
                  << std::setw(14) << Mped
                  << std::setw(16) << status
                  << std::endl;

        out << ch << ","
            << Von << ","
            << Voff << ","
            << Mon << ","
            << Moff << ","
            << gain << ","
            << Mped << ","
            << status << "\n";
    }

    out.close();

    std::cout << "\nSaved: moller_pedestal_from_calibrated_vqwk.csv\n";
    std::cout << "Done.\n";
}
