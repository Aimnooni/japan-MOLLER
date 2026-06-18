#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>

void check_bpm_saturation()
{
    TFile *f = TFile::Open("../RootFiles/isu_sample_1269.000.root");

    if (!f || f->IsZombie()) {
        std::cout << "ERROR: Could not open file" << std::endl;
        return;
    }

    TTree *evt = (TTree*) f->Get("evt");

    if (!evt) {
        std::cout << "ERROR: Could not find evt tree" << std::endl;
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

    std::cout << "\nBPM Saturation Check\n";
    std::cout << "Scale: channel * 4.0 / pow(2,17)\n";
    std::cout << "No cuts applied\n";
    std::cout << "Allowed range: [-4, 4]\n\n";

    std::cout << std::left
              << std::setw(18) << "Channel"
              << std::setw(12) << "Entries"
              << std::setw(14) << "MinScaled"
              << std::setw(14) << "MaxScaled"
              << std::setw(12) << "SatCount"
              << std::setw(12) << "Status"
              << std::endl;

    for (auto ch : channels)
    {
        TString expr = Form("(%s)*4.0/pow(2,17)", ch.Data());

        TString satCut = Form(
            "((%s)*4.0/pow(2,17) < -4 || (%s)*4.0/pow(2,17) > 4)",
            ch.Data(), ch.Data()
        );

        Long64_t nEntries = evt->GetEntries();
        Long64_t satCount = evt->GetEntries(satCut);

        Long64_t nDraw = evt->Draw(expr, "", "goff");

        if (nDraw <= 0) {
            std::cout << std::left
                      << std::setw(18) << ch
                      << std::setw(12) << 0
                      << std::setw(14) << "NA"
                      << std::setw(14) << "NA"
                      << std::setw(12) << "NA"
                      << std::setw(12) << "NO DATA"
                      << std::endl;
            continue;
        }

        double minVal = evt->GetV1()[0];
        double maxVal = evt->GetV1()[0];

        for (Long64_t i = 0; i < nDraw; i++) {
            double val = evt->GetV1()[i];
            if (val < minVal) minVal = val;
            if (val > maxVal) maxVal = val;
        }

        TString status = "OK";
        if (satCount > 0) status = "SATURATED";

        std::cout << std::left
                  << std::setw(18) << ch
                  << std::setw(12) << nEntries
                  << std::setw(14) << minVal
                  << std::setw(14) << maxVal
                  << std::setw(12) << satCount
                  << std::setw(12) << status
                  << std::endl;
    }

    std::cout << "\nDone.\n";
}
