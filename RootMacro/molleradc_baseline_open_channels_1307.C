#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TString.h>
#include <TStyle.h>

#include <iostream>
#include <fstream>
#include <vector>

void molleradc_baseline_open_channels_1307()
{
    TString mollerFile = "../RootFiles/isu_sample_1307.000.root";

    TFile *f = TFile::Open(mollerFile);

    if (!f || f->IsZombie()) {
        std::cout << "Cannot open file\n";
        return;
    }

    TTree *evt = (TTree*)f->Get("evt");

    if (!evt) {
        std::cout << "Cannot find evt tree\n";
        return;
    }

    gStyle->SetOptStat(1110);

    std::vector<TString> channels = {
        "tq05_r1","tq05_r2","tq05_r3","tq05_r4",
        "tq05_r5l","tq05_r5c","tq05_r5r","tq05_r6",

        "tq06_r1","tq06_r2","tq06_r3","tq06_r4",
        "tq06_r5l","tq06_r5c","tq06_r5r","tq06_r6"
    };

    TString pdfName =
        "molleradc_baseline_tq05_tq06.pdf";

    TString csvName =
        "molleradc_baseline_tq05_tq06.csv";

    std::ofstream out(csvName.Data());

    out << "Channel,"
        << "BaselineMean,"
        << "BaselineRMS,"
        << "RMSMean,"
        << "RMSRMS,"
        << "Entries\n";

    TCanvas *c =
        new TCanvas("c",
                    "Baseline Study",
                    1600,
                    900);

    bool firstPage = true;

    for (auto ch : channels)
    {
        if (!evt->GetBranch(ch))
        {
            std::cout << "Missing branch: "
                      << ch << std::endl;
            continue;
        }

        TString baselineExpr =
            ch + ".hw_sum";

        TString rmsExpr =
            ch + ".hw_sum_rms";

        c->Clear();
        c->Divide(2,1);

        //--------------------------------------------------
        // Baseline histogram
        //--------------------------------------------------

        c->cd(1);

        evt->Draw(baselineExpr);

        TH1 *hBase =
            (TH1*)gPad->GetPrimitive("htemp");

        double baseMean = 0;
        double baseRMS  = 0;
        double entries  = 0;

        if (hBase)
        {
            hBase->SetTitle(
                Form("%s Baseline",
                     ch.Data()));

            baseMean = hBase->GetMean();
            baseRMS  = hBase->GetRMS();
            entries  = hBase->GetEntries();
        }

        //--------------------------------------------------
        // RMS histogram
        //--------------------------------------------------

        c->cd(2);

        evt->Draw(rmsExpr);

        TH1 *hRMS =
            (TH1*)gPad->GetPrimitive("htemp");

        double rmsMean = 0;
        double rmsRMS  = 0;

        if (hRMS)
        {
            hRMS->SetTitle(
                Form("%s RMS",
                     ch.Data()));

            rmsMean = hRMS->GetMean();
            rmsRMS  = hRMS->GetRMS();
        }

        out << ch << ","
            << baseMean << ","
            << baseRMS << ","
            << rmsMean << ","
            << rmsRMS << ","
            << entries << "\n";

        if (firstPage)
        {
            c->Print(pdfName + "(");
            firstPage = false;
        }
        else
        {
            c->Print(pdfName);
        }

        std::cout
            << "Processed "
            << ch
            << std::endl;
    }

    c->Print(pdfName + ")");

    out.close();

    std::cout << "\nSaved PDF: "
              << pdfName
              << std::endl;

    std::cout << "Saved CSV: "
              << csvName
              << std::endl;
}
