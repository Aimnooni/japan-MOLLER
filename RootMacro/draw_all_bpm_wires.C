#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TString.h>
#include <TROOT.h>
#include <iostream>
#include <vector>
#include <cmath>

void draw_all_bpm_wires()
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

    TString pdfName = "all_bpm_wire_plots.pdf";

    for (size_t i = 0; i < channels.size(); i += 4)
    {
        TCanvas *c = new TCanvas(
            Form("c_bpm_%zu", i / 4),
            "BPM wire plots",
            1600,
            1000
        );

        c->Divide(2, 2);

        for (int j = 0; j < 4; j++)
        {
            if (i + j >= channels.size()) break;

            c->cd(j + 1);

            TString ch = channels[i + j];

            TString drawExpr = Form(
                "%s*4.0/pow(2,17):bcm0l02",
                ch.Data()
            );

            evt->Draw(drawExpr);
        }

        c->Update();

        if (i == 0)
            c->Print(pdfName + "(");
        else if (i + 4 >= channels.size())
            c->Print(pdfName + ")");
        else
            c->Print(pdfName);
    }

    std::cout << "Saved: " << pdfName << std::endl;
}
