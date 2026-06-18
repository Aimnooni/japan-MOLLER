#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TGraphErrors.h>
#include <TString.h>
#include <TStyle.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <cmath>

const char* FILE_PATTERN = "../RootFiles/isu_sample_%d.000.root";

struct ChannelInfo {
    std::string adcName;
    std::string channelName;
    int hwChan;
    double mean;
    double meanErr;
    double rms;
    double entries;
};

void baseline_by_adc(int run = 1307)
{
    gStyle->SetOptStat(1110);

    TString fileName = Form(FILE_PATTERN, run);
    TFile *f = TFile::Open(fileName, "READ");

    if (!f || f->IsZombie()) {
        std::cerr << "ERROR: Could not open " << fileName << std::endl;
        return;
    }

    TTree *evt = (TTree*) f->Get("evt");

    if (!evt) {
        std::cerr << "ERROR: Could not find evt tree." << std::endl;
        return;
    }

    TString pdfName = Form("baseline_by_adc_run%d.pdf", run);
    TString csvName = Form("baseline_by_adc_run%d.csv", run);

    std::ofstream csv(csvName.Data());
    csv << "Run,ADC,ChannelName,HWChannel,Mean,MeanError,RMS,Entries\n";

    // Example: TQ05_R1 ... TQ06_R6 are one ADC, channels 0-15.
    std::map<std::string, std::vector<ChannelInfo>> adcGroups;

    TObjArray *branches = evt->GetListOfBranches();

    for (int i = 0; i < branches->GetEntries(); i++) {

        TBranch *br = (TBranch*) branches->At(i);
        TString bname = br->GetName();

        TString lname = bname;
        lname.ToLower();

        // Keep only MOLLERADC TQ channels
        if (!lname.BeginsWith("tq")) continue;

        // Use hw_sum_rms branch only
        TString drawExpr = bname + ".hw_sum_rms";

        // Check if branch/member exists by trying Draw in graphics-off mode
        Long64_t n = evt->Draw(drawExpr, "", "goff");

        if (n <= 0) continue;

        // Parse tq number from name, e.g. tq05_r1 -> 5
        TString tqnumStr = lname(2, 2);
        int tqnum = tqnumStr.Atoi();

        if (tqnum <= 0) continue;

        // Group every two TQ groups as one ADC:
        // tq01/tq02 -> ADC01
        // tq03/tq04 -> ADC02
        // tq05/tq06 -> ADC03
        int adcIndex = (tqnum + 1) / 2;

        TString adcName = Form("ADC%02d", adcIndex);

        // Hardware channel number within one ADC
        // Odd tq group gives channels 0-7, even tq group gives channels 8-15
        int localBase = (tqnum % 2 == 1) ? 0 : 8;

        int localChan = -1;

        if      (lname.Contains("_r1"))  localChan = 0;
        else if (lname.Contains("_r2"))  localChan = 1;
        else if (lname.Contains("_r3"))  localChan = 2;
        else if (lname.Contains("_r4"))  localChan = 3;
        else if (lname.Contains("_r5l")) localChan = 4;
        else if (lname.Contains("_r5c")) localChan = 5;
        else if (lname.Contains("_r5r")) localChan = 6;
        else if (lname.Contains("_r6"))  localChan = 7;

        if (localChan < 0) continue;

        int hwChan = localBase + localChan;

        double minVal = evt->GetMinimum(drawExpr);
        double maxVal = evt->GetMaximum(drawExpr);

        if (minVal == maxVal) {
            minVal -= 1.0;
            maxVal += 1.0;
        }

        TH1D *h = new TH1D(Form("h_%s", lname.Data()),
                           Form("%s baseline;hw_sum_rms;Entries", bname.Data()),
                           200,
                           minVal,
                           maxVal);

        evt->Draw(drawExpr + Form(">>h_%s", lname.Data()), "", "goff");

        ChannelInfo info;
        info.adcName = adcName.Data();
        info.channelName = lname.Data();
        info.hwChan = hwChan;
        info.mean = h->GetMean();
        info.meanErr = h->GetMeanError();
        info.rms = h->GetRMS();
        info.entries = h->GetEntries();

        adcGroups[info.adcName].push_back(info);

        csv << run << ","
            << info.adcName << ","
            << info.channelName << ","
            << info.hwChan << ","
            << info.mean << ","
            << info.meanErr << ","
            << info.rms << ","
            << info.entries << "\n";

        delete h;
    }

    TCanvas *c = new TCanvas("c", "Baseline by ADC", 1200, 800);

    bool firstPage = true;

    for (auto &adc : adcGroups) {

        auto &channels = adc.second;

        TGraphErrors *g = new TGraphErrors();

        for (size_t i = 0; i < channels.size(); i++) {
            g->SetPoint(i, channels[i].hwChan, channels[i].mean);
            g->SetPointError(i, 0.0, channels[i].meanErr);
        }

        c->Clear();
        c->SetGrid(0,0);

        g->SetTitle(Form("Run %d Baseline: %s;ADC Channel Number;Mean Baseline hw_sum_rms",
                         run, adc.first.c_str()));

        g->SetMarkerStyle(20);
        g->SetMarkerSize(1.2);
        g->SetLineWidth(2);
        g->Draw("AP");

        g->GetXaxis()->SetLimits(-0.5, 15.5);
        g->GetYaxis()->SetTitleOffset(1.3);

        if (firstPage) {
            c->Print(pdfName + "(");
            firstPage = false;
        } else {
            c->Print(pdfName);
        }

        delete g;
    }

    c->Print(pdfName + ")");

    csv.close();

    std::cout << "Saved PDF: " << pdfName << std::endl;
    std::cout << "Saved CSV: " << csvName << std::endl;

    f->Close();
}
