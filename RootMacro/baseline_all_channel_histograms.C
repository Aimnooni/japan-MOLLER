#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TString.h>
#include <TStyle.h>
#include <TObjArray.h>
#include <TGraph.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>

const char* FILE_PATTERN = "../RootFiles/isu_sample_%d.000.root";

struct ChannelResult {
    TString name;
    double entries;
    double mean;
    double meanErr;
    double rms;
    double min;
    double max;
};

ChannelResult GetStats(TTree *evt, TString expr, TString ch)
{
    ChannelResult r;
    r.name = ch;
    r.entries = 0;
    r.mean = 0;
    r.meanErr = 0;
    r.rms = 0;
    r.min = 0;
    r.max = 0;

    Long64_t n = evt->Draw(expr, "", "goff");

    if (n <= 0) return r;

    double sum = 0.0;
    double sum2 = 0.0;
    double minVal =  1.0e30;
    double maxVal = -1.0e30;

    for (Long64_t i = 0; i < n; i++) {
        double y = evt->GetV1()[i];

        sum += y;
        sum2 += y * y;

        if (y < minVal) minVal = y;
        if (y > maxVal) maxVal = y;
    }

    r.entries = n;
    r.mean = sum / n;
    r.rms = std::sqrt(sum2 / n - r.mean * r.mean);
    r.meanErr = r.rms / std::sqrt((double)n);
    r.min = minVal;
    r.max = maxVal;

    return r;
}

void PrepareCanvas(TCanvas *c)
{
    c->Clear();
    c->Divide(4, 2, 0.01, 0.01);
}

void StylePad()
{
    gPad->SetLeftMargin(0.13);
    gPad->SetRightMargin(0.07);
    gPad->SetBottomMargin(0.13);
    gPad->SetTopMargin(0.10);
}

void baseline_all_channel_histograms(int run = 1307)
{
    gStyle->SetOptStat(1110);

    TString fileName = Form(FILE_PATTERN, run);

    TFile *f = TFile::Open(fileName, "READ");

    if (!f || f->IsZombie()) {
        std::cerr << "ERROR: Could not open " << fileName << std::endl;
        return;
    }

    TTree *evt = (TTree*)f->Get("evt");

    if (!evt) {
        std::cerr << "ERROR: Could not find evt tree." << std::endl;
        f->Close();
        return;
    }

    std::vector<TString> channels;

    TObjArray *branches = evt->GetListOfBranches();

    for (int i = 0; i < branches->GetEntries(); i++) {
        TBranch *br = (TBranch*)branches->At(i);
        TString bname = br->GetName();

        TString lname = bname;
        lname.ToLower();

        if (!lname.BeginsWith("tq")) continue;

        TString expr = bname + ".hw_sum";

        if (evt->Draw(expr, "", "goff") > 0)
            channels.push_back(bname);
    }

    TString baselinePDF = Form("baseline_hw_sum_histograms_run%d.pdf", run);
    TString rmsPDF      = Form("baseline_hw_sum_rms_histograms_run%d.pdf", run);
    TString eventPDF    = Form("baseline_hw_sum_vs_event_run%d.pdf", run);

    TString baselineCSV = Form("baseline_hw_sum_histograms_run%d.csv", run);
    TString rmsCSV      = Form("baseline_hw_sum_rms_histograms_run%d.csv", run);
    TString eventCSV    = Form("baseline_hw_sum_vs_event_run%d.csv", run);

    std::ofstream csvBase(baselineCSV.Data());
    std::ofstream csvRMS(rmsCSV.Data());
    std::ofstream csvEvent(eventCSV.Data());

    csvBase  << "Run,Channel,Entries,Mean,MeanError,RMS,Min,Max\n";
    csvRMS   << "Run,Channel,Entries,Mean,MeanError,RMS,Min,Max\n";
    csvEvent << "Run,Channel,Entries,Mean,MeanError,RMS,Min,Max\n";

    TCanvas *c = new TCanvas("c", "Baseline Checks", 1800, 1000);

    // ============================================================
    // 1. channel.hw_sum histograms
    // ============================================================

    bool first = true;
    int pad = 0;

    PrepareCanvas(c);

    for (size_t i = 0; i < channels.size(); i++) {

        TString ch = channels[i];
        TString expr = ch + ".hw_sum";

        ChannelResult r = GetStats(evt, expr, ch);

        csvBase << run << ","
                << ch << ","
                << r.entries << ","
                << r.mean << ","
                << r.meanErr << ","
                << r.rms << ","
                << r.min << ","
                << r.max << "\n";

        pad++;
        c->cd(pad);
        StylePad();

        evt->Draw(expr);

        TH1 *h = (TH1*)gPad->GetPrimitive("htemp");

        if (h) {
            h->SetTitle(Form("%s.hw_sum", ch.Data()));
            h->GetXaxis()->SetTitle("hw_sum");
            h->GetYaxis()->SetTitle("Entries");
            h->SetLineColor(kBlack);
            h->SetLineWidth(1);
            h->SetFillStyle(0);
        }

        if (pad == 8) {
            if (first) {
                c->Print(baselinePDF + "(");
                first = false;
            } else {
                c->Print(baselinePDF);
            }

            PrepareCanvas(c);
            pad = 0;
        }
    }

    if (pad > 0) {
        if (first) {
            c->Print(baselinePDF + "(");
            first = false;
        } else {
            c->Print(baselinePDF);
        }
    }

    if (!first) c->Print(baselinePDF + ")");

    // ============================================================
    // 2. channel.hw_sum_rms histograms
    // ============================================================

    first = true;
    pad = 0;

    PrepareCanvas(c);

    for (size_t i = 0; i < channels.size(); i++) {

        TString ch = channels[i];
        TString expr = ch + ".hw_sum_rms";

        if (evt->Draw(expr, "", "goff") <= 0) continue;

        ChannelResult r = GetStats(evt, expr, ch);

        csvRMS << run << ","
               << ch << ","
               << r.entries << ","
               << r.mean << ","
               << r.meanErr << ","
               << r.rms << ","
               << r.min << ","
               << r.max << "\n";

        pad++;
        c->cd(pad);
        StylePad();

        evt->Draw(expr);

        TH1 *h = (TH1*)gPad->GetPrimitive("htemp");

        if (h) {
            h->SetTitle(Form("%s.hw_sum_rms", ch.Data()));
            h->GetXaxis()->SetTitle("hw_sum_rms");
            h->GetYaxis()->SetTitle("Entries");
            h->SetLineColor(kBlack);
            h->SetLineWidth(1);
            h->SetFillStyle(0);
        }

        if (pad == 8) {
            if (first) {
                c->Print(rmsPDF + "(");
                first = false;
            } else {
                c->Print(rmsPDF);
            }

            PrepareCanvas(c);
            pad = 0;
        }
    }

    if (pad > 0) {
        if (first) {
            c->Print(rmsPDF + "(");
            first = false;
        } else {
            c->Print(rmsPDF);
        }
    }

    if (!first) c->Print(rmsPDF + ")");

    // ============================================================
    // 3. evt->Draw("channel.hw_sum:Entry$")
    // ============================================================

    first = true;
    pad = 0;

    PrepareCanvas(c);

    for (size_t i = 0; i < channels.size(); i++) {

        TString ch = channels[i];
        TString expr = ch + ".hw_sum";

        ChannelResult r = GetStats(evt, expr, ch);

        csvEvent << run << ","
                 << ch << ","
                 << r.entries << ","
                 << r.mean << ","
                 << r.meanErr << ","
                 << r.rms << ","
                 << r.min << ","
                 << r.max << "\n";

        pad++;
        c->cd(pad);
        StylePad();

        evt->Draw(Form("%s:Entry$", expr.Data()));

        TGraph *g = (TGraph*)gPad->GetPrimitive("Graph");

        if (g) {
            g->SetTitle(Form("%s.hw_sum vs Entry$;Entry$;hw_sum", ch.Data()));
            g->SetMarkerStyle(7);
            g->SetMarkerColor(kBlack);
            g->SetLineColor(kBlack);
        }

        if (pad == 8) {
            if (first) {
                c->Print(eventPDF + "(");
                first = false;
            } else {
                c->Print(eventPDF);
            }

            PrepareCanvas(c);
            pad = 0;
        }
    }

    if (pad > 0) {
        if (first) {
            c->Print(eventPDF + "(");
            first = false;
        } else {
            c->Print(eventPDF);
        }
    }

    if (!first) c->Print(eventPDF + ")");

    csvBase.close();
    csvRMS.close();
    csvEvent.close();

    std::cout << "\nProcessed channels: " << channels.size() << std::endl;

    std::cout << "Saved PDF: " << baselinePDF << std::endl;
    std::cout << "Saved CSV: " << baselineCSV << std::endl;

    std::cout << "Saved PDF: " << rmsPDF << std::endl;
    std::cout << "Saved CSV: " << rmsCSV << std::endl;

    std::cout << "Saved PDF: " << eventPDF << std::endl;
    std::cout << "Saved CSV: " << eventCSV << std::endl;

    f->Close();
}
