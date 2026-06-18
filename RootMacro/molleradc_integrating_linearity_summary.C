#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TLine.h>
#include <TStyle.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

struct RunInfo {
    int run;
    double vpp;
};

struct Result {
    double mean = 0.0;
    double rms = 0.0;
    double rmsAmp = 0.0;
    double rmsAmpErr = 0.0;
    double entries = 0.0;
};

Result AnalyzeRun(TTree *evt, TString ch)
{
    Result r;

    TString expr = ch + ".hw_sum";

    Long64_t n = evt->Draw(expr, "", "goff");

    if (n <= 0) return r;

    r.entries = n;

    double sum = 0.0;
    double sum2 = 0.0;

    for (Long64_t i = 0; i < n; i++) {
        double y = evt->GetV1()[i];
        sum += y;
        sum2 += y * y;
    }

    r.mean = sum / n;
    r.rms = std::sqrt(sum2 / n - r.mean * r.mean);

    r.rmsAmp = std::sqrt(2.0) * r.rms;
    r.rmsAmpErr = r.rmsAmp / std::sqrt(2.0 * n);

    return r;
}

void StyleGraph(TGraphErrors *g, int color, int marker)
{
    g->SetMarkerStyle(marker);
    g->SetMarkerSize(1.2);
    g->SetMarkerColor(color);
    g->SetLineColor(color);
    g->SetLineWidth(2);
}

void molleradc_integrating_linearity_summary()
{
    gStyle->SetOptStat(0);

    TString ch1 = "tq05_r1";
    TString ch2 = "tq05_r2";

    std::vector<RunInfo> runs = {
        {1379, 0.5},
        {1380, 1.0},
        {1381, 2.0},
        {1382, 3.0},
        {1385, 4.0},
        {1386, 5.0},
        {1387, 6.0},
        {1388, 7.0},
        {1389, 7.5},
        {1390, 8.0}
    };

    TString pdfName = "molleradc_integrating_linearity_summary.pdf";
    TString csvName = "molleradc_integrating_linearity_summary.csv";

    std::ofstream out(csvName.Data());

    out << "Run,Vpp,Channel,Entries,Mean_ADC,RMS_ADC,"
        << "RMSAmplitude_ADC,RMSAmplitudeError_ADC,"
        << "RMSAmplitudePerVpp_ADC\n";

    std::vector<double> vppList;
    std::vector<double> amp1List, amp2List;
    std::vector<double> err1List, err2List;

    for (auto r : runs) {
        TString fileName = Form("../RootFiles/isu_sample_%d.000.root", r.run);

        TFile *f = TFile::Open(fileName);

        if (!f || f->IsZombie()) {
            std::cout << "ERROR: Could not open " << fileName << std::endl;
            continue;
        }

        TTree *evt = (TTree*)f->Get("evt");

        if (!evt) {
            std::cout << "ERROR: Could not find evt tree in run "
                      << r.run << std::endl;
            f->Close();
            continue;
        }

        Result a = AnalyzeRun(evt, ch1);
        Result b = AnalyzeRun(evt, ch2);

        vppList.push_back(r.vpp);

        amp1List.push_back(a.rmsAmp);
        amp2List.push_back(b.rmsAmp);

        err1List.push_back(a.rmsAmpErr);
        err2List.push_back(b.rmsAmpErr);

        out << r.run << "," << r.vpp << "," << ch1 << ","
            << a.entries << "," << a.mean << "," << a.rms << ","
            << a.rmsAmp << "," << a.rmsAmpErr << ","
            << a.rmsAmp / r.vpp << "\n";

        out << r.run << "," << r.vpp << "," << ch2 << ","
            << b.entries << "," << b.mean << "," << b.rms << ","
            << b.rmsAmp << "," << b.rmsAmpErr << ","
            << b.rmsAmp / r.vpp << "\n";

        f->Close();
    }

    int n = vppList.size();

    TGraphErrors *g1 =
        new TGraphErrors(n, vppList.data(), amp1List.data(),
                         nullptr, err1List.data());

    TGraphErrors *g2 =
        new TGraphErrors(n, vppList.data(), amp2List.data(),
                         nullptr, err2List.data());

    TF1 *fit1 = new TF1("fit1", "[0]*x+[1]", 0.0, 8.2);
    TF1 *fit2 = new TF1("fit2", "[0]*x+[1]", 0.0, 8.2);

    fit1->SetLineColor(kBlue + 1);
    fit1->SetLineWidth(2);

    fit2->SetLineColor(kRed + 1);
    fit2->SetLineWidth(2);

    TCanvas *c =
        new TCanvas("c",
                    "Integrating Mode Linearity Summary",
                    1600,
                    1000);

    // ------------------------------------------------------------
    // Page 1: amplitude response
    // ------------------------------------------------------------

    c->Clear();
    c->SetGridx();
    c->SetGridy();

    StyleGraph(g1, kBlue + 1, 20);
    StyleGraph(g2, kRed + 1, 21);

    g1->SetTitle("Integrating Mode Amplitude Response;Input Amplitude (Vpp);Measured Amplitude #sqrt{2} RMS (ADC Counts)");
    g1->Draw("AP");
    g1->GetXaxis()->SetLimits(0.0, 8.2);

    g1->Fit(fit1, "Q");
    g2->Fit(fit2, "Q");

    g2->Draw("P SAME");
    fit1->Draw("same");
    fit2->Draw("same");

    TLegend *leg = new TLegend(0.58, 0.70, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(g1, "tq05_r1", "lp");
    leg->AddEntry(g2, "tq05_r2", "lp");
    leg->AddEntry(fit1, "Linear fit tq05_r1", "l");
    leg->AddEntry(fit2, "Linear fit tq05_r2", "l");
    leg->Draw();

    c->Print(pdfName + "(");

    // ------------------------------------------------------------
    // Page 2: residuals
    // ------------------------------------------------------------

    TGraphErrors *res1 = new TGraphErrors();
    TGraphErrors *res2 = new TGraphErrors();

    for (int i = 0; i < n; i++) {
        double x = vppList[i];

        double p1 = fit1->Eval(x);
        double p2 = fit2->Eval(x);

        double r1 = 100.0 * (amp1List[i] - p1) / p1;
        double r2 = 100.0 * (amp2List[i] - p2) / p2;

        res1->SetPoint(i, x, r1);
        res2->SetPoint(i, x, r2);
    }

    c->Clear();
    c->SetGridx();
    c->SetGridy();

    StyleGraph(res1, kBlue + 1, 20);
    StyleGraph(res2, kRed + 1, 21);

    res1->SetTitle("Integrating Mode Linearity Residuals;Input Amplitude (Vpp);Residual (%)");
    res1->Draw("APL");
    res1->GetXaxis()->SetLimits(0.0, 8.2);

    res2->Draw("PL SAME");

    TLine *zero = new TLine(0.0, 0.0, 8.2, 0.0);
    zero->SetLineStyle(2);
    zero->SetLineColor(kBlack);
    zero->Draw("same");

    TLegend *leg2 = new TLegend(0.62, 0.74, 0.88, 0.88);
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);
    leg2->AddEntry(res1, "tq05_r1", "lp");
    leg2->AddEntry(res2, "tq05_r2", "lp");
    leg2->AddEntry(zero, "Zero residual", "l");
    leg2->Draw();

    c->Print(pdfName);

    // ------------------------------------------------------------
    // Page 3: summary
    // ------------------------------------------------------------

    double slope1 = fit1->GetParameter(0);
    double int1   = fit1->GetParameter(1);
    double slope2 = fit2->GetParameter(0);
    double int2   = fit2->GetParameter(1);

    c->Clear();
    c->SetGridx(0);
    c->SetGridy(0);

    TPaveText *summary = new TPaveText(0.12, 0.20, 0.88, 0.82, "NDC");
    summary->SetFillColor(0);
    summary->SetBorderSize(1);
    summary->SetTextAlign(12);
    summary->SetTextSize(0.035);

    summary->AddText("Integrating Mode Amplitude Response Linearity");
    summary->AddText("----------------------------------------");
    summary->AddText("Analysis variable: hw_sum");
    summary->AddText("Amplitude estimator: sqrt(2) * RMS(hw_sum)");
    summary->AddText(Form("%s slope     : %.6f ADC counts / Vpp", ch1.Data(), slope1));
    summary->AddText(Form("%s intercept : %.6f ADC counts", ch1.Data(), int1));
    summary->AddText(Form("%s slope     : %.6f ADC counts / Vpp", ch2.Data(), slope2));
    summary->AddText(Form("%s intercept : %.6f ADC counts", ch2.Data(), int2));
    summary->AddText(Form("Slope ratio (%s/%s): %.8f",
                          ch2.Data(), ch1.Data(), slope2 / slope1));
    summary->AddText("Interpretation: system-level integrating-mode response test.");
    summary->Draw();

    c->Print(pdfName + ")");

    std::ofstream fitout("molleradc_integrating_linearity_fit_summary.csv");

    fitout << "Channel,Slope_ADC_per_Vpp,Intercept_ADC,SlopeError,InterceptError,Chi2,NDF\n";

    fitout << ch1 << ","
           << slope1 << ","
           << int1 << ","
           << fit1->GetParError(0) << ","
           << fit1->GetParError(1) << ","
           << fit1->GetChisquare() << ","
           << fit1->GetNDF() << "\n";

    fitout << ch2 << ","
           << slope2 << ","
           << int2 << ","
           << fit2->GetParError(0) << ","
           << fit2->GetParError(1) << ","
           << fit2->GetChisquare() << ","
           << fit2->GetNDF() << "\n";

    fitout.close();
    out.close();

    std::cout << "\nSaved PDF: " << pdfName << std::endl;
    std::cout << "Saved CSV: " << csvName << std::endl;
    std::cout << "Saved fit summary: molleradc_integrating_linearity_fit_summary.csv\n";
    std::cout << "Done.\n";
}
