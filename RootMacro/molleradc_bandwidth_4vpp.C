#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TF1.h>
#include <TH1D.h>
#include <TLine.h>
#include <TStyle.h>
#include <TMath.h>
#include <TColor.h>
#include <TPaveText.h>
#include <TLegend.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

struct RunInfo {
    int run;
    double freq;
};

struct FitResult {
    double amp = 0.0;
    double phase = 0.0;
    double offset = 0.0;
    double measuredRMS = 0.0;
    double expectedRMS = 0.0;
    double residualRMS = 0.0;
    double peakToPeak = 0.0;
    double nEventsUsed = 0.0;
};

void StyleGraph(TGraph *g, int color, int marker)
{
    g->SetMarkerStyle(marker);
    g->SetMarkerSize(1.2);
    g->SetMarkerColor(color);
    g->SetLineColor(color);
    g->SetLineWidth(2);
}

void DrawManualLegend(TGraph *g1,
                      TGraph *g2,
                      const char *label1,
                      const char *label2,
                      TLine *line = nullptr)
{
    TLegend *leg = new TLegend(0.60, 0.72, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.035);

    leg->AddEntry(g1, label1, "lp");
    leg->AddEntry(g2, label2, "lp");

    if (line)
        leg->AddEntry(line, "Unity response", "l");

    leg->Draw();
}

FitResult AnalyzeOneRunWindow(TTree *evt,
                              TString ch,
                              int run,
                              double freq,
                              double eventRate,
                              double nPeriods,
                              TCanvas *c,
                              TString pdfName,
                              bool firstPage)
{
    FitResult result;

    TString expr = ch + ".hw_sum";

    double periodEvents = eventRate / freq;
    double windowEvents = nPeriods * periodEvents;

    TString cut = Form("Entry$ < %.0f", windowEvents);

    Long64_t n = evt->Draw(Form("%s:Entry$", expr.Data()), cut, "goff");

    if (n <= 10) {
        std::cout << "ERROR: not enough entries for run "
                  << run << " channel " << ch << std::endl;
        return result;
    }

    result.nEventsUsed = n;

    TGraph *g = new TGraph();

    double sum = 0.0;
    double sum2 = 0.0;
    double ymin =  1e30;
    double ymax = -1e30;

    for (Long64_t i = 0; i < n; i++) {
        double y = evt->GetV1()[i];
        double x = evt->GetV2()[i];

        g->SetPoint(i, x, y);

        sum += y;
        sum2 += y * y;

        if (y < ymin) ymin = y;
        if (y > ymax) ymax = y;
    }

    double mean = sum / n;
    double rms = std::sqrt(sum2 / n - mean * mean);

    result.measuredRMS = rms;
    result.peakToPeak = ymax - ymin;

    double cyclesPerEvent = freq / eventRate;

    TF1 *fit = new TF1(Form("fit_%s_%d", ch.Data(), run),
                       "[0]*sin(2*TMath::Pi()*[1]*x+[2])+[3]",
                       0,
                       windowEvents);

    fit->SetParameters(rms * std::sqrt(2.0),
                       cyclesPerEvent,
                       0.0,
                       mean);

    fit->FixParameter(1, cyclesPerEvent);

    g->Fit(fit, "Q");

    result.amp = std::abs(fit->GetParameter(0));
    result.phase = fit->GetParameter(2);
    result.offset = fit->GetParameter(3);
    result.expectedRMS = result.amp / std::sqrt(2.0);

    TH1D *hResidual =
        new TH1D(Form("hResidual_%s_%d", ch.Data(), run),
                 Form("Run %d %s Residual;Residual (ADC Counts);Entries",
                      run, ch.Data()),
                 200,
                 -5.0 * rms,
                 5.0 * rms);

    for (Long64_t i = 0; i < n; i++) {
        double x, y;
        g->GetPoint(i, x, y);
        hResidual->Fill(y - fit->Eval(x));
    }

    result.residualRMS = hResidual->GetRMS();

    c->Clear();
    c->SetGridx(0);
    c->SetGridy(0);
    c->Divide(2, 2);

    c->cd(1);

    g->SetTitle(Form("Run %d, %.0f Hz, %s, first %.0f periods;Event Number;ADC Counts",
                     run, freq, ch.Data(), nPeriods));

    g->SetMarkerStyle(20);
    g->SetMarkerSize(0.7);
    g->SetMarkerColor(kBlack);
    g->SetLineColor(kBlack);
    g->Draw("AP");

    fit->SetLineColor(kRed);
    fit->SetLineWidth(3);
    fit->Draw("same");

    c->cd(2);

    hResidual->SetLineColor(kBlue + 1);
    hResidual->SetFillColor(0);
    hResidual->Draw();

    c->cd(3);

    TH1D *hData =
        new TH1D(Form("hData_%s_%d", ch.Data(), run),
                 Form("Run %d %s Data Distribution;ADC Counts;Entries",
                      run, ch.Data()),
                 200,
                 mean - 5.0 * rms,
                 mean + 5.0 * rms);

    for (Long64_t i = 0; i < n; i++) {
        double x, y;
        g->GetPoint(i, x, y);
        hData->Fill(y);
    }

    hData->SetLineColor(kBlack);
    hData->Draw();

    c->cd(4);

    TPaveText *summary = new TPaveText(0.10, 0.12, 0.90, 0.88, "NDC");
    summary->SetFillColor(0);
    summary->SetBorderSize(1);
    summary->SetTextAlign(12);
    summary->SetTextSize(0.035);

    double residualPercent = 0.0;
    if (result.measuredRMS > 0.0)
        residualPercent = 100.0 * result.residualRMS / result.measuredRMS;

    summary->AddText(Form("Run %d  |  %.0f Hz  |  %s", run, freq, ch.Data()));
    summary->AddText("----------------------------------------");
    summary->AddText(Form("Events used          : %.0f", result.nEventsUsed));
    summary->AddText(Form("Fit amplitude        : %.2f ADC counts", result.amp));
    summary->AddText(Form("Fit phase            : %.4f rad", result.phase));
    summary->AddText(Form("Fit offset           : %.2f ADC counts", result.offset));
    summary->AddText(Form("Measured RMS         : %.2f ADC counts", result.measuredRMS));
    summary->AddText(Form("A / sqrt(2)          : %.2f ADC counts", result.expectedRMS));
    summary->AddText(Form("Residual RMS         : %.2f ADC counts", result.residualRMS));
    summary->AddText(Form("Residual / Data RMS  : %.4f %%", residualPercent));
    summary->AddText(Form("Peak-to-peak         : %.2f ADC counts", result.peakToPeak));
    summary->Draw();

    if (firstPage)
        c->Print(pdfName + "(");
    else
        c->Print(pdfName);

    delete g;
    delete fit;
    delete hResidual;
    delete hData;
    delete summary;

    return result;
}

void molleradc_bandwidth_4vpp()
{
    gStyle->SetOptStat(1110);

    TString ch1 = "tq05_r1";
    TString ch2 = "tq05_r2";

    double eventRate = 960.0;
    double nPeriods = 10.0;

    double inputVpp = 4.0;
    double inputAmpVolts = inputVpp / 2.0;
    double expectedRMSVolts = inputAmpVolts / std::sqrt(2.0);
    double expectedPeakToPeakVolts = inputVpp;
    double expectedMeanVolts = 0.0;

    std::vector<RunInfo> runs = {
        {1371, 10.0},
        {1370, 15.0},
        {1369, 20.0},
        {1368, 50.0},
        {1367, 100.0}
    };

    TString pdfName = "molleradc_bandwidth_4vpp.pdf";
    TString csvName = "molleradc_bandwidth_4vpp_summary.csv";

    std::ofstream out(csvName.Data());

    out << "Run,Frequency_Hz,Channel,"
        << "EventsUsed,"
        << "FitAmplitude_ADC,"
        << "FitPhase,"
        << "FitOffset_ADC,"
        << "MeasuredRMS_ADC,"
        << "ExpectedRMS_FromFit_ADC,"
        << "ResidualRMS_ADC,"
        << "PeakToPeak_ADC,"
        << "AmplitudeRatioTo10Hz,"
        << "RMSAmplitude_Est_ADC,"
        << "RMSAmplitudeRatioTo10Hz,"
        << "ExpectedInputAmplitude_V,"
        << "ExpectedInputRMS_V,"
        << "ExpectedInputPeakToPeak_V,"
        << "ExpectedInputMean_V\n";

    std::cout << "\nMOLLERADC Bandwidth Test: Windowed Sine Fit\n";
    std::cout << "Channels: " << ch1 << ", " << ch2 << "\n";
    std::cout << "Input: 4 Vpp sine, zero offset\n";
    std::cout << "Event rate: " << eventRate << " Hz\n";
    std::cout << "Fit window: first " << nPeriods << " periods\n\n";

    TCanvas *c = new TCanvas("c",
                             "MOLLERADC Bandwidth Window Fit",
                             1800,
                             1100);

    std::vector<double> freqList;
    std::vector<double> amp1List;
    std::vector<double> amp2List;
    std::vector<double> ratio1List;
    std::vector<double> ratio2List;
    std::vector<double> rmsRatio1List;
    std::vector<double> rmsRatio2List;
    std::vector<double> residual1List;
    std::vector<double> residual2List;

    double refAmp1 = 0.0;
    double refAmp2 = 0.0;
    double refRmsAmp1 = 0.0;
    double refRmsAmp2 = 0.0;

    bool firstPage = true;

    for (size_t i = 0; i < runs.size(); i++) {
        int run = runs[i].run;
        double freq = runs[i].freq;

        TString fileName = Form("../RootFiles/isu_sample_%d.000.root", run);

        TFile *f = TFile::Open(fileName);

        if (!f || f->IsZombie()) {
            std::cout << "ERROR: Could not open " << fileName << std::endl;
            continue;
        }

        TTree *evt = (TTree*) f->Get("evt");

        if (!evt) {
            std::cout << "ERROR: Could not find evt tree in run " << run << std::endl;
            f->Close();
            continue;
        }

        if (!evt->GetBranch(ch1) || !evt->GetBranch(ch2)) {
            std::cout << "ERROR: Missing channel in run " << run << std::endl;
            f->Close();
            continue;
        }

        std::cout << "Analyzing run " << run
                  << " at " << freq << " Hz\n";

        FitResult r1 =
            AnalyzeOneRunWindow(evt, ch1, run, freq, eventRate,
                                nPeriods, c, pdfName, firstPage);

        firstPage = false;

        FitResult r2 =
            AnalyzeOneRunWindow(evt, ch2, run, freq, eventRate,
                                nPeriods, c, pdfName, false);

        double rmsAmp1 = r1.measuredRMS * std::sqrt(2.0);
        double rmsAmp2 = r2.measuredRMS * std::sqrt(2.0);

        if (freq == 10.0) {
            refAmp1 = r1.amp;
            refAmp2 = r2.amp;
            refRmsAmp1 = rmsAmp1;
            refRmsAmp2 = rmsAmp2;
        }

        double ratio1 = 0.0;
        double ratio2 = 0.0;
        double rmsRatio1 = 0.0;
        double rmsRatio2 = 0.0;

        if (refAmp1 > 0.0) ratio1 = r1.amp / refAmp1;
        if (refAmp2 > 0.0) ratio2 = r2.amp / refAmp2;
        if (refRmsAmp1 > 0.0) rmsRatio1 = rmsAmp1 / refRmsAmp1;
        if (refRmsAmp2 > 0.0) rmsRatio2 = rmsAmp2 / refRmsAmp2;

        freqList.push_back(freq);
        amp1List.push_back(r1.amp);
        amp2List.push_back(r2.amp);
        ratio1List.push_back(ratio1);
        ratio2List.push_back(ratio2);
        rmsRatio1List.push_back(rmsRatio1);
        rmsRatio2List.push_back(rmsRatio2);
        residual1List.push_back(r1.residualRMS);
        residual2List.push_back(r2.residualRMS);

        out << run << "," << freq << "," << ch1 << ","
            << r1.nEventsUsed << ","
            << r1.amp << ","
            << r1.phase << ","
            << r1.offset << ","
            << r1.measuredRMS << ","
            << r1.expectedRMS << ","
            << r1.residualRMS << ","
            << r1.peakToPeak << ","
            << ratio1 << ","
            << rmsAmp1 << ","
            << rmsRatio1 << ","
            << inputAmpVolts << ","
            << expectedRMSVolts << ","
            << expectedPeakToPeakVolts << ","
            << expectedMeanVolts << "\n";

        out << run << "," << freq << "," << ch2 << ","
            << r2.nEventsUsed << ","
            << r2.amp << ","
            << r2.phase << ","
            << r2.offset << ","
            << r2.measuredRMS << ","
            << r2.expectedRMS << ","
            << r2.residualRMS << ","
            << r2.peakToPeak << ","
            << ratio2 << ","
            << rmsAmp2 << ","
            << rmsRatio2 << ","
            << inputAmpVolts << ","
            << expectedRMSVolts << ","
            << expectedPeakToPeakVolts << ","
            << expectedMeanVolts << "\n";

        f->Close();
    }

    TGraph *gAmp1 = new TGraph();
    TGraph *gAmp2 = new TGraph();
    TGraph *gRatio1 = new TGraph();
    TGraph *gRatio2 = new TGraph();
    TGraph *gRmsRatio1 = new TGraph();
    TGraph *gRmsRatio2 = new TGraph();
    TGraph *gRes1 = new TGraph();
    TGraph *gRes2 = new TGraph();

    for (size_t i = 0; i < freqList.size(); i++) {
        gAmp1->SetPoint(i, freqList[i], amp1List[i]);
        gAmp2->SetPoint(i, freqList[i], amp2List[i]);
        gRatio1->SetPoint(i, freqList[i], ratio1List[i]);
        gRatio2->SetPoint(i, freqList[i], ratio2List[i]);
        gRmsRatio1->SetPoint(i, freqList[i], rmsRatio1List[i]);
        gRmsRatio2->SetPoint(i, freqList[i], rmsRatio2List[i]);
        gRes1->SetPoint(i, freqList[i], residual1List[i]);
        gRes2->SetPoint(i, freqList[i], residual2List[i]);
    }

    c->Clear();
    c->SetGridx();
    c->SetGridy();

    StyleGraph(gAmp1, kBlue + 1, 20);
    StyleGraph(gAmp2, kRed + 1, 21);

    gAmp1->SetTitle("Windowed Sine-Fit Amplitude vs Frequency;Frequency (Hz);Fit Amplitude (ADC Counts)");
    gAmp1->Draw("APL");
    gAmp1->GetXaxis()->SetLimits(0, 110);
    gAmp2->Draw("PL SAME");
    DrawManualLegend(gAmp1, gAmp2, "tq05_r1", "tq05_r2");
    c->Print(pdfName);

    c->Clear();
    c->SetGridx();
    c->SetGridy();

    StyleGraph(gRatio1, kBlue + 1, 20);
    StyleGraph(gRatio2, kRed + 1, 21);

    gRatio1->SetTitle("Windowed Sine-Fit Bandwidth Response;Frequency (Hz);Fit Amplitude / Fit Amplitude at 10 Hz");
    gRatio1->Draw("APL");
    gRatio1->GetXaxis()->SetLimits(0, 110);

    TLine *line1 = new TLine(0, 1.0, 110, 1.0);
    line1->SetLineStyle(2);
    line1->SetLineColor(kBlack);

    gRatio2->Draw("PL SAME");
    line1->Draw("same");
    DrawManualLegend(gRatio1, gRatio2, "tq05_r1", "tq05_r2", line1);
    c->Print(pdfName);

    c->Clear();
    c->SetGridx();
    c->SetGridy();

    StyleGraph(gRmsRatio1, kBlue + 1, 20);
    StyleGraph(gRmsRatio2, kRed + 1, 21);

    gRmsRatio1->SetTitle("RMS-Derived Bandwidth Response;Frequency (Hz);#sqrt{2} RMS / (#sqrt{2} RMS at 10 Hz)");
    gRmsRatio1->Draw("APL");
    gRmsRatio1->GetXaxis()->SetLimits(0, 110);

    TLine *line2 = new TLine(0, 1.0, 110, 1.0);
    line2->SetLineStyle(2);
    line2->SetLineColor(kBlack);

    gRmsRatio2->Draw("PL SAME");
    line2->Draw("same");
    DrawManualLegend(gRmsRatio1, gRmsRatio2, "tq05_r1", "tq05_r2", line2);
    c->Print(pdfName);

    c->Clear();
    c->SetGridx();
    c->SetGridy();

    StyleGraph(gRes1, kBlue + 1, 20);
    StyleGraph(gRes2, kRed + 1, 21);

    gRes1->SetTitle("Windowed Sine-Fit Residual RMS vs Frequency;Frequency (Hz);Residual RMS (ADC Counts)");
    gRes1->Draw("APL");
    gRes1->GetXaxis()->SetLimits(0, 110);
    gRes2->Draw("PL SAME");
    DrawManualLegend(gRes1, gRes2, "tq05_r1", "tq05_r2");
    c->Print(pdfName + ")");

    c->SetGridx(0);
    c->SetGridy(0);

    out.close();

    std::cout << "\nSaved PDF: " << pdfName << std::endl;
    std::cout << "Saved CSV: " << csvName << std::endl;
    std::cout << "Done.\n";
}
