#include <TRandom3.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TString.h>
#include <TStyle.h>
#include <TMath.h>

#include <iostream>
#include <vector>
#include <cmath>

void molleradc_trandom_minmax_cartoon()
{
    gStyle->SetOptStat(1110);

    TRandom3 r(0);

    TString pdfName = "molleradc_trandom_minmax_cartoon.pdf";

    const int nTrials = 10000;

    std::vector<int> sampleSizes = {100, 500, 1000, 2000, 4000};

    TGraph *gMeanMax = new TGraph();
    TGraph *gMeanMin = new TGraph();
    TGraph *gMeanAbsMax = new TGraph();

    std::cout << "\nTRandom Gaussian Min/Max Cartoon Model\n";
    std::cout << "Trials per sample size: " << nTrials << "\n\n";

    std::cout << "Nsamples"
              << "\tMeanMax"
              << "\tMeanMin"
              << "\tMeanLargestAbs"
              << std::endl;

    for (size_t is = 0; is < sampleSizes.size(); is++)
    {
        int nSamples = sampleSizes[is];

        double sumMax = 0.0;
        double sumMin = 0.0;
        double sumLargestAbs = 0.0;

        for (int t = 0; t < nTrials; t++)
        {
            double xmin =  1e9;
            double xmax = -1e9;

            for (int i = 0; i < nSamples; i++)
            {
                double x = r.Gaus(0.0, 1.0);

                if (x < xmin) xmin = x;
                if (x > xmax) xmax = x;
            }

            sumMax += xmax;
            sumMin += xmin;
            sumLargestAbs += std::max(std::abs(xmax), std::abs(xmin));
        }

        double meanMax = sumMax / nTrials;
        double meanMin = sumMin / nTrials;
        double meanLargestAbs = sumLargestAbs / nTrials;

        gMeanMax->SetPoint(is, nSamples, meanMax);
        gMeanMin->SetPoint(is, nSamples, meanMin);
        gMeanAbsMax->SetPoint(is, nSamples, meanLargestAbs);

        std::cout << nSamples
                  << "\t" << meanMax
                  << "\t" << meanMin
                  << "\t" << meanLargestAbs
                  << std::endl;
    }

    TCanvas *c = new TCanvas("c", "TRandom Min/Max Cartoon", 1600, 1000);

    // ------------------------------------------------------------
    // Page 1: Expected max/min vs number of samples
    // ------------------------------------------------------------
    c->Clear();

    gMeanMax->SetTitle("Gaussian Cartoon: Expected Min/Max vs Sample Count;Samples per Block;Distance from Mean in #sigma");
    gMeanMax->SetMarkerStyle(20);
    gMeanMax->SetMarkerSize(1.2);
    gMeanMax->SetLineWidth(2);

    gMeanMin->SetMarkerStyle(21);
    gMeanMin->SetMarkerSize(1.2);
    gMeanMin->SetLineWidth(2);

    gMeanMax->Draw("APL");
    gMeanMin->Draw("PL same");

    c->BuildLegend();
    c->Print(pdfName + "(");

    // ------------------------------------------------------------
    // Page 2: Largest absolute excursion vs number of samples
    // ------------------------------------------------------------
    c->Clear();

    gMeanAbsMax->SetTitle("Gaussian Cartoon: Largest Absolute Excursion vs Sample Count;Samples per Block;Largest |excursion| in #sigma");
    gMeanAbsMax->SetMarkerStyle(20);
    gMeanAbsMax->SetMarkerSize(1.2);
    gMeanAbsMax->SetLineWidth(2);
    gMeanAbsMax->Draw("APL");

    c->Print(pdfName);

    // ------------------------------------------------------------
    // Page 3: Distribution for real-like block size, N = 2000
    // ------------------------------------------------------------
    int nReal = 2000;

    TH1D *hMax2000 = new TH1D("hMax2000",
                              "N = 2000 Gaussian Samples: Max Distribution;Max value in #sigma;Entries",
                              120, 2.0, 5.5);

    TH1D *hMin2000 = new TH1D("hMin2000",
                              "N = 2000 Gaussian Samples: Min Distribution;Min value in #sigma;Entries",
                              120, -5.5, -2.0);

    TH1D *hLargest2000 = new TH1D("hLargest2000",
                                  "N = 2000 Gaussian Samples: Largest |Min/Max|;Largest |excursion| in #sigma;Entries",
                                  120, 2.0, 6.0);

    for (int t = 0; t < nTrials; t++)
    {
        double xmin =  1e9;
        double xmax = -1e9;

        for (int i = 0; i < nReal; i++)
        {
            double x = r.Gaus(0.0, 1.0);

            if (x < xmin) xmin = x;
            if (x > xmax) xmax = x;
        }

        hMax2000->Fill(xmax);
        hMin2000->Fill(xmin);
        hLargest2000->Fill(std::max(std::abs(xmax), std::abs(xmin)));
    }

    c->Clear();
    c->Divide(2, 2);

    c->cd(1);
    hMax2000->Draw();

    c->cd(2);
    hMin2000->Draw();

    c->cd(3);
    hLargest2000->Draw();

    c->Print(pdfName);

    // ------------------------------------------------------------
    // Page 4: Running maximum over samples
    // ------------------------------------------------------------
    int nRun = 4000;

    TGraph *gRunMax = new TGraph();
    TGraph *gRunMin = new TGraph();
    TGraph *gRunAbs = new TGraph();

    double runningMax = -1e9;
    double runningMin =  1e9;

    for (int i = 0; i < nRun; i++)
    {
        double x = r.Gaus(0.0, 1.0);

        if (x > runningMax) runningMax = x;
        if (x < runningMin) runningMin = x;

        gRunMax->SetPoint(i, i + 1, runningMax);
        gRunMin->SetPoint(i, i + 1, runningMin);
        gRunAbs->SetPoint(i, i + 1, std::max(std::abs(runningMax), std::abs(runningMin)));
    }

    c->Clear();

    gRunAbs->SetTitle("Gaussian Cartoon: Running Largest Excursion;Sample Number;Running largest |excursion| in #sigma");
    gRunAbs->SetMarkerStyle(20);
    gRunAbs->SetMarkerSize(0.4);
    gRunAbs->SetLineWidth(2);
    gRunAbs->Draw("AL");

    c->Print(pdfName);

    // ------------------------------------------------------------
    // Page 5: First block vs full window cartoon
    // ------------------------------------------------------------
    int firstBlockSamples = 2000;
    int fullWindowSamples = 8000;

    TH1D *hFirstBlock = new TH1D("hFirstBlock",
                                 "First Block: N = 2000;Largest |excursion| in #sigma;Entries",
                                 120, 2.0, 6.0);

    TH1D *hFullWindow = new TH1D("hFullWindow",
                                 "Full Window: N = 8000;Largest |excursion| in #sigma;Entries",
                                 120, 2.0, 6.5);

    for (int t = 0; t < nTrials; t++)
    {
        double minBlock =  1e9;
        double maxBlock = -1e9;

        for (int i = 0; i < firstBlockSamples; i++)
        {
            double x = r.Gaus(0.0, 1.0);

            if (x < minBlock) minBlock = x;
            if (x > maxBlock) maxBlock = x;
        }

        hFirstBlock->Fill(std::max(std::abs(maxBlock), std::abs(minBlock)));

        double minFull =  1e9;
        double maxFull = -1e9;

        for (int i = 0; i < fullWindowSamples; i++)
        {
            double x = r.Gaus(0.0, 1.0);

            if (x < minFull) minFull = x;
            if (x > maxFull) maxFull = x;
        }

        hFullWindow->Fill(std::max(std::abs(maxFull), std::abs(minFull)));
    }

    c->Clear();
    c->Divide(2, 1);

    c->cd(1);
    hFirstBlock->Draw();

    c->cd(2);
    hFullWindow->Draw();

    c->Print(pdfName + ")");

    std::cout << "\nFor N = 2000 samples per block:\n";
    std::cout << "  Mean max = " << hMax2000->GetMean() << " sigma\n";
    std::cout << "  Mean largest |min/max| = " << hLargest2000->GetMean() << " sigma\n";

    std::cout << "\nSaved PDF: " << pdfName << std::endl;
    std::cout << "Done.\n";
}
