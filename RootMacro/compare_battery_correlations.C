
// compare_battery_correlations.C
// Usage:
//   root -l
//   .x compare_battery_correlations.C("../RootFiles/isu_sample_1222.411.root")

#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH2.h>
#include <TH1.h>
#include <TF1.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TDirectory.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using namespace std;

struct CorrPlot {
  TString tree_name;
  TString xvar;
  TString yvar;
  TString label;
};

void compare_battery_correlations(
    TString filename = "../RootFiles/isu_sample_1222.411.root",
    TString output_prefix = "battery_correlations"
)
{
  gStyle->SetOptFit(1111);
  gStyle->SetOptStat(1111);

  TFile *file = TFile::Open(filename, "READ");
  if (!file || file->IsZombie()) {
    cout << "ERROR: Cannot open file " << filename << endl;
    return;
  }

  TTree *pr  = (TTree*) file->Get("pr");
  TTree *mul = (TTree*) file->Get("mul");

  if (!pr || !mul) {
    cout << "ERROR: Could not find pr or mul trees." << endl;
    file->Close();
    return;
  }

  TString pdfname  = output_prefix + ".pdf";
  TString rootname = output_prefix + ".root";
  TString csvname  = output_prefix + ".csv";

  TFile *outfile = TFile::Open(rootname, "RECREATE");
  ofstream csv(csvname.Data());

  csv << "tree,plot_name,entries,"
      << "mean_x,mean_y,rms_x,rms_y,corr_factor,"
      << "fit_p0,fit_p0_error,fit_p1,fit_p1_error,fit_chi2,fit_ndf,"
      << "residual_mean,residual_rms,residual_stddev,residual_entries"
      << endl;

  vector<CorrPlot> plots = {

    // Hall batteries, pair tree
    {"pr",  "asym_battery_hall_1.hw_sum", "asym_battery_hall_2.hw_sum", "pr_hall_batt1_vs_batt2"},
    {"pr",  "asym_battery_hall_1.hw_sum", "asym_battery_hall_3.hw_sum", "pr_hall_batt1_vs_batt3"},
    {"pr",  "asym_battery_hall_2.hw_sum", "asym_battery_hall_3.hw_sum", "pr_hall_batt2_vs_batt3"},

    // ISB batteries, pair tree
    {"pr",  "asym_rack_battery_isb.hw_sum", "asym_inj_laser_table_battery_3.hw_sum", "pr_isb_rack_vs_laser_table"},

    // Hall batteries, mul tree
    {"mul", "asym_battery_hall_1.hw_sum", "asym_battery_hall_2.hw_sum", "mul_hall_batt1_vs_batt2"},
    {"mul", "asym_battery_hall_1.hw_sum", "asym_battery_hall_3.hw_sum", "mul_hall_batt1_vs_batt3"},
    {"mul", "asym_battery_hall_2.hw_sum", "asym_battery_hall_3.hw_sum", "mul_hall_batt2_vs_batt3"},

    // ISB batteries, mul tree
    {"mul", "asym_rack_battery_isb.hw_sum", "asym_inj_laser_table_battery_3.hw_sum", "mul_isb_rack_vs_laser_table"}
  };

  TCanvas *c1 = new TCanvas("c1", "Battery Correlations", 1400, 700);

  c1->Print(pdfname + "[");

  for (size_t i = 0; i < plots.size(); i++) {

    TTree *tree = nullptr;
    if (plots[i].tree_name == "pr") tree = pr;
    if (plots[i].tree_name == "mul") tree = mul;

    if (!tree) continue;

    if (!tree->GetBranch(plots[i].xvar.Tokenize(".")->At(0)->GetName()) ||
        !tree->GetBranch(plots[i].yvar.Tokenize(".")->At(0)->GetName())) {
      cout << "WARNING: Missing branch for " << plots[i].label << endl;
      continue;
    }

    TString h2name = "h2_" + plots[i].label;
    TString fitname = "fit_" + plots[i].label;
    TString hresname = "hres_" + plots[i].label;
    TString cname = "canvas_" + plots[i].label;

    TString drawcmd =
      plots[i].yvar + ":" + plots[i].xvar +
      ">>" + h2name + "(300, -2e-5, 2e-5, 300, -2e-5, 2e-5)";

    cout << "Drawing correlation: " << drawcmd << endl;

    c1->Clear();
    file->cd();
    tree->Draw(drawcmd, "", "COLZ");

    TH2 *h2 = dynamic_cast<TH2*>(gROOT->FindObject(h2name));
    if (!h2) {
      h2 = dynamic_cast<TH2*>(gDirectory->Get(h2name));
    }

    if (!h2) {
      cout << "WARNING: Could not create histogram " << h2name << endl;
      continue;
    }

    h2->SetTitle(plots[i].label + ";" + plots[i].xvar + ";" + plots[i].yvar);

    TF1 *fit = new TF1(
      fitname,
      "[0] + [1]*x",
      h2->GetXaxis()->GetXmin(),
      h2->GetXaxis()->GetXmax()
    );

    h2->Fit(fit, "Q");

    double entries = h2->GetEntries();
    double meanx   = h2->GetMean(1);
    double meany   = h2->GetMean(2);
    double rmsx    = h2->GetRMS(1);
    double rmsy    = h2->GetRMS(2);
    double corr    = h2->GetCorrelationFactor();

    double p0      = fit->GetParameter(0);
    double p1      = fit->GetParameter(1);
    double p0err   = fit->GetParError(0);
    double p1err   = fit->GetParError(1);
    double chi2    = fit->GetChisquare();
    double ndf     = fit->GetNDF();

    TString residual_expr = Form(
      "%s - (%0.15g + %0.15g*%s)",
      plots[i].yvar.Data(), p0, p1, plots[i].xvar.Data()
    );

    TString hrescmd = residual_expr + ">>" + hresname + "(300)";

    cout << "Drawing residual: " << hrescmd << endl;

    tree->Draw(hrescmd, "", "goff");

    TH1 *hres = dynamic_cast<TH1*>(gROOT->FindObject(hresname));
    if (!hres) {
      hres = dynamic_cast<TH1*>(gDirectory->Get(hresname));
    }

    if (!hres) {
      cout << "WARNING: Could not create residual histogram "
           << hresname << endl;
      continue;
    }

    hres = dynamic_cast<TH1*>(hres->Clone(hresname + "_saved"));
    hres->SetDirectory(outfile);

    hres->SetTitle(plots[i].label + " residual;" +
                   plots[i].yvar + " - (p0 + p1*" + plots[i].xvar + ");Entries");

    double res_mean    = hres->GetMean();
    double res_rms     = hres->GetRMS();
    double res_stddev  = hres->GetStdDev();
    double res_entries = hres->GetEntries();

    c1->Clear();
    c1->Divide(2,1);

    c1->cd(1);
    h2->Draw("COLZ");
    fit->Draw("same");

    c1->cd(2);
    hres->Draw();

    c1->Print(pdfname);

    outfile->cd();

    TH2 *h2save = dynamic_cast<TH2*>(h2->Clone(h2name + "_saved"));
    h2save->SetDirectory(outfile);
    h2save->Write();

    fit->Write();
    hres->Write();

    c1->Write(cname);

    csv << plots[i].tree_name << ","
        << plots[i].label << ","
        << entries << ","
        << meanx << ","
        << meany << ","
        << rmsx << ","
        << rmsy << ","
        << corr << ","
        << p0 << ","
        << p0err << ","
        << p1 << ","
        << p1err << ","
        << chi2 << ","
        << ndf << ","
        << res_mean << ","
        << res_rms << ","
        << res_stddev << ","
        << res_entries
        << endl;

    cout << "Processed " << plots[i].label << endl;
    cout << "  Fit: y = " << p0 << " + " << p1 << "*x" << endl;
    cout << "  Residual mean = " << res_mean
         << ", RMS = " << res_rms << endl;
  }

  c1->Print(pdfname + "]");

  csv.close();
  outfile->Write();
  outfile->Close();
  file->Close();

  cout << endl;
  cout << "==========================================" << endl;
  cout << "Done." << endl;
  cout << "Saved PDF  : " << pdfname  << endl;
  cout << "Saved ROOT : " << rootname << endl;
  cout << "Saved CSV  : " << csvname  << endl;
  cout << "==========================================" << endl;
}
