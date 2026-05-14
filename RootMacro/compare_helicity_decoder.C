// compare_helicity_decoder.C
// Compare helicity decoder signals from two ADCs/ROCs.

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TH1.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TDirectory.h"

void compare_helicity_decoder(
    const char* input_root_file,
    const char* tree_name = "evt",
    const char* output_prefix = "helicity_decoder_comparison"
)
{
  gStyle->SetOptFit(1111);
  gStyle->SetOptStat(1111);

  TFile* infile = TFile::Open(input_root_file, "READ");
  if (!infile || infile->IsZombie()) {
    std::cerr << "ERROR: Cannot open input ROOT file: "
              << input_root_file << std::endl;
    return;
  }

  TTree* tree = dynamic_cast<TTree*>(infile->Get(tree_name));
  if (!tree) {
    std::cerr << "ERROR: Cannot find tree named "
              << tree_name << " in file." << std::endl;
    infile->Close();
    return;
  }

  std::vector<std::string> signals = {
    "pair_sync",
    "reported_helicity",
    "pattern_sync",
    "tsettle"
  };

  TString pdf_name  = TString(output_prefix) + ".pdf";
  TString csv_name  = TString(output_prefix) + ".csv";
  TString root_name = TString(output_prefix) + ".root";

  TFile* outfile = TFile::Open(root_name, "RECREATE");
  if (!outfile || outfile->IsZombie()) {
    std::cerr << "ERROR: Cannot create output ROOT file: "
              << root_name << std::endl;
    infile->Close();
    return;
  }

  std::ofstream csv(csv_name.Data());
  if (!csv.is_open()) {
    std::cerr << "ERROR: Cannot create CSV file: "
              << csv_name << std::endl;
    outfile->Close();
    infile->Close();
    return;
  }

  csv << "signal,"
      << "entries,"
      << "fit_p0,"
      << "fit_p0_error,"
      << "fit_p1,"
      << "fit_p1_error,"
      << "fit_chi2,"
      << "fit_ndf,"
      << "residual_mean,"
      << "residual_rms,"
      << "residual_stddev,"
      << "residual_entries"
      << "\n";

  TCanvas* c = new TCanvas("c", "Helicity Decoder Comparison", 1200, 900);

  c->Print(pdf_name + "[");

  for (const auto& sig : signals) {

    TString yvar = TString(sig) + "_1";
    TString xvar = TString(sig) + "_2";

    if (!tree->GetBranch(yvar) || !tree->GetBranch(xvar)) {
      std::cerr << "WARNING: Missing branch pair: "
                << yvar << " or " << xvar
                << ". Skipping." << std::endl;
      continue;
    }

    TString graph_name = "g_" + TString(sig);
    TString fit_name   = "fit_" + TString(sig);

    TString draw_expr = yvar + ":" + xvar;
    tree->Draw(draw_expr, "", "goff");

    Long64_t n = tree->GetSelectedRows();
    if (n <= 0) {
      std::cerr << "WARNING: No entries for " << sig
                << ". Skipping." << std::endl;
      continue;
    }

    TGraph* graph = new TGraph(n, tree->GetV2(), tree->GetV1());
    graph->SetName(graph_name);
    graph->SetTitle(TString(sig) + Form("; %s; %s", xvar.Data(), yvar.Data()));
    graph->SetMarkerStyle(20);
    graph->SetMarkerSize(0.45);

    double xmin = graph->GetXaxis()->GetXmin();
    double xmax = graph->GetXaxis()->GetXmax();

    TF1* fit = new TF1(fit_name, "[0] + [1]*x", xmin, xmax);
    graph->Fit(fit, "Q");

    double p0    = fit->GetParameter(0);
    double p1    = fit->GetParameter(1);
    double p0err = fit->GetParError(0);
    double p1err = fit->GetParError(1);
    double chi2  = fit->GetChisquare();
    int ndf      = fit->GetNDF();

    c->Clear();
    graph->Draw("AP");
    fit->Draw("same");
    c->Print(pdf_name);

    outfile->cd();
    graph->Write();
    fit->Write();

    TString residual_expr = Form(
      "%s - (%0.15g + %0.15g*%s)",
      yvar.Data(), p0, p1, xvar.Data()
    );

    TString hres_name = "hres_" + TString(sig);
    TString hres_cmd = residual_expr + ">>" + hres_name + "(200)";

    c->Clear();
    infile->cd();
    tree->Draw(hres_cmd, "", "");

    TH1* hres_tmp = dynamic_cast<TH1*>(gPad->GetPrimitive(hres_name));
    if (!hres_tmp) {
      hres_tmp = dynamic_cast<TH1*>(gDirectory->Get(hres_name));
    }

    if (!hres_tmp) {
      std::cerr << "WARNING: Could not create residual histogram for "
                << sig << std::endl;

      csv << sig << ","
          << n << ","
          << p0 << ","
          << p0err << ","
          << p1 << ","
          << p1err << ","
          << chi2 << ","
          << ndf << ","
          << "nan,nan,nan,0"
          << "\n";

      continue;
    }

    TH1* hres = dynamic_cast<TH1*>(hres_tmp->Clone(hres_name + "_saved"));
    hres->SetDirectory(outfile);

    hres->SetTitle(TString(sig) + " residual; "
                   + yvar + " - (p0 + p1*" + xvar + "); Entries");

    double res_mean    = hres->GetMean();
    double res_rms     = hres->GetRMS();
    double res_stddev  = hres->GetStdDev();
    double res_entries = hres->GetEntries();

    c->Clear();
    hres->Draw();
    c->Print(pdf_name);

    outfile->cd();
    hres->Write();

    csv << sig << ","
        << n << ","
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
        << "\n";

    std::cout << "Processed " << sig << std::endl;
    std::cout << "  Fit: " << yvar << " = "
              << p0 << " + " << p1 << "*" << xvar << std::endl;
    std::cout << "  Residual mean = " << res_mean
              << ", RMS = " << res_rms << std::endl;
  }

  c->Print(pdf_name + "]");

  csv.close();
  outfile->Write();
  outfile->Close();
  infile->Close();

  std::cout << "\nDone.\n";
  std::cout << "Saved PDF:  " << pdf_name << std::endl;
  std::cout << "Saved CSV:  " << csv_name << std::endl;
  std::cout << "Saved ROOT: " << root_name << std::endl;
}
