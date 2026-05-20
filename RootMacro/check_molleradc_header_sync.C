
//root -l
//.x check_molleradc_header_sync.C("../RootFiles/isu_sample_1212.000.root")

void check_molleradc_header_sync(
    const char* infile,
    const char* treename = "evt",
    const char* outbase = "molleradc_header_sync")
{
  TFile *fin = TFile::Open(infile, "READ");
  if (!fin || fin->IsZombie()) {
    std::cerr << "Cannot open input file: " << infile << std::endl;
    return;
  }

  TTree *evt = (TTree*) fin->Get(treename);
  if (!evt) {
    std::cerr << "Cannot find tree: " << treename << std::endl;
    return;
  }

  TString pdfname  = TString(outbase) + ".pdf";
  TString rootname = TString(outbase) + ".root";

  TFile *fout = TFile::Open(rootname, "RECREATE");
  TCanvas *c = new TCanvas("c", "MollerADC Header Sync Checks", 1000, 750);

  c->Print(pdfname + "[");

  auto draw_and_save = [&](const char* expr, const char* title) {
    c->Clear();
    evt->Draw(expr);
    TH1 *h = (TH1*) gPad->GetPrimitive("htemp");
    if (h) {
      h->SetTitle(title);
      h->Write(title);
    }
    c->Write(title);
    c->Print(pdfname);
  };

  draw_and_save("tq01_r1.region_timestamp - tq01_r2.region_timestamp",
                "BCM timestamp residual tq01_r1_minus_tq01_r2");

  draw_and_save("tq01_r1.header_packet_count - tq01_r2.header_packet_count",
                "BCM packet residual tq01_r1_minus_tq01_r2");

  draw_and_save("tq01_r1.region_timestamp:tq01_r2.region_timestamp",
                "BCM timestamp correlation tq01_r1_vs_tq01_r2");

  draw_and_save("tq01_r1.region_timestamp - tq02_r1.region_timestamp",
                "BCM timestamp residual tq01_r1_minus_tq02_r1");

  draw_and_save("tq01_r1.header_packet_count - tq02_r1.header_packet_count",
                "BCM packet residual tq01_r1_minus_tq02_r1");

  draw_and_save("bpm2i00XP.region_timestamp - bpm2i00XM.region_timestamp",
                "BPM timestamp residual bpm2i00XP_minus_bpm2i00XM");

  draw_and_save("bpm2i00XP.header_packet_count - bpm2i00XM.header_packet_count",
                "BPM packet residual bpm2i00XP_minus_bpm2i00XM");

  draw_and_save("bpm2i00XP.region_timestamp:bpm2i00XM.region_timestamp",
                "BPM timestamp correlation bpm2i00XP_vs_bpm2i00XM");

  c->Print(pdfname + "]");

  fout->Write();
  fout->Close();
  fin->Close();

  std::cout << "Saved PDF:  " << pdfname << std::endl;
  std::cout << "Saved ROOT: " << rootname << std::endl;
}
