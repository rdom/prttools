#include "PrtTools.h"

void proc_offset(TString in = "", int corrected = 0) {

  PrtTools t(in);
  int m = t.run()->getMomentum() + 0.1;
  int study = t.run()->getStudy();
  int fid = t.run()->getId();
  TString nid = t.run()->getName();
  bool bmc = t.run()->getMc();
  
  int h1a(0), h1b(50), h2a(0), h2b(50), hbin(1000); // h1a(200),h1b(400)

  if (corrected == 1) {
    h1a = 0;
    h1b = 50;
    h2a = 0;
    h2b = 50;
    hbin = 1000;
  }

  TH1F *hLeD = new TH1F("leD", "LE beam data ; LE [ns]; entries [#]", hbin, h1a, h1b);
  TH1F *hLeS = new TH1F("leS", "LE simulation; LE [ns]; entries [#]", hbin, h2a, h2b);
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  
  while (t.next() && t.i() < 1000000) {
    if (t.event()->getPid() != 2) continue;
    for (auto hit : t.event()->getHits()) {
      if (hit.getChannel() < t.maxdircch()) {
        double time = hit.getLeadTime();
	hLeD->Fill(time);
      }
    }
  }

  PrtTools ts(in.ReplaceAll("C.root","S.root"));    
  while (ts.next() && ts.i() < 1000000) {
    if (ts.event()->getPid() != 2) continue;
    for (auto hit : ts.event()->getHits()) {
      if (hit.getChannel() < ts.maxdircch()) {
        double time = hit.getLeadTime();
	time += gRandom->Gaus(0, 0.4);
        hLeS->Fill(time);
      }
    }
  }

  
  t.add_canvas("offset", 800, 400);

  hLeS->SetLineColor(kRed);
  TH1F *hh[] = {hLeS, hLeD};
  t.normalize_to(hh, 2, 1);
  hLeS->Draw("hist");
  hLeD->Draw("hist same");
  TLegend *leg = new TLegend(0.62, 0.7, 0.92, 0.9);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(hLeD, "beam data ", "lp");
  leg->AddEntry(hLeS, "simulation", "lp");
  leg->Draw();

  if (corrected == 0) {
    // double xmax1 = hLeD->GetXaxis()->GetBinCenter(hLeD->GetMaximumBin());
    // double xmax2 = hLeS->GetXaxis()->GetBinCenter(hLeS->GetMaximumBin());
    // double xmax1 = prt_fit(hLeD,0.5,50,1).X();
    // double xmax2 = prt_fit(hLeS,0.5,50,1).X();

    double threshold = hLeD->GetMaximum() * 0.6;
    int firstbin = hLeD->FindFirstBinAbove(threshold);
    double xmax1 = hLeD->GetXaxis()->GetBinCenter(firstbin);
    threshold = hLeS->GetMaximum() * 0.6;
    firstbin = hLeS->FindFirstBinAbove(threshold);
    double xmax2 = hLeS->GetXaxis()->GetBinCenter(firstbin);

    gPad->Update();
    TLine *gLine1 = new TLine(0, 0, 0, 1000);
    gLine1->SetX1(xmax1);
    gLine1->SetX2(xmax1);
    gLine1->SetY1(gPad->GetUymin());
    gLine1->SetY2(gPad->GetUymax());
    gLine1->SetLineColor(kBlue);
    gLine1->Draw();

    TLine *gLine2 = new TLine(0, 0, 0, 1000);
    gLine2->SetX1(xmax2);
    gLine2->SetX2(xmax2);
    gLine2->SetY1(gPad->GetUymin());
    gLine2->SetY2(gPad->GetUymax());
    gLine2->SetLineColor(kRed);
    gLine2->Draw();

    TString out = in.ReplaceAll(".root",".off.root");
        
    TFile efile(out, "RECREATE");
    TGraph *gr = new TGraph();
    gr->SetPoint(0, xmax1 - xmax2, xmax1 - xmax2);
    gr->SetName("off_" + nid);
    gr->Write();

    efile.Write();
    efile.Close();
    std::cout << "new offset " << xmax1 - xmax2 << std::endl;
  }

  t.save_canvas(t.dir(in) + Form("/%da/%d", study, fid), 0);
}
