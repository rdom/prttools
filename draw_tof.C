#include "PrtTools.h"

double walktheta(-4 * TMath::Pi() / 180.);
double tof1le(0), tof2le(0), tof1tot(0), tof2tot(0);

double fr11[11] = {0, 0.5, 0.5, 0.3, 0.3, 0.4, 0.3, 0.3, 0.2, 0.20, 0.15};
double fr12[11] = {0, 1.0, 1.0, 0.9, 0.9, 0.9, 0.9, 0.9, 0.8, 0.80, 0.70};
double fr21[11] = {0, 0.8, 0.8, 0.3, 0.3, 0.4, 0.3, 0.3, 0.2, 0.2, 0.2};
double fr22[11] = {0, 1.0, 1.0, 0.9, 0.9, 0.9, 0.9, 0.9, 0.8, 0.8, 0.8};

double c1y(0.5), c2y(0.5), c1x(0.9), c2x(0.9);

bool insideOfEllipce(double x, double y, double x0, double y0, double r1, double r2, double w = 0) {

  double xx = cos(w) * (x - x0) + sin(w) * (y - y0);
  double yy = sin(w) * (x - x0) - cos(w) * (y - y0);

  return xx * xx / (r1 * r1) + yy * yy / (r2 * r2) <= 1;
}

void draw_tof(TString infile = "hits.root", TString gcFile = "0") {
  if (gcFile != "0") {
    TFile f(gcFile);
    TIter nextkey(f.GetListOfKeys());
    TKey *key;

    while ((key = (TKey *)nextkey())) {
      TGraph *gr = (TGraph *)key->ReadObj();
      TString name = gr->GetName();
      if (name.Contains("tof")) {
        name.Remove(0, 4);
        if (infile.Contains(name)) {
          gr->GetPoint(0, tof1le, tof2le);
          gr->GetPoint(1, tof1tot, tof2tot);
        }
      }
    }
    f.Close();
  }

  PrtTools t(infile);
  
  int m = t.run()->getMomentum()+0.1;
  int study = t.run()->getStudy();

  // MCP-TOF
  // double le1(31), le2(36), l1(30), l2(35);
  // double t11(42), t12(70), t21(45), t22(52);

  // double le1(34), le2(40), l1(34), l2(40);
  // double t11(36), t12(45), t21(36), t22(45);
  double le1(32), le2(36), l1(32), l2(36);
  double t11(43), t12(55), t21(45), t22(55);

  // // SiTil
  // double le1(32), le2(36),l1(32), l2(36);
  // double t11(46), t12(50), t21(35), t22(39);

  if (m < 7) {
    le2 = 80;
    l2 = 80;
  }
  if (m == 2) {
    le2 = 85;
    l2 = 85;
  }
  std::cout<<"m "<<m<<std::endl;
  
  
  c1y = fr11[m];
  c2y = fr21[m];
  c1x = fr12[m];
  c2x = fr22[m];

  TH1F *hMult[4];

  TString names[] = {"trigger_1", "trigger_2", "tof_1", "tof_2"};
  int colors[] = {2, 4, 1, 8};
  for (int i = 0; i < 4; i++) {
    hMult[i] = new TH1F(names[i], ";multiplicity [#]; entries [#]", 10, 0, 10);
    hMult[i]->SetLineColor(colors[i]);
  }

  TH1F *hTof1 = new TH1F("tof1 ", "tof1;TOF2-TOF1 [ns]; entries [#]", 600, -1000, 1000);
  TH1F *hTof2 = new TH1F("tof2 ", "tof2;TOF2-TOF1 [ns]; entries [#]", 600, -1000, 1000);
  TH1F *hTof = new TH1F("tof ", "tof;TOF2-TOF1 [ns]; entries [#]", 600, le1, le2);
  TH1F *hTofC = new TH1F("tofC ", "tofC;TOF2-TOF1 [ns]; entries [#]", 600, le1, le2);
  TH1F *hTot = new TH1F("tot ", "tot;TOT1,TOT2 [ns]; entries [#]", 600, 0, 100);

  TH2F *hLeTot1 =
    new TH2F("letot1 ", "letot1;TOF2-TOF1 [ns]; TOT1 [ns]", 900, l1, l2, 125, t11, t12);
  TH2F *hLeTot2 =
    new TH2F("letot2 ", "letot2;TOF2-TOF1 [ns]; TOT2 [ns]", 500, l1, l2, 125, t21, t22);
  TH2F *hLeTotC =
    new TH2F("letotC ", "letotC;TOF2-TOF1 [ns]; TOT1 [ns]", 900, l1, l2, 125, t11, t12);
  TH2F *hLeTotC2 =
    new TH2F("letotC2 ", "letotC2;TOF2-TOF1 [ns]; TOT2 [ns]", 500, l1, l2, 125, t21, t22);

  gStyle->SetOptStat(1001111);
  gStyle->SetOptFit();

  while (t.next() && t.i() < 1000000) {
    
    bool btrig(false), bmcpout(false), btof1(false), btof2(false), str1b(false), stl1b(false),
      stl2b(false), str2b(false), str3b(false), stl3b(false);
    double tot1(0), tot2(0), tof1(0), tof2(0), str1l(0), stl1l(0), str3l(0), stl3l(0), str1t(0),
      stl1t(0), str3t(0), stl3t(0), str2l(0), stl2l(0), str2t(0), stl2t(0);
    int mult1(0), mult2(0), mult3(0), mult4(0);

    for (auto hit : t.event()->getHits()) {

      if (hit.getChannel() == 520) { // trigger 1
        btrig = true;
        mult1++;
      }
      if (hit.getChannel() == 513) { // trigger 2
        bmcpout = true;
        mult2++;
      }

      if (hit.getChannel() == 1136 && tof1 == 0) { // tof1
        btof1 = true;
        tof1 = hit.getLeadTime();
        tot1 = hit.getTotTime();
        mult3++;
      }
      if (hit.getChannel() == 1138 && tof2 == 0) { // tof2
        btof2 = true;
        tof2 = hit.getLeadTime();
        tot2 = hit.getTotTime();
        mult4++;
      }

      if (hit.getChannel() == 1140 && !str1b) {
        str1b = true;
        str1l = hit.getLeadTime();
        str1t = hit.getTotTime();
      }
      if (hit.getChannel() == 1142 && !stl1b) {
        stl1b = true;
        stl1l = hit.getLeadTime();
        stl1t = hit.getTotTime();
      }

      if (hit.getChannel() == 1144 && !str2b) {
        str2b = true;
        str2l = hit.getLeadTime();
        str2t = hit.getTotTime();
      }
      if (hit.getChannel() == 1146 && !stl2b) {
        stl2b = true;
        stl2l = hit.getLeadTime();
        stl2t = hit.getTotTime();
      }

      if (hit.getChannel() == 1148 && !str3b) {
        str3b = true;
        str3l = hit.getLeadTime();
        str3t = hit.getTotTime();
      }
      if (hit.getChannel() == 1150 && !stl3b) {
        stl3b = true;
        stl3l = hit.getLeadTime();
        stl3t = hit.getTotTime();
      }
    }

    if (!(btrig && btof1 && btof2)) continue;
    if (!(btof1 && btof2)) continue;

    // //    if(!(str1b && stl1b && str2b && stl2b)) continue;
    // if(!(str1b && stl1b && str2b && stl2b && str3b && stl3b)) continue;
    // tof1=(str1l+stl1l)/2.;
    // tof2=(str2l+stl2l)/2.;

    // tot1=str1t;
    // tot2=str3t;

    // if(fabs(tot1-tof1tot)>0.5 ||fabs(tot2-tof2tot)>0.5 )  continue;

    hMult[0]->Fill(mult1);
    hMult[1]->Fill(mult2);
    hMult[2]->Fill(mult3);
    hMult[3]->Fill(mult4);

    hTof1->Fill(tof1);
    hTof2->Fill(tof2);

    if (tof1 != 0 && tof2 != 0) {

      double time = tof2 - tof1;
      hTof->Fill(time);

      // //SiTil 1-3
      // time += (tot1-48.65)*tan(-6*TMath::Pi()/180.);
      // time += (tot2-36.83)*tan(5*TMath::Pi()/180.);

      // //SiTil 1-2
      // time += (tot1-48.65)*tan(-4*TMath::Pi()/180.);
      // time += (tot2-36.83)*tan(2*TMath::Pi()/180.);

      // MCP-TOF
      // time += (tot1-47.28)*tan(-2*TMath::Pi()/180.);
      // time += (tot2-49.59)*tan(-1*TMath::Pi()/180.);

      // time += (tot1-41.32)*tan(-4*TMath::Pi()/180.);
      // time += (tot2-40.75)*tan(2*TMath::Pi()/180.);

      if (gcFile != "0") {
        if (tof1tot < 55) time += (tot1 - tof1tot) * tan(-3 * TMath::DegToRad());
        if (tof2tot < 55) time += (tot2 - tof2tot) * tan(2 * TMath::DegToRad());
      }

      hTofC->Fill(time);

      // if(insideOfEllipce(time, tot1, tof1le, tof1tot, c1y, c1x) && insideOfEllipce(time, tot2,
      // tof1le, tof2tot, c1y, c1x)){
      hLeTotC->Fill(time, tot1);
      hLeTotC2->Fill(time, tot2);
      // }else if(insideOfEllipce(time, tot1, tof2le, tof1tot, c2y, c2x) && insideOfEllipce(time,
      // tot2, tof2le, tof2tot, c2y, c2x)){ 	hLeTotC->Fill(time,tot1); 	hLeTotC2->Fill(time,tot2);
      // }

      hLeTot1->Fill(tof2 - tof1, tot1);
      hLeTot2->Fill(tof2 - tof1, tot2);
    }
  }

  t.add_canvas("LeTot1", 800, 400);
  hLeTot1->Draw("colz");
  t.add_canvas("LeTot2", 800, 400);
  hLeTot2->Draw("colz");

  t.add_canvas("LeTotC", 800, 400);
  hLeTotC->Draw("colz");
  TEllipse *el1 = new TEllipse(tof1le, tof1tot, c1y, c1x);
  el1->SetLineColor(2);
  el1->SetLineWidth(2);
  el1->SetFillStyle(0);
  // el1->Draw();
  TEllipse *el2 = new TEllipse(tof2le, tof1tot, c2y, c2x);
  el2->SetLineColor(2);
  el2->SetLineWidth(2);
  el2->SetFillStyle(0);
  // el2->Draw();

  t.add_canvas("LeTotC2", 800, 400);
  hLeTotC2->Draw("colz");
  TEllipse *el3 = new TEllipse(tof1le, tof2tot, c1y, c1x);
  el3->SetLineColor(2);
  el3->SetLineWidth(2);
  el3->SetFillStyle(0);
  // el3->Draw();
  TEllipse *el4 = new TEllipse(tof2le, tof2tot, c2y, c2x);
  el4->SetLineColor(2);
  el4->SetLineWidth(2);
  el4->SetFillStyle(0);
  // el4->Draw();

  t.add_canvas("tof", 800, 400);
  t.fit(hTof, 10, 20, 2, 2);
  hTof->Draw();

  t.add_canvas("tofC", 800, 400);
  TVector3 r = t.fit(hTofC, 10, 20, 2, 2);
  hTofC->Draw();
  double tot1 = hLeTotC->GetMean(2);
  double tot2 = hLeTotC2->GetMean(2);
  std::cout << infile << "  LE1  " << r.X() << "   LE2 " << r.Z() << " TOT1 " << tot1 << "  TOT2 "
            << tot2 << std::endl;

  TFile efile(infile + ".tof.root", "RECREATE");
  TGraph *gr = new TGraph();
  gr->SetPoint(0, r.X(), r.Z());
  gr->SetPoint(1, tot1, tot2);
  gr->SetName("tof_" + t.run()->getName());
  gr->Write();
  efile.Write();
  efile.Close();

  t.add_canvas("mult", 800, 400);
  hMult[0]->Draw();
  for (int i = 1; i < 4; i++) {
    hMult[i]->Draw("same");
  }

  TLegend *leg = new TLegend(0.5, 0.6, 0.8, 0.85);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  for (int i = 0; i < 4; i++) leg->AddEntry(hMult[i], names[i], "lp");
  leg->Draw();

  gStyle->SetOptTitle(0);
  t.save_canvas("data/draw_tof",0);
}
