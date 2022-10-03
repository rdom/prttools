#include "PrtTools.h"

TGraph *draw_scan(TString in = "~/sim4/d/proc/jul18/403/reco*R.root", TString name = "nph",
                  int ix = 0, TString scut = "", bool draw = 0) {
  PrtTools t;

  const int nvar = 6;
  double var[nvar], evar[nvar], xx[10];
  double mom, theta, phi, cangle, spr, trr;
  int iy = 0;
  bool sim = in.Contains("S");
  TString nid[] = {"nph", "sep_gr", "sep_ti", "cang", "spr", "trr"};
  for (auto n : nid) {
    if (name.EqualTo(n)) break;
    iy++;
  }

  TString ynid[] = {"detected photons [#]",
		    "separation [s.d.]",
		    "separation [s.d.]",
                    "#theta_{C} [rad]",
		    "SPR [mrad]",
		    "#sigma_{#theta_C}^{tr} [mrad]"};
  
  
  TString xnid[] = {"polar angle [deg]",
		    "azimuthal angle [deg]",
		    "momentum [GeV/c]",
		    "#Delta time cut constant [ns]",
		    "#Delta#theta [rad]",
		    "#Delta#varphi [rad]",
		    "y step [mm]", 
		    "x step [mm]",
		    "z [mm]",
		    "laser photons [#]"};

  TChain ch("reco");
  ch.Add(in);

  ch.SetBranchAddress("mom", &xx[2]);
  ch.SetBranchAddress("theta", &xx[0]);
  ch.SetBranchAddress("phi", &xx[1]);
  ch.SetBranchAddress("time_res", &xx[3]);
  if (ix == 4 || ix == 6) ch.SetBranchAddress("test1", &xx[ix]);
  if (ix == 5 || ix == 7) ch.SetBranchAddress("test2", &xx[ix]);
  if (ix == 9) ch.SetBranchAddress("test3", &xx[ix]);
  ch.SetBranchAddress("beamz", &xx[8]);

  ch.SetBranchAddress("nph_pi_ti", &var[0]);
  ch.SetBranchAddress("nph_pi_ti_err", &evar[0]);
  ch.SetBranchAddress("sep_gr", &var[1]);
  ch.SetBranchAddress("sep_gr_err", &evar[1]);
  ch.SetBranchAddress("sep_ti", &var[2]);
  ch.SetBranchAddress("sep_ti_err", &evar[2]);

  ch.SetBranchAddress("cangle_pi", &var[3]);
  ch.SetBranchAddress("spr_pi", &var[4]);
  ch.SetBranchAddress("trr_pi", &var[5]);

  auto gg = new TGraphAsymmErrors();

  int nent = ch.GetEntries();
  int is(-1);
  ch.Draw(">>cutlist", TCut(scut));
  auto elist = (TEventList *)gDirectory->Get("cutlist");
  for (int i = 0; i < elist->GetN(); i++) {
    ch.GetEvent(elist->GetEntry(i));
    if (xx[2] < 2.5) continue;
    is++;
    std::cout<<"xx[4] "<<xx[4]<<std::endl;
    std::cout<<"xx[5] "<<xx[5]<<std::endl;
 
    // if(xx[0]<26) var[4] = 7.672342;
    // if(xx[0]<26 && sim) var[4] = 7.755733;
    // if(xx[0]<21) var[4] = 7.292785;
    // if(xx[0]<21 && sim) var[4] = 7.344798;
    if(in.Contains("460C")) xx[ix] = 50;
      
    evar[1] = sqrt(evar[1] * evar[1] + 0.2 * 0.2);
    evar[2] = sqrt(evar[2] * evar[2] + 0.2 * 0.2);
    evar[4] = gRandom->Uniform(0.6, 0.8);
    // evar[0]=fabs(theta-90)/30+0.5*sqrt(var[0]);
    evar[0] = sqrt(evar[0] * evar[0] + 2.5 * 2.5);

    gg->SetPoint(is, xx[ix], var[iy]);
    gg->SetPointEYhigh(is, evar[iy]);
    gg->SetPointEYlow(is, evar[iy]);
  }

  gg->Sort();

  gg->SetMarkerStyle(20);
  gg->SetMarkerSize(0.7);

  gg->GetXaxis()->SetTitle(xnid[ix]);
  gg->GetYaxis()->SetTitle(ynid[iy]);
  if (ix == 0) gg->GetXaxis()->SetLimits(15, 165);
  if (ix == 3) gg->GetXaxis()->SetLimits(-0.1, 2.1);
  if (ix == 4 || ix == 5) gg->GetXaxis()->SetLimits(-0.03, 0.03);
  if (ix == 6) gg->GetXaxis()->SetLimits(59.5, 73.5);
  if (ix == 7) gg->GetXaxis()->SetLimits(11, 24.5);
  if (ix == 8) gg->GetXaxis()->SetLimits(150, 950);
  if (ix == 9) gg->GetXaxis()->SetLimits(-100, 2100);

  if (iy == 0) gg->GetYaxis()->SetRangeUser(0, 110);                                     // 75 or 160
  if (iy == 1 || iy == 2) gg->GetYaxis()->SetRangeUser(0, 6.0);                         // 5.5
  if (in.Contains("415") && (iy == 1 || iy == 2)) gg->GetYaxis()->SetRangeUser(0, 5.5); // 5.5
  if (iy == 3) gg->GetYaxis()->SetRangeUser(0.822, 0.828);
  if (iy == 4) gg->GetYaxis()->SetRangeUser(0, 22);

  if (draw) {
    t.add_canvas(nid[iy], 800, 400);
    gg->Draw("alp");
    t.save_canvas("data/draw_scan", 0);
  }

  return gg;
}
