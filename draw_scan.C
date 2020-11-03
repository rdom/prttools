#include "prttools.C"

TGraph* draw_scan(TString in = "~/sim4/d/proc/jul18/403/reco*R.root", int iy=3, int ix=0, TString scut="",bool draw=0){

  const int nvar=5;
  double var[nvar], evar[nvar],xx[10];     
  double mom,theta,phi,cangle,spr,trr;

  TString nid[] = {"nph","sep","spr","cang","trr"};  

  TString ynid[] = {"detected photons [#]",
		    "separation [s.d.]",
		    "#theta_{C} [rad]",
		    "SPR [mrad]",
		    "#sigma_{#theta_C}^{tr} [mrad]"};
  
  TString xnid[] = {"polar angle [deg]",
		    "azimuthal angle [deg]",
		    "momentum [GeV/c]",
		    // "time presision [ns]",
		    "#Delta time cut constant [ns]",
		    "y step [mm]",
		    //		    "x step [mm]",
		    "#varphi [deg]",
		    "z [mm]"};
  
  TChain ch("reco"); ch.Add(in);
  
  ch.SetBranchAddress("mom",&xx[2]);
  ch.SetBranchAddress("theta",&xx[0]);
  ch.SetBranchAddress("phi",&xx[1]);
  ch.SetBranchAddress("time_res",&xx[3]);
  ch.SetBranchAddress("test1",&xx[4]);
  ch.SetBranchAddress("test2",&xx[5]);
  ch.SetBranchAddress("beamz",&xx[6]);
  
  ch.SetBranchAddress("nph_pi",&var[0]);
  ch.SetBranchAddress("nph_pi_err",&evar[0]);
  ch.SetBranchAddress("sep",&var[1]);
  ch.SetBranchAddress("sep_err",&evar[1]);

  ch.SetBranchAddress("cangle_pi",&var[2]);
  ch.SetBranchAddress("spr_pi",&var[3]);
  ch.SetBranchAddress("trr_pi",&var[4]);

  auto gg = new TGraphAsymmErrors();
    
  int nent = ch.GetEntries();
  int is(-1);
  ch.Draw(">>cutlist",TCut(scut));
  auto elist = (TEventList*)gDirectory->Get("cutlist"); 
  for (int i = 0; i < elist->GetN(); i++) {
    ch.GetEvent(elist->GetEntry(i));
    is++;

    evar[1]=sqrt(evar[1]*evar[1]+0.1*0.1);
    evar[4]=prt_rand.Uniform(0.8,1.2);
    //evar[0]=fabs(theta-90)/30+0.5*sqrt(var[0]);
    evar[0]=sqrt(evar[0]*evar[0]+0.8*0.8);

    std::cout<<"var[2] "<<var[2]<<std::endl;
    
    gg->SetPoint(is,xx[ix],var[iy]);
    gg->SetPointEYhigh(is,evar[iy]);
    gg->SetPointEYlow(is,evar[iy]);
  }

  gg->Sort();

  gg->SetMarkerStyle(20);
  gg->SetMarkerSize(0.7);
  
  gg->GetXaxis()->SetTitle(xnid[ix]);
  gg->GetYaxis()->SetTitle(ynid[iy]);
  if(ix==0) gg->GetXaxis()->SetLimits(15,145);
  if(ix==3) gg->GetXaxis()->SetLimits(-0.1,2.1);
  if(ix==4) gg->GetXaxis()->SetLimits(13.5,22);
  //if(ix==5) gg->GetXaxis()->SetLimits(61.5,75.5);
  if(ix==5) gg->GetXaxis()->SetLimits(-0.1,2.1);
  
  if(ix==6) gg->GetXaxis()->SetLimits(150,950);
  if(iy==0) gg->GetYaxis()->SetRangeUser(6,14);
  if(iy==1) gg->GetYaxis()->SetRangeUser(0,3.5);
  if(iy==2) gg->GetYaxis()->SetRangeUser(0.822,0.828);
  if(iy==3) gg->GetYaxis()->SetRangeUser(0,12);

  if(draw){
    prt_canvasAdd(nid[iy],800,500);
    gg->Draw("alp");
    prt_canvasSave("data/draw_scan",0);
  }
  
  return gg;
}
