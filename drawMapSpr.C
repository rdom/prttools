#include "prttools.C"

void drawMapSpr(TString in="tdata/rt.root"){
  prt_savepath = "data/drawMapSpr";

  TChain ch("dirc"); ch.Add(in);
  double cangle,spr,trr,nph,par1,par2,par3,par4,par5,par6,test1,test2,theta,phi; 
  
  ch.SetBranchAddress("spr",&spr);
  ch.SetBranchAddress("trr",&trr);
  ch.SetBranchAddress("nph",&nph);
  ch.SetBranchAddress("cangle",&cangle);
  ch.SetBranchAddress("par5",&par5);
  ch.SetBranchAddress("par6",&par6);
  ch.SetBranchAddress("test1",&test1);
  ch.SetBranchAddress("test2",&test2);
  ch.SetBranchAddress("theta",&theta);
  ch.SetBranchAddress("phi",&phi);


  gStyle->SetOptStat(0);
  TH2F *hSpr = new TH2F("hSpr",";#Delta#theta [mrad];#Delta#varphi [mrad]",40,-20,20,40,-10,10);

  for (auto i = 0; i < ch.GetEntries(); i++) {
    ch.GetEvent(i);
    hSpr->Fill(test1*1000,test2*1000,spr);
  }

  prt_canvasAdd("hspr",800,500);
  hSpr->Draw("colz");
  
}
