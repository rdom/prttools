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

  prt_canvasAdd(Form("hspr_%d",(int)theta),800,500);
  hSpr->Draw("colz");
  
  // int binx,biny,binz;  
  // hSpr->GetBinXYZ(hSpr->GetMinimumBin(),binx,biny,binz);
  // Double_t binxCenter = hSpr->GetXaxis()->GetBinCenter(binx);
  // Double_t binyCenter = hSpr->GetYaxis()->GetBinCenter(biny);
  // TEllipse *el = new TEllipse(binxCenter,binyCenter,.5,.5);
  // el->Draw();
  
  int n =15;
  double th[]  ={  20,   25,  30,  40,   50,   60,   70,   80,  90,  100,  110,  120,  130,  140,  150};
  double a[]   ={3.84,-2.38,3.98,1.60,-1.49,-2.60,-0.81,-2.98,0.56,-5.00, 4.59, 3.48, 7.25, 4.69, 1.88};
  double ph[]  ={4.45,6.11 ,5.36,6.03, 6.13, 3.31, 4.55, 3.14,0.00, -0.5,-5.46,-6.23,-5.29,-5.47,-6.03};

  TGraph *gadiff = new TGraph(n,th,a);
  TGraph *gphi = new TGraph(n,th,ph);
  gadiff->SetTitle(";#theta [degree]; #theta - #theta_{fit} [degree]");
  gadiff->SetMarkerStyle(20);
  gadiff->SetMarkerSize(0.8);
  gphi->SetTitle(";#theta [degree]; #varphi_{fit} [degree]");
  gphi->SetMarkerStyle(20);
  gphi->SetMarkerSize(0.8);

  prt_canvasAdd("adiff",800,400);
  gadiff->Draw("APL");
  prt_canvasAdd("phi",800,400);
  gphi->Draw("APL");
  
  prt_canvasSave();
  
}
