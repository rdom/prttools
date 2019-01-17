#include "prttools.C"
#include <algorithm>
#include <iterator>
#include <TGraphAsymmErrors.h>


void drawRes(TString in="data/recopdf_151/reco*root",TString in2=""){

  TChain ch("reco"); ch.Add(in);
  int studyId=prt_get3digit(in);
  std::cout<<"studyId "<<studyId <<std::endl;  
  prt_savepath = Form("data/drawRes_%d",studyId);
  
  Double_t sep,esep,sigma,mom,nph,enph;
  int theta;
  ch.SetBranchAddress("theta",&theta);
  ch.SetBranchAddress("sep",&sep);
  ch.SetBranchAddress("sigma",&sigma);
  ch.SetBranchAddress("mom",&mom);
  ch.SetBranchAddress("nph",&nph);
  ch.SetBranchAddress("enph",&enph);

  TString names[500]={""};
  names[403]="#varphi = 0 deg";
  names[406]="#varphi = 5 deg";
  names[407]="#varphi = 10 deg";
  names[408]="#varphi = 15 deg";
  names[409]="#varphi = 2.5 deg";

  int coll[]={kBlack,kRed+1,kGreen,kRed,4,kCyan-6,6,7,8,9,10};
  int colm[]={kBlack,kRed+1,kGreen+2,kRed+1,4,kCyan-6,6,7,8,9,10};

  TGraphAsymmErrors *g[20],*gc[20],*gn[20];
  int ng[20];
  for(int i=0; i<20; i++){
    ng[i]=0;
    g[i] = new TGraphAsymmErrors();
    g[i]->SetTitle(";polar angle [deg]; separation [s.d.]");
    g[i]->SetName("gr");
    g[i]->SetMarkerStyle(20);
    g[i]->SetMarkerSize(0.8);
    g[i]->SetLineColor(coll[i]);
    g[i]->SetMarkerColor(colm[i]);
    g[i]->SetLineStyle((i==0)? 1:2);
    
    gn[i] = new TGraphAsymmErrors();
    gn[i]->SetTitle(";polar angle [deg]; detected photons [#]");
    gn[i]->SetName("gr");
    gn[i]->SetMarkerStyle(20);
    gn[i]->SetMarkerSize(0.8);
    gn[i]->SetLineColor(coll[i]);
    gn[i]->SetMarkerColor(colm[i]);
    gn[i]->SetLineStyle((i==0)? 1:2);
  }
  
  prt_rand.SetSeed(445);
  std::vector<int> vec;

  for(int i=0; i<ch.GetEntries(); i++){
    ch.GetEvent(i);
    if (std::find(vec.begin(), vec.end(), sigma) == vec.end()) vec.push_back((int)sigma);
  }
  std::sort(vec.begin(), vec.end());
  int size=vec.size();
  bool beamdata = (vec[0]<0.1)? true: false;
  
  for(int i=0; i<ch.GetEntries(); i++){    
    ch.GetEvent(i);
    if(theta>140) continue;
    int sid = std::distance(vec.begin(),std::find(vec.begin(), vec.end(),sigma));
    
    if(sid==1 && mom==0) continue;
    g[sid]->SetPoint(ng[sid],theta,sep);
    gn[sid]->SetPoint(ng[sid],theta,nph);
    
    Double_t err= sqrt(esep*esep + 0.1*0.1);
    if(beamdata && sid==0) err = sqrt(esep*esep + 0.07+prt_rand.Gaus(0,0.02));    
    g[sid]->SetPointEYhigh(ng[sid], sqrt(esep*esep + 0.2*0.2));
    g[sid]->SetPointEYlow(ng[sid], sqrt(esep*esep + 0.2*0.2));
    
    gn[sid]->SetPointEYhigh(ng[sid], sqrt(nph));
    gn[sid]->SetPointEYlow(ng[sid], sqrt(nph));
    ng[sid]++;		  
    std::cout<<sid<<" " << mom<<" "<<theta<<" "<<sep <<" "<<sigma<<std::endl;    
  }
  
  for(int i=0; i<size; i++){
    g[i]->Sort();
    gn[i]->Sort();
    g[i]->GetXaxis()->SetLimits(15,145);
    g[i]->GetYaxis()->SetRangeUser(0,6);
    if(studyId>159) g[i]->GetYaxis()->SetRangeUser(0,6);
    // if(i==0) g[i]->Draw("apl");
    // g[i]->Draw("same pl");
  }

  prt_canvasAdd( Form("hSep_%d",studyId),800,500);

  for(int i=0; i<20; i++){
    gc[i] = (TGraphAsymmErrors*) g[i]->Clone();
    gc[i]->SetLineStyle(1);
  }
  
  // if(beamdata){
  //   g[0]->Draw("apl");
  //   gc[0]->Draw("same p");
  //   g[1]->Draw("same pl");
  //   gc[1]->Draw("same p");  
  // }else{
  //   g[1]->Draw("apl");
  //   gc[1]->Draw("same p");
  // }
     
  TLegend *leg = new TLegend(0.50,0.65,0.88,0.88);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);

  std::cout<<"size "<<size<<std::endl;
  
  for(int i=0; i<3; i++){
    if(i==1)continue;
    g[i]->Draw((i==0)? "apl":"same pl");

    if(beamdata && i==0) leg->AddEntry(g[0],"beam data "+names[studyId],"lp");
    else leg->AddEntry(g[i],Form("geant #sigma_{t} = %d ps ",vec[i]),"lp");
  }
  leg->Draw();

  prt_canvasAdd( Form("hNph_%d",studyId),800,500);
  for(int i=0; i<2; i++){
    gn[i]->GetXaxis()->SetLimits(15,145);
    gn[i]->GetYaxis()->SetRangeUser(0,140);
    gn[i]->Draw((i==0)? "apl":"same pl");
  }

  TLegend *leg1 = new TLegend(0.50,0.65,0.88,0.88);
  leg1->SetFillColor(0);
  leg1->SetFillStyle(0);
  leg1->SetBorderSize(0);
  leg1->SetFillStyle(0);
  leg1->AddEntry(gn[0],"beam data "+names[studyId],"lp");
  leg1->AddEntry(gn[1],"geant sim         ","lp");
  leg1->Draw();

  prt_canvasSave(0,0);
}
