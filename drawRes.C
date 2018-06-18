#include "prttools.C"
#include <algorithm>
#include <iterator>
#include <TGraphAsymmErrors.h>


void drawRes(TString in="data/recopdf_151/reco*root"){
  int studyId;
  sscanf(in.Data(),"%*[^0-9]%d{3}",&studyId);
  std::cout<<"studyId "<<studyId <<std::endl;  
  prt_savepath = "data/drawRes";
  //TChain ch("reco"); ch.Add("r152_ad1.root");
  TChain ch("reco"); ch.Add(in);
  
  Double_t sep,esep,sigma,mom,nph,enph;
  int theta;
  ch.SetBranchAddress("theta",&theta);
  ch.SetBranchAddress("sep",&sep);
  //ch.SetBranchAddress("esep",&esep);
  ch.SetBranchAddress("sigma",&sigma);
  ch.SetBranchAddress("mom",&mom);
  ch.SetBranchAddress("nph",&nph);
  ch.SetBranchAddress("enph",&enph);

  
  TGraphAsymmErrors *g[20],*gc[20],*gn[20];
  int ng[20];
  for(int i=0; i<20; i++){
    int col= (i<5)? i+1: i+5;
    ng[i]=0;
    g[i] = new TGraphAsymmErrors();
    g[i]->SetTitle(";polar angle [deg]; separation [s.d.]");
    g[i]->SetName("gr");
    g[i]->SetMarkerStyle(20);
    g[i]->SetMarkerSize(0.8);
    g[i]->SetLineColor(col);

    gn[i] = new TGraphAsymmErrors();
    gn[i]->SetTitle(";polar angle [deg]; detected photons [#]");
    gn[i]->SetName("gr");
    gn[i]->SetMarkerStyle(20);
    gn[i]->SetMarkerSize(0.8);
    gn[i]->SetLineColor(col);
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
    //    if(theta>155 || theta<20) continue;
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

  if(beamdata){
    g[0]->SetLineColor(1);
    g[0]->SetLineStyle(2);
  }

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
  
  int coll[]={kBlack,kGreen+1,kRed,kRed,4,kCyan-6,6,7,8,9,10};
  int colm[]={kBlack,kGreen+2,kRed+1,kRed+1,4,kCyan-6,6,7,8,9,10};
   
  TLegend *leg = new TLegend(0.55,0.65,0.75,0.88);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);

  std::cout<<"size "<<size<<std::endl;
  
  for(int i=0; i<size; i++){
    g[i]->SetLineColor(coll[i]);
    g[i]->SetMarkerColor(colm[i]);
    //if(vec[i] != 300 && vec[i] != 0 ) continue;
    if(i==0) g[i]->Draw("apl");
    else g[i]->Draw("same pl");

    if(beamdata && i==0) leg->AddEntry(g[0],"beam data","lp");
    else leg->AddEntry(g[i],Form("#sigma_{t} = %d ps ",vec[i]),"lp");
    //else leg->AddEntry(g[i],"geant sim","lp");
  }
  leg->Draw();

  prt_canvasAdd( Form("hNph_%d",studyId),800,500);
  gn[0]->Sort();
  if(beamdata) gn[0]->Draw("apl");
  else  gn[0]->Draw("apl");
  gn[0]->GetXaxis()->SetLimits(15,145);
  gn[0]->GetYaxis()->SetRangeUser(0,140);
  gn[1]->GetYaxis()->SetRangeUser(0,140);
  gn[1]->SetLineColor(kRed);
  gn[1]->SetMarkerColor(kRed+1);  
  gn[1]->Draw("same pl");

  TLegend *leg1 = new TLegend(0.55,0.65,0.80,0.85);
  leg1->SetFillColor(0);
  leg1->SetFillStyle(0);
  leg1->SetBorderSize(0);
  leg1->SetFillStyle(0);
  leg1->AddEntry(gn[0],"beam data","lp");
  leg1->AddEntry(gn[1],"geant sim","lp");
  leg1->Draw();

  prt_canvasSave(0,0);
}
