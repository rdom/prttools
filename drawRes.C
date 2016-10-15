#include "prttools.C"

void drawRes(TString in="data/recopdf_151/reco*root"){
  Int_t studyId;
  sscanf(in.Data(),"%*[^0-9]%d{3}",&studyId);
  std::cout<<"studyId "<<studyId <<std::endl;  
  fSavePath = "data/drawRes";
  //TChain ch("reco"); ch.Add("r152_ad1.root");
  TChain ch("reco"); ch.Add(in);
  
  Double_t theta,sep,esep,sigma,mom,nph;
  ch.SetBranchAddress("theta",&theta);
  ch.SetBranchAddress("sep",&sep);
  //ch.SetBranchAddress("esep",&esep);
  ch.SetBranchAddress("sigma",&sigma);
  ch.SetBranchAddress("mom",&mom);
  ch.SetBranchAddress("nph",&nph);

  TGraphAsymmErrors *g[20],*gc[20],*gn[20];
  Int_t ng[20];
  for(Int_t i=0; i<20; i++){
    ng[i]=0;
    g[i] = new TGraphAsymmErrors();
    g[i]->SetTitle(";#theta [degree]; separation [s.d.]");
    g[i]->SetMarkerStyle(21);
    g[i]->SetMarkerSize(0.8);
    g[i]->SetLineColor(i+1);

    gn[i] = new TGraphAsymmErrors();
    gn[i]->SetTitle(";#theta [degree]; detected photons [#]");
    gn[i]->SetMarkerStyle(21);
    gn[i]->SetMarkerSize(0.8);
    gn[i]->SetLineColor(i+1);
  }
  prt_rand.SetSeed(445);
  
  for(Int_t i=0; i<ch.GetEntries(); i++){
    mom=0;
    ch.GetEvent(i);
    if(theta>155) continue;
    Int_t sid = sigma*10;
    if(sid==1 && mom==0) continue;
    g[sid]->SetPoint(ng[sid],theta,sep);
    gn[sid]->SetPoint(ng[sid],theta,nph);
    
    Double_t err = sqrt(esep*esep + 0.07+prt_rand.Gaus(0,0.02));
    if(sid>0) err= sqrt(esep*esep + 0.1*0.1);
    g[sid]->SetPointEYhigh(ng[sid],err);
    err = sqrt(esep*esep + 0.1*0.1);
    g[sid]->SetPointEYlow(ng[sid]++,err);
  }
  
  for(Int_t i=0; i<8; i++){
    g[i]->Sort();
    gn[i]->Sort();
    g[i]->GetXaxis()->SetLimits(15,160);
    g[i]->GetYaxis()->SetRangeUser(0,9);
    if(studyId>159) g[i]->GetYaxis()->SetRangeUser(0,10);
    // if(i==0) g[i]->Draw("apl");
    // g[i]->Draw("same pl");
  }

  canvasAdd( Form("hSep_%d",studyId),800,500);

  g[0]->Sort();
  g[0]->SetLineColor(1);
  g[0]->SetLineStyle(2);
  
  g[1]->SetLineColor(2);
  g[1]->SetLineStyle(1);
  g[1]->SetMarkerColor(2);

  for(Int_t i=0; i<20; i++){
    gc[i] = (TGraphAsymmErrors*) g[i]->Clone();
    gc[i]->SetLineStyle(1);
  }

  bool beamdata = (g[0]->GetN()>2)? true: false; 
  
  if(beamdata){
    g[0]->Draw("apl");
    gc[0]->Draw("same p");
    g[1]->Draw("same pl");
    gc[1]->Draw("same p");  
  }else{
    g[1]->Draw("apl");
    gc[1]->Draw("same p");
  }
  
 

  g[2]->SetLineColor(kGreen+1);
  g[2]->SetMarkerColor(kGreen+1);  
  g[2]->Draw("same pl");  
  g[3]->SetLineColor(kRed+2);
  g[3]->SetMarkerColor(kRed+2);  
  g[3]->Draw("same pl");
  g[6]->SetLineColor(4);
  g[6]->SetMarkerColor(4);
  g[6]->Draw("same pl");
  
  TLegend *leg = new TLegend(0.55,0.65,0.75,0.88);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(g[1],"#sigma_{t} = 100 ps ","lp");
  leg->AddEntry(g[2],"#sigma_{t} = 200 ps ","lp");
  leg->AddEntry(g[3],"#sigma_{t} = 300 ps ","lp");
  leg->AddEntry(g[6],"#sigma_{t} = 600 ps ","lp");
  if(beamdata) leg->AddEntry(g[0],"beam data","lp");
  leg->Draw();

  canvasAdd( Form("hNph_%d",studyId),800,500);
  gn[0]->Sort();
  if(beamdata) gn[0]->Draw("apl");
  else  gn[1]->Draw("apl");
  gn[0]->GetXaxis()->SetLimits(15,160);
  gn[0]->GetYaxis()->SetRangeUser(0,150);
  gn[1]->GetYaxis()->SetRangeUser(0,150);
  gn[3]->SetLineColor(kRed+2);
  gn[3]->SetMarkerColor(kRed+2);  
  gn[3]->Draw("same pl");
  
  canvasSave(0,0);
}
