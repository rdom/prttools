#define prt__beam
#include "../prtdirc/src/PrtHit.h"
#include "../prtdirc/src/PrtEvent.h"
#include "prttools.C"
#include <TVirtualFitter.h>
#include <TKey.h>
#include <TRandom.h>

Int_t mcpdata[15][65];
Int_t cluster[15][65];
Int_t lneighbours[65];
Int_t lsize(0);

Int_t getneighbours(Int_t m, Int_t p){
  for(Int_t i=0; i<65; i++) if(p==lneighbours[i]) return -1;
  lneighbours[lsize]=p;
  lsize++;
  for(Int_t t=0; t<65; t++){
    if(mcpdata[m][t]){
      for(Int_t i=0; i<65; i++) if(t==lneighbours[i]) continue;
      if((t==p-1 && p%8!=0) || (t==p+1 && p%8!=7) ||
	 (t==p+8 && p<57) || (t==p-8 && p>8)) getneighbours(m,t);
    }
  }
  return lsize;
}

void getclusters(){
  for(Int_t m=0; m<15; m++){
    for(Int_t p=0; p<65; p++){
      if(mcpdata[m][p])  cluster[m][p] = getneighbours(m,p);
      lsize=0;
      for(Int_t i=0; i<65; i++) lneighbours[i]=0;
    }
  }
}



void recoPdf(TString path="$HOME/proc/152/beam_15183021251SF.root", TString pdf="$HOME/proc/152/beam_15183021251SF.root", Double_t sigma=1){
//void recoPdf(TString path="$HOME/proc/152/beam_15183013641SF.root", TString pdf="$HOME/proc/152/beam_15183013641SF.root", Double_t sigma=1){
  if(path=="") return;
  Int_t studyId;
  sscanf(path, "%*[^0-9]%d{3}",&studyId);
  fSavePath = Form("data/recopdf_%d",studyId);
  PrtInit(path,1);
  gStyle->SetOptStat(0);
  CreateMap();
  TGaxis::SetMaxDigits(4);
  
  //TCanvas *cc = new TCanvas("cc","cc");
  TH1F *hllf= new TH1F("hllf","hllf;ln L(p) - ln L(#pi); entries [#]",200,-50,50);
  TH1F *hlls= new TH1F("hlls","hlls;ln L(p) - ln L(#pi); entries [#]",200,-50,50);


  TH1F *hl1 = new TH1F("hl1","pdf;LE time [ns]; entries [#]", 1000,0,100);
  TH1F *hl2 = new TH1F("hl2","pdf;LE time [ns]; entries [#]", 1000,0,100);
  TH1F *hl3 = new TH1F("hl3","pdf;LE time [ns]; entries [#]", 1000,0,100);

  TRandom rand;
  TF1 *pdff[960],*pdfs[960];
  TH1F *hpdff[960],*hpdfs[960];
  if(path.Contains("F.root")) pdf.ReplaceAll("F.root","P.pdf.root");
  
  else  pdf.ReplaceAll(".root",".pdf.root");
  TFile f(pdf);
  
  hl3->Rebin((Int_t)sigma);
  Int_t integ1(0), integ2(0);
  for(Int_t i=0; i<960; i++){
    hpdff[i] = (TH1F*)f.Get(Form("hf_%d",i));
    hpdfs[i] = (TH1F*)f.Get(Form("hs_%d",i));
    hpdff[i]->Rebin((Int_t)sigma);
    hpdfs[i]->Rebin((Int_t)sigma);
    integ1+= hpdff[i]->Integral();
    integ2+= hpdfs[i]->Integral();
    hl3->Add(hpdff[i]);
    hl3->Add(hpdfs[i]);
  }
  if(path.Contains("C.root")) sigma=0;
  
  Double_t theta(0);
  TVirtualFitter *fitter;
  Double_t time,timeres(-1);
  PrtHit fHit;
  Int_t totalf(0),totals(0), ch, entries = 10000; //fCh->GetEntries();
  for (Int_t ievent=0; ievent<entries; ievent++){
    PrtNextEvent(ievent,1000);
    if(ievent==0){
      theta = prt_event->GetAngle();
    }
    timeres = prt_event->GetTimeRes();
    Double_t aminf,amins, sum(0),sumf(0),sums(0);

    Int_t nHits =prt_event->GetHitSize();

    
    // //clusters search
    // for(Int_t h=0; h<nHits; h++) {
    //   Int_t mid=prt_event->GetHit(h).GetMcpId();
    //   Int_t pid=prt_event->GetHit(h).GetPixelId()-1;
    //   mcpdata[mid][pid]=1;
    // }
    // getclusters();
    
    for(Int_t i=0; i<nHits; i++){
      fHit = prt_event->GetHit(i);
      ch=map_mpc[fHit.GetMcpId()][fHit.GetPixelId()-1];      
      time = fHit.GetLeadTime() + rand.Gaus(0,sigma/10.);
      
      Int_t mid=prt_event->GetHit(i).GetMcpId();
      Int_t pid=prt_event->GetHit(i).GetPixelId()-1;
      // if(cluster[mid][pid]>6) {
      // 	std::cout<<"cluster[mid][pid]  "<< cluster[mid][pid] <<std::endl;	
      // 	continue;
      // }
      
      // if(prt_event->GetParticle()==211) time += 0.2; //fix offset
      // if(prt_event->GetParticle()==2212) time -= 0.35;
      // if(time<10 || time >30) continue;
      // if(prt_event->GetParticle()==211) time -= 0.05; //fix offset
      // if(prt_event->GetParticle()==2212) time -= 0.05;

      
      if(time<5 || time>150) continue;
      aminf = hpdff[ch]->GetBinContent(hpdff[ch]->FindBin(time)); 
      amins = hpdfs[ch]->GetBinContent(hpdfs[ch]->FindBin(time));

      //      if(aminf==0 || amins==0) continue;
      Double_t noise = 1e-3; //1e-7;
      sumf+=TMath::Log((aminf+noise));
      sums+=TMath::Log((amins+noise));    


      // std::cout<< aminf<< " "<<amins <<std::endl;
      // //  if(aminf>amins)
      // //	if(aminf==0 || amins==0)
      // 	{
      // 	  prt_normalize(hpdff[ch],hpdfs[ch]);
      // 	hpdff[ch]->SetLineColor(2);
      // 	hpdff[ch]->Draw();
      // 	hpdfs[ch]->SetLineColor(4);
      // 	hpdfs[ch]->Draw("same");
      // 	TLine *line = new TLine(time,0,time,1000);
      // 	line->Draw("same");
      // 	cc->Update();
      // 	cc->WaitPrimitive();
      // }


      if(prt_event->GetParticle()==2212) hl1->Fill(time);
      if(prt_event->GetParticle()==211 || prt_event->GetParticle()==212) hl2->Fill(time);

    }
    if(fabs(sumf-sums)<0.3) continue;
    sum = sumf-sums;
    
    //std::cout<<"sum ===========  "<<sum  << "  "<< sumf<< "  "<< sums<<std::endl;
    
    if(prt_event->GetParticle()==2212) hllf->Fill(sum);
    if(prt_event->GetParticle()==211 || prt_event->GetParticle()==212)  hlls->Fill(sum);

    for(Int_t j=0; j<15; j++){
      for(Int_t i=0; i<65; i++){
	mcpdata[j][i]=0;
	cluster[j][i]=0;
      }
    }
     
  }

  TString name = Form("tis_%d_%1.1f_%1.1f.root",studyId,theta,sigma);
  if(path.Contains("C.root")) name = Form("tid_%d_%1.1f_%1.1f.root",studyId,theta,sigma);
  canvasAdd("ll_"+name,800,400);

  prt_normalize(hllf,hlls);
  
  TF1 *ff;
  Double_t m1,m2,s1,s2; 
  if(hllf->GetEntries()>10){
    hllf->Fit("gaus","Sq");
    ff = hllf->GetFunction("gaus");
    m1=ff->GetParameter(1);
    s1=ff->GetParameter(2);
  }
  if(hlls->GetEntries()>10){
    hlls->Fit("gaus","Sq");
    ff = hlls->GetFunction("gaus");
    m2=ff->GetParameter(1);
    s2=ff->GetParameter(2);
  }
  Double_t sep = (fabs(m1-m2))/(0.5*(s1+s2));
  std::cout<<path<<" separation "<< sep <<std::endl;
  hllf->SetTitle(Form("separation = %1.2f",sep));

  
  hllf->SetLineColor(2);
  hllf->Draw();
  hlls->SetLineColor(4);
  hlls->Draw("same");

  hl1->Scale(1/hl1->GetMaximum());
  hl2->Scale(1/hl2->GetMaximum());
  hl3->Scale(1/hl3->GetMaximum());

  prt_normalize(hl1,hl2);
  canvasAdd("hl_"+name,800,500);
  hl1->Draw();
  hl2->SetLineColor(4);
  hl2->Draw("same");
  hl3->SetLineColor(2);
  hl3->Draw("same");
  canvasSave(1,0);


  TFile fc("reco_"+name,"recreate");
  TTree *tc = new TTree("reco","reco");
  tc->Branch("theta",&theta,"theta/D");
  tc->Branch("sep",&sep,"sep/D");
  tc->Branch("sigma",&sigma,"sigma/D");
  sigma/=10.;
  tc->Fill();
  tc->Write();
}
