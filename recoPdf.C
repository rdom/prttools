#define prt__beam
#include "../prtdirc/src/PrtHit.h"
#include "../prtdirc/src/PrtEvent.h"
#include "prttools.C"
#include <TVirtualFitter.h>
#include <TKey.h>

void recoPdf(TString path="/data.local/data/jun15/beam_15177050804S"){

  //path="/data.local/data/jun15/beam_15183032354S"
  //  path="/data.local/data/jun15/hits_151_7";
  path="/data.local/data/jun15/hits_152_1s";
  //path="/data.local/data/jun15/beam_15183022858C";
  
  // path="/data.local/data/jun15/hits_153_1s_0ps";
  // path="/data.local/data/jun15/beam_15185011524C";
  
  fSavePath = "data/pdf";
  PrtInit(path+".root",1);
  gStyle->SetOptStat(0);
  CreateMap();

  TCanvas *cc = new TCanvas("cc","cc");
  TH1F *hllf= new TH1F("hllf","hllf",200,-100,100);
  TH1F *hlls= new TH1F("hlls","hlls",200,-100,100);
  
  TF1 *pdff[960],*pdfs[960];
  TH1F *hpdff[960],*hpdfs[960];
  //TFile f(path+".pdf.root");
  TFile f("/data.local/data/jun15/hits_152_1p.pdf.root");
  // TFile f("/data.local/data/jun15/hits_153_1_0ps.pdf.root");
  for(Int_t i=0; i<960; i++){
    hpdff[i] = (TH1F*)f.Get(Form("hf_%d",i));
    hpdfs[i] = (TH1F*)f.Get(Form("hs_%d",i));  
  }

  TVirtualFitter *fitter;
  Double_t time;
  PrtHit fHit;
  Int_t totalf(0),totals(0), ch, entries = fCh->GetEntries();
  for (Int_t ievent=0; ievent<entries; ievent++){
    PrtNextEvent(ievent,1000);

    Double_t aminf,amins, sum(0),sumf(0),sums(0);
    //  if(prt_event->GetHitSize()<5) continue; 
    for(Int_t i=0; i<prt_event->GetHitSize(); i++){
      fHit = prt_event->GetHit(i);
      ch=map_mpc[fHit.GetMcpId()][fHit.GetPixelId()-1];      
      time = fHit.GetLeadTime();
      
      //std::cout<<ch<<" "<< hpdff[ch]->FindBin(time)<<std::endl;
      aminf = hpdff[ch]->GetBinContent(hpdff[ch]->FindBin(time)); 
      amins = hpdfs[ch]->GetBinContent(hpdfs[ch]->FindBin(time));

      if(aminf>0) sumf+=TMath::Log(aminf);
      if(amins>0) sums+=TMath::Log(amins);

      // std::cout<< aminf<< " "<<amins <<std::endl;
      // // if(aminf>amins)
      // 	{
      // 	prt_normalize(hpdff[ch],hpdfs[ch]);
      // 	hpdff[ch]->SetLineColor(2);
      // 	hpdff[ch]->Draw();
      // 	hpdfs[ch]->SetLineColor(4);
      // 	hpdfs[ch]->Draw("same");
      // 	TLine *line = new TLine(time,0,time,1000);
      // 	line->Draw("same");
      // 	cc->Update();
      // 	cc->WaitPrimitive();
      // }

      
      // sumf+=aminf;
      // sums+=amins;
    }
    sum = sumf-sums;
    
    std::cout<<"sum ===========  "<<sum  << "  "<< sumf<< "  "<< sums<<std::endl;
    
    if(prt_event->GetParticle()==2212) hllf->Fill(sum);
    if(prt_event->GetParticle()==211 || prt_event->GetParticle()==212)  hlls->Fill(sum);
    
  }
  
  prt_normalize(hllf,hlls);
  canvasAdd("ll",800,500);
  hllf->SetLineColor(4);
  hllf->Draw();
  hlls->SetLineColor(2);
  hlls->Draw("same");

  canvasSave(1,0);



  TF1 *ff;
  Double_t m1,m2,s1,s2; 
  if(hllf->GetEntries()>10){
    hllf->Fit("gaus","S");
    ff = hllf->GetFunction("gaus");
    m1=ff->GetParameter(1);
    s1=ff->GetParameter(2);
  }
  if(hlls->GetEntries()>10){
    hlls->Fit("gaus","S");
    ff = hlls->GetFunction("gaus");
    m2=ff->GetParameter(1);
    s2=ff->GetParameter(2);
  }
  Double_t sep = (fabs(m1)+fabs(m2))/(0.5*(s1+s2));
  std::cout<<"separation "<< sep <<std::endl;
    
}
