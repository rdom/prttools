#define prt__beam
#include "../prtdirc/src/PrtHit.h"
#include "../prtdirc/src/PrtEvent.h"
#include "prttools.C"
#include <TVirtualFitter.h>
#include <TKey.h>

void recoPdf(TString path="/data.local/data/jun15/beam_15177050804S.root"){
  path.ReplaceAll(".root","");
  fSavePath = "data/recopdf_151";
  PrtInit(path+".root",1);
  gStyle->SetOptStat(0);
  CreateMap();
  TGaxis::SetMaxDigits(4);
  
  TCanvas *cc = new TCanvas("cc","cc");
  TH1F *hllf= new TH1F("hllf","hllf;ln L(p) - ln L(#pi); entries [#]",100,-100,100);
  TH1F *hlls= new TH1F("hlls","hlls;ln L(p) - ln L(#pi); entries [#]",100,-100,100);


  TH1F *hl1 = new TH1F("hl1","pdf;LE time [ns]; entries [#]", 500,0,50);
  TH1F *hl2 = new TH1F("hl2","pdf;LE time [ns]; entries [#]", 500,0,50);
  TH1F *hl3 = new TH1F("hl3","pdf;LE time [ns]; entries [#]", 500,0,50);
  
  TF1 *pdff[960],*pdfs[960];
  TH1F *hpdff[960],*hpdfs[960];
  //  TFile f(path+".pdf.root");
  // TFile f("/data.local/data/jun15/scan_tr/hh1521_04.pdf1.root");
  //TFile f("/data.local/data/jun15/scan_tr/hh1522_03M.pdf1.root");
  // TFile f("/data.local/data/jun15/hh1613_sa.pdf.root"); //plate

  //TFile f("/data.local/data/jun15/hh151_37_s.pdf.root");
  TFile f("/data.local/data/jun15/beam_15177135523C.pdf.root");
  
    Int_t integ1(0), integ2(0);
  for(Int_t i=0; i<960; i++){
    hpdff[i] = (TH1F*)f.Get(Form("hf_%d",i));
    hpdfs[i] = (TH1F*)f.Get(Form("hs_%d",i));
    integ1+= hpdff[i]->Integral();
    integ2+= hpdfs[i]->Integral();
    hl3->Add(hpdff[i]);
    hl3->Add(hpdfs[i]);
  }

  TVirtualFitter *fitter;
  Double_t time,timeres(-1);
  PrtHit fHit;
  Int_t totalf(0),totals(0), ch, entries = fCh->GetEntries();
  for (Int_t ievent=0; ievent<entries; ievent++){
    PrtNextEvent(ievent,1000);
    timeres = prt_event->GetTimeRes();
    Double_t aminf,amins, sum(0),sumf(0),sums(0);
    //  if(prt_event->GetHitSize()<5) continue; 
    for(Int_t i=0; i<prt_event->GetHitSize(); i++){
      fHit = prt_event->GetHit(i);
      ch=map_mpc[fHit.GetMcpId()][fHit.GetPixelId()-1];      
      time = fHit.GetLeadTime();

      // if(prt_event->GetParticle()==211) time += 0.2; //fix offset
      // if(prt_event->GetParticle()==2212) time -= 0.35;
      // if(time<10 || time >30) continue;

      if(prt_event->GetParticle()==211) time -= 0.05; //fix offset
      if(prt_event->GetParticle()==2212) time -= 0.05;
      //if(time>7) continue;
      
      
      //std::cout<<ch<<" "<< hpdff[ch]->FindBin(time)<<std::endl;
      aminf = hpdff[ch]->GetBinContent(hpdff[ch]->FindBin(time)); 
      amins = hpdfs[ch]->GetBinContent(hpdfs[ch]->FindBin(time));

      //      if(aminf==0 || amins==0) continue;
      Double_t noise = 1e-5; //1e-7;
      sumf+=-TMath::Log((aminf+noise));
      sums+=-TMath::Log((amins+noise));    


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
      if(prt_event->GetParticle()==211) hl2->Fill(time);

    }
    if(sumf==sums) continue;
    sum = sumf-sums;
    
    //std::cout<<"sum ===========  "<<sum  << "  "<< sumf<< "  "<< sums<<std::endl;
    
    if(prt_event->GetParticle()==2212) hllf->Fill(sum);
    if(prt_event->GetParticle()==211 || prt_event->GetParticle()==212)  hlls->Fill(sum);
    
  }
  
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
  canvasAdd("ll",800,400);
  hllf->SetLineColor(2);
  hllf->Draw();
  hlls->SetLineColor(4);
  hlls->Draw("same");

  hl1->Scale(1/hl1->GetMaximum());
  hl2->Scale(1/hl2->GetMaximum());
  hl3->Scale(1/hl3->GetMaximum());

  prt_normalize(hl1,hl2);
  canvasAdd("hl",800,500);
  hl1->Draw();
  hl2->SetLineColor(4);
  hl2->Draw("same");
  hl3->SetLineColor(2);
  hl3->Draw("same");
  canvasSave(1,0);
}
