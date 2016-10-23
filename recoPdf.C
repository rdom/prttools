#define prt__beam
#include "../prtdirc/src/PrtHit.h"
#include "../prtdirc/src/PrtEvent.h"
#include "prttools.C"
#include <TVirtualFitter.h>
#include <TKey.h>
#include <TRandom.h>

TLine *gLine = new TLine(0,0,0,1000);

void recoPdf(TString path="$HOME/simo/build/beam_15184203911SF.root", TString pdfEnding=".h2Zpdf.root", Double_t sigma=3,Bool_t debug=false){ 
  
  if(path=="") return;
  Int_t studyId;
  sscanf(path, "%*[^/]/%d{3}",&studyId);
  
  fSavePath = Form("data/recopdf_%d",studyId);
  PrtInit(path,1);
  // gStyle->SetOptStat(0);
  CreateMap();
  TGaxis::SetMaxDigits(4);
  
  TCanvas *cc = new TCanvas("cc","cc");

  TH1F *hpdff[maxch_dirc],*hpdfs[maxch_dirc], *hl[5],*hnph[5],*hll[5];
  for(Int_t i=0; i<5; i++){
    hl[i] = new TH1F(Form("hl_%d",i),"pdf;LE time [ns]; entries [#]", 1000,0,100);
    hnph[i] = new TH1F(Form("hnph_%d",i),";photon yield [#]; entries [#]", 250,0,250);
    hll[i] = new TH1F(Form("hll_%d",i),"hll;ln L(p) - ln L(#pi); entries [#]",110,-30,30);
  }  
  TH1F *hl3 = new TH1F("hl3","pdf;LE time [ns]; entries [#]", 1000,0,100);

  TRandom rand;
  TF1 *pdff[maxch_dirc],*pdfs[maxch_dirc];
  TString pdf = path;
  pdf.ReplaceAll(".root",pdfEnding);
  TFile f(pdf);
  
  if(sigma >0) hl3->Rebin((Int_t)sigma);
  Int_t integ1(0), integ2(0);
  for(Int_t i=0; i<maxch_dirc; i++){
    hpdff[i] = (TH1F*)f.Get(Form("hf_%d",i));
    hpdfs[i] = (TH1F*)f.Get(Form("hs_%d",i));
    hpdff[i]->SetLineColor(2);
    hpdfs[i]->SetLineColor(4);
    if(sigma >0) hpdff[i]->Rebin((Int_t)sigma);
    if(sigma >0) hpdfs[i]->Rebin((Int_t)sigma);
    integ1+= hpdff[i]->Integral();
    integ2+= hpdfs[i]->Integral();
    hl3->Add(hpdff[i]);
    hl3->Add(hpdfs[i]);
  }
  if(path.Contains("C.root")) sigma=0;
  if(path.Contains("Z.root")) sigma=0;

  TF1 *F1 = new TF1("gaus0","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2])",0,150);
  TF1 *F2 = new TF1("gaus1","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2])",0,150);
  F1->SetParameter(0,1);
  F2->SetParameter(0,1);
      
  F1->SetParameter(1,63);
  F2->SetParameter(1,50);
  F1->SetParameter(2,11);
  F2->SetParameter(2,9);
  
  Double_t theta(0);
  TVirtualFitter *fitter;
  Double_t mom,nph,time,timeres(-1);
  PrtHit fHit;
  Int_t totalf(0),totals(0),mcp,pix,ch, entries = fCh->GetEntries(); // [50000-rest] - is for pdf generation
  if(path.Contains("F.root")) entries = fCh->GetEntries();
  for (Int_t ievent=0; ievent<entries; ievent++){
    PrtNextEvent(ievent,1000);
    if(ievent==0){
      theta = prt_event->GetAngle();
      mom=prt_event->GetMomentum().Mag();
    }
    
    timeres = prt_event->GetTimeRes();
    Double_t aminf,amins, sum(0),sumf(0),sums(0);
    Int_t pid = prt_event->GetParticle();
    Int_t nHits =prt_event->GetHitSize();    
    if(nHits<5) continue;
    
    if(prt_event->GetType()==0){
      // if(fabs(prt_event->GetMomentum().Mag()-7)<0.1){
      // 	if( prt_event->GetParticle()==2212 && prt_event->GetTest1()<175.6 ) continue;
      // 	if( prt_event->GetParticle()==211  && prt_event->GetTest1()>175.1 ) continue;
      // }

      Bool_t tof1(false), tof2(false);
      Bool_t hodo1(false), hodo2(false);
      for(Int_t h=0; h<nHits; h++) {
  	fHit = prt_event->GetHit(h);
  	Int_t gch=fHit.GetChannel();
	
	//if(gch>1031 && gch<1034)
	  tof1=true;
	//if(gch>1060)
	  tof2=true;
	
	  //if(gch>776 && gch<780)
	  hodo1=true;
	  //if(gch>790 && gch<794)
	  hodo2=true;	   
      }
      
      if(!(tof1 && tof2 && hodo1 && hodo2)) continue;
    }

    hnph[prt_pid]->Fill(nHits);

    if(debug) std::cout<<"===================== event === "<< ievent <<std::endl;
    
    for(Int_t i=0; i<nHits; i++){
      fHit = prt_event->GetHit(i);
      mcp = fHit.GetMcpId();
      pix=fHit.GetPixelId()-1;
      ch = map_mpc[mcp][pix];
      time = fHit.GetLeadTime()+rand.Gaus(0,sigma/10.);
      if(mcp>6) continue;
      
      { //time cuts
      	// Double_t cut1(11);
      	// if(studyId==157 || studyId==155) cut1=8;
      	// if(theta<=80){
      	//   if(time<cut1 || time>75) continue; //45
      	// }else if(theta>94){
      	//   if(time<3 || time>40) continue; //40
      	// }
	if(time<9 || time >40) continue;				 
      }

      aminf = hpdff[ch]->GetBinContent(hpdff[ch]->FindBin(time-0.)); 
      amins = hpdfs[ch]->GetBinContent(hpdfs[ch]->FindBin(time+0.));   

      if(debug){
	TString x=(aminf>amins)? " <====== PROTON" : "";
	std::cout<<Form("f %1.6f s %1.6f mcp %d pix %d   pid %d",aminf,amins,mcp,pix  ,pid)<<"  "<<x <<std::endl;
	
	cc->cd();
	hpdff[ch]->Draw();
	hpdff[ch]->GetXaxis()->SetRangeUser(0,50);
	hpdfs[ch]->Draw("same");
	cc->Update();
	gLine->SetX1(time);
	gLine->SetX2(time);
	gLine->SetY1(cc->GetUymin());
	gLine->SetY2(cc->GetUymax());
	gLine->Draw();
	cc->Update();
	cc->WaitPrimitive();
      }
      // if(aminf==0 || amins==0) continue;
      Double_t noise = nHits * 1e-5; //1e-7;
      sumf+=TMath::Log((aminf+noise));
      sums+=TMath::Log((amins+noise));    

      hl[prt_pid]->Fill(time);
    }

    sum = sumf-sums; 
    // if(fabs(sum)<0.04) continue;
    
    // sumf += 10*TMath::Log(F2->Eval(nHits));
    // sums += 10*TMath::Log(F1->Eval(nHits));
    // sum = sumf-sums;
    
    hll[prt_pid]->Fill(sum);
  }

  TString name = Form("tis_%d_%1.1f_%1.1f_m%1.1f.root",studyId,theta,sigma,mom);
  if(path.Contains("C.root")) name = Form("tid_%d_%1.1f_%1.1f_m%1.1f.root",studyId,theta,sigma,mom);
  canvasAdd("ll_"+name,800,500);

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  prt_normalize(hll[4],hll[2]);
  hll[4]->GetYaxis()->SetNdivisions(9,5,0);
   
  TF1 *ff;
  Double_t m1,m2,s1,s2,dm1,dm2,ds1,ds2; 
  if(hll[4]->GetEntries()>10){
    hll[4]->Fit("gaus","Sq");
    ff = hll[4]->GetFunction("gaus");
    m1=ff->GetParameter(1);
    s1=ff->GetParameter(2);
    dm1=ff->GetParError(1);
    ds1=ff->GetParError(2);
  }
  if(hll[2]->GetEntries()>10){
    hll[2]->Fit("gaus","Sq");
    ff = hll[2]->GetFunction("gaus");
    m2=ff->GetParameter(1);
    s2=ff->GetParameter(2);
    dm2=ff->GetParError(1);
    ds2=ff->GetParError(2);
  }
  Double_t sep = (fabs(m1-m2))/(0.5*(s1+s2));
  
  hll[4]->SetTitle(Form("separation = %1.2f",sep));
  
  hll[4]->SetLineColor(2);
  hll[4]->Draw();
  hll[2]->SetLineColor(4);
  hll[2]->Draw("same");

  TLegend *leg = new TLegend(0.6,0.7,0.8,0.87);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(hll[2],"pions ","lp");
  leg->AddEntry(hll[4],"protons","lp");
  leg->Draw();
  
  canvasAdd("hl_"+name,800,500);


  // hl[4]->Scale(1/hl[4]->GetMaximum());
  // hl[2]->Scale(1/hl[2]->GetMaximum());
  // hl3->Scale(1/hl3->GetMaximum());

  prt_normalize(hl[4],hl[2]);
  
  hl[4]->Draw();
  hl[2]->SetLineColor(4);
  hl[2]->Draw("same");
  hl3->SetLineColor(2);
  hl3->Draw("same");
  
  canvasAdd("hnph_"+name,800,500);
  prt_normalize(hnph[4],hnph[2]);
  hnph[4]->Fit("gaus");
  ff = hnph[4]->GetFunction("gaus");
  nph=ff->GetParameter(1);
  
  hnph[4]->SetLineColor(4);
  hnph[4]->Draw();
  hnph[2]->SetLineColor(2);
  hnph[2]->Draw("same");
  
  canvasSave(0,0);

  
  Double_t e1,e2,e3,e4;

  std::cout<<dm1<<" "<<dm2<<" "<<ds1 <<" "<<ds2<<std::endl;
  
  
  e1=2/(s1+s2)*dm1;
  e2=2/(s1+s2)*dm2;
  e3=-((2*(m1 + m2))/((s1 + s2)*(s1 + s2)))*ds1;
  e4=-((2*(m1 + m2))/((s1 + s2)*(s1 + s2)))*ds2;

  Double_t esep=sqrt(e1*e1+e2*e2+e3*e3+e4*e4);
  std::cout<<path<<" separation "<< sep <<" +/- "<<esep <<std::endl;

  TFile fc(fSavePath+"/reco_"+name,"recreate");
  TTree *tc = new TTree("reco","reco");
  tc->Branch("theta",&theta,"theta/D");
  tc->Branch("sep",&sep,"sep/D");
  tc->Branch("esep",&esep,"esep/D");
  tc->Branch("sigma",&sigma,"sigma/D");
  tc->Branch("mom",&mom,"mom/D");
  tc->Branch("nph",&nph,"nph/D");
  sigma/=10.;
  tc->Fill();
  tc->Write();
}
