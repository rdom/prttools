#define prt__beam
#include "../prtdirc/src/PrtHit.h"
#include "../prtdirc/src/PrtEvent.h"
#include "prttools.C"
#include <TEllipse.h>


// Double_t walktheta(-16.5*TMath::Pi()/180.);
// Double_t deltatheta(2.5*TMath::Pi()/180.);
// Double_t walky0(42.23), walky1(44.00);
// Double_t walk1x0(174.85), walk2x0(175.64);

//7
// Double_t walktheta(-10.5*TMath::Pi()/180.);
// Double_t deltatheta(2.5*TMath::Pi()/180.);
// Double_t walky0(42.23), walky1(44.00);
// Double_t walk1x0(174.85), walk2x0(175.64);


//10
Double_t walktheta(-13*TMath::Pi()/180.);
Double_t fy1(42.2), fy2(44.00);

Double_t fx1[10]={175.0,175.2,175.2,175.1,175.18, 175.18,175.17,175.18,175.10,175.00};
Double_t fx2[10]={175.4,184.3,179.4,177.5,176.68, 176.22,175.95,175.76,175.58,175.44};


Double_t fr11[10]={0.5,0.5,0.3,0.3,0.4, 0.3,0.3,0.2,0.20,0.15};
Double_t fr12[10]={1.0,1.0,0.9,0.9,0.9, 0.9,0.9,0.8,0.80,0.70};

Double_t fr21[10]={0.8,0.8,0.3,0.3,0.4, 0.3,0.3,0.2,0.20,0.20};
Double_t fr22[10]={1.0,1.0,0.9,0.9,0.9, 0.9,0.9,0.8,0.80,0.70};
Int_t gmom(7);


Bool_t insideOfEllipce(Double_t x, Double_t y, Double_t x0, Double_t y0,  Double_t r1, Double_t r2, Double_t w=0){

  Double_t xx = cos(w)*(x-x0)+sin(w)*(y-y0);
  Double_t yy = sin(w)*(x-x0)-cos(w)*(y-y0);

  return xx*xx/(r1*r1)+yy*yy/(r2*r2)<=1;
}

void drawTof(TString infile="hits.root", Int_t momentum=7){
  gmom=momentum-1;
  TString fileid(infile);
  fileid.Remove(0,fileid.Last('/')-3);
  fileid.Remove(fileid.Last('.'));
  fSavePath = Form("tof/mom_%d",momentum)+fileid;
  PrtInit(infile,1);

  Int_t le1(170), le2(186);
  Int_t l1(174), l2(177);
  if(momentum<=5) {
    l2=186;
    fy1=42.3;
  }
  
  
  TH1F * hMult  = new TH1F("mult","mult",1000,-1000,1000);
  TH1F * hTof1  = new TH1F("tof1 ","tof1;TOT2-TOF1 [ns]; entries [#]",1000,-1000,1000);
  TH1F * hTof2  = new TH1F("tof2 ","tof2;TOT2-TOF1 [ns]; entries [#]",1000,-1000,1000);
  TH1F * hTof  = new TH1F("tof ","tof;TOT2-TOF1 [ns]; entries [#]",   1000,le1,le2);
  TH1F * hTofC  = new TH1F("tofC ","tofC;TOT2-TOF1 [ns]; entries [#]",1000,le1,le2);
  TH1F * hTot  = new TH1F("tot ","tot;TOT1,TOT2 [ns]; entries [#]",   1000,0,100);
 
  TH2F * hLeTot  = new TH2F("letot ","letot;TOT2-TOF1 [ns]; TOT1 [ns]",      500,l1,l2,200,40,44);
  TH2F * hLeTotW  = new TH2F("letotW ","letotW;TOT2-TOF1 [ns]; TOT1 [ns]",   500,l1,l2,200,40,44);
  TH2F * hLeTotC  = new TH2F("letotC ","letotC;TOT2-TOF1 [ns]; TOT1 [ns]",   500,l1,l2,200,40,44);
  TH2F * hLeTotC2  = new TH2F("letotC2 ","letotC2;TOT2-TOF1 [ns]; TOT2 [ns]",500,l1,l2,200,41,45);


  gStyle->SetOptStat(1001111);
  gStyle->SetOptFit();
  
  PrtHit fHit;
  for (Int_t ievent=0; ievent< fCh->GetEntries(); ievent++){ //fCh->GetEntries()
    PrtNextEvent(ievent,10000);
    
    Bool_t btrig(false),bmcpout(false),btof1(false),btof2(false);
    Double_t tot1(0),tot2(0),tof1(0),tof2(0);

    for(Int_t i=0; i<fEvent->GetHitSize(); i++){
      fHit = fEvent->GetHit(i);
      if(fHit.GetChannel()==1344) btrig = true;
      if(fHit.GetChannel()==1248) bmcpout = true;

      if(fHit.GetChannel()==960){
	btof1 = true;
	tof1 = fHit.GetLeadTime();
	tot1=fHit.GetTotTime();
      }
      if(fHit.GetChannel()==1104) {
	btof2 = true;
	tof2 = fHit.GetLeadTime();
	tot2=fHit.GetTotTime();
      }
    }
    
    if(!(btrig && btof1 && btof2 && bmcpout)) continue;
 
    hTof1->Fill(tof1);
    hTof2->Fill(tof2);

    if(tof1!=0 && tof2!=0){
      Double_t time = tof2-tof1;
      hTof->Fill(time);	
      hLeTotW->Fill(time,tot1);
      time += (tot1-fy1)*tan(walktheta);
      time += (tot2-fy2)*tan(-walktheta);
      hTofC->Fill(time);
      
      // if(insideOfEllipce(time, tot1, fx1[gmom], fy1, fr11[gmom], fr12[gmom]) && insideOfEllipce(time, tot2, fx1[gmom], fy2, fr11[gmom], fr12[gmom])){
	hLeTotC->Fill(time,tot1);
	hLeTotC2->Fill(time,tot2);
      // }else if(insideOfEllipce(time, tot1, fx2[gmom], fy1, fr21[gmom], fr22[gmom]) && insideOfEllipce(time, tot2, fx2[gmom], fy2, fr21[gmom], fr22[gmom])){
      // 	hLeTotC->Fill(time,tot1);
      // 	hLeTotC2->Fill(time,tot2);
      // }
	
      hLeTot->Fill(tof2-tof1,tot1);
    }
      
  }

  canvasAdd("LeTot",800,400);
  hLeTot->Draw("colz");

  // canvasAdd("LeTotW",800,400);
  // hLeTotW->Draw("colz");
  
  canvasAdd("LeTotC",800,400);
  hLeTotC->Draw("colz");
  TEllipse *el1 = new TEllipse(fx1[gmom],fy1,fr11[gmom],fr12[gmom]);
  el1->SetLineColor(2);
  el1->SetLineWidth(2);
  el1->SetFillStyle(0);
  el1->Draw();
  TEllipse *el2 = new TEllipse(fx2[gmom],fy1,fr21[gmom],fr22[gmom]);
  el2->SetLineColor(2);
  el2->SetLineWidth(2);
  el2->SetFillStyle(0);
  el2->Draw();
  
  canvasAdd("LeTotC2",800,400);
  hLeTotC2->Draw("colz");
  TEllipse *el3 = new TEllipse(fx1[gmom],fy2,fr11[gmom],fr12[gmom]);
  el3->SetLineColor(2);
  el3->SetLineWidth(2);
  el3->SetFillStyle(0);
  el3->Draw();
  TEllipse *el4 = new TEllipse(fx2[gmom],fy2,fr21[gmom],fr22[gmom]);
  el4->SetLineColor(2);
  el4->SetLineWidth(2);
  el4->SetFillStyle(0);
  el4->Draw();
 
  
  canvasAdd("tof",800,400);
  prt_fit(hTof,10,20,2,2);
  hTof->Draw();

  canvasAdd("tofC",800,400);
  prt_fit(hTofC,10,20,2,2);
  hTofC->Draw();

  // canvasAdd("tot",800,400);
  // hTot->Draw();
  
  canvasSave(1,0);
}
