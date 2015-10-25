#define prt__beam
#include "../prtdirc/src/PrtHit.h"
#include "../prtdirc/src/PrtEvent.h"
#include "prttools.C"
#include "datainfo.C"
#include <TEllipse.h>


Double_t walktheta(-13*TMath::Pi()/180.);
Double_t fy1(42.25), fy2(44.05);

Double_t fr11[10]={0.5,0.5,0.3,0.3,0.4, 0.3,0.3,0.2,0.20,0.15};
Double_t fr12[10]={1.0,1.0,0.9,0.9,0.9, 0.9,0.9,0.8,0.80,0.70};

Double_t fx1[200][10];
Double_t fx2[200][10];


Double_t fr21[10]={0.8,0.8,0.3,0.3,0.4, 0.3,0.3,0.2,0.2,0.2};
Double_t fr22[10]={1.0,1.0,0.9,0.9,0.9, 0.9,0.9,0.8,0.8,0.8};
Int_t gmom(7);


Bool_t insideOfEllipce(Double_t x, Double_t y, Double_t x0, Double_t y0,  Double_t r1, Double_t r2, Double_t w=0){

  Double_t xx = cos(w)*(x-x0)+sin(w)*(y-y0);
  Double_t yy = sin(w)*(x-x0)-cos(w)*(y-y0);

  return xx*xx/(r1*r1)+yy*yy/(r2*r2)<=1;
}

void drawTof(TString infile="hits.root"){

  TString fileid(infile);
  fileid.Remove(0,fileid.Last('/')+1);
  fileid.Remove(fileid.Last('.')-4);
  DataInfo di = getDataInfo(fileid);
  Int_t momentum = di.getMomentum();
  Int_t studyId = di.getStudyId();
  gmom=momentum-1;
  fSavePath = Form("tof/%d/%d_",studyId,momentum)+fileid;
  PrtInit(infile,1);

  Int_t le1(170), le2(186);
  Int_t l1(174), l2(177);

  // fx1[170][]={564,646,474};
  
  fx1[150][7]=174.95;
  fx2[150][7]=175.70;

  fx1[151][7]=174.95;
  fx2[151][7]=175.70;

  fx1[152][7]=175.0;  
  fx2[152][7]=175.8;  

  fx1[153][7]=175.05;  
  fx2[153][7]=175.85;
  
  fx1[154][7]=175.3;  
  fx2[154][7]=176.0;  
  fx1[155][7]=175.3;  
  fx2[155][7]=176.05; 
  fx1[156][7]=175.3;  
  fx2[156][7]=176.05;
  fx1[157][7]=175.2;  
  fx2[157][7]=175.95;
  fx1[158][7]=175.1;  
  fx2[158][7]=175.9;  
  fx1[159][7]=175.1;  
  fx2[159][7]=175.9;  

  fx1[160][5]=175.2;
  fx2[160][5]=176.7;


  fx1[170][2]=175.2;
  fx1[170][3]=175.2;
  fx1[170][4]=175.1;
  fx1[170][5]=175.18;
  fx1[170][6]=175.18;
  fx1[170][7]=175.17;
  fx1[170][8]=175.18;
  fx1[170][9]=175.10;
  fx1[170][10]=175.00;

  fx2[170][2]=184.3;
  fx2[170][3]=179.4;
  fx2[170][4]=177.5;
  fx2[170][5]=176.68;
  fx2[170][6]=176.22;
  fx2[170][7]=175.95;
  fx2[170][8]=175.76;
  fx2[170][9]=175.58;
  fx2[170][10]=175.44;

  fx1[171][3]=175.1;
  fx1[171][4]=175.1;
  fx1[171][5]=175.1;
  fx1[171][6]=174.9;
  fx1[171][7]=174.9;
  fx1[171][8]=174.9;
  fx1[171][9]=174.9;
  fx1[171][10]=174.90;
  
  fx2[171][3]=179.3;
  fx2[171][4]=177.4;
  fx2[171][5]=176.6;
  fx2[171][6]=176.1;
  fx2[171][7]=175.8;
  fx2[171][8]=175.6;
  fx2[171][9]=175.5;
  fx2[171][10]=175.35;
 
  fx1[172][3]=175.1;
  fx1[172][4]=175.1;
  fx1[172][5]=175.1;
  fx1[172][6]=175.0;
  fx1[172][7]=175.0;
  
  fx2[172][3]=179.3;
  fx2[172][4]=177.5;
  fx2[172][5]=176.6;
  fx2[172][6]=176.1;
  fx2[172][7]=175.8;
 
  fx1[173][3]=175.2;
  fx1[173][4]=175.1;
  fx1[173][5]=175.1;
  fx1[173][6]=175.1;
  fx1[173][7]=175.10;
  fx1[173][8]=175.0;
  fx1[173][9]=175.0;
  fx1[173][10]=174.90;
  fx2[173][3]=179.3;
  fx2[173][4]=177.5;
  fx2[173][5]=176.6;
  fx2[173][6]=176.1;
  fx2[173][7]=175.9;
  fx2[173][8]=175.6;
  fx2[173][9]=175.5;
  fx2[173][10]=175.40;

 
  fx1[174][3]=175.4;
  fx1[174][4]=175.3;
  fx1[174][5]=175.3;
  fx1[174][6]=175.3;
  fx1[174][7]=175.3;
  fx1[174][8]=175.15;
  fx1[174][9]=175.15;
  fx1[174][10]=175.15;
  fx2[174][3]=179.5;
  fx2[174][4]=177.7;
  fx2[174][5]=176.8;
  fx2[174][6]=176.3;
  fx2[174][7]=176.1;
  fx2[174][8]=175.80;
  fx2[174][9]=175.70;
  fx2[174][10]=175.60;

 
  fx1[175][4]=175.1;
  fx1[175][6]=175.18;
  fx2[175][4]=177.5;
  fx2[175][6]=176.22;

 
  fx1[176][4]=175.2;
  fx1[176][6]=175.1;
  fx2[176][4]=177.5;
  fx2[176][6]=176.1;

  
  fx1[177][4]=175.2;
  fx1[177][6]=175.1;
  fx2[177][4]=177.5;
  fx2[177][6]=176.2;

 
  fx1[178][1]=175.0;
  fx2[178][1]=175.4;
 
  fx1[179][2]=175.3;
  fx1[179][3]=175.2;
  fx1[179][4]=175.1;
  fx1[179][5]=175.1;
  fx1[179][6]=175.1;
  fx1[179][7]=175.1;
  fx1[179][8]=175.1;
  fx1[179][9]=175.1;
  fx1[179][10]=175.0;
  fx2[179][2]=184.4;
  fx2[179][3]=179.4;
  fx2[179][4]=177.5;
  fx2[179][5]=176.6;
  fx2[179][6]=176.1;
  fx2[179][7]=175.8;
  fx2[179][8]=175.6;
  fx2[179][9]=175.5;
  fx2[179][10]=175.4;

  fx1[180][7]=175.3;
  fx2[180][7]=176.0;


  if(momentum<=5) {
    l2=186;
    fy1=42.3;
  }

  Int_t id = studyId;

  if(studyId>=160 && studyId<170) id = 160;
  if(studyId>=180) id = 180;
  
  Double_t tpi =   fx1[id][momentum];
  Double_t tp =    fx2[id][momentum];
  Double_t cpiy = fr11[momentum];
  Double_t cpix = fr12[momentum];
  Double_t cpy =  fr21[momentum];
  Double_t cpx =  fr22[momentum];

  if(studyId==170) fy2=44.17;
  if(studyId==171) fy2=44.17;
  if(studyId==172) fy2=44.2;
  if(studyId==173) fy2=44.13;
  if(studyId==174) fy2=44.4;
  if(studyId==179) fy2=44.4;
  if(studyId>156 && studyId<160) {
    fy2=42.35;
    fy2=44.17;
  }

  if(studyId==151){
    if(infile.Contains("15177152") || infile.Contains("15177145")|| infile.Contains("15177143")){
      tpi = 175.1;
      tp = 175.9;
    }
    if(infile.Contains("15177141") || infile.Contains("15177135")){
      tpi = 175.0;
      tp = 175.8;
    }
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
      
      // if(insideOfEllipce(time, tot1, tpi, fy1, cpiy, cpix) && insideOfEllipce(time, tot2, tpi, fy2, cpiy, cpix)){
	hLeTotC->Fill(time,tot1);
	hLeTotC2->Fill(time,tot2);
      // }else if(insideOfEllipce(time, tot1, tp, fy1, cpy, cpx) && insideOfEllipce(time, tot2, tp, fy2, cpy, cpx)){
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
  TEllipse *el1 = new TEllipse(tpi,fy1,cpiy,cpix);
  el1->SetLineColor(2);
  el1->SetLineWidth(2);
  el1->SetFillStyle(0);
  el1->Draw();
  TEllipse *el2 = new TEllipse(tp,fy1,cpy,cpx);
  el2->SetLineColor(2);
  el2->SetLineWidth(2);
  el2->SetFillStyle(0);
  el2->Draw();
  
  canvasAdd("LeTotC2",800,400);
  hLeTotC2->Draw("colz");
  TEllipse *el3 = new TEllipse(tpi,fy2,cpiy,cpix);
  el3->SetLineColor(2);
  el3->SetLineWidth(2);
  el3->SetFillStyle(0);
  el3->Draw();
  TEllipse *el4 = new TEllipse(tp,fy2,cpy,cpx);
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
