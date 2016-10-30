#define prt__beam
#include "../prtdirc/src/PrtHit.h"
#include "../prtdirc/src/PrtEvent.h"
#include "prttools.C"
#include <TEllipse.h>
#include <TKey.h>

Double_t walktheta(-5*TMath::Pi()/180.);
Double_t tof1le(0),tof2le(0),tof1tot(0),tof2tot(0);

Double_t fr11[11]={0,0.5,0.5,0.3,0.3,0.4, 0.3,0.3,0.2,0.20,0.15};
Double_t fr12[11]={0,1.0,1.0,0.9,0.9,0.9, 0.9,0.9,0.8,0.80,0.70};
Double_t fr21[11]={0,0.8,0.8,0.3,0.3,0.4, 0.3,0.3,0.2,0.2,0.2};
Double_t fr22[11]={0,1.0,1.0,0.9,0.9,0.9, 0.9,0.9,0.8,0.8,0.8};

Double_t c1y(0.5),c2y(0.5),c1x(0.9),c2x(0.9);

Bool_t insideOfEllipce(Double_t x, Double_t y, Double_t x0, Double_t y0,  Double_t r1, Double_t r2, Double_t w=0){

  Double_t xx = cos(w)*(x-x0)+sin(w)*(y-y0);
  Double_t yy = sin(w)*(x-x0)-cos(w)*(y-y0);

  return xx*xx/(r1*r1)+yy*yy/(r2*r2)<=1;
}

void drawTof(TString infile="hits.root",TString gcFile="calib_2610.root"){
  if(gcFile!="0"){
    TFile f(gcFile);
    TIter nextkey(f.GetListOfKeys());
    TKey *key;
    
    while ((key = (TKey*)nextkey())) {
      TGraph *gr = (TGraph*)key->ReadObj();
      TString name = gr->GetName();
      if(name.Contains("tof")){
	name.Remove(0,4);
	if(infile.Contains(name)){
	  gr->GetPoint(0,tof1le,tof2le);
	  gr->GetPoint(1,tof1tot,tof2tot);
	}
      }
    }
    f.Close(); 
  }
  
  TString fileid(infile);
  fileid.Remove(0,fileid.Last('/')+1);
  fileid.Remove(fileid.Last('.')-1);
  prt_data_info = getDataInfo(fileid);
  Int_t momentum = prt_data_info.getMomentum();
  Int_t studyId = prt_data_info.getStudyId();
  if(studyId<0) studyId=1000;
  fSavePath = Form("tof/%d/%d",studyId,prt_data_info.getFileId());
  PrtInit(infile,1);
  std::cout<<"fSavePath  "<<fSavePath <<std::endl;
  if(studyId<0) return;
  
  Int_t le1(70), le2(75);
  Int_t l1(70), l2(75);

  if(studyId==171 && momentum==7){
    tof1le=174.98;
    tof2le=175.74;
  }
  c1y=fr11[momentum];
  c2y=fr21[momentum];
  c1x=fr12[momentum];
  c2x=fr22[momentum];

  TH1F * hMult[4];

  TString names[]={"trigger_1","trigger_2","tof_1","tof_2"};
  Int_t colors[]={2,4,1,8};
  for(Int_t i=0; i<4; i++){
    hMult[i]  = new TH1F(names[i],";multiplicity [#]; entries [#]",10,0,10);
    hMult[i]->SetLineColor(colors[i]);
  }
      
  TH1F * hTof1  = new TH1F("tof1 ","tof1;TOF2-TOF1 [ns]; entries [#]",600,-1000,1000);
  TH1F * hTof2  = new TH1F("tof2 ","tof2;TOF2-TOF1 [ns]; entries [#]",600,-1000,1000);
  TH1F * hTof  = new TH1F("tof ","tof;TOF2-TOF1 [ns]; entries [#]",   600,le1,le2);
  TH1F * hTofC  = new TH1F("tofC ","tofC;TOF2-TOF1 [ns]; entries [#]",600,le1,le2);
  TH1F * hTot  = new TH1F("tot ","tot;TOT1,TOT2 [ns]; entries [#]",   600,0,100);
 
  TH2F * hLeTot  = new TH2F("letot ","letot;TOF2-TOF1 [ns]; TOT1 [ns]",      500,l1,l2,125,43,47);
  TH2F * hLeTotW  = new TH2F("letotW ","letotW;TOF2-TOF1 [ns]; TOT1 [ns]",   500,l1,l2,125,43,47);
  TH2F * hLeTotC  = new TH2F("letotC ","letotC;TOF2-TOF1 [ns]; TOT1 [ns]",   500,l1,l2,125,43,47);
  TH2F * hLeTotC2  = new TH2F("letotC2 ","letotC2;TOF2-TOF1 [ns]; TOT2 [ns]",500,l1,l2,125,43,47);


  gStyle->SetOptStat(1001111);
  gStyle->SetOptFit();
  
  PrtHit fHit;
  for (Int_t ievent=0; ievent< fCh->GetEntries(); ievent++){ //fCh->GetEntries()
    PrtNextEvent(ievent,10000);
    
    Bool_t btrig(false),bmcpout(false),btof1(false),btof2(false);
    Double_t tot1(0),tot2(0),tof1(0),tof2(0);
    Int_t mult1(0), mult2(0), mult3(0), mult4(0);
    
    for(Int_t i=0; i<prt_event->GetHitSize(); i++){
      fHit = prt_event->GetHit(i);
      if(fHit.GetChannel()==818) {//trigger 1
	btrig = true;
	mult1++;
      }
      if(fHit.GetChannel()==821){ // trigger 2
	bmcpout = true;
	mult2++;
      }
      if(fHit.GetChannel()==720){ //tof 1
	btof1 = true;
	tof1=fHit.GetLeadTime();
	tot1=fHit.GetTotTime();
	mult3++;
      }
      if(fHit.GetChannel()==722){ //tof2
	btof2 = true;
	tof2 = fHit.GetLeadTime();
	tot2=fHit.GetTotTime();
	mult4++;
      }
    }
    
    if(!(btrig && bmcpout && btof1 && btof2)) continue;
    //if(fabs(tot1-tof1tot)>0.5 ||fabs(tot2-tof2tot)>0.5 )  continue;

    hMult[0]->Fill(mult1);
    hMult[1]->Fill(mult2);
    hMult[2]->Fill(mult3);
    hMult[3]->Fill(mult4);
    
    hTof1->Fill(tof1);
    hTof2->Fill(tof2);

    if(tof1!=0 && tof2!=0){
      Double_t time = tof2-tof1;
      hTof->Fill(time);	
      hLeTotW->Fill(time,tot1);
      time += (tot1-tof1tot)*tan(walktheta);
      time += (tot2-tof2tot)*tan(-walktheta);
      hTofC->Fill(time);
     
      // if(insideOfEllipce(time, tot1, tof1le, tof1tot, c1y, c1x) && insideOfEllipce(time, tot2, tof1le, tof2tot, c1y, c1x)){
   	hLeTotC->Fill(time,tot1);
	hLeTotC2->Fill(time,tot2);
      // }else if(insideOfEllipce(time, tot1, tof2le, tof1tot, c2y, c2x) && insideOfEllipce(time, tot2, tof2le, tof2tot, c2y, c2x)){
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
  TEllipse *el1 = new TEllipse(tof1le,tof1tot,c1y,c1x);
  el1->SetLineColor(2);
  el1->SetLineWidth(2);
  el1->SetFillStyle(0);
  el1->Draw();
  TEllipse *el2 = new TEllipse(tof2le,tof1tot,c2y,c2x);
  el2->SetLineColor(2);
  el2->SetLineWidth(2);
  el2->SetFillStyle(0);
  el2->Draw();
  
  canvasAdd("LeTotC2",800,400);
  hLeTotC2->Draw("colz");
  TEllipse *el3 = new TEllipse(tof1le,tof2tot,c1y,c1x);
  el3->SetLineColor(2);
  el3->SetLineWidth(2);
  el3->SetFillStyle(0);
  el3->Draw();
  TEllipse *el4 = new TEllipse(tof2le,tof2tot,c2y,c2x);
  el4->SetLineColor(2);
  el4->SetLineWidth(2);
  el4->SetFillStyle(0);
  el4->Draw();
 
  
  canvasAdd("tof",800,400);
  prt_fit(hTof,10,20,2,2);
  hTof->Draw();

  canvasAdd("tofC",800,400);
  TVector3 r=prt_fit(hTofC,10,20,2,2);
  hTofC->Draw();
  Double_t tot1=hLeTotC->GetMean(2);
  Double_t tot2=hLeTotC2->GetMean(2);
  std::cout<<infile  << "  LE1  "<<r.X() << "   LE2 "<< r.Z() << " TOT1 "<< tot1 <<"  TOT2 "<< tot2<<std::endl;

  
  fileid=infile;
  fileid.Remove(0,fileid.Last('_')+1);
  fileid.Remove(fileid.Last('C'));
  TFile efile(infile+".tof.root","RECREATE");
  TGraph *gr = new TGraph();
  gr->SetPoint(0,r.X(),r.Z());
  gr->SetPoint(1,tot1,tot2);
  gr->SetName("tof_"+fileid);
  gr->Write();
  efile.Write();
  efile.Close();
  
  
  canvasAdd("mult",800,400);
  hMult[0]->Draw();
  for(Int_t i=1; i<4; i++){
    hMult[i]->Draw("same");
  }
  
  TLegend *leg = new TLegend(0.5,0.6,0.8,0.85);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  for(Int_t i=0; i<4; i++) leg->AddEntry(hMult[i],names[i],"lp");
  leg->Draw();
  
  gStyle->SetOptTitle(0);
  canvasSave(0,0);
}
