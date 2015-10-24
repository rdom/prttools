#define prt__beam
#include "../prtdirc/src/PrtHit.h"
#include "../prtdirc/src/PrtEvent.h"
#include "prttools.C"

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
Double_t walktheta(-12*TMath::Pi()/180.);
Double_t walky0(42.23), walky1(44.00);
Double_t walk1x0(175), walk2x0(175.4);

Double_t fx1[10]={175.4,175.4,175.4,175.4,175.4,175.4,175.4,175.4,175.4,175.4};
Double_t fx2[10]={175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0,175.0};

Double_t fr11[10]={0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.20};
Double_t fr12[10]={0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.65};

Double_t fr21[10]={0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.20};
Double_t fr22[10]={0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.9,0.65};

Bool_t insideOfEllipce(Double_t x, Double_t y, Double_t x0, Double_t y0, Double_t w){
  Double_t r1(0.3), r2(0.9);

  Double_t xx = cos(w)*(x-x0)+sin(w)*(y-y0);
  Double_t yy = sin(w)*(x-x0)-cos(w)*(y-y0);

  return xx*xx/(r1*r1)+yy*yy/(r2*r2)<=1;
}

void drawTof(TString infile="hits.root"){
  fSavePath = "auto";
  PrtInit(infile,1);

  TH1F * hMult  = new TH1F("mult","mult",1000,-1000,1000);
  TH1F * hTof1  = new TH1F("tof1 ","tof1;TOT2-TOF1 [ns]; entries [#]",1000,-1000,1000);
  TH1F * hTof2  = new TH1F("tof2 ","tof2;TOT2-TOF1 [ns]; entries [#]",1000,-1000,1000);
  TH1F * hTof  = new TH1F("tof ","tof;TOT2-TOF1 [ns]; entries [#]",500,170,180);
  TH1F * hTofC  = new TH1F("tofC ","tofC;TOT2-TOF1 [ns]; entries [#]",1000,170,180);
  TH1F * hTot  = new TH1F("tot ","tot;TOT1,TOT2 [ns]; entries [#]",1000,0,100);
 
  TH2F * hLeTot  = new TH2F("letot ","letot;TOT2-TOF1 [ns]; TOT1 [ns]",500,174,176.5,500,40,45);
  TH2F * hLeTotW  = new TH2F("letotW ","letotW;TOT2-TOF1 [ns]; TOT1 [ns]",500,174,176.5,500,40,45);
  TH2F * hLeTotC  = new TH2F("letotC ","letotC;TOT2-TOF1 [ns]; TOT1 [ns]",500,174,176.5,500,40,45);
  TH2F * hLeTotC2  = new TH2F("letotC2 ","letotC2;TOT2-TOF1 [ns]; TOT2 [ns]",500,174,176.5,500,40,45);


  gStyle->SetOptStat(1001111);
  gStyle->SetOptFit();
  
  Double_t tof1, tof2;
  PrtHit fHit;
  for (Int_t ievent=0; ievent< fCh->GetEntries(); ievent++){ //fCh->GetEntries()
    PrtNextEvent(ievent,1000);
    
    Bool_t btrig(false),bmcpout(false),btof1(false),btof2(false);
    for(Int_t i=0; i<fEvent->GetHitSize(); i++){
      fHit = fEvent->GetHit(i);
      if(fHit.GetChannel()==1344) btrig = true;
      if(fHit.GetChannel()==960)  btof1 = true;
      if(fHit.GetChannel()==1104) btof2 = true;
      if(fHit.GetChannel()==1248) bmcpout = true;
    }
    
    if(!(btrig && btof1 && btof2 && bmcpout)) continue;

    tof1=0;
    tof2=0;
    Int_t counts(0);
    Double_t tot1(0),tot2(0);
    for(Int_t i=0; i<fEvent->GetHitSize(); i++){
      fHit = fEvent->GetHit(i);
      if(fHit.GetChannel()==960 ){
	tof1 = fHit.GetLeadTime();
	hTof1->Fill(fHit.GetLeadTime());
	tot1=fHit.GetTotTime();
	hTot->Fill(tot1);
      }
      if(fHit.GetChannel()==1104){
	tof2 = fHit.GetLeadTime();
	hTof2->Fill(fHit.GetLeadTime());
	//hTot->Fill(fHit.GetTotTime());	
	tot2=fHit.GetTotTime();
	hTot->Fill(tot2);
      }

      if(tof1!=0 && tof2!=0){
	Double_t time = (tof2-tof1);
	hTof->Fill(time);	
	//	if(insideOfEllipce(tof2-tof1, tot1, walk1x0, walky0,walktheta) && insideOfEllipce(tof2-tof1, tot2, walk1x0, walky1,-walktheta)){
	  hLeTotW->Fill(time,tot1);
	  time += (tot1-walky0)*tan(walktheta);
	  hLeTotC->Fill(time,tot1);
	  time += (tot2-walky1)*tan(-walktheta);
	  hLeTotC2->Fill(time,tot2);
	  hTofC->Fill(time);
	  //	}else if(insideOfEllipce(tof2-tof1, tot1, walk2x0, walky0,walktheta) && insideOfEllipce(tof2-tof1, tot2, walk2x0, walky1,-walktheta)){
	//   hLeTotW->Fill(time,tot1);
	//   time += (tot1-walky0)*tan(walktheta);
	//   hLeTotC->Fill(time,tot1);
	//   time += (tot2-walky1)*tan(-walktheta);
	//   hLeTotC2->Fill(time,tot2);
	//   hTofC->Fill(time);
	// }
	
	hLeTot->Fill(tof2-tof1,tot1);
      }
      if(fHit.GetLeadTime()<-270 || fHit.GetLeadTime() > -200)  continue;
      
      counts++;
    }

    hMult->Fill(counts);
  }

  canvasAdd("LeTot",800,500);
  hLeTot->Draw("colz");

  canvasAdd("LeTotW",800,500);
  hLeTotW->Draw("colz");
  
  canvasAdd("LeTotC",800,500);
  hLeTotC->Draw("colz");

  canvasAdd("LeTotC2",800,500);
  hLeTotC2->Draw("colz");
  
  canvasAdd("tof",800,500);
  prt_fit(hTof,3); 
  hTof->Draw();

  canvasAdd("tofC",800,500);
  prt_fit(hTofC,3);
  hTofC->Draw();

  // canvasAdd("tot",800,500);
  // hTot->Draw();
  
  //  canvasSave(1,0);
}
