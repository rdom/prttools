#define prt__beam
#include "../prtdirc/src/PrtHit.h"
#include "../prtdirc/src/PrtEvent.h"
#include "prttools.C"

#include <TSpectrum.h>

TSpectrum *spect = new TSpectrum(2);
TF1 *gaust;
TVector3 fit(TH1F *h, Double_t range = 3){
  int binmax = h->GetMaximumBin();
  double xmax = h->GetXaxis()->GetBinCenter(binmax);
  gaust = new TF1("gaust","gaus(0)",xmax-range,xmax+range);
  gaust->SetLineColor(1);
  Double_t integral = h->Integral(h->GetXaxis()->FindBin(xmax-0.6),h->GetXaxis()->FindBin(xmax+0.6));
  Double_t xxmin, xxmax, sigma1=0, mean1=0, sigma2, mean2;
  xxmax = xmax;
  xxmin = xxmax;
  Int_t nfound = 1, peakSearch = 2;
  if(integral>5){ 
    
    if(peakSearch == 1){
      gaust->SetParLimits(2,0.08,2);
      gaust->SetParameter(1,xmax);
      gaust->SetParameter(2,0.2);
    }
    
    if(peakSearch == 2){
      nfound = spect->Search(h,2,"",0.2);
      std::cout<<"nfound  "<<nfound <<std::endl;
      Float_t *xpeaks = spect->GetPositionX();
      if(nfound==1){
	gaust =new TF1("gaust","gaus(0)",xmax-range,xmax+range);
	gaust->SetParameter(1,xpeaks[0]);
      }else if(nfound==2) {
	if(xpeaks[0]>xpeaks[1]) {
	  xxmax = xpeaks[0];
	  xxmin = xpeaks[1];
	}else {
	  xxmax = xpeaks[1];
	  xxmin = xpeaks[0];
	}
	gaust =new TF1("gaust","gaus(0)+gaus(3)",xmax-range,xmax+range);
	gaust->SetParameter(1,xxmin);
	gaust->SetParameter(4,xxmax);
      }
    
      gaust->SetParameter(2,0.15);
      gaust->SetParameter(5,0.15);
    }

    h->Fit("gaust","","MQN",xxmin-range, xxmax+range);
    h->Fit("gaust","","MQN",xxmin-range, xxmax+range);
    mean1 = gaust->GetParameter(1);
    sigma1 = gaust->GetParameter(2);
    if(sigma1>10) sigma1=10;
    
    if(peakSearch == 2){ 
      mean2 = (nfound==1) ? gaust->GetParameter(1) : gaust->GetParameter(4);
      sigma2 = (nfound==1) ? gaust->GetParameter(2) : gaust->GetParameter(5);
    }
  }
  delete gaust;
  return TVector3(mean1,sigma1,0);
}
Double_t walktheta(-16.5*TMath::Pi()/180.);
Double_t deltatheta(2.5*TMath::Pi()/180.);
Double_t walky0(42.23), walky1(44.00);
Double_t walk1x0(174.85), walk2x0(175.64);

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

      if(fHit.GetChannel()==1104){ //1104
	tof2 = fHit.GetLeadTime();
	hTof2->Fill(fHit.GetLeadTime());
	//hTot->Fill(fHit.GetTotTime());	
	tot2=fHit.GetTotTime();
	hTot->Fill(tot2);
      }
      if(tof1!=0 && tof2!=0) {
	std::cout<<"tof2-tof1 "<<tof2-tof1 <<std::endl;
	
	hTof->Fill(tof2-tof1);
	if(insideOfEllipce(tof2-tof1, tot1, walk1x0, walky0,walktheta) && insideOfEllipce(tof2-tof1, tot2, walk1x0, walky1,-walktheta)){
	  Double_t time = (tof2-tof1);
	  hLeTotW->Fill(time,tot1);
	  time += (tot1-walky0)*tan(walktheta+deltatheta);
	  hLeTotC->Fill(time,tot1);
	  time += (tot2-walky1)*tan(-walktheta-deltatheta);
	  hLeTotC2->Fill(time,tot2);
	  hTofC->Fill(time);
	}
	if(insideOfEllipce(tof2-tof1, tot1, walk2x0, walky0,walktheta) && insideOfEllipce(tof2-tof1, tot2, walk2x0, walky1,-walktheta)){
	  Double_t time = (tof2-tof1);
	  hLeTotW->Fill(time,tot1);
	  time += (tot1-walky0)*tan(walktheta+deltatheta);
	  hLeTotC->Fill(time,tot1);
	  time += (tot2-walky1)*tan(-walktheta-deltatheta);
	  hLeTotC2->Fill(time,tot2);
	  hTofC->Fill(time);
	}
	
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
  fit(hTof,3); 
  hTof->Draw();

  canvasAdd("tofC",800,500);
  fit(hTofC,3);
  hTofC->Draw();

  // canvasAdd("tot",800,500);
  // hTot->Draw();
  
  //  canvasSave(1,0);
}
