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
  Int_t nfound = 1, peakSearch = 1;
  if(integral>5){ 
    
    if(peakSearch == 1){
      //gaust->SetParLimits(2,0.1,2);
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
    
      gaust->SetParameter(2,0.3);
      gaust->SetParameter(5,0.3);
    }

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

void procData(TString path="auto", TString infile="hits.root", Int_t studyid = 0, Int_t fileid=0, Double_t mom=0, Double_t angle=0, Double_t z =0, Double_t x= 0, Double_t step=0){
  fSavePath = path+Form("/%d/%d",studyid,fielid);
  PrtInit(infile,1);

  Double_t mult(0);
  TFile *file = new TFile(infile+".res.root","recreate");
  TTree *tree= new TTree("proc","proc");
  tree->Branch("studyid", &studyid,"studyid/I");
  tree->Branch("fileid", &fileid,"fileid/I");
  tree->Branch("mom", &mom,"mon/D");
  tree->Branch("angle", &angle,"angel/D");
  tree->Branch("z", &z,"z/D");
  tree->Branch("x", &x,"x/D");
  tree->Branch("step", &step,"step/D");
  tree->Branch("mult",&mult,"mult/D");
  
  TH1F * hMult  = new TH1F("mult","Mult;multiplicity [#];entries [#]",300,0,300);
  TH1F * hLe  = new TH1F("le","LE; LE [ns]; entries [#]",1000,-300,0);
  TH1F * hTot  = new TH1F("tot","TOT; TOT [#]; entries [#]",1000,20,80);
 
  gStyle->SetOptStat(1001111);

  Double_t tof1, tof2;
  PrtHit fHit;
  Int_t entries = fCh->GetEntries();
  for (Int_t ievent=0; ievent<entries; ievent++){
    PrtNextEvent(ievent,1000);

    Double_t time(0), tTime(0);
    for(Int_t i=0; i<fEvent->GetHitSize(); i++){
      fHit = fEvent->GetHit(i);
      if(fHit.GetChannel()==1344) tTime=fHit.GetLeadTime();
    }

    Int_t counts(0);
    Double_t tot(0);
    for(Int_t i=0; i<fEvent->GetHitSize(); i++){
      fHit = fEvent->GetHit(i);
      
      if(fHit.GetChannel()<960 ){
	time = fHit.GetLeadTime() - tTime;
	tot=fHit.GetTotTime();
	
	hLe->Fill(time);
	hTot->Fill(tot);

	if(tot>20 && tot<60 && time<-50 && time > -220)  counts++;
      } 
    }

    if(counts>5) hMult->Fill(counts);
  }

  canvasAdd("Le",800,500);
  hLe->Draw();
  
  canvasAdd("tot",800,500);
  //  fit(hTof,3); 
  hTot->Draw();

  canvasAdd("mult",800,500);
  mult = fit(hMult,30).X();

  hMult->Draw();

  tree->Fill();
  tree->Write();
  
  canvasSave(1,0);
}
