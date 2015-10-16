#define prt__beam
#include "../prtdirc/src/PrtHit.h"
#include "../prtdirc/src/PrtEvent.h"
#include "prttools.C"
#include <TSpectrum.h>

void procData(TString path="auto", TString infile="", Int_t studyId = 0, Int_t fileId=0, Double_t mom=0,Int_t radiatorId=0, Int_t lensId=0, Double_t angle=0, Double_t z =0, Double_t x= 0, Double_t xstep=0, Double_t ystep=0){

  if(infile=="") return;
  
  fSavePath = path+Form("/%d/%d",studyId,fileId);
  PrtInit(path+"/"+infile,1);

  Double_t mult(0);
  TFile *file = new TFile(path+"/"+infile+".res.root","recreate");
  TTree *tree= new TTree("proc","proc");
  tree->Branch("studyId", &studyId,"studyId/I");
  tree->Branch("fileId", &fileId,"fileId/I");
  tree->Branch("mom", &mom,"mon/D");
  tree->Branch("radiatorId", &radiatorId,"radiatorId/I");
  tree->Branch("lensId", &lensId,"lensId/I");
  tree->Branch("angle", &angle,"angel/D");
  tree->Branch("z", &z,"z/D");
  tree->Branch("x", &x,"x/D");
  tree->Branch("xstep", &xstep,"xstep/D");
  tree->Branch("ystep", &ystep,"ystep/D");
  tree->Branch("mult",&mult,"mult/D");
  
  TH1F * hMult  = new TH1F("mult","Mult;multiplicity [#];entries [#]",300,0,300);
  TH1F * hLe  = new TH1F("le","LE; LE [ns]; entries [#]",1000,-230,-100);
  TH1F * hTot  = new TH1F("tot","TOT; TOT [#]; entries [#]",1000,20,80);
 
  gStyle->SetOptStat(1001111);

  Double_t tof1, tof2;
  PrtHit fHit;
  Int_t entries = fCh->GetEntries();
  for (Int_t ievent=0; ievent<entries; ievent++){
    PrtNextEvent(ievent,1000);
    Int_t counts(0);
    Double_t tot(0);
    Bool_t btrig(false),bmcpout(false),btof1(false),btof2(false);

    Double_t time(0), tTime(0);
    for(Int_t i=0; i<fEvent->GetHitSize(); i++){
      fHit = fEvent->GetHit(i);
      if(fHit.GetChannel()==1344) tTime = fHit.GetLeadTime();
      if(fHit.GetChannel()==1344) btrig = true;
      if(fHit.GetChannel()==960)  btof1 = true;
      if(fHit.GetChannel()==1104) btof2 = true;
      if(fHit.GetChannel()==1248) bmcpout = true;
    }

    for(Int_t i=0; i<fEvent->GetHitSize(); i++){
      fHit = fEvent->GetHit(i);

      if(btrig && btof1 && btof2 && bmcpout){
	if(fHit.GetChannel()<960 ){

	  time = fHit.GetLeadTime() - tTime;
	  tot = fHit.GetTotTime();	
	  hLe->Fill(time);
	  hTot->Fill(tot);

	  if(tot>20 && tot<60 && time<-180 && time > -220)  {
	    Int_t mcpid = fHit.GetMcpId();
	    Int_t pixid = fHit.GetPixelId()-1;
	    fhDigi[mcpid]->Fill(pixid%8, pixid/8);
	    counts++;
	  }
	}
      }
    }

    if(counts>5) hMult->Fill(counts);
  }

  TString ext = Form("_%d_%d",studyId,fileId);
  
  canvasAdd("le"+ext,800,500);
  hLe->Draw();
  
  canvasAdd("tot"+ext,800,500);
  hTot->Draw();

  canvasAdd("mult"+ext,800,500);
  mult = prt_fit(hMult,30).X();
  hMult->Draw();
  
  drawDigi("m,p,v\n",2,-2,-2);
  cDigi->SetName("hits"+ext);
  canvasAdd(cDigi);  
  
  tree->Fill();
  tree->Write();
  
  canvasSave(1,0);
}
