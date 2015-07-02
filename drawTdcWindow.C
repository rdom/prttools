#define prt__sim
#include "../prtdirc/src/PrtHit.h"
#include "../prtdirc/src/PrtEvent.h"
#include "prttools.C"


void drawTdcWindow(TString infile="hits.root"){
  //fSavePath = "scan3";
  PrtInit(infile,1);

  TH1F * hTime[tdcnum];
  for(Int_t i=0; i<tdcnum; i++) hTime[i]   = new TH1F(Form("time%d",i),Form("time%d",i),1000,-1000,1000);
  
  PrtHit fHit;
  for (Int_t ievent=0; ievent<fCh->GetEntries(); ievent++){
    PrtNextEvent(ievent,1000);
    for(Int_t h=0; h<fEvent->GetHitSize(); h++){
      fHit = fEvent->GetHit(h);
      hTime[fHit.GetTdc()]->Fill(fHit.GetLeadTime());
      hTime[fHit.GetTdc()]->Fill(fHit.GetLeadTime()+fHit.GetTotTime());
    }
  }
  TFile *f = new TFile(infile+".time.root","recreate");
  for(Int_t i=0; i<tdcnum; i++) hTime[i]->Write();
  
  //canvasSave(1,0);
}
