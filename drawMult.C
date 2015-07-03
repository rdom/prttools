#define prt__sim
#include "../prtdirc/src/PrtHit.h"
#include "../prtdirc/src/PrtEvent.h"
#include "prttools.C"


void drawMult(TString infile="hits.root"){
  fSavePath = "auto";
  PrtInit(infile,1);

  TH1F * hMult  = new TH1F("mult","mult",1000,-1000,1000);
  TH1F * hTime  = new TH1F("time ","time",1000,-1000,1000);
  
  PrtHit fHit;
  for (Int_t ievent=0; ievent<10000; ievent++){ //fCh->GetEntries()
    PrtNextEvent(ievent,1000);

    Int_t counts=0;
    for(Int_t i=0; i<fEvent->GetHitSize(); i++){
      fHit = fEvent->GetHit(i);
      hTime->Fill(fHit.GetLeadTime());
      
      if(fHit.GetLeadTime()<-270 || fHit.GetLeadTime() > -150)  continue;
      counts++;
    }

    hMult->Fill(counts);
      
  }

  canvasAdd("Time",800,500);
  hTime->Draw();
  canvasAdd("Mult",800,500);
  hMult->Draw();
  
  canvasSave(1,0);
}
