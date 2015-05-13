#define prt__sim
#include "../prtdirc/src/PrtHit.h"
#include "../prtdirc/src/PrtEvent.h" 
#include "prttools.C"

void drawTimeDiff(TString infile="fileC0.root", Int_t triggerCh = 1776){
  TH1F* hTimeDiff = new TH1F("hTimeDiff","hTimeDiff;time [ns];entries [#]",100,-200,200);
  PrtInit(infile,0);
  Double_t triggerTime;
  PrtHit fHit;

  for (Int_t ievent=0; ievent<fCh->GetEntries(); ievent++){
    PrtNextEvent(ievent,1000);
    triggerTime=  0;

    Int_t nhits = fEvent->GetHitSize();
    for(Int_t h=0; h<fEvent->GetHitSize(); h++){
      fHit = fEvent->GetHit(h);
      if(triggerCh==fHit.GetChannel()) {
	triggerTime = fHit.GetLeadTime();
	break;
      }
    }

    for(Int_t h=0; h<fEvent->GetHitSize(); h++){
      fHit = fEvent->GetHit(h);
      Int_t ch = fHit.GetChannel();
      hTimeDiff->Fill(fHit.GetLeadTime() - triggerTime);
    }
  }
 
  hTimeDiff->Draw();
}
