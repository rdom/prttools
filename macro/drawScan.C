#define prt__sim
#include "../src/PrtHit.h"
#include "../src/PrtEvent.h"
#include "../../prttools/prttools.C"

void drawScan(TString infile="../../data/beam_15183032354C0.root"){
  fSavePath = "scan3";
  PrtInit(infile,1); //digi
  
  Int_t itest(0);
  PrtHit fHit;
  for (Int_t ievent=0; ievent<100000; ievent++){
    PrtNextEvent(ievent,1000);
    for(Int_t h=0; h<fEvent->GetHitSize(); h++){
      fHit = fEvent->GetHit(h);
      Int_t mcpid = fHit.GetMcpId();
      Int_t pixid = fHit.GetPixelId()-1;
      
      Double_t time = fHit.GetLeadTime();
      //      fhDigi[mcpid]->Fill(7-pixid/8, pixid%8);
      fhDigi[mcpid]->Fill(pixid%8, pixid/8);
    }
  }
  itest = fTest1+50;
  drawDigi("m,p,v\n",2,-2,-2);
  cDigi->SetName(Form("sc_%d_%d",fAngle,fMomentum/1000));
  canvasAdd(cDigi);  
  canvasSave(1,0);
  
}
