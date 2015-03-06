#define prt__sim
#include "../../prttools/prttools.C"

void drawScan(TString infile="../build/hits.root"){
  fSaveFlag = 2;
  fInfo = "drawScan.C \n";
  PrtInit(infile,1); //digi
  
  Int_t angle = 0;
  Int_t step = 0;
  for (Int_t ievent=0; ievent<fCh->GetEntries(); ievent++){
    PrtNextEvent(ievent,1000);
    if(ievent==0){
      angle = fEvent->GetAngle() + 0.01;
      step = fEvent->GetPrismStep();
      fInfo +=  fEvent->PrintInfo();
    }
    for(Int_t h=0; h<fEvent->GetHitSize(); h++){
      fHit = fEvent->GetHit(h);
      Int_t mcpid = fHit.GetMcpId();
      Int_t pixid = fHit.GetPixelId()-1;
      
      Double_t time = fHit.GetLeadTime();
      fhDigi[mcpid]->Fill(7-pixid/8, pixid%8);
    }
  }
  
  drawDigi("m,p,v\n",1);
  save(cDigi,"scan",Form("sc_%d_%d",angle,step),fInfo,fSaveFlag,1,-1);
}
