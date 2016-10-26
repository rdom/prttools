#define prt__beam
#include "../prtdirc/src/PrtHit.h"
#include "../prtdirc/src/PrtEvent.h"
#include "prttools.C"

void drawMult(TString infile="hits.root"){
  fSavePath = "auto";
  PrtInit(infile,1);

  TH1F * hMult  = new TH1F("mult","mult",1000,-1000,1000);
  TH1F * hTime  = new TH1F("time ","time",1000,-1000,1000);

  Int_t mcp, count[]={0,0,0,0};
  PrtHit fHit;
  for (Int_t ievent=0; ievent<50000; ievent++){ //fCh->GetEntries()
    PrtNextEvent(ievent,1000);
    Int_t pid = prt_event->GetParticle();

    Int_t counts=0;
    for(Int_t i=0; i<prt_event->GetHitSize(); i++){
      fHit = prt_event->GetHit(i);
      mcp = fHit.GetMcpId();
      hTime->Fill(fHit.GetLeadTime());
      
      //if(fHit.GetLeadTime()<-270 || fHit.GetLeadTime() > -150)  continue;
      //counts++;
      if(mcp==3){
	if(pid==211)count[0]++;
	if(pid==2212)count[1]++;
      }
      if(mcp==5){
	if(pid==211)count[2]++;
	if(pid==2212)count[3]++;
      }	          
    }

    hMult->Fill(counts);
  }

  std::cout<<"MCP3 "<< count[0]<<"/"<<count[1]<<" = "<<  count[0]/(Double_t)count[1]<<std::endl;
  std::cout<<"MCP5 "<< count[2]<<"/"<<count[3]<<" = "<<  count[2]/(Double_t)count[3]<<std::endl;
  
  
  canvasAdd("Time",800,500);
  hTime->Draw();
  canvasAdd("Mult",800,500);
  hMult->Draw();
  
  canvasSave(1,0);
}
