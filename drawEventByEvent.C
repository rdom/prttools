#define prt__beam
#include "../prtdirc/src/PrtHit.h"
#include "../prtdirc/src/PrtEvent.h"
#include "prttools.C"
#include <TCanvas.h>


void drawEventByEvent(TString infile="~/simo/laser/all_0910_2/th_800_16283152516C.root", Int_t nevents=100, bool mc=true){
  fSavePath = "data/draw_event_by_event";
  canvasSave(1,0);
  
  PrtInit(infile,1);

  Int_t ch, counts(0);
  Double_t tof1, tof2;
  PrtHit hit;
  Int_t entries = fCh->GetEntries();
  for (Int_t ievent=0; ievent<entries; ievent++){
    if(counts>nevents) continue;
    PrtNextEvent(ievent,1000);
    bool newhit(false);
    Double_t tot(0);
    Bool_t btrig(false),bmcpout(false),btof1(false),btof2(false);

    Double_t time(0);
    for(Int_t i=0; i<prt_event->GetHitSize(); i++){
      hit = prt_event->GetHit(i);

      if(mc || (btrig && btof1 && btof2 && bmcpout && hit.GetChannel()<960)){
	  time = hit.GetLeadTime();
	  tot = hit.GetTotTime();
	  ch  = hit.GetChannel();
	  if(time<9 || time >30) continue;
	  if(ch==-1) ch = map_mpc[mcp][ hit.GetPixelId()-1];
	  if(badcannel(ch) || time <0) continue;
	  if( mc || (tot>20 && tot<60 && time<-180 && time > -220)){
	    Int_t mcpid = hit.GetMcpId();
	    Int_t pixid = hit.GetPixelId()-1;
	    
	    if(pixid==0) std::cout<<"time  "<<time <<std::endl;
	    
	    fhDigi[mcpid]->Fill(pixid%8, pixid/8);
	    newhit = true;	    
	   }
      }
    }
    if(newhit){
      if(counts%1==0){
	drawDigi("m,p,v\n",7);
	cDigi->Print(fSavePath+Form("/hits_%d.png",counts));
	for(Int_t m=0; m<nmcp; m++) if(fhDigi[m]) fhDigi[m]->Reset();
      }
      counts++;
    }
    
  }

}
