#define prt__beam
#include "../prtdirc/src/PrtHit.h"
#include "../prtdirc/src/PrtEvent.h"
#include "prttools.C"
#include <TCanvas.h>


void drawEventByEvent(TString infile="~/simo/laser/all_0910_2/th_800_16283152516C.root", Int_t nevents=10, bool mc=0){

  if(!prt_init(infile,1,"data/draw_event_by_event_pilas")) return;
  prt_canvasSave();

  Int_t ch, counts(0);
  Double_t tof1, tof2;
  PrtHit hit;
  for (auto ievent=0; ievent< prt_entries; ievent++){
    if(counts>nevents) continue;
    prt_nextEvent(ievent,1000);

    bool newhit(false);
    Double_t tot(0);
    Bool_t btrig(false),bmcpout(false),btof1(false),btof2(false);

    Double_t time(0);
    for(auto i=0; i<prt_event->GetHitSize(); i++){
      hit = prt_event->GetHit(i);

      // if(mc || (btrig && btof1 && btof2 && bmcpout && hit.GetChannel()<960))
	{
	  time = hit.GetLeadTime();
	  tot = hit.GetTotTime();
	  ch  = hit.GetChannel();
	  if(time<20 || time >40) continue;
	  Int_t mcpid = hit.GetMcpId();
	  Int_t pixid = hit.GetPixelId()-1;
	  
	  if(ch==-1) ch = map_mpc[mcpid][ hit.GetPixelId()-1];
	  if(prt_isBadChannel(ch) || time <0) continue;
	  //	  if( mc || (tot>20 && tot<60 && time<-180 && time > -220))
	    {
	    if(pixid==0) std::cout<<"time  "<<time <<std::endl;
	    
	    prt_hdigi[mcpid]->Fill(pixid%8, pixid/8);
	    newhit = true;	    
	   }
      }
    }
    if(newhit){
      if(counts%1==0){
	prt_drawDigi("m,p,v\n",2017);
	prt_cdigi->Print(prt_savepath+Form("/hits_%d.png",counts));
	for(Int_t m=0; m<prt_nmcp; m++) if(prt_hdigi[m]) prt_hdigi[m]->Reset();
      }
      counts++;
    }
    
  }

}
