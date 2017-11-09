#define prt__beam
#include "../prtdirc/src/PrtHit.h"
#include "../prtdirc/src/PrtEvent.h"
#include "prttools.C"

void drawRatio(TString infile="../prtdirc/macro/hits.root"){
  
  if(!prt_init(infile,1,"data/drawRatio")) return;
 
  Double_t mult[5][prt_maxdircch];
  PrtHit hit;
  for (Int_t ievent=0; ievent< prt_entries; ievent++){
    prt_nextEvent(ievent,1000);
    for(Int_t h=0; h<prt_event->GetHitSize(); h++){
      hit = prt_event->GetHit(h);
      Int_t mcpid = hit.GetMcpId();
      Int_t pixid = hit.GetPixelId()-1;
      Double_t time = hit.GetLeadTime();
      Int_t ch = map_mpc[mcpid][pixid];
      if(prt_isBadChannel(ch)) continue;
      mult[prt_pid][ch]++;
    }
  }


  for(auto i=0; i<prt_maxdircch; i++){
    Int_t mcpid = map_mcp[i];
    Int_t pixid = map_pix[i];
    if(mult[2][i]>10)
      prt_hdigi[mcpid]->Fill(pixid%8, pixid/8,mult[2][i]/(double)mult[4][i]);
  }


  //i%3*4+i/3
 
  prt_drawDigi("m,p,v\n",2017,4,0);
  prt_cdigi->SetName(Form("hp_ratio_dat_%d_%2.1f",(Int_t)prt_theta,prt_phi));
  prt_canvasAdd(prt_cdigi);
  prt_canvasSave(1,0);
}
