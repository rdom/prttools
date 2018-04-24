#define prt__beam
#include "../prtdirc/src/PrtHit.h"
#include "../prtdirc/src/PrtEvent.h"
#include "prttools.C"
#include <TCanvas.h>

void drawTime(TString infile="~/sim4/d/proc/lab/pilas18110155722_100pFC.root"){

  TH1F * hLe  = new TH1F("le","LE;LE [ns];entries [#]",500,-51,-44);
  TH1F * hTot  = new TH1F("tot","TOT;TOT [ns];entries [#]",500,60,75);
  TH2F *hLeTot = new TH2F("hLeTot","hLeTot;LE [ns];TOT [ns]",500,-51,-44,500,60,75);

  
  if(!prt_init(infile,1,"data/drawTime")) return;
  PrtHit hit;
  for (int ievent=0; ievent< prt_entries; ievent++){
    prt_nextEvent(ievent,10000);
    for(int h=0; h<prt_event->GetHitSize(); h++){
      hit = prt_event->GetHit(h);
      int mcpid = hit.GetMcpId();
      int pixid = hit.GetPixelId()-1;
      double le = hit.GetLeadTime();
      double tot = hit.GetTotTime();
      int ch = map_mpc[mcpid][pixid];
      if(prt_isBadChannel(ch)) continue;

      hLe->Fill(le);
      hTot->Fill(tot);
      hLeTot->Fill(le,tot);
    }
  }

  prt_canvasAdd("hLe");
  hLe->Draw();
  prt_canvasAdd("hTot");
  hTot->Draw();
  prt_canvasAdd("hLeTot");
  hLeTot->Draw("colz");

  //prt_canvasSave(1,0);
}

