#define prt__beam
#include "../prtdirc/src/PrtHit.h"
#include "../prtdirc/src/PrtEvent.h"
#include "prttools.C"
#include <TCanvas.h>

void drawTime(TString infile="~/sim4/d/proc/lab/pilas18124170655C.root"){

  TH1F * hLe  = new TH1F("le","LE;LE [ns];entries [#]",5000,-55,40);
  TH1F * hTot  = new TH1F("tot","TOT;TOT [ns];entries [#]",1000,0,100);
  TH2F *hLeTot = new TH2F("hLeTot","hLeTot;LE [ns];TOT [ns]",5000,-55,40,2000,75,85);

  TH2F *hEvtLe = new TH2F("hEvtLe","hEvtLe;event [#];LE [ns]",500,0,500,100,-0.5,0.5);
  // // TH1F * hLe  = new TH1F("le","LE;LE [ns];entries [#]",500,-600,600);
  // // TH1F * hTot  = new TH1F("tot","TOT;TOT [ns];entries [#]",500,60,75);
  // // TH2F *hLeTot = new TH2F("hLeTot","hLeTot;LE [ns];TOT [ns]",500,-51,-44,500,60,75);

  
  // TH1F * hLe  = new TH1F("le","LE;LE [ns];entries [#]",500,-5,-2);
  // TH1F * hTot  = new TH1F("tot","TOT;TOT [ns];entries [#]",500,0,100);
  // TH2F *hLeTot = new TH2F("hLeTot","hLeTot;LE [ns];TOT [ns]",500,-3,-2,500,80,100);

  TGraph *gr = new TGraph();
  TF1 *fsin = new TF1("fsin","0.15*sin(x)-0.12",0,500);
  
  if(!prt_init(infile,1,"data/drawTime")) return;
  PrtHit hit;
  for (int ievent=0; ievent< prt_entries && ievent<1000; ievent++){
    prt_nextEvent(ievent,10000);

    int trig=816;
    //int trig=29;
    double trig_time=0;
    for(int h=0; h<prt_event->GetHitSize(); h++){
      hit = prt_event->GetHit(h);
      int ch = hit.GetChannel();
      if(ch==trig) trig_time = hit.GetLeadTime();
      //if(ch==16) trig_time16 = hit.GetLeadTime();
    }

    if(trig_time==0) continue;
    for(int h=0; h<prt_event->GetHitSize(); h++){
      hit = prt_event->GetHit(h);
      int ch = hit.GetChannel();
      if(ch!=16) continue;

      double le = hit.GetLeadTime()-trig_time;
      double tot = hit.GetTotTime();

      gr->SetPoint(ievent,ievent,le);
      hEvtLe->Fill(ievent,le);
      hLe->Fill(le);
      hTot->Fill(tot);
      hLeTot->Fill(le,tot);
    }
  }

  prt_canvasAdd("hLe");
  hLe->Fit("gaus");
  hLe->Draw();
  prt_canvasAdd("hTot");
  hTot->Draw();
  prt_canvasAdd("hLeTot");
  hLeTot->Draw("colz");

  prt_canvasAdd("hEvtLe");
  //  gr->Draw("APL");
  hEvtLe->Draw("colz");
  //  fsin->Draw("same");

  //prt_canvasSave(1,0);
}

