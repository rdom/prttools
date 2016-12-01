#define prt__beam
#include "../prtdirc/src/PrtHit.h"
#include "../prtdirc/src/PrtEvent.h"
#include "prttools.C"
#include <TCanvas.h>


void drawVsEvents(TString infile="~/simo/217n/pdf/beam_163*C.root", Int_t nevents=100000){
  fSavePath = "data/draw_vs_events_s217";
  canvasSave(1,0);
  gStyle->SetOptFit(1);
  
  PrtInit(infile,1);

  Int_t ch, setid(0),tevents(0);;
  Double_t tof1, tof2;
  PrtHit hit;
  Int_t entries = fCh->GetEntries();
  TH2F *hLeTot = new TH2F("hLeTot", ";LE [ns];TOT [ns]", 200,10,30, 200,1,8);
  TH2F *hLeTotTof1 = new TH2F("hLeTotTof1", ";LE [ns];TOT [ns]", 600,10,30, 1000,44.01,46);
  TH1F *hLe[5];
  for(Int_t i=0; i<5; i++) {
    hLe[i] = new TH1F("hLe_"+prt_name[i], ";LE [ns];entries [#]", 300,10,20);
    hLe[i]->SetLineColor(prt_color[i]);
  }
 
  TCanvas *cExport = new TCanvas("cExport","cExport",0,0,800,400);
  
  for (Int_t ievent=0; ievent<entries; ievent++){
    PrtNextEvent(ievent,1000);
    Double_t le(0),tot(0),tof1le(0),tof1tot(0);;

    for(Int_t i=0; i<prt_event->GetHitSize(); i++){
      hit = prt_event->GetHit(i);
      ch  = hit.GetChannel();
      
      if(ch==720){
	tof1le = hit.GetLeadTime();
	tof1tot = hit.GetTotTime();
      }
    }

    
    for(Int_t i=0; i<prt_event->GetHitSize(); i++){
      hit = prt_event->GetHit(i);
      ch  = hit.GetChannel();
      
      if(ch==295){
	le = hit.GetLeadTime();
	tot = hit.GetTotTime();

	le -= (tof1tot-44.9)/8.4; //7.1;
 
	Int_t mcpid = hit.GetMcpId();
	Int_t pixid = hit.GetPixelId()-1;
	hLeTotTof1->Fill(le,tof1tot);
	hLeTot->Fill(le,tot);
	hLe[prt_pid]->Fill(le);
	fhDigi[mcpid]->Fill(pixid%8, pixid/8);
      }
    }

    if(tevents++>nevents){
      // drawDigi("m,p,v\n",7);
      // cDigi->Print(fSavePath+Form("/hits_%d.png",setid));
      // for(Int_t m=0; m<nmcp; m++) if(fhDigi[m]) fhDigi[m]->Reset();
      cExport->SetName(Form("letot_%d",setid));
      canvasAdd(cExport);
      hLeTot->Draw("colz");
      canvasSave(1,0);
      cExport->SetName(Form("le_%d",setid));
      canvasAdd(cExport);
      prt_fit(hLe[2],0.5);
      hLe[2]->Draw();
      hLe[4]->Draw("same");
      canvasSave(1,0);
      
      hLeTot->Reset();
      for(Int_t i=0; i<5; i++) hLe[i]->Reset();
      tevents=0;
      setid++;	
    }
  }
  canvasAdd("cLeTotTof1",800,400);
  TGraph *g = prt_fitslices(hLeTotTof1,12,13,0.3,10);
  g->SetLineColor(2);
  g->SetMarkerStyle(21);
  g->SetMarkerSize(0.6);
  g->SetMarkerColor(2);
  
  hLeTotTof1->GetXaxis()->SetRangeUser(11,16);
  hLeTotTof1->SetMaximum(150);
  hLeTotTof1->Draw("colz");
  // g->Fit("pol1");//,"","",12.6,12.8);
  g->Draw("P same");
  canvasSave(0,0);
}
