#define prt__sim
#include "../src/PrtHit.h"
#include "../src/PrtEvent.h"
#include "../../prttools/prttools.C"

void createTimeMap(TString infile="../build/hits.root"){
  //infile="/u/rdzhigad/dirc/data/july14C.root";

  fSaveFlag = 2;
  fInfo = "drawScan.C \n";
  PrtInit(infile,1); //digi

  const int nmcp = 15, npix = 64;
  TH1F * hPTime[nmcp][npix];
  TH1F * hTime = new TH1F("time",";time, [ns];entries, [#]",500,0,150);
  for(Int_t m=0; m<nmcp; m++){
    for(Int_t p=0; p<npix; p++){
      hPTime[m][p]  = new TH1F(Form("hPTime_%d",m*100+p),Form("mcp %d, pixel %d",m, p),250,0,25); //800,1800
      axisTime800x500(hPTime[m][p]);
      hPTime[m][p]->SetStats(0);
      hPTime[m][p]->SetLineColor(1);
    }
  }
  
  PrtHit fHit;
  Int_t angle = 0;
  for (Int_t ievent=0; ievent<fCh->GetEntries(); ievent++){
    PrtNextEvent(ievent,1000);
    angle = fEvent->GetAngle() + 0.01;
    for(Int_t h=0; h<fEvent->GetHitSize(); h++){
      fHit = fEvent->GetHit(h);
      Int_t mcpid = fHit.GetMcpId();
      Int_t pixid = fHit.GetPixelId()-1;
      
      Double_t time = fHit.GetLeadTime();
      fhDigi[mcpid]->Fill(7-pixid/8, pixid%8);
      hPTime[mcpid][pixid]->Fill(time);
      hTime->Fill(time);
    }
  }

  TString path = createDir(Form("rdata/timemap/%d",angle), fInfo, fSaveFlag); 
  writeInfo("digi.csv", drawDigi("m,p,v\n",1), fSaveFlag);

  canvasAdd("time");
  hTime->Draw();
  canvasSave(2,"createTimeMap.C",1);

  // TCanvas* c1 = new TCanvas("c1","c1",0,0,800,400);
  // for(Int_t m=0; m<nmcp; m++){
  //   for(Int_t p=0; p<npix; p++){
  //     if(hPTime[m][p]->GetEntries()<1) continue;
  //     hPTime[m][p]->Draw();
  //     hPTime[m][p]->SetTitle(Form("Theta %d, mcp %d, pixel %d", angle, m, p));
  //     // c1->Modified(); c1->Update(); c1->WaitPrimitive();
  //     save(c1,path,Form("time_ang%dmcp%dpix%d",angle,m,p),fInfo,fSaveFlag,1);
  //   }
  // }

}
