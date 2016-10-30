#define prt__beam
#include "../prtdirc/src/PrtHit.h"
#include "../prtdirc/src/PrtEvent.h"
#include "prttools.C"

void drawMult(TString infile="hits.root"){
  fSavePath = "auto";
  PrtInit(infile,1);
  CreateMap();

  TH1F * hMult  = new TH1F("mult","mult",1000,-1000,1000);
  TH1F * hTime  = new TH1F("time ","time",1000,-1000,1000);

  Int_t count[5][9][64];

  for(Int_t i=0; i<5; i++){
    for (Int_t m=0; m <nmcp; m++) {
      for(Int_t p=0; p<npix; p++){
	count[i][m][p]=0;
      }
    }    
  }

  Int_t pid,mcp,pix;
  Double_t time;
  PrtHit fHit;
  for (Int_t ievent=0; ievent<fCh->GetEntries(); ievent++){ //fCh->GetEntries()
    PrtNextEvent(ievent,1000);
    pid = prt_event->GetParticle();

    Int_t counts=0;
    for(Int_t i=0; i<prt_event->GetHitSize(); i++){
      fHit = prt_event->GetHit(i);
      mcp = fHit.GetMcpId();
      pix=fHit.GetPixelId()-1;
      time = fHit.GetLeadTime();
      hTime->Fill(time);
      
      if(time<0 || time > 40)  continue;
      counts++;
      
      count[prt_pid][mcp][pix]++;
    }

    hMult->Fill(counts);
  }
    
  for (Int_t m=0; m <nmcp; m++) {
    Double_t t1(0),t2(0);
    for(Int_t p=0; p<npix; p++){
      t1+=count[4][m][p];
      t2+=count[2][m][p];
      fhDigi[m]->Fill(p%8,p/8,count[4][m][p]/(Double_t)count[2][m][p]);          
    }
    std::cout<<"t1   "<<t1 << "   t2 "<<t2<<std::endl;
    
    // for(Int_t p=0; p<npix; p++){
    //   fhDigi[m]->Fill(p%8,p/8,t1/t2);
    // }
  }

  

  drawDigi("m,p,v\n",7,1.6,0);
  cDigi->cd();
  (new TPaletteAxis(0.90,0.1,0.94,0.90,fhDigi[0]))->Draw();
  canvasAdd(cDigi);

  canvasAdd("Time",800,500);
  hTime->Draw();
  canvasAdd("Mult",800,500);
  hMult->Draw();
  
  canvasSave(1,0);
}
