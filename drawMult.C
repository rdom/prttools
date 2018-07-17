#define prt__beam
#include "../prtdirc/src/PrtHit.h"
#include "../prtdirc/src/PrtEvent.h"
#include "prttools.C"

void drawMult(TString infile="hits.root"){

  if(!prt_init(infile,1,"data/drawMult")) return;

  TH1F * hMult  = new TH1F("mult","mult",100,0,100);
  TH1F * hTime  = new TH1F("time ","time",1000,-200,200);

  int count[5][9][64]={0};

  int mcp,pix,ch;
  double time;
  PrtHit hit;
  for(int ievent=0; ievent<prt_entries && ievent<10000; ievent++){
    prt_nextEvent(ievent,1000);
    
    int counts=0;
    for(int h=0; h<prt_event->GetHitSize(); h++){
      hit = prt_event->GetHit(h);
      mcp = hit.GetMcpId();
      pix=hit.GetPixelId()-1;
      time = hit.GetLeadTime();
      ch = map_mpc[mcp][pix];
      
      hTime->Fill(time);
      
      if(time<20 || time > 40)  continue;
      //      if(ch==551) counts++;
      if(ch==423) counts++;
      count[prt_pid][mcp][pix]++;
    }

    hMult->Fill(counts);
  }
    
  for (int m=0; m <prt_nmcp; m++) {
    Double_t t1(0),t2(0);
    for(Int_t p=0; p<prt_npix; p++){
      t1+=count[4][m][p];
      t2+=count[2][m][p];
      prt_hdigi[m]->Fill(p%8,p/8,count[0][m][p]);          
    }
    std::cout<<"t1   "<<t1 << "   t2 "<<t2<<std::endl;
    
    // for(Int_t p=0; p<npix; p++){
    //   fhDigi[m]->Fill(p%8,p/8,t1/t2);
    // }
  }

  

  prt_drawDigi("m,p,v\n",2017,0,0);

  prt_canvasAdd("Time",800,500);
  hTime->Draw();
  prt_canvasAdd("Mult",800,500);
  hMult->Draw();
  
  prt_canvasSave();
}
