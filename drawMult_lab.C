#define prt__beam
#include "../prtdirc/src/PrtHit.h"
#include "../prtdirc/src/PrtEvent.h"
#include "prttools.C"

void drawMult_lab(TString infile="hits.root",double inten=0){

  if(!prt_init(infile,1,"data/drawMult")) return;

  TH1F * hMult  = new TH1F("mult","mult",100,0,100);
  TH1F * hTime  = new TH1F("time ","time",1000,-200,200);

  int un551(0),ln423(0),ua000(0),la000(0);
  double lpere(0);
  
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
      if(ch==551) un551++;
      if(mcp==8) ua000++;
      
      if(ch==423) ln423++;
      if(mcp==6) la000++;

      if(mcp==6) counts++;
      count[prt_pid][mcp][pix]++;
    }

    hMult->Fill(counts);
  }
    
  // for (int m=0; m <prt_nmcp; m++) {
  //   for(Int_t p=0; p<prt_npix; p++){
  //     prt_hdigi[m]->Fill(p%8,p/8,count[0][m][p]);          
  //   }
  // }
  

  prt_drawDigi("m,p,v\n",2017,0,0);
  
  prt_canvasAdd("Time",800,500);
  lpere = prt_fit(hMult,5,20,3).X();
 
  hTime->Draw();
  prt_canvasAdd("Mult",800,500);
  hMult->Draw();
  
  prt_canvasSave();

  infile.ReplaceAll("*C.root","R.root");
  TFile fc(infile,"recreate");
  TTree *tc = new TTree("lab","lab");
  tc->Branch("inten",&inten,"inten/D");
  tc->Branch("un551",&un551,"un551/I");
  tc->Branch("ln423",&ln423,"ln423/I");
  tc->Branch("ua000",&ua000,"ua000/I");
  tc->Branch("la000",&la000,"la000/I");
  tc->Branch("lpere",&lpere,"lpere/D");

  tc->Fill();
  tc->Write();
  fc.Write();
  fc.Close();

  
}
