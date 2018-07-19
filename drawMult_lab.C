#define prt__beam
#include "../prtdirc/src/PrtHit.h"
#include "../prtdirc/src/PrtEvent.h"
#include "prttools.C"

void drawMult_lab(TString infile="hits.root",double inten=0){

  if(!prt_init(infile,1,"data/drawMult")) return;

  TH1F * hMultL  = new TH1F("multL","multL",1000,0,5);
  TH1F * hMultR  = new TH1F("multR","multR",1000,0,5);
  TH1F * hTime  = new TH1F("time ","time",1000,-200,200);

  double un551(0),ln423(0),ua000(0),la000(0);
  double lpere(0);
  
  int count[5][12][64]={0};
  int mcp,pix,ch;
  double time;
  PrtHit hit;
  for(int ievent=0; ievent<prt_entries && ievent<10000; ievent++){
    prt_nextEvent(ievent,1000);
    
    int counts=0,lcounts=0,rcounts=0;
    for(int h=0; h<prt_event->GetHitSize(); h++){
      hit = prt_event->GetHit(h);
      mcp = hit.GetMcpId();
      pix=hit.GetPixelId()-1;
      time = hit.GetLeadTime();
      ch = map_mpc[mcp][pix];
      
      hTime->Fill(time);
      
      if(time<25 || time > 35)  continue;

      if(mcp==7 && pix%8<4) lcounts++;
      if(mcp==7) rcounts++;
      
      if(mcp==7 && pix%8<3)
	count[0][mcp][pix]++;
    }

    hMultL->Fill(lcounts/21.);
    hMultR->Fill(rcounts/21.);
  }
    
  for(int m=0; m<12; m++){
    for(int p=0; p<64; p++){
      prt_hdigi[m]->Fill(p%8,p/8,count[0][m][p]);          
    }
  }
  

  prt_drawDigi("m,p,v\n",2017,0,0);
  //prt_canvasAdd(prt_cdigi);
  
  prt_canvasAdd("Time",800,500);
  hTime->Draw();

  prt_canvasAdd("MultR",800,500);
  ua000 = hMultR->GetMean();// prt_fit(hMultR,5,20,3).X();
  hMultR->Draw();

  prt_canvasAdd("MultL",800,500);
  la000 = hMultL->GetMean();//prt_fit(hMultL,5,20,3).X();
  hMultL->Draw();
  
  prt_canvasSave();

  infile.ReplaceAll("C.root","R.root");
  TFile fc(infile,"recreate");
  TTree *tc = new TTree("lab","lab");
  tc->Branch("inten",&inten,"inten/D");
  tc->Branch("un551",&un551,"un551/D");
  tc->Branch("ln423",&ln423,"ln423/D");
  tc->Branch("ua000",&ua000,"ua000/D");
  tc->Branch("la000",&la000,"la000/D");
  tc->Branch("lpere",&lpere,"lpere/D");

  tc->Fill();
  tc->Write();
  fc.Write();
  fc.Close();

  
}
