#define prt__beam
#include "../prtdirc/src/PrtHit.h"
#include "../prtdirc/src/PrtEvent.h"
#include "prttools.C"

void drawMult_lab(TString infile="hits.root",double inten=0){

  TObjArray *sarr = infile.Tokenize("_");
  if(sarr->GetEntries()==3){
    TString soffset = ((TObjString *) sarr->At(1))->GetName();
    inten = soffset.Atof();
  }
  std::cout<<"inten "<<inten<<std::endl;
  
  if(!prt_init(infile,1,"data/drawMult_lab")) return;

  TH1F * hMult1  = new TH1F("mult1","mult1",50,0,50);
  TH1F * hMult2  = new TH1F("mult2","mult2",50,0,50);
  TH1F * hMult3  = new TH1F("mult3","mult3",50,0,50);
  TH1F * hMult4  = new TH1F("mult4","mult4",50,0,50);
  TH1F * hTime  = new TH1F("time ","time",1000,-200,200);

  double count1(0),count2(0), count3(0),count4(0);
  double lpere(0);
  
  int count[5][12][64]={0};
  int mcp,pix,ch;
  double time;
  PrtHit hit;
  for(int ievent=0; ievent<prt_entries && ievent<100000; ievent++){
    prt_nextEvent(ievent,1000);
    
    int counts=0,mcounts=0,ucounts=0,mcounts_t=0,ucounts_t=0;
    for(int h=0; h<prt_event->GetHitSize(); h++){
      hit = prt_event->GetHit(h);
      mcp = hit.GetMcpId();
      pix=hit.GetPixelId()-1;
      time = hit.GetLeadTime();
      ch = map_mpc[mcp][pix];
      
      hTime->Fill(time);

      if(mcp==0 && ch<31) ucounts++;
      else mcounts++;
      
      if(time<25 || time > 35)  continue;      

      if(mcp==0 && ch<31) ucounts_t++;
      else mcounts_t++;
    }
    hMult1->Fill(ucounts);
    hMult2->Fill(mcounts);
    
    hMult3->Fill(ucounts_t);
    hMult4->Fill(mcounts_t);
  }
    
  // for(int m=0; m<12; m++){
  //   for(int p=0; p<64; p++){
  //     prt_hdigi[m]->Fill(p%8,p/8,count[0][m][p]);          
  //   }
  // }
  

  prt_drawDigi("m,p,v\n",2018,0,0);
  //prt_canvasAdd(prt_cdigi);
  
  prt_canvasAdd("Time",800,500);
  hTime->Draw();

  prt_canvasAdd("Mult",800,500);
  count1 = hMult1->GetMean();// prt_fit(hMultR,5,20,3).X();
  count2 = hMult2->GetMean();
  count3 = hMult3->GetMean();
  count4 = hMult4->GetMean();
  hMult1->Draw();
  hMult2->SetLineColor(kRed);
  hMult2->Draw("same");
  hMult3->SetLineColor(kGreen);
  hMult3->Draw("same");
  hMult4->SetLineColor(kBlue);
  hMult4->Draw("same");
    
  prt_canvasSave();

  infile.ReplaceAll("C.root","R.root");
  TFile fc(infile,"recreate");
  TTree *tc = new TTree("lab","lab");
  tc->Branch("inten",&inten,"inten/D");
  tc->Branch("count1",&count1,"count1/D");
  tc->Branch("count2",&count2,"count2/D");
  tc->Branch("count3",&count3,"count3/D");
  tc->Branch("count4",&count4,"count4/D");
  int res=75;
  tc->Branch("res",&res,"res/D");
  
  tc->Fill();
  tc->Write();
  fc.Write();
  fc.Close();
}
