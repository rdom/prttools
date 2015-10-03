#define prt__sim
#include "../src/PrtHit.h"
#include "../src/PrtEvent.h"
#include "../../prttools/prttools.C"

const Double_t prismangle = 45;
Double_t prismSize[] = {50+300*tan(prismangle*TMath::Pi()/180.),170};
Double_t prismShift = (prismSize[0])/2. -50/2.;
void drawPrism(Double_t x, Double_t y){
  gPad->cd();
  TBox *box2 = new TBox(x-prismSize[0]/2.,y-prismSize[1]/2.,x+prismSize[0]/2.,y+prismSize[1]/2.);
  box2->SetFillStyle(0);
  box2->SetLineColor(4);
  box2->SetLineWidth(2);
  box2->Draw("same");
}

void drawLoad(TString infile="../build/hits.root"){
  gStyle->SetOptStat(0);
  fInfo = "drawScan.C \n";
  PrtInit(infile,0); //digi
  
  TH2F* hHits = new TH2F("hHits",";x, [mm];y, [mm]",500,-40,350,500,-100,100);
  Int_t angle(0), step(0);
  Double_t test(0);
  PrtHit fHit;
  for (Int_t ievent=0; ievent<fCh->GetEntries(); ievent++){
    PrtNextEvent(ievent,1000);
    if(ievent==0){
      angle = fEvent->GetAngle() + 0.01;
      step = fEvent->GetPrismStep()*2;
      test = fEvent->GetTest1();
      fInfo +=  fEvent->PrintInfo();
    }
    for(Int_t h=0; h<fEvent->GetHitSize(); h++){
      fHit = fEvent->GetHit(h);
      Int_t mcpid = fHit.GetMcpId();
      Int_t pixid = fHit.GetPixelId()-1;
      TVector3 pos = fHit.GetGlobalPos();
      
      Double_t time = fHit.GetLeadTime();
      hHits->Fill(pos.X(),pos.Y());
    }
  }
  
  canvasAdd(Form("load_%d",angle),800,500);
  hHits->SetStats(0);
  hHits->GetXaxis()->SetTitleOffset(0.85);
  hHits->GetYaxis()->SetTitleOffset(0.85);
  hHits->GetXaxis()->SetTitleSize(0.05);
  hHits->GetYaxis()->SetTitleSize(0.05);
  //hHits->SetTitle(Form("#theta_{track} = %d#circ",angle));
  hHits->Draw("colz");
  drawPrism(prismShift,0);
  canvasSave(0,"drawLoad.C",1,"data/load");
}
