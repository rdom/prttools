#define prt__beam
#include "../prtdirc/src/PrtHit.h"
#include "../prtdirc/src/PrtEvent.h"
#include "prttools.C"

void procData(TString path="/data.local/data/jun15", TString infile="", Int_t studyId = 0, Int_t fileId=0, Double_t mom=0,Int_t radiatorId=0, Int_t lensId=0, Double_t angle=0, Double_t z=0, Double_t x=0, Double_t xstep=0, Double_t ystep=0){
  
  if(infile=="") return;

  Double_t mult(0),le1(0),le2(50),offset(0);
  fSavePath = path+Form("/%ds/%d",studyId,fileId);
  
  if(infile.Contains("C.root")) { // beam data
    // offset=284.59;
    // le1=280;
    // le2=330;
    fSavePath = path+Form("/%dr/%d",studyId,fileId);
  }
  
  PrtInit(path+"/"+infile,1);
  CreateMap();
  
  TFile *file = new TFile(path+"/"+infile+".res.root","recreate");
  TTree *tree= new TTree("proc","proc");
  tree->Branch("studyId", &studyId,"studyId/I");
  tree->Branch("fileId", &fileId,"fileId/I");
  tree->Branch("mom", &mom,"mon/D");
  tree->Branch("radiatorId", &radiatorId,"radiatorId/I");
  tree->Branch("lensId", &lensId,"lensId/I");
  tree->Branch("angle", &angle,"angel/D");
  tree->Branch("z", &z,"z/D");
  tree->Branch("x", &x,"x/D");
  tree->Branch("xstep", &xstep,"xstep/D");
  tree->Branch("ystep", &ystep,"ystep/D");
  tree->Branch("mult",&mult,"mult/D");
  
  TH1F * hMult  = new TH1F("mult","Mult;multiplicity [#];entries [#]",300,0,300);
  TH1F * hLeA  = new TH1F("le","LE; LE [ns]; entries [#]",1000,le1,le2);
  TH1F * hTot  = new TH1F("tot","TOT; TOT [#]; entries [#]",1000,-2,15);

  Int_t colors[4]={2,4,2,4}
  TH1F * hLe1[maxch_dirc][4];
  for(Int_t t=0; t<4; t++){
    for(Int_t i=0; i<maxch_dirc; i++){
      hLe[i][t] = new TH1F(Form("le_ch_%d_%d",t,i),"LE; LE [ns]; entries [#]",1000,le1,le2);
      hLe[i][t]->SetLineColor(colors[t]);
    }    
  }
  
  gStyle->SetOptStat(1001111);
  gStyle->SetOptFit(1111);
 
  PrtHit fHit;
  Int_t entries = fCh->GetEntries();
  for (Int_t ievent=0; ievent<entries; ievent++){
    PrtNextEvent(ievent,1000);
    Int_t counts(0),pid(0);
    Double_t tot(0),time(0);
    if(prt_event->GetParticle()==211) pid=1;
    
    for(Int_t i=0; i<prt_event->GetHitSize(); i++){
      fHit = prt_event->GetHit(i);
      Int_t ch = fHit.GetChannel();
      if(ch==-1) ch = map_mpc[fHit.GetMcpId()][fHit.GetPixelId()-1];      
      
      if(ch<maxch_dirc && !badcannel(ch)){
	time = fHit.GetLeadTime()-offset;
	tot = fHit.GetTotTime();
	hLe[ch][pid]->Fill(time);
	hLeA->Fill(time);
	hTot->Fill(tot);

	if(time<le2 && time>le1){
	  Int_t mcpid = fHit.GetMcpId();
	  Int_t pixid = fHit.GetPixelId()-1;
	  fhDigi[mcpid]->Fill(pixid%8, pixid/8);
	  counts++;
	}
      }
    }

    if(counts>5) hMult->Fill(counts);
  }

  TString ext = Form("_%d_%d",studyId,fileId);
  TCanvas *cExport = new TCanvas("cExport","cExport",0,0,800,400);
  cExport->SetCanvasSize(800,400);
  for(Int_t i=0; i<maxch_dirc; i++){
    hLe[i][0]->Draw();
    hLe[i][1]->Draw("same");
    cExport->SetName(Form("hLe_%d",i));
    canvasAdd(cExport);
    canvasSave(1,0);
  }  
  
  
  canvasAdd("p_le"+ext,800,400);
  prt_fit(hLeA,0.3,100,100).X();
  hLeA->Draw();
  
  canvasAdd("p_tot"+ext,800,400);
  hTot->Draw();

  canvasAdd("p_mult"+ext,800,400);
  mult = prt_fit(hMult,20,20,100).X();
  hMult->Draw();
  
  drawDigi("m,p,v\n",7,-2,-2);
  cDigi->SetName("p_hits"+ext);
  canvasAdd(cDigi);  
  
  tree->Fill();
  tree->Write();
  
  canvasSave(1,0);
}
