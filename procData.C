#define prt__beam
#include "../prtdirc/src/PrtHit.h"
#include "../prtdirc/src/PrtEvent.h"
#include "datainfo.C"
#include "prttools.C"
#include "TLatex.h"


void procData(TString infile="", Int_t studyId = 0, Int_t fileId=0, Double_t mom=0,Int_t radiatorId=0, Int_t lensId=0, Double_t angle=0, Double_t z=0, Double_t x=0, Double_t xstep=0, Double_t ystep=0){
  
  if(infile=="") return;
  TString fileid(infile);
  fileid.Remove(0,fileid.Last('/')+1);
  fileid.ReplaceAll("C.root","");
  prt_data_info = getDataInfo(fileid);
  studyId=prt_data_info.getStudyId();

  TString savepath="data/procData";
  if(studyId>0) {
    savepath = infile;
    savepath.Remove(savepath.Last('/'));
    savepath += Form("/P%d",studyId);
  }
  Double_t mult(0),le1(0),le2(150),offset(0),timeres(0);

  if(infile.Contains("C.root")) { // beam data
    le1=0;
    le2=100;
  }
  if(infile.Contains("SP.root")) { // PDF data
    timeres=0.2;
  }
  
  if(!prt_init(infile,1,savepath)) return;   

  
  TFile *file = new TFile(infile.ReplaceAll(".root",".res.root"),"recreate");
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

  Int_t colors[4]={2,4,2,4};
  TH1F * hLe[prt_maxdircch][4];
  for(Int_t t=0; t<4; t++){
    for(Int_t i=0; i<prt_maxdircch; i++){
      hLe[i][t] = new TH1F(Form("le_ch_%d_%d",t,i),"LE; LE [ns]; entries [#]",500,le1,le2);
      hLe[i][t]->SetLineColor(colors[t]);
    }    
  }
  
  gStyle->SetOptStat(1001111);
  gStyle->SetOptFit(1111);
 
  PrtHit hit;
  for(int ievent=0; ievent<prt_entries; ievent++){
    prt_nextEvent(ievent,1000);
    Int_t counts(0);
    Double_t tot(0),time(0);
    for(int i=0; i<prt_event->GetHitSize(); i++){
      hit = prt_event->GetHit(i);
      Int_t mcpid = hit.GetMcpId();
      Int_t pixid = hit.GetPixelId()-1;
      Double_t time = hit.GetLeadTime();
      Int_t ch = map_mpc[mcpid][pixid];
      
      if(ch<prt_maxdircch){
	time = hit.GetLeadTime()-offset;
	//if(timeres>0) time=prt_rand.Gaus(time,timeres);
	tot = hit.GetTotTime();
	hLe[ch][prt_pid]->Fill(time);
	hLeA->Fill(time);
	hTot->Fill(tot);

	if(time<le2 && time>le1){
	  Int_t mcpid = hit.GetMcpId();
	  Int_t pixid = hit.GetPixelId()-1;
	  prt_hdigi[mcpid]->Fill(pixid%8, pixid/8);
	  counts++;
	}
      }
    }

    if(counts>5) hMult->Fill(counts);
  }

  angle = prt_theta;
  TString ext = Form("_%d_%d_%d",studyId,fileId,(Int_t)prt_theta);
  // TCanvas *cExport = new TCanvas("cExport","cExport",0,0,800,400);
  // cExport->SetCanvasSize(800,400);
  // for(Int_t i=0; i<maxch_dirc; i++){
  //   cExport->cd();
  //   prt_normalize( hLe[i][0], hLe[i][1]);
  //   hLe[i][0]->Draw();
  //   hLe[i][1]->Draw("same");
    
  //   cExport->SetName(Form("hLe_%d",i));
  //   canvasAdd(cExport);
  //   canvasSave(1,0);
  // }  
  
  prt_drawDigi("m,p,v\n",2018,0,0);
  prt_cdigi->cd();

  TLatex Tl;
  Tl.SetTextAlign(23);
  Tl.SetTextSize(0.05);
  Tl.DrawLatex(0.48,0.99,Form("%d deg",(Int_t)prt_theta));
  
  prt_cdigi->SetName(Form("hp_%d_%d",(Int_t)prt_theta,(Int_t)prt_test1)+ext);
  prt_canvasAdd(prt_cdigi);

  prt_canvasAdd("p_le"+ext,800,400);
  prt_fit(hLeA,0.3,100,100).X();
  hLeA->Draw();
  
  prt_canvasAdd("p_tot"+ext,800,400);
  hTot->Draw();

  prt_canvasAdd("p_mult"+ext,800,400);
  mult = prt_fit(hMult,20,20,100).X();
  hMult->Draw();

  tree->Fill();
  tree->Write();
  
  prt_canvasSave(1,0);
}
