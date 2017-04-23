#define prt__beam
#include "../prtdirc/src/PrtHit.h"
#include "../prtdirc/src/PrtEvent.h"
#include "prttools.C"

void procRaw(TString in="../prtdirc/build/hits.root"){

  if(!prt_init(in,1,"data/procRaw")) return;
  
  Int_t studyId(0),lensId(0),radiatorId(0);
  Double_t nph(0),z(0),x(0),xstep(0),ystep(0),t1(0),t2(0),
    le1(0),le2(50),offset(0),timeres(0);
  
  TFile *file = new TFile(in.ReplaceAll(".root","_proc.root"),"recreate");
  TTree *tree= new TTree("proc","proc");
  tree->Branch("studyId", &studyId,"studyId/I");
  tree->Branch("geometry", &prt_geometry,"prt_geometry/I");  
  tree->Branch("mom", &prt_mom,"prt_mom/D");
  tree->Branch("radiatorId", &radiatorId,"radiatorId/I");
  tree->Branch("lensId", &lensId,"lensId/I");
  tree->Branch("theta", &prt_theta,"prt_theta/I");
  tree->Branch("z", &z,"z/D");
  tree->Branch("x", &x,"x/D");
  tree->Branch("xstep", &xstep,"xstep/D");
  tree->Branch("ystep", &ystep,"ystep/D");
  tree->Branch("test1", &prt_test1,"prt_test1/D");
  tree->Branch("test2", &prt_test2,"prt_test2/D");
  tree->Branch("nph",&nph,"nph/D");
  
  
  TH1F * hNph = new TH1F("nph","Mult;multiplicity [#];entries [#]",300,0,300);
  TH1F * hLeA  = new TH1F("le","LE; LE [ns]; entries [#]",1000,le1,le2);
  TH1F * hTot  = new TH1F("tot","TOT; TOT [#]; entries [#]",1000,-2,15);

  TH1F * hLe[prt_maxch][5];
  for(Int_t t=0; t<5; t++){
    for(Int_t i=0; i<prt_maxch; i++){
      hLe[i][t] = new TH1F(Form("le_ch_%d_%d",t,i),"LE; LE [ns]; entries [#]",500,le1,le2);
      hLe[i][t]->SetLineColor(prt_color[t]);
    }    
  }
  
  gStyle->SetOptStat(1001111);
  gStyle->SetOptFit(1111);
 
  PrtHit hit;
  for (Int_t ievent=0; ievent<prt_entries; ievent++){
    prt_nextEvent(ievent,1000);
    Int_t counts(0);
    Double_t tot(0),time(0);
    
    for(Int_t i=0; i<prt_event->GetHitSize(); i++){
      hit = prt_event->GetHit(i);
      Int_t ch = hit.GetChannel();
      if(ch==-1) ch = map_mpc[hit.GetMcpId()][hit.GetPixelId()-1];      
      if(ch<prt_maxdircch && !prt_isBadChannel(ch)){
	time = hit.GetLeadTime()-offset;
	if(timeres>0) time=prt_rand.Gaus(time,timeres);
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
    if(counts>0) hNph->Fill(counts);
  }
std::cout<<"prt_theta "<<prt_theta <<std::endl;
 
  TString ext = Form("_%d",studyId);
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
  
  
  prt_canvasAdd("p_le"+ext,800,400);
  prt_fit(hLeA,0.3,100,100).X();
  hLeA->Draw();
  
  prt_canvasAdd("p_tot"+ext,800,400);
  hTot->Draw();

  prt_canvasAdd("p_nph"+ext,800,400);
  nph = prt_fit(hNph,20,20,100).X();
  hNph->Draw();
  
  prt_drawDigi("m,p,v\n",prt_geometry,-2,-2);
  prt_cdigi->SetName("p_hp"+ext);
  prt_canvasAdd(prt_cdigi);  
  
  tree->Fill();
  tree->Write();
  
  prt_canvasSave(1,0);
}
