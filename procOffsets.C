#define prt__beam
#include "../prtdirc/src/PrtHit.h"
#include "../prtdirc/src/PrtEvent.h"
#include "prttools.C"

void procOffsets(TString path="/data.local/data/jun15/beam_15180221900C.root",Int_t rawdata=1){
  
  if(path=="") return;

  Int_t h1a(200),h1b(400),h2a(0),h2b(100);

  if(rawdata==1){
    h1a=0;
    h1b=50;
    h2a=0;
    h2b=50;
  }
  
  TString outdir=path;outdir.Remove(outdir.Last('/'));
  TString sstudy=outdir; sstudy.Remove(0,sstudy.Last('/'));
  fSavePath = outdir+sstudy;

  
  TString insim = path; insim.ReplaceAll("C.root","S.root");
  PrtInit(path,1);
  fCh->Add(insim);
  fNEntries = fCh->GetEntries();

  TH1F * hLeD  = new TH1F("leD","LE beam data ; LE [ns]; entries [#]",2000,h1a,h1b);
  TH1F * hLeS  = new TH1F("leS","LE simulation; LE [ns]; entries [#]",2000,h2a,h2b);
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  PrtHit fHit;
  Int_t maxent(0);
  for (Int_t ievent=0; ievent<fNEntries; ievent++){
    PrtNextEvent(ievent,1000);
    if(fEvent->GetParticle()!=2212) continue;
    bool bsim(false);
    TString current_file_name  = fCh->GetCurrentFile()->GetName();
    if(current_file_name.Contains("S.root")) bsim = true;
    else maxent++;
    if(!bsim && maxent>10000) continue;
    
    Double_t time(0);
    for(Int_t i=0; i<fEvent->GetHitSize(); i++){
      fHit = fEvent->GetHit(i);

      if(fHit.GetChannel()<960 ){
	if(!bsim) hLeD->Fill(fHit.GetLeadTime());
	else hLeS->Fill(fHit.GetLeadTime());
      }
    }
  }
  canvasAdd("offset",800,400);
  hLeS->Draw();
  hLeD->SetLineColor(kRed);
  hLeD->Draw("same");
  TLegend *leg = new TLegend(0.62,0.7,0.92,0.9);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(hLeD,"beam data ","lp");
  leg->AddEntry(hLeS,"simulation","lp");
  leg->Draw();
    
  canvasSave(1,0);

  if(rawdata==0){
    TString fileid = path;
    fileid.Remove(0,fileid.Last('_')+1);
    fileid.Remove(fileid.Last('C'));
    double xmax1 = hLeD->GetXaxis()->GetBinCenter(hLeD->GetMaximumBin());
    double xmax2 = hLeS->GetXaxis()->GetBinCenter(hLeS->GetMaximumBin());
    TFile efile(path+ ".off.root","RECREATE");
    TGraph *gr = new TGraph();
    gr->SetPoint(0,xmax1-xmax2,  xmax1-xmax2);
    gr->SetName("off_"+fileid);
    gr->Write();
    efile.Write();
    efile.Close();
  }
}
