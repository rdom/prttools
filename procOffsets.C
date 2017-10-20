#define prt__beam
#include "../prtdirc/src/PrtHit.h"
#include "../prtdirc/src/PrtEvent.h"
#include "prttools.C"

void procOffsets(TString path="",Int_t corrected=1){
  
  if(path=="") return;
  TString fileid(path);
  fileid.Remove(0,fileid.Last('/')+1);
  fileid.Remove(fileid.Last('.')-1);
  
  prt_data_info = getDataInfo(fileid);
 
  Int_t h1a(0),h1b(400),h2a(0),h2b(400),hbin(5000);//h1a(200),h1b(400)

  if(corrected==1){
    h1a=0;
    h1b=50;
    h2a=0;
    h2b=50;
    hbin=2000;
  }
  
  TString outdir=path;outdir.Remove(outdir.Last('/'));
  TString sstudy=outdir; sstudy.Remove(0,sstudy.Last('/'));  
  TString insim = path; insim.ReplaceAll("C.root","S.root");
  
  if(!prt_init(path,1,outdir+Form("/%da/%d",prt_data_info.getStudyId(),prt_data_info.getFileId()))) return;
  prt_ch->Add(insim);

  TH1F * hLeD  = new TH1F("leD","LE beam data ; LE [ns]; entries [#]",hbin,h1a,h1b);
  TH1F * hLeS  = new TH1F("leS","LE simulation; LE [ns]; entries [#]",hbin,h2a,h2b);
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  PrtHit hit;
  Int_t maxent(0);
  for (Int_t ievent=0; ievent< prt_ch->GetEntries(); ievent++){
    prt_nextEvent(ievent,10000);
    if(prt_event->GetParticle()!=2212) continue;
    bool bsim(false);
    TString current_file_name  = prt_ch->GetCurrentFile()->GetName();
    if(current_file_name.Contains("S.root")) bsim = true;
    else maxent++;
    //if(!bsim && maxent>10000) continue;
    
    Double_t time(0);
    for(Int_t i=0; i<prt_event->GetHitSize(); i++){
      hit = prt_event->GetHit(i);

      if(hit.GetChannel()<960 ){
	if(!bsim) hLeD->Fill(hit.GetLeadTime());
	else hLeS->Fill(hit.GetLeadTime());
      }
    }
  }
  
  prt_canvasAdd("offset",800,400);
  hLeS->SetLineColor(kRed);
  hLeS->Draw();
  hLeD->Draw("same");
  TLegend *leg = new TLegend(0.62,0.7,0.92,0.9);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(hLeD,"beam data ","lp");
  leg->AddEntry(hLeS,"simulation","lp");
  leg->Draw();
    
  prt_canvasSave(1,0);

  if(corrected==0){
    double xmax1 = hLeD->GetXaxis()->GetBinCenter(hLeD->GetMaximumBin());
    double xmax2 = hLeS->GetXaxis()->GetBinCenter(hLeS->GetMaximumBin());
    TFile efile(path+ ".off.root","RECREATE");
    TGraph *gr = new TGraph();
    gr->SetPoint(0,xmax1-xmax2,  xmax1-xmax2);
    fileid=path;
    fileid.Remove(0,fileid.Last('_')+1);
    fileid.Remove(fileid.Last('C'));
    gr->SetName("off_"+fileid);    
    gr->Write();

    gr = new TGraph();
    gr->SetPoint(0,prt_theta, xmax1);
    gr->SetName("offbeam");    
    gr->Write();
    
    
    gr = new TGraph();
    gr->SetPoint(0,prt_theta, xmax2);
    gr->SetName("offsim");    
    gr->Write();    
    
    efile.Write();
    efile.Close();
    std::cout<<"new offset "<< xmax1-xmax2 <<std::endl;    
    
  }
}
