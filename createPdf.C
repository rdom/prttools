#define prt__beam
#include "../prtdirc/src/PrtHit.h"
#include "../prtdirc/src/PrtEvent.h"
#include "prttools.C"

TSpectrum *spect = new TSpectrum(2);
TF1 * fitpdf(TH1F *h){
  TF1 *gaust = NULL;
  Double_t rmin(0),rmax(60);
  Double_t integral = h->Integral(h->GetXaxis()->FindBin(rmin),h->GetXaxis()->FindBin(rmax));
  if(h->GetEntries()>50){
    Int_t nfound = spect->Search(h,2,"",0.1);
    std::cout<<"nfound  "<<nfound <<std::endl;
    if(nfound==1){
      gaust =new TF1("gaust","gaus(0)",rmin,rmax);
      gaust->SetNpx(500);
      gaust->FixParameter(1,spect->GetPositionX()[0]);
      gaust->SetParameter(2,0.3);
      gaust->SetParLimits(2,0.2,1);
    }
    if(nfound==2){
      gaust =new TF1("gaust","gaus(0)+gaus(3)",rmin,rmax);
      gaust->SetNpx(500);
      gaust->FixParameter(1,spect->GetPositionX()[0]);
      gaust->FixParameter(4,spect->GetPositionX()[1]);
      gaust->SetParameter(2,0.3);
      gaust->SetParameter(5,0.3);
      gaust->SetParLimits(2,0.2,1);
      gaust->SetParLimits(5,0.2,1);
    
      std::cout<<spect->GetPositionX()[0]<< " "<<spect->GetPositionX()[1] <<std::endl;
    
    }
    h->Fit("gaust","","MQN",0,60);
  }else{
    gaust =new TF1("gaust","pol0",rmin,rmax);
    gaust->FixParameter(0,0);
  }
	 
  return gaust;
}

void createPdf(TString path="/data.local/data/jun15/beam_15183022858C.root"){//beam_15177135523S.root
  //  path="/data.local/data/jun15/beam_15177135523S.root";
  //path="~/simo/build/beam_15184203911SP.root";
  if(path=="") return;
  
  fSavePath = "data/pdf3";
  PrtInit(path,1);
  gStyle->SetOptStat(0);
  CreateMap();

  TH1F *hlef[960], *hles[960];

  for(Int_t i=0; i<960; i++){
    hlef[i] = new TH1F(Form("lef_%d",i),"pdf;LE time [ns]; entries [#]", 1000,0,100);
    hles[i] = new TH1F(Form("les_%d",i),"pdf;LE time [ns]; entries [#]", 1000,0,100);
  }
  
  Double_t time;
  PrtHit fHit;
  Int_t totalf(0),totals(0), ch, entries = fCh->GetEntries();
  Int_t start = (path.Contains("C.root"))? 50000 : 0; 
  for (Int_t ievent=start; ievent<entries; ievent++){
    PrtNextEvent(ievent,1000);

    Int_t nHits =prt_event->GetHitSize();
    
    for(Int_t i=0; i<nHits; i++){
      fHit = prt_event->GetHit(i);
      Int_t mcpid = fHit.GetMcpId();
      Int_t pixid = fHit.GetPixelId()-1;
      ch=map_mpc[mcpid][pixid];      
      time=fHit.GetLeadTime();//+gRandom->Gaus(0,0.3);
      
      if(prt_event->GetParticle()==2212){
	//totalf++;
	hlef[ch]->Fill(time);
      }
      if(prt_event->GetParticle()==211 || prt_event->GetParticle()==212){
	//totals++;
	hles[ch]->Fill(time);
      }
      fhDigi[mcpid]->Fill(pixid%8, pixid/8);
    }
    if(prt_event->GetParticle()==2212) totalf++;
    if(prt_event->GetParticle()==211 || prt_event->GetParticle()==212) totals++;
  }

  std::cout<<"#1 "<< totalf <<"  #2 "<<totals <<std::endl;

  TCanvas *cExport = new TCanvas("cExport","cExport",0,0,800,400);
  
  if(totalf>0 && totals>0) {
    path.ReplaceAll("*","");
    path.ReplaceAll(".root",".pdf.root");
    TFile efile(path,"RECREATE");
    
    for(Int_t i=0; i<960; i++){
      hlef[i]->Scale(1/(Double_t)totalf);
      hles[i]->Scale(1/(Double_t)totals);

      // hlef[i]->Scale(1/(Double_t)(hlef[i]->GetEntries()));
      // hles[i]->Scale(1/(Double_t)(hlef[i]->GetEntries()));

      // Double_t mm = hles[i]->Integral();
      // if(hlef[i]->Integral()>hles[i]->Integral()) mm=hlef[i]->Integral();
      
      // hlef[i]-> Scale(1./mm, "width"); 
      // hles[i]-> Scale(1./mm, "width");
      
      // TF1 *f = fitpdf(hlef[i]);
      // TF1 *s = fitpdf(hles[i]);
      // for(Int_t p=0; p<f->GetNpar(); p++) f->FixParameter(p,f->GetParameter(p));
      // for(Int_t p=0; p<s->GetNpar(); p++) s->FixParameter(p,s->GetParameter(p));
      // f->SetName(Form("ff_%d",i));
      // s->SetName(Form("fs_%d",i));
      // f->Write();
      // s->Write();
      
      hlef[i]->SetName(Form("hf_%d",i));
      hles[i]->SetName(Form("hs_%d",i));
      hlef[i]->Write();
      hles[i]->Write();
      
      if(false){
      	cExport->cd();
      	//	canvasAdd(Form("pdf_%d",i),800,500);
      	cExport->SetName(Form("pdf_%d",i));
      	//canvasAdd(cExport);
      	//hlef[i]->GetYaxis()->SetRangeUser(0,1.5);
	prt_normalize(hlef[i],hles[i]);
	axisTime800x500(hlef[i]);
	axisTime800x500(hles[i]);
	hlef[i]->SetLineColor(2);
      	hlef[i]->Draw();
      	hles[i]->SetLineColor(4);
      	hles[i]->Draw("same");
      	// // f->Draw();
      	// s->SetLineColor(4);
      	// s->Draw("same");
      	cExport->Print(fSavePath+Form("/pdf_mcp%d_pix_%d.png",map_mcp[i],map_pix[i]));
      	//canvasSave(1,0);
      	//canvasDel(cExport->GetName());
      }
    }
    
    efile.Write();
    efile.Close();
  }

 
  canvasAdd("le",800,500);
  hlef[308]->SetLineColor(2);
  hlef[308]->Draw();
  hles[308]->SetLineColor(4);
  hles[308]->Draw("same");
  
  //  writeString(fSavePath+"/digi.csv", drawDigi("m,p,v\n",2,-2,-2));
  writeString(fSavePath+"/digi.csv", drawDigi("m,p,v\n",7,0,0));
  cDigi->SetName("hits");
  canvasAdd(cDigi);
  
  canvasSave(1,0);
}
