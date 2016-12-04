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

void createPdf(TString path="/data.local/data/jun15/beam_15183022858C.root", Int_t normtype=1 ,Bool_t save=false){
  // path="/data.local/data/jun15/beam_15177135523S.root";
  // path="~/simo/build/beam_15184203911SP.root";
  if(path=="") return;
  
  fSavePath = "data/pdf3";
  PrtInit(path,1);
  gStyle->SetOptStat(0);
  CreateMap();

  Int_t totalmcp[5][9], totalmcpr[5][9];
  for(Int_t i=0; i<5; i++){
    for(Int_t m=0; m<9; m++){
      totalmcp[i][m]=0;
      totalmcpr[i][m]=0;
    }
  }
  
  TH1F *hlef[maxch_dirc], *hles[maxch_dirc];

  for(Int_t i=0; i<maxch_dirc; i++){
    hlef[i] = new TH1F(Form("lef_%d",i),"pdf;LE time [ns]; entries [#]",500,0,50);
    hles[i] = new TH1F(Form("les_%d",i),"pdf;LE time [ns]; entries [#]",500,0,50);
  }
  
  Double_t time;
  PrtHit fHit;
  Int_t totalf(0),totals(0), ch, entries = fCh->GetEntries();
  Int_t start =0; //(path.Contains("C.root"))? 50000 : 0; 
  for (Int_t ievent=start; ievent<entries; ievent++){
    PrtNextEvent(ievent,1000);

    Int_t nHits =prt_event->GetHitSize();
    Double_t tof = prt_event->GetTest1();
    // if(tof>72 && tof<72.6) continue;
    
    if(prt_event->GetType()==0){
      Bool_t t1(false),t2(false),t3(false);
      Bool_t tof1(false), tof2(false);
      Bool_t hodo1(false), hodo2(false);
      for(Int_t h=0; h<nHits; h++) {
	fHit = prt_event->GetHit(h);
	Int_t gch=fHit.GetChannel();
 
       	if(gch==818)
	  t1=true;
	//if(gch==821)
	  t2=true;
	if(gch==819)
	  t3=true;

	//if(gch>1031 && gch<1034)
	  tof1=true;
	//if(gch>651 && gch<657)
	  tof2=true;

	//if(gch>776 && gch<=780) //4
	//if(gch>777 && gch<779) //2
	  hodo1=true;
	//if(gch>=790 && gch<794) //4
	//if(gch>791 && gch<793) //2
	if(gch>=790 && gch<=794)
	  hodo2=true;
      }

      if(!(t1 && t2 && t3 && tof1 && tof2 && hodo1 && hodo2)) continue;
    }

    Int_t mult[maxch];
    memset(mult, 0, sizeof(mult));

    for(Int_t i=0; i<nHits; i++){
      fHit = prt_event->GetHit(i);
      Int_t mcp = fHit.GetMcpId();
      Int_t pix = fHit.GetPixelId()-1;
      ch=map_mpc[mcp][pix];      
      time=fHit.GetLeadTime();//+gRandom->Gaus(0,0.3);
      // Double_t tot= fHit.GetTotTime();
      // if(tot<2 || tot>4) continue;
      

      if(++mult[ch]>1) continue;      
      if(time<5 || time>50) continue;
      
      if(prt_event->GetParticle()==2212){
	totalmcp[4][mcp]++;
	totalmcpr[4][mcp%3]++;
	totalf++;
	hlef[ch]->Fill(time);
      }
      if(prt_event->GetParticle()==211){
	totalmcp[2][mcp]++;
	totalmcpr[2][mcp%3]++;
	totals++;
	hles[ch]->Fill(time);
      }
      fhDigi[mcp]->Fill(pix%8, pix/8);
    }
    // if(prt_event->GetParticle()==2212) totalf++;
    // if(prt_event->GetParticle()==211 || prt_event->GetParticle()==212) totals++;
  }

  std::cout<<"#1 "<< totalf <<"  #2 "<<totals <<std::endl;

  TCanvas *cExport = new TCanvas("cExport","cExport",0,0,800,400);
  
  if(totalf>0 && totals>0) {
    path.ReplaceAll("*","");
    path.ReplaceAll(".root",Form(".pdf%d.root",normtype));
    
    TFile efile(path,"RECREATE");
    
    for(Int_t i=0; i<maxch_dirc; i++){
      Int_t mcp = i/64;

      if(normtype==1){
	hlef[i]->Scale(1/(Double_t)totalf);
	hles[i]->Scale(1/(Double_t)totals);
      }
      
      if(normtype==2){
	hlef[i]->Scale(1/(Double_t)totalmcp[4][mcp]);
	hles[i]->Scale(1/(Double_t)totalmcp[2][mcp]);
      }
	    
      if(normtype==3){
	hlef[i]->Scale(1/(Double_t)totalmcpr[4][mcp%3]);
	hles[i]->Scale(1/(Double_t)totalmcpr[2][mcp%3]);
      }

      
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
      
      if(save){
	cExport->cd();
      	hlef[i]->GetXaxis()->SetRangeUser(0,50);
	hles[i]->GetXaxis()->SetRangeUser(0,50);
	prt_normalize(hlef[i],hles[i]);
	hlef[i]->SetLineColor(2);
      	hles[i]->SetLineColor(4);
	hlef[i]->Draw();
      	hles[i]->Draw("same");
	cExport->SetName(Form("pdf_mcp%d_pix_%d",map_mcp[i],map_pix[i]));
	canvasAdd(cExport);
      	canvasSave(1,0);
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
