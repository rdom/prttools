#define prt__beam
#include "../prtdirc/src/PrtHit.h"
#include "../prtdirc/src/PrtEvent.h"
#include "prttools.C"
#include <TVirtualFitter.h>
#include <TKey.h>
#include <TRandom.h>

TLine *gLine = new TLine(0,0,0,1000);

void recoPdf(TString path="", TString pdfEnding=".pdf1.root", Double_t sigma=200,Bool_t debug=false, Double_t r1=0, Double_t r2=0, Int_t nforpdf=0){
  
  Int_t studyId;
  TString tpath=path;
  tpath.ReplaceAll("aug17","");
  sscanf(tpath, "%*[^0-9]%d{3}",&studyId);

  if(!prt_init(path,1,Form("data/recopdf_%d",studyId))) return;
  TGaxis::SetMaxDigits(4);
  
  TCanvas *cc = (debug)? new TCanvas("cc","cc",800,400): NULL ;

  TH1F *hpdff[prt_maxch],*hpdfs[prt_maxch], *hl[5],*hnph[5],*hll[5];
  TGraph *gpdff[prt_maxch],*gpdfs[prt_maxch];
  for(Int_t i=0; i<5; i++){
    hl[i] = new TH1F(Form("hl_%d",i),"pdf;LE time [ns]; entries [#]", 1000,0,50);
    hnph[i] = new TH1F(Form("hnph_%d",i),";detected photons [#]; entries [#]", 160,0,160);
    hll[i] = new TH1F(Form("hll_%d",i),";ln L(p) - ln L(#pi); entries",160,-40,40); //120,-60,60
  }  
  TH1F *hl3 = new TH1F("hl3","pdf;LE time [ns]; entries [#]", 1000,0,50);
  TH1F *hnphf =  new TH1F("hnphf","hnphf",200,0,200);
  TH1F *hnphs =  new TH1F("hnphs","hnphs",200,0,200);
  TH1F *hTof =  new TH1F("fTof",";TOF2-TOF1 [ns];entries [#]",400,30,34);
 
  Bool_t ismultnorm(false);
  //  if(pdfEnding.Contains("pdf0")) ismultnorm=true;
  TRandom rand;
  TF1 *pdff[prt_maxch],*pdfs[prt_maxch];
  TString pdf = path;
  pdf.ReplaceAll("*","");
  pdf.ReplaceAll(".root",pdfEnding);
  TFile f(pdf);

  Int_t binfactor = (Int_t)(sigma/50.+0.1);
  if(sigma >0) hl3->Rebin(binfactor);
  Double_t integ1(0), integ2(0);

  if(ismultnorm){
    hnphf = (TH1F*)f.Get("hnphf");
    hnphs = (TH1F*)f.Get("hnphs");
  }
  for(Int_t i=0; i<prt_maxdircch; i++){
    hpdff[i] = (TH1F*)f.Get(Form("hf_%d",i));
    hpdfs[i] = (TH1F*)f.Get(Form("hs_%d",i));
    hpdff[i]->SetLineColor(2);
    hpdfs[i]->SetLineColor(4);
    if(sigma >0) hpdff[i]->Rebin(binfactor);
    if(sigma >0) hpdfs[i]->Rebin(binfactor);
    //f hpdff[i]->Smooth(1);
    //f hpdfs[i]->Smooth(1);
    // hpdff[i]->Scale(1/hpdff[i]->Integral());
    // hpdfs[i]->Scale(1/hpdfs[i]->Integral());

    gpdff[i] = new TGraph(hpdff[i]);
    gpdfs[i] = new TGraph(hpdfs[i]); 
    
    integ1+= hpdff[i]->Integral();
    integ2+= hpdfs[i]->Integral();
    hl3->Add(hpdff[i]);
    hl3->Add(hpdfs[i]);
  }
  // for(Int_t i=0; i<prt_maxch; i++){
  //   hpdff[i]->Scale(1/integ1);
  //   hpdfs[i]->Scale(1/integ2);
  // }  
 
  if(path.Contains("C.root")) sigma=0;
  if(path.Contains("Z.root")) sigma=0;
  if(path.Contains("F.root")) sigma=400;

  TF1 *F1 = new TF1("gaus0","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2])",0,150);
  TF1 *F2 = new TF1("gaus1","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2])",0,150);
  F1->SetParameter(0,1);
  F2->SetParameter(0,1);      
  F1->SetParameter(1,63);
  F2->SetParameter(1,50);
  F1->SetParameter(2,11);
  F2->SetParameter(2,9);

  Int_t countall[prt_nmcp][64],countgood[prt_nmcp][64],countbad[prt_nmcp][64];
  Int_t mcpf[prt_nmcp], mcps[prt_nmcp];
      
  for (Int_t m=0; m<prt_nmcp; m++) {
    mcpf[m]=0;
    mcps[m]=0;
    for(Int_t p=0; p<prt_npix; p++){
      countall[m][p]=0;
      countgood[m][p]=0;
      countbad[m][p]=0;
    }
  }  
  
  Double_t theta(0);
  TVirtualFitter *fitter;
  Double_t nph,time,timeres(-1);
  PrtHit fHit;
  Int_t totalf(0),totals(0),mcp,pix,ch, entries = prt_entries; // [50000-rest] - is for pdf generation
  if(path.Contains("F.root")) entries = prt_entries;
  if(path.Contains("S.root")) entries = 4000;
  if(path.Contains("C.root")) entries = 100000;

  Int_t trigT1(816);
  Int_t trigT2(817);
  Int_t trigT3h(818);
  Int_t trigT3v(819);
  Int_t trigTof1(1392);
  Int_t trigTof2(1398);
  
  for (Int_t ievent=0; ievent<entries; ievent++){
    prt_nextEvent(ievent,1000);
    timeres = prt_event->GetTimeRes();
    Double_t aminf,amins, sum(0),sumf(0),sums(0);
    Int_t nGoodHits(0), nHits =prt_event->GetHitSize();
    if(prt_event->GetType()==0){


      // // 332
      // if(fabs(prt_event->GetMomentum().Mag()-7)<0.1){
      // 	if( prt_event->GetParticle()==2212 && prt_event->GetTest1()<32.65 ) continue;
      // 	if( prt_event->GetParticle()==211  && prt_event->GetTest1()>31.68 ) continue;
      // }

      // //plate
      // if(fabs(prt_event->GetMomentum().Mag()-7)<0.1){
      // 	if( prt_event->GetParticle()==2212 && prt_event->GetTest1()<32.5 ) continue;
      // 	if( prt_event->GetParticle()==211  && prt_event->GetTest1()>31.8 ) continue;
      // }
	
      hTof->Fill(prt_event->GetTest1());

      Bool_t t1(1),t2(0),t3h(0),t3v(0);
      Bool_t tof1(1), tof2(1);
      Bool_t hodo1(0), hodo2(0);
      
      for(Int_t h=0; h<nHits; h++) {
      	fHit = prt_event->GetHit(h);
      	Int_t gch=fHit.GetChannel();	
	
	if(gch==trigT2)
	  t2=true;
	if(gch==trigT3h)
	  t3h=true;
	if(gch==trigT3v)
	  t3v=true;

	//if(gch>=1350 && gch<=1351)
	//if(gch>=1351 && gch<=1352)
	if(gch>=1350 && gch<=1350)
	  hodo1=true;
	//if(gch>=1369 && gch<=1370)
	//if(gch>=1364 && gch<=1374)
	hodo2=true;	  
      }
      
      if(!( t1 && t2 && t3h && t3v && tof1 && tof2 && hodo1 && hodo2)) continue;
    }

    if(debug) std::cout<<"===================== event === "<< ievent <<std::endl;
     // if(prt_pid==2 && hll[2]->GetEntries()>2200)continue;
     // if(prt_pid==4 && hll[4]->GetEntries()>2200) continue;  
    
    

    Int_t mult[prt_maxch];
    memset(mult, 0, sizeof(mult));
    for(Int_t i=0; i<nHits; i++){
      fHit = prt_event->GetHit(i);
      mcp = fHit.GetMcpId();
      pix=fHit.GetPixelId()-1;
      ch = map_mpc[mcp][pix];
      time = fHit.GetLeadTime();

      // //cut-off from c2017
      // if(mcp%3==0 && pix<32) continue;
      // if(mcp%3==2 && pix>=32) continue; 

      if(ch>prt_maxdircch) continue;
      if(prt_event->GetType()!=0) time += rand.Gaus(0,sigma*0.001);
      //      if(++mult[ch]>1 || ch ==0) continue;
 
      prt_hdigi[mcp]->Fill(pix%8, pix/8);

      { //time cuts
      	// Double_t cut1(11);
      	// if(studyId==157 || studyId==155) cut1=8;
      	// if(theta<=80){
      	//   if(time<cut1 || time>75) continue; //45
      	// }else if(theta>94){
      	//   if(time<3 || time>40) continue; //40
      	// }
	//if(time<0 || time>50) continue;
      }
      nGoodHits++;
      // aminf = hpdff[ch]->GetBinContent(hpdff[ch]->FindBin(time-0.0)); 
      // amins = hpdfs[ch]->GetBinContent(hpdfs[ch]->FindBin(time-0.0));   

      aminf = gpdff[ch]->Eval(time);
      amins = gpdfs[ch]->Eval(time);
      
      countall[mcp][pix]++;
      if(prt_pid==4){
	mcpf[mcp]++;
	//	if(mcp ==7 || mcp ==10) continue;
	countgood [mcp][0]++;
	if(aminf>amins) countgood [mcp][pix]++;
	else countbad[mcp][pix]++;
      }else if (prt_pid==2){
	mcps[mcp]++;
	//	if(mcp ==1 || mcp ==3 || mcp ==5) continue;
	countbad[mcp][0]++;
	if(amins>aminf) countgood [mcp][pix]++;
        else countbad[mcp][pix]++;
      }
      // if(fabs(aminf-amins)/(aminf+amins)*0.5<0.01) continue;
      
      if(debug){
	TString x=(aminf>amins)? " <====== PROTON" : "";
	std::cout<<Form("f %1.6f s %1.6f mcp %d pix %d   pid %d",aminf,amins,mcp,pix  ,prt_particle)<<"  "<<x <<std::endl;
	
	cc->cd();	
	prt_axisTime800x500(hpdff[ch]);
	prt_axisTime800x500(hpdfs[ch]);
	prt_normalize(hpdff[ch],hpdfs[ch]);
	hpdff[ch]->SetLineColor(2);
	hpdfs[ch]->SetLineColor(4);
	hpdff[ch]->Draw();
	hpdff[ch]->SetTitle(Form("mcp=%d  pix=%d",mcp,pix));
	hpdff[ch]->GetXaxis()->SetTitle("LE time [ns]");
	hpdff[ch]->GetYaxis()->SetTitle("PDF value");
	hpdff[ch]->GetXaxis()->SetRangeUser(0,40);
	hpdfs[ch]->Draw("same");
	gpdff[ch]->Draw("PL same");
	gpdfs[ch]->Draw("PL same");
	cc->Update();
	gLine->SetLineWidth(2);
	gLine->SetX1(time);
	gLine->SetX2(time);
	gLine->SetY1(cc->GetUymin());
	gLine->SetY2(cc->GetUymax());
	gLine->Draw();
	cc->Update();
	cc->WaitPrimitive();
      }
      // if(aminf==0 || amins==0) continue;

      Double_t noise = 1e-6; //1e-7; // nHits //1e-5
      
      sumf+=TMath::Log((aminf+noise));
      sums+=TMath::Log((amins+noise));
      
      // Double_t res;
      // if(aminf>amins){
      // 	res=100*(aminf-amins)/amins;
      // 	sumf+=TMath::Log((res+noise));
      // }else{
      // 	res=100*(amins-aminf)/amins;
      // 	sums+=TMath::Log((res+noise));
      // }      
      hl[prt_pid]->Fill(time);
    }
    
    if(nGoodHits<5) continue;
    hnph[prt_pid]->Fill(nGoodHits);
    
    if(ismultnorm){
      Double_t lhmf = hnphf->GetBinContent(hnphf->FindBin(nGoodHits)); 
      Double_t lhms = hnphs->GetBinContent(hnphs->FindBin(nGoodHits)); 
      sum = sumf+TMath::Log(lhmf)-(sums+TMath::Log(lhms));
    }else sum = sumf-sums;    

    // if(fabs(sum)<0.04) continue;    
    // sumf += 10*TMath::Log(F2->Eval(nHits));
    // sums += 10*TMath::Log(F1->Eval(nHits));
    // sum = sumf-sums;
    
    hll[prt_pid]->Fill(sum);
  }

  for (Int_t m=0; m <prt_nmcp; m++) {
    std::cout<<mcpf[m]<< " "<< mcps[m]<<std::endl;    
    for(Int_t p=0; p<prt_npix; p++){
      //prt_hdigi[m]->Fill(p%8,p/8,countgood[m][0]/(Double_t)countbad[m][0]);
    }
  }


  gStyle->SetOptStat(0);
  //gStyle->SetOptTitle(0);
  
  //  prt_drawDigi("m,p,v\n",2017,1.3,0);
  prt_drawDigi("m,p,v\n",2017);
  prt_canvasAdd(prt_cdigi);

  TString name = Form("_%d_%d_%1.1f_m%1.1f_x%d_z%d.root",studyId,prt_theta,sigma,prt_mom,prt_beamx,prt_beamz);
  if(path.Contains("C.root")) name =  "tid"+ name;
  else name = "tis"+ name;
  
  prt_canvasAdd("ll_"+name,800,400);
  
  prt_normalize(hll[4],hll[2]);
  hll[4]->GetYaxis()->SetNdivisions(9,5,0);
  
  TF1 *ff;
  Double_t sep(0),esep(0),m1,m2,s1,s2,dm1,dm2,ds1,ds2; 
  if(hll[4]->GetEntries()>10 && hll[2]->GetEntries()>10){
    hll[4]->Fit("gaus","Sq","");
    ff = hll[4]->GetFunction("gaus");
    ff->SetLineColor(1);
    m1=ff->GetParameter(1);
    s1=ff->GetParameter(2);
    dm1=ff->GetParError(1);
    ds1=ff->GetParError(2);

    hll[2]->Fit("gaus","Sq");
    ff = hll[2]->GetFunction("gaus");
    ff->SetLineColor(1);
    m2=ff->GetParameter(1);
    s2=ff->GetParameter(2);
    dm2=ff->GetParError(1);
    ds2=ff->GetParError(2);

    sep = (fabs(m1-m2))/(0.5*(s1+s2));
    
    Double_t e1,e2,e3,e4;
    e1=2/(s1+s2)*dm1;
    e2=2/(s1+s2)*dm2;
    e3=-((2*(m1 + m2))/((s1 + s2)*(s1 + s2)))*ds1;
    e4=-((2*(m1 + m2))/((s1 + s2)*(s1 + s2)))*ds2;
    esep=sqrt(e1*e1+e2*e2+e3*e3+e4*e4);
  }
  
  hll[4]->SetTitle(Form("separation = %1.2f",sep));
  hll[4]->SetLineColor(2);
  hll[4]->Draw();
  hll[2]->SetLineColor(4);
  hll[2]->Draw("same");

  TLegend *leg = new TLegend(0.70,0.65,0.88,0.87);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(hll[2],"pions ","lp");
  leg->AddEntry(hll[4],"protons","lp");
  leg->Draw();
  
  // prt_canvasAdd("hl_"+name,800,500);
  // // hl[4]->Scale(1/hl[4]->GetMaximum());
  // // hl[2]->Scale(1/hl[2]->GetMaximum());
  // // hl3->Scale(1/hl3->GetMaximum());

  // prt_normalize(hl[4],hl[2]);
  
  // hl[4]->Draw();
  // hl[2]->SetLineColor(4);
  // hl[2]->Draw("same");
  // hl3->SetLineColor(2);
  // hl3->Draw("same");

  //  gStyle->SetOptFit(1);
  prt_canvasAdd("hnph_"+name,800,400);

  prt_normalize(hnph[4],hnph[2]);
  hnph[4]->Fit("gaus");
  ff = hnph[4]->GetFunction("gaus");
  ff->SetLineColor(1);
  nph=ff->GetParameter(1);  
  hnph[4]->SetLineColor(2);
  hnph[4]->Draw();
  hnph[2]->SetLineColor(4);
  hnph[2]->Draw("same");

  TLegend *leg1 = new TLegend(0.65,0.65,0.83,0.87);
  leg1->SetFillColor(0);
  leg1->SetFillStyle(0);
  leg1->SetBorderSize(0);
  leg1->SetFillStyle(0);
  leg1->AddEntry(hnph[2],"pions ","lp");
  leg1->AddEntry(hnph[4],"protons","lp");
  leg1->Draw();
  
  
  prt_canvasAdd("hTof_"+name,800,400);
  hTof->Draw();
  // gLine->SetLineWidth(2);
  // gLine->SetX1(32.65);
  // gLine->SetX2(32.65);
  // gLine->SetY1(gPad->GetUymin());
  // gLine->SetY2(gPad->GetUymax());
  //gLine->Draw();
  
  //prt_canvasSave(0);
  
  std::cout<<dm1<<" "<<dm2<<" "<<ds1 <<" "<<ds2<<std::endl; 
  std::cout<<path<<" separation "<< sep <<" +/- "<<esep <<std::endl;
  std::cout<<"entries:  pi "<<hll[2]->GetEntries()<<" p "<<hll[4]->GetEntries() <<std::endl;

  path.ReplaceAll("S.root","R.root");
  if(path.Contains("X.root")) path.ReplaceAll("S.root","R.root");
  else path=prt_savepath+"/reco_"+name;

  if(nforpdf!=0) path=prt_savepath+Form("/reco_%d.root",nforpdf);
  
  TFile fc(path,"recreate");
  TTree *tc = new TTree("reco","reco");
  tc->Branch("theta",&prt_theta,"prt_theta/I");
  tc->Branch("phi",&prt_phi,"prt_phi/D");
  tc->Branch("sep",&sep,"sep/D");
  tc->Branch("esep",&esep,"esep/D");
  tc->Branch("sigma",&sigma,"sigma/D");
  tc->Branch("mom",&prt_mom,"prt_mom/D");
  tc->Branch("nph",&nph,"nph/D");
  tc->Branch("r1",&r1,"r1/D");
  tc->Branch("r2",&r2,"r2/D");
  tc->Branch("beamz",&prt_beamz,"prt_beamz/I");
  tc->Branch("beamx",&prt_beamx,"prt_beamx/I");
  tc->Branch("nforpdf",&nforpdf,"nforpdf/I");
  tc->Fill();
  tc->Write();
  fc.Write();
  fc.Close();
}
