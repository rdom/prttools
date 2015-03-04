// mdisplay - tool to plot different quantities from the *M.root file
// original author: Roman Dzhygadlo - GSI Darmstadt 

#define prt__beam
#include "prttools.C"

#include "TStyle.h"
#include <TApplication.h>
#include <TGClient.h>
#include <TGButton.h>
#include <TGFrame.h>
#include <TFrame.h>
#include <TRootEmbeddedCanvas.h>
#include <TGStatusBar.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TF1.h>
#include "TSpectrum.h"
#include <TRandom.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TGNumberEntry.h>
#include <TGLabel.h>
#include <TGListBox.h>
#include <TGComboBox.h>
#include <TGButtonGroup.h>
#include <TGTextEntry.h>
#include <TGProgressBar.h>
#include <TThread.h>
#include <TProof.h>
#include <TGSplitter.h>
#include <TGSlider.h>

#include <sstream>

#include "../mz-unpacker-BarrelDirc/TPrtHit.h"
#include "../mz-unpacker-BarrelDirc/TPrtEvent.h"
#include "mdisplay.h"

class TPrtHit;
class TPrtEvent;

MyMainFrame *gMain;
const int nmcp = 15, npix = 64;
TH1F *hPTime[nmcp][npix];
TH1F *hPTot[nmcp][npix];
TH1F *hPMult[nmcp][npix];
TH1F *hEMult[nmcp+1];
TH2F *hLeTot[nmcp][npix];
TH2F *hShape[nmcp][npix];

MSelector            *fSelector;
const Int_t maxch =2000;
TGraph *gGrDiff[maxch];
Int_t chmap[nmcp][npix];

const Int_t maxMult = 30;
TH1F *hTotM[maxMult];
TH1F *hLeM[maxMult];

TH1F *hTot,*hLe,*hMult,*hCh;
TH1F *hMultEvtNum1,*hMultEvtNum2;

TCanvas *cTime;
Int_t gComboId=0, gTrigger = 0, gMode=0;
TGHProgressBar *pbar;
TString ginFile; 

Double_t gTimeCutMin=-10000, gTimeCutMax=10000;
Double_t gMultCutMin=0, gMultCutMax=0;  
TSpectrum *spect = new TSpectrum(2);;
Double_t gTimeCuts[nmcp][npix][2];
Double_t gTotMean[nmcp][npix];
TString gsTimeCuts = "0";
TString gsTotMean = "0";


void MSelector::Init(TTree *tree){
  fChain = tree; 
  fChain->SetBranchAddress("TPrtEvent", &fEvent);
}

void PrintStressProgress(Long64_t total, Long64_t processed, Float_t, Long64_t){
  pbar->SetPosition(100*processed/(Float_t)total);
}

void init(){
  SetRootPalette(1);
  fCh = new TChain("M");
  fCh->Add(ginFile);
  fNEntries = fCh->GetEntries();
  std::cout<<"Entries in chain:  "<< fNEntries<<std::endl;

  TString workers = "workers=4";
  if(gSystem->GetFromPipe("whoami")=="hadaq" && fNEntries>1000000) workers = "workers=12";

  TProof *proof = TProof::Open(workers);
  TString dir = gSystem->pwd();

  gProof->Exec("gSystem->Load(\""+ dir+"/../mz-unpacker-BarrelDirc/TPrtHit_cpp.so\")");
  gProof->Exec("gSystem->Load(\""+ dir+"/../mz-unpacker-BarrelDirc/TPrtEvent_cpp.so\")");
  proof->Load("mdisplay.C+");

  proof->SetPrintProgress(&PrintStressProgress);
  fCh->SetProof();
  fSelector = new MSelector();
  gStyle->SetOptStat(1001111);
  gStyle->SetOptFit(1111);

  // create channel - mcp/pixel map
  for(Int_t ch=0; ch<maxch; ch++){
    Int_t mcp = ch/128;
    Int_t pix = (ch - mcp*128)/2;
    Int_t col = pix/2 - 8*(pix/16);
    Int_t row = pix%2 + 2*(pix/16);
    pix = row*8+col;
    chmap[mcp][pix]=ch;
  }
}

void MSelector::SlaveBegin(TTree *){
  TString option = GetOption();
  std::istringstream source(option.Data());
  Int_t bins1 =100, bins2 = 100, min1=-5, max1=5, min2=-5, max2=5;
  source>>gMode>>gTrigger>>bins1>>min1>>max1>>bins2>>min2>>max2>>gTimeCutMin>>gTimeCutMax>>gMultCutMin>>gMultCutMax>>gsTimeCuts>>gsTotMean;

  TObjArray *sarr = gsTimeCuts.Tokenize(";");

  std::cout<<"sarr->GetEntries()  "<< sarr->GetEntries() << "    rest: "<< gsTimeCuts<<std::endl;
  if(sarr->GetEntries()==1920){
    for (Int_t m=0; m <nmcp; m++) {
      for(Int_t p=0; p<npix; p++){
	TString cut = ((TObjString *) sarr->At(2*(m*npix+p)))->GetName();
	gTimeCuts[m][p][0] = cut.Atof();
	cut = ((TObjString *) sarr->At(2*(m*npix+p)+1))->GetName();
	gTimeCuts[m][p][1] = cut.Atof();
      }
    } 
  }
  sarr = gsTotMean.Tokenize(";");
  std::cout<<"sarr2->GetEntries()  "<< sarr->GetEntries() << "    rest: "<< gsTotMean <<std::endl;
  if(sarr->GetEntries()==960){
    for (Int_t m=0; m <nmcp; m++) {
      for(Int_t p=0; p<npix; p++){
	TString totmean = ((TObjString *) sarr->At(m*npix+p))->GetName();
	gTotMean[m][p] = totmean.Atof();
	std::cout<<"gTotMean[m][p]  "<<gTotMean[m][p] <<std::endl;
	
      }
    } 
  }
  
  for(Int_t i=0; i<maxMult; i++){
    hTotM[i]=new TH1F(Form("hTot%d",i),Form("hTot%d",i),500,min2,max2);
    hLeM[i]=new TH1F(Form("hLe%d",i),Form("hLe%d",i),500,-600,600);
    fOutput->Add(hTotM[i]);
    fOutput->Add(hLeM[i]);
  }

  hMultEvtNum1=new TH1F("hMultEvtNum1","",1000001,-0.5,1000000.5);
  hMultEvtNum2=new TH1F("hMultEvtNum2","",1000001,-0.5,1000000.5);
  axisTime800x500(hMultEvtNum1,"event number , [#]");
  hMultEvtNum1->GetYaxis()->SetTitle("multiplicity, [#]");
  axisTime800x500(hMultEvtNum2,"event number , [#]");
  hMultEvtNum2->GetYaxis()->SetTitle("multiplicity, [#]");
  fOutput->Add(hMultEvtNum1);
  fOutput->Add(hMultEvtNum2);

  for(Int_t m=0; m<nmcp; m++){
    fhDigi[m] = new TH2F( Form("mcp%d", m),Form("mcp%d", m),8,0.,8.,8,0.,8.);
    fhDigi[m]->SetStats(0);
    fhDigi[m]->SetTitle(0);
    fhDigi[m]->GetXaxis()->SetNdivisions(10);
    fhDigi[m]->GetYaxis()->SetNdivisions(10);
    fhDigi[m]->GetXaxis()->SetLabelOffset(100);
    fhDigi[m]->GetYaxis()->SetLabelOffset(100);
    fhDigi[m]->GetXaxis()->SetTickLength(1);
    fhDigi[m]->GetYaxis()->SetTickLength(1);
    fhDigi[m]->GetXaxis()->SetAxisColor(15);
    fhDigi[m]->GetYaxis()->SetAxisColor(15);
    fOutput->Add(fhDigi[m]);

    hEMult[m]   = new TH1F(Form("emult_m%d",m),Form("mcp %d",m),  500,0,500);
    axisTime800x500(hEMult[m],"multiplicity per event, [#]");
    fOutput->Add(hEMult[m]);

    for(Int_t p=0; p<npix; p++){     
      hPTime[m][p]   = new TH1F(Form("le_mcp%dpix%d",m,p),Form("mcp %d, pixel %d",m, p),  bins1,min1,max1);
      hPTot[m][p]   = new TH1F(Form("tot_mcp%dpix%d",m,p),Form("mcp %d, pixel %d",m, p),  bins2,min2,max2);
      hPMult[m][p]   = new TH1F(Form("mult_mcp%dpix%d",m,p),Form("mcp %d, pixel %d",m, p),  50,0,50);
      
      if(gMode==1){
	hShape[m][p] = new TH2F(Form("hShape_mcp%dpix%d",m,p), Form("hShape_%d_%d",m,p) , 400,80,110,120,-30,30);
	hLeTot[m][p] = new TH2F(Form("hLeTot_mcp%dpix%d" ,m,p), Form("mcp %d, pixel %d",m, p), 200,min1,max1, bins2,min2,max2);
	axisTime800x500(hShape[m][p],"time, [ns]");
	hShape[m][p]->GetYaxis()->SetTitle("offset to the threshold, [mV]");
      
	fOutput->Add(hLeTot[m][p]);
	fOutput->Add(hShape[m][p]);
      }

      axisTime800x500(hPTime[m][p]);
      axisTime800x500(hPTot[m][p],"TOT time, [ns]");
      axisTime800x500(hPMult[m][p],"multiplicity, [#]");
  
      fOutput->Add(hPTime[m][p]);
      fOutput->Add(hPTot[m][p]);
      fOutput->Add(hPMult[m][p]);
    }
  }

  hTot=new TH1F("hTotA","",500,min2,max2);
  hLe=new TH1F("hLeA","",500,-600,600);
  hMult=new TH1F("hMultA","",50,0,50);
  hCh=new TH1F("hChA","",3000,0,3000);

  axisTime800x500(hTot,"TOT time, [ns]");
  axisTime800x500(hLe,"LE time, [ns]");
  axisTime800x500(hMult,"multiplicity, [#]");
 
  gStyle->SetOptStat(1001111);
  
  fOutput->Add(hLe);
  fOutput->Add(hTot);
  fOutput->Add(hMult);
  fOutput->Add(hEMult[nmcp]);
  fOutput->Add(hCh);
}

Bool_t MSelector::Process(Long64_t entry){
  GetEntry(entry);

  Double_t offset=0;
  if(gMode==1){
    TString current_file_name  = MSelector::fChain->GetCurrentFile()->GetName();
    TObjArray *sarr = current_file_name.Tokenize("_");
    if(sarr->GetEntries()==3){
      if(((TObjString *) sarr->At(0))->GetName()=="th");
      TString soffset = ((TObjString *) sarr->At(1))->GetName();
      offset = soffset.Atof();
    }
  }
 
  Double_t le,tot, refLe=-1;
  TPrtHit hit;
  Int_t mcp,pix, mult,col,row,ch;
  Int_t thitCount1=0, thitCount2=0, hitCount1=0, hitCount2=0;
  if(gTrigger>0){
    for(UInt_t h=0; h<fEvent->GetHitsSize(); h++){
      mult = fEvent->GetMultiplicity(h);
      for(UInt_t m=0; m<mult; m++){
  	hit = fEvent->GetHit(h,m);
	ch  = hit.GetChannel();

	// bad pixels july14
	if(ch == 379) continue;
	if(ch == 381) continue;
	if(ch == 1397) continue;
	if(ch == 1869) continue;

	if(ch == 1405) continue;
	if(ch == 1403) continue;
	if(ch == 1385) continue;
	if(ch == 1381) continue;
	if(ch == 1383) continue;
	if(ch == 1387) continue;

	if(hit.GetMcpId()>14) thitCount1++;
	else  thitCount2++;
  	if(ch == (UInt_t)gTrigger) {
  	  refLe = hit.GetLeadTime();
  	}
      }
    }
  }else{
    refLe=0;
  }

  for(UInt_t h=0; h<fEvent->GetHitsSize(); h++){
    mult = fEvent->GetMultiplicity(h);
    hMult->Fill(mult);
    // if(refLe==-1 && gTrigger!=0) continue; 
    for(UInt_t m=0; m<mult; m++){
      hit = fEvent->GetHit(h,m);
      if(m==0) hCh->Fill(hit.GetChannel());
      mcp = hit.GetMcpId();

      if(mcp>14){
	hitCount1++;
	continue;
      }
      if(gMultCutMin!=gMultCutMax && (thitCount2<gMultCutMin || thitCount2>gMultCutMax)) continue; 

      ch  = hit.GetChannel();
      pix = hit.GetPixelId()-1;
      col = 7-pix/8;
      row = pix%8;
      le = hit.GetLeadTime();
      tot = hit.GetTotTime();   

      pix = col*8+row;
      
      // bad pixels july14
      if(ch == 379) continue;
      if(ch == 381) continue;
      if(ch == 1397) continue;
      if(ch == 1869) continue;

      if(ch == 1405) continue;
      if(ch == 1403) continue;
      if(ch == 1385) continue;
      if(ch == 1381) continue;
      if(ch == 1383) continue;
      if(ch == 1387) continue;
 
      Double_t timeDiff = le-refLe;
      
      if(gsTimeCuts!="0" && (timeDiff<gTimeCuts[mcp][col+8*row][0] || timeDiff>gTimeCuts[mcp][col+8*row][1])) continue;
      else if(gTimeCutMin!=gTimeCutMax &&  (timeDiff<gTimeCutMin || timeDiff>gTimeCutMax)) continue;
      hitCount2++;
	
      if(gsTotMean!="0"){
	timeDiff += 0.3*(tot - gTotMean[mcp][col+8*row]);
      }

      if(refLe!=-1 || gTrigger==0) {
	fhDigi[mcp]->Fill(col, row);
	if(gMode==1){
	  Int_t tpix = col+8*row;
	  hLeTot[mcp][tpix]->SetTitle(Form("ch %d",hit.GetChannel()));
	  hLeTot[mcp][tpix]->Fill(timeDiff,tot);
	  hShape[mcp][tpix]->Fill(timeDiff,offset);
	  hShape[mcp][tpix]->Fill(timeDiff + tot,offset);
	}
	hPTime[mcp][col+8*row]->Fill(timeDiff);
	hPTime[mcp][col+8*row]->SetTitle(Form("%d " ,hit.GetChannel()));
      }

      hPTot[mcp][col+8*row]->Fill(tot);
      hPTot[mcp][col+8*row]->SetTitle(Form("mcp %d, pixel %d, channel %d",mcp, pix, hit.GetChannel()));
      if(m==0) {
	hPMult[mcp][col+8*row]->Fill(fEvent->GetMultiplicity(h));
	hPMult[mcp][col+8*row]->SetTitle(Form("mcp %d, pixel %d, channel %d",mcp, pix, hit.GetChannel()));
      }

      if(m<(UInt_t)maxMult){
	hTotM[m]->Fill(tot);
	hLeM[m]->Fill(le);
      }
      hTot->Fill(tot);
      hLe->Fill(le);
    }
  }
  if(hitCount1+hitCount2 !=0) hEMult[0]->Fill(hitCount1+hitCount2);
  if(hitCount1!=0) hEMult[1]->Fill(hitCount1);
  if(hitCount2!=0) hEMult[2]->Fill(hitCount2);
  hMultEvtNum1->Fill(entry,hitCount1);
  hMultEvtNum2->Fill(entry,hitCount2);

  return kTRUE;
}

TF1 *gaust;
TVector3 fit(TH1F *h, Double_t range = 3){
  int binmax = h->GetMaximumBin();
  double xmax = h->GetXaxis()->GetBinCenter(binmax);
  gaust = new TF1("gaust","gaus(0)",xmax-range,xmax+range);
  Double_t integral = h->Integral(h->GetXaxis()->FindBin(xmax-0.6),h->GetXaxis()->FindBin(xmax+0.6));
  Double_t xxmin, xxmax, sigma1=0, mean1=0, sigma2, mean2;
  xxmax = xmax;
  xxmin = xxmax;
  Int_t nfound = 1, peakSearch = 1;
  if(integral>5){ 
    
    if(peakSearch == 1){
      gaust->SetParLimits(2,0.1,2);
      gaust->SetParameter(1,xmax);
      gaust->SetParameter(2,0.2);
    }
    
    if(peakSearch == 2){
      nfound = spect->Search(h,2,"",0.2);
      std::cout<<"nfound  "<<nfound <<std::endl;
      Float_t *xpeaks = spect->GetPositionX();
      if(nfound==1){
	gaust =new TF1("gaust","gaus(0)",xmax-range,xmax+range);
	gaust->SetParameter(1,xpeaks[0]);
      }else if(nfound==2) {
	if(xpeaks[0]>xpeaks[1]) {
	  xxmax = xpeaks[0];
	  xxmin = xpeaks[1];
	}else {
	  xxmax = xpeaks[1];
	  xxmin = xpeaks[0];
	}
	gaust =new TF1("gaust","gaus(0)+gaus(3)",xmax-range,xmax+range);
	gaust->SetParameter(1,xxmin);
	gaust->SetParameter(4,xxmax);
      }
    
      gaust->SetParameter(2,0.3);
      gaust->SetParameter(5,0.3);
    }

    h->Fit("gaust","","MQN",xxmin-range, xxmax+range);
    mean1 = gaust->GetParameter(1);
    sigma1 = gaust->GetParameter(2);
    if(sigma1>10) sigma1=10;
    
    if(peakSearch == 2){ 
      mean2 = (nfound==1) ? gaust->GetParameter(1) : gaust->GetParameter(4);
      sigma2 = (nfound==1) ? gaust->GetParameter(2) : gaust->GetParameter(5);
    }
  }
  delete gaust;
  return TVector3(mean1,sigma1,0);
}

void calculateTimeCut(){
  TVector3 res;
  gsTimeCuts="";
  for (Int_t m=0; m <nmcp; m++) {
    for(Int_t p=0; p<npix; p++){
      res = fit(hPTime[m][p]);
      gTimeCuts[m][p][0]=res.X() - 3*res.Y();
      gTimeCuts[m][p][1]=res.X() + 3*res.Y();
      gsTimeCuts+=Form("%f;%f;",gTimeCuts[m][p][0],gTimeCuts[m][p][1]);
    }
  }
}

void getTimeOffset(){
  std::cout<<"Creating calibration"<<std::endl;
  TH1D* h;
  TH2F* hh;
  for (Int_t m=0; m <nmcp; m++) {
    for(Int_t p=0; p<npix; p++){
      Double_t mean = fit(hPTime[m][p],0.5).X();
      hh =(TH2F*) hLeTot[m][p]->Clone("hh");
      hh->RebinY(2);

      // TF1 *gaust = new TF1("gaust","gaus(0)",85,105);
      // gaust->SetParameter(1,90);
      // gaust->SetParameter(2,0.3);
      // hh->FitSlicesX(gaust,0,-1,10,"");
      // TGraph * gg = new TGraph((TH1D*)gDirectory->Get("hh_1")); 
      Int_t ch = chmap[m][p];
      gGrDiff[ch] = new TGraph();
      for (int i=0;i<100;i++){
	Double_t x = hh->GetYaxis()->GetBinCenter(i);
	h = hh->ProjectionX(Form("bin%d",i+1),i+1,i+2);
	Double_t vx = fit((TH1F*)h,0.5).X();
	if(vx==0) vx = mean;
	gGrDiff[ch]->SetPoint(i,x,vx);
      }

      gGrDiff[ch]->SetName(Form("gCalib_ch%d",ch));
      gGrDiff[ch]->GetXaxis()->SetTitle("fine bin, [#]");
      gGrDiff[ch]->GetYaxis()->SetTitle("fine time, [ns]");
    
    }
  }
}

void MyMainFrame::DoExportOffsets(){
  if(gMode==1) {
    getTimeOffset();

    TString filedir=ginFile;
    filedir.Remove(filedir.Last('/'));
    TFile efile(filedir+"/calibOffsets.root","RECREATE");
    Int_t c;
    for (Int_t m=0; m <nmcp; m++) {
      for(Int_t p=0; p<npix; p++){
	c = chmap[m][p];
	gGrDiff[c]->SetName(Form("%d_%d_%d",c,m,p));
	gGrDiff[c]->Write();
      }
    }
    efile.Write();
    efile.Close();
    std::cout<<"Exporting .. Done"<<std::endl;
  }else{
    std::cout<<"For exporting use -a1 flag"<<std::endl;
  }
}

Bool_t lock = false;
void exec3event(Int_t event, Int_t gx, Int_t gy, TObject *selected){
  if(gComboId==0 || gComboId==2 || gComboId==5 || gComboId==4 || gComboId==10 || gComboId==11){
    TCanvas *c = (TCanvas *) gTQSender;
    TPad *pad = (TPad *) c->GetSelectedPad();
    if (!pad) return;
    Float_t x = pad->AbsPixeltoX(gx);
    Float_t y = pad->AbsPixeltoY(gy);
    x = pad->PadtoX(x);
    y = pad->PadtoY(y);
    if(event ==1 && lock) lock = false;
    else if(event ==1) lock = true;
    if(lock) return;

    if (selected->InheritsFrom(TH2::Class())){
      TH2F *hDigi = (TH2F *) selected;
      Int_t binx = hDigi->GetXaxis()->FindBin(x);
      Int_t biny = hDigi->GetYaxis()->FindBin(y);
      TString smcp = selected->GetName();
      smcp = smcp(3,smcp.Sizeof());
      Int_t mcp = smcp.Atoi();
      Int_t pix = 8*(biny-1)+binx-1;
      //    printf("Canvas %s: event=%d, x=%d, y=%d, p=%d, selected=%d\n", smcp.Data(), event, binx, biny, pix,smcp.Atoi());
      cTime->cd();
      if(gComboId==0) {
	hPTime[mcp][pix]->Draw();
	fit(hPTime[mcp][pix],1);
	hPTime[mcp][pix]->Draw("same");
      }
      if(gComboId==2) hPTot[mcp][pix]->Draw();   
      if(gComboId==5) hPMult[mcp][pix]->Draw();      
      if(gComboId==4) hLeTot[mcp][pix]->Draw("colz");
      if(gComboId==10) hShape[mcp][pix]->Draw("colz");
      if(gComboId==11){
	Int_t ch = chmap[mcp][pix];
	hLeTot[mcp][pix]->Draw("colz");
	Double_t* xx = gGrDiff[ch]->GetX();
	Double_t* yy = gGrDiff[ch]->GetY();

	TGraph* gr = new TGraph(gGrDiff[ch]->GetN(),yy,xx);
	gr->SetMarkerStyle(7);
	gr->SetMarkerColor(2);
	gr->Draw("P same");
      }
      if(gMain->fCheckBtn2->GetState() == kButtonDown){
	gMain->fEdit3->SetText(Form("%2.2f %2.2f", gTimeCuts[mcp][pix][0], gTimeCuts[mcp][pix][1]));
      }
      cTime->Update(); 
    }
  }
}

void MSelector::Terminate(){
  for (Int_t m=0; m <nmcp; m++) {
    fhDigi[m] = dynamic_cast<TH2F *>(TProof::GetOutput(Form("mcp%d",m), fOutput));
    hEMult[m] = dynamic_cast<TH1F *>(TProof::GetOutput(Form("emult_m%d",m), fOutput)); 
    for(Int_t p=0; p<npix; p++){
      hPTime[m][p] = dynamic_cast<TH1F *>(TProof::GetOutput(Form("le_mcp%dpix%d",m,p), fOutput)); 
      hPTot[m][p] = dynamic_cast<TH1F *>(TProof::GetOutput(Form("tot_mcp%dpix%d",m,p), fOutput)); 
      hPMult[m][p] = dynamic_cast<TH1F *>(TProof::GetOutput(Form("mult_mcp%dpix%d",m,p), fOutput)); 
      if(gMode==1){
	hLeTot[m][p] = dynamic_cast<TH2F *>(TProof::GetOutput(Form("hLeTot_mcp%dpix%d",m,p), fOutput)); 
	hShape[m][p] = dynamic_cast<TH2F *>(TProof::GetOutput(Form("hShape_mcp%dpix%d",m,p), fOutput));
      }
    }
  }
  for(Int_t i=0; i<maxMult; i++){
    hLeM[i] = dynamic_cast<TH1F *>(TProof::GetOutput(Form("hLe%d",i), fOutput)); 
    hTotM[i] = dynamic_cast<TH1F *>(TProof::GetOutput(Form("hTot%d",i), fOutput)); 
  }

  hLe = dynamic_cast<TH1F *>(TProof::GetOutput("hLeA", fOutput));   
  hTot = dynamic_cast<TH1F *>(TProof::GetOutput("hTotA", fOutput)); 
  hMult = dynamic_cast<TH1F *>(TProof::GetOutput("hMultA", fOutput)); 
  hEMult[nmcp] = dynamic_cast<TH1F *>(TProof::GetOutput(Form("emult_m%d",nmcp), fOutput)); 
  hMultEvtNum1 = dynamic_cast<TH1F *>(TProof::GetOutput("hMultEvtNum1", fOutput));
  hMultEvtNum2 = dynamic_cast<TH1F *>(TProof::GetOutput("hMultEvtNum2", fOutput));

  hCh = dynamic_cast<TH1F *>(TProof::GetOutput("hChA", fOutput)); 
}

void MyMainFrame::DoDraw(){
  fCheckBtn2->SetText("3 sigma time cut                         ");
  fCheckBtn2->SetTextColor(Pixel_t(0x000000),kFALSE);   
  fCheckBtn3->SetText("Walk correction                         ");
  fCheckBtn3->SetTextColor(Pixel_t(0x000000),kFALSE);  
  if(gTrigger>-1) gTrigger = fNumber->GetIntNumber();

  for(Int_t m=0; m<nmcp; m++){
    if(fhDigi[m]) fhDigi[m]->Reset();
    if(hEMult[m]) hEMult[m]->Reset();
    for(Int_t p=0; p<npix; p++){
      if(hPTime[m][p]) hPTime[m][p]->Reset();
      if(hPTot[m][p]) hPTot[m][p]->Reset();
      if(hPMult[m][p]) hPMult[m][p]->Reset();
    }
  }
  if(hLe) hLe->Reset();
  if(hTot) hTot->Reset();
  if(hMult) hMult->Reset();
  if(hCh) hCh->Reset(); 

  fHProg3->Reset();
  //fNEntries = 100000;

  TString option = Form("%d %d %s %s %s %s %s %s",gMode,gTrigger,fEdit1->GetText(),fEdit2->GetText(),fEdit3->GetText(),fEdit4->GetText(),gsTimeCuts.Data(), gsTotMean.Data());

  gROOT->SetBatch(1);
  fCh->Process(fSelector,option,fNEntries);
  gROOT->SetBatch(0);
 
  if(gTrigger<0) gTrigger=0;

  Int_t tmax, max=0;
  for(Int_t p=0; p<15;p++){
    tmax = fhDigi[p]->GetMaximum();
    if(max<tmax) max = tmax;
  }
  max = (max>0)? max : 1;
  fHslider1->SetRange(0,max);
  fHslider1->SetPosition(max);

  fCheckBtn1->SetState(kButtonUp);

  drawDigi("m,p,v\n",1);
  updatePlot(gComboId);
}

TString MyMainFrame::updatePlot(Int_t id, TCanvas *cT){
  if(!cT) cT = cTime;
  TString histname="";
  gComboId = id;
  cT->cd();
  cT->SetLogy(0);
  Int_t max = 0;
  TLegend *leg;
  switch(id) {
  case 0: // LE 
    break;
  case 1: // LE All
    cT->SetLogy(1);
    hLe->Draw();
    for(Int_t i=0; i<maxMult; i++){
      Int_t colorid = i+1;
      if(colorid ==5) colorid=9;
      if(colorid ==3) colorid=8;
      if(hTotM[i]->GetEntries()>0){
	hLeM[i]->SetLineColor(colorid);
	if(i==0) hLeM[i]->SetLineWidth(2);
	hLeM[i]->Draw("same");
      }
    }
    histname=hLe->GetName();
    break;
  case 2: // ToT
    break;
  case 3: // ToT All
    hTot->Draw();
    for(Int_t i=0; i<maxMult; i++){
      Int_t colorid = i+1;
      if(colorid ==5) colorid=9;
      if(colorid ==3) colorid=8;
      if(hTotM[i]->GetEntries()>0){
	hTotM[i]->SetLineColor(colorid);
	if(i==0) hTotM[i]->SetLineWidth(2);
	hTotM[i]->Draw("same");
      }
    } 
    histname=hTot->GetName();
    break;
  case 4: // Le vs. ToT
    break;
  case 5: // Mult
    break;
  case 6: // Mult All
    cT->SetLogy(1);
    hMult->Draw();
    histname=hMult->GetName();
    break;
  case 7: // Channels
    hCh->Draw();
    histname=hCh->GetName();
    break; 
  case 8: // Event Mult
    if(hEMult[0]->GetMaximum()>max) max = hEMult[0]->GetMaximum();
    if(hEMult[1]->GetMaximum()>max) max = hEMult[1]->GetMaximum();
    if(hEMult[2]->GetMaximum()>max) max = hEMult[2]->GetMaximum();
    
    hEMult[0]->SetMaximum(max+0.05*max);
    hEMult[0]->Draw();
    hEMult[1]->SetLineColor(2);
    hEMult[1]->Draw("same");
    hEMult[2]->SetLineColor(4);
    hEMult[2]->Draw("same");
    leg = new TLegend(0.5,0.7,0.9,0.9);
    leg->SetFillColor(0);
    leg->AddEntry(hEMult[0],"Total","l");
    leg->AddEntry(hEMult[2],"Dirc MCPs","l");
    leg->AddEntry(hEMult[1],"Rest","l");
    leg->Draw();
    histname=hEMult[0]->GetName();
    break;
  case 9: // Mult vs. Event number
    if(hMultEvtNum1->GetMaximum()>max) max = hMultEvtNum1->GetMaximum();
    if(hMultEvtNum2->GetMaximum()>max) max = hMultEvtNum2->GetMaximum();
    
    hMultEvtNum1->SetMaximum(max+0.05*max);
    hMultEvtNum1->SetMinimum(0);
    hMultEvtNum1->SetLineColor(2);
    hMultEvtNum1->Draw();
    hMultEvtNum2->SetLineColor(4);
    hMultEvtNum2->Draw("same");
    leg = new TLegend(0.5,0.7,0.9,0.9);
    leg->SetFillColor(0);
    leg->AddEntry(hMultEvtNum2,"Dirc MCPs","l");
    leg->AddEntry(hMultEvtNum1,"Rest","l");
    leg->Draw();
    histname=hMultEvtNum1->GetName();
    break;
  }
  cT->Modified();
  cT->Update();
  return histname;
}

void MyMainFrame::DoMore(){
  if(GetState(fHm)) {
    HideFrame(fHm);
    fBtnMore->SetText("&More");
  }
  else{
    ShowFrame(fHm);
    fBtnMore->SetText("&Less");
  }
}

void MyMainFrame::DoExport(){
  gROOT->SetBatch(1);
  TCanvas *cExport = new TCanvas("cExport","cExport",0,0,800,400);
  cExport->SetCanvasSize(800,400);
  Int_t saveFlag = 1;
  TString histname="", filedir=ginFile;
  filedir.Remove(filedir.Last('/'));
  TString path = createDir(filedir+"/plots", "beam_data "+ginFile, saveFlag); 
  std::cout<<"Exporting into  "<<path <<std::endl;
  writeInfo("digi.csv", drawDigi("m,p,v\n",1), saveFlag);
  pbar->Reset();
  Float_t total = (nmcp-1)*(npix-1);
  if(gComboId==0 || gComboId==2 || gComboId==5 || gComboId==4 || gComboId==10 || gComboId==11){
    for(Int_t m=0; m<nmcp; m++){
      for(Int_t p=0; p<npix; p++){
	cExport->cd();
	if(gComboId==0) {
	  hPTime[m][p]->Draw();
	  fit(hPTime[m][p]);
	  hPTime[m][p]->Draw("same");
	  histname=hPTime[m][p]->GetName();
	}
	if(gComboId==2){
	  hPTot[m][p]->Draw(); 
	  histname=hPTot[m][p]->GetName();
	}  
	if(gComboId==5){
	  hPMult[m][p]->Draw();  
	  histname=hPMult[m][p]->GetName();    
	}
	if(gComboId==4){
	  hLeTot[m][p]->Draw("colz");
	  histname=hLeTot[m][p]->GetName();
	}
	if(gComboId==10){
	  hShape[m][p]->Draw("colz");
	  histname=hShape[m][p]->GetName();
	}
	save(cExport,path,histname,"mdisplay",saveFlag,1);
	pbar->SetPosition(100*(m*p)/total);
	gSystem->ProcessEvents();
      }
    }
  }else{
    histname = updatePlot(gComboId,cExport);
    save(cExport,path,histname,"mdisplay",saveFlag,1);
  }
  cExport = (TCanvas *) cDigi->DrawClone();
  cExport->SetCanvasSize(800,400);
  save(cExport,path,"digi","mdisplay",saveFlag,1);
  gROOT->SetBatch(0);
  std::cout<<"Exporting .. Done"<<std::endl;
  
}

TLine *gLine = new TLine(0,0,3000,0);
void MyMainFrame::DoSlider(Int_t pos){
  if(fCheckBtn1->GetState() != kButtonDown)
    drawDigi("m,p,v\n",1,pos);

  if(gComboId==7){
    cTime->cd();  
    gLine->SetY1(pos);
    gLine->SetY2(pos);
    gLine->SetLineColor(kRed);
    gLine->Draw();
    cTime->Update();
  }
  
}

void MyMainFrame::DoCheckBtnClecked1(){
  Int_t state = -1;
  if(fCheckBtn1->GetState() != kButtonDown){
    state = fHslider1->GetPosition();
    gMain->fHslider1->SetEnabled(kTRUE);
  }else{
    gMain->fHslider1->SetEnabled(kFALSE);
  }
  drawDigi("m,p,v\n",1,state);
}

void MyMainFrame::DoCheckBtnClecked2(){
  Int_t state = -1;
  fCheckBtn2->SetTextColor(Pixel_t(0xff0000),kFALSE);      
  fCheckBtn2->SetText("3 sigma time cut (Press Evaluate)");
  if(fCheckBtn2->GetState() == kButtonDown){
    calculateTimeCut();    
    gMain->fEdit3->SetEnabled(kFALSE);
  }
  if(fCheckBtn2->GetState() == kButtonUp){  
    gsTimeCuts="0";
    gMain->fEdit3->SetText("0 0");
    gMain->fEdit3->SetEnabled(kTRUE);
    for (Int_t m=0; m <nmcp; m++) {
      for(Int_t p=0; p<npix; p++){
	gTimeCuts[m][p][0]=0;
	gTimeCuts[m][p][1]=0;
      }
    }
  }
}

void MyMainFrame::DoCheckBtnClecked3(){
  Int_t state = -1;
  fCheckBtn3->SetTextColor(Pixel_t(0xff0000),kFALSE);       
  fCheckBtn3->SetText("Walk correction (Press Evaluate)");
  if(fCheckBtn3->GetState() == kButtonDown){
    // get mean TOT for each pixel
    gsTotMean="";
    TVector3 res;
    for (Int_t m=0; m <nmcp; m++) {
      for(Int_t p=0; p<npix; p++){
	res = fit(hPTot[m][p]);
	gsTotMean+=Form("%f;",res.X());
      }
    }  
  }
  if(fCheckBtn3->GetState() == kButtonUp){      
    gsTotMean="0";
  }
}

TH2F* fhDigi_temp[15];
void MyMainFrame::DoCheckBtnClecked4(){
  if(fCheckBtn4->GetState() == kButtonDown){
    TVector3 res;
    for(Int_t m=0; m<nmcp; m++){
      fhDigi_temp[m] = (TH2F*)fhDigi[m]->Clone();
      if(fhDigi[m]) fhDigi[m]->Reset();
    }
    for (Int_t m=0; m <nmcp; m++) {
      for(Int_t p=0; p<npix; p++){
	Int_t col = p/8;
	Int_t row = p%8;
	Double_t mean = fit(hPTime[m][p],0.5).X();
	if(mean>90) mean = 90; 
	fhDigi[m]->Fill(row,col,mean);
      }
    }
    drawDigi("m,p,v\n",1,-2,-2);
  }
  if(fCheckBtn4->GetState() == kButtonUp){
    for(Int_t m=0; m<nmcp; m++){
      fhDigi[m] = fhDigi_temp[m];
    }
    drawDigi("m,p,v\n",1);
  }
}

void MyMainFrame::DoExit(){
  gApplication->Terminate(0);
}

void MyMainFrame::SetStatusText(const char *txt, Int_t pi){
  fStatusBar->SetText(txt,pi);
}

void MyMainFrame::EventInfo(Int_t event, Int_t px, Int_t py, TObject *selected){
  const char *text0, *text1, *text3;
  char text2[50];
  text0 = selected->GetTitle();
  SetStatusText(text0,0);
  text1 = selected->GetName();
  SetStatusText(text1,1);
  if (event == kKeyPress)
    sprintf(text2, "%c", (char) px);
  else
    sprintf(text2, "%d,%d", px, py);
  SetStatusText(text2,2);
  text3 = selected->GetObjectInfo(px,py);
  SetStatusText(text3,3);
}

MyMainFrame::MyMainFrame(const TGWindow *p, UInt_t w, UInt_t h) : TGMainFrame(p, w, h){

  // Create the embedded canvas
  fEcan = new TRootEmbeddedCanvas(0,this,800,350);
  Int_t wid0 = fEcan->GetCanvasWindowId();
  cDigi = new TCanvas("cDigi",10,10,wid0);
  cDigi->SetMargin(0,0,0,0);
  fEcan->AdoptCanvas(cDigi);

  fTime = new TRootEmbeddedCanvas(0,this,800,350);
  Int_t wid1 = fTime->GetCanvasWindowId();
  cTime = new TCanvas("cTime",10,10,wid1);
  fTime->AdoptCanvas(cTime);

  AddFrame(fEcan, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 0, 0, 0, 0));
  AddFrame(fTime, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 0, 0, 0, 0));

  // =============== More frame ===================
  fHm = new TGVerticalFrame(this, 200, 120, kFixedSize);

  TGHorizontalFrame *fHm0 = new TGHorizontalFrame(fHm, 400, 40);
  TGLabel *fL1 = new TGLabel(fHm0, "LE's histogram parameters: ");
  fEdit1 = new TGTextEntry(fHm0, new TGTextBuffer(100));
  fEdit1->SetToolTipText("bins min max");
  fEdit1->Resize(80, fEdit1->GetDefaultHeight());
  fHm0->AddFrame(fL1, new TGLayoutHints(kLHintsTop | kLHintsLeft,5, 5, 5, 5));
  fHm0->AddFrame(fEdit1, new TGLayoutHints(kLHintsTop | kLHintsLeft,5, 5, 5, 5));

  TGLabel *fL2 = new TGLabel(fHm0, "Tot's histogram parameters: ");
  fEdit2 = new TGTextEntry(fHm0, new TGTextBuffer(100));
  fEdit2->SetToolTipText("bins min max");
  fEdit2->Resize(80, fEdit2->GetDefaultHeight());
  fHm0->AddFrame(fL2, new TGLayoutHints(kLHintsTop | kLHintsLeft,5, 5, 5, 5));
  fHm0->AddFrame(fEdit2, new TGLayoutHints(kLHintsTop | kLHintsLeft,5, 5, 5, 5));

  TGLabel *fL4 = new TGLabel(fHm0, "Time cut: ");
  fEdit3 = new TGTextEntry(fHm0, new TGTextBuffer(100));
  fEdit3->SetToolTipText("min max");
  fEdit3->Resize(80, fEdit3->GetDefaultHeight());
  fHm0->AddFrame(fL4, new TGLayoutHints(kLHintsTop | kLHintsLeft,5, 5, 5, 5));
  fHm0->AddFrame(fEdit3, new TGLayoutHints(kLHintsTop | kLHintsLeft,5, 5, 5, 5));


  TGLabel *fL5 = new TGLabel(fHm0, "Multiplicity cut: ");
  fEdit4 = new TGTextEntry(fHm0, new TGTextBuffer(100));
  fEdit4->SetToolTipText("min max");
  fEdit4->Resize(80, fEdit4->GetDefaultHeight());
  fHm0->AddFrame(fL5, new TGLayoutHints(kLHintsTop | kLHintsLeft,5, 5, 5, 5));
  fHm0->AddFrame(fEdit4, new TGLayoutHints(kLHintsTop | kLHintsLeft,5, 5, 5, 5));

  fHm->AddFrame(fHm0, new TGLayoutHints(kLHintsExpandX | kLHintsCenterX,5, 5, 5, 5));

  TGHorizontalFrame *fHm1 = new TGHorizontalFrame(fHm, 400, 40);

  TGLabel *fL3 = new TGLabel(fHm1, "Maximum in hit pattern: ");
  fHm1->AddFrame(fL3, new TGLayoutHints(kLHintsBottom | kLHintsLeft,5, 5, 5, 5));
  fHslider1 = new TGHSlider(fHm1, 200, kSlider1 | kScaleBoth, 0);
  fHslider1->Connect("PositionChanged(Int_t)", "MyMainFrame", this, "DoSlider(Int_t)");
  fHslider1->SetRange(0,50);
  fHm1->AddFrame(fHslider1, new TGLayoutHints(kLHintsBottom | kLHintsLeft,5, 5, 5, 5));

  fCheckBtn1  = new TGCheckButton(fHm1, new TGHotString("Local maximum"),        -1);
  fCheckBtn1->Connect("Clicked()", "MyMainFrame", this, "DoCheckBtnClecked1()");
  fHm1->AddFrame(fCheckBtn1, new TGLayoutHints(kLHintsBottom | kLHintsLeft,5, 5, 5, 5));

  fCheckBtn2  = new TGCheckButton(fHm1, new TGHotString("3 sigma time cut                         "),        -1);
  fCheckBtn2->Connect("Clicked()", "MyMainFrame", this, "DoCheckBtnClecked2()");
  fHm1->AddFrame(fCheckBtn2, new TGLayoutHints(kLHintsBottom | kLHintsLeft,5, 5, 5, 5));

  fCheckBtn4  = new TGCheckButton(fHm1, new TGHotString("show time map"),-1);
  fCheckBtn4->Connect("Clicked()", "MyMainFrame", this, "DoCheckBtnClecked4()");
  fHm1->AddFrame(fCheckBtn4, new TGLayoutHints(kLHintsBottom | kLHintsLeft,5, 5, 5, 5));

  TGTextButton * fBtnExport = new TGTextButton(fHm1, "Export for the &blog");
  fBtnExport->Connect("Clicked()", "MyMainFrame", this, "DoExport()");
  fHm1->AddFrame(fBtnExport, new TGLayoutHints(kLHintsBottom | kLHintsLeft,5, 5, 5, 5));

  fHm->AddFrame(fHm1, new TGLayoutHints(kLHintsExpandX | kLHintsCenterX,5, 5, 5, 5));

  TGHorizontalFrame *fHm2 = new TGHorizontalFrame(fHm, 400, 40);

  fCheckBtn3  = new TGCheckButton(fHm2, new TGHotString("Walk correction                         "),        -1);
  fCheckBtn3->Connect("Clicked()", "MyMainFrame", this, "DoCheckBtnClecked3()");
  fHm2->AddFrame(fCheckBtn3, new TGLayoutHints(kLHintsBottom | kLHintsLeft,5, 5, 5, 5));

  TGTextButton * fBtnExportOffsets = new TGTextButton(fHm2, "Export &offsets");
  fBtnExportOffsets->Connect("Clicked()", "MyMainFrame", this, "DoExportOffsets()");
  fHm2->AddFrame(fBtnExportOffsets, new TGLayoutHints(kLHintsBottom | kLHintsLeft,5, 5, 5, 5));

  fHm->AddFrame(fHm2, new TGLayoutHints(kLHintsExpandX | kLHintsCenterX,5, 5, 5, 5));


  AddFrame(fHm,   new TGLayoutHints(kLHintsExpandX | kLHintsCenterX, 2, 2, 2, 2));
  // =============== More frame ===================


  // status bar
  //Int_t parts[] = {45, 15, 10, 30};
  //fStatusBar = new TGStatusBar(this, 50, 10, kVerticalFrame);
  //fStatusBar->SetParts(parts, 4);
  //fStatusBar->Draw3DCorner(kFALSE);
  //AddFrame(fStatusBar, new TGLayoutHints(kLHintsExpandX, 0, 0, 10, 0));
   
  // Create a horizontal frame containing two buttons
  TGHorizontalFrame *hframe = new TGHorizontalFrame(this, 400, 40);

  fLabel = new TGLabel(hframe, "Refference channel: ");
  hframe->AddFrame(fLabel, new TGLayoutHints(kLHintsTop | kLHintsLeft,5, 5, 3, 4));

  fNumber = new TGNumberEntry(hframe, 0, 9,999, TGNumberFormat::kNESInteger,
			      TGNumberFormat::kNEANonNegative, 
			      TGNumberFormat::kNELLimitMinMax,
			      0, 99999);
  fNumber->SetNumber(0);
  hframe->AddFrame(fNumber, new TGLayoutHints(kLHintsTop | kLHintsLeft, 5, 5, 3, 4));
   
  TGTextButton *draw = new TGTextButton(hframe, "&Evaluate");
  draw->Connect("Clicked()", "MyMainFrame", this, "DoDraw()");
  hframe->AddFrame(draw, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));
   
  fHProg3 = new TGHProgressBar(hframe, TGProgressBar::kStandard, 100);
  fHProg3->SetFillType(TGProgressBar::kBlockFill);
  hframe->AddFrame(fHProg3, new TGLayoutHints(kLHintsTop | kLHintsLeft,5, 5, 5, 5));
  pbar = fHProg3;


  fLabel1 = new TGLabel(hframe, "   Show: ");
  hframe->AddFrame(fLabel1, new TGLayoutHints(kLHintsTop | kLHintsLeft,5, 5, 5, 5));


  fComboMode = new TGComboBox(hframe);
  hframe->AddFrame(fComboMode, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 5, 5));
  fComboMode->AddEntry("Leading Edge", 0);
  fComboMode->AddEntry("Leading Edge All", 1);
  fComboMode->AddEntry("ToT", 2);
  fComboMode->AddEntry("ToT All", 3);
  if(gMode==1){
    fComboMode->AddEntry("Le vs. ToT", 4);
    fComboMode->AddEntry("Signal shape",10);
    fComboMode->AddEntry("Offset graph",11);
  }
  fComboMode->AddEntry("Multiplicity",5);
  fComboMode->AddEntry("Multiplicity All",6);
  fComboMode->AddEntry("Event Mult",8);
  fComboMode->AddEntry("Mult vs. Evt#",9);
  fComboMode->AddEntry("Channels",7);

  fComboMode->Resize(300, 20);
  fComboMode->Connect("Selected(Int_t)", "MyMainFrame", this, "updatePlot(Int_t)");
  

  TGTextButton *exit = new TGTextButton(hframe, "&Exit ");
  exit->Connect("Pressed()", "MyMainFrame", this, "DoExit()");
  hframe->AddFrame(exit, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));

  fBtnMore = new TGTextButton(hframe, "&More ");
  fBtnMore->Connect("Pressed()", "MyMainFrame", this, "DoMore()");
  hframe->AddFrame(fBtnMore, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));

  AddFrame(hframe, new TGLayoutHints(kLHintsExpandX | kLHintsCenterX, 2, 2, 2, 2));

  // Set a name to the main frame   
  SetWindowName("mdisplay " + ginFile);
  MapSubwindows();

  // Initialize the layout algorithm via Resize()
  Resize(GetDefaultSize());
  fComboMode->Resize(300, 20);
  HideFrame(fHm);

  //if(gTrigger==0) fEdit1->SetText("400 -200 200");
  if(gTrigger==0) fEdit1->SetText("400 85 105");
  else if(gTrigger==1952 || gTrigger==1956 || gTrigger==1953 || gTrigger==1957) fEdit1->SetText("400 80 120");
  else if(gTrigger==1920 || gTrigger==1921) fEdit1->SetText("400 -100 -50");
  else if(gTrigger==2560 || gTrigger==2561) fEdit1->SetText("400 150 200");
  else fEdit1->SetText("300 0 60");
  if(ginFile.Contains("C.root"))  fEdit1->SetText("400 50 100");;

  fEdit2->SetText("200 -2 5");
  if(gTrigger==1952 || gTrigger==1956) fEdit2->SetText("200 -2 5");
  fEdit3->SetText("0 0");
  fEdit4->SetText("0 0");
  fNumber->SetIntNumber(gTrigger);
  fCheckBtn2->SetState(kButtonUp);

  init();
  DoDraw();
  fComboMode->Select(7);

  cDigi->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)", 0, 0,
		 "exec3event(Int_t,Int_t,Int_t,TObject*)");
  MapWindow();  // Map main frame
}

MyMainFrame::~MyMainFrame(){
  Cleanup();
  delete fEcan;
  delete fTime;
}

void mdisplay(TString inFile= "pilasM.root", Int_t trigger=0, Int_t mode=0){
  ginFile = inFile;
  gTrigger= trigger;
  gMode=mode;
  gMain = new MyMainFrame(gClient->GetRoot(), 800, 800);
}
