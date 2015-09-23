// prttools - useful functions for hld* 
// original author: Roman Dzhygadlo - GSI Darmstadt 

#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TVector3.h"
#include "TMath.h"
#include "TChain.h"
#include "TGaxis.h"
#include "TColor.h"
#include "TString.h"
#include "TArrayD.h"
#include "TSpectrum.h"
#include "TSpectrum2.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TRandom2.h"
#include "TError.h"
#include "TPaveStats.h"
#include "TObjString.h"
#include "TApplication.h"
#include <TLegend.h>
#include <TAxis.h>

#include <iostream>
#include <fstream>
#include <sstream>

#ifdef prt__sim
   class PrtEvent;
   class PrtHit;
   PrtEvent* fEvent = 0;
#endif 

#ifdef prt__beam
   class PrtEvent;
   class PrtHit;
   PrtEvent* fEvent = 0;
#endif  


TChain*  fCh = 0;
Int_t    fNEntries(0), fMomentum(0),fAngle(0),fParticle(0),fTest1(0),fTest2(0);
TString  fSavePath = "";
TString  fInfo = "", fPath;
TH2F*    fhDigi[15];
TPad*    fhPads[15];
TPad*    fhPglobal;
TCanvas* cDigi;

const Int_t nmcp = 15, npix = 64;
const Int_t maxmch(nmcp*npix);
const Int_t maxch =1500;
const Int_t tdcmax=10000;
const Int_t ctdc = 48; //41

Int_t map_tdc[tdcmax];
Int_t map_mpc[nmcp][npix];
Int_t map_mcp[maxch];
Int_t map_pix[maxch];
Int_t map_row[maxch];
Int_t map_col[maxch];

// const Int_t tdcnum=16;
// TString tdcsid[tdcnum] ={"10","11","12","13",
// 			 "20","21","22","23",
// 			 "780","781","782","783",
// 			 "840","841","842","843"
// };

// const Int_t tdcnum=41; //may2015
// TString tdcsid[tdcnum] ={"2000","2001","2002","2003","2004","2005","2006","2007","2008","2009",
// 			 "200a","200b","200c","200d","200e","200f","2010","2011","2012","2013",
// 			 "2014","2015","2016","2018","2019","201a","201c","2020","2023","2024",
// 			 "2025","2026","2027","2028","2029","202a","202b","202c","202d","202e","202f"
// };

const Int_t tdcnum=30;  //jun2015
TString tdcsid[tdcnum] ={"2000","2001","2002","2003","2004","2005","2006","2007","2008","2009",
			 "200a","200b","200c","200d","200e","200f","2010","2011","2012","2013",
			 "2014","2015","2016","2018","2019","201a","201c","201d","202c","202d"
};

void CreateMap(){
  Int_t seqid =-1;
  for(Int_t i=0; i<tdcmax; i++) map_tdc[i]=-1;
  for(Int_t i=0; i<tdcnum; i++){
    Int_t dec = TString::BaseConvert(tdcsid[i],16,10).Atoi();
    map_tdc[dec]=++seqid;
  }
  
  for(Int_t ch=0; ch<maxmch; ch++){
    Int_t mcp = ch/64;
    Int_t pix = ch%64;	
    Int_t col = pix/2 - 8*(pix/16);
    Int_t row = pix%2 + 2*(pix/16);
    pix = col+8*row;
      
    map_mpc[mcp][pix]=ch;
    map_mcp[ch] = mcp;
    map_pix[ch] = pix;
    map_row[ch] = row;
    map_col[ch] = col;
  }
}

Int_t GetChannelNumber(Int_t tdc, Int_t tdcChannel){
  Int_t ch = -1;
  if(tdc>=0) ch = 48*tdc+tdcChannel;
  return ch;
}

Int_t RemoveRefChannels(Int_t ch, Int_t tdcSeqId){
  return ch - tdcSeqId;
}

Int_t AddRefChannels(Int_t ch,Int_t tdcSeqId){
  return ch + tdcSeqId;
}

Bool_t badcannel(Int_t ch){
  if(ch<0) return true;
  
  // bad pixels july15

  if(ch==202) return true;
  if(ch==204) return true;
  if(ch==206) return true;
  if(ch==830) return true;
  if(ch==831) return true;
  if(ch==828) return true;
  if(ch==815) return true;
	
  return false;
}

TString drawDigi(TString digidata="", Int_t layoutId = 0, Double_t maxz = 0, Double_t minz = 0){
  if(!cDigi) cDigi = new TCanvas("cDigi","cDigi",0,0,800,400);
  cDigi->cd();
  // TPad * pp = new TPad("P","T",0.06,0.135,0.93,0.865);
  if(!fhPglobal){
    Double_t tt =(layoutId==3)? 0.88: 0.96; 
    fhPglobal = new TPad("P","T",0.04,0.04,tt,0.96);
    fhPglobal->SetFillStyle(0);
    fhPglobal->Draw();
  }
  fhPglobal->cd();
  
  Int_t nrow = 3, ncol = 5;
 
  if(layoutId > 1){
    float bw = 0.02, bh = 0.01, shift = 0,shiftw=0.02;
    float tbw = bw, tbh = bh;
    Int_t padi = 0;
    if(!fhPads[0]){
      for(int ii=0; ii<ncol; ii++){
	for(int j=0; j<nrow; j++){
	  if(j==1) shift = -0.028;
	  else shift = 0;
	  fhPads[padi] =  new TPad(Form("P%Double_t",ii*10+j),"T", ii/(Double_t)ncol+tbw+shift+shiftw , j/(Double_t)nrow+tbh, (ii+1)/(Double_t)ncol-tbw+shift+shiftw, (1+j)/(Double_t)nrow-tbh, 21);
	  fhPads[padi]->SetFillColor(kCyan-8);
	  fhPads[padi]->SetMargin(0.04,0.04,0.04,0.04);
	  fhPads[padi]->Draw();
	  padi++;
	}
      }
    }
  }else{
    float bw = 0.02, bh = 0.01, shift = 0,shiftw=-0.02;
    float tbw = bw, tbh = bh;
    Int_t padi = 0;
    if(!fhPads[0]){
      for(int ii=0; ii<ncol; ii++){
	for(int j=0; j<nrow; j++){
	  if(j==1) shift = 0.04;
	  else shift = 0;
	  fhPads[padi] =  new TPad(Form("P%d",ii*10+j),"T", ii/(Double_t)ncol+tbw+shift+shiftw , j/(Double_t)nrow+tbh, (ii+1)/(Double_t)ncol-tbw+shift+shiftw, (1+j)/(Double_t)nrow-tbh, 21);
	  fhPads[padi]->SetFillColor(kCyan-8);
	  fhPads[padi]->SetMargin(0.04,0.04,0.04,0.04);
	  fhPads[padi]->Draw(); 
	  padi++;
	}
      }
    }

  }
  
  Int_t np,tmax;
  Double_t max=0;
  if(maxz==0){
    for(Int_t p=0; p<nrow*ncol;p++){
      tmax = fhDigi[p]->GetBinContent(fhDigi[p]->GetMaximumBin());
      if(max<tmax) max = tmax;
    }
  }else{
    max = maxz;
  }

  if(maxz==-2 && minz==-2){ // optimize range
    for(Int_t p=0; p<nrow*ncol;p++){
      tmax = fhDigi[p]->GetMaximum();
      if(max<tmax) max = tmax;
    }
    if(max < 100) max = 100;
    Int_t tbins = 2000;
    TH1F *h = new TH1F("","",tbins,0,max);
    for(Int_t p=0; p<nrow*ncol;p++){
      for(Int_t i=0; i<64; i++){
	Double_t val = fhDigi[p]->GetBinContent(i);
	if(val!=0) h->Fill(val);
      }
    }
    Double_t integral;
    for(Int_t i=0; i<tbins; i++){
      integral = h->Integral(0,i);
      if(integral>5) {
	minz = h->GetBinCenter(i);
	break;
      } 
    }

    for(Int_t i=tbins; i>0; i--){
      integral = h->Integral(i,tbins);
      if(integral>5) {
	max = h->GetBinCenter(i);
	break;
      } 
    }
  }
 
  for(Int_t p=0; p<nrow*ncol;p++){
    if(layoutId == 1 || layoutId == 4)  np =p%3*5 + p/3;
    else np = p;
    
    fhPads[p]->cd();
    fhDigi[np]->Draw("col");
    if(maxz==-1)  max = fhDigi[np]->GetBinContent(fhDigi[np]->GetMaximumBin());
    fhDigi[np]->SetMaximum(max);
    fhDigi[np]->SetMinimum(minz);
    for(Int_t i=1; i<=8; i++){
      for(Int_t j=1; j<=8; j++){
  	Double_t weight = (double)(fhDigi[np]->GetBinContent(i,j))/(double)max *255;
  	if(weight > 0) digidata += Form("%d,%d,%d\n", np, (j-1)*8+i-1, (Int_t)weight);
      }
    }
  }
  cDigi->Modified();
  cDigi->Update();
  return digidata;
}

void initDigi(Int_t type=0){
  if(type == 0){
    for(Int_t m=0; m<3*5;m++){	
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
    }
  }
}

void resetDigi(){
    for(Int_t m=0; m<3*5;m++){	
      fhDigi[m]->Reset("M");
    }
}

void axisHits800x500(TH2 * hist){
  hist->SetStats(0);
  hist->SetTitle(Form("%d hits",(Int_t)hist->GetEntries()));
  hist->GetXaxis()->SetTitle("z, [cm]");
  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetXaxis()->SetTitleOffset(0.8);
  hist->GetYaxis()->SetTitle("y, [cm]");
  hist->GetYaxis()->SetTitleSize(0.05);
  hist->GetYaxis()->SetTitleOffset(0.7);
}

void axisAngle800x500(TH2 * hist){
  hist->SetStats(0);
  hist->SetTitle(Form("%d hits",(Int_t)hist->GetEntries()));
  hist->GetXaxis()->SetTitle("#theta, [degree]");
  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetXaxis()->SetTitleOffset(0.8);
  hist->GetYaxis()->SetTitle("photons per track, [#]");
  hist->GetYaxis()->SetTitleSize(0.05);
  hist->GetYaxis()->SetTitleOffset(0.7);
}

void axisAngle800x500(TH1 * hist){
  hist->SetStats(0);
  hist->SetTitle(Form("%d hits",(Int_t)hist->GetEntries()));
  hist->GetXaxis()->SetTitle("#theta, [degree]");
  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetXaxis()->SetTitleOffset(0.8);
  hist->GetYaxis()->SetTitle("photons per track, [#]");
  hist->GetYaxis()->SetTitleSize(0.05);
  hist->GetYaxis()->SetTitleOffset(0.7);
}

void axisTime800x500(TH2 * hist){
  hist->GetXaxis()->SetTitle("time, [ns]");
  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetXaxis()->SetTitleOffset(0.8);
  hist->GetYaxis()->SetTitle("entries, #");
  hist->GetYaxis()->SetTitleSize(0.05);
  hist->GetYaxis()->SetTitleOffset(0.7);
  hist->SetLineColor(1);
}

void axisTime800x500(TH1 * hist, TString xtitle = "time, [ns]"){
  TGaxis::SetMaxDigits(3);
  hist->GetXaxis()->SetTitle(xtitle);
  hist->GetXaxis()->SetTitleSize(0.05);
  hist->GetXaxis()->SetTitleOffset(0.8);
  hist->GetYaxis()->SetTitle("entries, [#]");
  hist->GetYaxis()->SetTitleSize(0.05);
  hist->GetYaxis()->SetTitleOffset(0.7);
  hist->SetLineColor(1);
}

void SetPrettyStyle(){
  // Canvas printing details: white bg, no borders.
  gStyle->SetCanvasColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasBorderSize(0);

  // Canvas frame printing details: white bg, no borders.
  gStyle->SetFrameFillColor(0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameBorderSize(0);

  // Plot title details: centered, no bg, no border, nice font.
  gStyle->SetTitleX(0.1);
  gStyle->SetTitleW(0.8);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleFillColor(0);

  // Font details for titles and labels.
  gStyle->SetTitleFont(42, "xyz");
  gStyle->SetTitleFont(42, "pad");
  gStyle->SetLabelFont(42, "xyz");
  gStyle->SetLabelFont(42, "pad");

  // Details for stat box.
  gStyle->SetStatColor(0);
  gStyle->SetStatFont(42);
  gStyle->SetStatBorderSize(1);
  gStyle->SetStatX(0.975);
  gStyle->SetStatY(0.9);

  // gStyle->SetOptStat(0);
}

void SetRootPalette(Int_t pal = 0){

 // pal =  1: rainbow\n"
 // pal =  2: reverse-rainbow\n"
 // pal =  3: amber\n"
 // pal =  4: reverse-amber\n"
 // pal =  5: blue/white\n"
 // pal =  6: white/blue\n"
 // pal =  7: red temperature\n"
 // pal =  8: reverse-red temperature\n"
 // pal =  9: green/white\n"
 // pal = 10: white/green\n"
 // pal = 11: orange/blue\n"
 // pal = 12: blue/orange\n"
 // pal = 13: white/black\n"
 // pal = 14: black/white\n"

  const Int_t NRGBs = 5;
  const Int_t NCont = 255;
  gStyle->SetNumberContours(NCont);

  if (pal < 1 && pal> 14) return;
  else pal--;

  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t red[14][NRGBs]   = {{ 0.00, 0.00, 0.87, 1.00, 0.51 },
			       { 0.51, 1.00, 0.87, 0.00, 0.00 },
			       { 0.17, 0.39, 0.62, 0.79, 1.00 },
			       { 1.00, 0.79, 0.62, 0.39, 0.17 },
			       { 0.00, 0.00, 0.00, 0.38, 1.00 },
			       { 1.00, 0.38, 0.00, 0.00, 0.00 },
			       { 0.00, 0.50, 0.89, 0.95, 1.00 },
			       { 1.00, 0.95, 0.89, 0.50, 0.00 },
			       { 0.00, 0.00, 0.38, 0.75, 1.00 },
			       { 0.00, 0.34, 0.61, 0.84, 1.00 },
			       { 0.75, 1.00, 0.24, 0.00, 0.00 },
			       { 0.00, 0.00, 0.24, 1.00, 0.75 },
			       { 0.00, 0.34, 0.61, 0.84, 1.00 },
			       { 1.00, 0.84, 0.61, 0.34, 0.00 }
  };
  Double_t green[14][NRGBs] = {{ 0.00, 0.81, 1.00, 0.20, 0.00 },		    
			       { 0.00, 0.20, 1.00, 0.81, 0.00 },
			       { 0.01, 0.02, 0.39, 0.68, 1.00 },
			       { 1.00, 0.68, 0.39, 0.02, 0.01 },
			       { 0.00, 0.00, 0.38, 0.76, 1.00 },
			       { 1.00, 0.76, 0.38, 0.00, 0.00 },
			       { 0.00, 0.00, 0.27, 0.71, 1.00 },
			       { 1.00, 0.71, 0.27, 0.00, 0.00 },
			       { 0.00, 0.35, 0.62, 0.85, 1.00 },
			       { 1.00, 0.75, 0.38, 0.00, 0.00 },
			       { 0.24, 1.00, 0.75, 0.18, 0.00 },
			       { 0.00, 0.18, 0.75, 1.00, 0.24 },
			       { 0.00, 0.34, 0.61, 0.84, 1.00 },
			       { 1.00, 0.84, 0.61, 0.34, 0.00 }
  };
  Double_t blue[14][NRGBs]  = {{ 0.51, 1.00, 0.12, 0.00, 0.00 },
			       { 0.00, 0.00, 0.12, 1.00, 0.51 },
			       { 0.00, 0.09, 0.18, 0.09, 0.00 },
			       { 0.00, 0.09, 0.18, 0.09, 0.00 },
			       { 0.00, 0.47, 0.83, 1.00, 1.00 },
			       { 1.00, 1.00, 0.83, 0.47, 0.00 },
			       { 0.00, 0.00, 0.00, 0.40, 1.00 },
			       { 1.00, 0.40, 0.00, 0.00, 0.00 },
			       { 0.00, 0.00, 0.00, 0.47, 1.00 },
			       { 1.00, 0.47, 0.00, 0.00, 0.00 },
			       { 0.00, 0.62, 1.00, 0.68, 0.12 },
			       { 0.12, 0.68, 1.00, 0.62, 0.00 },
			       { 0.00, 0.34, 0.61, 0.84, 1.00 },
			       { 1.00, 0.84, 0.61, 0.34, 0.00 }
  };


  TColor::CreateGradientColorTable(NRGBs, stops, red[pal], green[pal], blue[pal], NCont);
 
}

#ifdef prt__sim
void PrtInit(TString inFile="../build/hits.root", Int_t bdigi=0){

  SetRootPalette(1);
  gSystem->Load("../build/libprtdirclib.so");
  delete fCh;

  fCh = new TChain("data");

  fCh->Add(inFile);
  fCh->SetBranchAddress("PrtEvent", &fEvent);
  fNEntries = fCh->GetEntries();
  std::cout<<"Entries in chain:  "<<fNEntries <<std::endl;
  if(bdigi == 1) initDigi();
}

void PrtNextEvent(Int_t ievent, Int_t printstep){
  fCh->GetEntry(ievent);
  if(ievent%printstep==0 && ievent!=0) cout<<"Event # "<<ievent<< " # hits "<<fEvent->GetHitSize()<<endl;
  if(ievent == 0){
    if(gROOT->GetApplication()){
      TIter next(gROOT->GetApplication()->InputFiles());
      TObjString *os=0;
      while(os = (TObjString*)next()){
	fInfo += os->GetString()+" ";
      }
      fInfo += "\n";
    }
    fInfo += fEvent->PrintInfo();
    fMomentum = fEvent->GetMomentum().Mag() +0.01;
    fAngle = fEvent->GetAngle() + 0.01;
    fParticle =  fEvent->GetParticle();
    fTest1 = fEvent->GetTest1();
    fTest2 = fEvent->GetTest2();
  }
}
#endif

#ifdef prt__beam
void PrtInit(TString inFile="../build/hits.root", Int_t bdigi=0){

  SetRootPalette(1);
  delete fCh;

  fCh = new TChain("data");

  fCh->Add(inFile);
  fCh->SetBranchAddress("PrtEvent", &fEvent);
  fNEntries = fCh->GetEntries();
  std::cout<<"Entries in chain:  "<<fNEntries <<std::endl;
  if(bdigi == 1) initDigi();
}

void PrtNextEvent(Int_t ievent, Int_t printstep){
  fCh->GetEntry(ievent);
  if(ievent%printstep==0 && ievent!=0) cout<<"Event # "<<ievent<< " # hits "<<fEvent->GetHitSize()<<endl;
  if(ievent == 0){
    if(gROOT->GetApplication()){
      fInfo += "beam test";
      TIter next(gROOT->GetApplication()->InputFiles());
      TObjString *os=0;
      while(os = (TObjString*)next()){
	fInfo += os->GetString()+" ";
      }
      fInfo += "\n";
    }
    fInfo += fEvent->PrintInfo();
    fMomentum = fEvent->GetMomentum().Mag() +0.01;
    fAngle = fEvent->GetAngle() + 0.01;
    fParticle =  fEvent->GetParticle();
    fTest1 = fEvent->GetTest1();
    fTest2 = fEvent->GetTest2();
  }
}
#endif

TString randstr(Int_t len = 10){
  gSystem->Sleep(1500);
  srand (time(NULL));
  TString str = ""; 
  static const char alphanum[] =
    "0123456789"
    "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    "abcdefghijklmnopqrstuvwxyz";

  for (int i = 0; i < len; ++i) {
    str += alphanum[rand() % (sizeof(alphanum) - 1)];
  }
  return str;
}

Int_t getColorId(Int_t ind, Int_t style =0){
  Int_t cid = 1;
  if(style==0) {
    cid=ind+1;
    if(cid==5) cid =8;
    if(cid==3) cid =15;
  }
  if(style==1) cid=ind+300;
  return cid;
}

void waitPrimitive(TCanvas *c){
  c->Modified(); 
  c->Update(); 
  c->WaitPrimitive();
}

Int_t shiftHist(TH1F *hist, Double_t double_shift){
  Int_t bins=hist->GetXaxis()->GetNbins();
  Double_t xmin=hist->GetXaxis()->GetBinLowEdge(1);
  Double_t xmax=hist->GetXaxis()->GetBinUpEdge(bins);
  double_shift=double_shift*(bins/(xmax-xmin));
  Int_t shift=0;
  if(double_shift<0) shift=TMath::FloorNint(double_shift);
  if(double_shift>0) shift=TMath::CeilNint(double_shift);
  if(shift==0) return 0;
  if(shift>0){
    for(Int_t i=1; i<=bins; i++){
      if(i+shift<=bins) hist->SetBinContent(i,hist->GetBinContent(i+shift));
      if(i+shift>bins) hist->SetBinContent(i,0);
    }
    return 0;
  }
  if(shift<0){
    for(Int_t i=bins; i>0; i--){
      if(i+shift>0) hist->SetBinContent(i,hist->GetBinContent(i+shift));
      if(i+shift<=0) hist->SetBinContent(i,0);
    }    
    return 0;
  }
  return 1;
} 

void writeInfo(TString filename){
  ofstream myfile;
  myfile.open (filename);
  myfile << fInfo+"\n";
  myfile.close();
}

void writeString(TString filename, TString str){
  ofstream myfile;
  myfile.open (filename);
  myfile << str+"\n";
  myfile.close();
}

TString createDir(){
  TString finalpath = fSavePath;

  if(finalpath =="") return "";
  
  if(fSavePath == "auto") {
    TString dir = "data";
    gSystem->mkdir(dir);
    TDatime *time = new TDatime();
    TString path(""), stime = Form("%d.%d.%d", time->GetDay(),time->GetMonth(),time->GetYear()); 
    gSystem->mkdir(dir+"/"+stime);
    for(Int_t i=0; i<1000; i++){
      path = stime+"/"+Form("arid-%d",i);
      if(gSystem->mkdir(dir+"/"+path)==0) break;
    }
    gSystem->Unlink(dir+"/last");
    gSystem->Symlink(path, dir+"/last");
    finalpath = dir+"/"+path;
  }else{
    gSystem->mkdir(fSavePath,kTRUE);
  }
  writeInfo(finalpath+"/readme");
  return finalpath;
}

void save(TPad *c= NULL,TString path="", TString name="", Int_t what=0, Int_t style=0){
  if(c && path != "") {
    gROOT->SetBatch(1);
    Int_t w = 800, h = 400;
    if(style != -1){
      if(style == 1) {w = 800; h = 500;}
      if(style == 2) {w = 800; h = 600;}
      if(style == 3) {w = 800; h = 400;}
      if(style == 5) {w = 800; h = 900;} 
      if(style == 0){ 
	w = ((TCanvas*)c)->GetWindowWidth();
	h = ((TCanvas*)c)->GetWindowHeight();
      }
      	std::cout<<"w "<<w <<std::endl;
	std::cout<<"h "<<h <<std::endl;

      TCanvas *cc = new TCanvas(TString(c->GetName())+"exp","cExport",0,0,w,h);
      cc = (TCanvas*) c->DrawClone();
      cc->SetCanvasSize(w,h);
      if(style == 0) cc->SetBottomMargin(0.12);
      cc->Modified();
      cc->Update();
    
      cc->Print(path+"/"+name+".png");
      if(what==0) cc->Print(path+"/"+name+".pdf");
      if(what==0) cc->Print(path+"/"+name+".root");
    }else{
      c->Print(path+"/"+name+".png");
      if(what==0) c->Print(path+"/"+name+".pdf");
      if(what==0) c->Print(path+"/"+name+".root");
    }		    
    gROOT->SetBatch(0);
  }
}

TString createSubDir(TString dir="dir"){
  gSystem->mkdir(dir);
  return dir;
}

TList *gg_canvasList;
void canvasAdd(TString name="c",Int_t w=800, Int_t h=600){
  if(!gg_canvasList) gg_canvasList = new TList();
  TCanvas *c = new TCanvas(name,name,0,0,w,h); 
  gg_canvasList->Add(c);
}

void canvasAdd(TCanvas *c){
  if(!gg_canvasList) gg_canvasList = new TList();
  gg_canvasList->Add(c);
}

void canvasCd(TString name="c"){
  
}

void canvasDel(TString name="c"){
  TIter next(gg_canvasList);
  TCanvas *c=0;
  while((c = (TCanvas*) next())){
    if(c->GetName()==name) gg_canvasList->Remove(c);
  }
}

// style = 0 - for web blog
// style = 1 - for talk 
// what = 0 - save in png, pdf, root formats
// what = 1 - save in png format
void canvasSave(Int_t what=0, Int_t style=0){
  TIter next(gg_canvasList);
  TCanvas *c=0;
  TString path = createDir();
  while((c = (TCanvas*) next())){
    save(c, path, c->GetName(), what,style);
  }
}  

void normalize(TH1F* hists[],Int_t size){
  Double_t max = 0;
  Double_t min = 0;
  for(Int_t i=0; i<size; i++){
    Double_t tmax =  hists[i]->GetBinContent(hists[i]->GetMaximumBin());
    Double_t tmin = hists[i]->GetMinimum();
    if(tmax>max) max = tmax;
    if(tmin<min) min = tmin;
  }
  max += 0.05*max;
  for(Int_t i=0; i<size; i++){
    hists[i]->GetYaxis()->SetRangeUser(min,max);
  }
}
