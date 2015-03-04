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

#include <iostream>
#include <fstream>

#ifdef prt__sim
   class PrtEvent;
   class PrtHit;
   PrtEvent* fEvent = 0;
#endif 

#ifdef prt__beam
   class TPrtHit;
   class TPrtEvent;
   TPrtEvent* fEvent=0;
#endif  


TChain *fCh = 0;
Int_t fNEntries=0;
Int_t fSaveFlag = 1;
TString fInfo = "", fPath;
TH2F *fhDigi[15];
TPad* fhPads[15];
TPad * fhPglobal;
TCanvas* cDigi;

TString drawDigi(TString digidata="", Int_t layoutId = 0, Double_t maxz = 0, Double_t minz = 0){
  if(!cDigi) cDigi = new TCanvas("cDigi","cDigi",0,0,800,500);
  cDigi->cd();
  // TPad * pp = new TPad("P","T",0.06,0.135,0.93,0.865);
  if(!fhPglobal){
    fhPglobal = new TPad("P","T",0.01,0.01,0.99,0.99);
    fhPglobal->SetFillStyle(0);
    fhPglobal->Draw();
  }
  fhPglobal->cd();

  Int_t nrow = 3, ncol = 5;

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

  Int_t np,tmax;
  Double_t max=0;
  if(maxz==0){
    for(Int_t p=0; p<nrow*ncol;p++){
      tmax = fhDigi[p]->GetMaximum();
      if(max<tmax) max = tmax;
    }
  }else{
    max = maxz;
  }

  if(maxz==-2 && minz==-2){ // optimize range
    max = maxz;
    Int_t tbins = 2000;
    TH1F *h = new TH1F("","",tbins,0,200);
    for(Int_t p=0; p<nrow*ncol;p++){
      for(Int_t i=0; i<64; i++){
	Double_t val = fhDigi[p]->GetBinContent(i);
	if(val!=0) h->Fill(val);
      }
    }
   
    Double_t integral;
    for(Int_t i=0; i<tbins; i++){
      integral = h->Integral(0,i);
      if(integral>10) {
	minz = h->GetBinCenter(i);
	break;
      } 
    }

    for(Int_t i=tbins; i>0; i--){
      integral = h->Integral(i,tbins);
      if(integral>10) {
	max = h->GetBinCenter(i);
	break;
      } 
    }
  }

  for(Int_t p=0; p<nrow*ncol;p++){
    if(layoutId == 1)  np =p%3*5 + p/3;
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
  gSystem->Load("../build/libprtlib.so");
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
    fInfo += fEvent->PrintInfo();
  }
}
#endif  

#ifdef prt__beam
void PrtInit(TString inFile="../build/hits.root", Int_t bdigi=0){
  SetRootPalette(1);
  delete fCh;
  fCh = new TChain("M");
  fCh->Add(inFile);
  fCh->SetBranchAddress("TPrtEvent", &fEvent);
  fNEntries = fCh->GetEntries();
  std::cout<<"Entries in chain:  "<< fNEntries<<std::endl;
  if(bdigi == 1) initDigi();
}

void PrtNextEvent(Int_t ievent, Int_t printstep){
  fCh->GetEntry(ievent);
  if(ievent%(printstep/20)==0 && ievent!=0) std::cerr << ".";
  if(ievent%printstep==0 && ievent!=0) std::cout<<" Event # "<<ievent << std::endl;
  if(ievent == 0){
    fInfo += "beam test";
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
  if(style==0) cid=ind+1;
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

TString gg_path = "";

void writeInfo(TString filename, TString info,Int_t flag=1){
  if(flag<1) return;
  ofstream myfile;
  myfile.open (gg_path+"/"+filename);
  myfile << info+"\n";
  myfile.close();
}

// flag < 0 - do not save anything
// flag = 0 - save in a new folder based on date
// flag = 2 - save in a folder dir (second argument)  
void save(TPad *c= NULL, TString dir="rdata", TString name="", TString info="", Int_t flag=1, Int_t bfiles=0){
  if(flag<1) return;
  TString path = "";
  if(flag == 0 && flag!=2){
    TDatime *time = new TDatime();
    TString stime = Form("%d.%d.%d", time->GetDay(),time->GetMonth(),time->GetYear()); 
    gSystem->mkdir(dir+"/"+stime);
    for(Int_t i=0; i<1000; i++){
      path = dir+"/"+stime+"/"+Form("arid-%d",i);
      if(gSystem->mkdir(path)==0) break;
    }
    writeInfo("readme", info);
  }else{
    path = dir;
  }
  
  if(c) {
    c->Modified();
    c->Update();
    c->Print(path+"/"+name+".png");
    if(bfiles==0) c->Print(path+"/"+name+".pdf");
    if(bfiles==0) c->Print(path+"/"+name+".root");
  }
}

TString createDir(TString dir="rdata", TString info = "", Int_t flag=1){
  if(flag<1) return "";
  TString path = "";
  if(flag==2) {
    gSystem->mkdir(dir,kTRUE);
    path = dir;
  }else{
    gSystem->mkdir(dir);
    TDatime *time = new TDatime();
    TString stime = Form("%d.%d.%d", time->GetDay(),time->GetMonth(),time->GetYear()); 
    gSystem->mkdir(dir+"/"+stime);
    for(Int_t i=0; i<1000; i++){
      path = dir+"/"+stime+"/"+Form("arid-%d",i);
      if(gSystem->mkdir(path)==0) break;
    }
  }
  gSystem->Unlink(dir+"/last");
  gSystem->Symlink(path, dir+"/last");
  gg_path = path;
  writeInfo("readme", info);
  return path;
}

TString createSubDir(TString dir="dir", Int_t flag=1){
  if(flag<1) return "";
  TString path = "";
  gSystem->mkdir(dir);
  return dir;
}

