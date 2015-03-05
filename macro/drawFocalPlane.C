#define prt__sim

#include "../src/PrtEvent.h"
#include "../src/PrtHit.h"

#include "prttools.C"

double findVertex(TVector3 v1,TVector3 m1, TVector3 v2,TVector3 m2, TVector3* newvertex){
  TVector3 head0 = v1;
  TVector3 tail0 = v1+100*m1;

  TVector3 head1 = v2;
  TVector3 tail1 = v2+100*m2;

  TVector3 dir0 = head0 - tail0;
  TVector3 dir1 = head1 - tail1;

  TVector3 diff = tail0 - tail1;

  float a = dir0.Dot(dir0);//always >= 0
  float b = dir0.Dot(dir1);
  float c = dir1.Dot(dir1);//always >= 0
  float d = dir0.Dot(diff);
  float e = dir1.Dot(diff);
  float f = a * c - b * b;//always >= 0

  float sc;
  float tc;

  if(f<0.000001){//The lines are almost parallel
    sc = 0.0f;
    tc = b > c ? d / b : e / c;//Use the largest denominator
  }else{
    sc = (b * e - c * d) / f;
    tc = (a * e - b * d) / f;
  }
  
  TVector3 point0 = tail0 + dir0 * sc;
  TVector3 point1 = tail1 + dir1 * tc;
  newvertex->SetXYZ((point0.X()+point1.X())/2.,(point0.Y()+point1.Y())/2.,(point0.Z()+point1.Z())/2.);
  return (point1-point0).Mag();
}


void drawFocalPlane(TString infile="../build/hits.root", Double_t r1 = 48.8, Double_t r2 = 29.1, Int_t it1=0, Int_t it2=0, Double_t energy=-1){
  gSystem->Load("../build/libprtdirclib.so");
  PrtInit(infile,0);
  gStyle->SetOptStat(0);

  Double_t radiatorL(1250); //bar
  // Double_t radiatorL(1224.9); //plate
  TVector3 res;
  TH2F *hFp1 = new TH2F("hFp1",Form("r_{1}=%2.2f    r_{2}=%2.2f;x, [cm];y, [cm]",r1,r2),500,0,50,500,-30,30 );
  TH2F *hFp2 = new TH2F("hFp2",Form("r_{1}=%2.2f    r_{2}=%2.2f;z, [cm];y, [cm]",r1,r2),500,-30,30,500,-30,50 );

  PrtHit hit[2];
  Int_t nevents = fCh->GetEntries();
  Int_t count(0);
  for (Int_t ievent=0; ievent<nevents; ievent++){
    PrtNextEvent(ievent,1000);
    Int_t nhits = fEvent->GetHitSize();
    if(nhits!=2) continue;
    for(Int_t h=0; h<nhits; h++) hit[h] = fEvent->GetHit(h);
    Double_t d = findVertex(hit[0].GetGlobalPos(),hit[0].GetMomentum().Unit(),hit[1].GetGlobalPos(),hit[1].GetMomentum().Unit(), &res);
    if(d<1){
      Double_t x = -(res.X()+radiatorL/2.)/10.;
      if(x>6){
	hFp1->Fill(x,res.Z()/10.);
	hFp2->Fill(res.Y()/10.,res.Z()/10.);
      }
      count++;
    }
  } 
  Double_t eff = 100*count/(Double_t)nevents;

  TString senergy = "";
    if(energy>-1){
      Double_t m = 1, cm = 0.01, nanometer = 0.000000001, GeV = 1000000000;
      Double_t LambdaE = 2.0 * 3.14159265358979323846 * 1.973269602e-16 * m * GeV;
      Double_t lam = LambdaE/(energy*nanometer);

      senergy = Form("   E=%2.2f [eV] / #lambda=%2.0f [nm]",energy,lam);
    }
  hFp1->SetTitle(Form("r_{1}=%2.2f    r_{2}=%2.2f   #varepsilon=%2.0f",r1,r2,eff)+senergy);
  canvasAdd(Form("fp_%d_%d",it1,it2),800,500);
  hFp1->Draw("colz");
  canvasAdd(Form("fp2_%d_%d",it1,it2),600,800);
  hFp2->Draw("colz");
  canvasSave(0,"drawFocalPlane.C",1,"fp");
  
}

