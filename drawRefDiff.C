#define TTSelector_cxx
#include "prttools.C"
#include "drawRefDiff.h"


TCanvas * cTdc = new TCanvas("cTdc","cTdc",0,0,800,400);
TCanvas * cRef = new TCanvas("cRef","cRef",0,0,800,400);

TH2F *hRef = new TH2F("hRef","hRef;tdc [#]; tdc [#]",30,0,30,30,0,30);
TH1F *hTdc = new TH1F("hTdc","htdc;tdc [#]; entries [#]",50,0,50); 
TH1F *hDif[50][50];
Int_t gRange(20);


TLine *gLine1 = new TLine(0,0,0,1000);
TBox  *gBox = new TBox();
Int_t gState = 0;
Int_t gRef(0),gtRef(0);
Bool_t gDrawBox(1);
Bool_t gDrawLine(1);


void TTSelector::Begin(TTree *){
  SetRootPalette(1);
  fSavePath = "auto";
  gStyle->SetOptStat(1001111);
  CreateMap();
  
  
  for(Int_t i=0; i<50; i++){
    for(Int_t j=0; j<50; j++){
      hDif[i][j] = new TH1F(Form("hDif%d_%d",i,j),Form("hDif%d_%d;time [ns];entries [#]",i,j),5000,-gRange,gRange);
    }
  }
 
} 

Bool_t TTSelector::Process(Long64_t entry){
  if(entry%1000==0) std::cout<<"event # "<< entry <<std::endl;
  GetEntry(entry);
  Double_t times[50]={0};
  
  for(Int_t i=0; i<Hits_; i++){
    Int_t tdc = map_tdc[Hits_nTrbAddress[i]];
    Int_t lch = Hits_nTdcChannel[i];
    Int_t ch = 48*tdc+lch;
    if(lch==0){
      hTdc->Fill(tdc);
      times[tdc]=Hits_fTime[i];
    }
  }
  
  for(Int_t i=0; i<50; i++){
    for(Int_t j=0; j<50; j++){
      Double_t diff =times[i]-times[j];
      hDif[i][j]->Fill(diff);
      hRef->Fill(i,j,diff);
    }
  }
  
  return kTRUE;
}

void TTSelector::Terminate(){
  //  canvasAdd(Form("cChTime",0),800,5000);
  cRef->cd();
  hRef->Draw("colz");
  hDif[2][5]->Draw();
  // canvasAdd(Form("cTdc",0),800,5000);
  cTdc->cd();
  hTdc->Draw();
  //  canvasSave(1,0); //-1
} 

void drawRefDiff(TString inFile= "/data.local/may2015/ce15145205719.hld.root",Int_t range=1000, Int_t events=0){ //ce15144155257.hld.root
  gRange = range;
  
  TChain* ch = new TChain("T");
  ch->Add(inFile);
  
  Int_t entries = ch->GetEntries();
  std::cout<<"Entries in chain:  "<< entries<<std::endl;
  if(events==0) entries = events;
  //entries = 10000;
  TTSelector *selector = new TTSelector();
  ch->Process(selector,"",entries);

  cTdc->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)", 0, 0,
		"exec3event(Int_t,Int_t,Int_t,TObject*)");
}

void exec3event(Int_t event, Int_t gx, Int_t gy, TObject *selected){
  TCanvas *c = (TCanvas *) gTQSender;
  TPad *pad = (TPad *) c->GetSelectedPad();
  if (!pad) return;
  Float_t x = pad->AbsPixeltoX(gx);
  Int_t binx = hTdc->GetXaxis()->FindBin( pad->PadtoX(x))-1;

  if(binx<0 || binx >49) return;
    
  if(event==1 && gState==0){
    gState = 1;
    gDrawBox = false;
    gRef = binx;
    return;
  }

  if(event==1 && gState==1){
    gState = 2;
    gDrawLine = false;
    return;
  }

  if(event==1 && gState==2){
    gState = 0;
    gDrawLine = true;
    gDrawBox = true;
    return;
  }
  if(gState==2) return;

  if(gtRef==binx) return;
  else gtRef = binx;
  
  if(gDrawLine){
    cTdc->cd();
    gLine1->SetX1(binx+0.5);
    gLine1->SetX2(binx+0.5);
    gLine1->SetY1(cTdc->GetUymin());
    gLine1->SetY2(cTdc->GetUymax());
    gLine1->SetLineColor(kRed);
    gLine1->Draw();
    cTdc->Update(); 
  }
 
  if(gDrawBox){
    cTdc->cd();
    gBox->SetX1(binx);
    gBox->SetX2(binx+1);
    gBox->SetY1(cTdc->GetUymin());
    gBox->SetY2(cTdc->GetUymax());
    gBox->SetFillColor(kRed);
    gBox->SetFillStyle(3001);
    gBox->Draw();
    cTdc->Update(); 
  }
      
  cRef->cd();
  hDif[gRef][binx]->Draw();
  cRef->Update(); 
}
