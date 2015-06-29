// run as: root loadprtlibs.C tloop.C

#define TTSelector_cxx
#include "prttools.C"
#include "drawTot.h"

TH1F * hTimeL = new TH1F("hTime","hTime;LE [ns];counts [#]",500,80,150);
TH1F * hTimeT = new TH1F("hTime","hTime;TE [ns];counts [#]",500,80,150);
TH1F * hTot = new TH1F("hTot","hTot;TOT [ns];counts [#]",200,-60,-20);
TH1F * hCh = new TH1F("hCh","hCh;channel [#];counts [#]",1500,0,1500);

TH2F * hLeTot = new TH2F("hLeTot","hLeTot;LE [ns];TOT [ns]",200,88,94,200,-50,-40);
TH2F * hTeTot = new TH2F("hTeTot","hTeTot;TE [ns];TOT [ns]",200,131,140,200,-50,-40);
TH2F * hLeTe  = new TH2F("hLeTe","hLeTe;LE [ns];TE [ns]",200,88,94,200,131,140);

void TTSelector::Begin(TTree *){
  fSavePath = "auto";
  CreateMap();
}

Bool_t TTSelector::Process(Long64_t entry){
  if(entry%1000==0) std::cout<<"event # "<< entry <<std::endl;
  GetEntry(entry);
  
  Double_t reftimes[tdcnum]={0};
  Double_t trailtimes[maxch]={0};
  Double_t trigTime0(0),trigTime1(0);
  Int_t prech(0);
  
  for(Int_t i=0; i<Hits_; i++){
    Int_t tdc = map_tdc[Hits_nTrbAddress[i]];
    Int_t ch = GetChannelNumber(tdc,Hits_nTdcChannel[i]);
    if(badcannel(ch)) continue;
    
    if(Hits_nTdcChannel[i]==0) reftimes[tdc]=Hits_fTime[i];
    if(Hits_nSignalEdge[i]==0){
      trailtimes[prech]=Hits_fTime[i];
      continue;
    }
    prech = ch;
    
    if(ch ==1346) trigTime1 = Hits_fTime[i];
    if(ch ==1344) trigTime0 = Hits_fTime[i];
  }

  
  for(Int_t i=0; i<Hits_; i++){
    Int_t tdc = map_tdc[Hits_nTrbAddress[i]];
    Int_t ch = GetChannelNumber(tdc,Hits_nTdcChannel[i]);
    if(badcannel(ch)) continue;
    if(Hits_nSignalEdge[i]==0 || Hits_nTdcChannel[i]==0 ) continue;

    hCh->Fill(ch);	
    if(ch == 100){
      Double_t le = Hits_fTime[i]-reftimes[tdc]-(trigTime1-trigTime0);
      Double_t te = trailtimes[ch]-reftimes[tdc]-(trigTime1-trigTime0);
      //      if(le>90.5 && le<91.5)
	{
	  hTimeL->Fill(le);
	  hTimeT->Fill(te);
	  hTot->Fill(le-te);
	  
	  hLeTot->Fill(le,le-te);
	  hTeTot->Fill(te,le-te);
	  hLeTe->Fill(le,te);	  
	}
      //std::cout<<"LE   "<< le << " TE "<< te <<std::endl;
    }
  }
  
  return kTRUE;
}
 
void TTSelector::Terminate(){
  canvasAdd("LE",800,500);
  hTimeL->Draw();
  canvasAdd("TE",800,500);
  hTimeT->Draw();
  canvasAdd("TOT",800,500);
  hTot->Draw();
  
  canvasAdd("LeTOT",800,500);
  hLeTot->Draw("colz");
  canvasAdd("TeTOT",800,500);
  hTeTot->Draw("colz");
  canvasAdd("LeTe",800,500);
  hLeTe->Draw("colz");
  //hCh->Draw();
  canvasSave(1,0);
} 

void drawTot(TString inFile= "../data/pilas_15178162456.hld.root_calibrated.root", Int_t events = 100000){
  gStyle->SetOptStat(1001111);
  TChain* ch = new TChain("T");
  ch->Add(inFile);
  
  Int_t entries = ch->GetEntries();
  std::cout<<"Entries in chain:  "<< entries<<std::endl;
  if(events!=0) entries = events;
  
  TTSelector *selector = new TTSelector();
  ch->Process(selector,"",entries);
}
