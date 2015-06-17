// run as: root loadprtlibs.C testt.C

#define TTSelector_cxx
#include "prttools.C"
#include "tloop.h"

TH1F * hEvt1 = new TH1F("hEvt1","hEvt;event [#];counts [#]",1000,-3000,3000);
TH1F * hCh = new TH1F("hCh","hCh;event [#];counts [#]",2000,0,2000);

void TTSelector::Begin(TTree *){
  CreateMap();
}

Int_t g_1786 = 0;
Int_t g_1788 = 0;
Int_t g_1780 = 0;
Int_t g_1790 = 0;
Int_t g_g = 0;

Bool_t TTSelector::Process(Long64_t entry){
  if(entry%1000==0) std::cout<<"event # "<< entry <<std::endl;
  GetEntry(entry);
  Double_t trt=0;
  Bool_t ch1786 = false;
  Bool_t ch1788 = false;
  Bool_t ch1790 = false;
  Bool_t ch1780 = false;

  for(Int_t i=0; i<Hits_; i++){
    Int_t ch = 48*map_tdc[Hits_nTrbAddress[i]]+Hits_nTdcChannel[i]-1;
    hCh->Fill(ch);
    if(ch==1776) trt = Hits_fTime[i];
    if(ch==1786)  {
      ch1786 = true;
      g_1786++;
    }
    if(ch==1788){
      ch1788 = true;
      g_1788++;
    }
    if(ch==1780){
      ch1780 = true;
      g_1780++;
    }
    if(ch==1790){
      ch1790 = true;
      g_1790++;
    }
    
  }

  if(ch1790 && ch1780 && ch1788 && ch1786) g_g++;

  for(Int_t i=0; i<Hits_; i++){
    Int_t ch = 48*map_tdc[Hits_nTrbAddress[i]]+Hits_nTdcChannel[i]-1;
    if(ch == 1786 && ch1788 && ch1790 && ch1780){
      if(trt!=0) hEvt1->Fill(Hits_fTime[i]-trt);
    }
  }
  return kTRUE;
}
 
void TTSelector::Terminate(){
  std::cout<<"g_1786  "<< g_1786 <<"   g_1788  "<< g_1788<<std::endl;
  std::cout<<"g_1790  "<< g_1790<<std::endl;
  std::cout<<"g_1780  "<< g_1780<<std::endl;
  std::cout<<"g_g "<<g_g <<std::endl;
  
  
  hEvt1->Draw();
  // hCh->Draw();
} 

void testt(TString inFile= "/data.local/may2015/ce15144155257.hld.root"){
  TChain* ch = new TChain("T");
  ch->Add(inFile);
  
  Int_t entries = ch->GetEntries();
  entries = 100000;
  std::cout<<"Entries in chain:  "<< entries<<std::endl;
  
  TTSelector *selector = new TTSelector();
  ch->Process(selector,"",entries);
}
