#define TTSelector_cxx
#include "prttools.C"
#include "tloop.h"

TH1F * hEvt1 = new TH1F("hEvt1","hEvt;event [#];counts [#]",10000,0,10000);
TH1F * hEvt0 = new TH1F("hEvt0","hEvt;event [#];counts [#]",10000,0,10000);

void TTSelector::Begin(TTree *){
  CreateMap();
}

Bool_t TTSelector::Process(Long64_t entry){
  if(entry%1000==0) std::cout<<"event # "<< entry <<std::endl;
  GetEntry(entry);

  for(Int_t i=0; i<Hits_; i++){
    Int_t ch = 48*map_tdc[Hits_nTrbAddress[i]]+Hits_nTdcChannel[i];
    
    if(ch == 455){
      if(Hits_nSignalEdge[i]==1)  hEvt1->Fill(entry);
      if(Hits_nSignalEdge[i]==0)  hEvt0->Fill(entry);
    }
  }
  
  return kTRUE;
}
 
void TTSelector::Terminate(){
  hEvt1->Draw();
  hEvt0->SetLineColor(2);
  hEvt0->Draw("same");
} 

void tloop(TString inFile= "/data.local/may2015/pl15139141259.hld.root"){
  TChain* ch = new TChain("T");
  ch->Add(inFile);
  
  Int_t entries = ch->GetEntries();
  std::cout<<"Entries in chain:  "<< entries<<std::endl;
  
  TTSelector *selector = new TTSelector();
  ch->Process(selector,"",entries);
}
