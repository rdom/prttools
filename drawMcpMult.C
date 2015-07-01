// run as: root loadprtlibs.C tloop.C

#define TTSelector_cxx
#include "prttools.C"
#include "drawMcpMult.h"

TString m_fileid;
Int_t g_hv;
TGraph *mGr1,*mGr2;
TH1F *hm[15],*hp[60];


void TTSelector::Begin(TTree *){
  fSavePath = "data/hv_mult";

  for(Int_t i=0; i<15; i++) hm[i] = new TH1F(Form("h_m%d",i),Form("h_m%d",i),300,0,500);
  for(Int_t i=0; i<60; i++) hp[i] = new TH1F(Form("h_p%d",i),Form("h_p%d",i),300,0,500);
  
  mGr1 = new TGraph();
  mGr2 = new TGraph();
  CreateMap();
}

Bool_t TTSelector::Process(Long64_t entry){
  if(entry%1000==0) std::cout<<"event # "<< entry <<std::endl;
  GetEntry(entry);
  
  Int_t mmcp[15]={0};
  Int_t mpadiva[60]={0};
  for(Int_t i=0; i<15; i++) mmcp[i]=0;
  for(Int_t i=0; i<60; i++) mpadiva[i]=0;
	
  for(Int_t i=0; i<Hits_; i++){
    Int_t tdc = map_tdc[Hits_nTrbAddress[i]];
    Int_t ch = GetChannelNumber(tdc,Hits_nTdcChannel[i]);
    if(badcannel(ch) || ch >= 15*64) continue;
    mmcp[ch/64]++;
    mpadiva[ch/16]++;
  }
  for(Int_t i=0; i<15; i++)  hm[i]->Fill(mmcp[i]);
  for(Int_t i=0; i<15*4; i++)  hp[i]->Fill(mpadiva[i]);
  
  return kTRUE;
}
 
void TTSelector::Terminate(){
  for(Int_t i=0; i<15; i++){
    hm[i]->Fit("gaus");
    Double_t mean = hm[i]->GetFunction("gaus")->GetParameter(1);
    mGr1->SetPoint(i,i,mean);
  }
  for(Int_t i=0; i<15*4; i++){
    hp[i]->Fit("gaus");
    Double_t mean = hp[i]->GetFunction("gaus")->GetParameter(1);
    mGr2->SetPoint(i,i,mean);
  }
  
  mGr1->SetTitle(Form("%d",g_hv));
  mGr1->SetName(Form("mcp_%d",g_hv));
  mGr2->SetTitle(Form("%d",g_hv));
  mGr2->SetName(Form("pad_%d",g_hv));
  TFile f(Form("mult_mcp_%d.root",g_hv),"recreate"); 
  mGr1->Write();
  mGr2->Write();
} 


void drawMcpMult(TString inFile= "../data/pilas_15178162456.hld.root_calibrated.root",Int_t hv =1000 ,Int_t events = 10000){
  g_hv = hv;
  m_fileid =  inFile;
  m_fileid =  m_fileid.Remove(0,m_fileid.Last('/')+1);
  gStyle->SetOptStat(1001111);
  TChain* ch = new TChain("T");
  ch->Add(inFile);
  
  Int_t entries = ch->GetEntries();
  std::cout<<"Entries in chain:  "<< entries<<std::endl;
  if(events!=0) entries = events;
  
  TTSelector *selector = new TTSelector();
  ch->Process(selector,"",entries);
}
