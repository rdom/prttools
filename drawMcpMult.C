// run as: root loadprtlibs.C tloop.C

#define TTSelector_cxx
#include "prttools.C"
#include "drawMcpMult.h"

TH1F *hCh;
TString m_fileid;
TH1F *hm[15],*hp[60];
Double_t hv, mcp, mult;
Int_t detid;

TFile *f;
TTree *tree;

void TTSelector::Begin(TTree *){
  fSavePath = "data/hv_mult";
  hCh = new TH1F("hCh","hCh",1500,0,1500);
  for(Int_t i=0; i<15; i++) hm[i] = new TH1F(Form("h_m%d",i),Form("h_m%d",i),100,0,100);
  for(Int_t i=0; i<60; i++) hp[i] = new TH1F(Form("h_p%d",i),Form("h_p%d",i),100,0,100);

  f = new TFile(Form("rhv_%d.root",(Int_t)hv),"recreate");
  tree= new TTree("hv","hc scan");
  tree->Branch("hv", &hv,"hv/D");
  tree->Branch("mcp", &mcp,"mcp/D");
  tree->Branch("mult",&mult,"mult/D");

  CreateMap();
}

Double_t tdcRefTime[100];
Bool_t TTSelector::Process(Long64_t entry){
  if(entry%1000==0) std::cout<<"event # "<< entry <<std::endl;
  GetEntry(entry);
  
  Int_t mmcp[15],mpadiva[60];
  for(Int_t i=0; i<15; i++) mmcp[i]=0;
  for(Int_t i=0; i<60; i++) mpadiva[i]=0;

  for(Int_t i=0; i<Hits_; i++){
    Int_t tdc = map_tdc[Hits_nTrbAddress[i]];
    Int_t ch = GetChannelNumber(tdc,Hits_nTdcChannel[i]);
    if(badcannel(ch)) continue;

    Double_t coarseTime = 5*(Hits_nEpochCounter[i]*pow(2.0,11) + Hits_nCoarseTime[i]);
    if(Hits_nTdcChannel[i]==0)  tdcRefTime[tdc] = coarseTime-(Hits_nFineTime[i]-31)*0.0102; 
  }
	
  for(Int_t i=0; i<Hits_; i++){
    Int_t tdc = map_tdc[Hits_nTrbAddress[i]];
    Int_t ch = GetChannelNumber(tdc,Hits_nTdcChannel[i]);
    if(Hits_nSignalEdge[i]==0 || Hits_nTdcChannel[i]==0) continue;
    if(badcannel(ch) || ch >= 15*64) continue;
    
    Double_t coarseTime = 5*(Hits_nEpochCounter[i]*pow(2.0,11) + Hits_nCoarseTime[i]);
    Double_t time  = coarseTime-(Hits_nFineTime[i]-31)*0.0102 -tdcRefTime[tdc]; 
    if(time < 0 || time > 50 ) continue;
    mmcp[ch/64]++;
    mpadiva[ch/16]++;
    hCh->Fill(ch);
  }

  if(detid==0) for(Int_t i=0; i<15; i++)  if(mmcp[i]>0) hm[i]->Fill(mmcp[i]);
  if(detid==1) for(Int_t i=0; i<15*4; i++) if(mpadiva[i]>0)  hp[i]->Fill(mpadiva[i]);
  
  return kTRUE;
}

void TTSelector::Terminate(){
  if(detid==0){
    for(Int_t i=0; i<15; i++){
      mcp = i;
      mult = hm[i]->GetMean();
      tree->Fill();
    }
  }

  if(detid==1){
    for(Int_t i=0; i<15*4; i++){
      mcp = i;
      mult = hp[i]->GetMean();
      tree->Fill();
    }
  }

  tree->Write();
}

void drawMcpMult(TString inFile= "../data/pilas_15178162456.hld.root",Int_t events = 10000, Int_t detectorid=0){
  detid = detectorid;
  TString str = inFile;

  TString shv = str.Remove(0,str.Last('/')+1).Remove(0,4).Remove(4,9);
  hv = shv.Atof();
  
  gStyle->SetOptStat(1001111);
  TChain* ch = new TChain("T");
  ch->Add(inFile);
  
  Int_t entries = ch->GetEntries();
  std::cout<<"Entries in chain:  "<< entries<<std::endl;
  if(events!=0) entries = events;
  
  TTSelector *selector = new TTSelector();
  ch->Process(selector,"",entries);
}
