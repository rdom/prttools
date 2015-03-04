// tcalibration - routine for the prtdirc data calibration 
// original author: Roman Dzhygadlo - GSI Darmstad

#define TTSelector_cxx
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
#include <TChainElement.h>
#include <TKey.h>

#include "tcalibration.h"

const Int_t maxfiles(200);
const Int_t maxch =3000;
const Int_t nmcp(15), npix(64);
TString fileList[maxfiles];

TH1F *hCh;

TString ginFile(""), goutFile(""), gcFile(""), gtFile("");
Int_t gTrigger(0), gMode(0);;

const Int_t tdcnum(88);
const Int_t tdcmax(10000);
TString trbsid[tdcnum] = 
  {"0010","0011","0012","0013","0110","0111","0112","0113","0210","0211","0212","0213","0310","0311","0312","0313","0410","0411","0412","0413"
   ,"0510","0511","0512","0513","0610","0611","0612","0613","0710","0711","0712","0713","0810","0811","0812","0813","0910","0911","0912","0913"
   ,"1010","1011","1012","1013","1110","1111","1112","1113","1210","1211","1212","1213","1310","1311","1312","1313","1410","1411","1412","1413"
   ,"1510","1511","1512","1513","1610","1611","1612","1613","1710","1711","1712","1713","1810","1811","1812","1813","1910","1911","1912","1913"
   ,"2010","2011","2012","2013","2110","2111","2112","2113"};

Int_t tdcid[tdcnum];
Double_t trbRefTime[tdcnum];

Double_t timeTe0[tdcmax][50];
Int_t mult[tdcmax];

Int_t tdcmap[tdcmax];
Int_t mcpmap[tdcmax];
Int_t pixmap[tdcmax];
Int_t chmap[nmcp][npix];

Int_t gComboId=0;
TGraph *gGr[maxch];
TGraph *gGrDiff[maxch];

//TH1F  *hL = new TH1F("hL", "hL" , 500,150,200);

void CreateMap(){
  Int_t seqid =0;
  for(Int_t i=0; i<tdcmax; i++){
    tdcmap[i]=-1;
    mcpmap[i]=-1;
    for(Int_t j=0; j<tdcnum; j++){
      if(i==TString::BaseConvert(trbsid[j],16,10).Atoi()){
	tdcmap[i]=seqid++;
	mcpmap[i]=j/4;
	pixmap[i]=j-j/32;
	break;
      }
    }
  }

  for(Int_t ch=0; ch<maxch; ch++){
    Int_t mcp = ch/128;
    Int_t pix = (ch - mcp*128)/2;
    Int_t col = 7-(pix/2 - 8*(pix/16));
    Int_t row = pix%2 + 2*(pix/16);
    pix = col*8+row;
    //std::cout<<"ch  "<<ch <<"  m "<<mcp <<"  p  "<<pix <<std::endl;
    
    chmap[mcp][pix]=ch;
  }
}

void TTSelector::Begin(TTree *){
  TString option = GetOption();
  TObjArray *strobj = option.Tokenize(" ");
  gTrigger = ((TObjString*)strobj->At(0))->GetString().Atoi();
  
  gMode = ((TObjString*)strobj->At(1))->GetString().Atoi();
  CreateMap();
  TString filedir=ginFile;
  filedir.Remove(filedir.Last('.')-4);
  fFile = new TFile(goutFile,"RECREATE");
  fTree = new TTree("M","Tree for GSI Prt Analysis");  
  fEvent = new PrtEvent();
  fTree->Branch("PrtEvent", "PrtEvent", &fEvent, 64000, 2);

  TFile f(gcFile);
  TIter nextkey(f.GetListOfKeys());
  TKey *key;

  while ((key = (TKey*)nextkey())) {
    TGraph *gr = (TGraph*)key->ReadObj();
    TString name = gr->GetName();
    Int_t channel = name.Atoi();

    gGr[channel]= new TGraph(*gr);
  }
  f.Close();

  if(gMode == 1){
    TFile f2(gtFile);
    TIter nextkey2(f2.GetListOfKeys());
    TKey *key2;

    while ((key2 = (TKey*)nextkey2())) {
      TGraph *gr = (TGraph*)key2->ReadObj();
      TString name = gr->GetName();
      Int_t channel = name.Atoi();

      gGrDiff[channel]= new TGraph(*gr);
    }
    f2.Close();
  }
  std::cout<<"Initialization successful"<<std::endl;
  
}

Bool_t TTSelector::Process(Long64_t entry){
  Int_t trbSeqId,ch,mcp,pix,col,row;;
  Double_t timeTot(0), grTime0(0), grTime1(0),timeLe(0),coarseTime;
  Double_t time[50000];

  if(entry%1000==0) std::cout<<"event # "<< entry <<std::endl;
  GetEntry(entry);
  
  fEvent = new PrtEvent();
  // fEvent->SetReferenceChannel(gTrigger);
  for(Int_t i=0; i<Hits_; i++){
    trbSeqId = tdcmap[Hits_nTrbAddress[i]];
    ch = 32*trbSeqId+Hits_nTdcChannel[i];
    
    if(++mult[ch]>50) continue;
    coarseTime = 5*(Hits_nEpochCounter[i]*pow(2.0,11) + Hits_nCoarseTime[i]);
    time[i] = coarseTime-gGr[ch]->Eval(Hits_nFineTime[i]);

    timeTe0[ch][mult[ch]]=time[i];
    if(Hits_nTdcChannel[i]==0) {  // is ref channel
      trbRefTime[trbSeqId] = time[i];
      if(gTrigger!=0 && (gTrigger-ch)<=32 && (gTrigger-ch)>0) grTime0 = time[i];
    }
    if(gTrigger!=0 && ch==gTrigger) grTime1 = time[i];
  }
  PrtHit hit;
  Int_t nrhits=0;
  if((grTime0>0 && grTime1>0) || gTrigger==0){
    for(Int_t i=0; i<Hits_; i++){
      if(Hits_nTrbAddress[i]==0) continue;
      trbSeqId = tdcmap[Hits_nTrbAddress[i]];
      ch = 32*trbSeqId+Hits_nTdcChannel[i];
      mcp = ch/128;
      pix = (ch%128)/2;	
      col = pix/2 - 8*(pix/16);
      row = pix%2 + 2*(pix/16);
      pix = (7-col)*8+row;
      
      if(ch%2==0) continue; // go away trailing edge
      if(ch>3000) continue;

      timeLe = time[i]-trbRefTime[trbSeqId];
      timeLe = timeLe - (grTime1-grTime0);
      timeTot = timeTe0[ch+1][1] - timeTe0[ch][1];
      
      if(gMode == 1){
	if(ch>1920) continue;
	else timeLe -= gGrDiff[ch]->Eval(timeTot);
	if(abs(timeTot)>10) continue; 
	if(abs(timeLe)>300) continue; 
      }

      //PrtHit hit(Hits_nTrbAddress[i],Hits_nTdcChannel[i],ch,mcp,pix+1,timeLe,timeTot);
      hit.SetTdc(Hits_nTdcChannel[i]);
      hit.SetChannel(ch);
      hit.SetMcpId(mcp);
      hit.SetPixelId(pix+1);
      hit.SetLeadTime(timeLe);
      hit.SetTotTime(timeTot);
      fEvent->AddHit(hit);
      nrhits++;
    }
  }

  for(Int_t i=0; i<Hits_; i++){
    trbSeqId = tdcmap[Hits_nTrbAddress[i]];
    ch = 32*trbSeqId+Hits_nTdcChannel[i];
    mult[ch]=0;
    for(Int_t j=0; j<50; j++){
      timeTe0[ch][j]=0; 
    }
  }
  if(nrhits!=0) fTree->Fill();
  fEvent->Clear();
  delete fEvent;

  return kTRUE;
}

void TTSelector::Terminate(){
  fFile->Write();
  fFile->Close();
  // hL->Draw();
  // hL->Fit("gaus","V","E1",175,185);
}

void tcalibration(TString inFile= "../../data/cj.hld.root", TString outFile= "outFileC.root", TString cFile= "calib.root", TString tFile= "calibOffsets.root", Int_t trigger=2560,  Int_t sEvent =0, Int_t eEvent=0, Int_t build = 0){ //1920
  if(build==1) return;
  ginFile = inFile;
  goutFile = outFile;
  gcFile = cFile; // fine time calibration
  gtFile = tFile; // pilas offsets + walk corrections
  gTrigger = trigger;
  if(gtFile=="") gMode = 0;
  else  gMode = 1;

  TChain* ch = new TChain("T");
  ch->Add(ginFile);
  
  Int_t entries = ch->GetEntries();
  TTSelector *selector = new TTSelector();
  TString option = Form("%d %d",gTrigger,gMode);
  
  if(eEvent==0){
    std::cout<<"Entries in chain:  "<< entries<<std::endl;
    ch->Process(selector,option,entries);
  }else{
    ch->Process(selector,option,eEvent-sEvent,sEvent); 
  }
}
