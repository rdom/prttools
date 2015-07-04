// tcalibration - routine for the prtdirc data calibration 
// original author: Roman Dzhygadlo - GSI Darmstad

#define TTSelector_cxx
#include "prttools.C"
#include "tcalibration.h"
 
TString ginFile(""), goutFile(""), gcFile(""), gtFile("");
Int_t gSetup=2015, gTrigger(0), gMode(0);
Double_t tdcRefTime[100];
Int_t gComboId=0;
TGraph *gGr[maxch];
TGraph *gGrDiff[maxch];

void TTSelector::Begin(TTree *){
  TString option = GetOption();
  TObjArray *strobj = option.Tokenize(" ");
  gTrigger = ((TObjString*)strobj->At(0))->GetString().Atoi();
  gMode = ((TObjString*)strobj->At(1))->GetString().Atoi();
  CreateMap();
  TString filedir=ginFile;
  filedir.Remove(filedir.Last('.')-4);
  fFile = new TFile(goutFile,"RECREATE");
  fTree = new TTree("data","Tree for GSI Prt Analysis");  
  fEvent = new PrtEvent();
  fTree->Branch("PrtEvent", "PrtEvent", &fEvent, 64000, 2);
  
  if(gcFile!=""){
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
  }

  if(gtFile!=""){
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
  Int_t tdc,ch;
  Double_t timeTot(0), grTime0(0), grTime1(0),timeLe(0),coarseTime(0),offset(0);
  Double_t time[10000], timeT[10000];

  TString current_file_name  = TTSelector::fChain->GetCurrentFile()->GetName();
  TObjArray *sarr = current_file_name.Tokenize("_");

  if(entry%1000==0) std::cout<<"event # "<< entry <<std::endl;
  GetEntry(entry);
  
  fEvent = new PrtEvent();
  // fEvent->SetReferenceChannel(gTrigger);
  
  for(Int_t i=0; i<Hits_ && i<10000; i++){
    tdc = map_tdc[Hits_nTrbAddress[i]];
    ch = GetChannelNumber(tdc,Hits_nTdcChannel[i]);
    if(badcannel(ch)) continue;

    coarseTime = 5*(Hits_nEpochCounter[i]*pow(2.0,11) + Hits_nCoarseTime[i]);

    if(gcFile!="") time[i] = coarseTime-gGr[AddRefChannels(ch)]->Eval(Hits_nFineTime[i]);
    else time[i] = Hits_fTime[i]; //coarseTime-(Hits_nFineTime[i]-31)*0.0102; //Hits_fTime[i];//
    
    if(Hits_nSignalEdge[i]==0){
      timeT[i]=time[i];
      continue;
    }
    
    if(Hits_nTdcChannel[i]==0) {  // is ref channel
      tdcRefTime[tdc] = time[i];
      if(gTrigger!=0 && (gTrigger-ch)<=ctdc && (gTrigger-ch)>0) grTime0 = time[i];
    }
    if(gTrigger!=0 && ch==gTrigger) grTime1 = time[i];
  }
  
  PrtHit hit;
  Int_t nrhits=0;
  if((grTime0>0 && grTime1>0) || gTrigger==0){
    for(Int_t i=0; i<Hits_ && i<10000; i++){
      if(Hits_nTdcErrCode[i]!=0) continue;
      if(Hits_nTdcChannel[i]==0) continue; // ref channel
      if(Hits_nSignalEdge[i]==0) continue; // tailing edge
      
      tdc = map_tdc[Hits_nTrbAddress[i]];
      ch = GetChannelNumber(tdc,Hits_nTdcChannel[i])-1;
      if(badcannel(ch)) continue;

      if(gMode==1)timeLe = time[i]-tdcRefTime[tdc];
      else timeLe = time[i];

      if(gTrigger!=0) timeLe = timeLe - (grTime1-grTime0);
      
      timeTot = timeT[i+1] - time[i];

      
      // if(timeLe>-300 && timeLe <-100) std::cout<<"timeLe "<<timeLe << " - "<<timeTot  <<std::endl;      
           
      if(gtFile!=""){
	if(ch<960) timeLe -= gGrDiff[ch]->Eval(timeTot);	
	//if(abs(timeTot)>10) continue; 
	//if(abs(timeLe)>300) continue; 
      }

      hit.SetTdc(tdc);
      hit.SetChannel(ch);
      hit.SetMcpId(map_mcp[ch]);
      hit.SetPixelId(map_pix[ch]+1);
      hit.SetLeadTime(timeLe);
      hit.SetTotTime(timeTot);
      fEvent->AddHit(hit);
      nrhits++;
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
}

void tcalibration(TString inFile= "../../data/cj.hld.root", TString outFile= "outFileC.root", TString cFile= "calib.root", TString tFile= "calibOffsets.root", Int_t trigger=0,  Int_t sEvent =0, Int_t eEvent=0, Int_t mode=1, Int_t build=0){
  if(build==1) return;
  ginFile = inFile;
  goutFile = outFile;
  gcFile = cFile; // fine time calibration
  gtFile = tFile; // pilas offsets + walk corrections
  gTrigger = trigger;
  gMode = mode;
  
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
