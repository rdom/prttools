// tcalibration - routine for the prtdirc data calibration 
// original author: Roman Dzhygadlo - GSI Darmstad

#define TTSelector_cxx
#include "prttools.C"

#include "tcalibration.h"

const Int_t maxch =3000;
const Int_t nmcp(15), npix(64);

TString ginFile(""), goutFile(""), gcFile(""), gtFile("");
Int_t gSetup=2015, gTrigger(0), gMode(0), ctdc;

Int_t tdcnum(100);
const Int_t tdcmax(10000);

TString tdcsid2014[88] = 
  {"0010","0011","0012","0013","0110","0111","0112","0113","0210","0211","0212","0213","0310","0311","0312","0313","0410","0411","0412","0413"
   ,"0510","0511","0512","0513","0610","0611","0612","0613","0710","0711","0712","0713","0810","0811","0812","0813","0910","0911","0912","0913"
   ,"1010","1011","1012","1013","1110","1111","1112","1113","1210","1211","1212","1213","1310","1311","1312","1313","1410","1411","1412","1413"
   ,"1510","1511","1512","1513","1610","1611","1612","1613","1710","1711","1712","1713","1810","1811","1812","1813","1910","1911","1912","1913"
   ,"2010","2011","2012","2013","2110","2111","2112","2113"};


TString tdcsid2015[41] ={"2000","2001","2002","2003","2004","2005","2006","2007","2008","2009",
			 "200a","200b","200c","200d","200e","200f","2010","2011","2012","2013",
			 "2014","2015","2016","2018","2019","201a","201c","2020","2023","2024",
			 "2025","2026","2027","2028","2029","202a","202b","202c","202d","202e","202f"
};

TString tdcsid[100];
Int_t tdcid[100];
Double_t trbRefTime[100];

Double_t timeTe0[tdcmax][50];
Int_t mult[tdcmax];

Int_t tdcmap[tdcmax];
Int_t map_mpc[nmcp][npix];
Int_t map_mcp[maxch];
Int_t map_pix[maxch];
Int_t map_row[maxch];
Int_t map_col[maxch];

Int_t gComboId=0;
TGraph *gGr[maxch];
TGraph *gGrDiff[maxch];

void CreateMap(){
  Int_t seqid =0;
  for(Int_t i=0; i<tdcmax; i++){
    tdcmap[i]=-1;
    for(Int_t j=0; j<tdcnum; j++){
      if(i==TString::BaseConvert(tdcsid[j],16,10).Atoi()){
	tdcmap[i]=seqid++;
	break;
      }
    }
  }
  
  for(Int_t ch=0; ch<maxch; ch++){
    Int_t mcp = ch/128;
    Int_t pix = (ch - mcp*128)/2;
    Int_t col = pix/2 - 8*(pix/16);
    Int_t row = pix%2 + 2*(pix/16);
    pix = col*8+row;
    
    if(gSetup==2015){
      mcp = ch/64;
      pix = ch%64;	
      col = pix/2 - 8*(pix/16);
      row = pix%2 + 2*(pix/16);
      pix = col+8*row;
    }
    
    map_mpc[mcp][pix]=ch;
    map_mcp[ch] = mcp;
    map_pix[ch] = pix;
    map_row[ch] = row;
    map_col[ch] = col;
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
}

Bool_t badcannel(Int_t ch){

  // bad pixels july14
  if(ch ==  379) return true;
  if(ch ==  381) return true;
  if(ch == 1397) return true;
  if(ch == 1869) return true;

  if(ch == 1405) return true;
  if(ch == 1403) return true;
  if(ch == 1385) return true;
  if(ch == 1381) return true;
  if(ch == 1383) return true;
  if(ch == 1387) return true;

  return false;
}

Bool_t TTSelector::Process(Long64_t entry){
  Int_t trbSeqId,ch;
  Double_t timeTot(0), grTime0(0), grTime1(0),timeLe(0),coarseTime;
  Double_t time[50000];

  if(entry%1000==0) std::cout<<"event # "<< entry <<std::endl;
  GetEntry(entry);
  
  fEvent = new PrtEvent();
  // fEvent->SetReferenceChannel(gTrigger);
  for(Int_t i=0; i<Hits_; i++){
    trbSeqId = tdcmap[Hits_nTrbAddress[i]];
    ch = ctdc*trbSeqId+Hits_nTdcChannel[i];
    if(badcannel(ch)) continue; 

    if(++mult[ch]>50) continue;
    coarseTime = 5*(Hits_nEpochCounter[i]*pow(2.0,11) + Hits_nCoarseTime[i]);
    if(gcFile!="") time[i] = coarseTime-gGr[ch]->Eval(Hits_nFineTime[i]);
    else time[i] = coarseTime-(Hits_nFineTime[i]-31)*0.0102;
      
    timeTe0[ch][mult[ch]]=time[i];
    if(Hits_nTdcChannel[i]==0) {  // is ref channel
      trbRefTime[trbSeqId] = time[i];
      if(gTrigger!=0 && (gTrigger-ch)<=ctdc && (gTrigger-ch)>0) grTime0 = time[i];
    }
    if(gTrigger!=0 && ch==gTrigger) grTime1 = time[i];
  }
  PrtHit hit;
  Int_t nrhits=0;
  if((grTime0>0 && grTime1>0) || gTrigger==0){
    for(Int_t i=0; i<Hits_; i++){
      if(Hits_nTrbAddress[i]==0) continue;
      trbSeqId = tdcmap[Hits_nTrbAddress[i]];
      ch = ctdc*trbSeqId+Hits_nTdcChannel[i];
      if(badcannel(ch)) continue; 
      
      if(gSetup==2014 && ch%2==0 || Hits_nTdcChannel[i]==0) continue; // go away trailing edge
      if(gSetup==2015 && Hits_nTdcChannel[i]==0) continue; // go away ref channel
      
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

      hit.SetTdc(Hits_nTdcChannel[i]);
      hit.SetChannel(ch);
      hit.SetMcpId(map_mcp[ch]);
      hit.SetPixelId(map_pix[ch]+1);
      hit.SetLeadTime(timeLe);
      hit.SetTotTime(timeTot);
      fEvent->AddHit(hit);
      nrhits++;
    }
  }

  for(Int_t i=0; i<Hits_; i++){
    trbSeqId = tdcmap[Hits_nTrbAddress[i]];
    ch = ctdc*trbSeqId+Hits_nTdcChannel[i];
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
}

void tcalibration(TString inFile= "../../data/cj.hld.root", TString outFile= "outFileC.root", TString cFile= "calib.root", TString tFile= "calibOffsets.root", Int_t trigger=0,  Int_t sEvent =0, Int_t eEvent=0, Int_t build = 0){ //1920
  if(build==1) return;
  ginFile = inFile;
  goutFile = outFile;
  gcFile = cFile; // fine time calibration
  gtFile = tFile; // pilas offsets + walk corrections
  gTrigger = trigger;
  if(gtFile=="") gMode = 0;
  else  gMode = 1;

  gSetup = 2015;
  if(gSetup==2014){
    ctdc = 32;
    tdcnum = 88;
    for(Int_t i=0; i<tdcnum; i++){
      tdcsid[i]=tdcsid2014[i];
    }
  }else{
    ctdc = 48;
    tdcnum = 41;
    for(Int_t i=0; i<tdcnum; i++){
      tdcsid[i]=tdcsid2015[i];
    }
  }
  
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
