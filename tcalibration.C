// tcalibration - routine for the prtdirc data calibration 
// original author: Roman Dzhygadlo - GSI Darmstad

#define TTSelector_cxx
#include "prttools.C"
#include "tcalibration.h"
 
TString ginFile(""), goutFile(""), gcFile("");
Int_t gSetup=2015, gTrigger(0), gMode(0), gComboId(0),  gMaxIn[maxch];
Double_t tdcRefTime[maxtdc],gTotO[maxch], gTotP[960][10],gLeOffArr[960],gEvtOffset(0);
TGraph *gGrIn[maxch], *gLeO[maxch], *gGrDiff[maxch];

Double_t walktheta(-16.5*TMath::Pi()/180.);
Double_t deltatheta(2.5*TMath::Pi()/180.);
Double_t walky0(42.23), walky1(44.00);
Double_t walk1x0(174.85), walk2x0(175.64);

Bool_t insideOfEllipce(Double_t x, Double_t y, Double_t x0, Double_t y0, Double_t w){
  Double_t r1(0.3), r2(0.9);

  Double_t xx = cos(w)*(x-x0)+sin(w)*(y-y0);
  Double_t yy = sin(w)*(x-x0)-cos(w)*(y-y0);

  return xx*xx/(r1*r1)+yy*yy/(r2*r2)<=1;
}

Double_t getTotWalk(Double_t tot,Int_t ch, Int_t type=0){ 
  Double_t minp(0), walk(0), d(0), min(100);

  if(type==0){
    if(ch<960){
      for(Int_t i=0; i<9; i++){
	if(gTotP[ch][i]<0.00000001) continue;
	d = gTotP[ch][i]-tot;
	if(fabs(d)<fabs(min)){
	  minp = gTotP[ch][i];
	  min = d;
	}
      }
    }
    Double_t wcorr(10);
    // if(ch/48==1) wcorr=0;
    // if(ch/48==5) wcorr=0;
    // if(ch/48==9) wcorr=0;
      
    if(fabs(min)<0.8) walk=-min*tan(wcorr*TMath::Pi()/180.);
    if(tot<10) walk-=(10-tot)*tan(6*TMath::Pi()/180.);
  }

  if(type==1){ //walk of the 1345
    walk += (38.85-tot)*tan(25*TMath::Pi()/180.); 
  }
  
  return walk;
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

  if(gcFile!="0"){
    TFile f(gcFile);
    TIter nextkey(f.GetListOfKeys());
    TKey *key;
    
    while ((key = (TKey*)nextkey())) {
      TGraph *gr = (TGraph*)key->ReadObj();
      TString name = gr->GetName();
      long long  ch = name.Atoll();
      std::cout<<name<<"  ch  "<<ch <<std::endl;
      
      Double_t x,y;
      if(ch <10000){ // spline calibration
	gGrIn[ch]= new TGraph(*gr);
      }else if(ch == 10000){ // line calibration
	for(Int_t i=0; i<maxch; i++){
	  gr->GetPoint(i,x,y);
	  gMaxIn[i] = (Int_t)(y+0.01);
	  //std::cout<<"ch  "<<i<< "  FT max"<<  gMaxIn[i]<<std::endl;	  
	}
      }else if(ch == 10001){ // read tot offsets
	for(Int_t i=0; i<maxch; i++){
	  gr->GetPoint(i,gTotO[i],y);
	  //std::cout<<"ch  "<<i<< " TOT off "<<  gTotO[i]<<std::endl;
	}
      }else if(ch == 10002){ // read tot peaks
	for(Int_t i=0; i<960*10; i++){
	  gr->GetPoint(i,x,y);
	  gTotP[i/10][i%10] = y;
	  //std::cout<<"ch  "<<i/10<< " peak "<< i%10<< " = " <<y<<std::endl;
	}
      }else if(ch == 10003){ // read LE offsets 1
	for(Int_t i=0; i<960; i++){
	  gr->GetPoint(i,gLeOffArr[i],y);
	}
      }else if(ch >= 20000 && ch < 30000){ // read LE offsets 3
	gLeO[ch-20000] = new TGraph(*gr);
      }else if(ch > 100000){ // read event offsets
	if(ginFile.Contains(name)){
	  gr->GetPoint(0,x,gEvtOffset);
	}
      }
    }
    f.Close();
  }

  std::cout<<"Initialization successful"<<std::endl;
}  

Bool_t TTSelector::Process(Long64_t entry){
  //  if(entry >10000 ) return kTRUE;
  Int_t tdc,ch,tofpid(0);
  Double_t grTime0(0), grTime1(0),grTime2(0),coarseTime(0),offset(0),triggerLe(0),triggerTot(0);
  Double_t time[10000], timeLe(0),timeT[10000],timeTot(0);
  Bool_t btrig(false),bmcpout(false),btof1(false),btof2(false);
 
  TString current_file_name  = TTSelector::fChain->GetCurrentFile()->GetName();
  Bool_t trbdata = current_file_name.Contains("trb");
  TObjArray *sarr = current_file_name.Tokenize("_");

  if(entry%1000==0) std::cout<<"event # "<< entry <<std::endl;
  GetEntry(entry);
  
  fEvent = new PrtEvent();
  //fEvent->SetReferenceChannel(gTrigger);
  
  for(Int_t i=0; i<Hits_ && i<10000; i++){
    tdc = map_tdc[Hits_nTrbAddress[i]];
    ch = GetChannelNumber(tdc,Hits_nTdcChannel[i])-1;
    
    if(ch==1344) btrig = true;
    if(ch==960)  btof1 = true;
    if(ch==1104) btof2 = true;
    if(ch==1248) bmcpout = true;
    
    time[i] =  5*(Hits_nEpochCounter[i]*pow(2.0,11) + Hits_nCoarseTime[i]); //coarsetime
    if(gcFile!="0") {
      //spline calib
      time[i] -= gGrIn[AddRefChannels(ch+1,tdc)]->Eval(Hits_nFineTime[i]+1);
      //linear calib
      // Double_t max = (Double_t) gMaxIn[AddRefChannels(ch+1,tdc)]-2;
      // time[i] = coarseTime-5*(Hits_nFineTime[i]-31)/(max-31);
    } // else time[i] -= (Hits_nFineTime[i]-31)*0.0102;

    if(Hits_nSignalEdge[i]==1){
      if(ch==gTrigger) grTime1 = time[i];
      if(Hits_nTdcChannel[i]==0) { //ref channel
	tdcRefTime[tdc] = time[i];
	if(gTrigger/48==tdc) grTime0 = time[i];
      }
    }else{
      timeT[i]=time[i];
      grTime2=time[i];
    }
  }

  Double_t tof1(0),tof2(0),tot1(0),tot2(0),toftime(0),mass(0);
  if(gMode==5){
    if(!(btrig && btof1 && btof2 && bmcpout)){
      fEvent->Clear();
      delete fEvent;
      return kTRUE;
    }

    for(Int_t i=0; i<Hits_ && i<10000; i++){
      if(Hits_nTdcErrCode[i]!=0) continue;
      if(Hits_nTdcChannel[i]==0) continue; // ref channel
      if(Hits_nSignalEdge[i]==0) continue; // tailing edge
      
      tdc = map_tdc[Hits_nTrbAddress[i]];
      ch = GetChannelNumber(tdc,Hits_nTdcChannel[i])-1;
      if(ch==960){
    	tof1 = time[i]-tdcRefTime[tdc];
    	tot1 = timeT[i+1] - time[i];
      }
      if(ch==1104){
    	tof2 = time[i]-tdcRefTime[tdc];
    	tot2 = timeT[i+1] - time[i];
      }
    }
    
    if(tof1!=0 && tof2!=0 && gMode==5) {
      toftime = tof2-tof1;
      if(insideOfEllipce(toftime, tot1, walk1x0, walky0,walktheta) && insideOfEllipce(toftime, tot2, walk1x0, walky1,-walktheta)){
    	toftime += (tot1-walky0)*tan(walktheta+deltatheta);
    	toftime += (tot2-walky1)*tan(-walktheta-deltatheta);
    	tofpid=2212;
    	mass = 0.938272046;
      }
      toftime = tof2-tof1;
      if(insideOfEllipce(toftime, tot1, walk2x0, walky0,walktheta) && insideOfEllipce(toftime, tot2, walk2x0, walky1,-walktheta)){
    	toftime = (tof2-tof1);
    	toftime += (tot1-walky0)*tan(walktheta+deltatheta);
    	toftime += (tot2-walky1)*tan(-walktheta-deltatheta);
    	tofpid=212;
    	mass=0.13957018;
      }
    }
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
      if(!trbdata && badcannel(ch)) continue;

      if(gMode>0){
	timeLe = time[i]-tdcRefTime[tdc];
	if(gTrigger!=0 && ch<960) timeLe = timeLe - (grTime1-grTime0);
      }else {
	timeLe = time[i];
	if(gTrigger!=0 && ch<960) timeLe = timeLe - grTime1;
      }
      if(gTrigger!=0) {
	triggerLe = grTime1 - grTime0;
	triggerTot=grTime2-grTime1;
      }
	
      timeTot = timeT[i+1] - time[i];

      if(ch<960) {
	//if(timeTot<0 || timeLe<20 || timeLe>40) continue;
	timeTot += 30-gTotO[ch];
	timeLe += getTotWalk(timeTot,ch);
	//timeLe += getTotWalk(triggerTot,ch,1);
	//if(gLeO[ch]) timeLe -=  gLeO[ch]->Eval(tot)-30;
	timeLe -= gLeOffArr[ch];
	if(tofpid>0){
	  Double_t mom = 7;
	  timeLe += 24.109/((mom/sqrt(mass*mass+mom*mom)*299792458))*1E9;
	}
      }

      if(gMode!=5 || tofpid!=0){
	hit.SetTdc(tdc);
	hit.SetChannel(ch);
	hit.SetMcpId(map_mcp[ch]);
	hit.SetPixelId(map_pix[ch]+1);
	hit.SetLeadTime(timeLe-gEvtOffset);
	hit.SetTotTime(timeTot);
	fEvent->AddHit(hit);
	nrhits++;
      }
    }
  }
  
  if(nrhits!=0) {
    fEvent->SetParticle(tofpid);
    fEvent->SetTest1(toftime);
    fTree->Fill();
  }
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
  gcFile = (cFile!="")? cFile: "0"; // calibration
  gTrigger = trigger;
  gMode = mode;
  if(gMode == 5) gTrigger=960;
  
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
