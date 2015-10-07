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
Int_t gMax[maxch];
Double_t gTotO[maxch], gTotP[960][10];;

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

Double_t getTotWalk(Double_t tot,Int_t ch){
  Double_t minp(0), walk(0), d(0), min(100);
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
  
  if(fabs(min)<0.8) walk=-min*tan(12*TMath::Pi()/180.);
  if(tot<10) walk-=(10-tot)*tan(6*TMath::Pi()/180.);

  // std::cout<<ch<<"  tot "<< tot << "  min "<<min << " minp "<< minp << "      corr "<< walk<<std::endl;
  // for(Int_t i=0; i<10; i++){
  //   std::cout<<i<<" "<<  gTotP[ch][i] <<std::endl;
  // }
  // std::cout<<" " <<std::endl;
  
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
  
  if(gcFile!=""){
    TFile f(gcFile);
    TIter nextkey(f.GetListOfKeys());
    TKey *key;

    while ((key = (TKey*)nextkey())) {
      TGraph *gr = (TGraph*)key->ReadObj();
      TString name = gr->GetName();
      Int_t channel = name.Atoi();
      Double_t x,y;
      if(channel == 10000){ // line calibration
	for(Int_t i=0; i<maxch; i++){
	  gr->GetPoint(i,x,y);
	  gMax[i] = (Int_t)(y+0.01);
	  std::cout<<"i  "<<i<< "  "<<  gMax[i]<<std::endl;
	  
	}
      }else if(channel == 20000){ // read tot offsets
	for(Int_t i=0; i<maxch; i++){
	  gr->GetPoint(i,x,y);
	  gTotO[i] = y;
	  std::cout<<"ch  "<<i<< " tot off "<<  gTotO[i]<<std::endl;
	}
      }else if(channel == 30000){ // read tot peaks
	for(Int_t i=0; i<960*10; i++){
	  gr->GetPoint(i,x,y);
	  gTotP[i/10][i%10] = y;
	  std::cout<<"ch  "<<i/10<< " peak "<< i%10<< " = " <<y<<std::endl;
	}
      }else{                      // spline calibration
	gGr[channel]= new TGraph(*gr);
      }
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
  //  if(entry >10000 ) return kTRUE;
  Int_t tdc,ch,tofpid(0);
  Double_t timeTot(0), grTime0(0), grTime1(0),timeLe(0),coarseTime(0),offset(0);
  Double_t time[10000], timeT[10000];
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
    
    coarseTime = 5*(Hits_nEpochCounter[i]*pow(2.0,11) + Hits_nCoarseTime[i]);
    if(gcFile!="") {
      //spline calib
      //time[i] = coarseTime-gGr[AddRefChannels(ch+1,tdc)]->Eval(Hits_nFineTime[i]);

      //linear calib
      Double_t max = (Double_t) gMax[AddRefChannels(ch+1,tdc)]-5;
      time[i] = coarseTime-5*(Hits_nFineTime[i]-31)/(max-31);
    }
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
	//	std::cout<<gTrigger << "  TRIG "<<(grTime1-grTime0)<<std::endl;
	
      }else {
	timeLe = time[i];
	if(gTrigger!=0 && ch<960) timeLe = timeLe - grTime1;
      }
      
      timeTot = timeT[i+1] - time[i]-gTotO[ch]+30;
      if(ch<960) timeLe += getTotWalk(timeTot,ch);
      
      if(gtFile!=""){
	if(ch<960) timeLe -= gGrDiff[ch]->Eval(timeTot);	
      }
      
      if(tofpid>0 && ch < 960){
	Double_t mom = 7;
	//	std::cout<<"m  "<<mass <<"  s  "<< 24.109/((mom/sqrt(mass*mass+mom*mom)*299792458)) *1E9 <<std::endl;
	timeLe += 24.109/((mom/sqrt(mass*mass+mom*mom)*299792458))*1E9;
      }

      if(gMode!=5 || tofpid!=0){
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
  gcFile = cFile; // fine time calibration
  gtFile = tFile; // pilas offsets + walk corrections
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
