// tcalibration - routine for the prtdirc data calibration 
// original author: Roman Dzhygadlo - GSI Darmstad

#include "datainfo.C"
#define TTSelector_cxx
#include "prttools.C"
#include "tcalibration.h"

DataInfo prt_data_info;
TString ginFile(""), goutFile(""), gcFile("");
Int_t gSetup=2015, gTrigger(0), gMode(0), gComboId(0),  gMaxIn[prt_maxch];
Double_t tdcRefTime[prt_ntdc],gTotO[prt_maxch], gTotP[prt_maxdircch][10],gLeOffArr[prt_maxdircch],gEvtOffset(0);
TGraph *gGrIn[prt_maxch], *gWalk[prt_maxch], *gGrDiff[prt_maxch];

Double_t walktheta(-5*TMath::Pi()/180.);

Double_t tof1le(0),tof2le(0),tof1tot(0),tof2tot(0);
Double_t fr11[11]={0,0.5,0.5,0.3,0.3,0.4, 0.3,0.3,0.2,0.20,0.15};
Double_t fr12[11]={0,1.0,1.0,0.9,0.9,0.9, 0.9,0.9,0.8,0.80,0.70};
Double_t fr21[11]={0,0.8,0.8,0.3,0.3,0.4, 0.3,0.3,0.2,0.2,0.2};
Double_t fr22[11]={0,1.0,1.0,0.9,0.9,0.9, 0.9,0.9,0.8,0.8,0.8};
Double_t c1y(0.5),c2y(0.5),c1x(0.9),c2x(0.9);

Double_t tof1lea[]= {0,0,81.84,76.48,74.51,73.58,73.06, 71.17, 71.91};
Double_t tof2lea[]= {0,0,72.12,72.03,71.99,71.96,71.94, 72.00, 72.55};
Double_t tof1tota[]={0,0,45.14,45.11,45.05,45.03,44.97, 47.10, 44.91};
Double_t tof2tota[]={0,0,45.26,45.31,45.32,45.34,45.36, 46.50, 45.38};

//aug 2017
Double_t tofpi1[]={0,0,71.50,71.50,71.60, 71.50,71.55, 70.00,71.60};
Double_t tofpi2[]={0,0,72.50,72.50,72.40, 72.20,72.20, 71.40,72.10};

Double_t tofp1[] ={0,0,80.80,75.80,73.80, 73.20,72.65, 71.80,72.30};
Double_t tofp2[] ={0,0,83.20,77.20,75.00, 73.90,73.40, 73.20,72.90};
   
Bool_t IsPion(Double_t tof, Int_t mom){
  return tofpi1[mom]<tof && tof<tofpi2[mom];
}

Bool_t IsProton(Double_t tof, Int_t mom){
  return tofp1[mom]<tof && tof<tofp2[mom];
}

Bool_t insideOfEllipce(Double_t x, Double_t y, Double_t x0, Double_t y0,  Double_t r1, Double_t r2, Double_t w=0){

  Double_t xx = cos(w)*(x-x0)+sin(w)*(y-y0);
  Double_t yy = sin(w)*(x-x0)-cos(w)*(y-y0);

  return xx*xx/(r1*r1)+yy*yy/(r2*r2)<=1;
}

Double_t getTotWalk(Double_t tot,Int_t ch, Int_t type=0){
  Double_t minp(0), walk(0), d(0), min(100);

  if(type==0){
    if(ch<prt_maxdircch){
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
    if(ch/48==1) wcorr=5;
    if(ch/48==3) wcorr=15;
    if(ch/48==5) wcorr=12;
    // if(ch/48==9) wcorr=0;
      
    if(fabs(min)<0.8) walk=-min*tan(wcorr*TMath::Pi()/180.);
    if(tot<8) walk-=(4-tot)*tan(10*TMath::Pi()/180.);
  }

  if(type==1){ //walk of the xxx
    walk += (38.85-tot)*tan(25*TMath::Pi()/180.); 
  }
  
  return walk;
}

void TTSelector::Begin(TTree *){
  TString option = GetOption();
  TObjArray *strobj = option.Tokenize(" ");
  gTrigger = ((TObjString*)strobj->At(0))->GetString().Atoi();
  gMode = ((TObjString*)strobj->At(1))->GetString().Atoi();
  prt_createMap();
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
      Double_t x,y;
      if(name.Contains("tof")){
	name.Remove(0,4);
	//	if(ginFile.Contains(name)){
	  gr->GetPoint(0,tof1le,tof2le);
	  gr->GetPoint(1,tof1tot,tof2tot);	  
	  //	}
	continue;
      }
      if(name.Contains("off")){ // read event offsets
	name.Remove(0,4);
	if(ginFile.Contains(name)) gr->GetPoint(0,x,gEvtOffset);
	continue;
      }

      if(name.Contains("walk")){ // read walk corrections
	name.Remove(0,5);
	Int_t ch = name.Atoi();
	gWalk[ch]= new TGraph(*gr);
	continue;
      }
      
      long long  ch = name.Atoll();
      if(ch <10000){ // spline calibration
	gGrIn[ch]= new TGraph(*gr);
      }else if(ch == 10000){ // line calibration
	for(Int_t i=0; i<prt_maxch; i++){
	  gr->GetPoint(i,x,y);
	  gMaxIn[i] = (Int_t)(y+0.01);
	  //std::cout<<"ch  "<<i<< "  FT max"<<  gMaxIn[i]<<std::endl;	  
	}
      }else if(ch == 10001){ // read tot offsets
	for(Int_t i=0; i<prt_maxch; i++){
	  gr->GetPoint(i,gTotO[i],y);
	  //std::cout<<"ch  "<<i<< " TOT off "<<  gTotO[i]<<std::endl;
	}
      }else if(ch == 10002){ // read tot peaks
	for(Int_t i=0; i<prt_maxdircch*10; i++){
	  gr->GetPoint(i,x,y);
	  gTotP[i/10][i%10] = y;
	  //std::cout<<"ch  "<<i/10<< " peak "<< i%10<< " = " <<y<<std::endl;
	}
      }else if(ch == 10003){ // read LE offsets 1
	for(Int_t i=0; i<prt_maxdircch; i++){
	  gr->GetPoint(i,gLeOffArr[i],y);
	}
      }
    }
    f.Close();
  }
  
  TString fileid(ginFile);
  fileid.Remove(0,fileid.Last('/')+1);
  //fileid.Remove(fileid.Last('.')-4);
  fileid.ReplaceAll("C.hld.root","");
  fileid.ReplaceAll("P.hld.root","");
  fileid.ReplaceAll(".hld.root","");
  prt_data_info = getDataInfo(fileid);
  Int_t momentum = prt_data_info.getMomentum();
  std::cout<<fileid<<" study id "<<prt_data_info.getStudyId() << " mom "<<momentum <<std::endl;
   
  if(prt_data_info.getStudyId()<0) momentum=7;
  c1y=fr11[momentum];
  c2y=fr21[momentum];
  c1x=fr12[momentum];
  c2x=fr22[momentum];
   
  std::cout<<"Initialization successful"<<std::endl;
}  

Bool_t TTSelector::Process(Long64_t entry){
  // if(entry >1000 ) return kTRUE;
  Int_t tdc,ch,tofpid(0);
  Double_t grTime0(0), grTime1(0),grTime2(0),coarseTime(0),offset(0),triggerLe(0),triggerTot(0);
  Double_t time[10000], timeLe(0),timeT[10000],timeTot(0),mom(7),simOffset(12.59-59);
  Int_t multT1(0), multT2(0), multT3v(0), multT3h(0), multTof1(0), multTof2(0);
  
  TString current_file_name  = TTSelector::fChain->GetCurrentFile()->GetName();
  Bool_t trbdata = current_file_name.Contains("trb");
  Bool_t laser = current_file_name.Contains("pilas") || current_file_name.Contains("pico");
    
  TObjArray *sarr = current_file_name.Tokenize("_");

  if(entry%10000==0) std::cout<<"event # "<< entry <<std::endl;
  GetEntry(entry);

  Int_t trigT1(816);
  Int_t trigT2(817);
  Int_t trigT3h(818);
  Int_t trigT3v(819);  
  Int_t trigTof1(1392);
  Int_t trigTof2(1398);
    
  fEvent = new PrtEvent();
  if(gMode==5){
    Int_t studyId=prt_data_info.getStudyId();
    if(studyId>0) {
      mom = prt_data_info.getMomentum();
      simOffset = prt_data_info.getSimTO();
      fEvent->SetAngle(prt_data_info.getAngle());
      fEvent->SetMomentum(TVector3(0,0,mom));
      fEvent->SetTrigger(816);
      fEvent->SetGeometry(studyId);
      fEvent->SetLens(prt_data_info.getLensId());
      fEvent->SetPrismStepX(prt_data_info.getXstep());
      fEvent->SetPrismStepY(prt_data_info.getYstep());
      fEvent->SetBeamX(prt_data_info.getX());
      fEvent->SetBeamZ(prt_data_info.getZ());
    }
  }
  
  for(Int_t i=0; i<Hits_ && i<10000; i++){
    tdc = map_tdc[Hits_nTrbAddress[i]];
    if(tdc<0) continue;
    ch = prt_getChannelNumber(tdc,Hits_nTdcChannel[i])-1;
    time[i] = 5*(Hits_nEpochCounter[i]*pow(2.0,11) + Hits_nCoarseTime[i]); //coarsetime
    if(gcFile!="0") {
      //spline calib
      //time[i] -= gGrIn[prt_addRefChannels(ch+1,tdc)]->Eval(Hits_nFineTime[i]+1); //slow
      Double_t xx,yy;
      gGrIn[prt_addRefChannels(ch+1,tdc)]->GetPoint(Hits_nFineTime[i],xx,yy); time[i] -=yy;//fast

      //linear calib
      // Double_t max = (Double_t) gMaxIn[prt_addRefChannels(ch+1,tdc)]-2;
      // time[i] = coarseTime-5*(Hits_nFineTime[i]-31)/(max-31);
    } // else time[i] -= (Hits_nFineTime[i]-31)*0.0102;

    if(Hits_nSignalEdge[i]==1){
      if(ch==gTrigger && grTime1==0) grTime1 = time[i];
      if(Hits_nTdcChannel[i]==0) { //ref channel
	tdcRefTime[tdc] = time[i];
	if(gTrigger/48==tdc) grTime0 = time[i];
      }
      if(ch==trigT1) multT1++; //trigger1
      if(ch==trigT2) multT2++; //trigger2
      if(ch==trigTof1) multTof1++; //tof1
      if(ch==trigTof2) multTof2++; //tof2
      if(ch==trigT3h) multT3h++; //trigger3h
      if(ch==trigT3v) multT3v++; //trigger3v
    }else{
      timeT[i]=time[i];
      if(ch==gTrigger && grTime2==0) grTime2=time[i];
    }
  }
  
  Double_t tof1(0),tof2(0),tot1(0),tot2(0),toftime(0),mass(0);

  
  if(gMode==5){
    if(multT1!=1 || multTof1<1 || multTof2<1 || multT3h<1 || multT3v<1){ //  || mult2!=1 || mult5!=1
      fEvent->Clear();
      delete fEvent;
      return kTRUE;
    }
 
    for(Int_t i=0; i<Hits_ && i<10000; i++){
      if(Hits_nTdcErrCode[i]!=0) continue;
      if(Hits_nTdcChannel[i]==0) continue; // ref channel
      if(Hits_nSignalEdge[i]==0) continue; // tailing edge
      
      tdc = map_tdc[Hits_nTrbAddress[i]];
      ch = prt_getChannelNumber(tdc,Hits_nTdcChannel[i])-1;
      if(ch==trigTof1 && tof1==0){
    	tof1 = time[i]-tdcRefTime[tdc];
    	tot1 = timeT[i+1] - time[i];
      }
      if(ch==trigTof2 && tof2==0){
    	tof2 = time[i]-tdcRefTime[tdc];
    	tot2 = timeT[i+1] - time[i];
      }      
    }

    if(tof1!=0 && tof2!=0) {
      Double_t time = tof2-tof1;
      time += (tot1-tof1tot)*tan(walktheta);
      time += (tot2-tof2tot)*tan(-walktheta);

      // time += (tot1-tof1tot)*tan(-6*TMath::Pi()/180.);
      // time += (tot2-tof2tot)*tan(0.5*TMath::Pi()/180.);      
      
      toftime = time;
      Int_t m = (Double_t) (mom+0.1);

      // if(insideOfEllipce(time, tot1, tof1lea[m], tof1tota[m], c1y, c1x) && insideOfEllipce(time, tot2, tof1lea[m], tof2tota[m], c1y, c1x)){
      if(IsPion(time,m)){
	tofpid=211;
	mass=0.13957018;
	//}else if(insideOfEllipce(time, tot1, tof2lea[m], tof1tota[m], c2y, c2x) && insideOfEllipce(time, tot2, tof2lea[m], tof2tota[m], c2y, c2x)){
      }else if(IsProton(time,m)){ 
	tofpid=2212;
	mass = 0.938272046;
      }else{
	fEvent->Clear();
	delete fEvent;
	return kTRUE;
      }
    }
  }
  
  PrtHit hit;
  Int_t nrhits=0;
  if((grTime0>0 && grTime1>0) || gTrigger==0){
    if(gTrigger!=0) {
      triggerLe = grTime1 - grTime0;
      triggerTot= grTime2 - grTime1;
    }
    
    for(Int_t i=0; i<Hits_ && i<10000; i++){
      if(Hits_nTdcErrCode[i]!=0) continue;
      if(Hits_nTdcChannel[i]==0) continue; // ref channel
      if(Hits_nSignalEdge[i]==0) continue; // tailing edge
      
      tdc = map_tdc[Hits_nTrbAddress[i]];
      ch = prt_getChannelNumber(tdc,Hits_nTdcChannel[i])-1;
      if(!trbdata && prt_isBadChannel(ch)) continue;
      
      if(gMode>0){
	timeLe = time[i]-tdcRefTime[tdc];
	if(gTrigger!=0 && ch<prt_maxdircch) timeLe = timeLe - (grTime1-grTime0);
      }else {
	timeLe = time[i];
	if(gTrigger!=0 && ch<prt_maxdircch) timeLe = timeLe - grTime1;
      }
      
      timeTot = timeT[i+1] - time[i];

      if(ch<prt_maxdircch) {
	//if(timeTot<0 || timeLe<20 || timeLe>40) continue;
	timeTot += 30-gTotO[ch];

	// tmp commented! 
	// timeLe += getTotWalk(timeTot,ch);
	// if(gTrigger==trigT1 && fabs(triggerTot-tof1tot)<1) timeLe -= (triggerTot-tof1tot)*tan(5*TMath::Pi()/180.);

        if(timeTot>0.5 && timeTot<9 && gWalk[ch]) timeLe -=  gWalk[ch]->Eval(timeTot);
		
	timeLe -= gLeOffArr[ch];

	if(!laser && gMode==5){
	  if(gTrigger==trigT1) timeLe -= (7.829 +0.39)/((mom/sqrt(mass*mass+mom*mom)*299792458))*1E9; //25 degree trig1	
	  if(gTrigger==trigTof1) timeLe -= (24.460+0.39)/((mom/sqrt(mass*mass+mom*mom)*299792458))*1E9; //25 degree tof1
	  if(gTrigger==trigTof2) timeLe += ( 3.998-0.39)/((mom/sqrt(mass*mass+mom*mom)*299792458))*1E9; //25 degree tof2

	  timeLe += simOffset;
	}else{
	  if(fabs(tof1tot-44.9)<1) timeLe -= (tof1tot-44.9)/8.4; //7.1;
	}
      }   
      
      if(gMode==5){
	//timeLe-=gEvtOffset;
	//if(ch>prt_maxdircch && ch != 1104 && ch != 1344 && ch != 1248) continue;
	if(ch<prt_maxdircch && (timeLe<0 || timeLe>100)) continue;
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
  gcFile = (cFile!="")? cFile: "0"; // calibration
  gTrigger = trigger;
  gMode = mode;
  if(gMode == 5) gTrigger=816; //720;
  
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
