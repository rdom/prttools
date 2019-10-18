// tcalibration - routine for the prtdirc data calibration 
// original author: Roman Dzhygadlo - GSI Darmstad

#include "datainfo.C"
#define TTSelector_cxx
#include "prttools.C"
#include "tcalibration.h"

DataInfo prt_data_info;
TString ginFile(""), goutFile(""), gcFile("");
Int_t gTrigger(0), gSetup(2019), gMode(0), gComboId(0),  gMaxIn[prt_maxchm];
Double_t tdcRefTime[prt_ntdcm],gTotO[prt_maxchm], gTotP[prt_maxdircch][10],gLeOffArr[prt_maxdircch],gEvtOffset(0);
TGraph *gGrIn[prt_maxchm], *gWalk[prt_maxchm], *gGrDiff[prt_maxchm];

Double_t walktheta(-5*TMath::Pi()/180.);

Double_t tof1le(0),tof2le(0),tof1tot(0),tof2tot(0);
Double_t fr11[11]={0,0.5,0.5,0.3,0.3,0.4, 0.3,0.3,0.2,0.20,0.15};
Double_t fr12[11]={0,1.0,1.0,0.9,0.9,0.9, 0.9,0.9,0.8,0.80,0.70};
Double_t fr21[11]={0,0.8,0.8,0.3,0.3,0.4, 0.3,0.3,0.2,0.2,0.2};
Double_t fr22[11]={0,1.0,1.0,0.9,0.9,0.9, 0.9,0.9,0.8,0.8,0.8};
Double_t c1y(0.5),c2y(0.5),c1x(0.9),c2x(0.9);

Double_t tof1lea[]= {0,0,81.84,76.48,74.51,73.58,73.06, 31.74, 71.91};
Double_t tof2lea[]= {0,0,72.12,72.03,71.99,71.96,71.94, 32.58, 72.55};
Double_t tof1tota[]={0,0,45.14,45.11,45.05,45.03,44.97, 47.24, 44.91};
Double_t tof2tota[]={0,0,45.26,45.31,45.32,45.34,45.36, 47.09, 45.38};

// //aug 2017
// Double_t tofpi1[]={0,0,71.50,31.50,0.00, 31.00,31.00, 31.00,31.00, 0, 31.0};
// Double_t tofpi2[]={0,0,72.50,32.20,0.00, 32.20,32.00, 31.95,31.80, 0 ,31.8};

// Double_t tofp1[] ={0,0,80.80,35.70,0.00, 33.00,32.60, 32.35,32.40, 0, 32.2};
// Double_t tofp2[] ={0,0,83.20,37.00,0.00, 34.00,33.50, 33.50,33.00, 0, 33.0};

//jul 2018
Double_t tofpi1a[]={0,0,71.50,31.50,0.00, 31.00,31.00, 32.00,31.00, 0, 31.0};
Double_t tofpi2a[]={0,0,72.50,32.20,0.00, 32.20,32.00, 33.50,31.80, 0 ,31.8};

Double_t tofp1a[] ={0,0,80.80,35.70,0.00, 33.00,32.60, 33.90,32.40, 0, 32.2};
Double_t tofp2a[] ={0,0,83.20,37.00,0.00, 34.00,33.50, 35.50,33.00, 0, 33.0};

Double_t tofpi1[]={0,0,35.00,35.00,35.00, 35.00,35.00, 35.00,35.00, 35.00, 35.0};
Double_t tofpi2[]={0,0,36.60,36.60,36.20, 36.30,36.20, 36.00,35.80, 35.80 ,35.6};

Double_t tofp1[] ={0,0,44.20,39.50,37.80, 37.00,36.60, 36.45,36.30, 36.25, 36.2};
Double_t tofp2[] ={0,0,47.50,41.50,39.20, 38.20,37.60, 37.20,37.00, 37.00, 36.7};

Int_t gg_nevents(0);

Bool_t IsPion(Double_t tof, Int_t mom){
  
  if(prt_data_info.getStudyId()<415) return tofpi1a[mom]<tof && tof<tofpi2a[mom];
  else return tofpi1[mom]<tof && tof<tofpi2[mom];
}

Bool_t IsProton(Double_t tof, Int_t mom){
  if(prt_data_info.getStudyId()<415) return tofp1a[mom]<tof && tof<tofp2a[mom];
  else return tofp1[mom]<tof && tof<tofp2[mom];
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
  gSetup = ((TObjString*)strobj->At(2))->GetString().Atoi();
  std::cout<<"gSetup "<<gSetup<<std::endl;
  
  prt_createMap(gSetup);
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

  {
    if(prt_data_info.getStudyId()==400
       || fileid.Contains("401_20C")
       || fileid.Contains("401_25C")
       || fileid.Contains("401_30C")){
      tofpi1[7]=31.00;
      tofpi2[7]=32.00;
      tofp1[7]=32.60;
      tofp2[7]=33.50;
    }
    if(fileid.Contains("401_90C")){
      tofpi1[7]=31.00;
      tofpi2[7]=31.9;
      tofp1[7]=32.30;
      tofp2[7]=33.50;
    }
  }
  
  if(prt_data_info.getStudyId()<0) momentum=7;
  c1y=fr11[momentum];
  c2y=fr21[momentum];
  c1x=fr12[momentum];
  c2x=fr22[momentum];
   
  std::cout<<"Initialization successful"<<std::endl;
}  

Bool_t TTSelector::Process(Long64_t entry){
  // if(gg_nevents >= 50000 ) return kTRUE;
  Int_t tdc,ch,tofpid(0);
  Double_t grTime0(0), grTime1(0),grTime2(0),coarseTime(0),offset(0),triggerLe(0),triggerTot(0);
  Double_t time[10000], timeLe(0),timeT[10000],timeTot(0),mom(7),simOffset(74.20);
  Int_t multT1(0), multT2(0), multT3v(0), multT3h(0), multTof1(0), multTof2(0),
    multStr1(0),multStl1(0),multStr2(0),multStl2(0);
  
  TString current_file_name  = TTSelector::fChain->GetCurrentFile()->GetName();
  Bool_t trbdata = current_file_name.Contains("trb");
  Bool_t laser = current_file_name.Contains("pilas") || current_file_name.Contains("pico");
    
  TObjArray *sarr = current_file_name.Tokenize("_");

  if(entry%10000==0) std::cout<<"event # "<< entry <<std::endl;
  GetEntry(entry);

  // Int_t trigT1(816);
  // Int_t trigT2(817);
  // Int_t trigT3h(818);
  // Int_t trigT3v(819);  
  // Int_t trigTof1(1392);
  // Int_t trigTof2(1398);

  Int_t trigT1(520);
  Int_t trigT2(513);
  Int_t trigT3h(514);
  Int_t trigT3v(515);  
  Int_t trigTof1(1136);
  Int_t trigTof2(1138);

  Int_t trigStr1(1140);
  Int_t trigStl1(1142);
  Int_t trigStr2(1144);
  Int_t trigStl2(1146);

    
  fEvent = new PrtEvent();
  if(gMode==5){
    Int_t studyId=prt_data_info.getStudyId();
    
    if(studyId>0) {
      mom = prt_data_info.getMomentum();
      simOffset = prt_data_info.getSimTO();
      fEvent->SetAngle(prt_data_info.getAngle());
      fEvent->SetPhi(prt_data_info.getPhi());
      fEvent->SetMomentum(TVector3(0,0,mom));
      fEvent->SetTrigger(816);
      fEvent->SetGeometry(studyId);
      fEvent->SetLens(prt_data_info.getLensId());
      fEvent->SetPrismStepX(prt_data_info.getXstep());
      fEvent->SetPrismStepY(prt_data_info.getYstep());
      fEvent->SetBeamX(prt_data_info.getX());
      fEvent->SetBeamZ(prt_data_info.getZ());

      Double_t rad = TMath::Pi()/180.0,
	zrot=146,
	xrot=100,
	prtangle=fEvent->GetAngle(),
	z=fEvent->GetBeamZ(),
	b = xrot*tan(0.5*(prtangle-90)*rad),
	lenz = (z-zrot+b)/cos((prtangle-90)*rad)+b+zrot;

      if(fEvent->GetLens()==6) lenz-=12;
      if(fEvent->GetLens()==3) lenz-=15;
      if(fEvent->GetLens()==2) lenz-=14.4;
      fEvent->SetPosition(TVector3(0,0,lenz));
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
	if(prt_getTdcId(gTrigger)==tdc) grTime0 = time[i];
      }
      if(ch==trigT1) multT1++; //trigger1
      if(ch==trigT2) multT2++; //trigger2
      if(ch==trigTof1) multTof1++; //tof1
      if(ch==trigTof2) multTof2++; //tof2
      if(ch==trigT3h) multT3h++; //trigger3h
      if(ch==trigT3v) multT3v++; //trigger3v
      
      if(ch==trigStr1) multStr1++;
      if(ch==trigStl1) multStl1++;
      if(ch==trigStr2) multStr2++;
      if(ch==trigStl2) multStl2++;
    }else{
      timeT[i]=time[i];
      if(ch==gTrigger && grTime2==0) grTime2=time[i];
    }
  }
  
  Double_t tof1(0),tof2(0),tot1(0),tot2(0),toftime(0),mass(0);
  Double_t tofstr1(0),tofstr2(0),totstr1(0),totstr2(0),
		   tofstl1(0),tofstl2(0),totstl1(0),totstl2(0);
  
  if(gMode==5){
    if(multT1<1 || multTof1<1 || multTof2<1 || multT3h<1 || multT3v<1){ //  || mult2!=1 || mult5!=1
    //if(multT1<1 || multTof1<1 || multTof2<1 || multStr1<1 ||  multStl1<1|| multStr2<1 ||  multStl2<1 || multT3h<1 || multT3v<1){
    //if(multT1<1 || multTof1<1 || multTof2<1){ //  || mult2!=1 || mult5!=1
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
      if(ch==trigStr1 && tofstr1==0){
    	tofstr1 = time[i]-tdcRefTime[tdc];
    	totstr1 = timeT[i+1] - time[i];
      }
      if(ch==trigStr2 && tofstr2==0){
    	tofstr2 = time[i]-tdcRefTime[tdc];
    	totstr2 = timeT[i+1] - time[i];
      }
      if(ch==trigStl1 && tofstl1==0){
    	tofstl1 = time[i]-tdcRefTime[tdc];
    	totstl1 = timeT[i+1] - time[i];
      }
      if(ch==trigStl2 && tofstl2==0){
    	tofstl2 = time[i]-tdcRefTime[tdc];
    	totstl2 = timeT[i+1] - time[i];
      }
    }

    if(tof1!=0 && tof2!=0) {
      Double_t time = tof2-tof1;
      // time += (tot1-tof1tot)*tan(walktheta);
      // time += (tot2-tof2tot)*tan(-walktheta);

      // time += (tot1-tof1tot)*tan(-3*TMath::Pi()/180.);
      // time += (tot2-tof2tot)*tan(2*TMath::Pi()/180.);

      // time += (tot1-41.32)*tan(-4*TMath::Pi()/180.);
      // time += (tot2-40.75)*tan(2*TMath::Pi()/180.);

      toftime = time;
      Int_t m = (Double_t) (mom+0.1);

      if(IsPion(time,m)){
    	tofpid=211;
    	mass=0.13957018;
      }else if(IsProton(time,m)){ 
    	tofpid=2212;
    	mass = 0.938272046;
      }else{
      	fEvent->Clear();
      	delete fEvent;
      	return kTRUE;
      }
    }

    // if(tofstr1!=0 && tofstr2!=0 && tofstl1!=0 && tofstl2!=0 ) { 
    // tofpi1[7]=32.0;
    // tofpi2[7]=33.3;      
    // tofp1[7] =33.7;
    // tofp2[7] =35;	
    // Double_t time =(tofstr2+tofstl2)/2.-(tofstr1+tofstl1)/2.;

    // toftime = time;
    // Int_t m = (Double_t) (mom+0.1);

    // if(IsPion(time,m)){
    // 	tofpid=211;
    // 	mass=0.13957018;
    // }else if(IsProton(time,m)){ 
    // 	tofpid=2212;
    // 	mass = 0.938272046;
    // }else{
    // 	fEvent->Clear();
    // 	delete fEvent;
    // 	return kTRUE;
    // }
    // }
  }
  
  PrtHit hit;
  Int_t nrhits=0;
  
  if((grTime0>0 && grTime1>0) || gTrigger==0){    
    if(gTrigger!=0) {
      triggerLe = grTime1 - grTime0;
      if(gTrigger==1140) triggerLe = (tofstr1+tofstl1)/2.;
      if(gTrigger==1144) triggerLe = (tofstr2+tofstl2)/2.;
      triggerTot= grTime2 - grTime1;
    }
    
    for(Int_t i=0; i<Hits_ && i<10000; i++){      
      //if(Hits_nTdcErrCode[i]!=0) continue;
      if(Hits_nTdcChannel[i]==0) continue; // ref channel
      if(Hits_nSignalEdge[i]==0) continue; // tailing edge 
      
      tdc = map_tdc[Hits_nTrbAddress[i]];
      ch = prt_getChannelNumber(tdc,Hits_nTdcChannel[i])-1;
      
      //if(!trbdata && prt_isBadChannel(ch)) continue;

      if(gMode>0){
	timeLe = time[i]-tdcRefTime[tdc];
	if(gTrigger!=0 && ch<prt_maxdircch) timeLe = timeLe - triggerLe;
      }else {
	timeLe = time[i];
	if(gTrigger!=0 && ch<prt_maxdircch) timeLe = timeLe - grTime1;
      }
      
      timeTot = timeT[i+1] - time[i];
      if(ch<prt_maxdircch){
	timeTot += 30-gTotO[ch];

	//timeLe += getTotWalk(timeTot,ch);
	//if(gTrigger==trigT1 && fabs(triggerTot-tof1tot)<1) timeLe -= (triggerTot-tof1tot)*tan(5*TMath::Pi()/180.);

	if(fabs(tot2-tof2tot)<1.2) timeLe -= (tot2-tof2tot)*tan(2*TMath::Pi()/180.); //7.1;
	
        if(timeTot>0.5 && timeTot<9 && gWalk[ch]) timeLe -=  gWalk[ch]->Eval(timeTot);	
	timeLe -= gLeOffArr[ch];

	if(!laser && gMode==5){
	  Double_t rad = TMath::Pi()/180.,
	    zrot=155,
	    xrot=98,
	    prtangle= fEvent->GetAngle(),
	    z = fEvent->GetBeamZ(),
	    rot_dist = ((z-zrot)-xrot/cos(prtangle*rad))*tan((90-prtangle)*rad)/1000.;

	  if(gTrigger==trigT1)   timeLe -= ( 2.719+rot_dist)/((mom/sqrt(mass*mass+mom*mom)*299792458))*1E9; // trig1	
	  if(gTrigger==trigTof1) timeLe -= (24.490+rot_dist)/((mom/sqrt(mass*mass+mom*mom)*299792458))*1E9; // tof1
	  if(gTrigger==1140)     timeLe -= (24.490+0.030+rot_dist)/((mom/sqrt(mass*mass+mom*mom)*299792458))*1E9; // tof1
	  if(gTrigger==trigTof2) timeLe += ( 4.151-rot_dist)/((mom/sqrt(mass*mass+mom*mom)*299792458))*1E9; // tof2
	  if(gTrigger==1144)     timeLe += ( 4.151-rot_dist)/((mom/sqrt(mass*mass+mom*mom)*299792458))*1E9; // scitil2
	  
	  timeLe += simOffset;
	  if(gTrigger==trigTof2) timeLe -= 59;
	}	
	if(!laser && gMode!=5) timeLe += -2.5; //simOffset;
      }
      
      if(gMode==5){
	if(fabs(timeLe)>600) continue;
	// if(ch==trigTof1) timeLe -= (tot1-tof1tot)*tan(-3*TMath::Pi()/180.);
	// if(ch==trigTof2) timeLe -= (tot2-tof2tot)*tan(2*TMath::Pi()/180.);
	
	////timeLe-=gEvtOffset;

	// if(ch>820 && ch<1340) continue;
        // if(ch<prt_maxdircch && (timeLe<0 || timeLe>100)) continue;
      }

      if(gMode!=5 || tofpid!=0){

	if(prt_geometry==2023){
	  if(timeTot<0 || timeTot>20) continue;
	  if(timeLe<-100 || timeLe>100) continue;
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
  }
  
  if(nrhits!=0) {
    gg_nevents++;
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

void tcalibration(TString inFile= "../../data/cj.hld.root", TString outFile= "outFileC.root", TString cFile= "calib.root", TString tFile= "calibOffsets.root", Int_t trigger=0,  Int_t sEvent =0, Int_t eEvent=0, Int_t mode=1, Int_t build=0, Int_t setupid=2019){
  if(build==1) return;
  ginFile = inFile;
  goutFile = outFile;
  gcFile = (cFile!="")? cFile: "0"; // calibration
  gTrigger = trigger;
  gMode = mode;  
  gSetup = setupid;
  if(gMode == 5) gTrigger=1136; //1136
  
  TChain* ch = new TChain("T");
  ch->Add(ginFile);
  
  Int_t entries = ch->GetEntries();
  TTSelector *selector = new TTSelector();  
  TString option = Form("%d %d %d",gTrigger,gMode,gSetup);
  
  if(eEvent==0){
    std::cout<<"Entries in chain:  "<< entries<<std::endl;
    ch->Process(selector,option,entries);
  }else{
    ch->Process(selector,option,eEvent-sEvent,sEvent); 
  }
}
