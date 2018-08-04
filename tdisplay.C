// tdisplay - tool to plot different quantities from the .hld file
// original author: Roman Dzhygadlo - GSI Darmstadt 

#define TTSelector_cxx
#include "prttools.C"
#include "tdisplay.h"

const Int_t maxfiles = 150;
const Int_t geometry=2018;
Int_t gSetup=2015, gTrigger,gEntries=0, gMode=0, gComboId=0, gWorkers=4, nfiles = 10;
TString ginFile(""),gcFile(""), fileList[maxfiles];
TH1F *hCh, *hRefDiff, *hFine[maxfiles][prt_maxch], *hTot[maxfiles][prt_maxch],
  *hTimeL[prt_nmcp][prt_npix], *hTimeT[prt_nmcp][prt_npix], *hSigma[prt_ntdc];
TH2F *hLeTot[prt_maxch], *hShape[prt_nmcp][prt_npix];
TGraph *gMaxFine, *gLeOff, *gTotOff, *gTotPeaks;
Int_t gMaxIn[prt_maxch];
Double_t  tdcRefTime[prt_ntdc], gTotO[prt_maxch], gTotP[prt_maxdircch][10],gLeOffArr[prt_maxdircch],gEvtOffset(0);
TGraph *gGrIn[prt_maxch], *gLeO[prt_maxch], *gGr[maxfiles][prt_maxch], *gGrDiff[prt_maxch];
TCanvas *cTime;

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
    
    Double_t wcorr(8);
    if(ch/48==1) wcorr=0;
    if(ch/48==5) wcorr=0;
    if(ch/48==9) wcorr=0;
    if(ch/48==13) wcorr=0;
      
    if(fabs(min)<0.8) walk=-min*tan(wcorr*TMath::Pi()/180.);
    if(tot<10) walk-=(10-tot)*tan(6*TMath::Pi()/180.); //10
  }

  if(type==1){ //walk of the 1345
    walk += (38.85-tot)*tan(25*TMath::Pi()/180.);
  }
  
  return walk;
}

void TTSelector::SlaveBegin(TTree *){
  std::cout<<"init starts "<<std::endl;
  
  TString option = GetOption();
  TObjArray *strobj = option.Tokenize(" ");
  nfiles = ((TObjString*)strobj->At(0))->GetString().Atoi();
  gTrigger = ((TObjString*)strobj->At(1))->GetString().Atoi();
  gMode = ((TObjString*)strobj->At(2))->GetString().Atoi();
  gSetup = ((TObjString*)strobj->At(3))->GetString().Atoi();
  gcFile = ((TObjString*)strobj->At(4))->GetString();
  
  for(Int_t i=5; i<nfiles+5; i++){
    fileList[i-5]=((TObjString*)strobj->At(i))->GetString();
    std::cout<<" fileList[i]  "<<fileList[i-5] <<std::endl;
  }

  Int_t totb(200), totl(0), toth(100); //12
  Int_t leb(400), le1(20), le2(40);
  if(fileList[0].Contains("trb")){
    gTrigger=0;
    //totb=4000; totl=-1000; toth=1000;
    // leb=4000, le1=-5, le2=100;

    totb=4000; totl=50; toth=80;
    leb=4000, le1=-5, le2=5;
  }
  if(fileList[0].Contains("pilas") || fileList[0].Contains("th_") ){
    if(!gTrigger) gTrigger=820;
    if(geometry==2018) gTrigger=520;
    // leb = 2000; le1=40; le2=80;
    leb = 2000; le1=55; le2=85;
    //  leb = 2000; le1=75; le2=100;
    totb=240; totl=0; toth=12;
  }
  if(fileList[0].Contains("pico")){
   if(!gTrigger) gTrigger=820;
    totb=240; totl=0; toth=12;
    leb = 400; le1=10; le2=40; //20 40
  }
  
  for(Int_t j=0; j<nfiles; j++){
    for(Int_t c=0; c<prt_maxch; c++){
      hFine[j][c] = new TH1F(Form("hFine_%d_ch%d",j,c),Form("hFine_%d_ch%d;bin [#];LE entries [#]",j,c) , 600,0,600);
      hTot[j][c] = new TH1F(Form("hTot_%d_ch%d",j,c), Form("hTot_%d_ch%d;TOT [ns];entries [#]",j,c) , totb,totl,toth); //2000  50 80 //400 -2 2
      hTot[j][c]->SetLineColor(j+1);
      fOutput->Add(hFine[j][c]);
      fOutput->Add(hTot[j][c]);
    }
  }
  
  for(Int_t c=0; c<prt_maxch; c++){
    hLeTot[c] = new TH2F(Form("hLeTot_ch%d",c), Form("hLeTot_ch%d",c) ,1000,le1,le2, 100,0,12); //35 40 for 1345
    fOutput->Add(hLeTot[c]);
  }

  prt_initDigi();
  for(Int_t m=0; m<prt_nmcp; m++){
    for(Int_t p=0; p<prt_npix; p++){

      hTimeL[m][p] = new TH1F(Form("hTimeL_mcp%dpix%d",m,p), Form("hTimeL_%d_%d;LE time [ns];entries [#]",m,p) , leb,le1,le2);
      hTimeT[m][p] = new TH1F(Form("hTimeT_mcp%dpix%d",m,p), Form("hTimeT_%d_%d;TOT [ns];entries [#]",m,p) , 400,le1,le2);
      hShape[m][p] = new TH2F(Form("hShape_mcp%dpix%d",m,p), Form("hShape_%d_%d;time [ns];offset to the threshold [mV]",m,p) , 400,le1,le2,130,-5,35);

      fOutput->Add(hTimeL[m][p]);
      fOutput->Add(hTimeT[m][p]);
      fOutput->Add(hShape[m][p]);
    }
    fOutput->Add(prt_hdigi[m]);
  }

  hCh = new TH1F("hCh","hCh;channel [#];entries [#]",prt_maxch,0,prt_maxch);
  fOutput->Add(hCh);
  hRefDiff = new TH1F("hRefDiff","ch-ref. resolution;sigma [ns];entries [#]",200,0,1);//0.05
  for(Int_t i=0; i<prt_ntdc; i++){
    hSigma[i] = new TH1F(Form("hSigma%d",i),"ch-ref. resolution;sigma [ns];entries [#]",200,0,0.1);
    fOutput->Add(hSigma[i]);
  }
  fOutput->Add(hRefDiff);
  prt_createMap();
  
  if(gcFile!="0"){
    TFile f(gcFile);
    TIter nextkey(f.GetListOfKeys());
    TKey *key;
    
    while ((key = (TKey*)nextkey())) {
      TGraph *gr = (TGraph*)key->ReadObj();
      TString name = gr->GetName();
      Double_t x,y;
      if(name.Contains("tof")){
	continue;
      }
      if(name.Contains("off")){ // read event offsets
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
      }else if(ch >= 20000 && ch < 30000){ // read LE offsets 3
	gLeO[ch-20000] = new TGraph(*gr);
      }
    }
    f.Close();
  }

  std::cout<<"init done " <<std::endl;
}

Int_t mult[prt_maxch]={0};
Bool_t TTSelector::Process(Long64_t entry){
  if(entry%1000==0) std::cout<<"event # "<< entry <<std::endl;
  Int_t tdc,ch,mcp,pix,col,row;
  Double_t triggerLe(0),triggerTot(0), grTime0(0), grTime1(0),grTime2(0),timeLe(0), timeTe(0), offset(0);
  const Int_t maxhits(100000);
  Double_t timeL[maxhits], timeT[maxhits];

  TString current_file_name  = TTSelector::fChain->GetCurrentFile()->GetName();
  TObjArray *sarr = current_file_name.Tokenize("_");
  if(sarr->At(1)){
    TString soffset = ((TObjString *) sarr->At(1))->GetName();
    offset = soffset.Atof();
    offset = offset/400.;
  }
  
  Int_t fileid=0;
  if(gMode<3){
    for(fileid=0; fileid<nfiles; fileid++){
      if(current_file_name.Contains(fileList[fileid])) break;
    }
  }
  if(gMode==4) if(current_file_name.Contains("cc"))  fileid=1;
  Bool_t calib = current_file_name.Contains("trb");

  GetEntry(entry);
  memset(mult, 0, sizeof(mult));
  Int_t hh(0);
  for(auto i=0; i<Hits_&& i<maxhits; i++){
    tdc = map_tdc[Hits_nTrbAddress[i]];

    if(tdc<0) continue;
    ch = prt_getChannelNumber(tdc,Hits_nTdcChannel[i])-1;
    //if(Hits_nTdcErrCode[i]!=0) continue;
    timeL[i] =  5*(Hits_nEpochCounter[i]*pow(2.0,11) + Hits_nCoarseTime[i]); //coarsetime
    if(gcFile!="0"){
      //spline calib
      //timeL[i]-= gGrIn[prt_addRefChannels(ch+1,tdc)]->Eval(Hits_nFineTime[i]+1);
      Double_t xx,yy;
      gGrIn[prt_addRefChannels(ch+1,tdc)]->GetPoint(Hits_nFineTime[i],xx,yy); timeL[i] -=yy;//fast
      
      //linear calib
      // Double_t max = (Double_t) gMaxIn[prt_addRefChannels(ch+1,tdc)]-2;
      // timeL[i] -= 5*(Hits_nFineTime[i]-31)/(max-31);
    }

    if(Hits_nSignalEdge[i]==1){
      if(ch>prt_maxdircch) hh++;
      mult[ch]++;
      if(ch==gTrigger) grTime1 = timeL[i];
      if(Hits_nTdcChannel[i]==0) { //ref channel
	tdcRefTime[tdc] = timeL[i];
	if(prt_getTdcId(gTrigger)==tdc) grTime0 = timeL[i];
      }
    }else{
      timeT[i]=timeL[i];
      grTime2=timeL[i];
    }
  }

  //if(hh<20) return kTRUE;
  // if((grTime0>0 && grTime1>0) || gTrigger==0)
    {
    for(auto i=0; i<Hits_ && i<maxhits; i++){
      
      //if(Hits_nTdcErrCode[i]!=0) continue;

      tdc = map_tdc[Hits_nTrbAddress[i]];
      //if (tdc!=7)continue;
      
      if(tdc<0) continue;
      if(Hits_nSignalEdge[i]==0) continue; // tailing edge
      
      ch = prt_getChannelNumber(tdc,Hits_nTdcChannel[i])-1;      
      hFine[fileid][prt_addRefChannels(ch+1,tdc)]->Fill(Hits_nFineTime[i]);	  
      //if(mult[ch]>1) continue;
      
      if(Hits_nTdcChannel[i]==0) continue; // ref channel
      //      if(ch<prt_maxdircch)
      {
	mcp = map_mcp[ch];
	pix = map_pix[ch];	
	if(gTrigger!=0) {
	  triggerLe = grTime1 - grTime0;
	  triggerTot=grTime2-grTime1;
	  //if(triggerTot<44 || triggerTot>44.05) continue;  // pilas
	  //if(triggerTot<38.7 || triggerTot>38.8) continue; // pico
	  hTimeL[0][16]->Fill(triggerTot);
	}
	hCh->Fill(ch);
	
	if(ch<prt_maxch){
	  timeLe = timeL[i]-tdcRefTime[tdc] - triggerLe;	  
	  //timeLe = tdcRefTime[tdc] - triggerLe;
	  timeTe = timeT[i+1]-tdcRefTime[tdc]-triggerLe;
	  Double_t tot = timeT[i+1]-timeL[i];
	  
	  if(!calib){
	    //if(tot<0 || timeLe<20 || timeLe>40) continue;
	    tot += 30 - gTotO[ch];
	    //timeLe += getTotWalk(tot,ch);
	    //timeLe += getTotWalk(triggerTot,ch,1);
	    ////if(gLeO[ch]) timeLe -=  gLeO[ch]->Eval(tot)-30;
	    //timeLe -= gLeOffArr[ch];
	  }

	  
	  if(mcp< prt_nmcp) {
	    prt_hdigi[mcp]->Fill(map_col[ch],map_row[ch]);
	  
	    TString tdchex = TString::BaseConvert(Form("%d",Hits_nTrbAddress[i]),10,16);
	    hFine[fileid][ch]->SetTitle(Form("tdc 0x%s, chain %d, lch %d = ch %d, m%dp%d ",tdchex.Data() ,(ch/16)%3,ch%16,ch, mcp, pix));
	    hTimeL[mcp][pix]->Fill(timeLe); 
	    hTimeT[mcp][pix]->Fill(timeTe);
	    hTot[fileid][ch]->Fill(tot);
	  
	    hLeTot[ch]->Fill(timeLe,tot);
	    hShape[mcp][pix]->Fill(timeLe,offset);
	    hShape[mcp][pix]->Fill(timeLe+tot,offset);
	  }
	}
      }
    }
  }

  for(Int_t i=0; i<prt_ntdc; i++) tdcRefTime[i]=0;
  
  return kTRUE;
}

TLine *gLine1 = new TLine(0,0,0,1000);
TLine *gLine2 = new TLine(0,0,0,1000);
TString drawHist(Int_t m, Int_t p){
  TString histname="";
  Int_t ch = map_mpc[m][p];
  
  if(gComboId==0){
    TLegend *leg = new TLegend(0.5,0.7,0.9,0.9);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    Int_t num=0;
    ch = prt_addRefChannels(ch,ch/48)+1;
    for(Int_t j=0; j<nfiles; j++){
      if(gGr[j][ch]->GetN()<1){
	gPad->Clear();
	(new TText(0.5,0.1,gGr[j][ch]->GetTitle()))->Draw();
	continue;
      }
      num++;
      if(j==0) gGr[j][ch]->Draw("AL");  
      else gGr[j][ch]->Draw("L");
      leg->AddEntry(gGr[j][ch],gGr[j][ch]->GetTitle(),"l");
    }
    Double_t y1 = 0.9-0.05*num;
    y1 = (y1>0)? y1 : 0; 
    leg->SetY1(y1);
    leg->Draw();
    histname=Form("gCalib_mcp%dpix%d",m,p);
  }
  if(gComboId==1){
    ch = prt_addRefChannels(ch,ch/48)+1;
    for(Int_t j=0; j<nfiles; j++){
      if(j==0) hFine[j][ch]->Draw();
      else hFine[j][ch]->Draw("same");
    }
    histname=hFine[0][ch]->GetName();
  }
  if(gComboId==2){
    hTimeL[m][p]->Draw();
    //hTimeT[m][p]->Draw("same");
    histname=hTimeL[m][p]->GetName();
    if(!ginFile.Contains("trb")){
      Double_t a,b;
      TGraph* gr = new TGraph(); 
      gLeOff->GetPoint(ch,a,b);

      gLine1->SetX1(a+30);
      gLine1->SetX2(a+30);
      gLine1->SetY1(cTime->GetUymin());
      gLine1->SetY2(cTime->GetUymax());
      gLine1->SetLineColor(kRed);
      gLine1->Draw();

      gLine2->SetX1(b+30);
      gLine2->SetX2(b+30);
      gLine2->SetY1(cTime->GetUymin());
      gLine2->SetY2(cTime->GetUymax());
      gLine2->SetLineColor(kBlue);
      gLine2->Draw();
    }
	
  } 
  if(gComboId==3){
    hShape[m][p]->Draw("colz");
    histname=hShape[m][p]->GetName();
  }
  if(gComboId==4){
    hCh->Draw();
    histname=hCh->GetName();
  }
  if(gComboId==5){
    for(Int_t j=0; j<nfiles; j++){
      if(j==0) hTot[j][ch]->Draw();
      else hTot[j][ch]->Draw("same");
    }

    if(ginFile.Contains("trb")){
      Double_t a,b;
      TGraph* gr = new TGraph(); 
      gTotOff->GetPoint(ch,a,b);

      gLine1->SetX1(a);
      gLine1->SetX2(a);
      gLine1->SetY1(cTime->GetUymin());
      gLine1->SetY2(cTime->GetUymax());
      gLine1->SetLineColor(kRed);
      gLine1->Draw();

      gLine2->SetX1(b);
      gLine2->SetX2(b);
      gLine2->SetY1(cTime->GetUymin());
      gLine2->SetY2(cTime->GetUymax());
      gLine2->SetLineColor(kBlue);
      gLine2->Draw();
    }
    histname=hTot[0][ch]->GetName();
    //hTot[0][1953]->Draw();
  }
  if(gComboId==6){
    hLeTot[ch]->Draw("colz");
    histname=hLeTot[ch]->GetName();
    if(gGrDiff[ch]){
      Double_t* xx = gGrDiff[ch]->GetX();
      Double_t* yy = gGrDiff[ch]->GetY();

      TGraph* gr = new TGraph(gGrDiff[ch]->GetN(),yy,xx);
      gr->SetMarkerStyle(7);
      gr->SetMarkerColor(2);
      gr->Draw("P same");
    }
  }
  if(gComboId==7){
    hRefDiff->Draw();
    histname=hRefDiff->GetName();
    // for(Int_t i=0; i<4; i++){
    //   hSigma[i]->SetLineColor(i+1);
    //   hSigma[i]->Draw("same");
    // }
  }

  return histname;
}

Bool_t tlock = false;
void exec3event(Int_t event, Int_t gx, Int_t gy, TObject *selected){
  // if(gComboId==0 || gComboId==1)
  {
    TCanvas *c = (TCanvas *) gTQSender;
    TPad *pad = (TPad *) c->GetSelectedPad();
    if (!pad) return;
    Float_t x = pad->AbsPixeltoX(gx);
    Float_t y = pad->AbsPixeltoY(gy);
    x = pad->PadtoX(x);
    y = pad->PadtoY(y);
    if(event ==1 && tlock) tlock = false;
    else if(event ==1) tlock = true;
    if(tlock) return;

    if (selected->InheritsFrom(TH2::Class())){
      TH2F *hDigi = (TH2F *) selected;
      Int_t binx = hDigi->GetXaxis()->FindBin(x);
      Int_t biny = hDigi->GetYaxis()->FindBin(y);
      TString smcp = selected->GetName();
      smcp = smcp(3,smcp.Sizeof());
      Int_t mcp = smcp.Atoi();
      //      Int_t pix = 8*(binx-1)+biny-1;
      Int_t pix = 8*(biny-1)+binx-1;
  
      cTime->cd();
      drawHist(mcp,pix);
      cTime->Update(); 
    }
  }
}

void MyMainFrame::DoExport(Int_t type){
  if(type==0){
    gROOT->SetBatch(1);
    TCanvas *cExport = new TCanvas("cExport","cExport",0,0,800,400);
    cExport->SetCanvasSize(800,400);
    TString histname="", filedir=ginFile;
    filedir.Remove(filedir.Last('/'));
    prt_savepath = filedir+"/plots";
  
    std::cout<<"Exporting into  "<<prt_savepath <<std::endl;
    prt_writeString(prt_savepath+"/digi.csv", prt_drawDigi("m,p,v\n",geometry));
    Float_t total = (prt_nmcp-1)*(prt_npix-1);
    if(gComboId==0 || gComboId==1 || gComboId==2 || gComboId==3){
      for(Int_t m=0; m<prt_nmcp; m++){
	for(Int_t p=0; p<prt_npix; p++){
	  cExport->cd();
	  histname=drawHist(m,p);

	  cExport->SetName(histname);
	  prt_canvasAdd(cExport);
	  prt_canvasSave(1,1);
	  prt_canvasDel(cExport->GetName());
    
	  gSystem->ProcessEvents();
	}
      }
    }else{
      histname = updatePlot(gComboId,cExport);
      cExport->SetName(histname);
      prt_canvasAdd(cExport);
      prt_canvasSave(1,1);
      prt_canvasDel(cExport->GetName());
    }

    cExport = (TCanvas *) prt_cdigi->DrawClone();
    cExport->SetCanvasSize(800,400);

    cExport->SetName("digi");
    prt_canvasAdd(cExport);
    prt_canvasSave(1,1);
    prt_canvasDel(cExport->GetName());
  
    gROOT->SetBatch(0);
    std::cout<<"Exporting .. Done"<<std::endl;
    return;
  }

  TString name("FT");
  if(type==2) name="TO";
  if(type==3) name="TP";
  if(type==4) name="LO"; 
  TString filedir=ginFile;
  filedir.Remove(filedir.Last('/'));
  TFile efile(filedir+"/calib"+name+ ".root","RECREATE");
   
  if(type==1){
    for(Int_t c=0; c<prt_maxch; c++){
      gGr[0][c]->SetName(Form("%d",c));
      gGr[0][c]->Write();
    }
    gMaxFine->SetName("10000");
    gMaxFine->Write();
  }
  
  if(type==2){
    gTotOff->SetName("10001");
    gTotOff->Write();
  }

  if(type==3){
    gTotPeaks->SetName("10002");
    gTotPeaks->Write();
  }

  if(type==4){
    gLeOff->SetName("10003");
    gLeOff->Write();
    for (Int_t m=0; m <prt_nmcp; m++) {
      for(Int_t p=0; p<prt_npix; p++){
	Int_t ch = map_mpc[m][p];
	if(gGrDiff[ch]){
	  gGrDiff[ch]->SetName(Form("%d",20000+ch));
	  gGrDiff[ch]->Write();
	}
      }
    }
  }
    
  efile.Write();
  efile.Close();
  std::cout<<"Exporting"+name+" .. Done"<<std::endl;
}

TGraph * getGarph(TH1F *hist){
  TGraph * gr = new TGraph();
  Int_t nbins=600;
  TAxis *axis = hist->GetXaxis();
  Double_t step = 600/(Double_t)nbins;
  Double_t total = hist->GetEntries();
  Int_t point=0;
  for(Int_t b=0; b<nbins; b++){
    Double_t xmin  = 0;//b*step;
    Double_t xmax = (b+1)*step;
    int bmin = axis->FindBin(xmin);
    int bmax = axis->FindBin(xmax);
    Double_t integral = hist->Integral(bmin,bmax);
    integral -= hist->GetBinContent(bmin)*(xmin-axis->GetBinLowEdge(bmin))/axis->GetBinWidth(bmin);
    integral -= hist->GetBinContent(bmax)*(axis->GetBinUpEdge(bmax)-xmax)/axis->GetBinWidth(bmax);
    if(integral>=0){
      gr->SetPoint(point,xmax+step/2., 5*(integral/total));
      point++;
    }
  }
  return gr;
}

TSpectrum *spect = new TSpectrum(10);
void Calibrate(){
  std::cout<<"Creating calibration"<<std::endl;
  gMaxFine = new TGraph();
  gTotOff = new TGraph();
  gLeOff = new TGraph();
  gTotPeaks = new TGraph();

  if(gMode==5) for(Int_t i=0; i<9; i++) prt_hdigi[i]->Reset();	
  
  for(Int_t j=0; j<nfiles; j++){
    for(Int_t c=0; c<prt_maxch; c++){
      TString title = Form("%s  %d", hFine[j][c]->GetTitle(), (Int_t)hFine[j][c]->GetEntries());
      if(gMode==3) title = Form("All  %d", (Int_t)hFine[j][c]->GetEntries());
      if(gMode==4 && j==0) title = Form("All beam  %d", (Int_t)hFine[j][c]->GetEntries());
      if(gMode==4 && j==1) title = Form("All pilas  %d", (Int_t)hFine[j][c]->GetEntries());
      hFine[j][c]->SetLineColor(prt_getColorId(j));

      gGr[j][c] = getGarph(hFine[j][c]);
      gGr[j][c]->SetName(Form("gCalib_%d_ch%d",j,c));
      gGr[j][c]->SetTitle(title);
      gGr[j][c]->GetXaxis()->SetTitle("fine bin, [#]");
      gGr[j][c]->GetYaxis()->SetTitle("fine time, [ns]");
      gGr[j][c]->SetLineColor(prt_getColorId(j));

      Int_t firstbin = hFine[j][c]->FindFirstBinAbove(0);
      Int_t lastbin = hFine[j][c]->FindLastBinAbove(0);
      //      std::cout<<c<<"   "<<firstbin << "  "<<lastbin <<std::endl;
      gMaxFine->SetPoint(c,firstbin,lastbin);

      Int_t threshold =  hTot[j][c]->GetMaximum()*0.2;
      firstbin = hTot[j][c]->FindFirstBinAbove(threshold);
      Double_t xle = hTot[j][c]->GetXaxis()->GetBinCenter(firstbin);
      gTotOff->SetPoint(c, xle, hTot[j][c]->GetMean());

      if(gMode==5){
	gGr[j][c]->Fit("pol1","","",50,400);
	Double_t chi=gGr[j][c]->GetFunction("pol1")->GetChisquare();
	Int_t ch =c-1;
	if(ch%49==0) continue;
	ch= 48*(ch/49)+ch%49;
	Int_t m = ch/64;
	if(m<9) prt_hdigi[m]->Fill(map_col[ch],map_row[ch],chi);
      }
    }
  }

  for(Int_t c=0; c<prt_maxdircch; c++){
    if(hTot[0][c]->Integral()>100){
      Int_t nfound = spect->Search(hTot[0][c],3,"",0.1);
      std::cout<<"nfound  "<<nfound <<std::endl;
      Double_t peak(0);
      for(Int_t i=0; i<10; i++){
	if(i<nfound) peak = spect->GetPositionX()[i];
	else peak = 0;
	std::cout<<i<<" xpeaks "<< peak <<std::endl;
	gTotPeaks->SetPoint(c*10+i, 1, peak);
      }
    }
  }

  for(Int_t m=0; m<prt_nmcp; m++){
    for(Int_t p=0; p<prt_npix; p++){
      if(hTimeL[m][p]->Integral()>100){
	double sigma = prt_fit(hTimeL[m][p],0.25,20,2,20).Y();
	hRefDiff->Fill(sigma);
	Int_t ch = map_mpc[m][p];
	hSigma[(ch/48)%4]->Fill(sigma);
      }
    }
  }

  if(true){
    std::cout<<"Getting LE offsets"<<std::endl;
    for (Int_t m=0; m <prt_nmcp; m++) {
      for(Int_t p=0; p<prt_npix; p++){
	Int_t ch = map_mpc[m][p];
	Int_t threshold =  hTimeL[m][p]->GetMaximum()*0.3;      
	Int_t firstbin = hTimeL[m][p]->FindFirstBinAbove(threshold);
	Double_t xmax = hTimeL[m][p]->GetXaxis()->GetBinCenter(hTimeL[m][p]->GetMaximumBin());
	Double_t xle = hTimeL[m][p]->GetXaxis()->GetBinCenter(firstbin);
	gLeOff->SetPoint(ch, xle-30, xmax-30);
      }
    }
  }

  if(false){
    std::cout<<"Getting LE offsets"<<std::endl;
    TH1D* h;
    TH2F* hh;
    for (Int_t m=0; m <prt_nmcp; m++) {
      for(Int_t p=0; p<prt_npix; p++){
	Int_t ch = map_mpc[m][p];
	Double_t mean = prt_fit(hTimeL[m][p],0.2).X();
	hh =(TH2F*) hLeTot[ch]->Clone("hh");
	//hh->RebinY(2);

	gGrDiff[ch] = new TGraph();
	for (int i=0;i<400;i++){
	  Double_t x = hh->GetYaxis()->GetBinCenter(i);
	  h = hh->ProjectionX(Form("bin%d",i+1),i+1,i+2);
	  Double_t vx = prt_fit((TH1F*)h,0.2,200).X();
	  if(vx==0) vx = mean;
	  if(vx-mean>1) vx=mean+1;
	  if(vx-mean<-1) vx=mean-1;
	  gGrDiff[ch]->SetPoint(i,x,vx);
	}
	gGrDiff[ch]->SetPoint(400,-1,mean);

	gGrDiff[ch]->SetName(Form("gCalib_ch%d",ch));
	gGrDiff[ch]->GetXaxis()->SetTitle("fine bin, [#]");
	gGrDiff[ch]->GetYaxis()->SetTitle("fine time, [ns]");
    
      }
    }
  }
  std::cout<<"Calibration .. Done"<<std::endl;
  
}

void TTSelector::Terminate(){
  std::cout<<"terminate start "<<std::endl;
  for (Int_t m=0; m < prt_nmcp; m++) {
    for (Int_t p=0; p < prt_npix; p++) {
      hTimeL[m][p] = dynamic_cast<TH1F *>(TProof::GetOutput(Form("hTimeL_mcp%dpix%d",m,p), fOutput));
      hTimeT[m][p] = dynamic_cast<TH1F *>(TProof::GetOutput(Form("hTimeT_mcp%dpix%d",m,p), fOutput));
      hShape[m][p] = dynamic_cast<TH2F *>(TProof::GetOutput(Form("hShape_mcp%dpix%d",m,p), fOutput));
    }
    prt_hdigi[m] = dynamic_cast<TH2F *>(TProof::GetOutput(Form("mcp%d", m), fOutput));
  }

  for(Int_t j=0; j<nfiles; j++){
    for(Int_t c=0; c<prt_maxch; c++){
      hFine[j][c] = dynamic_cast<TH1F *>(TProof::GetOutput(Form("hFine_%d_ch%d",j,c), fOutput));
      hTot[j][c] = dynamic_cast<TH1F *>(TProof::GetOutput(Form("hTot_%d_ch%d",j,c), fOutput));
    }
  }

  for(Int_t c=0; c<prt_maxch; c++) hLeTot[c] = dynamic_cast<TH2F *>(TProof::GetOutput(Form("hLeTot_ch%d",c), fOutput));
  for(Int_t c=0; c<prt_ntdc; c++) hSigma[c] = dynamic_cast<TH1F *>(TProof::GetOutput(Form("hSigma%d",c), fOutput));
  
  hCh = dynamic_cast<TH1F *>(TProof::GetOutput("hCh", fOutput));
  hRefDiff = dynamic_cast<TH1F *>(TProof::GetOutput("hRefDiff", fOutput));
  Calibrate();
  std::cout<<"terminate done "<<std::endl;
  
}

TString MyMainFrame::updatePlot(Int_t id, TCanvas *cT){
  if(!cT) cT = cTime;
  gComboId = id;
  cTime->cd();
  if(gComboId==7) drawHist(1,1);
  cTime->Modified();
  cTime->Update();
  return "";
}

MyMainFrame::MyMainFrame(const TGWindow *p, UInt_t w, UInt_t h) : TGMainFrame(p, w, h){

  // Create the embedded canvas
  fEcan = new TRootEmbeddedCanvas(0,this,800,350);
  Int_t wid0 = fEcan->GetCanvasWindowId();
  prt_cdigi = new TCanvas("prt_cdigi",10,10,wid0);
  prt_cdigi->SetMargin(0,0,0,0);
  fEcan->AdoptCanvas(prt_cdigi);

  fTime = new TRootEmbeddedCanvas(0,this,800,350);
  Int_t wid1 = fTime->GetCanvasWindowId();
  cTime = new TCanvas("cTime",10,10,wid1);
  fTime->AdoptCanvas(cTime);

  AddFrame(fEcan, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 0, 0, 0, 0));
  AddFrame(fTime, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 0, 0, 0, 0));
   
  TGHorizontalFrame *hframe = new TGHorizontalFrame(this, 400, 40);

  TGLabel *label1 = new TGLabel(hframe, "   Show: ");
  hframe->AddFrame(label1, new TGLayoutHints(kLHintsTop | kLHintsLeft,5, 5, 5, 5));
  fComboMode = new TGComboBox(hframe);
  hframe->AddFrame(fComboMode, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 5, 5));

  fComboMode->AddEntry("Calibration graph", 0);
  fComboMode->AddEntry("Fine time", 1);
  fComboMode->AddEntry("Le time", 2);
  fComboMode->AddEntry("TOT time", 5);
  fComboMode->AddEntry("Le vs.TOT time", 6);
  fComboMode->AddEntry("Signal shape", 3);
  fComboMode->AddEntry("Channels", 4);
  fComboMode->AddEntry("Ref time diff", 7);
  fComboMode->Connect("Selected(Int_t)", "MyMainFrame", this, "updatePlot(Int_t)");

  TGLayoutHints * layout = new TGLayoutHints(kLHintsBottom | kLHintsLeft,5, 5, 5, 5);
  TString names[]={"&Blog", "Fine &Calib","TOT &offset","TOT &peaks","&LE offset"};
  for(Int_t i=0; i<5; i++){
    TGTextButton * btn = new TGTextButton(hframe,names[i]);
    btn->Connect("Clicked()", "MyMainFrame", this, Form("DoExport(=%d)",i));
    hframe->AddFrame(btn, layout);
  }
  
  TGTextButton *exit = new TGTextButton(hframe, "&Exit ");
  exit->Connect("Pressed()", "MyMainFrame", this, "DoExit()");
  hframe->AddFrame(exit, layout);

  AddFrame(hframe, new TGLayoutHints(kLHintsExpandX | kLHintsCenterX, 2, 2, 2, 2));

  // Set a name to the main frame   
  SetWindowName("tdisplay " + ginFile);
  MapSubwindows();

  // Initialize the layout algorithm via Resize()
  Resize(GetDefaultSize());
  fComboMode->Resize(300, 20);
  fComboMode->Select(0);

 
  prt_setRootPalette(1);
  gStyle->SetOptStat(1001111);
  gStyle->SetOptFit();
  
  TChain* ch = new TChain("T");
  ch->Add(ginFile);

  nfiles=0;
  TObjArray *fileElements=ch->GetListOfFiles();
  TIter next(fileElements);
  TChainElement *chEl=0; 
  TString strfiles="";

  while (( chEl=(TChainElement*)next() )) {
    fileList[nfiles]=chEl->GetTitle();
    strfiles += fileList[nfiles] +" ";
    nfiles++;
  }
  
  if(gMode==3) nfiles=1;
  if(gMode==4) nfiles=2;
  
  TString option = Form("%d %d %d %d ",nfiles,gTrigger,gMode,gSetup)+gcFile+" "+strfiles;
  
  std::cout<<"nfiles "<<nfiles <<std::endl;

  Int_t entries = ch->GetEntries();
  std::cout<<"Entries in chain: "<< entries<<std::endl;
  if(gEntries>0) {
    entries=gEntries;
    std::cout<<"Proccesing "<<entries<<" entries"<<std::endl;
  }
  
  TString workers = Form("workers=%d",gWorkers);
  TProof *proof;
  if(gWorkers>1){
    proof = TProof::Open(workers);
    proof->SetProgressDialog(0);
    proof->Load("tdisplay.C+");
    ch->SetProof();
  }
  
  TTSelector *selector = new TTSelector();
  prt_createMap();
  ch->SetCacheSize(10000000);
  ch->AddBranchToCache("*");
  
  ch->Process(selector,option,entries);

  prt_drawDigi("m,p,v\n",geometry,-2,0);
  updatePlot(0); //gComboId

  prt_cdigi->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)", 0, 0,
		 "exec3event(Int_t,Int_t,Int_t,TObject*)");
  MapWindow();
}

MyMainFrame::~MyMainFrame(){
  delete fEcan;
  delete fTime;
}

void tdisplay(TString inFile= "file.hld.root", Int_t trigger=0, Int_t mode=0, Int_t workers = 4, TString cFile= "calib.root",Int_t entries=0){
  //inFile= "data/dirc/scan1/th_1*.hld.root";
  ginFile = inFile;
  gTrigger = trigger;
  gEntries=entries;
  gcFile = (cFile!="")? cFile: "0"; // fine time calibration
  gMode=mode;
  gSetup = 2015;
  gWorkers = workers;

  new MyMainFrame(gClient->GetRoot(), 800, 800);
}
