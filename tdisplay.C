// tdisplay - tool to plot different quantities from the .hld file
// original author: Roman Dzhygadlo - GSI Darmstadt 

#define TTSelector_cxx
#include "prttools.C"
#include "tdisplay.h"

Int_t gSetup=2015, gTrigger, gMode=0, gComboId=0, gWorkers=4, nfiles = 10;
TString ginFile="";
const Int_t maxfiles = 150;
TString fileList[maxfiles];

TH1F *hCh, *hRefDiff;
TH1F *hFine[maxfiles][maxch];
TH1F *hTot[maxfiles][maxch];
TH1F *hTimeL[nmcp][npix];
TH1F *hTimeT[nmcp][npix];

TH2F *hLeTot[maxch];
TH2F *hShape[nmcp][npix];
TGraph *gGr[maxfiles][maxch];

Double_t tdcRefTime[100];
Double_t timeTe0[tdcmax][50];
Int_t mult[tdcmax];
TCanvas *cTime;

void TTSelector::SlaveBegin(TTree *){
  std::cout<<"init starts "<<std::endl;
  TString option = GetOption();
  TObjArray *strobj = option.Tokenize(" ");
  nfiles = ((TObjString*)strobj->At(0))->GetString().Atoi();
  gTrigger = ((TObjString*)strobj->At(1))->GetString().Atoi();
  gMode = ((TObjString*)strobj->At(2))->GetString().Atoi();
  gSetup = ((TObjString*)strobj->At(3))->GetString().Atoi();
  for(Int_t i=4; i<nfiles+4; i++){
    fileList[i-4]=((TObjString*)strobj->At(i))->GetString();
    std::cout<<" fileList[i]  "<<fileList[i-4] <<std::endl;
  } 
  
  for(Int_t j=0; j<nfiles; j++){
    for(Int_t c=0; c<maxch; c++){
      hFine[j][c] = new TH1F(Form("hFine_%d_ch%d",j,c),Form("hFine_%d_ch%d;bin [#];LE entries [#]",j,c) , 600,1,600);
      hTot[j][c] = new TH1F(Form("hTot_%d_ch%d",j,c), Form("hTot_%d_ch%d;TOT [ns];entries [#]",j,c) , 500,0,70);
      hTot[j][c]->SetLineColor(j+1);
      fOutput->Add(hFine[j][c]);
      fOutput->Add(hTot[j][c]);
    }
  }
  
  const Long_t lb = 0, hb = 100;
  for(Int_t c=0; c<maxch; c++){
    hLeTot[c] = new TH2F(Form("hLeTot_ch%d",c), Form("hLeTot_ch%d",c) ,200,lb,hb, 100,-2,5);
    fOutput->Add(hLeTot[c]);
  }

  initDigi(0);
  for(Int_t m=0; m<nmcp; m++){
    for(Int_t p=0; p<npix; p++){

      hTimeL[m][p] = new TH1F(Form("hTimeL_mcp%dpix%d",m,p), Form("hTimeL_%d_%d;LE time [ns];entries [#]",m,p) , 500,lb,hb);
      hTimeT[m][p] = new TH1F(Form("hTimeT_mcp%dpix%d",m,p), Form("hTimeT_%d_%d;TOT [ns];entries [#]",m,p) , 500,lb,hb);
      hShape[m][p] = new TH2F(Form("hShape_mcp%dpix%d",m,p), Form("hShape_%d_%d;time [ns];offset to the threshold [mV]",m,p) , 400,lb,hb,130,-5,35);

      fOutput->Add(hTimeL[m][p]);
      fOutput->Add(hTimeT[m][p]);
      fOutput->Add(hShape[m][p]);
    }
    fOutput->Add(fhDigi[m]);
  }

  hCh = new TH1F("hCh","hCh;channel [#];entries [#]",3000,0,3000);
  fOutput->Add(hCh);
  hRefDiff = new TH1F("hRefDiff","hRefDiff;time [ns];entries [#]",500,-10,10); 
  fOutput->Add(hRefDiff);
  CreateMap();
  std::cout<<"init done " <<std::endl;
  
}

Bool_t TTSelector::Process(Long64_t entry){
  if(entry%1000==0) std::cout<<"event # "<< entry <<std::endl;
  Int_t tdcSeqId,ch,mcp,pix,col,row;
  Double_t triggerTime(0), grTime0(0), grTime1(0),timeLe(0), timeTe(0), offset(0);
  Double_t timeT[50000];

 
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

  GetEntry(entry);
  for(Int_t i=0; i<Hits_; i++){
    tdcSeqId = map_tdc[Hits_nTrbAddress[i]];
    if(tdcSeqId<0) continue;
    ch = ctdc*tdcSeqId+Hits_nTdcChannel[i];
    //if(Hits_nTdcErrCode[i]!=0) continue

    if(ch==gTrigger) grTime1 = Hits_fTime[i];
    if(Hits_nTdcChannel[i]==0) { //ref channel
      tdcRefTime[tdcSeqId] = Hits_fTime[i];
      if(gTrigger/48==tdcSeqId) grTime0 = Hits_fTime[i];
      continue;
    }

    if(Hits_nSignalEdge[i]==0){
      timeT[ch]=Hits_fTime[i];
      continue;
    }
  }
 
  if((grTime0>0 && grTime1>0) || gTrigger==0){
    for(Int_t i=0; i<Hits_; i++){
      if(Hits_nSignalEdge[i]==0) continue; //tailing edge
      //if(Hits_nTdcErrCode[i]!=0) continue;
      
      tdcSeqId = map_tdc[Hits_nTrbAddress[i]];
      if(tdcSeqId<0) continue;
      ch = ctdc*tdcSeqId+Hits_nTdcChannel[i];
      hFine[fileid][AddRefChannels(ch)]->Fill(Hits_nFineTime[i]);
      if(Hits_nTdcChannel[i]==0) continue; // ref channel
      ch -=1;
      if(ch<3000) {
	mcp = map_mcp[ch];
	pix = map_pix[ch];	
	triggerTime = grTime1 - grTime0;

	hCh->Fill(ch);
	timeLe = Hits_fTime[i]-tdcRefTime[tdcSeqId] - triggerTime;
	timeTe = timeT[ch]-tdcRefTime[tdcSeqId]-triggerTime;
	  
	if(ch<maxmch){
	  fhDigi[mcp]->Fill(map_col[ch],map_row[ch]);
	  TString tdchex = TString::BaseConvert(Form("%d",Hits_nTrbAddress[i]),10,16);
	  hFine[fileid][ch]->SetTitle(Form("tdc 0x%s, chain %d, lch %d = ch %d, m%dp%d ",tdchex.Data() ,(ch/16)%3,ch%16,ch, mcp, pix));
	  hTimeL[mcp][pix]->Fill(timeLe); 
	  hTimeT[mcp][pix]->Fill(timeTe);
	  hShape[mcp][pix]->Fill(timeLe,offset);
	  hShape[mcp][pix]->Fill(timeTe,offset);
	}
	hTot[fileid][ch]->Fill(timeTe-timeLe);
	hLeTot[ch]->Fill(timeLe,timeTe-timeLe);
      }
    }
  }

  for(Int_t i=0; i<tdcnum; i++){
    for(Int_t j=i+1; j<tdcnum; j++){
      hRefDiff->Fill(tdcRefTime[j]-tdcRefTime[i]);
    }
  }

  for(Int_t i=0; i<10; i++){
    tdcRefTime[i]=0;
  }
  
  return kTRUE;
}

TString drawHist(Int_t m, Int_t p){
  TString histname="";
  Int_t ch = map_mpc[m][p];
  ch = AddRefChannels(ch)+1;
  
  if(gComboId==0){
    TLegend *leg = new TLegend(0.5,0.7,0.9,0.9);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    Int_t num=0;
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
    histname=hTot[0][ch]->GetName();
    //hTot[0][1953]->Draw();
  }
  if(gComboId==6){
    hLeTot[ch]->Draw("colz");
    histname=hLeTot[ch]->GetName();
  }

  if(gComboId==7){
    hRefDiff->Draw();
    histname=hRefDiff->GetName();
  }

  return histname;
}

Bool_t lock = false;
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
    if(event ==1 && lock) lock = false;
    else if(event ==1) lock = true;
    if(lock) return;

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

void MyMainFrame::DoExport(){
  gROOT->SetBatch(1);
  TCanvas *cExport = new TCanvas("cExport","cExport",0,0,800,400);
  cExport->SetCanvasSize(800,400);
  TString histname="", filedir=ginFile;
  filedir.Remove(filedir.Last('/'));
  fSavePath = filedir+"/plots";
  
  std::cout<<"Exporting into  "<<fSavePath <<std::endl;
  writeString(fSavePath+"/digi.csv", drawDigi("m,p,v\n",2));
  Float_t total = (nmcp-1)*(npix-1);
  if(gComboId==0 || gComboId==1 || gComboId==2 || gComboId==3){
    for(Int_t m=0; m<nmcp; m++){
      for(Int_t p=0; p<npix; p++){
	cExport->cd();
	histname=drawHist(m,p);

	cExport->SetName(histname);
	canvasAdd(cExport);
	canvasSave(1,1);
	canvasDel(cExport->GetName());
    
	gSystem->ProcessEvents();
      }
    }
  }else{
    histname = updatePlot(gComboId,cExport);
    cExport->SetName(histname);
    canvasAdd(cExport);
    canvasSave(1,1);
    canvasDel(cExport->GetName());
  }

  cExport = (TCanvas *) cDigi->DrawClone();
  cExport->SetCanvasSize(800,400);

  cExport->SetName("digi");
  canvasAdd(cExport);
  canvasSave(1,1);
  canvasDel(cExport->GetName());
  
  gROOT->SetBatch(0);
  std::cout<<"Exporting .. Done"<<std::endl;
}

void MyMainFrame::DoExportGr(){
  TString filedir=ginFile;
  filedir.Remove(filedir.Last('/'));
  TFile efile(filedir+"/calib1.root","RECREATE");
  
  for(Int_t c=0; c<maxch; c++){
    gGr[0][c]->SetName(Form("%d",c));
    gGr[0][c]->Write();
  }
  efile.Write();
  efile.Close();
  std::cout<<"Exporting .. Done"<<std::endl;
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
    if(integral>0){
      gr->SetPoint(point,xmax+step/2., 5*(integral/total));
      point++;
    }
  }
  return gr;
}

void Calibrate(){
  std::cout<<"Creating calibration"<<std::endl;

  for(Int_t j=0; j<nfiles; j++){
    for(Int_t c=0; c<maxch; c++){
      TString title = Form("%s  %d", hFine[j][c]->GetTitle(), (Int_t)hFine[j][c]->GetEntries());
      if(gMode==3) title = Form("All  %d", (Int_t)hFine[j][c]->GetEntries());
      if(gMode==4 && j==0) title = Form("All beam  %d", (Int_t)hFine[j][c]->GetEntries());
      if(gMode==4 && j==1) title = Form("All pilas  %d", (Int_t)hFine[j][c]->GetEntries());
      hFine[j][c]->SetLineColor(getColorId(j));

      gGr[j][c] = getGarph(hFine[j][c]);
      gGr[j][c]->SetName(Form("gCalib_%d_ch%d",j,c));
      gGr[j][c]->SetTitle(title);
      gGr[j][c]->GetXaxis()->SetTitle("fine bin, [#]");
      gGr[j][c]->GetYaxis()->SetTitle("fine time, [ns]");
      gGr[j][c]->SetLineColor(getColorId(j));
    }
  }
}

void TTSelector::Terminate(){
  std::cout<<"terminate start "<<std::endl;
  for (Int_t m=0; m < nmcp; m++) {
    for (Int_t p=0; p < npix; p++) {
      hTimeL[m][p] = dynamic_cast<TH1F *>(TProof::GetOutput(Form("hTimeL_mcp%dpix%d",m,p), fOutput));
      hTimeT[m][p] = dynamic_cast<TH1F *>(TProof::GetOutput(Form("hTimeT_mcp%dpix%d",m,p), fOutput));
      hShape[m][p] = dynamic_cast<TH2F *>(TProof::GetOutput(Form("hShape_mcp%dpix%d",m,p), fOutput));
    }
    fhDigi[m] = dynamic_cast<TH2F *>(TProof::GetOutput(Form("mcp%d", m), fOutput));
  }

  for(Int_t j=0; j<nfiles; j++){
    for(Int_t c=0; c<maxch; c++){
      hFine[j][c] = dynamic_cast<TH1F *>(TProof::GetOutput(Form("hFine_%d_ch%d",j,c), fOutput));
      hTot[j][c] = dynamic_cast<TH1F *>(TProof::GetOutput(Form("hTot_%d_ch%d",j,c), fOutput));
    }
  }

  for(Int_t c=0; c<maxch; c++){
    hLeTot[c] = dynamic_cast<TH2F *>(TProof::GetOutput(Form("hLeTot_ch%d",c), fOutput));
  }

  hCh = dynamic_cast<TH1F *>(TProof::GetOutput("hCh", fOutput));
  hRefDiff = dynamic_cast<TH1F *>(TProof::GetOutput("hRefDiff", fOutput));
  Calibrate();
  std::cout<<"terminate done "<<std::endl;
  
}

TString MyMainFrame::updatePlot(Int_t id, TCanvas *cT){
  if(!cT) cT = cTime;
  gComboId = id;
  cTime->cd();
  cTime->Modified();
  cTime->Update();
  return "";
}

MyMainFrame::MyMainFrame(const TGWindow *p, UInt_t w, UInt_t h) : TGMainFrame(p, w, h){

  // Create the embedded canvas
  fEcan = new TRootEmbeddedCanvas(0,this,800,350);
  Int_t wid0 = fEcan->GetCanvasWindowId();
  cDigi = new TCanvas("cDigi",10,10,wid0);
  cDigi->SetMargin(0,0,0,0);
  fEcan->AdoptCanvas(cDigi);

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

  TGTextButton * fBtnExport = new TGTextButton(hframe, "Export for the &blog");
  fBtnExport->Connect("Clicked()", "MyMainFrame", this, "DoExport()");
  hframe->AddFrame(fBtnExport, new TGLayoutHints(kLHintsBottom | kLHintsLeft,5, 5, 5, 5));

  TGTextButton * fBtnExportGr = new TGTextButton(hframe, "Export curves");
  fBtnExportGr->Connect("Clicked()", "MyMainFrame", this, "DoExportGr()");
  hframe->AddFrame(fBtnExportGr, new TGLayoutHints(kLHintsBottom | kLHintsLeft,5, 5, 5, 5));


  TGTextButton *exit = new TGTextButton(hframe, "&Exit ");
  exit->Connect("Pressed()", "MyMainFrame", this, "DoExit()");
  hframe->AddFrame(exit, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));

  AddFrame(hframe, new TGLayoutHints(kLHintsExpandX | kLHintsCenterX, 2, 2, 2, 2));

  // Set a name to the main frame   
  SetWindowName("tdisplay " + ginFile);
  MapSubwindows();

  // Initialize the layout algorithm via Resize()
  Resize(GetDefaultSize());
  fComboMode->Resize(300, 20);
  fComboMode->Select(0);

 
  SetRootPalette(1);
  gStyle->SetOptStat(1001111);

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
  
  TString option = Form("%d %d %d %d ",nfiles,gTrigger,gMode,gSetup)+strfiles;
  
  std::cout<<"nfiles "<<nfiles <<std::endl;

  Int_t entries = ch->GetEntries();
  std::cout<<"Entries in chain:  "<< entries<<std::endl;
 
  TString workers = Form("workers=%d",gWorkers);
  TProof *proof;
  if(gWorkers>1){
    proof = TProof::Open(workers);
    proof->SetProgressDialog(0);
    proof->Load("tdisplay.C+");
    ch->SetProof();
  }
  
  TTSelector *selector = new TTSelector();
  CreateMap();
  ch->SetCacheSize(10000000);
  ch->AddBranchToCache("*");
  
  ch->Process(selector,option,entries);

  drawDigi("m,p,v\n",2);  
  updatePlot(0); //gComboId

  cDigi->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)", 0, 0,
		 "exec3event(Int_t,Int_t,Int_t,TObject*)");
  MapWindow();
}

MyMainFrame::~MyMainFrame(){
  delete fEcan;
  delete fTime;
}

void tdisplay(TString inFile= "file.hld.root", Int_t trigger=0, Int_t mode=0, Int_t workers = 4){
  //inFile= "data/dirc/scan1/th_1*.hld.root";
  ginFile = inFile;
  gTrigger = trigger;
  gMode=mode;
  gSetup = 2015;
  gWorkers = workers;

  new MyMainFrame(gClient->GetRoot(), 800, 800);
}
