// tdisplay - tool to plot different quantities from the .hld file
// original author: Roman Dzhygadlo - GSI Darmstadt 

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

#include "tdisplay.h"

Int_t nfiles = 10;
const Int_t maxfiles = 200;
const Int_t maxch =3000;
const Int_t nmcp = 15, npix = 64;
TString fileList[maxfiles];

TH1F *hFine[maxfiles][maxch];
TH1F *hTot[maxfiles][maxch];
TH2F *hLeTot[maxch];

TH1F *hTimeL[nmcp][npix];
TH1F *hTimeT[nmcp][npix];
TH2F *hShape[nmcp][npix];

TH1F *hCh;

TString ginFile="";
Int_t gTrigger, gMode=0;;

const Int_t tdcnum=88;
const Int_t tdcmax=10000;
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
TGraph *gGr[maxfiles][maxch];

TCanvas *cTime;

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
    Int_t col = pix/2 - 8*(pix/16);
    Int_t row = pix%2 + 2*(pix/16);
    pix = col*8+row;
    chmap[mcp][pix]=ch;
  }
}

void TTSelector::SlaveBegin(TTree *){
  TString option = GetOption();
  TObjArray *strobj = option.Tokenize(" ");
  nfiles = ((TObjString*)strobj->At(0))->GetString().Atoi();
  gTrigger = ((TObjString*)strobj->At(1))->GetString().Atoi();
  gMode = ((TObjString*)strobj->At(2))->GetString().Atoi();
  for(Int_t i=3; i<nfiles+2; i++){
    fileList[i-3]=((TObjString*)strobj->At(i))->GetString();
    std::cout<<" fileList[i]  "<<fileList[i-3] <<std::endl;
  }

  for(Int_t j=0; j<nfiles; j++){
    for(Int_t c=0; c<maxch; c++){
      hFine[j][c] = new TH1F(Form("hFine_%d_ch%d",j,c),Form("hFine_%d_ch%d",j,c) , 600,1,600);
      hTot[j][c] = new TH1F(Form("hTot_%d_ch%d",j,c), Form("hTot_%d_ch%d",j,c) , 500,-2,6);
      hTot[j][c]->SetLineColor(j+1);
      fOutput->Add(hFine[j][c]);
      fOutput->Add(hTot[j][c]);
    }
  }

  const Int_t lb = 85, hb = 105;
  for(Int_t c=0; c<maxch; c++){
    hLeTot[c] = new TH2F(Form("hLeTot_ch%d",c), Form("hLeTot_ch%d",c) ,200,lb,hb, 100,-2,5);
    fOutput->Add(hLeTot[c]);
  }

  for(Int_t m=0; m<nmcp; m++){
    for(Int_t p=0; p<npix; p++){
     
      hTimeL[m][p] = new TH1F(Form("hTimeL_mcp%dpix%d",m,p), Form("hTimeL_%d_%d",m,p) , 500,lb,hb);
      hTimeT[m][p] = new TH1F(Form("hTimeT_mcp%dpix%d",m,p), Form("hTimeT_%d_%d",m,p) , 500,lb,hb);
      hShape[m][p] = new TH2F(Form("hShape_mcp%dpix%d",m,p), Form("hShape_%d_%d",m,p) , 400,lb,hb,130,-5,35);

      fOutput->Add(hTimeL[m][p]);
      fOutput->Add(hTimeT[m][p]);
      fOutput->Add(hShape[m][p]);
    }

    fhDigi[m] = new TH2F( Form("mcp%d", m),Form("mcp%d", m),8,0.,8.,8,0.,8.);
    fhDigi[m]->SetStats(0);
    fhDigi[m]->SetTitle(0);
    fhDigi[m]->GetXaxis()->SetNdivisions(10);
    fhDigi[m]->GetYaxis()->SetNdivisions(10);
    fhDigi[m]->GetXaxis()->SetLabelOffset(100);
    fhDigi[m]->GetYaxis()->SetLabelOffset(100);
    fhDigi[m]->GetXaxis()->SetTickLength(1);
    fhDigi[m]->GetYaxis()->SetTickLength(1);
    fhDigi[m]->GetXaxis()->SetAxisColor(15);
    fhDigi[m]->GetYaxis()->SetAxisColor(15);
    fOutput->Add(fhDigi[m]);
  }

  hCh = new TH1F("hCh","hCh",3000,0,3000); 
  fOutput->Add(hCh);
  CreateMap();
}

Bool_t TTSelector::Process(Long64_t entry){
  Int_t trbSeqId,ch,mcp,pix,col,row;
  Double_t grTime0(0), grTime1(0),timeLe(0), timeTe(0), offset(0);

  TString current_file_name  = TTSelector::fChain->GetCurrentFile()->GetName();
  TObjArray *sarr = current_file_name.Tokenize("_");
  if(sarr->At(1)){
    TString soffset = ((TObjString *) sarr->At(1))->GetName();
    offset = soffset.Atof();
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
    trbSeqId = tdcmap[Hits_nTrbAddress[i]];
    ch = 32*trbSeqId+Hits_nTdcChannel[i];
    if(++mult[ch]>50) continue;
    timeTe0[ch][mult[ch]]=Hits_fTime[i];
    if(Hits_nTdcChannel[i]==0) { //ref channel
      trbRefTime[trbSeqId] = Hits_fTime[i];
      if((gTrigger-ch)<=32 && (gTrigger-ch)>0) grTime0 = Hits_fTime[i];
    }
    if(ch==gTrigger) grTime1 = Hits_fTime[i];
  }

  if((grTime0>0 && grTime1>0) || gTrigger==0){
    for(Int_t i=0; i<Hits_; i++){
      trbSeqId = tdcmap[Hits_nTrbAddress[i]];
      ch = 32*trbSeqId+Hits_nTdcChannel[i];
   
      hFine[fileid][ch]->Fill(Hits_nFineTime[i]);

      if(ch%2==0 || Hits_nTrbAddress[i]==0) continue; // go away trailing edge and ref channel
      hCh->Fill(ch);
      if(ch<3000) {
	mcp = ch/128;
	pix = (ch%128)/2;	
	col = pix/2 - 8*(pix/16);
	row = pix%2 + 2*(pix/16);
	pix = col+8*row;

	if(gMode!=3){
	  // noisy pixels
	  if(ch==1397) continue;
	  if(mcp==2  && pix==55) continue;
	  if(mcp==2  && pix==62) continue;
	  if(mcp==14 && pix==35) continue;
	}
       
	timeLe = Hits_fTime[i]-trbRefTime[trbSeqId];
	timeTe = timeTe0[ch][0]-trbRefTime[trbSeqId];
	Double_t triggerTime = grTime1-grTime0;
	if(mcp<15){
	  fhDigi[mcp]->Fill(col,row);
	  	 
	  hFine[fileid][ch]->SetTitle(Form("ch %d m%dp%d ",ch, mcp, pix));
	  hTimeL[mcp][pix]->Fill(timeLe - triggerTime); 
	  hTimeT[mcp][pix]->Fill(timeTe - triggerTime);
	  hShape[mcp][pix]->Fill(timeLe - triggerTime,offset);
	  hShape[mcp][pix]->Fill(timeTe - triggerTime,offset);
	}
	//if(ch==gTrigger) 
	hTot[fileid][ch]->Fill(timeTe0[ch+1][0] - timeTe0[ch][0]);
	hLeTot[ch]->Fill(timeLe - triggerTime,timeTe0[ch+1][0] - timeTe0[ch][0]);
      }
    }
  }

  for(Int_t i=0; i<Hits_; i++){
    trbSeqId = tdcmap[Hits_nTrbAddress[i]];
    ch = 32*trbSeqId+Hits_nTdcChannel[i];
    mult[ch]=-1;
    for(Int_t j=0; j<50; j++){
      timeTe0[ch][j]=0; 
    }
  }
  return kTRUE;
}

TString drawHist(Int_t m, Int_t p){
  TString histname="";
  Int_t ch = chmap[m][p];
  
  if(gComboId==0){
    TLegend *leg = new TLegend(0.5,0.7,0.9,0.9);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    Int_t num=0;
    for(Int_t j=0; j<nfiles; j++){
      if(gGr[j][ch]->GetN()<1) continue;
      num++;
      if(j==0) gGr[j][ch]->Draw("AL");  
      else gGr[j][ch]->Draw("L");
      leg->AddEntry(gGr[j][ch],gGr[j][ch]->GetTitle(),"l");
    }
    Double_t y1 = 0.9-0.05*num;
    y1 = (y1>0)? y1 : 0; 
    leg->SetY1(y1);
    leg->Draw();
    histname=gGr[0][ch]->GetName();
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
    hTimeT[m][p]->Draw("same");
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
      Int_t pix = 8*(binx-1)+biny-1;
  
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
  Int_t saveFlag = 1;
  TString histname="", filedir=ginFile;
  filedir.Remove(filedir.Last('/'));
  TString path = createDir(filedir+"/plots", "beam_data "+ginFile, saveFlag); 
  std::cout<<"Exporting into  "<<path <<std::endl;
  writeInfo("digi.csv", drawDigi("m,p,v\n",1), saveFlag);
  //pbar->Reset();
  Float_t total = (nmcp-1)*(npix-1);
  if(gComboId==0 || gComboId==1 || gComboId==2 || gComboId==3){
    for(Int_t m=0; m<nmcp; m++){
      for(Int_t p=0; p<npix; p++){
	cExport->cd();
	histname=drawHist(m,p);
	save(cExport,path,histname,"tdisplay",saveFlag,1);
	//pbar->SetPosition(100*(m*p)/total);
	gSystem->ProcessEvents();
      }
    }
  }else{
    histname = updatePlot(gComboId,cExport);
    save(cExport,path,histname,"tdisplay",saveFlag,1);
  }
  cExport = (TCanvas *) cDigi->DrawClone();
  cExport->SetCanvasSize(800,400);
  save(cExport,path,"digi","tdisplay",saveFlag,1);
  gROOT->SetBatch(0);
  std::cout<<"Exporting .. Done"<<std::endl;
}

void MyMainFrame::DoExportGr(){
  TString filedir=ginFile;
  filedir.Remove(filedir.Last('/'));
  TFile efile(filedir+"/calib.root","RECREATE");
  
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
  Calibrate();
}

TString MyMainFrame::updatePlot(Int_t id, TCanvas *cT){
  if(!cT) cT = cTime;
  gComboId = id;
  cTime->cd();
  cTime->Modified();
  cTime->Update();
  return "";
}

void MyMainFrame::DoExit(){
  gApplication->Terminate(0);
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
  
  TString option = Form("%d %d %d ",nfiles,gTrigger,gMode)+strfiles;
  
  std::cout<<"nfiles "<<nfiles <<std::endl;

  Int_t entries = ch->GetEntries();
  std::cout<<"Entries in chain:  "<< entries<<std::endl;
 
  TString workers = "workers=4";
  if(gSystem->GetFromPipe("whoami")=="hadaq" && entries>1000000) workers = "workers=12";

  TProof *proof = TProof::Open(workers);
  proof->SetProgressDialog(0);
  proof->Load("tdisplay.C+");
  ch->SetProof();
  TTSelector *selector = new TTSelector();
  CreateMap();
  ch->Process(selector,option,entries);
    
  drawDigi("m,p,v\n",1);
  updatePlot(0); //gComboId

  cDigi->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)", 0, 0,
		 "exec3event(Int_t,Int_t,Int_t,TObject*)");
  MapWindow();
}

MyMainFrame::~MyMainFrame(){
  Cleanup();
  delete fEcan;
  delete fTime;
}

void tdisplay(TString inFile= "file.hld.root", Int_t trigger=1920, Int_t mode =0){ //1952
  //inFile= "data/dirc/scan1/th_1*.hld.root";
  ginFile = inFile;
  gTrigger = trigger;
  gMode=mode;
  new MyMainFrame(gClient->GetRoot(), 800, 800);
}
