#include "prttools.C"

Int_t tharr[tdcmax];

Double_t min(1000000),max(0);
TH1F* hLow = new TH1F("hLow",";channel [#];threshold [?]",1000,0,1000);
TH1F* hHigh = new TH1F("hHigh",";channel [#];threshold [?]",1000,0,1000);


void show_thresholds(TString infile1="thresholds.thr", TString infile2="",Int_t range = 0){
  fSavePath = "auto";
  CreateMap();
  initDigi();

  TString s, stdc,sth,schain;
  Int_t ch, tdc, th, chain;
  ifstream in;

  if(infile2 != ""){
    in.open(infile2.Data());
  
    while (1) {
      in >> s >> s >> s >> stdc >> s >> schain >> s >> ch >> s >>sth >> s >> s;
      if (!in.good()) break;
    
      stdc = stdc.Strip(TString::kTrailing,',');
      tdc = TString::BaseConvert(stdc,16,10).Atoi();

      schain = schain.Strip(TString::kLeading,'0');
      schain = schain.Strip(TString::kTrailing,',');
      chain = schain.Atoi();
    
      sth = sth.Strip(TString::kTrailing,',');
      th = TString::BaseConvert(sth,16,10).Atoi();
      if(th>max) max = th;
      if(th<min) min = th;
    
      Int_t tdcSeq = map_tdc[tdc];
      Int_t channel= tdcSeq*48 + chain*16 + ch;
      tharr[channel]=th;
      hLow->Fill(channel,th);
    }
    in.close();
  }

  if(infile1 != ""){
    if(infile2 != ""){
      max=-10000000;
      min=10000000;
    }
    in.open(infile1.Data());  
    while (1) {
      in >> s >> s >> s >> stdc >> s >> schain >> s >> ch >> s >>sth >> s >> s;
      if (!in.good()) break;
    
      stdc = stdc.Strip(TString::kTrailing,',');
      tdc = TString::BaseConvert(stdc,16,10).Atoi();

      schain = schain.Strip(TString::kLeading,'0');
      schain = schain.Strip(TString::kTrailing,',');
      chain = schain.Atoi();
    
      sth = sth.Strip(TString::kTrailing,',');
      th = TString::BaseConvert(sth,16,10).Atoi();
      if(th>max) max = th;
      if(th<min) min = th;
    
      Int_t tdcSeq = map_tdc[tdc];
      std::cout<<"tdcSeq "<< tdc << "  " <<  tdcSeq<<std::endl;
      Int_t channel= tdcSeq*48 + chain*16 + ch;
      Int_t mcp = channel/64;
      Int_t pix = channel%64;
      Int_t col = pix/2 - 8*(pix/16);
      Int_t row = pix%2 + 2*(pix/16);
      fhDigi[mcp]->Fill(col,row,th-tharr[channel]);
      hHigh->Fill(channel,th);
      
      if(infile2 != ""){
	th = th-tharr[channel];
	if(th>max) max = th;
	if(th<min) min = th;
      }
      //std::cout<<"mcp  "<<mcp << " tdcSeq  " <<tdcSeq << " col  "<< col << "  row "<< row<<std::endl;
    }
    in.close();
  }

  if(range==1) gStyle->SetOptLogz();
  if(range>1){
    min = -range;
    max = range;
  }
  drawDigi("m,p,v\n",3,max,min);

  TPaletteAxis *pal = new TPaletteAxis(0.90,0.1,0.94,0.90,fhDigi[0]);   
  cDigi->cd();
  pal->Draw();

  cTh = new TCanvas("cTh","cTh",0,0,800,400);

  hLow->SetLineColor(2);
  hHigh->SetLineColor(4);

  if(infile2 != ""){
    hLow->Draw();
    hHigh->Draw("same");
  }else{
    hHigh->Draw();
    hLow->Draw("same");
  }
  
  // TClass *cl = cDigi->IsA();
  // l = cl->GetMenuList();
  // n = new TClassMenuItem(TClassMenuItem::kPopupUserFunction,cl, "Draw thresholds. Mode 1","mode1",0,"int,int");
  // l->AddFirst(n);
  // n = new TClassMenuItem(TClassMenuItem::kPopupUserFunction,cl, "Draw thresholds. Mode 2","mode2",0,"int,int");
  // l->AddFirst(n);
  

  
  // TClass *cl = gClient->GetRoot()->IsA();
  // l = cl->GetMenuList();
  
  // n = new TClassMenuItem(TClassMenuItem::kPopupUserFunction,cl,
  //  			 "test no 1","poptest1",0,"int,int");
  // l->AddFirst(n);
   
  // TGMenuBar *fMenuBar = new TGMenuBar(gClient->GetRoot(),100,20,kHorizontalFrame);
  // fMenuBar->AddPopup("UUUUZ");

  // fMenuBar = new TGMenuBar(gClient->GetRoot(),100,20,kHorizontalFrame);
  // fMenuBar->AddPopup("&ZZZZZZZZZZZZZZZZZZZ");
  // EnableEntry(fMenuBar);
    // gClient->SetRoot();
}

// void mode1(){
//   cDigi->cd();
//   hLow->SetLineColor(2);
//   hLow->Draw();
//   hHigh->SetLineColor(4);
//   hHigh->Draw("same");
//   cDigi->Modified();
//   cDigi->Update();
// }

// void mode2(){
//   cDigi->cd();
//   drawDigi("m,p,v\n",3,max,min);
// }
