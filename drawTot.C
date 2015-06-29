// run as: root loadprtlibs.C tloop.C

#define TTSelector_cxx
#include "prttools.C"
#include "drawTot.h"

TString m_fileid;
Int_t ttr1 = -50, ttr2 = -35;

TH1F *hTimeL,*hTimeT,*hTot,*hCh;
TH2F *hLeTot,*hTeTot,*hLeTe;

void TTSelector::Begin(TTree *){
  fSavePath = "auto";

  hTimeL = new TH1F("hTime",m_fileid+";LE [ns];counts [#]",500,-135,-74);
  hTimeT = new TH1F("hTime",m_fileid+";TE [ns];counts [#]",500,-135,-74);
  hTot = new TH1F("hTot",m_fileid+";TOT [ns];counts [#]",200,ttr1,ttr2);
  hCh = new TH1F("hCh",m_fileid+";channel [#];counts [#]",1500,0,1500);

  // TH2F * hLeTot = new TH2F("hLeTot","hLeTot;LE [ns];TOT [ns]",200,88,94,200,-50,-40);
  // TH2F * hTeTot = new TH2F("hTeTot","hTeTot;TE [ns];TOT [ns]",200,131,140,200,-50,-40);
  // TH2F * hLeTe  = new TH2F("hLeTe","hLeTe;LE [ns];TE [ns]",200,88,94,200,131,140);

  hLeTot = new TH2F("hLeTot",m_fileid+";LE [ns];TOT [ns]",200,-135,-120,200,ttr1,ttr2);
  hTeTot = new TH2F("hTeTot",m_fileid+";TE [ns];TOT [ns]",200,-95,-75,200,ttr1,ttr2);
  hLeTe  = new TH2F("hLeTe",m_fileid+";LE [ns];TE [ns]",200,-135,-120,200,-95,-75);
 
  CreateMap();
}

Bool_t TTSelector::Process(Long64_t entry){
  if(entry%1000==0) std::cout<<"event # "<< entry <<std::endl;
  GetEntry(entry);
  
  Double_t reftimes[tdcnum]={0};
  Double_t trailtimes[maxch]={0};
  Double_t trigTime0(0),trigTime1(0);
  Int_t prech(0);
  
  for(Int_t i=0; i<Hits_; i++){
    Int_t tdc = map_tdc[Hits_nTrbAddress[i]];
    Int_t ch = GetChannelNumber(tdc,Hits_nTdcChannel[i]);
    if(badcannel(ch)) continue;
    
    if(Hits_nTdcChannel[i]==0) reftimes[tdc]=Hits_fTime[i];
    if(Hits_nSignalEdge[i]==0){
      trailtimes[prech]=Hits_fTime[i];
      continue;
    }
    prech = ch;
    
    if(ch ==1345) trigTime1 = Hits_fTime[i];
    if(ch ==1344) trigTime0 = Hits_fTime[i];
  }

  
  for(Int_t i=0; i<Hits_; i++){
    Int_t tdc = map_tdc[Hits_nTrbAddress[i]];
    Int_t ch = GetChannelNumber(tdc,Hits_nTdcChannel[i]);
    if(badcannel(ch)) continue;
    if(Hits_nSignalEdge[i]==0 || Hits_nTdcChannel[i]==0 ) continue;

    hCh->Fill(ch);	
    if(ch == 39){
      Double_t le = Hits_fTime[i]-reftimes[tdc]-(trigTime1-trigTime0);
      Double_t te = trailtimes[ch]-reftimes[tdc]-(trigTime1-trigTime0);
      //      if(le>90.5 && le<91.5)
	{
	  hTimeL->Fill(le);
	  hTimeT->Fill(te);
	  hTot->Fill(le-te);
	  
	  hLeTot->Fill(le,le-te);
	  hTeTot->Fill(te,le-te);
	  hLeTe->Fill(le,te);	  
	}
      //std::cout<<"LE   "<< le << " TE "<< te <<std::endl;
    }
  }
  
  return kTRUE;
}
 
void TTSelector::Terminate(){
  canvasAdd("LE",800,500);
  hTimeL->Draw();
  canvasAdd("TE",800,500);
  hTimeT->Draw();
  canvasAdd("TOT",800,500);
  hTot->Draw();
  
  canvasAdd("LeTOT",800,500);
  hLeTot->Draw("colz");
  canvasAdd("TeTOT",800,500);
  hTeTot->Draw("colz");
  canvasAdd("LeTe",800,500);
  hLeTe->Draw("colz");
  //hCh->Draw();
  canvasSave(0,0);
} 


void drawTot(TString inFile= "../data/pilas_15178162456.hld.root_calibrated.root", Int_t events = 100000){
  m_fileid =  inFile;
  m_fileid =  m_fileid.Remove(0,m_fileid.Last('/')+1);
  gStyle->SetOptStat(1001111);
  TChain* ch = new TChain("T");
  ch->Add(inFile);
  
  Int_t entries = ch->GetEntries();
  std::cout<<"Entries in chain:  "<< entries<<std::endl;
  if(events!=0) entries = events;
  
  TTSelector *selector = new TTSelector();
  ch->Process(selector,"",entries);
}
