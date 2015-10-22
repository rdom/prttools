#define prt__beam
#include "../prtdirc/src/PrtHit.h"
#include "../prtdirc/src/PrtEvent.h"
#include "prttools.C"

void procOffsets(TString path="/data.local/data/jun15/beam_15180221900C.root"){
  
  if(path=="") return;
 
  TString insim = path; insim.ReplaceAll("C.root","S.root");
  PrtInit(path,1);
  fCh->Add(insim);
  fNEntries = fCh->GetEntries();

  TH1F * hLeD  = new TH1F("leD","LE; LE [ns]; entries [#]",5000,200,400);
  TH1F * hLeS  = new TH1F("leS","LE; LE [ns]; entries [#]",5000,0,100);

  PrtHit fHit;
  Int_t entries = fCh->GetEntries();
  for (Int_t ievent=0; ievent<entries; ievent++){
    PrtNextEvent(ievent,1000);
    bool bsim(false);
    TString current_file_name  = fCh->GetCurrentFile()->GetName();
    if(current_file_name.Contains("S.root")) bsim = true;

    Double_t time(0);
    if(fEvent->GetParticle()!=2212) continue;

    for(Int_t i=0; i<fEvent->GetHitSize(); i++){
      fHit = fEvent->GetHit(i);

      if(fHit.GetChannel()<960 ){
	if(!bsim) hLeD->Fill(fHit.GetLeadTime());
	else hLeS->Fill(fHit.GetLeadTime());
      }
    }
  }
  TString fileid = path;
  fileid.Remove(0,fileid.Last('_')+1);
  fileid.Remove(fileid.Last('C'));
  double xmax1 = hLeD->GetXaxis()->GetBinCenter(hLeD->GetMaximumBin());
  double xmax2 = hLeS->GetXaxis()->GetBinCenter(hLeS->GetMaximumBin());
  TFile efile(path+ ".off.root","RECREATE");
  TGraph *gr = new TGraph();
  gr->SetPoint(0,xmax1-xmax2,  xmax1-xmax2);
  gr->SetName(fileid);
  gr->Write();
  efile.Write();
  efile.Close();
}
