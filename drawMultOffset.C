#define prt__sim
#include "../prtdirc/src/PrtHit.h"
#include "../prtdirc/src/PrtEvent.h"
#include "../../prttools/datainfo.C"
#include "../../prttools/prttools.C"

void drawMultOffset(TString infile="../build/hits.root"){

  if(!prt_init(infile,1,"data/drawMultOffset")) return;

  Int_t studyId, fileId(0), radiatorId(0), lensId(0), offset(0);
  Double_t mom(0), angle(0), z(0), x(0), xstep(0), ystep(0),mult(0);

  TFile *file = new TFile(infile.ReplaceAll(".root",".ro.root"),"recreate");
  TTree *tree= new TTree("proc","proc");
  tree->Branch("studyId", &studyId,"studyId/I");
  tree->Branch("fileId", &fileId,"fileId/I");
  tree->Branch("mom", &mom,"mon/D");
  tree->Branch("radiatorId", &radiatorId,"radiatorId/I");
  tree->Branch("lensId", &lensId,"lensId/I");
  tree->Branch("angle", &angle,"angel/D");
  tree->Branch("z", &z,"z/D");
  tree->Branch("x", &x,"x/D");
  tree->Branch("xstep", &xstep,"xstep/D");
  tree->Branch("ystep", &ystep,"ystep/D");
  tree->Branch("mult",&mult,"mult/D");
  tree->Branch("offset",&offset,"offset/I");

  TH1F * hMult[prt_nmcp];
  Int_t multa[prt_nmcp];
  for(auto i=0; i<prt_nmcp; i++){
    hMult[i]  = new TH1F(Form("mult_%i",i),"Mult;multiplicity [#];entries [#]",300,0,300);
  }

  TH1F * hTime  = new TH1F("le","LE; LE [ns]; entries [#]",1000,-50,50);

  infile.Remove(0, infile.Last('/')+1);  
  TObjArray *sarr = infile.Tokenize("_");
  if(sarr->GetEntries()==3){
    if(infile.Contains("th_")){      
      TString soffset = ((TObjString *) sarr->At(1))->GetName();
      offset = soffset.Atof();
    }
  }

  std::cout<<"offset "<<offset<<std::endl;
    
  PrtHit hit;
  for (auto ievent=0; ievent< prt_entries; ievent++){
    prt_nextEvent(ievent,1000);
    memset(multa, 0, sizeof(multa));    
    for(auto h=0; h<prt_event->GetHitSize(); h++){
      hit = prt_event->GetHit(h);
      Int_t mcpid = hit.GetMcpId();
      Double_t time = hit.GetLeadTime();
      if(time>-30 && time<0) multa[mcpid]++;
      hTime->Fill(time);
    }
    for(auto i=0; i<prt_nmcp; i++){
      hMult[i]->Fill(multa[i]);
    }
  }

  for(auto i=0; i<prt_nmcp; i++){
    mult = hMult[i]->GetMean();
    tree->Fill();
  }
  
  std::cout<<"hMult[i]->GetMean() "<<hMult[2]->GetMean()<<std::endl;
  
  
  prt_canvasAdd("time");
  hTime->Draw();
  
  prt_canvasAdd("mult_off");
  hMult[2]->Draw();
  prt_canvasSave(1,0);
  
  tree->Write();    
}
