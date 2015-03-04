#define prt__beam
#include "prttools.C"

#include "../mz-unpacker-BarrelDirc/TPrtHit.h"
#include "../mz-unpacker-BarrelDirc/TPrtEvent.h"


void mentries(TString inFile = "inM.root",TString outFile = "out.root", Int_t fileid=0 ){
 
  Int_t totalHits =0, channelHits=0;
  TFile *oFile = new TFile(outFile,"RECREATE");
  TTree *outTree = new TTree("T","CERN 1988 staff data");
  outTree->Branch("fileid",&fileid,"fileid/I");
  outTree->Branch("channelHits",&channelHits,"channelHits/I");
  outTree->Branch("totalHits",&totalHits,"totalHits/I");
 
  TPrtHit hit;
  PrtInit(inFile,0);

  for (Int_t ievent=0; ievent<fNEntries; ievent++){
    PrtNextEvent(ievent,100000);
    for(Int_t h=0; h<fEvent->GetHitsSize(); h++){
      for(Int_t m=0; m<fEvent->GetMultiplicity(h); m++){
    	hit = fEvent->GetHit(h,m);
    	if(hit.GetChannel()==1928) channelHits++;;	
    	totalHits++;
      }
    }
  }

  std::cout<<"inFile "<< inFile << "  ch1928  "<< channelHits << "   total  "<<totalHits <<std::endl;
  
  outTree->Fill();
  outTree->Write();

}


