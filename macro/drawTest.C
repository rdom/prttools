#define prt__sim
#include "../../prttools/prttools.C"

void drawTest(TString infile="../build/hits.root", Int_t flag = 0){
  std::cout<<"new scan"<<std::endl;
  
  fSaveFlag = 2; //2
  fInfo = "drawTest.C \n";
  PrtInit(infile,1);

  TH2F *hHits = new TH2F(infile+"1","hHist",500,-3,30,300,-10,10 );
  TH1F *hTime = new TH1F("hTime","hTime",500,0,50);

  PrtHit fHit;
  Double_t z = 0;
  Int_t tVal, angle = 0;
  TVector3 mom;
  for (Int_t ievent=0; ievent<fCh->GetEntries(); ievent++){
    PrtNextEvent(ievent,1000);
    if(ievent == 0) fPath = createDir("data/stepL0",fInfo, fSaveFlag);  
    tVal = fEvent->GetTest();
    angle = fEvent->GetAngle() + 0.01;
    for(Int_t h=0; h<fEvent->GetHitSize(); h++){
      fHit = fEvent->GetHit(h);
      Double_t time = fHit.GetLeadTime();
      if(time<50){
	int mcpid = fHit.GetMcpId();
	int pixid = fHit.GetPixelId();
	fhDigi[mcpid]->Fill(pixid/8, pixid%8);

	z = fHit.GetGlobalPos().Z();
	hHits->Fill(fHit.GetGlobalPos().X()/10.,fHit.GetGlobalPos().Y()/10.);
	hTime->Fill(time); 	  
      }	
    }
  } 

  Int_t tIVal = tVal;
  TString title = Form(" a=%d^{o} d=%f",angle,tIVal);
  if(flag==0) title = Form(" a=%d^{o}",angle);
  TCanvas* c2 = new TCanvas(infile+"c2","c2",800,500);
  axisHits800x500(hHits);
  hHits->SetTitle(Form("%d hits",(Int_t)hHits->GetEntries())+title);
  hHits->Draw("colz");
  save(c2,fPath,Form("a_hit0_a%d_s%d",angle,tIVal),fInfo,fSaveFlag,1);

  drawDigi();
  save(cDigi,fPath,Form("a_hit_a%d_s%d",angle,tIVal),fInfo,fSaveFlag,1);

  TCanvas* c3 = new TCanvas(infile+"c3","c3",800,500);
  axisTime800x500(hTime);
  hTime->SetTitle(Form("%d hits",(Int_t)hTime->GetEntries())+title);
  hTime->Draw();
  save(c3,fPath,Form("a_time_a%d_s%d",angle,tIVal),fInfo,fSaveFlag,1);
}


