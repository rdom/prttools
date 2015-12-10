#define prt__beam
#include "../prtdirc/src/PrtHit.h"
#include "../prtdirc/src/PrtEvent.h"
#include "prttools.C"
void drawPid(TString inFile = "../build/reco_spr.root", TString outFile="c_pid.root"){

  //fSavePath = "data/temp/pid";

  TString fileid(inFile);
  fileid.Remove(0,fileid.Last('/')+1);
  fileid.Remove(fileid.Last('.')-1);
  prt_data_info = getDataInfo(fileid);
  TString outdir=inFile;outdir.Remove(outdir.Last('/'));
  fSavePath = outdir+Form("/%da/%d",prt_data_info.getStudyId(),prt_data_info.getFileId());

  TChain ch("dirc"); ch.Add(inFile);
  Int_t tofPid;
  Double_t cangle,spr,trr,nph,par1,par2,par3,par4,par5,par6,test1,test2,theta,phi; 
  
  TH1F *hPi = new TH1F("hPi","hPi;#theta_{c} [rad]; entries [#]",500,0.6,0.9);
  TH1F *hP  = new TH1F("hP ","hP ;#theta_{c} [rad]; entries [#]",500,0.6,0.9);

  
  ch.SetBranchAddress("tofPid",&tofPid);
  ch.SetBranchAddress("spr",&spr);
  ch.SetBranchAddress("trr",&trr);
  ch.SetBranchAddress("nph",&nph);
  ch.SetBranchAddress("cangle",&cangle);
  ch.SetBranchAddress("par4",&par4);
  ch.SetBranchAddress("par5",&par5);
  ch.SetBranchAddress("par6",&par6);
  ch.SetBranchAddress("test1",&test1);
  ch.SetBranchAddress("test2",&test2);
  ch.SetBranchAddress("theta",&theta);
  ch.SetBranchAddress("phi",&phi);
  
  Int_t nent = ch.GetEntries();
  std::cout<<"# entries  "<< nent <<std::endl;
  std::cout<<"infor  "<< ch.GetTree()->GetTitle()<<std::endl;
  
  for (Int_t i = 0; i < nent; i++) {
    ch.GetEvent(i);
    if(nph<10) continue;
    if(tofPid==211) hPi->Fill(cangle);
    if(tofPid==2212) hP->Fill(cangle);    
  }
 
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  canvasAdd("mix_canglepid",800,400);
  hP->SetLineColor(2);
  hP->Draw();
  hPi->SetLineColor(4);
  hPi->Draw("same");

  canvasSave(1,0);
}
