#include "prttools.C"
void drawProcRes(TString inFile = "../data/res151.root"){
  fSavePath = "mult";
  TChain ch("proc"); ch.Add(inFile);

  Int_t studyId = 0, fileId=0,radiatorId, lensId=0;
  Double_t mom=0, angle=0, z =0, x= 0, xstep=0, ystep =0, mult=0;
   
  TGraph *gNph = new TGraph();

  ch.SetBranchAddress("studyId",&studyId);
  ch.SetBranchAddress("fileId",&fileId);
  //  ch.SetBranchAddress("radiatorId",&radiatorId);
  ch.SetBranchAddress("lensId",&lensId);
  ch.SetBranchAddress("mom",&mom);
  ch.SetBranchAddress("angle",&angle);
  ch.SetBranchAddress("z",&z);
  ch.SetBranchAddress("x",&x);
  ch.SetBranchAddress("xstep",&xstep);
  ch.SetBranchAddress("ystep",&ystep);
  ch.SetBranchAddress("mult",&mult);
  
  Int_t nent = ch.GetEntries();
  std::cout<<"# entries  "<< nent <<std::endl;
  Int_t point(0);
  for (Int_t i = 0; i < nent; i++) {
    ch.GetEvent(i);
    if(angle<160)
      gNph->SetPoint(point++,angle,mult);
  }
  gNph->Sort();

  gNph->SetLineColor(38);
  gNph->SetMarkerStyle(20);
  gNph->SetMarkerSize(0.7);
  gNph->GetYaxis()->SetRangeUser(0,170);
  gNph->GetXaxis()->SetRangeUser(10,170);
  gNph->GetYaxis()->SetTitle("multiplicity [#]");


  gNph->GetXaxis()->SetLabelSize(0.05);
  gNph->GetXaxis()->SetTitleSize(0.06);
  gNph->GetXaxis()->SetTitleOffset(0.84);

  gNph->GetYaxis()->SetLabelSize(0.05);
  gNph->GetYaxis()->SetTitleSize(0.06);
  gNph->GetYaxis()->SetTitleOffset(0.7);

  gNph->GetXaxis()->SetTitle("#theta_{track} [#circ]");

  TCanvas* c2 = new TCanvas(Form("mult_%d",studyId),"c2",800,500);c2->SetBottomMargin(0.12);
  gNph->Draw("APL");
  canvasAdd(c2);
 
  canvasSave(1,1);
}
