#include "prttools.C"
void drawProcRes(TString inFile = "../data/res151.root"){
  prt_savepath = "data/drawProcRes";
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
    std::cout<<"mult "<<mult<< " t "<< angle<<std::endl;
    
    if(angle<160)
      gNph->SetPoint(point++,angle,mult);
  }
  gNph->Sort();

  TString names[500];
  names[150]="bar 2LS @ 7 GeV/c";
  names[151]="bar 3LS @ 7 GeV/c";
  names[152]="plate WL @ 7 GeV/c";
  names[153]="plate 2LC @ 7 GeV/c";
  names[154]="bar 2LAG @ 7 GeV/c";
  names[160]="bar 3LS @ 5 GeV/c";
  names[161]="plate WL @ 5 GeV/c";
  names[162]="plate 2LC @ 5 GeV/c";

  names[201]="sim, plate w/o lens @ 7 GeV/c";
  names[202]="sim, plate with 2LCL @ 7 GeV/c";
  names[219]="data, plate with 2LCL @ 7 GeV/c";
  names[221]="data, plate w/o lens @ 7 GeV/c";

  names[400]="Fast angle scan, bar, 3LS lens  @ 7 GeV/c";
  
  
  
  gNph->SetTitle(names[studyId]);
  
  gNph->SetLineColor(38);
  gNph->SetMarkerStyle(20);
  gNph->SetMarkerSize(0.7);
  gNph->GetYaxis()->SetTitle("multiplicity [#]");


  gNph->GetXaxis()->SetLabelSize(0.05);
  gNph->GetXaxis()->SetTitleSize(0.06);
  gNph->GetXaxis()->SetTitleOffset(0.84);

  gNph->GetYaxis()->SetLabelSize(0.05);
  gNph->GetYaxis()->SetTitleSize(0.06);
  gNph->GetYaxis()->SetTitleOffset(0.7);

  gNph->GetXaxis()->SetTitle("#theta_{track} [deg]");

  TCanvas* c2 = new TCanvas(Form("mult_%d",studyId),"c2",800,500);c2->SetBottomMargin(0.12);
  gNph->Draw("APL");
  gNph->GetYaxis()->SetRangeUser(0,100);
  gNph->GetXaxis()->SetRangeUser(10,170);
  
  prt_canvasAdd(c2);

  // TLegend *leg = new TLegend(0.2,0.7,0.5,0.9);
  // leg->SetFillColor(0);
  // leg->SetFillStyle(0);
  // leg->SetBorderSize(0);
  // leg->SetFillStyle(0);
  // leg->AddEntry("S221",names[221],"lp");
  // leg->AddEntry("S219",names[219],"lp");
  // leg->AddEntry("S201",names[201],"lp");
  // leg->AddEntry("S202",names[202],"lp");   
  // leg->Draw();
  
  prt_canvasSave(0,1);
}
