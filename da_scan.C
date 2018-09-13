#include "prttools.C"

void setGStyle(TGraph *g, int id){
  int c[]={kBlack,kRed+1,kGreen+1,kBlue,kOrange+1,kCyan+1,kBlue,kOrange+1,1,1,1,1,1,1,1,1};
  
  g->SetLineColor(c[id]);
  g->SetMarkerColor(id>1 ? ++c[id]:c[id]);
  g->SetMarkerStyle(20);
  g->SetMarkerSize(0.5);
  g->SetLineStyle(id>5 ? 2:1);

  g->SetTitle("shape a");
}

void da_js(){
  prt_savepath = "data/joinScan";
  TGraph *gret;

  TString gpath="~/sim4/d/proc/jul18";

  const int max(2);
  TString wid[max][3]{
    {gpath+"/403/reco*R_spr.root","(fR1==1)&&(fR2==0.00)","beam data"},
    {gpath+"/403_00/reco*S_spr.root","(fR1==1)&&(fR2==0.00)","geant4 sim"}       
  };
  
  TLegend *leg = new TLegend(0.65,0.61,0.85,0.89);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  
  TString nid = "403_";

  prt_canvasAdd(nid+"nph",800,500);
  for(int j=0; j<max; j++){
    gret = da_scan(wid[j][0],1,wid[j][1]);
    setGStyle(gret,j);
    gret->Draw((j==0)?"APL":"PL same");
    leg->AddEntry(gret,wid[j][2],"lp");
  }
  leg->Draw();

  prt_canvasAdd(nid+"spr",800,500);
  for(int j=0; j<max; j++){
    gret = da_scan(wid[j][0],2,wid[j][1]);
    setGStyle(gret,j);
    gret->Draw((j==0)?"APL":"PL same");
    leg->AddEntry(gret,wid[j][2],"lp");
  }
  leg->Draw();

  prt_canvasAdd(nid+"trr",800,500);
  for(int j=0; j<max; j++){
    gret = da_scan(wid[j][0],3,wid[j][1]);
    setGStyle(gret,j);
    gret->Draw((j==0)?"APL":"PL same");
    leg->AddEntry(gret,wid[j][2],"lp");
  }
  leg->Draw();
  prt_canvasSave(2,0);
}

TGraph* da_scan(TString inFile = "r_spr.root", const Int_t func=0, TString scut=""){
  
  TString outdir=inFile;outdir.Remove(outdir.Last('/'));
  TString sfile=inFile; sfile.Remove(0,sfile.Last('/')+1);
  TString sstudy=outdir; sstudy.Remove(0,sstudy.Last('/'));
  prt_savepath = outdir+sstudy+"r";
  if(inFile.Contains("S.root"))  prt_savepath = outdir+sstudy+"s";
  outdir.ReplaceAll("*", "");
  TString outFile=outdir+"/j"+sfile;

  TChain ch("dirc"); ch.Add(inFile);
  Double_t cangle,spr,trr,nph,par1,par2,par3,par4,par5,par6,test1,test2,theta,phi; 
  
  TGraph *gSpr = new TGraph();
  TGraph *gNph = new TGraph();
  TGraph *gTrr = new TGraph();

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

  Int_t it(0);
  for (Int_t i = 0; i < nent; i++) {
    ch.GetEvent(i);
    if(spr==0 || theta>156) continue;
    gSpr->SetPoint(it,theta,TMath::Abs(spr));
    gNph->SetPoint(it,theta,nph);
    gTrr->SetPoint(it,theta,TMath::Abs(trr));
    it++;
  }
  gSpr->Sort();
  gNph->Sort();
  gTrr->Sort();

  gSpr->SetLineColor(38);
  gNph->SetLineColor(38);
  gTrr->SetLineColor(38);
  gSpr->SetMarkerStyle(20);
  gNph->SetMarkerStyle(20);
  gTrr->SetMarkerStyle(20);
  gSpr->SetMarkerSize(0.7);
  gNph->SetMarkerSize(0.7);
  gTrr->SetMarkerSize(0.7);
  gNph->GetYaxis()->SetRangeUser(0,100);
  gSpr->GetYaxis()->SetRangeUser(0,20);
  gTrr->GetYaxis()->SetRangeUser(0,5);

  gNph->GetXaxis()->SetLimits(15,145);
  gSpr->GetXaxis()->SetLimits(15,145);
  gTrr->GetXaxis()->SetLimits(15,145);
  
  gNph->GetXaxis()->SetRangeUser(15,145);
  gSpr->GetXaxis()->SetRangeUser(15,145);
  gTrr->GetXaxis()->SetRangeUser(15,145);

  gSpr->GetYaxis()->SetTitle("SPR [mrad]");
  gNph->GetYaxis()->SetTitle("multiplicity [#]");
  gTrr->GetYaxis()->SetTitle("#sigma_{#theta_{C} tr} [mrad]");
  
  gSpr->GetXaxis()->SetLabelSize(0.05);
  gSpr->GetXaxis()->SetTitleSize(0.06);
  gSpr->GetXaxis()->SetTitleOffset(0.84);

  gTrr->GetXaxis()->SetLabelSize(0.05);
  gTrr->GetXaxis()->SetTitleSize(0.06);
  gTrr->GetXaxis()->SetTitleOffset(0.84);

  gNph->GetXaxis()->SetLabelSize(0.05);
  gNph->GetXaxis()->SetTitleSize(0.06);
  gNph->GetXaxis()->SetTitleOffset(0.84);

  gSpr->GetYaxis()->SetLabelSize(0.05);
  gSpr->GetYaxis()->SetTitleSize(0.06);
  gSpr->GetYaxis()->SetTitleOffset(0.7);

  gTrr->GetYaxis()->SetLabelSize(0.05);
  gTrr->GetYaxis()->SetTitleSize(0.06);
  gTrr->GetYaxis()->SetTitleOffset(0.7);

  gNph->GetYaxis()->SetLabelSize(0.05);
  gNph->GetYaxis()->SetTitleSize(0.06);
  gNph->GetYaxis()->SetTitleOffset(0.7);


  gSpr->GetXaxis()->SetTitle("#theta_{track} [deg]");
  gNph->GetXaxis()->SetTitle("#theta_{track} [deg]");
  gTrr->GetXaxis()->SetTitle("#theta_{track} [deg]");

     
  if(func){
    TGraph *gret;
    if(func==1) gret=gNph;
    if(func==2) gret=gSpr;
    if(func==3) gret=gTrr;    
    return gret;
  }else{  
    TFile *file = new TFile(outFile,"RECREATE");
    TCanvas* c1 = new TCanvas("spr","spr",800,400);c1->SetBottomMargin(0.12);
    gSpr->Draw("APL");
    prt_canvasAdd(c1);
    TCanvas* c2 = new TCanvas("nph","nph",800,400);c2->SetBottomMargin(0.12);
    gNph->Draw("APL");
    prt_canvasAdd(c2);
    TCanvas* c3 = new TCanvas("trr","trr",800,400);c3->SetBottomMargin(0.12);
    gTrr->Draw("APL");
    prt_canvasAdd(c3);

    prt_canvasSave(1,0);
   
    file->cd();
    c1->Write();
    c2->Write();
    c3->Write();
    file->Close();
  }
  return 0;
}
