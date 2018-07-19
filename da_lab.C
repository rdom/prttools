#include "prttools.C"

void setGStyle(TGraph *g, int id){
  int c[]={kBlack,kRed+1,kGreen+1,kBlue,kOrange+1,kCyan+1,kBlue,kOrange+1,1,1,1,1,1,1,1,1};
  
  g->SetLineColor(c[id]);
  g->SetMarkerColor(id>1 ? ++c[id]:c[id]);
  g->SetMarkerStyle(20);
  g->SetMarkerSize(0.5);
  g->SetLineStyle(id>5 ? 2:1);
  g->SetTitle(";PiLas intencity [%];detected photons [#]");
}

void da_lab(TString infile = ""){
  
  TChain ch("lab"); ch.Add(infile);

  

  const int max=6;
  TGraph *gr[max];
  int it[max]={0};
  double var[max]={0};
  double inten, lpere(0);
  for(int v=0; v<5; v++){
    gr[v] = new TGraph();
    setGStyle(gr[v],v);
  }
  
  TString nid[] = {"ch 551 (masked)","ch 423 (open)","open pixels","masked pixels","per event for open mcp"};
  
  ch.SetBranchAddress("inten",&inten);

  ch.SetBranchAddress("un551",&var[0]);
  ch.SetBranchAddress("ln423",&var[1]);
  ch.SetBranchAddress("ua000",&var[2]);
  ch.SetBranchAddress("la000",&var[3]);
  ch.SetBranchAddress("lpere",&lpere);
  
  for(int i = 0; i < ch.GetEntries(); i++) {
    ch.GetEvent(i);
    std::cout<<"inten "<<inten<<" "<< var[3]<<std::endl;
    
    for(int v=0; v<4; v++) gr[v]->SetPoint(it[v]++,inten,var[v]*64);
    gr[4]->SetPoint(it[4]++,inten,lpere*64);
  }

  for(int v=0; v<5; v++) gr[v]->Sort();

  prt_canvasAdd("lab_nph_10ns_mcp7_norm",800,500);

  //TLegend *leg = new TLegend(0.50,0.57,0.85,0.83);
  TLegend *leg = new TLegend(0.20,0.57,0.5,0.8);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  
  for(int v=2; v<4; v++){
    gr[v]->Draw(v>2? "PL":"APL");
    leg->AddEntry(gr[v],nid[v],"lp");

    // gr[v]->SetMinimum(1);
    // gr[v]->SetMaximum(1000000);
    gr[v]->GetYaxis()->SetRangeUser(0, 55);

  }

  //gPad->SetLogy();
  leg->Draw();

  prt_savepath ="data/da_lab";
  prt_canvasSave(2);
  
  
}
