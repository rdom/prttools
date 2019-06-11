#include "prttools.C"

void setGStyle(TGraph *g, int id){
  int c[]={kBlack,kRed+1,kGreen+1,kBlue,kOrange+1 ,kBlack,kRed+2,kGreen+2,kBlue+1,kOrange+2};
  
  g->SetLineColor(c[id%5]);
  g->SetMarkerColor(c[id%5+5]);
  g->SetMarkerStyle(20);
  g->SetMarkerSize(0.5);
  g->SetLineStyle(id>4 ? 2:1);
  g->SetTitle(";PiLas intensity [%];detected photons [#]");
}

void da_lab(TString infile = ""){
  
  TChain ch("lab"); ch.Add(infile);

  const int max=5;
  TGraph *gr[max][max];
  int it[max][max]={0};
  double var[max]={0};
  double inten,thresh,sigma;
  int res(0),ind(0);
  for(int v=0; v<max; v++){
    for(int h=0; h<max; h++){
      gr[v][h] = new TGraph();
      setGStyle(gr[v][h],ind++);
    }
  }
  
  TString nid[] = {"1 mV","2 mV","3 mV","4 mV","5 mV"};
  
  ch.SetBranchAddress("inten",&inten);
  ch.SetBranchAddress("thresh",&thresh);

  ch.SetBranchAddress("count1",&var[0]);
  ch.SetBranchAddress("count2",&var[1]);
  ch.SetBranchAddress("count3",&var[2]);
  ch.SetBranchAddress("count4",&var[3]);
  ch.SetBranchAddress("sigma",&var[4]);
  ch.SetBranchAddress("res",&res);
 
  int h =0;
  for(int i = 0; i < ch.GetEntries(); i++) {
    ch.GetEvent(i);
    if(thresh==600) h=0;
    if(thresh==1200) h=1;
    if(thresh==1800) h=2;
    if(thresh==2400) h=3;
    if(thresh==3000) h=4;    
    for(int v=0; v<5; v++) {
      gr[v][h]->SetPoint(it[v][h]++,inten,var[v]*(32/24.));
    }
  }

  for(int v=0; v<max; v++) for(int h=0; h<max; h++) gr[v][h]->Sort();

  prt_canvasAdd("mutl_simple",800,500);

  TLegend *leg = new TLegend(0.13,0.45,0.39,0.89);
  //TLegend *leg = new TLegend(0.53,0.15,0.86,0.51);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry("","opened:","");
  
  int v=0;
  for(int h=0; h<max; h++){    
    gr[v][h]->Draw(h>0? "PL same":"APL");
    leg->AddEntry(gr[v][h],nid[h],"lp");    
    gr[v][h]->GetYaxis()->SetRangeUser(0.001, 35);
    gr[v][h]->GetXaxis()->SetRangeUser(5,60);
  }
  leg->AddEntry("","masked:","");
  v=1;
  for(int h=0; h<max; h++){
    gr[v][h]->SetMarkerStyle(22);
    gr[v][h]->SetLineStyle(7);
    gr[v][h]->SetMarkerSize(0.8);
    gr[v][h]->Draw("PL same");
    leg->AddEntry(gr[v][h],nid[h],"lp");    
    gr[v][h]->GetYaxis()->SetRangeUser(0.001, 35);
    gr[v][h]->GetXaxis()->SetRangeUser(5,60);
  }

  
  // //gPad->SetLogy();
  leg->Draw();

  // TLegend *leg1 = new TLegend(0.64,0.23,0.97,0.38);
  // leg1->SetFillColor(0);
  // leg1->SetFillStyle(0);
  // leg1->SetBorderSize(0);
  // leg1->SetFillStyle(0);
  // prt_canvasAdd("resolution",800,500);
  // gr[4]->SetTitle(";PiLas intensity [%];time precision [ns]");
  // gr[4]->GetYaxis()->SetRangeUser(0.1, 0.3);
  // leg1->AddEntry(gr[4],nid[4],"lp");
  // leg1->AddEntry(gr[9],nid[9],"lp");
  // gr[4]->Draw();
  // gr[9]->Draw("same PL");
  // leg1->Draw();
  
  prt_savepath ="data/da_lab_19";
  prt_canvasSave(2);  
  
}
