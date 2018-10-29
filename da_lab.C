#include "prttools.C"

void setGStyle(TGraph *g, int id){
  //int c[]={kBlack,kRed+1,kGreen+1,kBlue,kOrange+1,kCyan+1,kBlue,kOrange+1,1,1,1,1,1,1,1,1};
  int c[]={kBlack,kRed+1,kGreen+1,kBlue,kBlack,kBlack,kRed+1,kGreen+1,kBlue,kRed+1,1,1,1,1,1,1,1,1};
  
  g->SetLineColor(c[id]);
  g->SetMarkerColor(id>1 && id!=4 && id!=5 ? ++c[id]:c[id]);
  g->SetMarkerStyle(20);
  g->SetMarkerSize(0.5);
  g->SetLineStyle(id>4 ? 2:1);
  g->SetTitle(";PiLas intensity [%];detected photons [#]");
}

void da_lab(TString infile = ""){
  
  TChain ch("lab"); ch.Add(infile);
  
  const int max=10;
  TGraph *gr[max];
  int it[max]={0};
  double var[max]={0};
  double inten,sigma;
  int res(0);
  for(int v=0; v<max; v++){
    gr[v] = new TGraph();
    setGStyle(gr[v],v);
  }
  
  TString nid[] = {"0 #Omega; opened","0 #Omega; masked","0 #Omega; opened + time cut","0 #Omega; masked + time cut", "0 #Omega", 
		   "75 #Omega; opened","75 #Omega; masked","75 #Omega; opened + time cut","75 #Omega; masked + time cut","75 #Omega"};
  
  ch.SetBranchAddress("inten",&inten);

  ch.SetBranchAddress("count1",&var[0]);
  ch.SetBranchAddress("count2",&var[1]);
  ch.SetBranchAddress("count3",&var[2]);
  ch.SetBranchAddress("count4",&var[3]);
  ch.SetBranchAddress("sigma",&var[4]);
  ch.SetBranchAddress("res",&res);
  
  for(int i = 0; i < ch.GetEntries(); i++) {
    ch.GetEvent(i);
    std::cout<<"res "<<res<<" inten "<<inten<<" "<< var[3]<<std::endl;
    
    for(int v=0; v<5; v++) {
      if(res==0) gr[v]->SetPoint(it[v]++,inten,var[v]);
      if(res==75) gr[v+5]->SetPoint(it[v+5]++,inten,var[v]);
    }
  }

  for(int v=0; v<max; v++) gr[v]->Sort();

  prt_canvasAdd("mutl_simple",800,500);

  TLegend *leg = new TLegend(0.13,0.56,0.39,0.89);
  //TLegend *leg = new TLegend(0.53,0.15,0.86,0.51);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  
  for(int v=0; v<max; v++){
    if(v==2 || v==3 || v==4 || v>6)continue;
    gr[v]->Draw(v>0? "PL same":"APL");
    leg->AddEntry(gr[v],nid[v],"lp");
    
    gr[v]->GetYaxis()->SetRangeUser(0.001, 35);
    gr[v]->GetXaxis()->SetRangeUser(5,60);
  }

  
  //gPad->SetLogy();
  leg->Draw();

  TLegend *leg1 = new TLegend(0.64,0.23,0.97,0.38);
  leg1->SetFillColor(0);
  leg1->SetFillStyle(0);
  leg1->SetBorderSize(0);
  leg1->SetFillStyle(0);
  prt_canvasAdd("resolution",800,500);
  gr[4]->SetTitle(";PiLas intensity [%];time precision [ns]");
  gr[4]->GetYaxis()->SetRangeUser(0.1, 0.3);
  leg1->AddEntry(gr[4],nid[4],"lp");
  leg1->AddEntry(gr[9],nid[9],"lp");
  gr[4]->Draw();
  gr[9]->Draw("same PL");
  leg1->Draw();
  
  prt_savepath ="data/da_lab_2";
  prt_canvasSave(2);  
  
}
