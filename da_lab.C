#include "prttools.C"

void setGStyle(TGraph *g, int id){
  //int c[]={kBlack,kRed+1,kGreen+1,kBlue,kOrange+1,kCyan+1,kBlue,kOrange+1,1,1,1,1,1,1,1,1};
  int c[]={kBlack,kRed+1,kGreen+1,kBlue,kBlack,kRed+1,kGreen+1,kBlue,1,1,1,1,1,1,1,1};
  
  g->SetLineColor(c[id]);
  g->SetMarkerColor(id>1 && id!=4 ? ++c[id]:c[id]);
  g->SetMarkerStyle(20);
  g->SetMarkerSize(0.5);
  g->SetLineStyle(id>3 ? 2:1);
  g->SetTitle(";PiLas intensity [%];detected photons [#]");
}

void da_lab(TString infile = ""){
  
  TChain ch("lab"); ch.Add(infile);

  

  const int max=8;
  TGraph *gr[max];
  int it[max]={0};
  double var[max]={0};
  double inten;
  int res(0);
  for(int v=0; v<max; v++){
    gr[v] = new TGraph();
    setGStyle(gr[v],v);
  }
  
  TString nid[] = {"0 Om; opened","0 Om; masked","0 Om; opened + time cut","0 Om; masked + time cut",
		   "75 Om; opened","75 Om; masked","75 Om; opened + time cut","75 Om; masked + time cut"};
  
  ch.SetBranchAddress("inten",&inten);

  ch.SetBranchAddress("count1",&var[0]);
  ch.SetBranchAddress("count2",&var[1]);
  ch.SetBranchAddress("count3",&var[2]);
  ch.SetBranchAddress("count4",&var[3]);
  ch.SetBranchAddress("res",&res);
  
  for(int i = 0; i < ch.GetEntries(); i++) {
    ch.GetEvent(i);
    std::cout<<"res "<<res<<" inten "<<inten<<" "<< var[3]<<std::endl;
    
    for(int v=0; v<4; v++) {
      if(res==0) gr[v]->SetPoint(it[v]++,inten,var[v]);
      if(res==75) gr[v+4]->SetPoint(it[v+4]++,inten,var[v]);
    }
  }

  for(int v=0; v<max; v++) gr[v]->Sort();

  prt_canvasAdd("mutl_log",800,500);

  //TLegend *leg = new TLegend(0.13,0.56,0.39,0.89);
  TLegend *leg = new TLegend(0.53,0.15,0.86,0.51);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  
  for(int v=0; v<max; v++){
    gr[v]->Draw(v>0? "PL same":"APL");
    leg->AddEntry(gr[v],nid[v],"lp");

    // gr[v]->SetMinimum(10);
    // gr[v]->SetMaximum(1000000);
    //gr[v]->GetYaxis()->SetRangeUser(0, 22);
    gr[v]->GetYaxis()->SetRangeUser(0.001, 30);
    gr[v]->GetXaxis()->SetRangeUser(5,60);
  }
  
  
  gPad->SetLogy();
  leg->Draw();

  prt_savepath ="data/da_lab_11";
  prt_canvasSave(2);  
  
}
