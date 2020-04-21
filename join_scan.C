#include "prttools.C"
#include "draw_scan.C"

void setGStyle(TGraph *g, int id){
  int c[]={kBlack,kRed+1,kBlue,kGreen+1,kOrange+1,kCyan+1,kBlue,kOrange+1,1,1,1,1,1,1,1,1};
  
  g->SetLineColor(c[id]);
  g->SetMarkerColor(id>1 ? c[id]+1:c[id]);
  g->SetMarkerStyle(20);
  g->SetMarkerSize(0.75);
  g->SetLineStyle(id>5 ? 2:1);

  g->SetTitle(" ");
}

void join_scan(){
  prt_savepath = "data/join_scan";

  const int max(2);
  TString wid[max][3]{
    {"~/sim4/d/proc/jul18/403/reco*R.root","",   "beam data"},
    {"~/sim4/d/proc/jul18/403/reco*S.root","",   "geant4 simulation"}   
  };
  
  TLegend *leg = new TLegend(0.4,0.65,0.9,0.8);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  
  TString nid[] = {"nph","sep"};

  TGraph *gret;
  for(int i=0; i<2; i++){
    prt_canvasAdd(nid[i],800,500);
    for(int j=0; j<max; j++){
      gret = draw_scan(wid[j][0],i,0,wid[j][1]);
      setGStyle(gret,j);
      gret->Draw((j==0)?"APL":"PL same");
      if(i==0) leg->AddEntry(gret,wid[j][2],"lp");
    }
    leg->Draw();
  }
  
  prt_canvasSave(1);
}

