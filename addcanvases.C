#include "prttools.C"
void addcanvases(TString f1n="cspr_150S.root", TString f2n="spr_150R.root", Int_t data=-1){
  
  TString outdir=f1n;outdir.Remove(outdir.Last('/'));
  prt_savepath = outdir;
  if(data==-1) {
    f2n = f1n; f2n.ReplaceAll("S.root","R.root");
    TString sstudy=outdir; sstudy.Remove(0,sstudy.Last('/'));
    prt_savepath = outdir+sstudy+"a";
  }
  
  std::cout<<"reading  "<<f1n <<std::endl;
  std::cout<<"reading  "<<f2n <<std::endl;
  
  const Int_t narr = 20;
  gStyle->SetOptStat(0); 
  gStyle->SetOptTitle(0); 

  TFile *f1 = TFile::Open(f1n);
  TIter next1(f1->GetListOfKeys());
  TKey *key1;
  Int_t it1 = 0;
  TCanvas *carr1[narr];
  while((key1 = (TKey*)next1())) {
    TClass *cl = gROOT->GetClass(key1->GetClassName());
    if (!cl->InheritsFrom("TCanvas")) continue;
    carr1[it1] = (TCanvas*)key1->ReadObj();
    it1++;
  }

  TFile *f2 = TFile::Open(f2n);
  TIter next2(f2->GetListOfKeys());
  TKey *key2;
  Int_t it2 = 0;
  TCanvas *carr2[narr];
  while ((key2 = (TKey*)next2())) {
    TClass *cl = gROOT->GetClass(key2->GetClassName());
    if (!cl->InheritsFrom("TCanvas")) continue;
    carr2[it2] = (TCanvas*)key2->ReadObj();
    it2++;
  }

  for(Int_t i=0; i<it2; i++){
    TLegend *leg = new TLegend(0.5,0.7,0.8,0.9);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    carr2[i]->Draw();
    TH1F *tt = new TH1F(); tt->SetMarkerStyle(20);tt->SetMarkerSize(0.8);
    leg->AddEntry(tt,"beam data","lp");
    carr2[i]->SetName(Form("mix_%s",carr2[i]->GetName()));
    prt_canvasAdd(carr2[i]);

    TIter next(carr1[i]->GetListOfPrimitives());
    TObject *obj;
    while((obj = next()) ){
      // if(obj->InheritsFrom("TH1F")){
      // 	TH1F *h = (TH1F*)obj;
      // 	std::cout<<"name "<< h->GetName() <<std::endl;      
      // 	h->SetLineStyle(7);
      // 	h->SetLineWidth(2);
      //   h->Draw("same");
      // }
      if(obj->InheritsFrom("TGraph")){
	TGraph *h = (TGraph*)obj;
	h->SetName("g2");
	std::cout<<"name "<< h->GetName() <<std::endl;      
	h->SetLineColor(32);
	h->SetMarkerColor(2);
	h->SetMarkerSize(0.8);
	//	h->SetLineWidth(2);
	h->Draw("same PL");
	leg->AddEntry(h,"simulation","lp");
	leg->Draw();
      }

    }
  }
  std::cout<<"save all  " <<std::endl;
  prt_canvasSave(0,3);
}
