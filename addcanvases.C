#include "prttools.C"
void addcanvases(TString f1n="cspr_150S.root", TString f2n="spr_150R.root"){

  f2n = f1n; f2n.ReplaceAll("S.root","R.root");

  TString outdir=f1n;outdir.Remove(outdir.Last('/'));
  TString sstudy=outdir; sstudy.Remove(0,sstudy.Last('/'));
  fSavePath = outdir+sstudy;
  
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
    TLegend *leg = new TLegend(0.2,0.7,0.5,0.9);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    carr1[i]->Draw();
    leg->AddEntry(carr1[i],"sim","lp");
    leg->Draw();
    carr1[i]->SetName("mix_"+carr1[i]->GetName());
    canvasAdd(carr1[i]);
    TIter next(carr2[i]->GetListOfPrimitives());
    TObject *obj;

    while((obj = next())){
      // if(obj->InheritsFrom("TH1F")){
      // 	TH1F *h = (TH1F*)obj;
      // 	std::cout<<"name "<< h->GetName() <<std::endl;      
      // 	h->SetLineStyle(7);
      // 	h->SetLineWidth(2);
      //   h->Draw("same");
      // }
      if(obj->InheritsFrom("TGraph")){
	TGraph *h = (TGraph*)obj;
	std::cout<<"name "<< h->GetName() <<std::endl;      
	h->SetLineColor(32);
	h->SetMarkerColor(4);
	//	h->SetLineWidth(2);
        h->Draw("same PL");
	leg->AddEntry(h,"beam data","lp");
	leg->Draw();
      }

    }
  }
  std::cout<<"save all  " <<std::endl;

  canvasSave(0,1);
}
