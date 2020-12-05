#include "prttools.C"

void drawMapSpr(TString in="tdata/rt.root"){
  prt_setRootPalette(1);
  TChain ch("reco"); ch.Add(in);
  double spr_pi,nph,test1,test2,theta,phi,sep_gr,sep_ti; 

  ch.SetBranchAddress("sep_gr",&sep_gr);
  ch.SetBranchAddress("sep_ti",&sep_ti);
  ch.SetBranchAddress("spr_pi",&spr_pi);
  ch.SetBranchAddress("test1",&test1);
  ch.SetBranchAddress("test2",&test2);
  ch.SetBranchAddress("theta",&theta);
  ch.SetBranchAddress("phi",&phi);

  gStyle->SetOptStat(0);
  TH2F *hSpr = new TH2F("hSpr",";#Delta#theta [mrad];#Delta#varphi [mrad]",21,-21,21,21,-21,21);
  TH2F *hSep_gr = new TH2F("hSep_gr",";#Delta#theta [mrad];#Delta#varphi [mrad]",21,-21,21,21,-21,21);
  TH2F *hSep_ti = new TH2F("hSep_ti",";#Delta#theta [mrad];#Delta#varphi [mrad]",21,-21,21,21,-21,21);

  for (auto i = 0; i < ch.GetEntries(); i++) {
    ch.GetEvent(i);
    hSpr->Fill(test1*1000,test2*1000,spr_pi);
    hSep_gr->Fill(test1*1000,test2*1000,sep_gr);
    hSep_ti->Fill(test1*1000,test2*1000,sep_ti);
  }

  TString nid = Form("#theta = %1.2f deg",theta);
  hSpr->SetTitle(nid);
  hSep_gr->SetTitle(nid);
  hSep_ti->SetTitle(nid);
  
  prt_canvasAdd(Form("hspr_%d",(int)(theta+0.6)),800,500);
  hSpr->Draw("colz");

  
  prt_canvasAdd(Form("hsep_ti_%d",(int)(theta+0.6)),800,500);
  hSep_ti->Draw("colz");


  prt_canvasAdd(Form("hsep_gr_%d",(int)(theta+0.6)),800,500);
  hSep_gr->Draw("colz");
  
  hSep_gr->Smooth();
  int maxbin = hSep_gr->GetMaximumBin();
  int binx,biny,binz;
  hSep_gr->GetBinXYZ(maxbin, binx, biny, binz);
    
  double maxx = hSep_gr->GetXaxis()->GetBinCenter(binx);
  double maxy = hSep_gr->GetYaxis()->GetBinCenter(biny);

  
  TEllipse *el3 = new TEllipse(maxx,maxy,0.6,1);
  el3->SetFillColor(1);
  el3->Draw();
  std::cout<<"theta "<<theta<<"  "<<Form("%1.2f",0.001*maxx*TMath::RadToDeg()) <<"  "<<Form("%1.2f",0.001*maxy*TMath::RadToDeg())<<std::endl;
  
  
  // int binx,biny,binz;  
  // hSpr->GetBinXYZ(hSpr->GetMinimumBin(),binx,biny,binz);
  // Double_t binxCenter = hSpr->GetXaxis()->GetBinCenter(binx);
  // Double_t binyCenter = hSpr->GetYaxis()->GetBinCenter(biny);
  // TEllipse *el = new TEllipse(binxCenter,binyCenter,.5,.5);
  // el->Draw();


  // int n =25;
  // double th[]  ={  20,  25,  30,  35,   40,   45,   50,   55,   60,   65,   70,   75,   80,   85,   90,   95,  100,  105,  110, 115,   120,  125,  130,  135,  140};
  // //p
  // double a[]   ={5.00,9.00,7.00,4.00, 5.00, 3.00, 4.00, 4.00, -4.0,-7.00, -4.0, -2.00, -5.0, -6.0, 5.00, 5.00, 12.0, 16.0, 7.00, 10.0, 10.0, 5.00, 0.00, 4.00, 8.00};
  // double ph[]  ={20.0,15.0,14.0,13.0, 11.0, 9.00, 11.0, 11.0, 8.00, 6.00, 8.00,  6.00, 8.00, 5.00, 8.00, 7.00, 6.00, 3.00, 5.00, 8.00, 6.00, 10.0, 10.0, 7.00, 8.00};
  // for(int i=0; i<n; i++){
  //   std::cout<<" "<<ph[i]/TMath::Pi()*0.180<<std::endl;
    
  // }
  
  // // int n =15;
  // // double th[]  ={  20,   25,  30,  40,   50,   60,   70,   80,  90,  100,  110,  120,  130,  140,  150};
  // // //p
  // // double a[]   ={3.84,-2.38,3.98,1.60,-1.49,-2.60,-0.81,-2.98,0.56,-5.00, 4.59, 3.48, 7.25, 4.69, 1.88};
  // // double ph[]  ={4.45,6.11 ,5.36,6.03, 6.13, 3.31, 4.55, 3.14,0.00, -0.5,-5.46,-6.23,-5.29,-5.47,-6.03};

  // // //pi
  // // double a_pi[]   ={-2.23,0.43,0.18,1.60,-3.71,-1.01, 2.94, 1.10,0.00,-5.50, 1.30, 0.3, 7.00, 5.90, 1.9};
  // // double ph_pi[]  ={3.11 ,5.21,4.34,5.82, 2.25, 2.31, 1.51, 1.90,0.00,-1.30,-3.70,-4.2,-5.10,-5.20,-5.1};

  // TGraph *gadiff = new TGraph(n,th,a);
  // TGraph *gphi = new TGraph(n,th,ph);
  // gadiff->SetTitle(";#theta [degree]; #theta - #theta_{fit} [mrad]");
  // gadiff->SetMarkerStyle(20);
  // gadiff->SetMarkerSize(0.8);
  // gadiff->SetLineColor(kBlack);
  // gadiff->SetMarkerColor(kBlack);
  // gphi->SetTitle(";#theta [degree]; #varphi_{fit} [mrad]");
  // gphi->SetMarkerStyle(20);
  // gphi->SetMarkerSize(0.8);
  // gphi->SetLineColor(kBlack);
  // gphi->SetMarkerColor(kBlack);
  
  // // TGraph *gadiff_pi = new TGraph(n,th,a_pi);
  // // TGraph *gphi_pi = new TGraph(n,th,ph_pi);
  // // gadiff_pi->SetTitle(";#theta [degree]; #theta - #theta_{fit} [mrad]");
  // // gadiff_pi->SetMarkerStyle(20);
  // // gadiff_pi->SetMarkerSize(0.8);
  // // gadiff_pi->SetLineColor(kRed+1);
  // // gadiff_pi->SetMarkerColor(kRed+1);
  // // gphi_pi->SetTitle(";#theta [degree]; #varphi_{fit} [mrad]");
  // // gphi_pi->SetMarkerStyle(20);
  // // gphi_pi->SetMarkerSize(0.8);
  // // gphi_pi->SetLineColor(kRed+1);
  // // gphi_pi->SetMarkerColor(kRed+1);


  // TLegend *leg_p = new TLegend(0.2,0.7,0.5,0.87);
  // leg_p->SetFillColor(0);
  // leg_p->SetFillStyle(0);
  // leg_p->SetBorderSize(0);
  // leg_p->SetFillStyle(0);
  // //leg_p->AddEntry(gphi,"protons ","lp");
  // //leg_p->AddEntry(gphi_pi,"pions","lp");
  
  // // TLegend *leg_pi = new TLegend(0.6,0.7,0.9,0.87);
  // // leg_pi->SetFillColor(0);
  // // leg_pi->SetFillStyle(0);
  // // leg_pi->SetBorderSize(0);
  // // leg_pi->SetFillStyle(0);
  // // leg_pi->AddEntry(gadiff,"protons ","lp");
  // // leg_pi->AddEntry(gadiff_pi,"pions","lp");

  // prt_canvasAdd("adiff",800,400);
  // gadiff->Draw("APL");
  // //gadiff_pi->Draw("same PL");
  // leg_p->Draw();
  // prt_canvasAdd("phi",800,400);
  // gphi->Draw("APL");
  // // gphi_pi->Draw("same PL");
  // // leg_pi->Draw();
  
  prt_canvasSave("data/drawMapSep_401_nocorr",0);
  
}
