#include "prttools.C"

const Int_t tdcmax=10000;
Int_t tdcnum=20;
Int_t tdcmap[tdcmax];
TString trbsid[20] ={"2000","2001","2002","2003","2004","2005","2006","2007","2008","2009",
		     "2010","2011","2012","2013","2014","2015",
		     "2016","2017","2018","2019"};


void CreateMap(){
  Int_t seqid =0;
  for(Int_t i=0; i<tdcmax; i++){
    tdcmap[i]=-1;
    for(Int_t j=0; j<tdcnum; j++){
      if(i==TString::BaseConvert(trbsid[j],16,10).Atoi()){
	tdcmap[i]=seqid++;
	break;
      }
    }
  }
}


void show_thresholds(TString infile1="thresholds.thr", TString infile2=""){
  fSavePath = "auto";
  CreateMap();
  initDigi();
  
  ifstream in;
  in.open(infile1.Data());
  TString s, stdc,sth,schain;
  Int_t ch, tdc, th, chain;
  Double_t min(1000000),max(0);
  
  while (1) {
    in >> s >> s >> s >> stdc >> s >> schain >> s >> ch >> s >>sth >> s >> s;
    if (!in.good()) break;
    
    stdc = stdc.Strip(TString::kTrailing,',');
    tdc = TString::BaseConvert(stdc,16,10).Atoi();

    schain = schain.Strip(TString::kLeading,'0');
    schain = schain.Strip(TString::kTrailing,',');
    chain = schain.Atoi();
    
    sth = sth.Strip(TString::kTrailing,',');
    th = TString::BaseConvert(sth,16,10).Atoi();
    if(th>max) max = th;
    if(th<min) min = th;
    
    Int_t tdcSeq = tdcmap[tdc];
    Int_t channel= tdcSeq*48 + chain*16 + ch;
    Int_t mcp = channel/64;
    Int_t pix = channel%64;
    Int_t col = pix/2 - 8*(pix/16);
    Int_t row = pix%2 + 2*(pix/16);
    fhDigi[mcp]->Fill(col,row,th);

    std::cout<<"mcp  "<<mcp << " tdcSeq  " <<tdcSeq << " col  "<< col << "  row "<< row<<std::endl;
  }
  
  drawDigi("m,p,v\n",3,max,min);
    
  TPaletteAxis *pal = new TPaletteAxis(0.90,0.1,0.94,0.90,fhDigi[0]);   
  cDigi->cd();
  pal->Draw();

}
