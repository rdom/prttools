#define prt__beam
#include "../prtdirc/src/PrtHit.h"
#include "../prtdirc/src/PrtEvent.h"
#include "prttools.C"

void procOffsets(TString path="",int corrected=1){
  
  if(path=="") return;
  TString fileid(path);
  fileid.Remove(0,fileid.Last('/')+1);
  fileid.Remove(fileid.Last('.')-1);
  
  prt_data_info = getDataInfo(fileid);
  int studyid = prt_data_info.getStudyId();
 
  int h1a(0),h1b(50),h2a(0),h2b(50),hbin(1000);//h1a(200),h1b(400)

  if(corrected==1){
    h1a=0;
    h1b=50;
    h2a=0;
    h2b=50;
    hbin=1000;
  }
  
  TString outdir=path;outdir.Remove(outdir.Last('/'));
  TString sstudy=outdir; sstudy.Remove(0,sstudy.Last('/'));  
  TString insim = path; insim.ReplaceAll("C.root","S.root");
  
  if(!prt_init(insim,1,outdir+Form("/%da/%d",studyid,prt_data_info.getFileId()))) return;
  prt_ch->Add(path);

  TH1F * hLeD  = new TH1F("leD","LE beam data ; LE [ns]; entries [#]",hbin,h1a,h1b);
  TH1F * hLeS  = new TH1F("leS","LE simulation; LE [ns]; entries [#]",hbin,h2a,h2b);
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);

  int maxent(0);
  for (int ievent=0; ievent< prt_ch->GetEntries(); ievent++){
    if(maxent>20000) continue;
    prt_nextEvent(ievent,10000);

    if(prt_event->GetType()==0){
      
      int gch, ndirc(0), t2(0), t3h(0), t3v(0);
      int hodo1(0), hodo2(0);
      for(auto h=0; h<prt_event->GetHitSize(); h++){
	gch = prt_event->GetHit(h).GetChannel();
      
	if(gch<prt_maxdircch) ndirc++;
      
	if(gch==818) t3h++;
	if(gch==819) t3v++;
	if(gch>=1350 && gch<=1352) hodo1++;
	if(gch>=1367 && gch<=1372) hodo2++;
      }
      // if(ndirc<5) continue;
      // if(!(t3h && t3v)) continue;
      // if(!(t3h && t3v && hodo1 && hodo2)) continue;
    }
    
    if(prt_event->GetParticle()!=211) continue;
    bool bsim(false);
    TString current_file_name  = prt_ch->GetCurrentFile()->GetName();
    if(current_file_name.Contains("S.root")) bsim = true;
    else maxent++;
    
    for(auto hit : prt_event->GetHits()){
      int ch = hit.GetChannel();
      if(ch<512 ){
	double time = hit.GetLeadTime();
	if(bsim){
	  time += prt_rand.Gaus(0,0.38);
	  hLeS->Fill(time);
	}else{
	  if(studyid == 401){
	    double o, prtangle = prt_theta;
	    if(fabs(prtangle-20)<1) o = -5.65;
	    if(fabs(prtangle-25)<1) o = -5.65;
	    if(fabs(prtangle-30)<1) o = -5.65;
	    if(fabs(prtangle-60)<1) o = -4.3;
	    if(fabs(prtangle-90)<1) o = -5.4;
	    if(fabs(prtangle-140)<1) o = -4.1;

	    
	    time += o;
	  }
	  hLeD->Fill(time);	  
	}
      }
    }
  }
  
  prt_canvasAdd("offset",800,400);
  
  hLeS->SetLineColor(kRed);
  TH1F *hh[] = {hLeS,hLeD};
  prt_normalizeto(hh,2,1);
  hLeS->Draw("hist");
  hLeD->Draw("hist same");
  TLegend *leg = new TLegend(0.62,0.7,0.92,0.9);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(hLeD,"beam data ","lp");
  leg->AddEntry(hLeS,"simulation","lp");
  leg->Draw();    

  if(corrected==0){
    // double xmax1 = hLeD->GetXaxis()->GetBinCenter(hLeD->GetMaximumBin());
    // double xmax2 = hLeS->GetXaxis()->GetBinCenter(hLeS->GetMaximumBin());
    // double xmax1 = prt_fit(hLeD,0.5,50,1).X();
    // double xmax2 = prt_fit(hLeS,0.5,50,1).X();

    double threshold =  hLeD->GetMaximum()*0.3;      
    int firstbin = hLeD->FindFirstBinAbove(threshold);
    double xmax1 = hLeD->GetXaxis()->GetBinCenter(firstbin);
    threshold =  hLeS->GetMaximum()*0.3;      
    firstbin = hLeS->FindFirstBinAbove(threshold);
    double xmax2 = hLeS->GetXaxis()->GetBinCenter(firstbin);

    gPad->Update();
    TLine *gLine1 = new TLine(0,0,0,1000);
    gLine1->SetX1(xmax1);
    gLine1->SetX2(xmax1);
    gLine1->SetY1(gPad->GetUymin());
    gLine1->SetY2(gPad->GetUymax());
    gLine1->SetLineColor(kBlue);
    gLine1->Draw();

    TLine *gLine2 = new TLine(0,0,0,1000);
    gLine2->SetX1(xmax2);
    gLine2->SetX2(xmax2);
    gLine2->SetY1(gPad->GetUymin());
    gLine2->SetY2(gPad->GetUymax());
    gLine2->SetLineColor(kRed);
    gLine2->Draw();
        
    TFile efile(path+ ".off.root","RECREATE");
    TGraph *gr = new TGraph();
    gr->SetPoint(0,xmax1-xmax2,  xmax1-xmax2);
    fileid=path;
    fileid.Remove(0,fileid.Last('_')+1);
    fileid.Remove(fileid.Last('C'));
    gr->SetName("off_"+fileid);    
    gr->Write();

    gr = new TGraph();
    gr->SetPoint(0,prt_theta, xmax1);
    gr->SetName("offbeam");    
    gr->Write();
    
    
    gr = new TGraph();
    gr->SetPoint(0,prt_theta, xmax2);
    gr->SetName("offsim");    
    gr->Write();

    gr = new TGraph();
    gr->SetPoint(0,prt_theta, xmax1-xmax2);
    gr->SetName("offdiff");    
    gr->Write();
    
    efile.Write();
    efile.Close();
    std::cout<<"new offset "<< xmax1-xmax2 <<std::endl;    
    
  }

  prt_canvasSave(1,0);
}
