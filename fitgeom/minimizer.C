#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TRandom2.h"
#include "TError.h"

#define prt__beam
#include "../../prtdirc/src/PrtHit.h"
#include "../../prtdirc/src/PrtEvent.h"
#include "../prttools.C"

double tpar,chisq;
int status=0, iter;
ifstream gStatusFileR;
ofstream gStatusFileW;

double ratiosim[prt_maxdircch],ratiodat[prt_maxdircch];

void getRatioArray(TString infile="hits.root",bool sim=false){
  if(!prt_init(infile,1,"data/fitgeom")) return;
 
  Double_t mult[5][prt_maxdircch];
  for(auto i=0; i<5; i++){
    for(auto c=0; c<prt_maxdircch; c++){
      mult[i][c]=0;;
    }
  }
  
  PrtHit hit;
  for (auto ievent=0; ievent< prt_entries; ievent++){
    prt_nextEvent(ievent,1000);
    for(auto h=0; h<prt_event->GetHitSize(); h++){
      hit = prt_event->GetHit(h);
      Int_t mcpid = hit.GetMcpId();
      Int_t pixid = hit.GetPixelId()-1;
      Int_t ch = map_mpc[mcpid][pixid];
      if(prt_isBadChannel(ch)) continue;
      mult[prt_pid][ch]++;
    }
  }

  for(auto i=0; i<prt_maxdircch; i++){
    Int_t mcpid = map_mcp[i];
    Int_t pixid = map_pix[i];

    if(sim)ratiosim[i]=0;
    else ratiodat[i]=0;

    if(mult[2][i]<10 || mult[4][i]<10) continue;

    if(sim) ratiosim[i]=mult[2][i]/(double)mult[4][i];
    else ratiodat[i]=mult[2][i]/(double)mult[4][i];
    
    prt_hdigi[mcpid]->Fill(pixid%8, pixid/8,mult[2][i]/(double)mult[4][i]);
  }

  prt_drawDigi("m,p,v\n",2017,4,0);
  if(sim) prt_cdigi->Print(Form("hp_sfit_%d.png",iter));
  else prt_cdigi->Print(Form("hp_dfit_%d.png",iter));
}

double getChiSq(const double *xx){

  iter++;
  TString command = Form("( ./botfit %d %f %f %f  %f %f )", iter, xx[0], xx[1], xx[2], xx[3],chisq);
  gSystem->Exec(command.Data());

  while(status<1){
    gSystem->Sleep(1000);
    gStatusFileR.open("status.dat");
    gStatusFileR >> status >> iter >> tpar >> tpar >> tpar >> tpar >>tpar;
    std::cout<<"status "<<status<<std::endl;
    
    gStatusFileR.close();
  }

  getRatioArray(Form("hits_%d.root",iter));

  double chi;
  chisq = 0;
  for(auto i=0; i<prt_maxdircch; i++){
    chi = fabs(ratiodat[i]-ratiosim[i])/0.1;
    chisq += chi*chi;
  }
  std::cout<<"chisq "<<chisq <<" "<<  xx[0]<<" " <<xx[1]<<" "<< xx[2]<<" "<< xx[3]<<std::endl;
  
  status=0;
  return chisq;
}
 
int minimizer(){
  
  getRatioArray("/d/proc/aug17/332/beam_s332_50S.root",true);

  std::cout<<"test "<<ratiosim[10]<<" "<< ratiosim[100]<<std::endl;
  
  
  iter = 0;
  // algoName Migrad, Simplex,Combined,Scan  (default is Migrad)
  ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit", "Migrad");

  // set tolerance , etc...
  min->SetMaxFunctionCalls(10000); // for Minuit/Minuit2 
  min->SetMaxIterations(10000);  // for GSL 
  min->SetTolerance(0.001);
  min->SetPrintLevel(1);

  // create funciton wrapper for minmizer. a IMultiGenFunction type 

  const int npar=4;
  ROOT::Math::Functor f(&getChiSq,npar); 
  double step[npar] ={ 0.1  0.1,  0.5,  0.5};
  double par[npar] = {50.0, 0.0, 70.0, 16.5};
  min->SetFunction(f);
   
  min->SetVariable(0,"theta",par[0], step[0]);
  min->SetVariable(1,"phi",par[1], step[1]);
  min->SetVariable(2,"gsx",par[2], step[2]);
  min->SetVariable(3,"gsy",par[3], step[3]);

  min->SetVariableLimits(0, par[0]-0.4, par[0]+0.4);
  min->SetVariableLimits(1, par[1]-0.8, par[1]+0.8);
  min->SetVariableLimits(2, par[2]-10, par[2]+10);
  min->SetVariableLimits(3, par[3]-5, par[3]+5);

  min->Minimize(); 
    
  const double *xs = min->X();
  std::cout << "Minimum: f(" << xs[0] << "," << xs[1] << "): " 
	    << min->MinValue()  << std::endl;

  if( min->MinValue() < 8 ) 
    std::cout << "Converged" << std::endl;
  else {
    Error("NumericalMinimization","fail to converge");
  }
 
  return 0;
}
