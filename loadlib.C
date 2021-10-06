void loadlib(TString path=""){
  gROOT->ProcessLine(".L "+path+"../prtdirc/src/PrtHit.cxx+");
  gROOT->ProcessLine(".L "+path+"../prtdirc/src/PrtEvent.cxx+");
  gROOT->ProcessLine(".L "+path+"../prtdirc/src/PrtRun.cxx+");
  gROOT->ProcessLine(".L PrtTools.cxx+");
  

  // gSystem->Load(path+"../prtdirc/src/PrtHit_cxx.so");
  // gSystem->Load(path+"../prtdirc/src/PrtEvent_cxx.so");
}
