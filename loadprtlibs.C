void loadprtlibs(TString path=""){
  gROOT->ProcessLine(".L "+path+"../prtdirc/src/PrtHit.cxx+");
  gROOT->ProcessLine(".L "+path+"../prtdirc/src/PrtEvent.cxx+");

  gSystem->Load(path+"../prtdirc/src/PrtHit_cxx.so");
  gSystem->Load(path+"../prtdirc/src/PrtEvent_cxx.so");
}
