void loadprtlibs(){
  gROOT->ProcessLine(".L ../prtdirc/src/PrtHit.cxx+");
  gROOT->ProcessLine(".L ../prtdirc/src/PrtEvent.cxx+");

  gSystem->Load("../prtdirc/src/PrtHit_cxx.so");
  gSystem->Load("../prtdirc/src/PrtEvent_cxx.so");
}
