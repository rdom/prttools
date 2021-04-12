#include "PrtTools.h"

void wait(int n) {
  TString c = "ps aux | grep \"[r]oot.exe\\|[p]rtdirc\" | wc -l";
  while (true) {
    gSystem->Sleep(1000);
    int active = atoi(gSystem->GetFromPipe(c)) - 1;
    if (active <= n) break;
    std::cout << "A " << active << std::endl;
  }
}

void loop(int study = 409, int level = 0, int fid = -1) {
  PrtTools t;
  int threads = 25;
  int events = 45000;
  int nfiles = 0;
  int eventsj = events / threads;

  TString path = Form("$HOME/data/jul18/%d/", study);
  TString exe = "../prtdirc/build/prtdirc ";
  TString lut, sim, simj, rec, rlist = "";

  for (auto run : t.get_runs(study)) {
    TString nid = path + run->getName();
    int id = run->getId();
    TString end = Form("-b 1 -v 1 > %s.%d.log", nid.Data(),level);

    sim = Form("-r 0 -o %sS.root -study %d -fid %d -e %d ", nid.Data(), study, id, events);
    lut = Form("-r 1 -o %sS.lut.root -study %d -fid %d -e 10000000 ", nid.Data(), study, id);
    rec = Form("-r 2 -i %sS.root -o %s.rec.root -u %s.lut.cs_avr.root -e 2000 -tr 0.5 ", nid.Data(),
               nid.Data(), nid.Data());

    rlist += nid + ".rec.root ";

    lut += end + Form(" && cd ~/dirc/prtdirc/macro > /dev/null && root -q -b loadlib.C "
                      "lutmean_cs.C'(\"%sS.lut.root\")'",
                      nid.Data());

    if (level == 0) {
      for (int i = 0; i < threads; i++) {
        simj = sim + Form(" -o %sSJ%d.root -e %d -seed %d ", nid.Data(), i, eventsj, i);
        if (level == 0) gSystem->Exec(exe + simj + end + "&");
      }

      wait(0);
      gSystem->Exec(Form("hadd -f %sS.root %sSJ*.root ", nid.Data(), nid.Data()));
      gSystem->Exec(Form("rm %sSJ*.root ", nid.Data()));
    }
    if (level == 1) gSystem->Exec(exe + lut + "&");
    if (level == 2) gSystem->Exec(exe + rec + end + "&");
  }

  wait(0);

  if (level == 2) {
    gSystem->Exec(Form("hadd -f %srec_%d.root ", path.Data(), study) + rlist);
  }
}
