#include "draw_scan.C"
#if defined(__ACLIC__)
#include "PrtTools.h"
#else
R__LOAD_LIBRARY(../prtdirc/build/libPrt.so)
#endif

void wait(int n, int id) {
  TString c = Form("ps aux | grep \"[r]oot.exe\\|[p]rtdirc\" | grep %d | wc -l", id);
  while (true) {
    int active = atoi(gSystem->GetFromPipe(c)) - 1;
    if (active <= n) break;
    else gSystem->Sleep(3000);
    std::cout << "A " << active << std::endl;
  }
}

// level
// 0 simulation
// 1 create lut
// 2 reconstruction

// mc
// 0 process beam data only
// 1 process simulation data only
// 2 process both

void loop(int study = 409, int level = 0, int mc = 2, int fid = -1) {
  PrtTools t(1);
  int threads = 25;
  int events = 20000;
  int nfiles = 0;
  int eventsj = events / threads;

  TString path = Form("/d/proc/jul18/%d/", study);
  if (study < 400) path = Form("$HOME/data/aug17/%d/", study);
  TString exe = "../prtdirc/build/prtdirc ";

  TString smc[2] = {"C", "S"};
  TString lut, sim, simj, rec, pdf;

  for (int imc = 0; imc < 2; imc++) {
    if (imc == 0 && level < 2) continue;
    if (mc == 0 && imc == 1) continue;
    if (mc == 1 && imc == 0) continue;

    TString  rlist = "";
    for (auto run : t.get_runs(study)) {
      int id = run->getId();
      if (fid > -1 && id != fid) continue;      
      TString nid = path + run->getName() + smc[imc];
      TString end = Form("-b 1 -v 1 > %s.%d.log", nid.Data(), level);
      
      int theta = run->getTheta() + 0.1;
      // if (theta % 5 == 0 && theta % 10 != 0) continue;

      sim = Form("-r 0 -o %s.root -study %d -fid %d -e %d ", nid.Data(), study, id, events);
      lut = Form("-r 1 -o %s.lut.root -study %d -fid %d -e 10000000 ", nid.Data(), study, id);
      rec = Form("-r 2 -i %s.root -o %s.rec.root -e 1500 -tr 0.5 ", nid.Data(), nid.Data());
      pdf = Form("-r 4 -i %s.root -o %s.rec.root -e 1500 -tr 0.5 ", nid.Data(), nid.Data());

      rlist += nid + ".rec.root ";

      lut += end + Form(" && cd ~/dirc/prtdirc/macro > /dev/null && root -q -b loadlib.C "
                        "lutmean_cs.C'(\"%s.lut.root\")'",
                        nid.Data());

      if (level == 0) {
        for (int i = 0; i < threads; i++) {
          simj = sim + Form(" -o %sJ%d.root -e %d -seed %d ", nid.Data(), i, eventsj, i);
          if (level == 0) gSystem->Exec(exe + simj + end + "&");
        }

        wait(0,study);
        gSystem->Exec(Form("hadd -f %s.root %sJ*.root ", nid.Data(), nid.Data()));
        gSystem->Exec(Form("rm %sJ*.root ", nid.Data()));
      }
      if (level == 1) gSystem->Exec(exe + lut + "&");
      if (level == 2) gSystem->Exec(exe + rec + end + "&");
      if (level == 4) gSystem->Exec(exe + pdf + end + "&");
    }

    wait(0,study);

    if (level == 2) {
      gSystem->Exec(Form("hadd -f %srec_%d%s.root ", path.Data(), study, smc[imc].Data()) + rlist);
    }
  }
}
