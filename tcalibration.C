// tcalibration - routine for the prtdirc data calibration
// original author: Roman Dzhygadlo - GSI Darmstad

#define TTSelector_cxx
#include "tcalibration.h"

TString ginFile(""), goutFile(""), gcFile("");

const int maxch = 2000;
int gTrigger(0), gSetup(2019), gMode(0), gComboId(0), gMaxIn[maxch];
double tdcRefTime[maxch], gTotO[maxch], gTotP[maxch][10], gLeOffArr[maxch],
  gEvtOffset(0), gPilasOffset[maxch];
TGraph *gGrIn[maxch], *gWalk[maxch], *gGrDiff[maxch];

double walktheta(-5 * TMath::Pi() / 180.);

double tof1le(0), tof2le(0), tof1tot(0), tof2tot(0);
double fr11[11] = {0, 0.5, 0.5, 0.3, 0.3, 0.4, 0.3, 0.3, 0.2, 0.20, 0.15};
double fr12[11] = {0, 1.0, 1.0, 0.9, 0.9, 0.9, 0.9, 0.9, 0.8, 0.80, 0.70};
double fr21[11] = {0, 0.8, 0.8, 0.3, 0.3, 0.4, 0.3, 0.3, 0.2, 0.2, 0.2};
double fr22[11] = {0, 1.0, 1.0, 0.9, 0.9, 0.9, 0.9, 0.9, 0.8, 0.8, 0.8};

double tof1lea[] = {0, 0, 81.84, 76.48, 74.51, 73.58, 73.06, 31.74, 71.91};
double tof2lea[] = {0, 0, 72.12, 72.03, 71.99, 71.96, 71.94, 32.58, 72.55};
double tof1tota[] = {0, 0, 45.14, 45.11, 45.05, 45.03, 44.97, 47.24, 44.91};
double tof2tota[] = {0, 0, 45.26, 45.31, 45.32, 45.34, 45.36, 47.09, 45.38};

// //aug 2017
// double tofpi1[]={0,0,71.50,31.50,0.00, 31.00,31.00, 31.00,31.00, 0, 31.0};
// double tofpi2[]={0,0,72.50,32.20,0.00, 32.20,32.00, 31.95,31.80, 0 ,31.8};

// double tofp1[] ={0,0,80.80,35.70,0.00, 33.00,32.60, 32.35,32.40, 0, 32.2};
// double tofp2[] ={0,0,83.20,37.00,0.00, 34.00,33.50, 33.50,33.00, 0, 33.0};

// jul 2018
double tofpi1a[] = {0, 0, 71.50, 31.50, 0.00, 31.00, 31.00, 32.00, 31.00, 0, 31.0};
double tofpi2a[] = {0, 0, 72.50, 32.20, 0.00, 32.20, 32.00, 33.50, 31.80, 0, 31.8};

double tofp1a[] = {0, 0, 80.80, 35.70, 0.00, 33.00, 32.60, 33.90, 32.40, 0, 32.2};
double tofp2a[] = {0, 0, 83.20, 37.00, 0.00, 34.00, 33.50, 35.50, 33.00, 0, 33.0};

double tofpi1[] = {0, 0, 35.00, 35.00, 35.00, 35.00, 35.00, 35.00, 35.00, 35.00, 35.0};
double tofpi2[] = {0, 0, 36.60, 36.60, 36.20, 36.30, 36.20, 36.00, 35.80, 35.80, 35.6};

double tofp1[] = {0, 0, 44.20, 39.50, 37.80, 37.00, 36.60, 36.45, 36.30, 36.25, 36.2};
double tofp2[] = {0, 0, 47.50, 41.50, 39.20, 38.20, 37.60, 37.20, 37.00, 37.00, 36.7};

int gg_nevents(0);

bool TTSelector::isPion(double tof, int mom) {
  // std::cout<<"tof1le "<<tof1le<<std::endl;
  // std::cout<<"tof2le "<<tof2le<<std::endl;

  return tof1le - 5 * 0.2 < tof && tof < 0.5 * (tof1le + tof2le);
}

bool TTSelector::isProton(double tof, int mom) {
  return 0.5 * (tof1le + tof2le) < tof && tof < tof2le + 5 * 0.2;
}

double TTSelector::getTotWalk(double tot, int ch, int type) {
  double minp(0), walk(0), d(0), min(100);

  if (type == 0) {
    if (ch < t.maxdircch()) {
      for (int i = 0; i < 9; i++) {
        if (gTotP[ch][i] < 0.00000001) continue;
        d = gTotP[ch][i] - tot;
        if (fabs(d) < fabs(min)) {
          minp = gTotP[ch][i];
          min = d;
        }
      }
    }
    double wcorr(10);
    if (ch / 48 == 1) wcorr = 5;
    if (ch / 48 == 3) wcorr = 15;
    if (ch / 48 == 5) wcorr = 12;
    // if(ch/48==9) wcorr=0;

    if (fabs(min) < 0.8) walk = -min * tan(wcorr * TMath::Pi() / 180.);
    if (tot < 8) walk -= (4 - tot) * tan(10 * TMath::Pi() / 180.);
  }

  if (type == 1) { // walk of the xxx
    walk += (38.85 - tot) * tan(25 * TMath::Pi() / 180.);
  }

  return walk;
}

void TTSelector::Begin(TTree *) {
  TString option = GetOption();
  TObjArray *strobj = option.Tokenize(" ");
  gTrigger = ((TObjString *)strobj->At(0))->GetString().Atoi();
  gMode = ((TObjString *)strobj->At(1))->GetString().Atoi();
  gSetup = ((TObjString *)strobj->At(2))->GetString().Atoi();

  t.read_db("data_db.dat");
  run = t.find_run(ginFile);
  std::cout << "gSetup " << gSetup << std::endl;
  
  t.create_maps(gSetup);

  int study = run->getStudy();
  int fileid = run->getId();

  int radiator = run->getRadiator();
  double radiatorL = (radiator == 2) ? 1224.9 : 1200; // plate : bar
  double radiatorW = (radiator == 2) ? 174.8 : 34.9;  // plate : bar
  double radiatorH = (radiator == 2) ? 1224.9 : 17.1; // plate : bar
  run->setRadiatorL(radiatorL);
  run->setRadiatorW(radiatorW);
  run->setRadiatorH(radiatorH);

  TString filedir = ginFile;
  filedir.Remove(filedir.Last('.') - 4);
  fFile = new TFile(goutFile, "RECREATE");
  fTree = new TTree("data", "Tree for GSI Prt Analysis");
  fEvent = new PrtEvent();
  fTree->Branch("PrtEvent", "PrtEvent", &fEvent, 64000, 2);

  fRunTree = new TTree("header", "run info");
  fRunTree->Branch("PrtRun", "PrtRun", &run, 64000, 2);
  fRunTree->Fill();
  
  if (gcFile != "0") {
    TFile f(gcFile);
    TIter nextkey(f.GetListOfKeys());
    TKey *key;

    while ((key = (TKey *)nextkey())) {
      TGraph *gr = (TGraph *)key->ReadObj();
      TString name = gr->GetName();

      double x, y;
      if (name.Contains("tof_")) {
        name.Remove(0, 4);

        if (ginFile.Contains(name)) {
          gr->GetPoint(0, tof1le, tof2le);
          gr->GetPoint(1, tof1tot, tof2tot);
        }
        continue;
      }
      if (name.Contains("off_")) { // read event offsets
        name.Remove(0, 4);
        if (ginFile.Contains(name)) gr->GetPoint(0, x, gEvtOffset);
        continue;
      }

      if (name.Contains("walk_")) { // read walk corrections
        name.Remove(0, 5);
        int ch = name.Atoi();
        gWalk[ch] = new TGraph(*gr);
        continue;
      }

      if (name.Contains("evl")) {
        for (int i = 0; i < t.maxdircch(); i++) {
          gr->GetPoint(i, x, gPilasOffset[i]);
        }
        continue;
      }

      long long ch = name.Atoll();
      if (ch < 10000) { // spline calibration
        gGrIn[ch] = new TGraph(*gr);
      } else if (ch == 10000) { // line calibration
        for (int i = 0; i < maxch; i++) {
          gr->GetPoint(i, x, y);
          gMaxIn[i] = (int)(y + 0.01);
          // std::cout<<"ch  "<<i<< "  FT max"<<  gMaxIn[i]<<std::endl;
        }
      } else if (ch == 10001) { // read tot offsets
        for (int i = 0; i < maxch; i++) {
          gr->GetPoint(i, gTotO[i], y);
          // std::cout<<"ch  "<<i<< " TOT off "<<  gTotO[i]<<std::endl;
        }
      } else if (ch == 10002) { // read tot peaks
        for (int i = 0; i < t.maxdircch() * 10; i++) {
          gr->GetPoint(i, x, y);
          gTotP[i / 10][i % 10] = y;
          // std::cout<<"ch  "<<i/10<< " peak "<< i%10<< " = " <<y<<std::endl;
        }
      } else if (ch == 10003) { // read LE offsets 1
	std::cout<<"t.maxdircch() "<<t.maxdircch()<<std::endl;
	
        for (int i = 0; i < t.maxdircch(); i++) {
          gr->GetPoint(i, gLeOffArr[i], y);
        }
      }
    }
    f.Close();
  }

  
  std::cout << "study id " << study << " file id " << fileid << std::endl;

  std::cout << "Initialization successful" << std::endl;
}

bool TTSelector::Process(Long64_t entry) {

  int tdc, ch, tofpid(0);
  double grTime0(0), grTime1(0), grTime2(0), coarseTime(0), offset(0), triggerLe(0), triggerTot(0);
  double time[10000], timeLe(0), timeT[10000], timeTot(0), mom(7), simOffset(74.20);
  int multT1(0), multT2(0), multT3v(0), multT3h(0), multTof1(0), multTof2(0), multStr1(0),
    multStl1(0), multStr2(0), multStl2(0);

  TString current_file_name = TTSelector::fChain->GetCurrentFile()->GetName();
  bool trbdata = current_file_name.Contains("trb");
  bool laser = current_file_name.Contains("pilas") || current_file_name.Contains("pico");

  TObjArray *sarr = current_file_name.Tokenize("_");

  if (entry % 10000 == 0) std::cout << "event # " << entry << std::endl;
  GetEntry(entry);

  int trigT1(520);
  int trigT2(513);
  int trigT3h(514);
  int trigT3v(515);
  int trigTof1(1136);
  int trigTof2(1138);

  int trigStr1(1140);
  int trigStl1(1142);
  int trigStr2(1144);
  int trigStl2(1146);

  if (gSetup == 2017) {
    trigT1 = 816;
    trigT2 = 817;
    trigT3h = 818;
    trigT3v = 819;
    trigTof1 = 1392;
    trigTof2 = 1398;
    simOffset += 15;
  }

  fEvent = new PrtEvent();
  if (gMode == 5) {
    int study = run->getStudy();
    if (study == 401) gTrigger = 1138;
    if (study == 403) gTrigger = 1136;
    if (study > 0) {
      mom = run->getMomentum();
      fEvent->setMomentum(TVector3(0, 0, mom));

      // if(fEvent->GetLens()==6) lenw = 12;
      // if(fEvent->GetLens()==3) lenw = 15;
      // if(fEvent->GetLens()==2) lenw = 14.4;
      double zrot = 155, xrot = 98;

      if (study < 400) {
	zrot = 146;
	xrot = 97.8;
      }

      double rad = TMath::Pi() / 180.0, prtangle = run->getTheta(), z = run->getBeamZ(),
             b = xrot * tan(0.5 * (prtangle - 90) * rad),
             lenz = (z - zrot + b) / cos((prtangle - 90) * rad) + b + zrot;

      fEvent->setPosition(TVector3(0, 0, lenz));
    }
  }

  int tdc_trig = t.get_tdcid(gTrigger);

  for (int i = 0; i < Hits_ && i < 10000; i++) {
    
    tdc = t.map_tdc[Hits_nTrbAddress[i]];
    if (tdc < 0) continue;
    ch = t.get_channel(tdc, Hits_nTdcChannel[i]) - 1;
    
    time[i] = 5 * (Hits_nEpochCounter[i] * pow(2.0, 11) + Hits_nCoarseTime[i]); // coarsetime    
    if (gcFile != "0") {
      int rch = ch + 1 + tdc;
      
      // spline calib
      // time[i] -= gGrIn[rch]->Eval(Hits_nFineTime[i]+1); //slow
      double xx, yy;      
      gGrIn[rch]->GetPoint(Hits_nFineTime[i], xx, yy);
      time[i] -= yy; // fast

      // linear calib
      // double max = (double) gMaxIn[rch]-2;
      // time[i] = coarseTime-5*(Hits_nFineTime[i]-31)/(max-31);
    } // else time[i] -= (Hits_nFineTime[i]-31)*0.0102;

    if (Hits_nSignalEdge[i] == 1) {
      if (ch == gTrigger && grTime1 == 0) grTime1 = time[i];
      if (Hits_nTdcChannel[i] == 0) { // ref channel
        tdcRefTime[tdc] = time[i];
        if (tdc_trig == tdc) grTime0 = time[i];
      }
      if (ch == trigT1) multT1++;     // trigger1
      if (ch == trigT2) multT2++;     // trigger2
      if (ch == trigTof1) multTof1++; // tof1
      if (ch == trigTof2) multTof2++; // tof2
      if (ch == trigT3h) multT3h++;   // trigger3h
      if (ch == trigT3v) multT3v++;   // trigger3v

      if (ch == trigStr1) multStr1++;
      if (ch == trigStl1) multStl1++;
      if (ch == trigStr2) multStr2++;
      if (ch == trigStl2) multStl2++;
    } else {
      timeT[i] = time[i];
      if (ch == gTrigger && grTime2 == 0) grTime2 = time[i];
    }
  }

  double tof1(0), tof2(0), tot1(0), tot2(0), toftime(0), mass(0);
  double tofstr1(0), tofstr2(0), totstr1(0), totstr2(0), tofstl1(0), tofstl2(0), totstl1(0),
    totstl2(0);

  if (gMode == 5) {
    if (gSetup == 2018)
      if (multT1 < 1 || multTof1 < 1 || multTof2 < 1 || multT3h < 1 ||
          multT3v < 1) { //  || mult2!=1 || mult5!=1
        // if(multT1<1 || multTof1<1 || multTof2<1 || multStr1<1 ||  multStl1<1|| multStr2<1 ||
        // multStl2<1 || multT3h<1 || multT3v<1){ if(multT1<1 || multTof1<1 || multTof2<1){ //  ||
        // mult2!=1 || mult5!=1
        fEvent->Clear();
        delete fEvent;
        return kTRUE;
      }

    for (int i = 0; i < Hits_ && i < 10000; i++) {
      if (Hits_nTdcErrCode[i] != 0) continue;
      if (Hits_nTdcChannel[i] == 0) continue; // ref channel
      if (Hits_nSignalEdge[i] == 0) continue; // tailing edge

      tdc = t.map_tdc[Hits_nTrbAddress[i]];
      ch = t.get_channel(tdc, Hits_nTdcChannel[i]) - 1;
      if (ch == trigTof1 && tof1 == 0) {
        tof1 = time[i] - tdcRefTime[tdc];
        tot1 = timeT[i + 1] - time[i];
      }
      if (ch == trigTof2 && tof2 == 0) {
        tof2 = time[i] - tdcRefTime[tdc];
        tot2 = timeT[i + 1] - time[i];
      }
      if (ch == trigStr1 && tofstr1 == 0) {
        tofstr1 = time[i] - tdcRefTime[tdc];
        totstr1 = timeT[i + 1] - time[i];
      }
      if (ch == trigStr2 && tofstr2 == 0) {
        tofstr2 = time[i] - tdcRefTime[tdc];
        totstr2 = timeT[i + 1] - time[i];
      }
      if (ch == trigStl1 && tofstl1 == 0) {
        tofstl1 = time[i] - tdcRefTime[tdc];
        totstl1 = timeT[i + 1] - time[i];
      }
      if (ch == trigStl2 && tofstl2 == 0) {
        tofstl2 = time[i] - tdcRefTime[tdc];
        totstl2 = timeT[i + 1] - time[i];
      }
    }

    if (tof1 != 0 && tof2 != 0) {
      double time = tof2 - tof1;
      // time += (tot1-tof1tot)*tan(walktheta);
      // time += (tot2-tof2tot)*tan(-walktheta);

      if (tof1tot < 55 && tof1tot > 10) time += (tot1 - tof1tot) * tan(-3 * TMath::DegToRad());
      if (tof2tot < 55 && tof2tot > 10) time += (tot2 - tof2tot) * tan(2 * TMath::DegToRad());

      // time += (tot1-41.32)*tan(-4*TMath::Pi()/180.);
      // time += (tot2-40.75)*tan(2*TMath::Pi()/180.);

      toftime = time;
      int m = (double)(mom + 0.1);
      if (isPion(time, m)) {
        tofpid = 2; // 211;
        mass = 0.13957018;
      } else if (isProton(time, m)) {
        tofpid = 4; // 2212;
        mass = 0.938272046;
      } else {
        fEvent->Clear();
        delete fEvent;
        return kTRUE;
      }
    }

    // if(tofstr1!=0 && tofstr2!=0 && tofstl1!=0 && tofstl2!=0 ) {
    // tofpi1[7]=32.0;
    // tofpi2[7]=33.3;
    // tofp1[7] =33.7;
    // tofp2[7] =35;
    // double time =(tofstr2+tofstl2)/2.-(tofstr1+tofstl1)/2.;

    // toftime = time;
    // int m = (double) (mom+0.1);

    // if(isPion(time,m)){
    // 	tofpid=2;
    // 	mass=0.13957018;
    // }else if(isProton(time,m)){
    // 	tofpid=4;
    // 	mass = 0.938272046;
    // }else{
    // 	fEvent->Clear();
    // 	delete fEvent;
    // 	return kTRUE;
    // }
    // }
  }

  PrtHit hit;
  int nrhits = 0;

  if ((grTime0 > 0 && grTime1 > 0) || gTrigger == 0) {
    if (gTrigger != 0) {
      triggerLe = grTime1 - grTime0;
      if (gTrigger == 1140) triggerLe = (tofstr1 + tofstl1) / 2.;
      if (gTrigger == 1144) triggerLe = (tofstr2 + tofstl2) / 2.;
      triggerTot = grTime2 - grTime1;
    }

    for (int i = 0; i < Hits_ && i < 10000; i++) {
      // if(Hits_nTdcErrCode[i]!=0) continue;
      if (Hits_nTdcChannel[i] == 0) continue; // ref channel
      if (Hits_nSignalEdge[i] == 0) continue; // tailing edge

      tdc = t.map_tdc[Hits_nTrbAddress[i]];
      ch = t.get_channel(tdc, Hits_nTdcChannel[i]) - 1;

      // if(!trbdata && prt_isBadChannel(ch)) continue;

      if (gMode > 0) {
        timeLe = time[i] - tdcRefTime[tdc];
        if (gTrigger != 0 && ch < t.maxdircch()) timeLe = timeLe - triggerLe;
      } else {
        timeLe = time[i];
        if (gTrigger != 0 && ch < t.maxdircch()) timeLe = timeLe - grTime1;
      }

      timeTot = timeT[i + 1] - time[i];
      if (ch < t.maxdircch()) {
        timeTot += 30 - gTotO[ch];

        // timeLe += getTotWalk(timeTot,ch);
        // if(gTrigger==trigT1 && fabs(triggerTot-tof1tot)<1) timeLe -=
        // (triggerTot-tof1tot)*tan(5*TMath::Pi()/180.);

        if (gTrigger == trigTof2 && fabs(tot2 - tof2tot) < 1.2)
          timeLe += (tot2 - tof2tot) * tan(2 * TMath::Pi() / 180.); // 7.1;
        if (timeTot > 0.5 && timeTot < 9 && gWalk[ch]) timeLe -= gWalk[ch]->Eval(timeTot);
        timeLe -= (gLeOffArr[ch] - gPilasOffset[ch]);

        if (!laser && gMode == 5) {
          double rad = TMath::Pi() / 180., zrot = 155, xrot = 98, prtangle = run->getTheta(),
                 z = run->getBeamZ(),
                 rot_dist =
                   ((z - zrot) - xrot / cos(prtangle * rad)) * tan((90 - prtangle) * rad) / 1000.;

          if (gTrigger == trigT1)
            timeLe -= (2.719 + rot_dist) / ((mom / sqrt(mass * mass + mom * mom) * 299792458)) *
                      1E9; // trig1
          if (gTrigger == trigTof1)
            timeLe -= (24.490 + rot_dist) / ((mom / sqrt(mass * mass + mom * mom) * 299792458)) *
                      1E9; // tof1
          if (gTrigger == 1140)
            timeLe -= (24.490 + 0.030 + rot_dist) /
                      ((mom / sqrt(mass * mass + mom * mom) * 299792458)) * 1E9; // tof1
          if (gTrigger == trigTof2)
            timeLe += (4.151 - rot_dist) / ((mom / sqrt(mass * mass + mom * mom) * 299792458)) *
                      1E9; // tof2
          if (gTrigger == 1144)
            timeLe += (4.151 - rot_dist) / ((mom / sqrt(mass * mass + mom * mom) * 299792458)) *
                      1E9; // scitil2

          timeLe += simOffset;
          if (gTrigger == trigTof2) timeLe -= 59;
        }
        if (!laser && gMode != 5) timeLe += -2.5; // simOffset;
      }

      if (gMode == 5) {
        if (fabs(timeLe) > 600) continue;
        // if(ch==trigTof1) timeLe -= (tot1-tof1tot)*tan(-3*TMath::Pi()/180.);
        // if(ch==trigTof2) timeLe -= (tot2-tof2tot)*tan(2*TMath::Pi()/180.);

        timeLe -= gEvtOffset;

        // if(ch>820 && ch<1340) continue;
        // if(ch<t.maxdircch() && (timeLe<0 || timeLe>100)) continue;
      }

      if (gMode != 5 || tofpid != 0) {

        if (run->getGeometry() == 2023) {
          if (timeTot < 0 || timeTot > 20) continue;
          if (timeLe < -100 || timeLe > 100) continue;
        }

        hit.setChannel(ch);
        hit.setPmt(t.map_pmt[ch]);
        hit.setPixel(t.map_pix[ch]);
        hit.setLeadTime(timeLe);
        hit.setTotTime(timeTot);
        fEvent->addHit(hit);
        nrhits++;
      }
    }
  }

  if (nrhits != 0) {
    gg_nevents++;
    fEvent->setPid(tofpid);
    fEvent->setTof(toftime);
    fEvent->setTofPi(tof1le);
    fEvent->setTofP(tof2le);
    fTree->Fill();
  }
  fEvent->Clear();
  delete fEvent;

  return kTRUE;
}

void TTSelector::Terminate() {
  fFile->Write();
  fFile->Close();
}

void tcalibration(TString inFile = "../../data/cj.hld.root", TString outFile = "outFileC.root",
                  TString cFile = "calib.root", TString tFile = "calibOffsets.root",
                  int trigger = 0, int sEvent = 0, int eEvent = 0, int mode = 5, int build = 0,
                  int setupid = 2018) {
  if (build == 1) return;
  ginFile = inFile;
  goutFile = outFile;
  gcFile = (cFile != "") ? cFile : "0"; // calibration
  gTrigger = trigger;
  gMode = mode;
  gSetup = setupid;
  
  if (gMode == 5) gTrigger = 1136;
  if (setupid == 2017 && gTrigger != 820) gTrigger = 1392;
  TChain *ch = new TChain("T");
  ch->Add(ginFile);

  int entries = ch->GetEntries();
  TTSelector *selector = new TTSelector();
  TString option = Form("%d %d %d", gTrigger, gMode, gSetup);

  if (eEvent == 0) {
    std::cout << "Entries in chain:  " << entries << std::endl;
    ch->Process(selector, option, entries);
  } else {
    ch->Process(selector, option, eEvent - sEvent, sEvent);
  }
}
