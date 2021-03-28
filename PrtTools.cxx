#include "PrtTools.h"

PrtTools::PrtTools() {
  init();
}

PrtTools::PrtTools(PrtRun *run) {
  init();
  _npmt = run->getNpmt();
  _npix = run->getNpix();
  _maxdircch = _npmt * _npix;
  _pmtlayout = run->getPmtLayout();
  _info += run->getInfo();
}

void PrtTools::init() {
  _maxch = 5000;
  _npix = 64;
  _npmt = 8;
  _maxdircch = _npmt * _npix;
  _pmtlayout = 2018;
  map_pid = {0};

  _spectrum = new TSpectrum(2);
  _info = "";

  if (gROOT->GetApplication()) {
    TIter next(gROOT->GetApplication()->InputFiles());
    TObjString *os = 0;
    while ((os = (TObjString *)next())) {
      _info += os->GetString() + " ";
    }
    _info += "\n";
  }
}

bool PrtTools::init_run(TString in, int bdigi, TString savepath, int setupid) {

  if (in == "" || gSystem->AccessPathName(in)) {
    std::cout << "file not found " << in << std::endl;
    return false;
  }

  if (savepath != "") _savepath = savepath;
  TGaxis::SetMaxDigits(4);
  set_palette(1);
  create_maps(setupid);
  delete _chain;

  _chain = new TChain("data");
  _chain->Add(in);
  _chain->SetBranchAddress("PrtEvent", &_event);

  _entries = _chain->GetEntries();
  std::cout << "Entries in chain:  " << _entries << std::endl;
  if (bdigi) init_digi();
  return true;
}

void PrtTools::init_digi() {
  int nrow = sqrt(_npix);
  for (int m = 0; m < _npmt; m++) {
    if (_hdigi[m])
      _hdigi[m]->Reset("M");
    else {
      _hdigi[m] = new TH2F(Form("pmt%d", m), Form("pmt%d", m), nrow, 0, nrow, nrow, 0, nrow);
      _hdigi[m]->SetStats(0);
      _hdigi[m]->SetTitle(0);
      _hdigi[m]->GetXaxis()->SetNdivisions(((nrow > 10) ? 20 : 10));
      _hdigi[m]->GetYaxis()->SetNdivisions(((nrow > 10) ? 20 : 10));
      _hdigi[m]->GetXaxis()->SetLabelOffset(100);
      _hdigi[m]->GetYaxis()->SetLabelOffset(100);
      _hdigi[m]->GetXaxis()->SetTickLength(1);
      _hdigi[m]->GetYaxis()->SetTickLength(1);
      _hdigi[m]->GetXaxis()->SetAxisColor(15);
      _hdigi[m]->GetYaxis()->SetAxisColor(15);
    }
  }
}

// _pmtlayout == 5    - 5 row's design for the PANDA Barrel DIRC
// _pmtlayout == 2015 - cern 2015
// _pmtlayout == 2016 - cern 2016
// _pmtlayout == 2017 - cern 2017
// _pmtlayout == 2018 - cern 2018
// _pmtlayout == 2021 - new 3.6 row's design for the PANDA Barrel DIRC
// _pmtlayout == 2023 - new 2x4 layout for the PANDA Barrel DIRC
// _pmtlayout == 2031 - EIC DIRC beam test
// _pmtlayout == 2030 - EIC DIRC prism
// _pmtlayout == 2032 - EIC DIRC focusing prism
TCanvas *PrtTools::draw_digi(double maxz, double minz, TCanvas *cdigi) {

  _last_maxz = maxz;
  _last_minz = minz;

  TString sid = rand_str(3);
  if (cdigi)
    cdigi->cd();
  else
    cdigi = new TCanvas("hp=" + sid, "hp_" + sid, 800, 400);

  TPad *pads[_nmaxpmt];
  TPad *toppad;

  if (_pmtlayout == 2015 || _pmtlayout == 5)
    toppad = new TPad(sid, "T", 0.04, 0.04, 0.88, 0.96);
  else if (_pmtlayout == 2021)
    toppad = new TPad(sid, "T", 0.12, 0.02, 0.78, 0.98);
  else if (_pmtlayout == 2016)
    toppad = new TPad(sid, "T", 0.2, 0.02, 0.75, 0.98);
  else if (_pmtlayout == 2017)
    toppad = new TPad(sid, "T", 0.15, 0.02, 0.80, 0.98);
  else if (_pmtlayout == 2018)
    toppad = new TPad(sid, "T", 0.05, 0.07, 0.9, 0.93);
  else if (_pmtlayout == 2023)
    toppad = new TPad(sid, "T", 0.073, 0.02, 0.877, 0.98);
  else if (_pmtlayout == 2030)
    toppad = new TPad(sid, "T", 0.10, 0.025, 0.82, 0.975);
  else if (_pmtlayout == 2032)
    toppad = new TPad(sid, "T", 0.04, 0.025, 0.91, 0.975);
  else if (_pmtlayout == 2031)
    toppad = new TPad(sid, "T", 0.12, 0.01, 0.80, 0.99);
  else
    toppad = new TPad(sid, "T", 0.04, 0.04, 0.96, 0.96);

  toppad->SetFillStyle(0);
  toppad->Draw();
  toppad->cd();

  int nrow = 3, ncol = 5;
  if (_pmtlayout == 2016) ncol = 3;
  if (_pmtlayout == 2017) ncol = 4;
  if (_pmtlayout == 2018 || _pmtlayout == 2023) {
    nrow = 2;
    ncol = 4;
  }
  if (_pmtlayout == 2021) ncol = 4;
  if (_pmtlayout == 2030) {
    nrow = 4;
    ncol = 6;
  }
  if (_pmtlayout == 2032) {
    nrow = 4;
    ncol = 7;
  }
  if (_pmtlayout == 2031) {
    nrow = 3;
    ncol = 4;
  }

  if (_pmtlayout > 1) {
    float tbw(0.02), tbh(0.01), shift(0), shiftw(0.02), shifth(0), margin(0.01);
    int padi(0);

    for (int i = 0; i < ncol; i++) {
      for (int j = 0; j < nrow; j++) {
        if (j == 1)
          shift = -0.028;
        else
          shift = 0;
        shifth = 0;
        if (_pmtlayout == 5) {
          shift = 0;
          shiftw = 0.001;
          tbw = 0.001;
          tbh = 0.001;
        }
        if (_pmtlayout == 2021) {
          if (i == 0 && j == nrow - 1) continue;
          shift = 0;
          shiftw = 0.001;
          tbw = 0.001;
          tbh = 0.001;
          if (i == 0) shifth = 0.167;
        }
        if (_pmtlayout == 2016) {
          shift = -0.01;
          shiftw = 0.01;
          tbw = 0.03;
          tbh = 0.006;
          if (j == 1) shift += 0.015;
        }
        if (_pmtlayout == 2017) {
          margin = 0.1;
          shift = 0;
          shiftw = 0.01;
          tbw = 0.005;
          tbh = 0.006;
        }
        if (_pmtlayout == 2018) {
          margin = 0.1;
          shift = 0;
          shiftw = 0.01;
          tbw = 0.005;
          tbh = 0.006;
        }
        if (_pmtlayout == 2023) {
          margin = 0.1;
          shift = 0;
          shiftw = 0.01;
          tbw = 0.0015;
          tbh = 0.042;
        }
        if (_pmtlayout == 2030) {
          margin = 0.1;
          shift = 0;
          shiftw = 0.01;
          tbw = 0.001;
          tbh = 0.001;
          padi = j * ncol + i;
        }
        if (_pmtlayout == 2032) {
          margin = 0.1;
          shift = 0;
          shiftw = 0.01;
          tbw = 0.001;
          tbh = 0.001;
          padi = j * ncol + i;
        }
        if (_pmtlayout == 2031) {
          margin = 0.1;
          shift = 0;
          shiftw = 0.01;
          tbw = 0.001;
          tbh = 0.001;
          padi = i * nrow + j;
        }

        pads[padi] = new TPad(
          sid + Form("P%d", i * 10 + j), "T", i / (ncol + 2 * margin) + tbw + shift + shiftw,
          j / (double)nrow + tbh + shifth, (i + 1) / (ncol + 2 * margin) - tbw + shift + shiftw,
          (1 + j) / (double)nrow - tbh + shifth, 21);
        pads[padi]->SetFillColor(kCyan - 8);
        pads[padi]->SetMargin(0.055, 0.055, 0.055, 0.055);
        pads[padi]->Draw();
        padi++;
      }
    }
  } else {
    float tbw(0.02), tbh(0.01), shift(0), shiftw(-0.02);
    int padi(0);
    for (int ii = 0; ii < ncol; ii++) {
      for (int j = 0; j < nrow; j++) {
        if (j == 1)
          shift = 0.04;
        else
          shift = 0;
        pads[padi] =
          new TPad(Form("P%d", ii * 10 + j), "T", ii / (double)ncol + tbw + shift + shiftw,
                   j / (double)nrow + tbh, (ii + 1) / (double)ncol - tbw + shift + shiftw,
                   (1 + j) / (double)nrow - tbh, 21);
        pads[padi]->SetFillColor(kCyan - 8);
        pads[padi]->SetMargin(0.04, 0.04, 0.04, 0.04);
        pads[padi]->Draw();
        padi++;
      }
    }
  }

  int np;
  double max = 0;

  {
    double tmax;
    if (maxz == 0) {
      for (int p = 0; p < nrow * ncol; p++) {
        tmax = _hdigi[p]->GetBinContent(_hdigi[p]->GetMaximumBin());
        if (max < tmax) max = tmax;
      }
    } else {
      max = maxz;
    }

    if (maxz == -2 || minz == -2) { // optimize range
      for (int p = 0; p < nrow * ncol; p++) {
        tmax = _hdigi[p]->GetMaximum();
        if (max < tmax) max = tmax;
      }
      int tbins = 2000;
      TH1F *h = new TH1F("", "", tbins, 0, max);
      for (int p = 0; p < nrow * ncol; p++) {
        for (int i = 0; i < 8; i++) {
          for (int j = 0; j < 8; j++) {
            double val = _hdigi[p]->GetBinContent(i + 1, j + 1);
            if (val != 0) h->Fill(val);
          }
        }
      }
      double integral;
      for (int i = 0; i < tbins; i++) {
        integral = h->Integral(0, i);
        if (integral > 0) {
          if (minz == -2) minz = h->GetBinCenter(i);
          break;
        }
      }

      for (int i = tbins; i > 0; i--) {
        integral = h->Integral(i, tbins);
        if (integral > 10) {
          if (maxz == -2) max = h->GetBinCenter(i);
          break;
        }
      }
    }
  }

  _last_max = max;
  int nnmax(0);
  for (int p = 0; p < nrow * ncol; p++) {
    if (_pmtlayout == 1 || _pmtlayout == 4)
      np = p % nrow * ncol + p / 3;
    else
      np = p;

    if (_pmtlayout == 6 && p > 10) continue;

    pads[p]->cd();
    _hdigi[np]->Draw("col"); //"col+text"
    if (maxz == -1) max = _hdigi[np]->GetBinContent(_hdigi[np]->GetMaximumBin());
    if (nnmax < _hdigi[np]->GetEntries()) nnmax = np;
    _hdigi[np]->SetMaximum(max);
    _hdigi[np]->SetMinimum(minz);
  }

  cdigi->cd();
  TPaletteAxis *palette;
  if (_pmtlayout == 2018 || _pmtlayout == 2023)
    palette = new TPaletteAxis(0.89, 0.1, 0.93, 0.90, (TH1 *)_hdigi[nnmax]);
  else if (_pmtlayout == 2032)
    palette = new TPaletteAxis(0.91, 0.1, 0.94, 0.90, (TH1 *)_hdigi[nnmax]);
  else
    palette = new TPaletteAxis(0.82, 0.1, 0.86, 0.90, (TH1 *)_hdigi[nnmax]);

  palette->SetName("prt_palette");
  palette->Draw();

  cdigi->Modified();
  cdigi->Update();

  return cdigi;
}

TString PrtTools::pix_digi(TString s) {
  int nrow = 3, ncol = 5, np, nmax = 0, npix = 8;
  if (_pmtlayout == 2016) ncol = 3;
  if (_pmtlayout == 2017) ncol = 4;
  if (_pmtlayout == 2018 || _pmtlayout == 2023) {
    nrow = 2;
    ncol = 4;
  }
  if (_pmtlayout == 2021) ncol = 4;
  if (_pmtlayout == 2030) {
    nrow = 4;
    ncol = 6;
    npix = 16;
  }
  if (_pmtlayout == 2032) {
    nrow = 4;
    ncol = 7;
    npix = 16;
  }
  if (_pmtlayout == 2031) {
    nrow = 3;
    ncol = 4;
    npix = 16;
  }

  int nnmax(0);
  for (int p = 0; p < nrow * ncol; p++) {
    if (_pmtlayout == 1 || _pmtlayout == 4)
      np = p % nrow * ncol + p / 3;
    else if (_pmtlayout == 2030)
      np = p % ncol * nrow + p / ncol;
    else
      np = p;

    if (_last_maxz == -1)
      _last_max = _hdigi[p]->GetBinContent(_hdigi[p]->GetMaximumBin());
    if (nnmax < _hdigi[p]->GetEntries()) nnmax = p;
    _hdigi[p]->SetMaximum(_last_max);
    _hdigi[p]->SetMinimum(_last_minz);
    for (int i = 1; i <= npix; i++) {
      for (int j = 1; j <= npix; j++) {
        double weight = (double)(_hdigi[p]->GetBinContent(i, j)) / (double)_last_max * 255;
        if (weight > 255) weight = 255;
        if (weight > 0) s += Form("%d,%d,%d\n", np, (i - 1) * npix + j - 1, (int)weight);
      }
    }
  }
  return s;
}

bool PrtTools::next(int i, int printstep) {

  _chain->GetEntry(i);

  if (i % printstep == 0 && i != 0)
    std::cout << "Event # " << i << " # hits " << _event->getHits().size() << std::endl;

  int pid = _event->getPid();
  if (pid < 5000 && pid > 0) {
    _pid = map_pid[pid];
  }

  return true;
}

bool PrtTools::is_bad_channel(int ch) {
  if (ch < 0 || ch >= _maxdircch) return true;
  return false;
}

bool PrtTools::read_db(TString in) {

  std::cout << "Parsing " << in << "=============" << std::endl;

  std::ifstream ins(in);
  double dx = 0, vx = 0, dy = 0, vy = 0;
  TString info, name;
  std::string line, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, s13;
  int study = 0, fileid = 0, radiatorid = 0, lensid = 0;
  double theta = 0, phi = 0, z = 0, x = 0, sx = 0, sy = 0, mom = 0, beamsize = 0, simo = 0;

  while (std::getline(ins, line)) {
    std::istringstream iss(line);

    if (line.rfind("S", 0) == 0) { // parse header
      line.erase(0, 1);
      std::istringstream iss(line);
      if (iss >> s1) {
        s2 = iss.str();
        study = atoi(s1.c_str());
        info = s2.c_str();
      } else {
        std::cout << "Error parsing data_db: " << line << std::endl;
      }
    } else {
      if (iss >> s1 >> s2 >> s3 >> s4 >> s5 >> s6 >> s7 >> s8 >> s9 >> s10 >> s11 >> s12 >> s13) {
        name = s1.c_str();
        fileid = atoi(s2.c_str());
        radiatorid = atoi(s3.c_str());
        lensid = atoi(s4.c_str());
        theta = atof(s5.c_str());
        phi = atof(s6.c_str());
        z = atof(s7.c_str());
        x = atof(s8.c_str());
        sx = atof(s9.c_str());
        sy = atof(s10.c_str());
        mom = atof(s11.c_str());
        beamsize = atof(s12.c_str());
        simo = atof(s13.c_str());
      } else {
        // std::cout<<"Error parsing data_db: "<<line<<std::endl;
      }
    }

    PrtRun *r = new PrtRun();
    r->setInfo(info);
    r->setName(name);
    r->setInfo(info);
    r->setId(fileid);
    r->setRadiator(radiatorid);
    r->setLens(lensid);
    r->setTheta(theta);
    r->setPhi(phi);
    r->setBeamX(x);
    r->setBeamZ(z);
    r->setPrismStepX(sx);
    r->setPrismStepY(sy);
    r->setMomentum(mom);
    r->setBeamSize(beamsize);

    _runs.push_back(r);
  }

  std::cout << "Parsed " << _runs.size() << " runs ================" << std::endl;

  return 1;
}

void PrtTools::set_palette(int pal) {

  // pal =  1: rainbow\n"
  // pal =  2: reverse-rainbow\n"
  // pal =  3: amber\n"
  // pal =  4: reverse-amber\n"
  // pal =  5: blue/white\n"
  // pal =  6: white/blue\n"
  // pal =  7: red temperature\n"
  // pal =  8: reverse-red temperature\n"
  // pal =  9: green/white\n"
  // pal = 10: white/green\n"
  // pal = 11: orange/blue\n"
  // pal = 12: blue/orange\n"
  // pal = 13: white/black\n"
  // pal = 14: black/white\n"

  const int NRGBs = 5;
  const int NCont = 255;
  gStyle->SetNumberContours(NCont);

  if (pal < 1 && pal > 15)
    return;
  else
    pal--;

  double stops[NRGBs] = {0.00, 0.34, 0.61, 0.84, 1.00};
  double red[15][NRGBs] = {
    {0.00, 0.00, 0.87, 1.00, 0.51}, {0.51, 1.00, 0.87, 0.00, 0.00}, {0.17, 0.39, 0.62, 0.79, 1.00},
    {1.00, 0.79, 0.62, 0.39, 0.17}, {0.00, 0.00, 0.00, 0.38, 1.00}, {1.00, 0.38, 0.00, 0.00, 0.00},
    {0.00, 0.50, 0.89, 0.95, 1.00}, {1.00, 0.95, 0.89, 0.50, 0.00}, {0.00, 0.00, 0.38, 0.75, 1.00},
    {0.00, 0.34, 0.61, 0.84, 1.00}, {0.75, 1.00, 0.24, 0.00, 0.00}, {0.00, 0.00, 0.24, 1.00, 0.75},
    {0.00, 0.34, 0.61, 0.84, 1.00}, {1.00, 0.84, 0.61, 0.34, 0.00}, {0.00, 0.00, 0.80, 1.00, 0.80}};
  double green[15][NRGBs] = {
    {0.00, 0.81, 1.00, 0.20, 0.00}, {0.00, 0.20, 1.00, 0.81, 0.00}, {0.01, 0.02, 0.39, 0.68, 1.00},
    {1.00, 0.68, 0.39, 0.02, 0.01}, {0.00, 0.00, 0.38, 0.76, 1.00}, {1.00, 0.76, 0.38, 0.00, 0.00},
    {0.00, 0.00, 0.27, 0.71, 1.00}, {1.00, 0.71, 0.27, 0.00, 0.00}, {0.00, 0.35, 0.62, 0.85, 1.00},
    {1.00, 0.75, 0.38, 0.00, 0.00}, {0.24, 1.00, 0.75, 0.18, 0.00}, {0.00, 0.18, 0.75, 1.00, 0.24},
    {0.00, 0.34, 0.61, 0.84, 1.00}, {1.00, 0.84, 0.61, 0.34, 0.00}, {0.00, 0.85, 1.00, 0.30, 0.00}};
  double blue[15][NRGBs] = {
    {0.51, 1.00, 0.12, 0.00, 0.00}, {0.00, 0.00, 0.12, 1.00, 0.51}, {0.00, 0.09, 0.18, 0.09, 0.00},
    {0.00, 0.09, 0.18, 0.09, 0.00}, {0.00, 0.47, 0.83, 1.00, 1.00}, {1.00, 1.00, 0.83, 0.47, 0.00},
    {0.00, 0.00, 0.00, 0.40, 1.00}, {1.00, 0.40, 0.00, 0.00, 0.00}, {0.00, 0.00, 0.00, 0.47, 1.00},
    {1.00, 0.47, 0.00, 0.00, 0.00}, {0.00, 0.62, 1.00, 0.68, 0.12}, {0.12, 0.68, 1.00, 0.62, 0.00},
    {0.00, 0.34, 0.61, 0.84, 1.00}, {1.00, 0.84, 0.61, 0.34, 0.00}, {0.60, 1.00, 0.10, 0.00, 0.00}};

  TColor::CreateGradientColorTable(NRGBs, stops, red[pal], green[pal], blue[pal], NCont);
}

void PrtTools::create_maps(int pmtlayout) {

  if (pmtlayout == 2019) {
    for (size_t i = 0; i < _tdcsid_jul2019.size(); i++) {
      int dec = TString::BaseConvert(_tdcsid_jul2019[i], 16, 10).Atoi();
      map_tdc[dec] = i;
    }
  } else {
    for (size_t i = 0; i < _tdcsid_jul2018.size(); i++) {
      int dec = TString::BaseConvert(_tdcsid_jul2018[i], 16, 10).Atoi();
      map_tdc[dec] = i;
    }
  }

  for (int ch = 0; ch < _maxch; ch++) {
    int pmt = ch / _npix;
    int pix = ch % _npix;
    int col = pix / 2 - 8 * (pix / 8);
    int row = pix % 2 + 2 * (pix / 8);
    pix = col + sqrt(_npix) * row;

    map_pmt[ch] = pmt;
    map_pix[ch] = pix;
    map_row[ch] = row;
    map_col[ch] = col;
  }

  for (int i = 0; i < 5; i++) {
    map_pid[_pdg[i]] = i;
  }
}

TString PrtTools::rand_str(int len) {
  TString str = "";
  static const char alphanum[] = "0123456789"
                                 "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                                 "abcdefghijklmnopqrstuvwxyz";

  for (int i = 0; i < len; ++i) {
    str += alphanum[rand() % (sizeof(alphanum) - 1)];
  }
  return str;
}

TVector3 PrtTools::fit(TH1 *h, double range, double threshold, double limit, int peakSearch,
                       int bkg, TString opt) {

  int binmax = h->GetMaximumBin();
  auto ax = h->GetXaxis();
  double xmax = ax->GetBinCenter(binmax);
  TString sfun = "[0]*exp(-0.5*((x-[1])/[2])^2)";
  if (bkg == 1) sfun = "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]+x*[4]";

  _fgaus = new TF1("_fgaus", sfun, xmax - range, xmax + range);
  _fgaus->SetNpx(500);
  _fgaus->SetParNames("const", "mean", "sigma");
  _fgaus->SetLineColor(2);

  double integral = h->Integral(ax->FindBin(xmax - range), ax->FindBin(xmax + range));
  double xxmin, xxmax, sigma1(0), mean1(0), mean2(0);
  xxmax = xmax;
  xxmin = xxmax;
  int nfound(1);
  if (integral > threshold) {

    if (peakSearch == 1) {
      _fgaus->SetParameter(1, xmax);
      _fgaus->SetParameter(2, 0.005);
      _fgaus->SetParLimits(2, 0.003, limit);
      h->Fit("_fgaus", opt, "", xxmin - range, xxmax + range);
    }

    if (peakSearch > 1) {
      nfound = _spectrum->Search(h, 4, "goff", 0.1);
      std::cout << "nfound  " << nfound << std::endl;
      if (nfound == 1) {
        _fgaus = new TF1("_fgaus", "gaus(0)", xmax - range, xmax + range);
        _fgaus->SetNpx(500);
        _fgaus->SetParameter(1, _spectrum->GetPositionX()[0]);
      } else if (nfound >= 2) {
        double p1 = _spectrum->GetPositionX()[0];
        double p2 = _spectrum->GetPositionX()[1];
        if (p1 > p2) {
          xxmax = p1;
          xxmin = p2;
        } else {
          xxmax = p2;
          xxmin = p1;
        }
        if (peakSearch == 20) {
          xxmax = xxmin;
          _fgaus = new TF1("_fgaus", "gaus(0)", xxmin - range, xxmin + range);
          _fgaus->SetNpx(500);
          _fgaus->SetParameter(1, _spectrum->GetPositionX()[0]);
        } else {
          _fgaus = new TF1("_fgaus", "gaus(0)+gaus(3)", xmax - range, xmax + range);
          _fgaus->SetNpx(500);
          _fgaus->SetParameter(0, 1000);
          _fgaus->SetParameter(3, 1000);

          _fgaus->FixParameter(1, xxmin);
          _fgaus->FixParameter(4, xxmax);
          _fgaus->SetParameter(2, 0.1);
          _fgaus->SetParameter(5, 0.1);
          h->Fit("_fgaus", opt, "", xxmin - range, xxmax + range);
          _fgaus->ReleaseParameter(1);
          _fgaus->ReleaseParameter(4);
        }
      }

      _fgaus->SetParameter(2, 0.2);
      _fgaus->SetParameter(5, 0.2);
    }

    h->Fit("_fgaus", opt, "", xxmin - range, xxmax + range);
    mean1 = _fgaus->GetParameter(1);
    sigma1 = _fgaus->GetParameter(2);
    if (sigma1 > 10) sigma1 = 10;

    if (peakSearch == 2) {
      mean2 = (nfound == 1) ? _fgaus->GetParameter(1) : _fgaus->GetParameter(4);
    }
  }
  delete _fgaus;
  return TVector3(mean1, sigma1, mean2);
}

TGraph *PrtTools::fit_slices(TH2F *h, double minrange, double maxrange, double fitrange, int rebin,
                             int ret) {
  TH2F *ht = (TH2F *)h->Clone("ht");
  ht->RebinY(rebin);
  int point(0);
  TGraph *gres = new TGraph();
  for (int i = 1; i < ht->GetNbinsY(); i++) {
    double x = ht->GetYaxis()->GetBinCenter(i);
    TH1D *hp;
    if (minrange != maxrange) {
      TCutG *cut = new TCutG("onepeakcut", 5);
      cut->SetVarX("y");
      cut->SetVarY("x");
      cut->SetPoint(0, minrange, -1E6);
      cut->SetPoint(1, minrange, 1E6);
      cut->SetPoint(2, maxrange, 1E6);
      cut->SetPoint(3, maxrange, -1E6);
      cut->SetPoint(4, minrange, -1E6);

      hp = ht->ProjectionX(Form("bin%d", i), i, i, "[onepeakcut]");
    } else {
      hp = ht->ProjectionX(Form("bin%d", i), i, i);
    }

    TVector3 res = fit((TH1F *)hp, fitrange, 100, 2, 1, 1);
    double y = 0;
    if (ret == 0) y = res.X();
    if (ret == 1) y = res.Y();
    if (ret == 2) y = res.X() + 0.5 * res.Y();
    if (ret == 3) y = res.X() - 0.5 * res.Y();
    if (y == 0 || y < minrange || y > maxrange) continue;

    gres->SetPoint(point, y, x);
    gres->SetLineWidth(2);
    gres->SetLineColor(kRed);
    point++;
  }
  return gres;
}

void style_graph(TGraph *g, int id) {
  int coll[] = {kBlack, kRed + 1, kGreen, kBlue, 4, kCyan - 6, kOrange, 7, 8, 9, 10};
  int colm[] = {kBlack, kRed + 1, kGreen + 2, kBlue + 1, 4, kCyan - 6, kOrange + 1, 7, 8, 9, 10};

  int cl = (id < 10) ? coll[id] : id;
  int cm = (id < 10) ? colm[id] : id;
  g->SetLineColor(cl);
  g->SetMarkerColor(cm);
  g->SetMarkerStyle(20);
  g->SetMarkerSize(0.8);
  g->SetName(Form("gr_%d", id));
}

double PrtTools::integral(TH1F *h, double xmin, double xmax) {
  TAxis *axis = h->GetXaxis();
  int bmin = axis->FindBin(xmin);
  int bmax = axis->FindBin(xmax);
  double integral = h->Integral(bmin, bmax);
  integral -= h->GetBinContent(bmin) * (xmin - axis->GetBinLowEdge(bmin)) / axis->GetBinWidth(bmin);
  integral -= h->GetBinContent(bmax) * (axis->GetBinUpEdge(bmax) - xmax) / axis->GetBinWidth(bmax);
  return integral;
}

void PrtTools::normalize(TH1F *hists[], int size) {

  // for(int i=0; i<size; i++){
  //   hists[i]->Scale(1/hists[i]->Integral(), "width");
  // }

  double max = 0, min = 0;
  for (int i = 0; i < size; i++) {
    double tmax = hists[i]->GetBinContent(hists[i]->GetMaximumBin());
    double tmin = hists[i]->GetMinimum();
    if (tmax > max) max = tmax;
    if (tmin < min) min = tmin;
  }
  max += 0.05 * max;
  for (int i = 0; i < size; i++) {
    hists[i]->GetYaxis()->SetRangeUser(min, max);
  }
}

void PrtTools::normalize_to(TH1F *hists[], int size, double max) {

  for (int i = 0; i < size; i++) {
    double tmax = hists[i]->GetBinContent(hists[i]->GetMaximumBin());
    if (tmax > 0) hists[i]->Scale(max / tmax);
  }
}

void PrtTools::normalize(TH1F *h1, TH1F *h2) {
  double max = (h1->GetMaximum() > h2->GetMaximum()) ? h1->GetMaximum() : h2->GetMaximum();
  max += max * 0.1;
  h1->GetYaxis()->SetRangeUser(0, max);
  h2->GetYaxis()->SetRangeUser(0, max);
}

// just x for now
TGraph *PrtTools::smooth(TGraph *g, int smoothness) {
  double x, y;
  int n = g->GetN();
  TH1F *h = new TH1F("h", "h", g->GetN(), 0, n);
  TGraph *gr = new TGraph();
  gr->SetName(g->GetName());
  for (auto i = 0; i < n; i++) {
    g->GetPoint(i, x, y);
    h->Fill(i, x);
  }

  h->Smooth(smoothness);

  for (auto i = 0; i < n; i++) {
    g->GetPoint(i, x, y);
    gr->SetPoint(i, h->GetBinContent(i), y);
  }
  return gr;
}

int PrtTools::shift_hist(TH1 *hist, double double_shift) {

  int bins = hist->GetXaxis()->GetNbins();
  double xmin = hist->GetXaxis()->GetBinLowEdge(1);
  double xmax = hist->GetXaxis()->GetBinUpEdge(bins);
  double_shift = double_shift * (bins / (xmax - xmin));
  int shift = 0;
  if (double_shift < 0) shift = TMath::FloorNint(double_shift);
  if (double_shift > 0) shift = TMath::CeilNint(double_shift);
  if (shift == 0) return 0;
  if (shift > 0) {
    for (int i = 1; i <= bins; i++) {
      if (i + shift <= bins) hist->SetBinContent(i, hist->GetBinContent(i + shift));
      if (i + shift > bins) hist->SetBinContent(i, 0);
    }
    return 0;
  }
  if (shift < 0) {
    for (int i = bins; i > 0; i--) {
      if (i + shift > 0) hist->SetBinContent(i, hist->GetBinContent(i + shift));
      if (i + shift <= 0) hist->SetBinContent(i, 0);
    }
    return 0;
  }
  return 1;
}

void PrtTools::save_canvas(int what, int style, bool rm){
  TIter next(_canvaslist);
  TCanvas *c=0;
  TString path = create_dir();
  while((c = (TCanvas*) next())){
    set_style(c);
    save(c, path, what,style);
    if(rm){
      _canvaslist->Remove(c);
      c->Close();
    }
  }
}

// path - folder for saving
void PrtTools::save_canvas(TString path, int what, int style, bool rm) {
  _savepath = path;
  save_canvas(what, style, rm);
}

void PrtTools::add_canvas(TString name,int w, int h){
  if(!get_canvas(name)){
    if(!_canvaslist) _canvaslist = new TList();
    TCanvas *c = new TCanvas(name,name,0,0,w,h);
    _canvaslist->Add(c);
  }
}

void PrtTools::add_canvas(TCanvas *c){
  if(!_canvaslist) _canvaslist = new TList();
  c->cd();
  _canvaslist->Add(c);
}

TCanvas *PrtTools::get_canvas(TString name){
  TIter next(_canvaslist);
  TCanvas *c = 0;
  while((c = (TCanvas*) next())){
    if(c->GetName()==name || name=="*") break;
  }
  return c;
}

void PrtTools::del_canvas(TString name){
  TIter next(_canvaslist);
  TCanvas *c=0;
  while((c = (TCanvas*) next())){
    if(c->GetName()==name || name=="*") _canvaslist->Remove(c);
    c->Delete();
  }
}

void PrtTools::wait_primitive(TString name, TString prim) {
  TIter next(_canvaslist);
  TCanvas *c = 0;
  while ((c = (TCanvas *)next())) {
    if (TString(c->GetName()) == name) {
      c->Modified();
      c->Update();
      c->WaitPrimitive(prim);
    }
  }
}

void PrtTools::set_style() {
  TIter next(_canvaslist);
  TCanvas *c = 0;
  TString path = create_dir();
  while ((c = (TCanvas *)next())) {
    set_style(c);
  }
}

void PrtTools::set_style(TCanvas *c) {

  if(fabs(c->GetBottomMargin()-0.1)<0.001) c->SetBottomMargin(0.12);
  TIter next(c->GetListOfPrimitives());
  TObject *obj;

  while((obj = next())){
    if(obj->InheritsFrom("TH1")){
      TH1F *hh = (TH1F*)obj;
      hh->GetXaxis()->SetTitleSize(0.06);
      hh->GetYaxis()->SetTitleSize(0.06);

      hh->GetXaxis()->SetLabelSize(0.05);
      hh->GetYaxis()->SetLabelSize(0.05);

      hh->GetXaxis()->SetTitleOffset(0.85);
      hh->GetYaxis()->SetTitleOffset(0.76);

      if(c->GetWindowHeight()>700){
	hh->GetXaxis()->SetTitleSize(0.045);
	hh->GetYaxis()->SetTitleSize(0.045);

	hh->GetXaxis()->SetLabelSize(0.035);
	hh->GetYaxis()->SetLabelSize(0.035);

	hh->GetXaxis()->SetTitleOffset(0.85);
	hh->GetYaxis()->SetTitleOffset(0.98);
	c->SetRightMargin(0.12);
      }

      if(c->GetWindowWidth()<=600){
      	hh->GetXaxis()->SetTitleSize(0.05);
      	hh->GetYaxis()->SetTitleSize(0.05);

      	hh->GetXaxis()->SetLabelSize(0.04);
      	hh->GetYaxis()->SetLabelSize(0.04);

      	hh->GetXaxis()->SetTitleOffset(0.85);
      	hh->GetYaxis()->SetTitleOffset(0.92);
      	c->SetRightMargin(0.12);
      }

      if(fabs(c->GetBottomMargin()-0.12)<0.001){
	TPaletteAxis *palette = (TPaletteAxis*)hh->GetListOfFunctions()->FindObject("palette");
	if(palette) {
	  palette->SetY1NDC(0.12);
	  c->Modified();
	}
      }
      c->Modified();
      c->Update();
    }

    if(obj->InheritsFrom("TGraph")){
      TGraph *gg = (TGraph*)obj;
      gg->GetXaxis()->SetLabelSize(0.05);
      gg->GetXaxis()->SetTitleSize(0.06);
      gg->GetXaxis()->SetTitleOffset(0.84);

      gg->GetYaxis()->SetLabelSize(0.05);
      gg->GetYaxis()->SetTitleSize(0.06);
      gg->GetYaxis()->SetTitleOffset(0.8);
    }

    if(obj->InheritsFrom("TMultiGraph")){
      TMultiGraph *gg = (TMultiGraph*)obj;
      gg->GetXaxis()->SetLabelSize(0.05);
      gg->GetXaxis()->SetTitleSize(0.06);
      gg->GetXaxis()->SetTitleOffset(0.84);

      gg->GetYaxis()->SetLabelSize(0.05);
      gg->GetYaxis()->SetTitleSize(0.06);
      gg->GetYaxis()->SetTitleOffset(0.8);
    }

    if(obj->InheritsFrom("TF1")){
      TF1 *f = (TF1*)obj;
      f->SetNpx(500);
    }
  }
}

void PrtTools::save(TPad *c,TString path, int what, int style){
  TString name = c->GetName();
  bool batch = gROOT->IsBatch();
  gROOT->SetBatch(1);

  if(c && path != "") {
    int w = 800, h = 400;
    if(style != -1){
      if(style == 1) {w = 800; h = 500;}
      if(style == 2) {w = 800; h = 600;}
      if(style == 3) {w = 800; h = 400;}
      if(style == 5) {w = 800; h = 900;}
      if(style == 0){
    	w = ((TCanvas*)c)->GetWindowWidth();
    	h = ((TCanvas*)c)->GetWindowHeight();
      }

      TCanvas *cc;
      if(TString(c->GetName()).Contains("hp") || TString(c->GetName()).Contains("cdigi")) {
	cc = draw_digi(_last_maxz,_last_minz);
	cc->SetCanvasSize(800,400);
	if(name.Contains("=")) name =  name.Tokenize('=')->First()->GetName();
      }else{
      	cc = new TCanvas(TString(c->GetName())+"exp","cExport",0,0,w,h);
      	cc = (TCanvas*) c->DrawClone();
	cc->SetCanvasSize(w,h);
      }

      if(style == 0) set_style(cc);

      print_canvas(cc,name,path,what);
    }else{
      c->SetCanvasSize(w,h);
      print_canvas(c,name,path,what);
    }
  }

  gROOT->SetBatch(batch);
}

void PrtTools::print_canvas(TPad *c, TString name, TString path, int what) {
  c->Modified();
  c->Update();
  c->Print(path + "/" + name + ".png");
  if (what > 0) c->Print(path + "/" + name + ".C");
  if (what > 1) c->Print(path + "/" + name + ".pdf");
  if (what > 2) c->Print(path + "/" + name + ".eps");
}

TString PrtTools::create_dir(TString inpath) {
  if (inpath != "") _savepath = inpath;
  TString finalpath = _savepath;

  if (finalpath == "") return "";

  if (_savepath.EndsWith("auto")) {
    TString dir = _savepath.ReplaceAll("auto", "data");
    gSystem->mkdir(dir);
    TDatime *time = new TDatime();
    TString path(""), stime = Form("%d.%d.%d", time->GetDay(), time->GetMonth(), time->GetYear());
    gSystem->mkdir(dir + "/" + stime);
    for (int i = 0; i < 1000; i++) {
      path = stime + "/" + Form("arid-%d", i);
      if (gSystem->mkdir(dir + "/" + path) == 0) break;
    }
    gSystem->Unlink(dir + "/last");
    gSystem->Symlink(path, dir + "/last");
    finalpath = dir + "/" + path;
    _savepath = finalpath;
  } else {
    gSystem->mkdir(_savepath, kTRUE);
  }

  write_info(finalpath + "/readme");
  return finalpath;
}

void PrtTools::write_info(TString filename) {
  std::ofstream myfile;
  myfile.open(filename);
  myfile << _info + "\n";
  myfile.close();
}

void PrtTools::write_string(TString filename, TString str) {
  std::ofstream myfile;
  myfile.open(filename);
  myfile << str + "\n";
  myfile.close();
  std::cout << "output: " << filename << std::endl;
}

int PrtTools::get_pid(int pdg) {
  int pid = 0;
  if (pdg == 11) pid = 0;   // e
  if (pdg == 13) pid = 1;   // mu
  if (pdg == 211) pid = 2;  // pi
  if (pdg == 321) pid = 3;  // K
  if (pdg == 2212) pid = 4; // p
  return pid;
}
