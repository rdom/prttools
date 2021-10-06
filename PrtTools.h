// prttools - useful functions for hld*, prt*
// original author: Roman Dzhygadlo - GSI Darmstadt

#ifndef prttools_h
#define prttools_h 1

#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TSpline.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TVector3.h"
#include "TMath.h"
#include "TChain.h"
#include "TGaxis.h"
#include "TColor.h"
#include "TString.h"
#include "TArrayD.h"
#include "TSpectrum.h"
#include "TSpectrum2.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TRandom2.h"
#include "TError.h"
#include "TPaveStats.h"
#include "TObjString.h"
#include "TApplication.h"
#include <TLegend.h>
#include <TAxis.h>
#include <TPaletteAxis.h>
#include <TRandom.h>
#include <TCutG.h>
#include <TKey.h>
#include "TPRegexp.h"
#include "TFitResult.h"

#include <iostream>
#include <fstream>
#include <sstream>

#include "../prtdirc/src/PrtHit.h"
#include "../prtdirc/src/PrtEvent.h"
#include "../prtdirc/src/PrtRun.h"

using std::array;

class PrtTools {

 public:
  PrtTools();
  PrtTools(int study);
  PrtTools(TString in);
  PrtTools(PrtRun *run);
  ~PrtTools() {}
  void init(int study = 400);

  bool init_run(TString in = "", int bdigi = 0, TString savepath = "", int setupid = 2019);
  void init_digi();
  TCanvas *draw_digi(double maxz = 0, double minz = 0, TCanvas *cdigi = NULL);
  TString pix_digi(TString header = "m,p,v\n");
  bool next(int i, int printstep = 1000);
  bool next();
  bool is_bad_channel(int ch);
  int get_pid(int pdg);
  bool read_db(TString in = "data_db.dat");
  bool write_db(TString in = "data_db.dat");
  PrtRun *set_run();
  PrtRun *get_run(TString in);
  std::vector<PrtRun *> get_runs(int study);
  void modify_run(PrtRun *run);
  PrtRun *find_run(int sid, int fid = 0);
  PrtRun *find_run(TString path);
  void fill_digi(int pmt, int pix);
  
  void set_palette(int p = 1);
  void set_pmtlayout(int v) { _pmtlayout = v; }
  void create_maps(int setupid = 2019);
  TString get_inpath();
  TString get_outpath();
  TString get_lutpath();
  TString get_pdfpath();  
  int get_channel(int tdc, int tdcch);
  int get_tdcid(int ch);
  TString rand_str(int len = 10);
  TVector3 fit(TH1 *h, double range = 3, double threshold = 20, double limit = 2,
               int peakSearch = 1, int bkg = 0, TString opt = "MQ");
  TGraph *fit_slices(TH2F *h, double minrange = 0, double maxrange = 0, double fitrange = 1,
                     int rebin = 1, int ret = 0);
  void style_graph(TGraph *g, int id);
  double integral(TH1F *h, double xmin, double xmax);
  void normalize(TH1F *h1, TH1F *h2);
  void normalize(TH1F *hists[], int size);
  void normalize_to(TH1F *hists[], int size, double max = 1);
  TGraph *smooth(TGraph *g, int smoothness = 1);
  int shift_hist(TH1 *hist, double double_shift);

  TString dir(TString path);
  TString create_dir(TString inpath = "");
  void write_info(TString filename);
  void write_string(TString filename, TString str);
  
  void set_style();
  void set_style(TCanvas *c);
  void save(TPad *c= NULL,TString path="", int what=0, int style=0);


  void add_canvas(TString name="c",int w=800, int h=400);
  void add_canvas(TCanvas *c);
  TCanvas *get_canvas(TString name="c");
  void del_canvas(TString name="c");

  // style = 0 - for web blog
  // style = 1 - for talk
  // what = 0 - save in png
  // what = 1 - save in png, C
  // what = 2 - save in png, C, pdf
  // what = 3 - save in png, C, pdf, eps
  void save_canvas(int what=1, int style=0, bool rm=false);
  void save_canvas(TString path, int what=1, int style=0, bool rm=false);
  void print_canvas(TPad *c, TString name = "", TString path = "", int what = 0);

  void wait_primitive(TString name, TString prim="");
  
  // accessors
  PrtEvent *event() const { return _event; }
  PrtRun *run() { return _run; }
  TChain *chain() { return _chain; }
  int pdg(int v) { return _pdg[v]; }
  int pid() { return _pid; }
  int i( ) { return _iter; }
  int entries(){ return _entries;}
  double mass(int v) { return _mass[v]; }
  int color(int v) { return _color[v]; }
  TString name(int v) { return _name[v]; }
  TString lname(int v) { return _lname[v]; }
  int maxdircch(){ return _maxdircch;}
  TH2F *get_digi(int v){ return _hdigi[v];}

  // mutators
  void set_path(TString v) { _savepath = v; }

  array<int, 10000> map_pmt{}, map_pix{}, map_row{}, map_col{}, map_tdc{};
  array<array<int, 1000>, 100> map_pmtpix{};

 private:
  TChain *_chain;
  PrtEvent *_event;
  TString _dbpath, _savepath;
  TString _info;
  array<int, 5> _pdg = {11, 13, 211, 321, 2212};
  array<double, 5> _mass = {0.000511, 0.1056584, 0.139570, 0.49368, 0.9382723};
  array<TString, 5> _name = {"e", "muon", "pion", "kaon", "proton"};
  array<TString, 5> _lname = {"e", "#mu", "#pi", "K", "p"};
  array<int, 5> _color = {kOrange + 6, kCyan + 1, kBlue + 1, kRed + 1, kRed + 1};

  array<TString, 32> _tdcsid_jul2018 = {
    "2000", "2001", "2002", "2003", "2004", "2005", "2006", "2007", "2008", "2009", "200a",
    "200b", "200c", "200d", "200e", "200f", "2010", "2011", "2012", "2013", "2014", "2015",
    "2016", "2017", "2018", "2019", "201a", "201b", "201c", "201d", "201e", "201f"};

  array<TString, 21> _tdcsid_jul2019 = {"2014", "2015", "2016", "2017", "2000", "2001", "2002",
                                        "2003", "2004", "2005", "2006", "2007", "2008", "2009",
                                        "200a", "200b", "200c", "200d", "200e", "200f", "2010"};

  std::vector<PrtRun *> _runs;
  PrtRun* _run;
  array<TH2F *, 21> _hdigi{};
  int _entries;
  int _iter = -1;
  int _printstep = 1000;
  int _pid;
  int _maxch, _maxdircch;
  int _npmt, _npix;
  int _pmtlayout;
  int _nmaxpmt = 28;
  int _last_max, _last_maxz, _last_minz;
  TF1 *_fgaus;
  TSpectrum *_spectrum;
  TList *_canvaslist;
};


#endif



// const int prt_ntdcm = 80; //41
// int prt_ntdc = prt_ntdcm;
// int prt_maxch = prt_ntdc*48;



// TString prt_getTdcName(int ch){
//   return prt_tdcsid[prt_getTdcId(ch)];
// }

// int prt_getTdcChannel(int ch){
//   int tch=0,tdcc=0;
//   if(prt_geometry==2018){
//     for(int i=0; i<=prt_ntdc; i++){
//       tdcc=ch-tch+1;
//       if(i==2 || i==6 || i==10 || i==14 || i==1 || i==2 || i==5 || i==6 || i==9 || i==10 || i==13
//       || i==14) tch += 16; else tch += 48; if(tch>ch){
// 	break;
//       };
//     }
//   }else if(prt_geometry==2023){
//     tdcc=ch%32+1;
//   }else{
//     tdcc=ch%48+1;
//   }
//   return tdcc;
// }

// int prt_getColorId(int ind, int style =0){
//   int cid = 1;
//   if(style==0) {
//     cid=ind+1;
//     if(cid==5) cid =8;
//     if(cid==3) cid =kOrange+2;
//   }
//   if(style==1) cid=ind+300;
//   return cid;
// }


// void prt_addInfo(TString str){
//   prt_info += str+"\n";
// }



// TString prt_createSubDir(TString dir="dir"){
//   gSystem->mkdir(dir);
//   return dir;
// }


// double prt_get_momentum_from_tof(double dist,double dtof){
//   double s = dtof*0.299792458/dist;
//   double x = s*s;
//   double a = prt_mass[2]*prt_mass[2]; //pi
//   double b = prt_mass[4]*prt_mass[4]; //p
//   double p = sqrt((a - 2*sqrt(a*a+a*b*x-2*a*b+b*b)/s + b)/(x - 4));

//   return p;
// }

// // return TOF difference [ns] for mominum p [GeV] and flight path l [m]
// double prt_get_tof_diff(int pid1=211, int pid2=321, double p=1, double l=2){
//   double c = 299792458;
//   double m1 = prt_mass[prt_get_pid(pid1)];
//   double m2 = prt_mass[prt_get_pid(pid2)];
//   double td = l*(sqrt(p*p+m1*m1)-sqrt(p*p+m2*m2))/(p*c)*1E9;

//   // relativistic
//   // td = l*(m1*m1-m2*m2)/(2*p*p*c)*1E9;

//   return td;
// }

// bool prt_ispath(TString path){
//   Long_t *id(0),*size(0),*flags(0),*modtime(0);
//   return !gSystem->GetPathInfo(path,id,size,flags,modtime);
// }

// int prt_get3digit(TString str){
//   TPRegexp e("[0-9][0-9][0-9]");
//   return ((TString)str(e)).Atoi();
// }
