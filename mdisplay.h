// mdisplay - tool to plot different quantities from the *M.root file
// original author: Roman Dzhygadlo - GSI Darmstadt 

#ifndef mdisplay_h
#define mdisplay_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TH2.h>
#include <vector>

class MyMainFrame : public TGMainFrame {
private:
  TRootEmbeddedCanvas  *fEcan;
  TRootEmbeddedCanvas  *fTime;
  TGStatusBar          *fStatusBar;
  TGNumberEntry        *fNumber;
  TGLabel              *fLabel;
  TGLabel              *fLabel1;
  TGComboBox           *fComboMode;
  TGCheckButton        *fCheckBox;
  TGHProgressBar       *fHProg3;
  TGVerticalFrame      *fHm;
  TGTextEntry          *fEdit1;
  TGTextEntry          *fEdit2;
  TGTextEntry          *fEdit4;
  TGCheckButton        *fCheckBtn1;
  TGCheckButton        *fCheckBtn2;
  TGCheckButton        *fCheckBtn3;
  TGCheckButton        *fCheckBtn4;
  TGTextButton         *fBtnMore;

public:
  MyMainFrame(const TGWindow *p, UInt_t w, UInt_t h);
  virtual ~MyMainFrame();
  void DoExit();
  void DoDraw();
  void DoMore();
  void DoExport();
  void DoExportOffsets();
  void DoSlider(Int_t pos =0);
  void DoCheckBtnClecked1();
  void DoCheckBtnClecked2();
  void DoCheckBtnClecked3();
  void DoCheckBtnClecked4();
  TString updatePlot(Int_t id=0, TCanvas *cT=0);
  void SetStatusText(const char *txt, Int_t pi);
  void EventInfo(Int_t event, Int_t px, Int_t py, TObject *selected);

  TGTextEntry          *fEdit3;
  TGHSlider            *fHslider1;
  
  ClassDef(MyMainFrame, 0);
};

class MSelector : public TSelector {
public :
  TTree          *fChain; 
  MSelector(TTree *  =0) : fChain(0) { }
  virtual ~MSelector() { }
  virtual Int_t   Version() const { return 2; }
  virtual void    SlaveBegin(TTree *tree);
  virtual void    Init(TTree *tree);
  virtual Bool_t  Notify(){ return kTRUE; };
  virtual Bool_t  Process(Long64_t entry);
  virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
  virtual void    SetOption(const char *option) { fOption = option; }
  virtual void    SetObject(TObject *obj) { fObject = obj; }
  virtual void    SetInputList(TList *input) { fInput = input; }
  virtual TList  *GetOutputList() const { return fOutput; }
  virtual void    Terminate();

  ClassDef(MSelector,0);
};

#endif
