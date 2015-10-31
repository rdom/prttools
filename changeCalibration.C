void changeCalibration(TString ifile="calib.root", TString ofile="newcalib.root"){
  TFile inf(ifile);
  TIter nextkey(inf.GetListOfKeys());
  TKey *key;

  TFile ouf(ofile,"RECREATE");
  
  while ((key = (TKey*)nextkey())) {
    TGraph *gr = (TGraph*)key->ReadObj();
    TString name = gr->GetName();
    if(name.BeginsWith("15") && name.Length()>6) continue;
    gr->Write();
  }
  inf.Close();
  ouf.Write();
  ouf.Close();
}
