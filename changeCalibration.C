void changeCalibration(TString ifile="calib.root", TString ofile="newcalib.root"){
  TFile inf(ifile);
  TIter nextkey(inf.GetListOfKeys());
  TKey *key;

  TFile ouf(ofile,"RECREATE");
  
  while ((key = (TKey*)nextkey())) {
    TGraph *gr = (TGraph*)key->ReadObj();
    TString name = gr->GetName();
    //if(name.BeginsWith("15") && name.Length()>6) continue;
    gr->Write();
    if(name.Contains("off_15189043006")){
      gr->SetName("off_15189043007");
      gr->Write();
    }
    if(name.Contains("tof_15189043006")){
      gr->SetName("tof_15189043007");
      gr->Write();
    }
    if(name.Contains("off_15186160145")){
      gr->SetName("off_15186160146");
      gr->Write();
    }
    if(name.Contains("tof_15186160145")){
      gr->SetName("tof_15186160146");
      gr->Write();
    }
  }
  inf.Close();
  ouf.Write();
  ouf.Close();
}
