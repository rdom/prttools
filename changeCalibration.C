void changeCalibration(TString ifile="calib.root", TString ofile="newcalib.root"){
  TFile inf(ifile);
  TIter nextkey(inf.GetListOfKeys());
  TKey *key;

  TFile ouf(ofile,"RECREATE");
  
  while ((key = (TKey*)nextkey())) {
    TGraph *gr = (TGraph*)key->ReadObj();
    TString name = gr->GetName();
    
    // if(name.Contains("15189005042"))continue;
    // if(name.Contains("15189010958"))continue;
    // if(name.Contains("15189003008"))continue;
    // if(name.Contains("15189013001"))continue;
    // if(name.Contains("15189001041"))continue;
    // if(name.Contains("15189015019"))continue;
    // if(name.Contains("15188235006"))continue;
    // if(name.Contains("15189020959"))continue;
    // if(name.Contains("15188233003"))continue;
    // if(name.Contains("15189023009"))continue;
    // if(name.Contains("15188231010"))continue;
    // if(name.Contains("15189025012"))continue;
    // if(name.Contains("15188225123"))continue;
    // if(name.Contains("15189031006"))continue;
    // if(name.Contains("15188203142"))continue;
    // if(name.Contains("15189033937"))continue;
    // if(name.Contains("15188223014"))continue;
    // if(name.Contains("15189040011"))continue;
    // if(name.Contains("15188220905"))continue;
    // if(name.Contains("15189041530"))continue;
    // if(name.Contains("15188214935"))continue;
    // if(name.Contains("15189043006"))continue;
    // if(name.Contains("15188213501"))continue;
    // if(name.Contains("15189045215"))continue;
    // if(name.Contains("15188211923"))continue;
    // if(name.Contains("15189050257"))continue;
    // if(name.Contains("15188205503"))continue;
    // if(name.Contains("15189051255"))continue;
   
    if(name.Contains("off_")) continue;
    
    // if(name.Contains("off_15189043006")){
    //   gr->SetName("off_15189043007");
    //   gr->Write();
    // }
    // if(name.Contains("tof_15189043006")){
    //   gr->SetName("tof_15189043007");
    //   gr->Write();
    // }
    // if(name.Contains("off_15186160145")){
    //   gr->SetName("off_15186160146");
    //   gr->Write();
    // }
    // if(name.Contains("tof_15186160145")){
    //   gr->SetName("tof_15186160146");
    //   gr->Write();
    // }

    gr->Write();
  }
  inf.Close();
  ouf.Write();
  ouf.Close();
}
