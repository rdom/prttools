#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include "TString.h"

class DataInfo {

  TString _runId;
  Int_t _studyId;
  Int_t _lensId;
  Double_t _angle;
  Double_t _z;
  Double_t _x;
  Double_t _step;

  Int_t _fileId;
  Int_t _zId;
  Int_t _xId;
  Int_t _stepId;
  Int_t _nchildren;
  TString _aliasId;
  std::vector<TString> _childRuns;

public:
  DataInfo(){}; 	//the default constructor
  DataInfo(Int_t studyId, TString r, Int_t l, Double_t a, Double_t z,Double_t x,Double_t s):
    _studyId(studyId),_runId(r),_lensId(l),_angle(a),_z(z),_x(x),_step(s),_aliasId(""),_nchildren(0),_fileId(0){
  };

  ~DataInfo() {}
  
  friend std::ostream& operator<<(std::ostream& os, const DataInfo& f){
    os<<f._runId<<"  "<<f._angle;
    return os;
  };

  bool operator == (const DataInfo& d) const{
    return _lensId == d._lensId && _angle == d._angle && _z == d._z && _x == d._x && _step == d._step;
  }

  bool operator < (const DataInfo& d) const{
    if(_studyId==0 && _angle<d._angle) return true; //angle
    if(_studyId==1 && _z<d._z) return true; //z
    if(_studyId==4 && _angle<d._angle) return true; //angle
    if(_studyId==5 && _z>d._z) return true; //x
    return false; 
  }

  void addChildRunId(TString runid){
    bool exists = false;
    for(Int_t i=0; i<_nchildren;i++ ){
      if(_childRuns[i]==runid) exists=true;
    }

    if(!exists){
      _childRuns.push_back(runid);
      _nchildren++;
    }
  }

  TString info(){
    TString info = Form("Study Id = %d;",_studyId);
    info += "Alias Id = "+_aliasId+";";
    for(Int_t i=0; i<_nchildren;i++){
      info += Form("Child[%d] Id = ",i)+_childRuns[i] +";";
    }
    info += Form("Lens Id = %d;",_lensId);
    info += Form("Angle = %f;",_angle);
    info += Form("X = %f [mm];",_x);
    info += Form("Z = %f [mm];",_z);
    info += Form("Step = %f [mm];",_step);
    info += Form("Uniq file id for current study = %d;",_fileId);
    info += Form("Folder path = %d/%d;",_studyId,_fileId);
    return info;
  }

  /* Accessors */
  Int_t getStudyId() const { return _studyId; }
  TString getRunId() const { return _runId; }
  Int_t getLensId() const {return _lensId; }
  Double_t getAngle() const { return _angle; }
  Double_t getZ() const { return _z; }
  Double_t getX() const { return _x; }
  Double_t getStep() const { return _step; }
  TString getAliasId() const {return _aliasId; }
  Int_t getFileId(){return _fileId;}
  TString getChildRunId(Int_t ind){return _childRuns[ind];}
  Int_t getNChildren(){return _nchildren;}
  
  /* Mutators */
  void setRunId(TString var) { _runId = var; }
  void setAliasId(TString var) { _aliasId = var; }
  void setFileId(Int_t var) { _fileId = var; }
  
};	// end of class DataInfo

Int_t gg_alias=0;
std::vector<DataInfo> dataArray;
std::vector<DataInfo> aliasArray;
const Int_t gg_nstudies = 6;
TString study[gg_nstudies];
void init(){
  study[0]="Angle scan. July 14.";
  study[1]="Z scan. July 14.";
  study[2]="Uniq Angle and Z. July 14.";
  study[3]="Uniq Angle and Z. July 14.";
  study[4]="Angle scan with lens. July 14.";
  study[5]="X scan with lens. July 14.";

  //= July 14 ===================================

  dataArray.push_back(DataInfo(0,"14190005337",0,122.00,442.46,-1.00,30.6));
  dataArray.push_back(DataInfo(0,"14191201517",0,122.00,442.46,-1.00,30.6));
  dataArray.push_back(DataInfo(0,"14191211929",0,122.00,442.46,-1.00,30.6));
  dataArray.push_back(DataInfo(0,"14192221441",0,122.00,442.46,-1.00,30.6));
  dataArray.push_back(DataInfo(0,"14192224720",0,122.00,442.46,-1.00,30.6));
  dataArray.push_back(DataInfo(0,"14192230122",0,120.00,436.57,-1.00,30.6));
  dataArray.push_back(DataInfo(0,"14193200214",0,122.00,442.46,-1.00,30.6));
  dataArray.push_back(DataInfo(0,"14193201620",0,122.00,442.46,-1.00,30.6));
  dataArray.push_back(DataInfo(0,"14193202529",0,122.00,442.46,-1.00,30.6));
  dataArray.push_back(DataInfo(0,"14193203751",0,120.00,436.57,-1.00,30.6));
  dataArray.push_back(DataInfo(0,"14193204418",0,119.50,435.14,-1.00,30.6));
  dataArray.push_back(DataInfo(0,"14193205040",0,119.00,433.72,-1.00,30.6));
  dataArray.push_back(DataInfo(0,"14193205754",0,118.50,432.32,-1.00,30.6));
  dataArray.push_back(DataInfo(0,"14193212008",0,118.00,430.94,-1.00,30.6));
  dataArray.push_back(DataInfo(0,"14193212908",0,116.00,425.55,-1.00,30.6));
  dataArray.push_back(DataInfo(0,"14193213740",0,114.00,420.38,-1.00,30.6));
  dataArray.push_back(DataInfo(0,"14193214623",0,112.00,415.40,-1.00,30.6));
  dataArray.push_back(DataInfo(0,"14193215453",0,110.00,411.36,-1.00,30.6));
  dataArray.push_back(DataInfo(2,"14194193450",0,55.700,1047.4,-1.00,30.6));
  dataArray.push_back(DataInfo(2,"14194193623",0,55.700,1047.4,-1.00,30.6));
  dataArray.push_back(DataInfo(2,"14194193901",0,55.700,1047.4,-1.00,30.6));
  dataArray.push_back(DataInfo(1,"14194200623",0,56.920,1042.1,-1.00,30.6));
  dataArray.push_back(DataInfo(1,"14194200712",0,56.920,1042.1,-1.00,30.6));
  dataArray.push_back(DataInfo(1,"14194202800",0,56.920,677.10,-1.00,30.6));
  dataArray.push_back(DataInfo(1,"14194204509",0,56.920,392.10,-1.00,30.6));
  dataArray.push_back(DataInfo(1,"14194210453",0,56.920,162.10,-1.00,30.6));
  dataArray.push_back(DataInfo(3,"14194212652",0,120.13,180.12,-1.00,30.6));
                                  
  dataArray.push_back(DataInfo(4,"14196211606",1,127.00,291.69,-1.00,18.0));
  dataArray.push_back(DataInfo(4,"14196212144",1,127.00,291.69,-1.00,18.0));
  dataArray.push_back(DataInfo(4,"14196213127",1,126.00,290.21,-1.00,18.0));
  dataArray.push_back(DataInfo(4,"14196214102",1,125.00,288.73,-1.00,18.0));
  dataArray.push_back(DataInfo(4,"14196215015",1,124.00,287.24,-1.00,18.0));
  dataArray.push_back(DataInfo(4,"14196220032",1,123.00,285.74,-1.00,18.0));
  dataArray.push_back(DataInfo(4,"14196220932",1,123.00,285.74,-1.00,18.0));
  dataArray.push_back(DataInfo(4,"14196221241",1,122.00,284.23,-1.00,18.0));
  dataArray.push_back(DataInfo(4,"14196222402",1,121.00,282.71,-1.00,18.0));
  dataArray.push_back(DataInfo(4,"14196222743",1,121.00,282.71,-1.00,18.0));
  dataArray.push_back(DataInfo(4,"14196223412",1,120.00,281.18,-1.00,18.0));
  dataArray.push_back(DataInfo(4,"14196224618",1,119.00,279.65,-1.00,18.0));
  dataArray.push_back(DataInfo(4,"14196225437",1,119.00,279.65,-1.00,18.0));
  dataArray.push_back(DataInfo(4,"14196230221",1,113.00,270.18,-1.00,18.0));
  dataArray.push_back(DataInfo(4,"14196231032",1,112.00,268.56,-1.00,18.0));
  dataArray.push_back(DataInfo(5,"14197203207",1,125.00,318.32,-1.00,18.0));
  dataArray.push_back(DataInfo(5,"14197203627",1,125.00,318.32,-1.00,18.0));
  dataArray.push_back(DataInfo(5,"14197204025",1,125.00,318.32,-1.00,18.0));
  dataArray.push_back(DataInfo(5,"14197204224",1,125.00,318.32,-1.00,18.0));
  dataArray.push_back(DataInfo(5,"14197205430",1,125.00,318.32,-1.00,18.0));
  dataArray.push_back(DataInfo(5,"14197205528",1,125.00,318.32,-21.0,18.0));
  dataArray.push_back(DataInfo(5,"14197210005",1,125.00,318.32,-21.0,18.0));
  dataArray.push_back(DataInfo(5,"14197211125",1,125.00,318.32,-41.0,18.0));
  dataArray.push_back(DataInfo(5,"14197211314",1,125.00,318.32,-41.0,18.0));
  dataArray.push_back(DataInfo(5,"14197211528",1,125.00,318.32,-41.0,18.0));
  dataArray.push_back(DataInfo(5,"14197211648",1,125.00,318.32,-41.0,18.0));
  dataArray.push_back(DataInfo(5,"14197212345",1,125.00,318.32,-61.0,18.0));
  dataArray.push_back(DataInfo(5,"14197212900",1,125.00,318.32,-61.0,18.0));



  //= August 14 =================================
}

std::vector<DataInfo> getStudy(Int_t id){
  std::vector<DataInfo> newset;
  for(UInt_t i = 0; i != aliasArray.size(); i++) {
    if(aliasArray[i].getStudyId()==id){
      newset.push_back(aliasArray[i]);
    }
  }
  std::sort(newset.begin(), newset.end());
  return  newset;
}

void createAliases(){
  TString alias="a";
  for(UInt_t i = 0; i != dataArray.size(); i++) {
    TString same = dataArray[i].getRunId();
    if(dataArray[i].getAliasId()==""){
      gg_alias++;
      dataArray[i].setAliasId(Form("%011d",gg_alias));
    }
    
    for(UInt_t j = 0; j != dataArray.size(); j++) {
      if(dataArray[i] == dataArray[j] && dataArray[i].getRunId() != dataArray[j].getRunId()){
	dataArray[j].setAliasId(dataArray[i].getAliasId());
       	same += " "+dataArray[j].getRunId();
      }
    }
  }

  for(UInt_t i = 0; i != dataArray.size(); i++) {
    TString aid = dataArray[i].getAliasId();
    for(UInt_t j = 0; j != dataArray.size(); j++) {
      if(aid==dataArray[j].getAliasId()){
	bool found = false;
	for(UInt_t k = 0; k != aliasArray.size(); k++) {
	  if(aliasArray[k].getAliasId()==aid){
	    aliasArray[k].addChildRunId(dataArray[j].getRunId());
	    found = true;
	    break;
	  }
	}
	if(!found){
	  DataInfo newdi = dataArray[i];
	  newdi.addChildRunId(dataArray[j].getRunId());
	  aliasArray.push_back(newdi);
	}
      }
    }
  }
  std::sort(aliasArray.begin(), aliasArray.end());

  for(Int_t i=0; i<gg_nstudies; i++){
    std::vector<DataInfo> newset= getStudy(i);
    for(UInt_t j = 0; j != newset.size(); j++) {
      for(UInt_t k = 0; k != aliasArray.size(); k++) {
	if(aliasArray[k].getAliasId()==newset[j].getAliasId()){
	  aliasArray[k].setFileId(j);
	  break;
	}
      }
    }
  }


}

void p_hadd(){
  for(UInt_t i = 0; i != aliasArray.size(); i++) {
    std::cout<<"hadd  cc"<<aliasArray[i].getAliasId()<< ".hld.root  ";
     for(Int_t j=0; j<aliasArray[i].getNChildren();j++ ){
       cout<<"cc"<<aliasArray[i].getChildRunId(j)<<".hld.root ";
     }
    std::cout<<std::endl;
  }
}

void p_print(std::vector<DataInfo> newset, Int_t format){

  if(format==0){ // file name 
    for(UInt_t i = 0; i != newset.size(); i++) {
      std::cout<<"cc"<<newset[i].getAliasId()<< "C.root"<<std::endl;
    }
  }

  if(format==1){ // path 
    for(UInt_t i = 0; i != newset.size(); i++) {
      std::cout<< newset[i].getStudyId()<<"/"<<i<<std::endl; 
   }
  }

  if(format==2){ // info 
    for(UInt_t i = 0; i != newset.size(); i++) {
      std::cout<<newset[i].info()<<std::endl;
    }
  }

  if(format==3){ // prtdirc
    for(UInt_t i = 0; i != newset.size(); i++) {
      std::cout<<"cc"<<newset[i].getAliasId()<< "S.root  " <<newset[i].getStudyId()<<"/"<<i<<std::endl;
    }
  }
  if(format==10){ // cp
    for(UInt_t i = 0; i != dataArray.size(); i++) {
      std::cout<<"cc"<<dataArray[i].getRunId()<< ".hld  ";
    }
  }

}

void datainfo(Int_t studyId=0, Int_t format = 0){
  init();
  createAliases();

  std::vector<DataInfo> newset= getStudy(studyId);

  p_print(newset, format);


  // std::cout<<"ST"<<studyId<<std::endl;
  // for(UInt_t i = 0; i != newset.size(); i++) {
  //   std::cout<< newset[i].getAngle()<<" ";
  // }  

  //for(UInt_t i = 0; i != aliasArray.size(); i++) {
  //   std::cout<<"A  "<<aliasArray[i].getAliasId()<< "  "<< aliasArray[i].getStudyId() << "  "<< aliasArray[i].getAngle() <<std::endl;
  // }
  // p_hadd();
}

