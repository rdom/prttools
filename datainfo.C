#ifndef PrtDataInfo
#define PrtDataInfo

#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include "TString.h"

class DataInfo {

  TString _runId;
  Int_t _studyId;
  Int_t _radiatorId;
  Int_t _lensId;
  Double_t _angle;
  Double_t _z;
  Double_t _x;
  Double_t _xstep;
  Double_t _ystep;
  Double_t _momentum;
  Double_t _simToffset;

  Int_t _fileId;
  Int_t _zId;
  Int_t _xId;
  Int_t _stepId;
  Int_t _nchildren;
  TString _aliasId;
  std::vector<TString> _childRuns;

public:
  DataInfo(){_studyId=-1;}; 	//the default constructor
  DataInfo(Int_t studyId, TString r, Int_t radiator, Int_t l, Double_t a, Double_t z,Double_t x,Double_t xs,Double_t ys, Double_t m):
    _studyId(studyId),_runId(r),_radiatorId(radiator),_lensId(l),_angle(a),_z(z),_x(x),_xstep(xs),_ystep(ys),_momentum(m),_aliasId(""),_nchildren(0),_fileId(0){
  };

  ~DataInfo() {}
  
  friend std::ostream& operator<<(std::ostream& os, const DataInfo& f){
    os<<f._runId<<"  "<<f._angle;
    return os;
  };

  bool operator == (const DataInfo& d) const{
    return _studyId == d._studyId && _radiatorId == d._radiatorId && _lensId == d._lensId && _angle == d._angle && _z == d._z && _x == d._x && _xstep == d._xstep && _ystep == d._ystep && _momentum == d._momentum;
  }

  bool operator < (const DataInfo& d) const{
    if(_studyId==0 && _angle<d._angle) return true; //angle
    if(_studyId==1 && _z<d._z) return true; //z
    if(_studyId==4 && _angle<d._angle) return true; //angle
    if(_studyId==5 && _z>d._z) return true; //x
    if(_studyId>100 && _studyId<170 &&_angle<d._angle) return true; //angle
    if(_studyId>169 && _studyId<180 && _momentum < d._momentum) return true; //momentum
    if(_studyId>179 && _studyId<190 && _z < d._z) return true; //z
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
    info += Form("Radiator Id = %d;",_radiatorId);
    info += Form("Lens Id = %d;",_lensId);
    info += Form("Angle = %f;",_angle);
    info += Form("X = %f [mm];",_x);
    info += Form("Z = %f [mm];",_z);
    info += Form("X Step = %f [mm];",_xstep);
    info += Form("Y Step = %f [mm];",_ystep);
    info += Form("Momentum = %f [mm];",_momentum);
    info += Form("Uniq file id for current study = %d;",_fileId);
    info += Form("Folder path = %d/%d;",_studyId,_fileId);
    return info;
  }

  /* Accessors */
  Int_t getStudyId() const { return _studyId; }
  TString getRunId() const { return _runId; }
  Int_t getRadiatorId() const {return _radiatorId; }
  Int_t getLensId() const {return _lensId; }
  Double_t getAngle() const { return _angle; }
  Double_t getZ() const { return _z; }
  Double_t getX() const { return _x; }
  Double_t getXstep() const { return _xstep; }
  Double_t getYstep() const { return _ystep; }
  Double_t getMomentum() const { return _momentum; }
  TString getAliasId() const {return _aliasId; }
  Int_t getFileId(){return _fileId;}
  TString getChildRunId(Int_t ind){return _childRuns[ind];}
  Int_t getNChildren(){return _nchildren;}
  Double_t getSimTO() const { return _simToffset; }
  
  /* Mutators */
  void setRunId(TString var) { _runId = var; }
  void setAliasId(TString var) { _aliasId = var; }
  void setFileId(Int_t var) { _fileId = var; }
  void setSimTO(Double_t var) { _simToffset = var; }
  
};	// end of class DataInfo

Int_t gg_alias=0;
std::vector<DataInfo> dataArray;
std::vector<DataInfo> aliasArray;
const Int_t gg_nstudies = 200;
Int_t gg_studyArray[gg_nstudies];
TString study[gg_nstudies];
void datainfo_init(){

  for(Int_t i=0; i<gg_nstudies; i++) gg_studyArray[i]=0;
  
  study[0]="Angle scan. July 14.";
  study[1]="Z scan. July 14.";
  study[2]="Uniq Angle and Z. July 14.";
  study[3]="Uniq Angle and Z. July 14.";
  study[4]="Angle scan with lens. July 14.";
  study[5]="X scan with lens. July 14.";

  //= July 14 ===================================
  {
    dataArray.push_back(DataInfo(0,"14190005337",1,0,122.00,442.46,85.0,-1.00,30.6,1.9));
    dataArray.push_back(DataInfo(0,"14191201517",1,0,122.00,442.46,85.0,-1.00,30.6,1.9));
    dataArray.push_back(DataInfo(0,"14191211929",1,0,122.00,442.46,85.0,-1.00,30.6,1.9));
    dataArray.push_back(DataInfo(0,"14192221441",1,0,122.00,442.46,85.0,-1.00,30.6,1.9));
    dataArray.push_back(DataInfo(0,"14192224720",1,0,122.00,442.46,85.0,-1.00,30.6,1.9));
    dataArray.push_back(DataInfo(0,"14192230122",1,0,120.00,436.57,85.0,-1.00,30.6,1.9));
    dataArray.push_back(DataInfo(0,"14193200214",1,0,122.00,442.46,85.0,-1.00,30.6,1.9));
    dataArray.push_back(DataInfo(0,"14193201620",1,0,122.00,442.46,85.0,-1.00,30.6,1.9));
    dataArray.push_back(DataInfo(0,"14193202529",1,0,122.00,442.46,85.0,-1.00,30.6,1.9));
    dataArray.push_back(DataInfo(0,"14193203751",1,0,120.00,436.57,85.0,-1.00,30.6,1.9));
    dataArray.push_back(DataInfo(0,"14193204418",1,0,119.50,435.14,85.0,-1.00,30.6,1.9));
    dataArray.push_back(DataInfo(0,"14193205040",1,0,119.00,433.72,85.0,-1.00,30.6,1.9));
    dataArray.push_back(DataInfo(0,"14193205754",1,0,118.50,432.32,85.0,-1.00,30.6,1.9));
    dataArray.push_back(DataInfo(0,"14193212008",1,0,118.00,430.94,85.0,-1.00,30.6,1.9));
    dataArray.push_back(DataInfo(0,"14193212908",1,0,116.00,425.55,85.0,-1.00,30.6,1.9));
    dataArray.push_back(DataInfo(0,"14193213740",1,0,114.00,420.38,85.0,-1.00,30.6,1.9));
    dataArray.push_back(DataInfo(0,"14193214623",1,0,112.00,415.40,85.0,-1.00,30.6,1.9));
    dataArray.push_back(DataInfo(0,"14193215453",1,0,110.00,411.36,85.0,-1.00,30.6,1.9));
    dataArray.push_back(DataInfo(2,"14194193450",1,0,55.700,1047.4,85.0,-1.00,30.6,1.9));
    dataArray.push_back(DataInfo(2,"14194193623",1,0,55.700,1047.4,85.0,-1.00,30.6,1.9));
    dataArray.push_back(DataInfo(2,"14194193901",1,0,55.700,1047.4,85.0,-1.00,30.6,1.9));
    dataArray.push_back(DataInfo(1,"14194200623",1,0,56.920,1042.1,85.0,-1.00,30.6,1.9));
    dataArray.push_back(DataInfo(1,"14194200712",1,0,56.920,1042.1,85.0,-1.00,30.6,1.9));
    dataArray.push_back(DataInfo(1,"14194202800",1,0,56.920,677.10,85.0,-1.00,30.6,1.9));
    dataArray.push_back(DataInfo(1,"14194204509",1,0,56.920,392.10,85.0,-1.00,30.6,1.9));
    dataArray.push_back(DataInfo(1,"14194210453",1,0,56.920,162.10,85.0,-1.00,30.6,1.9));
    dataArray.push_back(DataInfo(3,"14194212652",1,0,120.13,180.12,85.0,-1.00,30.6,1.9));
                                  
    dataArray.push_back(DataInfo(4,"14196211606",1,2,127.00,291.69,85.0,-1.00,18.0,1.9));
    dataArray.push_back(DataInfo(4,"14196212144",1,2,127.00,291.69,85.0,-1.00,18.0,1.9));
    dataArray.push_back(DataInfo(4,"14196213127",1,2,126.00,290.21,85.0,-1.00,18.0,1.9));
    dataArray.push_back(DataInfo(4,"14196214102",1,2,125.00,288.73,85.0,-1.00,18.0,1.9));
    dataArray.push_back(DataInfo(4,"14196215015",1,2,124.00,287.24,85.0,-1.00,18.0,1.9));
    dataArray.push_back(DataInfo(4,"14196220032",1,2,123.00,285.74,85.0,-1.00,18.0,1.9));
    dataArray.push_back(DataInfo(4,"14196220932",1,2,123.00,285.74,85.0,-1.00,18.0,1.9));
    dataArray.push_back(DataInfo(4,"14196221241",1,2,122.00,284.23,85.0,-1.00,18.0,1.9));
    dataArray.push_back(DataInfo(4,"14196222402",1,2,121.00,282.71,85.0,-1.00,18.0,1.9));
    dataArray.push_back(DataInfo(4,"14196222743",1,2,121.00,282.71,85.0,-1.00,18.0,1.9));
    dataArray.push_back(DataInfo(4,"14196223412",1,2,120.00,281.18,85.0,-1.00,18.0,1.9));
    dataArray.push_back(DataInfo(4,"14196224618",1,2,119.00,279.65,85.0,-1.00,18.0,1.9));
    dataArray.push_back(DataInfo(4,"14196225437",1,2,119.00,279.65,85.0,-1.00,18.0,1.9));
    dataArray.push_back(DataInfo(4,"14196230221",1,2,113.00,270.18,85.0,-1.00,18.0,1.9));
    dataArray.push_back(DataInfo(4,"14196231032",1,2,112.00,268.56,85.0,-1.00,18.0,1.9));
    dataArray.push_back(DataInfo(5,"14197203207",1,2,125.00,318.32,85.0,-1.00,18.0,1.9));
    dataArray.push_back(DataInfo(5,"14197203627",1,2,125.00,318.32,85.0,-1.00,18.0,1.9));
    dataArray.push_back(DataInfo(5,"14197204025",1,2,125.00,318.32,85.0,-1.00,18.0,1.9));
    dataArray.push_back(DataInfo(5,"14197204224",1,2,125.00,318.32,85.0,-1.00,18.0,1.9));
    dataArray.push_back(DataInfo(5,"14197205430",1,2,125.00,318.32,85.0,-1.00,18.0,1.9));
    dataArray.push_back(DataInfo(5,"14197205528",1,2,125.00,318.32,85.0,-21.0,18.0,1.9));
    dataArray.push_back(DataInfo(5,"14197210005",1,2,125.00,318.32,85.0,-21.0,18.0,1.9));
    dataArray.push_back(DataInfo(5,"14197211125",1,2,125.00,318.32,85.0,-41.0,18.0,1.9));
    dataArray.push_back(DataInfo(5,"14197211314",1,2,125.00,318.32,85.0,-41.0,18.0,1.9));
    dataArray.push_back(DataInfo(5,"14197211528",1,2,125.00,318.32,85.0,-41.0,18.0,1.9));
    dataArray.push_back(DataInfo(5,"14197211648",1,2,125.00,318.32,85.0,-41.0,18.0,1.9));
    dataArray.push_back(DataInfo(5,"14197212345",1,2,125.00,318.32,85.0,-61.0,18.0,1.9));
    dataArray.push_back(DataInfo(5,"14197212900",1,2,125.00,318.32,85.0,-61.0,18.0,1.9));
  }

  //= August 14 =================================
  {
    study[50]="X scan with lens. August 14.";
    
    //= May 15    =================================
    study[100]="X scan with lens. May 15.";
  }
  
  //= June 15   =================================
  {
    // lenses
    // 0 - no lens
    // 1 - 2-componet sph. 2CS
    // 2 - 2-componet cyl. 2CC
    // 3 - 3-componet sph. 3CS
    // 4 - 1-componet sph. with air gap 1CS
    // 5 - 1-componet cyl. with air gap 1CC

    // radiators
    // 1 - bar
    // 2 - plate  

    // ======= Angle scans ========================================
    study[150]="Angle scan with 7 GeV/c, bar, 2-layer spherical lens";
    {
      // study id | run name | radiator | lens | angle | z pos | x pos | x step | y step | momentum
      dataArray.push_back(DataInfo(150,"beam_15181014802",1,1,20.00,378,85.0,67.5,0,7.0)); // 11
      dataArray.push_back(DataInfo(150,"beam_15181020846",1,1,25.00,378,85.0,67.5,0,7.0));
      dataArray.push_back(DataInfo(150,"beam_15181012738",1,1,30.00,378,85.0,67.5,0,7.0));
      dataArray.push_back(DataInfo(150,"beam_15181022734",1,1,35.00,378,85.0,67.5,0,7.0));
      dataArray.push_back(DataInfo(150,"beam_15181011013",1,1,40.00,378,85.0,67.5,0,7.0));
      dataArray.push_back(DataInfo(150,"beam_15181024437",1,1,45.00,378,85.0,67.5,0,7.0));
      dataArray.push_back(DataInfo(150,"beam_15181004916",1,1,50.00,378,85.0,67.5,0,7.0));
      dataArray.push_back(DataInfo(150,"beam_15181030400",1,1,55.00,378,85.0,67.5,0,7.0));
      dataArray.push_back(DataInfo(150,"beam_15181002807",1,1,60.00,378,85.0,67.5,0,7.0));
      dataArray.push_back(DataInfo(150,"beam_15181000647",1,1,70.00,378,85.0,67.5,0,7.0));
      dataArray.push_back(DataInfo(150,"beam_15181034220",1,1,75.00,378,85.0,67.5,0,7.0));
      dataArray.push_back(DataInfo(150,"beam_15180233642",1,1,80.00,378,85.0,67.5,0,7.0));
      dataArray.push_back(DataInfo(150,"beam_15181040321",1,1,85.00,378,85.0,67.5,0,7.0));
      dataArray.push_back(DataInfo(150,"beam_15181042444",1,1,89.00,378,85.0,67.5,0,7.0));
      dataArray.push_back(DataInfo(150,"beam_15181044329",1,1,89.50,378,85.0,67.5,0,7.0));
      dataArray.push_back(DataInfo(150,"beam_15181050434",1,1,90.50,378,85.0,67.5,0,7.0));
      // dataArray.push_back(DataInfo(150,"beam_15180231556",1,1,90.50,378,85.0,67.5,0,7.0));
      dataArray.push_back(DataInfo(150,"beam_15181052518",1,1,91.00,378,85.0,67.5,0,7.0));
      dataArray.push_back(DataInfo(150,"beam_15181054550",1,1,95.00,378,85.0,67.5,0,7.0));
      dataArray.push_back(DataInfo(150,"beam_15180225646",1,1,100.0,378,85.0,67.5,0,7.0));
      dataArray.push_back(DataInfo(150,"beam_15181060653",1,1,105.0,378,85.0,67.5,0,7.0));
      dataArray.push_back(DataInfo(150,"beam_15180223816",1,1,110.0,378,85.0,67.5,0,7.0));
      dataArray.push_back(DataInfo(150,"beam_15181062614",1,1,115.0,378,85.0,67.5,0,7.0));
      dataArray.push_back(DataInfo(150,"beam_15180221900",1,1,120.0,378,85.0,67.5,0,7.0));
      dataArray.push_back(DataInfo(150,"beam_15181064133",1,1,124.0,378,85.0,67.5,0,7.0));
      dataArray.push_back(DataInfo(150,"beam_15181065632",1,1,124.5,378,85.0,67.5,0,7.0));
      dataArray.push_back(DataInfo(150,"beam_15181071156",1,1,125.0,378,85.0,67.5,0,7.0));
      dataArray.push_back(DataInfo(150,"beam_15181072725",1,1,125.5,378,85.0,67.5,0,7.0));
      dataArray.push_back(DataInfo(150,"beam_15181074753",1,1,126.0,378,85.0,67.5,0,7.0));
      dataArray.push_back(DataInfo(150,"beam_15180215928",1,1,130.0,378,85.0,67.5,0,7.0));
      dataArray.push_back(DataInfo(150,"beam_15181080832",1,1,135.0,378,85.0,67.5,0,7.0));
      dataArray.push_back(DataInfo(150,"beam_15181082901",1,1,145.0,378,85.0,67.5,0,7.0));
      dataArray.push_back(DataInfo(150,"beam_15180214439",1,1,140.0,378,85.0,67.5,0,7.0));
      dataArray.push_back(DataInfo(150,"beam_15180212654",1,1,150.0,378,85.0,67.5,0,7.0));
      dataArray.push_back(DataInfo(150,"beam_15180211216",1,1,160.0,378,85.0,67.5,0,7.0));
    }
  
    study[151]="Angle scan with 7 GeV/c, bar, 3-layer spherical lens";
    {
      // study id | run name | radiator | lens | angle | z pos | x pos | x step | y step | momentum
      dataArray.push_back(DataInfo(151,"beam_15176224521",1,3,160.0,378,85.0,67.5,0,7.0));  // 69.0 ???
      dataArray.push_back(DataInfo(151,"beam_15176230610",1,3,150.0,378,85.0,67.5,0,7.0));
      dataArray.push_back(DataInfo(151,"beam_15176232840",1,3,140.0,378,85.0,67.5,0,7.0));
      dataArray.push_back(DataInfo(151,"beam_15176234947",1,3,130.0,378,85.0,67.5,0,7.0));
      dataArray.push_back(DataInfo(151,"beam_15177001109",1,3,120.0,378,85.0,67.5,0,7.0));
      dataArray.push_back(DataInfo(151,"beam_15177003141",1,3,110.0,378,85.0,67.5,0,7.0));
      dataArray.push_back(DataInfo(151,"beam_15177005149",1,3,100.0,378,85.0,67.5,0,7.0));
      dataArray.push_back(DataInfo(151,"beam_15177011203",1,3,90.00,378,85.0,67.5,0,7.0));
      dataArray.push_back(DataInfo(151,"beam_15177013222",1,3,80.00,378,85.0,67.5,0,7.0));
      dataArray.push_back(DataInfo(151,"beam_15177015318",1,3,70.00,378,85.0,67.5,0,7.0));
      dataArray.push_back(DataInfo(151,"beam_15177021502",1,3,60.00,378,85.0,67.5,0,7.0));
      dataArray.push_back(DataInfo(151,"beam_15177023715",1,3,50.00,378,85.0,67.5,0,7.0));
      dataArray.push_back(DataInfo(151,"beam_15177025918",1,3,40.00,378,85.0,67.5,0,7.0));
      dataArray.push_back(DataInfo(151,"beam_15177032004",1,3,30.00,378,85.0,67.5,0,7.0));
      dataArray.push_back(DataInfo(151,"beam_15177034013",1,3,20.00,378,85.0,67.5,0,7.0));
      dataArray.push_back(DataInfo(151,"beam_15177040631",1,3,25.00,378,85.0,67.5,0,7.0));
      dataArray.push_back(DataInfo(151,"beam_15177042740",1,3,35.00,378,85.0,67.5,0,7.0));
      dataArray.push_back(DataInfo(151,"beam_15177044731",1,3,45.00,378,85.0,67.5,0,7.0));
      dataArray.push_back(DataInfo(151,"beam_15177050804",1,3,55.00,378,85.0,67.5,0,7.0));
      dataArray.push_back(DataInfo(151,"beam_15177052734",1,3,65.00,378,85.0,67.5,0,7.0));
      dataArray.push_back(DataInfo(151,"beam_15177054709",1,3,75.00,378,85.0,67.5,0,7.0));
      dataArray.push_back(DataInfo(151,"beam_15177060655",1,3,85.00,378,85.0,67.5,0,7.0));
      dataArray.push_back(DataInfo(151,"beam_15177062636",1,3,89.00,378,85.0,67.5,0,7.0));
      dataArray.push_back(DataInfo(151,"beam_15177064354",1,3,89.50,378,85.0,67.5,0,7.0));
      dataArray.push_back(DataInfo(151,"beam_15177065950",1,3,90.50,378,85.0,67.5,0,7.0));
      dataArray.push_back(DataInfo(151,"beam_15177071541",1,3,91.00,378,85.0,67.5,0,7.0));
      dataArray.push_back(DataInfo(151,"beam_15177073150",1,3,95.00,378,85.0,67.5,0,7.0));
      dataArray.push_back(DataInfo(151,"beam_15177075210",1,3,101.0,378,85.0,67.5,0,7.0));
      dataArray.push_back(DataInfo(151,"beam_15177080959",1,3,101.5,378,85.0,67.5,0,7.0));
      dataArray.push_back(DataInfo(151,"beam_15177082744",1,3,102.0,378,85.0,67.5,0,7.0));
      dataArray.push_back(DataInfo(151,"beam_15177085207",1,3,102.5,378,85.0,67.5,0,7.0));
      dataArray.push_back(DataInfo(151,"beam_15177091448",1,3,103.0,378,85.0,67.5,0,7.0));
      // dataArray.push_back(DataInfo(151,"beam_15177092348",1,3,103.0,378,85.0,67.5,0,7.0));
      dataArray.push_back(DataInfo(151,"beam_15177094547",1,3,105.0,378,85.0,67.5,0,7.0));
      dataArray.push_back(DataInfo(151,"beam_15177101200",1,3,112.0,378,85.0,67.5,0,7.0));
      dataArray.push_back(DataInfo(151,"beam_15177103432",1,3,112.5,378,85.0,67.5,0,7.0));
      dataArray.push_back(DataInfo(151,"beam_15177105939",1,3,113.0,378,85.0,67.5,0,7.0));
      dataArray.push_back(DataInfo(151,"beam_15177112552",1,3,114.0,378,85.0,67.5,0,7.0));
      dataArray.push_back(DataInfo(151,"beam_15177114933",1,3,113.5,378,85.0,67.5,0,7.0));
      dataArray.push_back(DataInfo(151,"beam_15177123229",1,3,115.0,378,85.0,67.5,0,7.0));
      dataArray.push_back(DataInfo(151,"beam_15177130959",1,3,124.0,378,85.0,67.5,0,7.0));
      dataArray.push_back(DataInfo(151,"beam_15177133159",1,3,124.5,378,85.0,67.5,0,7.0));
      //dataArray.push_back(DataInfo(151,"beam_15177133258",1,3,124.5,378,85.0,67.5,0,7.0));
      dataArray.push_back(DataInfo(151,"beam_15177135523",1,3,125.0,378,85.0,67.5,0,7.0));
      dataArray.push_back(DataInfo(151,"beam_15177141530",1,3,125.5,378,85.0,67.5,0,7.0));
      dataArray.push_back(DataInfo(151,"beam_15177143523",1,3,126.0,378,85.0,67.5,0,7.0));
      dataArray.push_back(DataInfo(151,"beam_15177145503",1,3,135.0,378,85.0,67.5,0,7.0));
      dataArray.push_back(DataInfo(151,"beam_15177152058",1,3,145.0,378,85.0,67.5,0,7.0));
    }
    
    study[152]="Angle scan with 7 GeV/c, plate, no lens";
    {
      // study id | run name | radiator | lens | angle | z pos | x pos | x step | y step | momentum
      dataArray.push_back(DataInfo(152,"beam_15182215112",2,0,160.0,378,85.0,0.00,0,7.0)); // 11->0 ??? 378 ???
      dataArray.push_back(DataInfo(152,"beam_15182221114",2,0,150.0,378,85.0,0.00,0,7.0));
      dataArray.push_back(DataInfo(152,"beam_15182223629",2,0,140.0,378,85.0,0.00,0,7.0));
      dataArray.push_back(DataInfo(152,"beam_15182225656",2,0,130.0,378,85.0,0.00,0,7.0));
      dataArray.push_back(DataInfo(152,"beam_15182231700",2,0,120.0,378,85.0,0.00,0,7.0));
      dataArray.push_back(DataInfo(152,"beam_15182233553",2,0,110.0,378,85.0,0.00,0,7.0));
      dataArray.push_back(DataInfo(152,"beam_15182235455",2,0,100.0,378,85.0,0.00,0,7.0));
      dataArray.push_back(DataInfo(152,"beam_15183001420",2,0,90.00,378,85.0,0.00,0,7.0));
      dataArray.push_back(DataInfo(152,"beam_15183003118",2,0,80.00,378,85.0,0.00,0,7.0));
      dataArray.push_back(DataInfo(152,"beam_15183004751",2,0,70.00,378,85.0,0.00,0,7.0));
      dataArray.push_back(DataInfo(152,"beam_15183010429",2,0,60.00,378,85.0,0.00,0,7.0));
      dataArray.push_back(DataInfo(152,"beam_15183012049",2,0,50.00,378,85.0,0.00,0,7.0));
      dataArray.push_back(DataInfo(152,"beam_15183013641",2,0,40.00,378,85.0,0.00,0,7.0));
      dataArray.push_back(DataInfo(152,"beam_15183015512",2,0,30.00,378,85.0,0.00,0,7.0));
      dataArray.push_back(DataInfo(152,"beam_15183021251",2,0,20.00,378,85.0,0.00,0,7.0));
      dataArray.push_back(DataInfo(152,"beam_15183022858",2,0,25.00,378,85.0,0.00,0,7.0));
      dataArray.push_back(DataInfo(152,"beam_15183024513",2,0,35.00,378,85.0,0.00,0,7.0));
      dataArray.push_back(DataInfo(152,"beam_15183030053",2,0,45.00,378,85.0,0.00,0,7.0));
      dataArray.push_back(DataInfo(152,"beam_15183032354",2,0,55.00,378,85.0,0.00,0,7.0));
      dataArray.push_back(DataInfo(152,"beam_15183034625",2,0,65.00,378,85.0,0.00,0,7.0));
      // dataArray.push_back(DataInfo(152,"beam_15183042908",2,0,65.00,378,85.0,0.00,0,7.0));
      dataArray.push_back(DataInfo(152,"beam_15183050104",2,0,75.00,378,85.0,0.00,0,7.0));
      dataArray.push_back(DataInfo(152,"beam_15183052450",2,0,85.00,378,85.0,0.00,0,7.0));
      dataArray.push_back(DataInfo(152,"beam_15183054332",2,0,95.00,378,85.0,0.00,0,7.0));
      dataArray.push_back(DataInfo(152,"beam_15183060212",2,0,105.0,378,85.0,0.00,0,7.0));
      dataArray.push_back(DataInfo(152,"beam_15183062001",2,0,115.0,378,85.0,0.00,0,7.0));
      dataArray.push_back(DataInfo(152,"beam_15183063920",2,0,135.0,378,85.0,0.00,0,7.0));
      dataArray.push_back(DataInfo(152,"beam_15183065721",2,0,145.0,378,85.0,0.00,0,7.0));
      dataArray.push_back(DataInfo(152,"beam_15183071509",2,0,155.0,378,85.0,0.00,0,7.0));
    }
  
    study[153]="Angle scan with 7 GeV/c, plate, 2-layer cylindracal lens";
    {
      // study id | run name | radiator | lens | angle | z pos | x pos | x step | y step | momentum
      dataArray.push_back(DataInfo(153,"beam_15184195338",2,2,125.0,378,85.0,0.00,0,7.0)); // 11->0 ??? 378 ???
      dataArray.push_back(DataInfo(153,"beam_15184201540",2,2,160.0,378,85.0,0.00,0,7.0));
      //dataArray.push_back(DataInfo(153,"beam_15184203031",2,2,160.0,378,85.0,0.00,0,7.0));
      dataArray.push_back(DataInfo(153,"beam_15184203911",2,2,150.0,378,85.0,0.00,0,7.0));
      dataArray.push_back(DataInfo(153,"beam_15184205832",2,2,140.0,378,85.0,0.00,0,7.0));
      dataArray.push_back(DataInfo(153,"beam_15184211622",2,2,130.0,378,85.0,0.00,0,7.0));
      dataArray.push_back(DataInfo(153,"beam_15184213201",2,2,120.0,378,85.0,0.00,0,7.0));
      dataArray.push_back(DataInfo(153,"beam_15184214917",2,2,110.0,378,85.0,0.00,0,7.0));
      dataArray.push_back(DataInfo(153,"beam_15184220500",2,2,100.0,378,85.0,0.00,0,7.0));
      dataArray.push_back(DataInfo(153,"beam_15184222053",2,2,90.00,378,85.0,0.00,0,7.0));
      dataArray.push_back(DataInfo(153,"beam_15184223652",2,2,80.00,378,85.0,0.00,0,7.0));
      //dataArray.push_back(DataInfo(153,"beam_15184224645",2,2,80.00,378,85.0,0.00,0,7.0));
      dataArray.push_back(DataInfo(153,"beam_15184230203",2,2,70.00,378,85.0,0.00,0,7.0));
      dataArray.push_back(DataInfo(153,"beam_15184231800",2,2,60.00,378,85.0,0.00,0,7.0));
      dataArray.push_back(DataInfo(153,"beam_15184233347",2,2,50.00,378,85.0,0.00,0,7.0));
      dataArray.push_back(DataInfo(153,"beam_15184234926",2,2,40.00,378,85.0,0.00,0,7.0));
      dataArray.push_back(DataInfo(153,"beam_15185000456",2,2,30.00,378,85.0,0.00,0,7.0));
      // dataArray.push_back(DataInfo(153,"beam_15185003038",2,2,30.00,378,85.0,0.00,0,7.0));
      dataArray.push_back(DataInfo(153,"beam_15185003934",2,2,20.00,378,85.0,0.00,0,7.0));
      // dataArray.push_back(DataInfo(153,"beam_15185005640",2,2,20.00,378,85.0,0.00,0,7.0));
      dataArray.push_back(DataInfo(153,"beam_15185011524",2,2,25.00,378,85.0,0.00,0,7.0));
      dataArray.push_back(DataInfo(153,"beam_15185013725",2,2,35.00,378,85.0,0.00,0,7.0));
      dataArray.push_back(DataInfo(153,"beam_15185015251",2,2,45.00,378,85.0,0.00,0,7.0));
      dataArray.push_back(DataInfo(153,"beam_15185020751",2,2,55.00,378,85.0,0.00,0,7.0));
      dataArray.push_back(DataInfo(153,"beam_15185022252",2,2,65.00,378,85.0,0.00,0,7.0));
      //dataArray.push_back(DataInfo(153,"beam_15185023414",2,2,65.00,378,85.0,0.00,0,7.0));
      dataArray.push_back(DataInfo(153,"beam_15185024910",2,2,75.00,378,85.0,0.00,0,7.0));
      dataArray.push_back(DataInfo(153,"beam_15185030412",2,2,85.00,378,85.0,0.00,0,7.0));
      dataArray.push_back(DataInfo(153,"beam_15185031948",2,2,95.00,378,85.0,0.00,0,7.0));
      dataArray.push_back(DataInfo(153,"beam_15185033448",2,2,105.0,378,85.0,0.00,0,7.0));
      dataArray.push_back(DataInfo(153,"beam_15185035126",2,2,115.0,378,85.0,0.00,0,7.0));
      dataArray.push_back(DataInfo(153,"beam_15185041543",2,2,135.0,378,85.0,0.00,0,7.0));
      dataArray.push_back(DataInfo(153,"beam_15185043932",2,2,145.0,378,85.0,0.00,0,7.0));
      dataArray.push_back(DataInfo(153,"beam_15185050435",2,2,155.0,378,85.0,0.00,0,7.0));
      // dataArray.push_back(DataInfo(153,"beam_15185052036",2,2,155.0,378,85.0,0.00,0,7.0));
    }

    study[154]="Angle scan with 7 GeV/c, bar, 1-layer air gap lens";
    {
      // study id | run name | radiator | lens | angle | z pos | x pos | x step | y step | momentum
      dataArray.push_back(DataInfo(154,"beam_15186160145",1,4,125.0,378,85.0,70.5,15.5,7.0)); //125->125-
      dataArray.push_back(DataInfo(154,"beam_15186162443",1,4,160.0,378,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(154,"beam_15186164323",1,4,150.0,378,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(154,"beam_15186170314",1,4,140.0,378,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(154,"beam_15186172520",1,4,130.0,378,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(154,"beam_15186174759",1,4,120.0,378,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(154,"beam_15186181245",1,4,126.0,378,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(154,"beam_15186183415",1,4,125.5,378,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(154,"beam_15186185602",1,4,124.5,378,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(154,"beam_15186192506",1,4,124.0,378,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(154,"beam_15186194720",1,4,114.0,378,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(154,"beam_15186201432",1,4,114.5,378,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(154,"beam_15186204658",1,4,113.5,378,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(154,"beam_15186211818",1,4,113.0,378,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(154,"beam_15186214322",1,4,112.5,378,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(154,"beam_15186220107",1,4,112.0,378,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(154,"beam_15186221834",1,4,110.0,378,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(154,"beam_15186223441",1,4,103.0,378,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(154,"beam_15186224606",1,4,102.5,378,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(154,"beam_15186230736",1,4,102.0,378,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(154,"beam_15186232426",1,4,101.5,378,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(154,"beam_15186234219",1,4,101.0,378,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(154,"beam_15186235922",1,4,100.0,378,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(154,"beam_15187001242",1,4,90.00,378,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(154,"beam_15187002713",1,4,80.00,378,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(154,"beam_15187004116",1,4,70.00,378,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(154,"beam_15187005522",1,4,60.00,378,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(154,"beam_15187011300",1,4,50.00,378,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(154,"beam_15187012709",1,4,40.00,378,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(154,"beam_15187014015",1,4,30.00,378,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(154,"beam_15187015523",1,4,20.00,378,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(154,"beam_15187021317",1,4,25.00,378,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(154,"beam_15187022833",1,4,35.00,378,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(154,"beam_15187024237",1,4,45.00,378,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(154,"beam_15187025751",1,4,55.00,378,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(154,"beam_15187031101",1,4,65.00,378,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(154,"beam_15187032616",1,4,75.00,378,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(154,"beam_15187034209",1,4,85.00,378,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(154,"beam_15187042652",1,4,95.00,378,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(154,"beam_15187044406",1,4,105.0,378,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(154,"beam_15187050255",1,4,115.0,378,85.0,70.5,15.5,7.0));
      //dataArray.push_back(DataInfo(154,"beam_15187052346",1,4,125.0,378,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(154,"beam_15187054808",1,4,135.0,378,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(154,"beam_15187060856",1,4,145.0,378,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(154,"beam_15187062311",1,4,155.0,378,85.0,70.5,15.5,7.0));
	
    }

    study[155]="Angle scan with 7 GeV/c, bar, 1-layer air gap lens";
    {
      // study id | run name | radiator | lens | angle | z pos | x pos | x step | y step | momentum
      dataArray.push_back(DataInfo(155,"beam_15187131418",1,4,125.0,636,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(155,"beam_15187140809",1,4,90.0,636,85.0,70.5,15.5,7.0));
      //dataArray.push_back(DataInfo(155,"beam_15187142755",1,4,90.0,636,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(155,"beam_15187143608",1,4,20.0,636,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(155,"beam_15187143706",1,4,30.0,636,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(155,"beam_15187143755",1,4,40.0,636,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(155,"beam_15187143930",1,4,50.0,636,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(155,"beam_15187144033",1,4,60.0,636,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(155,"beam_15187144120",1,4,70.0,636,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(155,"beam_15187144247",1,4,80.0,636,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(155,"beam_15187144421",1,4,85.0,636,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(155,"beam_15187144616",1,4,95.0,636,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(155,"beam_15187144659",1,4,100.0,636,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(155,"beam_15187144815",1,4,110.0,636,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(155,"beam_15187144846",1,4,120.0,636,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(155,"beam_15187144937",1,4,130.0,636,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(155,"beam_15187145029",1,4,140.0,636,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(155,"beam_15187145215",1,4,150.0,636,85.0,70.5,15.5,7.0));
    }

    // study[156]="Angle scan with 7 GeV/c, no lens, 2mm air gap";
    // {
    //   // study id | run name | radiator | lens | angle | z pos | x pos | x step | y step | momentum
    //   dataArray.push_back(DataInfo(156,"beam_15187154308",1,0,90.0,636,85.0,70.5,15.5,7.0));
    //   dataArray.push_back(DataInfo(156,"beam_15187155731",1,0,125.0,636,85.0,70.5,15.5,7.0));
    //   dataArray.push_back(DataInfo(156,"beam_15187162106",1,0,155.0,636,85.0,70.5,15.5,7.0));
    // }
    
    study[156]="Angle scan with 7 GeV/c, bar, no lens";
    {
      // study id | run name | radiator | lens | angle | z pos | x pos | x step | y step | momentum
      dataArray.push_back(DataInfo(156,"beam_15187173355",1,0,90.0,636,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(156,"beam_15187170906",1,0,125.0,636,85.0,70.5,15.5,7.0));
      //dataArray.push_back(DataInfo(156,"beam_15187171441",1,0,125.0,636,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(156,"beam_15187165102",1,0,155.0,636,85.0,70.5,15.5,7.0));
    }
      
    study[157]="Angle scan with 7 GeV/c, bar, 2-layer cyl. lens";
    {
      // study id | run name | radiator | lens | angle | z pos | x pos | x step | y step | momentum
      dataArray.push_back(DataInfo(157,"beam_15187185303",1,2,125.0,636,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(157,"beam_15187191230",1,2,150.0,636,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(157,"beam_15187193150",1,2,140.0,636,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(157,"beam_15187194944",1,2,130.0,636,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(157,"beam_15187200758",1,2,126.0,636,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(157,"beam_15187201920",1,2,125.5,636,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(157,"beam_15187202744",1,2,124.5,636,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(157,"beam_15187220536",1,2,124.0,636,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(157,"beam_15187221949",1,2,120.0,636,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(157,"beam_15187223501",1,2,114.5,636,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(157,"beam_15187224841",1,2,114.0,636,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(157,"beam_15187230235",1,2,113.5,636,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(157,"beam_15187231721",1,2,113.0,636,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(157,"beam_15187233559",1,2,112.5,636,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(157,"beam_15187235252",1,2,112.0,636,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(157,"beam_15188000948",1,2,110.0,636,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(157,"beam_15188002737",1,2,103.0,636,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(157,"beam_15188004551",1,2,102.5,636,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(157,"beam_15188010235",1,2,102.0,636,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(157,"beam_15188011929",1,2,101.5,636,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(157,"beam_15188013606",1,2,101.0,636,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(157,"beam_15188015333",1,2,100.0,636,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(157,"beam_15188021123",1,2, 90.0,636,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(157,"beam_15188023039",1,2, 80.0,636,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(157,"beam_15188024955",1,2, 70.0,636,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(157,"beam_15188030927",1,2, 60.0,636,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(157,"beam_15188032753",1,2, 50.0,636,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(157,"beam_15188034628",1,2, 40.0,636,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(157,"beam_15188040609",1,2, 30.0,636,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(157,"beam_15188042903",1,2, 20.0,636,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(157,"beam_15188044801",1,2, 25.0,636,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(157,"beam_15188050525",1,2, 35.0,636,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(157,"beam_15188052044",1,2, 45.0,636,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(157,"beam_15188053521",1,2, 55.0,636,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(157,"beam_15188054911",1,2, 65.0,636,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(157,"beam_15188060326",1,2, 75.0,636,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(157,"beam_15188061830",1,2, 85.0,636,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(157,"beam_15188063302",1,2, 95.0,636,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(157,"beam_15188064753",1,2,105.0,636,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(157,"beam_15188070224",1,2,115.0,636,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(157,"beam_15188071942",1,2,135.0,636,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(157,"beam_15188073431",1,2,145.0,636,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(157,"beam_15188074918",1,2,155.0,636,85.0,70.5,15.5,7.0));

    }
    
    study[158]="Angle scan with 7 GeV/c, bar, 3-layer spherical lens";
    {
      // study id | run name | radiator | lens | angle | z pos | x pos | x step | y step | momentum
      dataArray.push_back(DataInfo(158,"beam_15188161517",1,3, 90,380,85.0,67.0,17,7.0));
      dataArray.push_back(DataInfo(158,"beam_15188163015",1,3,100,380,85.0,67.0,17,7.0));
      dataArray.push_back(DataInfo(158,"beam_15188164521",1,3,110,380,85.0,67.0,17,7.0));
      dataArray.push_back(DataInfo(158,"beam_15188170049",1,3,120,380,85.0,67.0,17,7.0));
      dataArray.push_back(DataInfo(158,"beam_15188171757",1,3,125,380,85.0,67.0,17,7.0));
      dataArray.push_back(DataInfo(158,"beam_15188173339",1,3,130,380,85.0,67.0,17,7.0));
      dataArray.push_back(DataInfo(158,"beam_15188174820",1,3,140,380,85.0,67.0,17,7.0));
      dataArray.push_back(DataInfo(158,"beam_15188175938",1,3,150,380,85.0,67.0,17,7.0));
      dataArray.push_back(DataInfo(158,"beam_15188181744",1,3, 80,380,85.0,67.0,17,7.0));
      dataArray.push_back(DataInfo(158,"beam_15188182922",1,3, 70,380,85.0,67.0,17,7.0));
      dataArray.push_back(DataInfo(158,"beam_15188184242",1,3, 60,380,85.0,67.0,17,7.0));
      dataArray.push_back(DataInfo(158,"beam_15188185804",1,3, 50,380,85.0,67.0,17,7.0));
      dataArray.push_back(DataInfo(158,"beam_15188191251",1,3, 40,380,85.0,67.0,17,7.0));
      dataArray.push_back(DataInfo(158,"beam_15188192749",1,3, 30,380,85.0,67.0,17,7.0));
      dataArray.push_back(DataInfo(158,"beam_15188194347",1,3, 20,380,85.0,67.0,17,7.0));
    }

    study[159]="Angle scan with 7 GeV/c, bar, no lens";
    {
      // study id | run name | radiator | lens | angle | z pos | x pos | x step | y step | momentum
      dataArray.push_back(DataInfo(159,"beam_15188203142",1,0, 90,380,85.0,67.0,17,7.0));
      dataArray.push_back(DataInfo(159,"beam_15188205503",1,0,150,380,85.0,67.0,17,7.0));
      dataArray.push_back(DataInfo(159,"beam_15188211923",1,0,140,380,85.0,67.0,17,7.0));
      dataArray.push_back(DataInfo(159,"beam_15188213501",1,0,130,380,85.0,67.0,17,7.0));
      dataArray.push_back(DataInfo(159,"beam_15188214935",1,0,120,380,85.0,67.0,17,7.0));
      dataArray.push_back(DataInfo(159,"beam_15188220905",1,0,110,380,85.0,67.0,17,7.0));
      dataArray.push_back(DataInfo(159,"beam_15188223014",1,0,100,380,85.0,67.0,17,7.0));
      dataArray.push_back(DataInfo(159,"beam_15188225123",1,0, 80,380,85.0,67.0,17,7.0));
      dataArray.push_back(DataInfo(159,"beam_15188231010",1,0, 70,380,85.0,67.0,17,7.0));
      dataArray.push_back(DataInfo(159,"beam_15188233003",1,0, 60,380,85.0,67.0,17,7.0));
      dataArray.push_back(DataInfo(159,"beam_15188235006",1,0, 50,380,85.0,67.0,17,7.0));
      dataArray.push_back(DataInfo(159,"beam_15189001041",1,0, 40,380,85.0,67.0,17,7.0));
      dataArray.push_back(DataInfo(159,"beam_15189003008",1,0, 30,380,85.0,67.0,17,7.0));
      dataArray.push_back(DataInfo(159,"beam_15189005042",1,0, 20,380,85.0,67.0,17,7.0));
      dataArray.push_back(DataInfo(159,"beam_15189010958",1,0, 25,380,85.0,67.0,17,7.0));
      dataArray.push_back(DataInfo(159,"beam_15189013001",1,0, 35,380,85.0,67.0,17,7.0));
      dataArray.push_back(DataInfo(159,"beam_15189015019",1,0, 45,380,85.0,67.0,17,7.0));
      dataArray.push_back(DataInfo(159,"beam_15189020959",1,0, 55,380,85.0,67.0,17,7.0));
      dataArray.push_back(DataInfo(159,"beam_15189023009",1,0, 65,380,85.0,67.0,17,7.0));
      dataArray.push_back(DataInfo(159,"beam_15189025012",1,0, 75,380,85.0,67.0,17,7.0));
      dataArray.push_back(DataInfo(159,"beam_15189031006",1,0, 85,380,85.0,67.0,17,7.0));
      dataArray.push_back(DataInfo(159,"beam_15189033937",1,0, 95,380,85.0,67.0,17,7.0));
      dataArray.push_back(DataInfo(159,"beam_15189040011",1,0,105,380,85.0,67.0,17,7.0));
      dataArray.push_back(DataInfo(159,"beam_15189041530",1,0,115,380,85.0,67.0,17,7.0));
      dataArray.push_back(DataInfo(159,"beam_15189043006",1,0,125,380,85.0,67.0,17,7.0));
      dataArray.push_back(DataInfo(159,"beam_15189045215",1,0,135,380,85.0,67.0,17,7.0));
      dataArray.push_back(DataInfo(159,"beam_15189050257",1,0,145,380,85.0,67.0,17,7.0));
      dataArray.push_back(DataInfo(159,"beam_15189051255",1,0,155,380,85.0,67.0,17,7.0));
    }
    
    study[160]="Angle scan with 5 GeV/c, bar, 3-layer spherical lens";
    {
      // study id | run name | radiator | lens | angle | z pos | x pos | x step | y step | momentum
      dataArray.push_back(DataInfo(160,"beam_15177231650",1,3,160.0,378,85.0,67.5,16.5,5.0)); //378 ???  67.6 ??? 
      dataArray.push_back(DataInfo(160,"beam_15177233702",1,3,150.0,378,85.0,67.5,16.5,5.0));
      dataArray.push_back(DataInfo(160,"beam_15177235309",1,3,140.0,378,85.0,67.5,16.5,5.0));
      dataArray.push_back(DataInfo(160,"beam_15178000918",1,3,130.0,378,85.0,67.5,16.5,5.0));
      dataArray.push_back(DataInfo(160,"beam_15178002528",1,3,120.0,378,85.0,67.5,16.5,5.0));
      dataArray.push_back(DataInfo(160,"beam_15178004106",1,3,110.0,378,85.0,67.5,16.5,5.0));
      dataArray.push_back(DataInfo(160,"beam_15178005649",1,3,100.0,378,85.0,67.5,16.5,5.0));
      dataArray.push_back(DataInfo(160,"beam_15178011248",1,3,90.00,378,85.0,67.5,16.5,5.0));
      dataArray.push_back(DataInfo(160,"beam_15178012836",1,3,80.00,378,85.0,67.5,16.5,5.0));
      dataArray.push_back(DataInfo(160,"beam_15178014443",1,3,70.00,378,85.0,67.5,16.5,5.0));
      dataArray.push_back(DataInfo(160,"beam_15178020305",1,3,60.00,378,85.0,67.5,16.5,5.0));
      dataArray.push_back(DataInfo(160,"beam_15178021844",1,3,50.00,378,85.0,67.5,16.5,5.0));
      dataArray.push_back(DataInfo(160,"beam_15178023505",1,3,40.00,378,85.0,67.5,16.5,5.0));
      dataArray.push_back(DataInfo(160,"beam_15178025013",1,3,30.00,378,85.0,67.5,16.5,5.0));
      dataArray.push_back(DataInfo(160,"beam_15178030551",1,3,20.00,378,85.0,67.5,16.5,5.0));
      dataArray.push_back(DataInfo(160,"beam_15178032122",1,3,25.00,378,85.0,67.5,16.5,5.0));
      dataArray.push_back(DataInfo(160,"beam_15178033716",1,3,35.00,378,85.0,67.5,16.5,5.0));
      dataArray.push_back(DataInfo(160,"beam_15178035350",1,3,45.00,378,85.0,67.5,16.5,5.0));
      dataArray.push_back(DataInfo(160,"beam_15178041202",1,3,55.00,378,85.0,67.5,16.5,5.0));
      dataArray.push_back(DataInfo(160,"beam_15178042803",1,3,65.00,378,85.0,67.5,16.5,5.0));
      dataArray.push_back(DataInfo(160,"beam_15178042803",1,3,75.00,378,85.0,67.5,16.5,5.0));
      dataArray.push_back(DataInfo(160,"beam_15178050916",1,3,85.00,378,85.0,67.5,16.5,5.0));
      dataArray.push_back(DataInfo(160,"beam_15178053420",1,3,89.00,378,85.0,67.5,16.5,5.0));
      dataArray.push_back(DataInfo(160,"beam_15178055857",1,3,89.50,378,85.0,67.5,16.5,5.0));
      dataArray.push_back(DataInfo(160,"beam_15178061943",1,3,90.50,378,85.0,67.5,16.5,5.0));
      dataArray.push_back(DataInfo(160,"beam_15178064148",1,3,91.00,378,85.0,67.5,16.5,5.0));
      //dataArray.push_back(DataInfo(160,"beam_15178072643",1,3,125.0,378,85.0,67.5,16.5,5.0));
      dataArray.push_back(DataInfo(160,"beam_15178091807",1,3,95.00,378,85.0,67.5,16.5,5.0));
      dataArray.push_back(DataInfo(160,"beam_15178095254",1,3,101.0,378,85.0,67.5,16.5,5.0));
      dataArray.push_back(DataInfo(160,"beam_15178095254",1,3,112.5,378,85.0,67.5,16.5,5.0));
      dataArray.push_back(DataInfo(160,"beam_15178102306",1,3,101.5,378,85.0,67.5,16.5,5.0));
      dataArray.push_back(DataInfo(160,"beam_15178104722",1,3,102.0,378,85.0,67.5,16.5,5.0));
      dataArray.push_back(DataInfo(160,"beam_15178111130",1,3,102.5,378,85.0,67.5,16.5,5.0));
      dataArray.push_back(DataInfo(160,"beam_15178113513",1,3,103.0,378,85.0,67.5,16.5,5.0));
      dataArray.push_back(DataInfo(160,"beam_15178115721",1,3,105.0,378,85.0,67.5,16.5,5.0));
      dataArray.push_back(DataInfo(160,"beam_15178170915",1,3,112.0,378,85.0,67.5,16.5,5.0));
      dataArray.push_back(DataInfo(160,"beam_15178173130",1,3,112.5,378,85.0,67.5,16.5,5.0));
      dataArray.push_back(DataInfo(160,"beam_15178175341",1,3,113.0,378,85.0,67.5,16.5,5.0));
      dataArray.push_back(DataInfo(160,"beam_15178181502",1,3,114.0,378,85.0,67.5,16.5,5.0));
      dataArray.push_back(DataInfo(160,"beam_15178183715",1,3,113.5,378,85.0,67.5,16.5,5.0));
      dataArray.push_back(DataInfo(160,"beam_15178185737",1,3,115.0,378,85.0,67.5,16.5,5.0));
      dataArray.push_back(DataInfo(160,"beam_15178191905",1,3,124.0,378,85.0,67.5,16.5,5.0));
      dataArray.push_back(DataInfo(160,"beam_15178193924",1,3,124.5,378,85.0,67.5,16.5,5.0));
      dataArray.push_back(DataInfo(160,"beam_15178195923",1,3,125.0,378,85.0,67.5,16.5,5.0));
      dataArray.push_back(DataInfo(160,"beam_15178202056",1,3,125.5,378,85.0,67.5,16.5,5.0));
      dataArray.push_back(DataInfo(160,"beam_15178204222",1,3,126.0,378,85.0,67.5,16.5,5.0));
      dataArray.push_back(DataInfo(160,"beam_15178210315",1,3,135.0,378,85.0,67.5,16.5,5.0));
      dataArray.push_back(DataInfo(160,"beam_15178212446",1,3,145.0,378,85.0,67.5,16.5,5.0));
    }

    study[161]="Angle scan with 5 GeV/c, plate, no lens";
    {
      // study id | run name | radiator | lens | angle | z pos | x pos | x step | y step | momentum
      dataArray.push_back(DataInfo(161,"beam_15183203900",2,0,125.0,378,85.0,0.00,11,5.0)); // 11 ??? 378 (636) ???
      dataArray.push_back(DataInfo(161,"beam_15183210433",2,0,160.0,378,85.0,0.00,11,5.0));
      dataArray.push_back(DataInfo(161,"beam_15183212442",2,0,150.0,378,85.0,0.00,11,5.0));
      dataArray.push_back(DataInfo(161,"beam_15183214525",2,0,140.0,378,85.0,0.00,11,5.0));
      dataArray.push_back(DataInfo(161,"beam_15183221343",2,0,130.0,378,85.0,0.00,11,5.0));
      dataArray.push_back(DataInfo(161,"beam_15183223548",2,0,120.0,378,85.0,0.00,11,5.0));
      dataArray.push_back(DataInfo(161,"beam_15183225356",2,0,110.0,378,85.0,0.00,11,5.0));
      dataArray.push_back(DataInfo(161,"beam_15183232854",2,0,100.0,378,85.0,0.00,11,5.0));
      dataArray.push_back(DataInfo(161,"beam_15183235052",2,0,90.00,378,85.0,0.00,11,5.0));
      dataArray.push_back(DataInfo(161,"beam_15184001051",2,0,80.00,378,85.0,0.00,11,5.0));
      dataArray.push_back(DataInfo(161,"beam_15184003113",2,0,70.00,378,85.0,0.00,11,5.0));
      //dataArray.push_back(DataInfo(161,"beam_15184004359",2,0,70.00,378,85.0,0.00,11,5.0));
      dataArray.push_back(DataInfo(161,"beam_15184010045",2,0,60.00,378,85.0,0.00,11,5.0));
      dataArray.push_back(DataInfo(161,"beam_15184012521",2,0,50.00,378,85.0,0.00,11,5.0));
      dataArray.push_back(DataInfo(161,"beam_15184020111",2,0,40.00,378,85.0,0.00,11,5.0));
      dataArray.push_back(DataInfo(161,"beam_15184021921",2,0,30.00,378,85.0,0.00,11,5.0));
      dataArray.push_back(DataInfo(161,"beam_15184023642",2,0,20.00,378,85.0,0.00,11,5.0));
      dataArray.push_back(DataInfo(161,"beam_15184025508",2,0,25.00,378,85.0,0.00,11,5.0));
      dataArray.push_back(DataInfo(161,"beam_15184031335",2,0,35.00,378,85.0,0.00,11,5.0));
      dataArray.push_back(DataInfo(161,"beam_15184033142",2,0,45.00,378,85.0,0.00,11,5.0));
      dataArray.push_back(DataInfo(161,"beam_15184035005",2,0,55.00,378,85.0,0.00,11,5.0));
      dataArray.push_back(DataInfo(161,"beam_15184040554",2,0,65.00,378,85.0,0.00,11,5.0));
      dataArray.push_back(DataInfo(161,"beam_15184042150",2,0,75.00,378,85.0,0.00,11,5.0));
      dataArray.push_back(DataInfo(161,"beam_15184043751",2,0,85.00,378,85.0,0.00,11,5.0));
      dataArray.push_back(DataInfo(161,"beam_15184045615",2,0,95.00,378,85.0,0.00,11,5.0));
      dataArray.push_back(DataInfo(161,"beam_15184051433",2,0,105.0,378,85.0,0.00,11,5.0));
      dataArray.push_back(DataInfo(161,"beam_15184053229",2,0,115.0,378,85.0,0.00,11,5.0));
      dataArray.push_back(DataInfo(161,"beam_15184055137",2,0,135.0,378,85.0,0.00,11,5.0));
      dataArray.push_back(DataInfo(161,"beam_15184060948",2,0,145.0,378,85.0,0.00,11,5.0));
      dataArray.push_back(DataInfo(161,"beam_15184062758",2,0,155.0,378,85.0,0.00,11,5.0));
    }

    study[162]="Angle scan with 5 GeV/c, plate, 2-layer cylindrical lens";
    {
      // study id | run name | radiator | lens | angle | z pos | x pos | x step | y step | momentum
      dataArray.push_back(DataInfo(162,"beam_15185181752",2,2,125.0,378,85.0,0.00,11,5.0)); // 11 ??? 378 (636) ???
      dataArray.push_back(DataInfo(162,"beam_15185190050",2,2,160.0,378,85.0,0.00,11,5.0));
      dataArray.push_back(DataInfo(162,"beam_15185192720",2,2,150.0,378,85.0,0.00,11,5.0));
      // dataArray.push_back(DataInfo(162,"beam_15185194933",2,2,150.0,378,85.0,0.00,11,5.0));
      dataArray.push_back(DataInfo(162,"beam_15185204824",2,2,140.0,378,85.0,0.00,11,5.0));
      dataArray.push_back(DataInfo(162,"beam_15185210513",2,2,130.0,378,85.0,0.00,11,5.0));
      dataArray.push_back(DataInfo(162,"beam_15185213404",2,2,120.0,378,85.0,0.00,11,5.0));
      dataArray.push_back(DataInfo(162,"beam_15185215050",2,2,110.0,378,85.0,0.00,11,5.0));
      // dataArray.push_back(DataInfo(162,"beam_15185220043",2,2,110.0,378,85.0,0.00,11,5.0));
      dataArray.push_back(DataInfo(162,"beam_15185222028",2,2,100.0,378,85.0,0.00,11,5.0));
      dataArray.push_back(DataInfo(162,"beam_15185224240",2,2,90.00,378,85.0,0.00,11,5.0));
      dataArray.push_back(DataInfo(162,"beam_15185230343",2,2,80.00,378,85.0,0.00,11,5.0));
      dataArray.push_back(DataInfo(162,"beam_15185232501",2,2,70.00,378,85.0,0.00,11,5.0));
      //dataArray.push_back(DataInfo(162,"beam_15185234151",2,2,70.00,378,85.0,0.00,11,5.0));
      dataArray.push_back(DataInfo(162,"beam_15185235741",2,2,60.00,378,85.0,0.00,11,5.0));
      dataArray.push_back(DataInfo(162,"beam_15186001336",2,2,50.00,378,85.0,0.00,11,5.0));
      dataArray.push_back(DataInfo(162,"beam_15186002750",2,2,40.00,378,85.0,0.00,11,5.0));
      // dataArray.push_back(DataInfo(162,"beam_15186004614",2,2,40.00,378,85.0,0.00,11,5.0));
      dataArray.push_back(DataInfo(162,"beam_15186005241",2,2,30.00,378,85.0,0.00,11,5.0));
      // dataArray.push_back(DataInfo(162,"beam_15186010424",2,2,30.00,378,85.0,0.00,11,5.0));
      // dataArray.push_back(DataInfo(162,"beam_15186012740",2,2,30.00,378,85.0,0.00,11,5.0));
      dataArray.push_back(DataInfo(162,"beam_15186013446",2,2,20.00,378,85.0,0.00,11,5.0));
      dataArray.push_back(DataInfo(162,"beam_15186015715",2,2,25.00,378,85.0,0.00,11,5.0));
      // dataArray.push_back(DataInfo(162,"beam_15186021939",2,2,25.00,378,85.0,0.00,11,5.0));
      dataArray.push_back(DataInfo(162,"beam_15186022437",2,2,35.00,378,85.0,0.00,11,5.0));
      // dataArray.push_back(DataInfo(162,"beam_15186024205",2,2,35.00,378,85.0,0.00,11,5.0));
      dataArray.push_back(DataInfo(162,"beam_15186025002",2,2,45.00,378,85.0,0.00,11,5.0));
      // dataArray.push_back(DataInfo(162,"beam_15186025950",2,2,45.00,378,85.0,0.00,11,5.0));
      dataArray.push_back(DataInfo(162,"beam_15186031225",2,2,55.00,378,85.0,0.00,11,5.0));
      dataArray.push_back(DataInfo(162,"beam_15186033935",2,2,65.00,378,85.0,0.00,11,5.0));
      dataArray.push_back(DataInfo(162,"beam_15186042438",2,2,75.00,378,85.0,0.00,11,5.0));
      dataArray.push_back(DataInfo(162,"beam_15186044621",2,2,85.00,378,85.0,0.00,11,5.0));
      dataArray.push_back(DataInfo(162,"beam_15186050706",2,2,95.00,378,85.0,0.00,11,5.0));
      dataArray.push_back(DataInfo(162,"beam_15186052828",2,2,105.00,378,85.0,0.00,11,5.0));
      dataArray.push_back(DataInfo(162,"beam_15186054643",2,2,115.00,378,85.0,0.00,11,5.0));
      dataArray.push_back(DataInfo(162,"beam_15186060212",2,2,135.00,378,85.0,0.00,11,5.0));
      dataArray.push_back(DataInfo(162,"beam_15186061853",2,2,145.00,378,85.0,0.00,11,5.0));
      dataArray.push_back(DataInfo(162,"beam_15186063542",2,2,155.00,378,85.0,0.00,11,5.0));
    }

    // ======= Momentum scans ======================================
    study[170]="Momentum scan, bar,  3-layer sph. lens, 125 degree";
    {
      // study id | run name | radiator | lens | angle | z pos | x pos | x step | y step | momentum
      // 16.5 ??? 378 (600) ???
         
      dataArray.push_back(DataInfo(170,"beam_15179034945",1,3,125.0,378,85.0,67.5,16.5,10.));
      dataArray.push_back(DataInfo(170,"beam_15179041521",1,3,125.0,378,85.0,67.5,16.5,9.0));
      dataArray.push_back(DataInfo(170,"beam_15179043010",1,3,125.0,378,85.0,67.5,16.5,8.0));
      dataArray.push_back(DataInfo(170,"beam_15179044540",1,3,125.0,378,85.0,67.5,16.5,7.0));
      dataArray.push_back(DataInfo(170,"beam_15179050147",1,3,125.0,378,85.0,67.5,16.5,6.0));
      dataArray.push_back(DataInfo(170,"beam_15179051932",1,3,125.0,378,85.0,67.5,16.5,5.0));
      dataArray.push_back(DataInfo(170,"beam_15179054114",1,3,125.0,378,85.0,67.5,16.5,4.0));
      dataArray.push_back(DataInfo(170,"beam_15179023137",1,3,125.0,378,85.0,67.5,16.5,3.0));
      dataArray.push_back(DataInfo(170,"beam_15179061046",1,3,125.0,378,85.0,67.5,16.5,2.0));
    }

    study[171]="Momentum scan, plate, no lens, 125 degree";
    {
      // study id | run name | radiator | lens | angle | z pos | x pos | x step | y step | momentum
      dataArray.push_back(DataInfo(171,"beam_15184102124",2,0,125.0,378,85.0,0.0,11,10.)); //11 ???  378 (636) ???
      dataArray.push_back(DataInfo(171,"beam_15184112613",2,0,125.0,378,85.0,0.0,11,9.0));
      dataArray.push_back(DataInfo(171,"beam_15184110440",2,0,125.0,378,85.0,0.0,11,8.0));
      dataArray.push_back(DataInfo(171,"beam_15184114505",2,0,125.0,378,85.0,0.0,11,7.0));
      dataArray.push_back(DataInfo(171,"beam_15184121116",2,0,125.0,378,85.0,0.0,11,6.0));
      dataArray.push_back(DataInfo(171,"beam_15184123830",2,0,125.0,378,85.0,0.0,11,5.0));
      dataArray.push_back(DataInfo(171,"beam_15184131342",2,0,125.0,378,85.0,0.0,11,4.0));
      dataArray.push_back(DataInfo(171,"beam_15184065128",2,0,125.0,378,85.0,0.0,11,3.0));
    }

    study[172]="Momentum scan, plate, no lens, 55 degree";
    {
      // study id | run name | radiator | lens | angle | z pos | x pos | x step | y step | momentum
      dataArray.push_back(DataInfo(172,"beam_15184141820",2,0,55.0,378,85.0,0.0,11,3.0));//11 ???  378 (636) ???
      dataArray.push_back(DataInfo(172,"beam_15184143916",2,0,55.0,378,85.0,0.0,11,4.0));
      dataArray.push_back(DataInfo(172,"beam_15184150700",2,0,55.0,378,85.0,0.0,11,5.0));
      dataArray.push_back(DataInfo(172,"beam_15184154153",2,0,55.0,378,85.0,0.0,11,6.0));
      dataArray.push_back(DataInfo(172,"beam_15184160844",2,0,55.0,378,85.0,0.0,11,7.0)); 
    }
    
    study[173]="Momentum scan, plate, 2-layer cyl. lens, 125 degree";
    {
      // study id | run name | radiator | lens | angle | z pos | x pos | x step | y step | momentum
      dataArray.push_back(DataInfo(173,"beam_15185054335",2,2,125.0,378,85.0,0.0,11,10.)); //11 ???  378 (636) ???
      dataArray.push_back(DataInfo(173,"beam_15185055820",2,2,125.0,378,85.0,0.0,11,9.0));
      dataArray.push_back(DataInfo(173,"beam_15185061258",2,2,125.0,378,85.0,0.0,11,8.0));
      dataArray.push_back(DataInfo(173,"beam_15184195338",2,2,125.0,378,85.0,0.0,11,7.0)); // form 153
      dataArray.push_back(DataInfo(173,"beam_15185062806",2,2,125.0,378,85.0,0.0,11,6.0));
      dataArray.push_back(DataInfo(173,"beam_15185064716",2,2,125.0,378,85.0,0.0,11,5.0));
      dataArray.push_back(DataInfo(173,"beam_15185071041",2,2,125.0,378,85.0,0.0,11,4.0));
      dataArray.push_back(DataInfo(173,"beam_15185075232",2,2,125.0,378,85.0,0.0,11,3.0));
    }

    study[174]="Momentum scan, bar, 1-layer air gap lens, 125 degree";
    {
      // study id | run name | radiator | lens | angle | z pos | x pos | x step | y step | momentum
      dataArray.push_back(DataInfo(174,"beam_15187064331",1,4,125.0,378,85.0,70.5,15.5,10.0));
      dataArray.push_back(DataInfo(174,"beam_15187071510",1,4,125.0,378,85.0,70.5,15.5,9.00));
      dataArray.push_back(DataInfo(174,"beam_15187075151",1,4,125.0,378,85.0,70.5,15.5,8.00));
      dataArray.push_back(DataInfo(174,"beam_15186160146",1,4,125.0,378,85.0,70.5,15.5,7.00)); // from 154
      dataArray.push_back(DataInfo(174,"beam_15187081410",1,4,125.0,378,85.0,70.5,15.5,6.00));
      dataArray.push_back(DataInfo(174,"beam_15187085029",1,4,125.0,378,85.0,70.5,15.5,5.00));
      dataArray.push_back(DataInfo(174,"beam_15187093204",1,4,125.0,378,85.0,70.5,15.5,4.00));
      dataArray.push_back(DataInfo(174,"beam_15187101807",1,4,125.0,378,85.0,70.5,15.5,3.00));
    }

    study[175]="Momentum scan, plate,  2-layer cyl. lens, 22 degree";
    {
      // study id | run name | radiator | lens | angle | z pos | x pos | x step | y step | momentum
      dataArray.push_back(DataInfo(175,"beam_15186070709",2,2,22.0,378,85.0,0.0,11,4.0));
      dataArray.push_back(DataInfo(175,"beam_15186065352",2,2,22.0,378,85.0,0.0,11,6.0));
    }

    study[176]="Momentum scan, plate, 2-layer cyl. lens, 40 degree";
    {
      // study id | run name | radiator | lens | angle | z pos | x pos | x step | y step | momentum
      dataArray.push_back(DataInfo(176,"beam_15186074729",2,2,40.0,378,85.0,0.0,11,4.0));
      dataArray.push_back(DataInfo(176,"beam_15186073113",2,2,40.0,378,85.0,0.0,11,6.0));
    }

    study[177]="Momentum scan, plate, 2-layer cyl. lens, 140 degree";
    {
      // study id | run name | radiator | lens | angle | z pos | x pos | x step | y step | momentum
      dataArray.push_back(DataInfo(177,"beam_15186082623",2,2,140.0,378,85.0,0.0,11,4.0));
      dataArray.push_back(DataInfo(177,"beam_15186085834",2,2,140.0,378,85.0,0.0,11,6.0));
    }

    study[178]="Momentum scan, plate, 2-layer cyl. lens, 27 degree";
    {
      // study id | run name | radiator | lens | angle | z pos | x pos | x step | y step | momentum
      dataArray.push_back(DataInfo(178,"beam_15186100326",2,2,27.0,378,85.0,0.0,11,1.0));
    }

    study[179]="Momentum scan, bar, no lens, 125 degree, +6.3% BHZ3";
    {
      dataArray.push_back(DataInfo(179,"beam_15189053808",1,0,125.0,378,85.0,67.0,17,10.));
      dataArray.push_back(DataInfo(179,"beam_15189054854",1,0,125.0,378,85.0,67.0,17,9.0));
      dataArray.push_back(DataInfo(179,"beam_15189060129",1,0,125.0,378,85.0,67.0,17,8.0));
      dataArray.push_back(DataInfo(179,"beam_15189043007",1,0,125.0,380,85.0,67.0,17,7.0)); // from 159
      dataArray.push_back(DataInfo(179,"beam_15189063013",1,0,125.0,378,85.0,67.0,17,6.0));
      dataArray.push_back(DataInfo(179,"beam_15189064836",1,0,125.0,378,85.0,67.0,17,5.0));
      dataArray.push_back(DataInfo(179,"beam_15189071238",1,0,125.0,378,85.0,67.0,17,4.0));
      dataArray.push_back(DataInfo(179,"beam_15189073600",1,0,125.0,378,85.0,67.0,17,3.0));
      dataArray.push_back(DataInfo(179,"beam_15189080038",1,0,125.0,378,85.0,67.0,17,2.0));
    }

    // ======= Z-X scans ===========================================
    study[180]="Z scan, bar, 1-layer air gap lens, 90 degree";
    {
      // study id | run name | radiator | lens | angle | z pos | x pos | x step | y step | momentum
      dataArray.push_back(DataInfo(180,"beam_15187110254",1,4,90.0,61 ,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(180,"beam_15187112512",1,4,90.0,858,85.0,70.5,15.5,7.0));
      dataArray.push_back(DataInfo(180,"beam_15187114612",1,4,90.0,618,85.0,70.5,15.5,7.0));
    }

    study[181]="Z scan, bar, 1-layer air gap lens, 125 degree";
    {
      // study id | run name | radiator | lens | angle | z pos | x pos | x step | y step | momentum
      dataArray.push_back(DataInfo(181,"beam_15187121519",1,4,125,930 ,85.0,70.5,15.5,7.0));
    }
  }
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
      Int_t st = dataArray[i].getStudyId();
      gg_studyArray[st]++;
      dataArray[i].setAliasId(Form("%011d",st*10000000 +gg_studyArray[st]));
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

  for(UInt_t k = 0; k != aliasArray.size(); k++) {
    Double_t offset =0;
    if("00000000001" == aliasArray[k].getAliasId()) offset = 72.080000;
    if("00000000002" == aliasArray[k].getAliasId()) offset = 72.240000;
    if("00000000003" == aliasArray[k].getAliasId()) offset = 72.240000;
    if("00000000004" == aliasArray[k].getAliasId()) offset = 72.360000;
    if("00000000005" == aliasArray[k].getAliasId()) offset = 72.400000;
    if("00000000006" == aliasArray[k].getAliasId()) offset = 72.480000;
    if("00000000007" == aliasArray[k].getAliasId()) offset = 72.200000;
    if("00000000008" == aliasArray[k].getAliasId()) offset = 72.240000;
    if("00000000009" == aliasArray[k].getAliasId()) offset = 72.160000;
    if("00000000010" == aliasArray[k].getAliasId()) offset = 72.160000;
    if("00000000011" == aliasArray[k].getAliasId()) offset = 73.560000;
    if("00000000012" == aliasArray[k].getAliasId()) offset = 73.560000;
    if("00000000013" == aliasArray[k].getAliasId()) offset = 72.960000;
    if("00000000014" == aliasArray[k].getAliasId()) offset = 72.480000;
    if("00000000015" == aliasArray[k].getAliasId()) offset = 72.120000;
    if("00000000016" == aliasArray[k].getAliasId()) offset = 72.680000;
    if("00000000017" == aliasArray[k].getAliasId()) offset = 72.040000;
    if("00000000018" == aliasArray[k].getAliasId()) offset = 72.120000;
    if("00000000019" == aliasArray[k].getAliasId()) offset = 72.040000;
    if("00000000020" == aliasArray[k].getAliasId()) offset = 72.080000;
    if("00000000021" == aliasArray[k].getAliasId()) offset = 72.040000;
    if("00000000022" == aliasArray[k].getAliasId()) offset = 72.080000;
    if("00000000023" == aliasArray[k].getAliasId()) offset = 72.080000;
    if("00000000024" == aliasArray[k].getAliasId()) offset = 72.040000;
    if("00000000025" == aliasArray[k].getAliasId()) offset = 72.080000;
    if("00000000026" == aliasArray[k].getAliasId()) offset = 72.000000;
    if("00000000027" == aliasArray[k].getAliasId()) offset = 72.040000;
    if("00000000028" == aliasArray[k].getAliasId()) offset = 71.960000;
    if("00000000029" == aliasArray[k].getAliasId()) offset = 71.920000;
    if("00000000030" == aliasArray[k].getAliasId()) offset = 71.920000;
    if("00000000031" == aliasArray[k].getAliasId()) offset = 71.960000;
    if("00000000014" == aliasArray[k].getAliasId()) offset = 72.480000;
    if("00000000002" == aliasArray[k].getAliasId()) offset = 72.240000;
    if("00000000022" == aliasArray[k].getAliasId()) offset = 72.080000;
    aliasArray[k].setSimTO(offset);

  }

  

}

void p_hadd(){
  // for(UInt_t i = 0; i != aliasArray.size(); i++) {
  //   std::cout<<"hadd "<<aliasArray[i].getAliasId()<< ".hld.root  ";
  //   for(Int_t j=0; j<aliasArray[i].getNChildren();j++ ){
  //     std::cout<<aliasArray[i].getChildRunId(j)<<".hld.root ";
  //   }
  //   std::cout<<std::endl;
  // }
  for(UInt_t i = 0; i < aliasArray.size(); i++) {
    Int_t sid =aliasArray[i].getStudyId();
    if(aliasArray[i].getNChildren()>1 && sid>=150){
      std::cout<<"hadd "<<sid<<"/"<<aliasArray[i].getAliasId()<< ".hld.root  ";
      for(Int_t j=0; j<aliasArray[i].getNChildren();j++ ){
	std::cout<<sid<<"/"<<aliasArray[i].getChildRunId(j)<<".hld.root ";
      }
      std::cout<<std::endl;
      std::cout<<"mv "<<sid<<"/"<<aliasArray[i].getAliasId()<< ".hld.root  " <<sid<<"/"<<aliasArray[i].getChildRunId(0)<<".hld.root && rm "<<sid<<"/"<<aliasArray[i].getChildRunId(1)<< ".hld.root  "  <<std::endl;
      
    }
  }
  
}

void p_print(std::vector<DataInfo> newset, Int_t format){

  if(format==0){ // file name 
    for(UInt_t i = 0; i != newset.size(); i++) {
      for(Int_t j=0; j<newset[i].getNChildren();j++ ){
	if(j<1) std::cout<<newset[i].getChildRunId(j)<<std::endl;
      }
    }
  }

  if(format==1){ // path 
    for(UInt_t i = 0; i != newset.size(); i++) {
      std::cout<< newset[i].getStudyId()<<"/"<<newset[i].getFileId()<<std::endl; 
    }
  }

  if(format==2){ // info 
    for(UInt_t i = 0; i != newset.size(); i++) {
      std::cout<<newset[i].info()<<std::endl;
    }
  }

  if(format==3){ // prtdirc
    for(UInt_t i = 0; i != newset.size(); i++) {
      std::cout//<<newset[i].getChildRunId(0)<<":  "
	       <<" -p "<< newset[i].getMomentum() <<" -h "<<newset[i].getRadiatorId()
	       <<" -l "<<newset[i].getLensId()
	       <<" -a "<<newset[i].getAngle()<<" -gz "<<newset[i].getZ()
	       <<" -gx "<<newset[i].getX()<<" -gsx "<<newset[i].getXstep()<<" -gsy "<<newset[i].getYstep()
	       <<" -z 10 ";

      if(newset[i].getStudyId()>=150) std::cout<<" -g 2015 -c 2015 "<<std::endl;
      else  std::cout<<" "<<std::endl;
    }
  }

  if(format==4){ // sim out name  
    for(UInt_t i = 0; i != newset.size(); i++) {
      std::cout<<newset[i].getAliasId()<< "S.root"<<std::endl;
    }
  }
  
  if(format==5){ // alias  name 
    for(UInt_t i = 0; i != newset.size(); i++) {
      std::cout<<newset[i].getAliasId()<< "C.root"<<std::endl;
    }
  }
  
  if(format==6){ // reco
    for(UInt_t i = 0; i != newset.size(); i++) {
      std::cout<<"\"/d/proc/jun15/"<<newset[i].getStudyId()<<"\",\""<<newset[i].getChildRunId(0)<<"C.root\","
	       << newset[i].getStudyId() <<","<<i<<","<<newset[i].getMomentum() << ","<<newset[i].getRadiatorId()
	       <<","<<newset[i].getLensId()
	       <<","<<newset[i].getAngle()<<","<<newset[i].getZ()
	       <<","<<newset[i].getX()<<","<<newset[i].getXstep()
	       <<","<<newset[i].getYstep()<<std::endl;
    }
  }
  
  if(format==7){ // reco sim
    for(UInt_t i = 0; i != newset.size(); i++) {
      std::cout<<"\"/d/proc/jun15/"<<newset[i].getStudyId()<<"\",\""<<newset[i].getChildRunId(0)<<"S.root\","
	       << newset[i].getStudyId() <<","<<i<<","<<newset[i].getMomentum() << ","<<newset[i].getRadiatorId()
	       <<","<<newset[i].getLensId()
	       <<","<<newset[i].getAngle()<<","<<newset[i].getZ()
	       <<","<<newset[i].getX()<<","<<newset[i].getXstep()
	       <<","<<newset[i].getYstep()<<std::endl;
    }
  }

  if(format==8){ // js array
    if(newset.size()>0){
      std::cout<<"p["<< newset[0].getStudyId() <<"]=["<<newset.size()<<",";
      for(UInt_t i = 0; i < newset.size(); i++) {
	std::cout<<i;
	if(i<newset.size()-1) std::cout<< ",";
      }
      std::cout<<"]; " <<std::endl;
      
      std::cout<<"m["<< newset[0].getStudyId() <<"]=[";
      for(UInt_t i = 0; i < newset.size(); i++) {
	std::cout<<newset[i].getAngle();
	if(i<newset.size()-1) std::cout<< ",";
      }
      std::cout<<"]; " <<std::endl;
    }
  }
  
  if(format==10){ // cp
    for(UInt_t i = 0; i != newset.size(); i++) {
      for(Int_t j=0; j<newset[i].getNChildren();j++ ){
	std::cout<<newset[i].getChildRunId(j)<< ".hld.root  ";
      }
    }
  }

}

void p_exportinfo(TString name="alias.html"){
  std::ofstream out;
  out.open(name);

  for(Int_t i=0; i<gg_nstudies; i++){
    std::vector<DataInfo> newset= getStudy(i);
    if(i<100 || newset.size()==0) continue;
    
    out<<"<h3>Study id # "<<i<< "  "<<study[i] <<"</h3>\n";

    out<<"<table style=\"width:100%\">\n";
    out<<"<tr>\n";
    //  out<<"<td>Files alias </td>"; 
    out<<"<td>Files *.hld </td>";
    out<<"<td>Study Id  </td>";
    out<<"<td>Rad- iator </td>";
    out<<"<td>Lens      </td>";
    out<<"<td>Angle     </td>";
    out<<"<td>Z         </td>";
    out<<"<td>X         </td>";
    out<<"<td>X step    </td>";
    out<<"<td>Y step    </td>";
    out<<"<td>p [GeV/c]  </td>";
    out<<"<td>File Id   </td>";
    //  out<<"<td>Sim off set</td>";
    out<<"</tr>";

    for(UInt_t k = 0; k != newset.size(); k++) {  
      out<<"<tr>\n";
      // out<<"<td>"<<newset[k].getAliasId()<<"</td>\n"; 
      out<<"<td>\n";
      for(Int_t j=0; j<newset[k].getNChildren();j++ ){
	out<<newset[k].getChildRunId(j)<<"<br>";
      }
      out<<"</td>";

      out<<"<td>"<<newset[k].getStudyId()<<"</td>";
      out<<"<td>"<<newset[k].getRadiatorId()<<"</td>";
      out<<"<td>"<<newset[k].getLensId()<<"</td>";
      out<<"<td>"<<newset[k].getAngle()<<"</td>";
      out<<"<td>"<<newset[k].getZ()<<"</td>";
      out<<"<td>"<<newset[k].getX()<<"</td>";
      out<<"<td>"<<newset[k].getXstep()<<"</td>";
      out<<"<td>"<<newset[k].getYstep()<<"</td>";
      out<<"<td>"<<newset[k].getMomentum()<<"</td>";
      out<<"<td>"<<newset[k].getFileId()<<"</td>";
      //  out<<"<td>"<<newset[k].getSimTO()<<"</td>";
      out<<"</tr>";
    }
    out<<"</table>\n";
  }

  out.close();
  std::cout<<"p_exportinfo  done  " <<std::endl;
  
}

void p_export(TString name="data.info"){
  std::ofstream out;
  out.open(name);

  for(UInt_t k = 0; k != dataArray.size(); k++) {  

    out<<dataArray[k].getRunId()<<"  "
       <<dataArray[k].getAliasId()<<"  "
       <<dataArray[k].getStudyId()<<"  "
       <<dataArray[k].getRadiatorId()<<"  "
       <<dataArray[k].getLensId()<<"  "
       <<dataArray[k].getAngle()<<"  "
       <<dataArray[k].getZ()<<"  "
       <<dataArray[k].getX()<<"  "
       <<dataArray[k].getXstep()<<"  "
       <<dataArray[k].getYstep()<<"  "
       <<dataArray[k].getFileId()<<"  "
       <<dataArray[k].getSimTO()<<"\n";
  }
  out.close();
}

void p_import(TString name="data.info"){
  // ifstream in;
  // in.open(name);
  // while (1) {
  //   //      in >> x >> y >> z;
  //   if (!in.good()) break;
      
  // }
}

DataInfo getDataInfo(TString name){
  datainfo_init();
  createAliases();
  
  for(UInt_t i = 0; i < aliasArray.size(); i++) {
    TString sname = aliasArray[i].getRunId();
    if(sname.Contains(name)) return aliasArray[i];
  }
  return DataInfo();
}

void datainfo(Int_t studyId=0, Int_t format = 0){
  datainfo_init();
  createAliases();

  std::vector<DataInfo> newset= getStudy(studyId);

  p_print(newset, format);

  if(format==8){
    for(Int_t i=studyId; i<200; i++){
      std::vector<DataInfo> newset= getStudy(i);
      if(newset.size()>0){
	p_print(newset, format);
	std::cout<<Form("desc[%d]=\"",i)<<study[i]<<"\";" <<std::endl;
      }
    }
    for(Int_t i=studyId; i<200; i++){
      std::vector<DataInfo> newset= getStudy(i);
      if(newset.size()>0){
	std::cout<<Form("<option value=%d> %d ",i,i)<<study[i]<<"</option>" <<std::endl;
      }
    }
  }

  if(format==10){ //cp all
    for(Int_t i=0; i<200; i++){
      std::vector<DataInfo> newset= getStudy(i);
      if(newset.size()>0){
	std::cout<<"mkdir ../proc/"<<i <<std::endl;
	std::cout<<"cp ";
	p_print(newset, format);
	std::cout<<" ../proc/"<<i<<std::endl;
      }
    }
  }

  if(format==11){ //print all ID starting from studyId
    for(Int_t i=studyId; i<200; i++){
      std::vector<DataInfo> newset= getStudy(i);
      if(newset.size()>0){
	std::cout<<i<<std::endl;
      }
    }
  }
  
  if(format==12) p_exportinfo(); // html

  //p_export();

  // std::cout<<"ST"<<studyId<<std::endl;
  // for(UInt_t i = 0; i != newset.size(); i++) {
  //   std::cout<< newset[i].getX()<<" ";
  // }  

  //for(UInt_t i = 0; i != aliasArray.size(); i++) {
  //   std::cout<<"A  "<<aliasArray[i].getAliasId()<< "  "<< aliasArray[i].getStudyId() << "  "<< aliasArray[i].getAngle() <<std::endl;
  // }
  // p_hadd();
}


#endif
