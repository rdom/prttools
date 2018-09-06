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
  Double_t _phi;
  Double_t _test;
  Double_t _z;
  Double_t _x;
  Double_t _xstep;
  Double_t _ystep;
  Double_t _momentum;
  Double_t _simToffset;
  Double_t _beamDimension;

  Int_t _fileId;
  Int_t _zId;
  Int_t _xId;
  Int_t _stepId;
  Int_t _nchildren;
  TString _aliasId;
  std::vector<TString> _childRuns;

public:
  DataInfo(){_studyId=-1;}; 	//the default constructor
  DataInfo(Int_t studyId, TString r, Int_t radiator, Int_t l, Double_t a, Double_t z,Double_t x,Double_t xs,Double_t ys, Double_t m, Double_t beamDimension=10, Double_t simToffset=0, Double_t phi=0, Double_t test=0):
    _studyId(studyId),_runId(r),_radiatorId(radiator),_lensId(l),_angle(a),_z(z),_x(x),_xstep(xs),_ystep(ys),_momentum(m),_aliasId(""),_nchildren(0),_fileId(0),_beamDimension(beamDimension),_simToffset(simToffset),_phi(phi),_test(test){
  };

  ~DataInfo() {}
  
  friend std::ostream& operator<<(std::ostream& os, const DataInfo& f){
    os<<f._runId<<"  "<<f._angle;
    return os;
  };

  bool operator == (const DataInfo& d) const{
    return _studyId == d._studyId && _radiatorId == d._radiatorId && _lensId == d._lensId && _angle == d._angle && _z == d._z && _x == d._x && _xstep == d._xstep && _ystep == d._ystep && _momentum == d._momentum &&_beamDimension == d._beamDimension && _phi == d._phi && _test == d._test;
  }

  bool operator < (const DataInfo& d) const{
    if(_studyId==0 && _angle<d._angle) return true; //angle
    if(_studyId==1 && _z<d._z) return true; //z
    if(_studyId==4 && _angle<d._angle) return true; //angle
    if(_studyId==5 && _z>d._z) return true; //x
    if(_studyId>100 && _studyId<170 &&_angle<d._angle) return true; //angle
    if(_studyId>169 && _studyId<180 && _momentum < d._momentum) return true; //momentum
    if((_studyId==313 || _studyId==322) && _momentum < d._momentum) return true; //momentum
    if(_studyId>179 && _studyId<190 && _z < d._z) return true; //z
    if(_studyId>=200 && _angle<d._angle) return true; //angle
    if(_studyId==314 && _phi<d._phi) return true; //angle
    if(_studyId==404 && _test<d._test) return true; //test value
    if(_studyId==415 && _momentum < d._momentum) return true; //momentum
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
    info += Form("Phi = %f [mm];",_phi);
    info += Form("X = %f [mm];",_x);
    info += Form("Z = %f [mm];",_z);
    info += Form("X Step = %f [mm];",_xstep);
    info += Form("Y Step = %f [mm];",_ystep);
    info += Form("Momentum = %f [mm];",_momentum);
    info += Form("Test = %f [mm];",_test);
    info += Form("Uniq file id for current study = %d;",_fileId);
    info += Form("Folder path = %d/%d;",_studyId,_fileId);
    return info;
  }

  
  TString getOpt(){
    
    TString sopt =  Form(" -h %d -l %d -p %2.2f -a %2.2f -phi %2.2f -gz %2.2f -gx %2.2f -gsx %2.2f -gsy %2.2f -z %2.2f",
			 _radiatorId,_lensId,_momentum,_angle,_phi,_z,_x,_xstep,_ystep,_beamDimension);

    if(_studyId<200) sopt +=" -g 2015 -c 2015 ";
    else if(_studyId==252) sopt += " -g 2021 -c 2021 ";
    else if(_studyId>399) sopt += " -g 2018 -c 2018 ";
    else if(_studyId<300) sopt += " -g 2016 -c 2016 ";
    else if(_studyId>299) sopt += " -g 2017 -c 2017 ";

    if(_studyId>399) sopt += Form("-study %d ",_studyId);
    
    return sopt;
  }

  TString getAlias(){
    return Form("S%d_A%2.2f",_studyId,_angle);
  }

  /* Accessors */
  Int_t getStudyId() const { return _studyId; }
  TString getRunId() const { return _runId; }
  Int_t getRadiatorId() const {return _radiatorId; }
  Int_t getLensId() const {return _lensId; }
  Double_t getAngle() const { return _angle; }
  Double_t getPhi() const { return _phi; }
  Double_t getTest() const { return _test; }
  Double_t getZ() const { return _z; }
  Double_t getX() const { return _x; }
  Double_t getXstep() const { return _xstep; }
  Double_t getYstep() const { return _ystep; }
  Double_t getMomentum() const { return _momentum; }
  Double_t getBeamDimension() const { return _beamDimension; }
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
  
};

Int_t gg_alias=0;
std::vector<DataInfo> dataArray;
std::vector<DataInfo> aliasArray;
const Int_t gg_nstudies = 500;
Int_t gg_studyArray[gg_nstudies];
TString study[gg_nstudies];

void datainfo_init(){
  
  for(Int_t i=0; i<gg_nstudies; i++) gg_studyArray[i]=0;
  
  //= Aug 2017 =================================
  {
    study[300]="Angle scan, bar, 3LS lens, cookies for lens, air gap for the prizm-MCP";
    {
      Double_t o =  28.81;
      // study id | run name | radiator | lens | angle | z pos | x pos | x step | y step | momentum | beam dim | sim offset
      dataArray.push_back(DataInfo(300,"beam_17240151314",1,3,140.0,447,85.0,70.00,5,7.0,10,o));
      dataArray.push_back(DataInfo(300,"beam_17240181058",1,3,130.0,447,85.0,70.00,5,7.0,10,o));
      dataArray.push_back(DataInfo(300,"beam_17240182858",1,3,120.0,447,85.0,70.00,5,7.0,10,o));
      dataArray.push_back(DataInfo(300,"beam_17240185318",1,3,110.0,447,85.0,70.00,5,7.0,10,o));
      dataArray.push_back(DataInfo(300,"beam_17240191244",1,3,100.0,447,85.0,70.00,5,7.0,10,o));
      dataArray.push_back(DataInfo(300,"beam_17240201655",1,3, 90.0,447,85.0,70.00,5,7.0,10,o));
      dataArray.push_back(DataInfo(300,"beam_17240203001",1,3, 80.0,447,85.0,70.00,5,7.0,10,o));
      dataArray.push_back(DataInfo(300,"beam_17240204445",1,3, 70.0,447,85.0,70.00,5,7.0,10,o));
      dataArray.push_back(DataInfo(300,"beam_17240205753",1,3, 60.0,447,85.0,70.00,5,7.0,10,o));
      dataArray.push_back(DataInfo(300,"beam_17240211057",1,3, 50.0,447,85.0,70.00,5,7.0,10,o));
      dataArray.push_back(DataInfo(300,"beam_17240212322",1,3, 40.0,447,85.0,70.00,5,7.0,10,o));
      dataArray.push_back(DataInfo(300,"beam_17240213623",1,3, 30.0,447,85.0,70.00,5,7.0,10,o));
      dataArray.push_back(DataInfo(300,"beam_17240214847",1,3, 20.0,447,85.0,70.00,5,7.0,10,o));
    }

    study[301]="Angle scan, bar, 3LS lens, cookies for lens, air gap for the prizm-MCP; HV +150V; th offset 500";
    {
      Double_t o =  28.81;
      // study id | run name | radiator | lens | angle | z pos | x pos | x step | y step | momentum | beam dim | sim offset
      dataArray.push_back(DataInfo(301,"beam_17242114304",1,3,140.0,447,85.0,70.00,5,7.0,10,o));
      dataArray.push_back(DataInfo(301,"beam_17242113826",1,3,130.0,447,85.0,70.00,5,7.0,10,o));
      dataArray.push_back(DataInfo(301,"beam_17242105452",1,3,120.0,447,85.0,70.00,5,7.0,10,o));
      dataArray.push_back(DataInfo(301,"beam_17242122658",1,3,110.0,447,85.0,70.00,5,7.0,10,o));
      dataArray.push_back(DataInfo(301,"beam_17242121608",1,3,100.0,447,85.0,70.00,5,7.0,10,o));
      dataArray.push_back(DataInfo(301,"beam_17242120250",1,3, 90.0,447,85.0,70.00,5,7.0,10,o));
      dataArray.push_back(DataInfo(301,"beam_17242115300",1,3, 80.0,447,85.0,70.00,5,7.0,10,o));
      dataArray.push_back(DataInfo(301,"beam_17242070844",1,3, 70.0,447,85.0,70.00,5,7.0,10,o));
      dataArray.push_back(DataInfo(301,"beam_17241235836",1,3, 60.0,447,85.0,70.00,5,7.0,10,o));
      dataArray.push_back(DataInfo(301,"beam_17241232404",1,3, 50.0,447,85.0,70.00,5,7.0,10,o));
      dataArray.push_back(DataInfo(301,"beam_17241230615",1,3, 40.0,447,85.0,70.00,5,7.0,10,o));
      dataArray.push_back(DataInfo(301,"beam_17241224638",1,3, 30.0,447,85.0,70.00,5,7.0,10,o));
      dataArray.push_back(DataInfo(301,"beam_17242053631",1,3, 25.0,447,85.0,70.00,5,7.0,10,o));
      dataArray.push_back(DataInfo(301,"beam_17241222625",1,3, 20.0,447,85.0,70.00,5,7.0,10,o));
    }
    
    study[310]="Angle scan, bar, 3LS lens, all with grease; HV +150V; th offset 600";
    {
      Double_t o =  28.81;
      // study id | run name | radiator | lens | angle | z pos | x pos | x step | y step | momentum | beam dim | sim offset
      dataArray.push_back(DataInfo(310,"beam_17244170154",1,3,20.0,447,85.0,70.00,5,7.0,10,o));
      dataArray.push_back(DataInfo(310,"beam_17244171356",1,3,30.0,447,85.0,70.00,5,7.0,10,o));
      dataArray.push_back(DataInfo(310,"beam_17244171854",1,3,40.0,447,85.0,70.00,5,7.0,10,o));
      dataArray.push_back(DataInfo(310,"beam_17244172340",1,3,50.0,447,85.0,70.00,5,7.0,10,o));
      dataArray.push_back(DataInfo(310,"beam_17244172834",1,3,60.0,447,85.0,70.00,5,7.0,10,o));
      dataArray.push_back(DataInfo(310,"beam_17244173452",1,3,70.0,447,85.0,70.00,5,7.0,10,o));
      dataArray.push_back(DataInfo(310,"beam_17244174209",1,3,80.0,447,85.0,70.00,5,7.0,10,o));
      dataArray.push_back(DataInfo(310,"beam_17244175018",1,3,90.0,447,85.0,70.00,5,7.0,10,o));
      dataArray.push_back(DataInfo(310,"beam_17244180555",1,3,100.0,447,85.0,70.00,5,7.0,10,o));
      dataArray.push_back(DataInfo(310,"beam_17244182110",1,3,110.0,447,85.0,70.00,5,7.0,10,o));
      dataArray.push_back(DataInfo(310,"beam_17244182818",1,3,120.0,447,85.0,70.00,5,7.0,10,o));
      dataArray.push_back(DataInfo(310,"beam_17244183924",1,3,130.0,447,85.0,70.00,5,7.0,10,o));
      dataArray.push_back(DataInfo(310,"beam_17244184701",1,3,140.0,447,85.0,70.00,5,7.0,10,o));
      dataArray.push_back(DataInfo(310,"beam_17244185431",1,3,150.0,447,85.0,70.00,5,7.0,10,o));
    }

    // improved TOF
    study[311]="Theta scan, bar, 3LS lens, all with grease; HV +150V; th offset 600" ;
    {
      Double_t o = 28.81;
      // study id | run name | radiator | lens | angle | z pos | x pos | x step | y step | momentum | beam dim | sim offset
      dataArray.push_back(DataInfo(311,"beam_s311_20",1,3,20.0,447,85.0,70.00,5,7.0,10,o));
      dataArray.push_back(DataInfo(311,"beam_s311_25",1,3,25.0,447,85.0,70.00,5,7.0,10,o));
      dataArray.push_back(DataInfo(311,"beam_s311_30",1,3,30.0,447,85.0,70.00,5,7.0,10,o));
      dataArray.push_back(DataInfo(311,"beam_s311_50",1,3,50.0,447,85.0,70.00,5,7.0,10,o));
      dataArray.push_back(DataInfo(311,"beam_s311_70",1,3,70.0,447,85.0,70.00,5,7.0,10,o));
      dataArray.push_back(DataInfo(311,"beam_s311_90",1,3,90.0,447,85.0,70.00,5,7.0,10,o));
      dataArray.push_back(DataInfo(311,"beam_s311_110",1,3,110.0,447,85.0,70.00,5,7.0,10,o));
      dataArray.push_back(DataInfo(311,"beam_s311_130",1,3,130.0,447,85.0,70.00,5,7.0,10,o));

      dataArray.push_back(DataInfo(311,"beam_s311_35",1,3,35.0,447,85.0,70.00,5,7.0,10,o));
      dataArray.push_back(DataInfo(311,"beam_s311_40",1,3,40.0,447,85.0,70.00,5,7.0,10,o));
      dataArray.push_back(DataInfo(311,"beam_s311_45",1,3,45.0,447,85.0,70.00,5,7.0,10,o));
      dataArray.push_back(DataInfo(311,"beam_s311_55",1,3,55.0,447,85.0,70.00,5,7.0,10,o));
      dataArray.push_back(DataInfo(311,"beam_s311_60",1,3,60.0,447,85.0,70.00,5,7.0,10,o));
      dataArray.push_back(DataInfo(311,"beam_s311_75",1,3,75.0,447,85.0,70.00,5,7.0,10,o));
      dataArray.push_back(DataInfo(311,"beam_s311_80",1,3,80.0,447,85.0,70.00,5,7.0,10,o));
      dataArray.push_back(DataInfo(311,"beam_s311_85",1,3,85.0,447,85.0,70.00,5,7.0,10,o));
      dataArray.push_back(DataInfo(311,"beam_s311_95",1,3,95.0,447,85.0,70.00,5,7.0,10,o));
      dataArray.push_back(DataInfo(311,"beam_s311_100",1,3,100.0,447,85.0,70.00,5,7.0,10,o));
      dataArray.push_back(DataInfo(311,"beam_s311_105",1,3,105.0,447,85.0,70.00,5,7.0,10,o));
      dataArray.push_back(DataInfo(311,"beam_s311_115",1,3,115.0,447,85.0,70.00,5,7.0,10,o));
      dataArray.push_back(DataInfo(311,"beam_s311_120",1,3,120.0,447,85.0,70.00,5,7.0,10,o));
      dataArray.push_back(DataInfo(311,"beam_s311_125",1,3,125.0,447,85.0,70.00,5,7.0,10,o));
      dataArray.push_back(DataInfo(311,"beam_s311_135",1,3,135.0,447,85.0,70.00,5,7.0,10,o));
      dataArray.push_back(DataInfo(311,"beam_s311_140",1,3,140.0,447,85.0,70.00,5,7.0,10,o));
    }

    study[312]="Theta scan, phi = 10; bar + 3LS lens, all with grease; HV +150V; flipped edges" ;
    {
      Double_t o = 28.81;
      // study id | run name | radiator | lens | angle | z pos | x pos | x step | y step | momentum | beam dim | sim offset | phi
      dataArray.push_back(DataInfo(312,"beam_s312_20",1,3,20.0,447,85.0,70.00,5,7.0,10,o,10));
      dataArray.push_back(DataInfo(312,"beam_s312_25",1,3,25.0,447,85.0,70.00,5,7.0,10,o,10));
      dataArray.push_back(DataInfo(312,"beam_s312_40",1,3,40.0,447,85.0,70.00,5,7.0,10,o,10));
      dataArray.push_back(DataInfo(312,"beam_s312_50",1,3,50.0,447,85.0,70.00,5,7.0,10,o,10));
      dataArray.push_back(DataInfo(312,"beam_s312_60",1,3,60.0,447,85.0,70.00,5,7.0,10,o,10));
      dataArray.push_back(DataInfo(312,"beam_s312_70",1,3,70.0,447,85.0,70.00,5,7.0,10,o,10));
      dataArray.push_back(DataInfo(312,"beam_s312_80",1,3,80.0,447,85.0,70.00,5,7.0,10,o,10));
      dataArray.push_back(DataInfo(312,"beam_s312_90",1,3,90.0,447,85.0,70.00,5,7.0,10,o,10));
      dataArray.push_back(DataInfo(312,"beam_s312_100",1,3,100.0,447,85.0,70.00,5,7.0,10,o,10));
      dataArray.push_back(DataInfo(312,"beam_s312_110",1,3,110.0,447,85.0,70.00,5,7.0,10,o,10));
      dataArray.push_back(DataInfo(312,"beam_s312_120",1,3,120.0,447,85.0,70.00,5,7.0,10,o,10));
      dataArray.push_back(DataInfo(312,"beam_s312_130",1,3,130.0,447,85.0,70.00,5,7.0,10,o,10));
      dataArray.push_back(DataInfo(312,"beam_s312_140",1,3,140.0,447,85.0,70.00,5,7.0,10,o,10));
    }

    study[313]="Mom scan, theta = 25 ; phi = 0; bar + 3LS lens; HV +150V;" ;
    {
      Double_t o = 28.81;
      // study id | run name | radiator | lens | angle | z pos | x pos | x step | y step | momentum | beam dim | sim offset | phi
      dataArray.push_back(DataInfo(313,"beam_s313_25_3",1,3,20.0,447,85.0,70.00,5,3.0,10,o,0));
      dataArray.push_back(DataInfo(313,"beam_s313_25_5",1,3,20.0,447,85.0,70.00,5,5.0,10,o,0));
      dataArray.push_back(DataInfo(313,"beam_s313_25_6",1,3,20.0,447,85.0,70.00,5,6.0,10,o,0));
      dataArray.push_back(DataInfo(313,"beam_s313_25_7",1,3,20.0,447,85.0,70.00,5,7.0,10,o,0));
      dataArray.push_back(DataInfo(313,"beam_s313_25_8",1,3,20.0,447,85.0,70.00,5,8.0,10,o,0));
      dataArray.push_back(DataInfo(313,"beam_s313_25_10",1,3,20.0,447,85.0,70.00,5,10.0,10,o,0));
    }
    
    study[314]="Phi scan, theta = 25 ; bar + 3LS lens; HV +150V;";
    {
      Double_t o = 28.81;
      // study id | run name | radiator | lens | angle | z pos | x pos | x step | y step | momentum | beam dim | sim offset | phi
      dataArray.push_back(DataInfo(314,"beam_s314_25_1",   1,3,25.0,447,85.0,70.00,14.8,7.0,3,o,1));
      dataArray.push_back(DataInfo(314,"beam_s314_25_2.5", 1,3,25.0,447,85.0,70.00,14.8,7.0,3,o,2.5));
      dataArray.push_back(DataInfo(314,"beam_s314_25_5.0", 1,3,25.0,447,85.0,70.00,14.8,7.0,3,o,5.0));
      dataArray.push_back(DataInfo(314,"beam_s314_25_7.5", 1,3,25.0,447,85.0,70.00,14.8,7.0,3,o,7.5));
      dataArray.push_back(DataInfo(314,"beam_s314_25_10",  1,3,25.0,447,85.0,70.00,14.8,7.0,3,o,10));
      dataArray.push_back(DataInfo(314,"beam_s314_25_12.5",1,3,25.0,447,85.0,70.00,14.8,7.0,3,o,12.5));
      dataArray.push_back(DataInfo(314,"beam_s314_25_15",  1,3,25.0,447,85.0,70.00,14.8,7.0,3,o,15));
    }

    study[315]="Theta scan, phi = 5 ; bar + 3LS lens; HV +150V;";
    {
      Double_t o = 28.81;
      // study id | run name | radiator | lens | angle | z pos | x pos | x step | y step | momentum | beam dim | sim offset | phi
      dataArray.push_back(DataInfo(315,"beam_s315_20_5",1,3,20.0,447,85.0,70.00,14.8,7.0,3,o,5));
      dataArray.push_back(DataInfo(315,"beam_s315_25_5",1,3,25.0,447,85.0,70.00,14.8,7.0,3,o,5));
      dataArray.push_back(DataInfo(315,"beam_s315_30_5",1,3,30.0,447,85.0,70.00,14.8,7.0,3,o,5));
      dataArray.push_back(DataInfo(315,"beam_s315_40_5",1,3,40.0,447,85.0,70.00,14.8,7.0,3,o,5));
      dataArray.push_back(DataInfo(315,"beam_s315_50_5",1,3,50.0,447,85.0,70.00,14.8,7.0,3,o,5));
      dataArray.push_back(DataInfo(315,"beam_s315_60_5",1,3,60.0,447,85.0,70.00,14.8,7.0,3,o,5));
      dataArray.push_back(DataInfo(315,"beam_s315_70_5",1,3,70.0,447,85.0,70.00,14.8,7.0,3,o,5));
      dataArray.push_back(DataInfo(315,"beam_s315_80_5",1,3,80.0,447,85.0,70.00,14.8,7.0,3,o,5));
      dataArray.push_back(DataInfo(315,"beam_s315_85_5",1,3,85.0,447,85.0,70.00,14.8,7.0,3,o,5));
      dataArray.push_back(DataInfo(315,"beam_s315_90_5",1,3,90.0,447,85.0,70.00,14.8,7.0,3,o,5));
      dataArray.push_back(DataInfo(315,"beam_s315_95_5",1,3,95.0,447,85.0,70.00,14.8,7.0,3,o,5));
      dataArray.push_back(DataInfo(315,"beam_s315_100_5",1,3,100.0,447,85.0,70.00,14.8,7.0,3,o,5));
      dataArray.push_back(DataInfo(315,"beam_s315_110_5",1,3,110.0,447,85.0,70.00,14.8,7.0,3,o,5));
      dataArray.push_back(DataInfo(315,"beam_s315_120_5",1,3,120.0,447,85.0,70.00,14.8,7.0,3,o,5));
      dataArray.push_back(DataInfo(315,"beam_s315_130_5",1,3,130.0,447,85.0,70.00,14.8,7.0,3,o,5));
      dataArray.push_back(DataInfo(315,"beam_s315_140_5",1,3,140.0,447,85.0,70.00,14.8,7.0,3,o,5));
      dataArray.push_back(DataInfo(315,"beam_s315_145_5",1,3,145.0,447,85.0,70.00,14.8,7.0,3,o,5));      
      dataArray.push_back(DataInfo(315,"beam_s315_150_5",1,3,150.0,447,85.0,70.00,14.8,7.0,3,o,5));
    }
    
    study[316]="Theta scan, phi = 0 ; plate + 2LC lens; HV +150V;";
    {
      Double_t o = 28.81;
      // study id | run name | radiator | lens | angle | z pos | x pos | x step | y step | momentum | beam dim | sim offset | phi
      dataArray.push_back(DataInfo(316,"beam_s316_15",2,2,15.0,447,85.0,70.00,17,7.0,3,o,0));
      dataArray.push_back(DataInfo(316,"beam_s316_20",2,2,20.0,447,85.0,70.00,17,7.0,3,o,0));
      dataArray.push_back(DataInfo(316,"beam_s316_25",2,2,25.0,447,85.0,70.00,17,7.0,3,o,0));
      dataArray.push_back(DataInfo(316,"beam_s316_30",2,2,30.0,447,85.0,70.00,17,7.0,3,o,0));
      dataArray.push_back(DataInfo(316,"beam_s316_40",2,2,40.0,447,85.0,70.00,17,7.0,3,o,0));
      dataArray.push_back(DataInfo(316,"beam_s316_50",2,2,50.0,447,85.0,70.00,17,7.0,3,o,0));
      dataArray.push_back(DataInfo(316,"beam_s316_60",2,2,60.0,447,85.0,70.00,17,7.0,3,o,0));
      dataArray.push_back(DataInfo(316,"beam_s316_70",2,2,70.0,447,85.0,70.00,17,7.0,3,o,0));
      dataArray.push_back(DataInfo(316,"beam_s316_80",2,2,80.0,447,85.0,70.00,17,7.0,3,o,0));
      dataArray.push_back(DataInfo(316,"beam_s316_90",2,2,90.0,447,85.0,70.00,17,7.0,3,o,0));
      dataArray.push_back(DataInfo(316,"beam_s316_100",2,2,100.0,447,85.0,70.00,17,7.0,3,o,0));
      dataArray.push_back(DataInfo(316,"beam_s316_110",2,2,110.0,447,85.0,70.00,17,7.0,3,o,0));
      dataArray.push_back(DataInfo(316,"beam_s316_120",2,2,120.0,447,85.0,70.00,17,7.0,3,o,0));
      dataArray.push_back(DataInfo(316,"beam_s316_130",2,2,130.0,447,85.0,70.00,17,7.0,3,o,0));
      dataArray.push_back(DataInfo(316,"beam_s316_140",2,2,140.0,447,85.0,70.00,17,7.0,3,o,0));

      dataArray.push_back(DataInfo(316,"beam_s316_125",2,2,125.0,447,85.0,70.00,17,7.0,3,o,0));
      dataArray.push_back(DataInfo(316,"beam_s316_55",2,2,55.0,447,85.0,70.00,17,7.0,3,o,0));
      dataArray.push_back(DataInfo(316,"beam_s316_45",2,2,45.0,447,85.0,70.00,17,7.0,3,o,0));
      dataArray.push_back(DataInfo(316,"beam_s316_35",2,2,35.0,447,85.0,70.00,17,7.0,3,o,0));
      dataArray.push_back(DataInfo(316,"beam_s316_65",2,2,65.0,447,85.0,70.00,17,7.0,3,o,0));
      dataArray.push_back(DataInfo(316,"beam_s316_75",2,2,75.0,447,85.0,70.00,17,7.0,3,o,0));
      dataArray.push_back(DataInfo(316,"beam_s316_85",2,2,85.0,447,85.0,70.00,17,7.0,3,o,0));
      dataArray.push_back(DataInfo(316,"beam_s316_95",2,2,95.0,447,85.0,70.00,17,7.0,3,o,0));
      dataArray.push_back(DataInfo(316,"beam_s316_105",2,2,105.0,447,85.0,70.00,17,7.0,3,o,0));
      dataArray.push_back(DataInfo(316,"beam_s316_115",2,2,115.0,447,85.0,70.00,17,7.0,3,o,0));
      dataArray.push_back(DataInfo(316,"beam_s316_135",2,2,135.0,447,85.0,70.00,17,7.0,3,o,0));
    }

    study[317]="Theta scan, phi = 0 ; plate + 3LC lens";
    {
      Double_t o = 28.81;
      // study id | run name | radiator | lens | angle | z pos | x pos | x step | y step | momentum | beam dim | sim offset | phi
      dataArray.push_back(DataInfo(317,"beam_s317_20",2,6,20.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(317,"beam_s317_25",2,6,25.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(317,"beam_s317_30",2,6,30.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(317,"beam_s317_40",2,6,40.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(317,"beam_s317_50",2,6,50.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(317,"beam_s317_55",2,6,55.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(317,"beam_s317_60",2,6,60.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(317,"beam_s317_70",2,6,70.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(317,"beam_s317_80",2,6,80.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(317,"beam_s317_90",2,6,90.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(317,"beam_s317_100",2,6,100.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(317,"beam_s317_110",2,6,110.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(317,"beam_s317_120",2,6,120.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(317,"beam_s317_130",2,6,130.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(317,"beam_s317_140",2,6,140.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(317,"beam_s317_125",2,6,125.0,447,85.0,70.00,0,7.0,3,o,0));
      
    }      

    study[318]="Theta scan, phi = 5 ; plate + 3LC lens; bobbles"; 
    {
      Double_t o = 28.81;
      // study id | run name | radiator | lens | angle | z pos | x pos | x step | y step | momentum | beam dim | sim offset | phi
      dataArray.push_back(DataInfo(316,"beam_s318_20",2,6,20.0,447,85.0,70.00,0,7.0,3,o,5));
    }

    study[319]="z scan, theta = 25; phi = 0 ; plate + 3LC lens";
    {
      Double_t o = 28.81;
      // study id | run name | radiator | lens | angle | z pos | x pos | x step | y step | momentum | beam dim | sim offset | phi
      dataArray.push_back(DataInfo(319,"beam_s319_262",2,6,20.0,262,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(319,"beam_s319_513",2,6,20.0,513,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(319,"beam_s319_668",2,6,20.0,668,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(319,"beam_s319_xm50",2,6,20.0,447,135.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(319,"beam_s319_x25" ,2,6,20.0,447,60.0,70.00,0,7.0,3,o,0));
    }

    study[320]="Theta scan; phi = 5 ; plate + 3LC lens";
    {
      Double_t o = 28.81;
      // study id | run name | radiator | lens | angle | z pos | x pos | x step | y step | momentum | beam dim | sim offset | phi
      dataArray.push_back(DataInfo(319,"beam_s320_20",2,6,20.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(319,"beam_s320_25",2,6,25.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(319,"beam_s320_30",2,6,30.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(319,"beam_s320_40",2,6,40.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(319,"beam_s320_50",2,6,50.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(319,"beam_s320_60",2,6,60.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(319,"beam_s320_70",2,6,70.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(319,"beam_s320_80",2,6,80.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(319,"beam_s320_90",2,6,90.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(319,"beam_s320_100",2,6,100.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(319,"beam_s320_110",2,6,110.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(319,"beam_s320_120",2,6,120.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(319,"beam_s320_130",2,6,130.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(319,"beam_s320_140",2,6,140.0,447,85.0,70.00,0,7.0,3,o,0));
    }

    study[321]="Offset scan; plate + 3LC lens";
    {
      Double_t o = 28.81;
      // study id | run name | radiator | lens | angle | z pos | x pos | x step | y step | momentum | beam dim | sim offset | phi
      dataArray.push_back(DataInfo(321,"beam_s321_25_900",2,6,25.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(321,"beam_s321_25_1200",2,6,25.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(321,"beam_s321_25_2400",2,6,25.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(321,"beam_s321_125_900",2,6,125.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(321,"beam_s321_125_1200",2,6,125.0,447,85.0,70.00,0,7.0,3,o,0));
    }

    study[322]="Momentum scan; theta=25; phi=0; plate + 3LC lens";
    {
      Double_t o = 28.81;
      // study id | run name | radiator | lens | angle | z pos | x pos | x step | y step | momentum | beam dim | sim offset | phi
      dataArray.push_back(DataInfo(322,"beam_s322_10",2,6,25.0,447,85.0,70.00,0,10.0,3,o,0));
      dataArray.push_back(DataInfo(322,"beam_s322_9" ,2,6,25.0,447,85.0,70.00,0, 9.0,3,o,0));
      dataArray.push_back(DataInfo(322,"beam_s322_8" ,2,6,25.0,447,85.0,70.00,0, 8.0,3,o,0));
      dataArray.push_back(DataInfo(322,"beam_s322_7" ,2,6,25.0,447,85.0,70.00,0, 7.0,3,o,0));
      dataArray.push_back(DataInfo(322,"beam_s322_6" ,2,6,25.0,447,85.0,70.00,0, 6.0,3,o,0));
      dataArray.push_back(DataInfo(322,"beam_s322_5" ,2,6,25.0,447,85.0,70.00,0, 5.0,3,o,0));
      dataArray.push_back(DataInfo(322,"beam_s322_4" ,2,6,25.0,447,85.0,70.00,0, 4.0,3,o,0));
      dataArray.push_back(DataInfo(322,"beam_s322_3" ,2,6,25.0,447,85.0,70.00,0, 3.0,3,o,0));
      dataArray.push_back(DataInfo(322,"beam_s322_2" ,2,6,25.0,447,85.0,70.00,0, 2.0,3,o,0));
    }

    study[323]="Theta scan; phi=0; plate + 3LC lens";
    {
      Double_t o = 28.81;
      // study id | run name | radiator | lens | angle | z pos | x pos | x step | y step | momentum | beam dim | sim offset | phi
      dataArray.push_back(DataInfo(323,"beam_17253115922",2,6,20.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(323,"beam_17253114621",2,6,25.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(323,"beam_17253120320",2,6,30.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(323,"beam_17253120717",2,6,35.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(323,"beam_17253121145",2,6,40.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(323,"beam_17253121546",2,6,45.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(323,"beam_17253121937",2,6,50.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(323,"beam_17253122256",2,6,55.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(323,"beam_17253122621",2,6,60.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(323,"beam_17253122927",2,6,65.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(323,"beam_17253123238",2,6,70.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(323,"beam_17253123541",2,6,75.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(323,"beam_17253123831",2,6,80.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(323,"beam_17253124041",2,6,85.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(323,"beam_17253124232",2,6,90.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(323,"beam_17253124413",2,6,95.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(323,"beam_17253124713",2,6,100.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(323,"beam_17253124856",2,6,105.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(323,"beam_17253125035",2,6,110.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(323,"beam_17253125158",2,6,115.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(323,"beam_17253125333",2,6,120.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(323,"beam_17253125511",2,6,125.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(323,"beam_17253125649",2,6,130.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(323,"beam_17253125831",2,6,135.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(323,"beam_17253130011",2,6,140.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(323,"beam_17253130206",2,6,145.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(323,"beam_17253130341",2,6,150.0,447,85.0,70.00,0,7.0,3,o,0));
    }

    study[324]="Colimator scan; theta=25; phi=0; plate + 3LC lens";
    {
      Double_t o = 28.81;
      // study id | run name | radiator | lens | angle | z pos | x pos | x step | y step | momentum | beam dim | sim offset | phi
      dataArray.push_back(DataInfo(324,"beam_s324_5x5",2,6,20.0,447,85.0,70.00,0,7.0,3,o,0));
    }

    study[325]="X edge scan; theta=25; phi=0; plate + cookies + 3LC lens;";
    {
      Double_t o = 28.81;
      // study id | run name | radiator | lens | angle | z pos | x pos | x step | y step | momentum | beam dim | sim offset | phi
      dataArray.push_back(DataInfo(325,"beam_s325_x105",2,6,25.0,447,190.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(325,"beam_s325_x95" ,2,6,25.0,447,180.0,70.00,0,7.0,3,o,0));
    }
    
    study[326]="Theta scan; phi=0; plate + cookies + 3LC lens";
    {
      Double_t o = 28.81;
      // study id | run name | radiator | lens | angle | z pos | x pos | x step | y step | momentum | beam dim | sim offset | phi
      dataArray.push_back(DataInfo(326,"beam_17253142508",2,6,20.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(326,"beam_s326_25"    ,2,6,25.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(326,"beam_17253134837",2,6,30.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(326,"beam_17253135121",2,6,40.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(326,"beam_17253135327",2,6,50.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(326,"beam_17253135538",2,6,60.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(326,"beam_17253135758",2,6,70.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(326,"beam_17253140013",2,6,80.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(326,"beam_17253140220",2,6,90.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(326,"beam_17253140439",2,6,100.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(326,"beam_17253140658",2,6,110.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(326,"beam_17253140950",2,6,120.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(326,"beam_17253141258",2,6,130.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(326,"beam_17253141532",2,6,140.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(326,"beam_17253141820",2,6,150.0,447,85.0,70.00,0,7.0,3,o,0));      
    }

    study[330]="Theta scan; phi=0; bar + grease + 3LC lens"; //Zygo bar #2.
    {
      Double_t o = 28.81;
      // study id | run name | radiator | lens | angle | z pos | x pos | x step | y step | momentum | beam dim | sim offset | phi
      dataArray.push_back(DataInfo(330,"beam_s330_20",1,6,20.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(330,"beam_s330_25",1,6,25.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(330,"beam_s330_30",1,6,30.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(330,"beam_s330_40",1,6,40.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(330,"beam_s330_50",1,6,50.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(330,"beam_s330_60",1,6,60.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(330,"beam_s330_70",1,6,70.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(330,"beam_s330_80",1,6,80.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(330,"beam_s330_90",1,6,90.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(330,"beam_s330_100",1,6,100.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(330,"beam_s330_110",1,6,110.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(330,"beam_s330_120",1,6,120.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(330,"beam_s330_130",1,6,130.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(330,"beam_s330_140",1,6,140.0,447,85.0,70.00,0,7.0,3,o,0));

      dataArray.push_back(DataInfo(330,"beam_s330_35",1,6,35.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(330,"beam_s330_45",1,6,45.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(330,"beam_s330_55",1,6,55.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(330,"beam_s330_65",1,6,65.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(330,"beam_s330_75",1,6,75.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(330,"beam_s330_85",1,6,85.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(330,"beam_s330_95",1,6,95.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(330,"beam_s330_105",1,6,105.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(330,"beam_s330_115",1,6,115.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(330,"beam_s330_125",1,6,125.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(330,"beam_s330_135",1,6,135.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(330,"beam_s330_145",1,6,145.0,447,85.0,70.00,0,7.0,3,o,0));      
    }

    study[331]="Theta scan; phi=0; bar + grease + 3LC lens";  //Zygo bar #2.
    {
      Double_t o = 28.81;      
      // study id | run name | radiator | lens | angle | z pos | x pos | x step | y step | momentum | beam dim | sim offset | phi
      dataArray.push_back(DataInfo(331,"beam_17253233933",1,6,20.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(331,"beam_17254031211",1,6,25.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(331,"beam_17253192643",1,6,30.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(331,"beam_17254094201",1,6,35.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(331,"beam_17253195712",1,6,40.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(331,"beam_17254092959",1,6,45.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(331,"beam_17253202748",1,6,50.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(331,"beam_17254090322",1,6,55.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(331,"beam_17253205819",1,6,60.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(331,"beam_17254095328",1,6,65.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(331,"beam_17253213338",1,6,70.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(331,"beam_17254100438",1,6,75.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(331,"beam_17253220547",1,6,80.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(331,"beam_17254101544",1,6,85.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(331,"beam_17253223758",1,6,90.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(331,"beam_17254102621",1,6,95.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(331,"beam_17253230713",1,6,100.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(331,"beam_17254103710",1,6,105.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(331,"beam_17254063100",1,6,110.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(331,"beam_17254104756",1,6,115.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(331,"beam_17254070617",1,6,120.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(331,"beam_17254083734",1,6,125.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(331,"beam_17254074316",1,6,130.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(331,"beam_17254110232",1,6,135.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(331,"beam_17254081240",1,6,140.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(331,"beam_17254111419",1,6,145.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(331,"beam_17254112735",1,6,150.0,447,85.0,70.00,0,7.0,3,o,0));
    }

    study[332]="Theta scan; phi=0; bar + grease + 3LS lens";
    {
      Double_t o = 28.81;//   12.59 + 62 + 20-63-2.95;
      // study id | run name | radiator | lens | angle | z pos | x pos | x step | y step | momentum | beam dim | sim offset | phi
      dataArray.push_back(DataInfo(332,"beam_s332_20",1,3, 20.0 +0.11,447,85.0,66.80,14.8,7.0,3,o, 0.23));
      dataArray.push_back(DataInfo(332,"beam_s332_25",1,3, 25.0 +0.15,447,85.0,66.80,14.8,7.0,3,o, 0.32));
      dataArray.push_back(DataInfo(332,"beam_s332_30",1,3, 30.0 +0.17,447,85.0,66.80,14.8,7.0,3,o, 0.29));
      dataArray.push_back(DataInfo(332,"beam_s332_40",1,3, 40.0 +0.00,447,85.0,66.80,14.8,7.0,3,o, 0.32));
      dataArray.push_back(DataInfo(332,"beam_s332_50",1,3, 50.0 -0.06,447,85.0,66.80,14.8,7.0,3,o, 0.23));
      dataArray.push_back(DataInfo(332,"beam_s332_60",1,3, 60.0 +0.00,447,85.0,66.80,14.8,7.0,3,o, 0.17));
      dataArray.push_back(DataInfo(332,"beam_s332_70",1,3, 70.0 +0.11,447,85.0,66.80,14.8,7.0,3,o, 0.14));
      dataArray.push_back(DataInfo(332,"beam_s332_80",1,3, 80.0 +0.11,447,85.0,66.80,14.8,7.0,3,o, 0.14));
      dataArray.push_back(DataInfo(332,"beam_s332_90",1,3, 90.0 +0.03,447,85.0,66.80,14.8,7.0,3,o, 0.00));
      dataArray.push_back(DataInfo(332,"beam_s332_100",1,3,100.0-0.28,447,85.0,66.80,14.8,7.0,3,o,-0.06));
      dataArray.push_back(DataInfo(332,"beam_s332_110",1,3,110.0+0.11,447,85.0,66.80,14.8,7.0,3,o,-0.26));
      dataArray.push_back(DataInfo(332,"beam_s332_120",1,3,120.0+0.11,447,85.0,66.80,14.8,7.0,3,o,-0.29));
      dataArray.push_back(DataInfo(332,"beam_s332_130",1,3,130.0+0.34,447,85.0,66.80,14.8,7.0,3,o,-0.29));
      dataArray.push_back(DataInfo(332,"beam_s332_140",1,3,140.0+0.30,447,85.0,66.80,14.8,7.0,3,o,-0.30));
      dataArray.push_back(DataInfo(332,"beam_s332_150",1,3,150.0+0.11,447,85.0,66.80,14.8,7.0,3,o,-0.31));
    }

    study[333]="X scan; theta=25; phi=0; bar + grease + 3LS lens";
    {
      Double_t o = 28.81;
      // study id | run name | radiator | lens | angle | z pos | x pos | x step | y step | momentum | beam dim | sim offset | phi
      dataArray.push_back(DataInfo(333,"beam_17255152659",1,3,25.0,447,125.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(333,"beam_17255152233",1,3,25.0,447,121.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(333,"beam_17255151854",1,3,25.0,447,117.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(333,"beam_17255151437",1,3,25.0,447,113.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(333,"beam_17255151040",1,3,25.0,447,109.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(333,"beam_17255150625",1,3,25.0,447,105.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(333,"beam_17255150223",1,3,25.0,447,101.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(333,"beam_17255145801",1,3,25.0,447,97.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(333,"beam_17255145344",1,3,25.0,447,93.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(333,"beam_17255144946",1,3,25.0,447,89.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(333,"beam_17255144431",1,3,25.0,447,85.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(333,"beam_17255143352",1,3,25.0,447,81.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(333,"beam_17255141534",1,3,25.0,447,77.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(333,"beam_17255141056",1,3,25.0,447,73.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(333,"beam_17255140703",1,3,25.0,447,69.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(333,"beam_17255134047",1,3,25.0,447,65.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(333,"beam_17255134435",1,3,25.0,447,61.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(333,"beam_17255134834",1,3,25.0,447,57.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(333,"beam_17255135222",1,3,25.0,447,53.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(333,"beam_17255135819",1,3,25.0,447,49,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(333,"beam_17255140213",1,3,25.0,447,45.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(333,"beam_17255153345",1,3,25.0,447,41.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(333,"beam_17255153904",1,3,25.0,447,37.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(333,"beam_17255154411",1,3,25.0,447,33.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(333,"beam_17255154850",1,3,25.0,447,29.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(333,"beam_17255155308",1,3,25.0,447,25.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(333,"beam_17255155835",1,3,25.0,447,-5.0,70.00,0,7.0,3,o,0));
      dataArray.push_back(DataInfo(333,"beam_17255160454",1,3,25.0,447,-15.0,70.00,0,7.0,3,o,0));
    }

    study[334]="Theta scan; phi=0; bar + grease + 3LS lens; 10 GeV/c";
    {
      Double_t o = 28.81;
      // study id | run name | radiator | lens | angle | z pos | x pos | x step | y step | momentum | beam dim | sim offset | phi
      dataArray.push_back(DataInfo(334,"beam_s334_20" ,1,3,20.0, 447,85.0,70.00,0,10.0,3,o,0));
      dataArray.push_back(DataInfo(334,"beam_s334_25" ,1,3,25.0, 447,85.0,70.00,0,10.0,3,o,0));
      dataArray.push_back(DataInfo(334,"beam_s334_30" ,1,3,30.0, 447,85.0,70.00,0,10.0,3,o,0));
      dataArray.push_back(DataInfo(334,"beam_s334_60" ,1,3,60.0, 447,85.0,70.00,0,10.0,3,o,0));
      dataArray.push_back(DataInfo(334,"beam_s334_120",1,3,120.0,447,85.0,70.00,0,10.0,3,o,0));
    }
    
  }
  

  //= Jul 2018 =================================
  {
    Double_t o =  74.20;
    study[400]="Fast angle scan. Bar, 3LS lens, grease everywere";
    {
      // study id | run name | radiator | lens | angle | z pos | x pos | x step | y step | momentum | beam dim | sim offset
      dataArray.push_back(DataInfo(400,"beam_18212214716",1,3,150.0,442,85.0,70.00,16.45,7.0,10,o));
      dataArray.push_back(DataInfo(400,"beam_18212214000",1,3,140.0,442,85.0,70.00,16.45,7.0,10,o));
      dataArray.push_back(DataInfo(400,"beam_18212212157",1,3,130.0,442,85.0,70.00,16.45,7.0,10,o));
      dataArray.push_back(DataInfo(400,"beam_18212211420",1,3,120.0,442,85.0,70.00,16.45,7.0,10,o));
      dataArray.push_back(DataInfo(400,"beam_18212210407",1,3,110.0,442,85.0,70.00,16.45,7.0,10,o));
      dataArray.push_back(DataInfo(400,"beam_18212205205",1,3,100.0,442,85.0,70.00,16.45,7.0,10,o));
      dataArray.push_back(DataInfo(400,"beam_18212203705",1,3, 90.0,442,85.0,70.00,16.45,7.0,10,o));
      dataArray.push_back(DataInfo(400,"beam_18212202744",1,3, 80.0,442,85.0,70.00,16.45,7.0,10,o));
      dataArray.push_back(DataInfo(400,"beam_18212201933",1,3, 70.0,442,85.0,70.00,16.45,7.0,10,o));
      dataArray.push_back(DataInfo(400,"beam_18212200803",1,3, 60.0,442,85.0,70.00,16.45,7.0,10,o));
      dataArray.push_back(DataInfo(400,"beam_18212195303",1,3, 50.0,442,85.0,70.00,16.45,7.0,10,o));
      dataArray.push_back(DataInfo(400,"beam_18212194136",1,3, 40.0,442,85.0,70.00,16.45,7.0,10,o));
      dataArray.push_back(DataInfo(400,"beam_18212193358",1,3, 30.0,442,85.0,70.00,16.45,7.0,10,o));
      dataArray.push_back(DataInfo(400,"beam_18212192707",1,3, 25.0,442,85.0,70.00,16.45,7.0,10,o));
      dataArray.push_back(DataInfo(400,"beam_18212192129",1,3, 20.0,442,85.0,70.00,16.45,7.0,10,o));
    }
      
    study[401]="High stat scan. Bar, 3LS lens, grease everywere";
    {
      // study id | run name | radiator | lens | angle | z pos | x pos | x step | y step | momentum | beam dim | sim offset
      dataArray.push_back(DataInfo(401,"beam_401_20",1,3, 20.0,442,85.0,70.00,16.45,7.0,10,o));
      dataArray.push_back(DataInfo(401,"beam_401_25",1,3, 25.0,442,85.0,70.00,16.45,7.0,10,o));
      dataArray.push_back(DataInfo(401,"beam_401_30",1,3, 30.0,442,85.0,70.00,16.45,7.0,10,o));
      dataArray.push_back(DataInfo(401,"beam_401_60",1,3, 60.0,442,85.0,70.00,16.45,7.0,10,o));
      dataArray.push_back(DataInfo(401,"beam_401_90",1,3, 90.0,442,85.0,70.00,16.45,7.0,10,o));
      dataArray.push_back(DataInfo(401,"beam_401_140",1,3, 140.0,442,85.0,70.00,16.45,7.0,10,o));
    }

    study[402]="Offset study. Bar, 3LS lens, grease everywere, 20 deg, 2 mV";
    {
      // study id | run name | radiator | lens | angle | z pos | x pos | x step | y step | momentum | beam dim | sim offset
      dataArray.push_back(DataInfo(402,"beam_402_20",1,3, 20.0,442,85.0,70.00,16.45,7.0,10,o));
    }
    
    // study[403]="Theta scan @ phi=0 deg. Bar, 3LS lens, grease everywere";
    // {
    //   // study id | run name | radiator | lens | angle | z pos | x pos | x step | y step | momentum | beam dim | sim offset
    //   dataArray.push_back(DataInfo(403,"beam_18215163732",1,3,20.0  -0.286479 ,442,85.0,70.00,16.45,7.0,10,o, 1.14592  ));
    //   dataArray.push_back(DataInfo(403,"beam_18215164946",1,3,25.0  -0.515662 ,442,85.0,70.00,16.45,7.0,10,o, 0.859437 ));
    //   dataArray.push_back(DataInfo(403,"beam_18215170338",1,3,30.0  -0.40107  ,442,85.0,70.00,16.45,7.0,10,o, 0.802141 ));
    //   dataArray.push_back(DataInfo(403,"beam_18215172752",1,3,35.0  -0.229183 ,442,85.0,70.00,16.45,7.0,10,o, 0.744845 ));
    //   dataArray.push_back(DataInfo(403,"beam_18215174559",1,3,40.0  -0.286479 ,442,85.0,70.00,16.45,7.0,10,o, 0.630254 ));
    //   dataArray.push_back(DataInfo(403,"beam_18215180516",1,3,45.0  -0.171887 ,442,85.0,70.00,16.45,7.0,10,o, 0.515662 ));
    //   dataArray.push_back(DataInfo(403,"beam_18215182758",1,3,50.0  -0.229183 ,442,85.0,70.00,16.45,7.0,10,o, 0.630254 ));
    //   dataArray.push_back(DataInfo(403,"beam_18215184548",1,3,55.0  -0.229183 ,442,85.0,70.00,16.45,7.0,10,o, 0.630254 ));
    //   dataArray.push_back(DataInfo(403,"beam_18215190401",1,3,60.0  +0.229183 ,442,85.0,70.00,16.45,7.0,10,o, 0.458366 ));
    //   dataArray.push_back(DataInfo(403,"beam_18215192247",1,3,65.0  +0.40107  ,442,85.0,70.00,16.45,7.0,10,o, 0.343775 ));
    //   dataArray.push_back(DataInfo(403,"beam_18215194145",1,3,70.0  +0.229183 ,442,85.0,70.00,16.45,7.0,10,o, 0.458366 ));
    //   dataArray.push_back(DataInfo(403,"beam_18215200044",1,3,75.0  +0.114592 ,442,85.0,70.00,16.45,7.0,10,o, 0.343775 ));
    //   dataArray.push_back(DataInfo(403,"beam_18215201936",1,3,80.0  +0.286479 ,442,85.0,70.00,16.45,7.0,10,o, 0.458366 ));
    //   dataArray.push_back(DataInfo(403,"beam_18215203823",1,3,85.0  +0.343775 ,442,85.0,70.00,16.45,7.0,10,o, 0.286479 ));
    //   dataArray.push_back(DataInfo(403,"beam_18215210112",1,3,90.0  -0.286479 ,442,85.0,70.00,16.45,7.0,10,o, 0.458366 ));
    //   dataArray.push_back(DataInfo(403,"beam_18215212118",1,3,95.0  -0.286479 ,442,85.0,70.00,16.45,7.0,10,o, 0.40107  ));
    //   dataArray.push_back(DataInfo(403,"beam_18215213958",1,3,100.0 -0.687549 ,442,85.0,70.00,16.45,7.0,10,o, 0.343775 ));
    //   dataArray.push_back(DataInfo(403,"beam_18215220217",1,3,105.0 -0.916732 ,442,85.0,70.00,16.45,7.0,10,o, 0.171887 ));
    //   dataArray.push_back(DataInfo(403,"beam_18215222213",1,3,110.0 -0.40107  ,442,85.0,70.00,16.45,7.0,10,o, 0.286479 ));
    //   dataArray.push_back(DataInfo(403,"beam_18215224129",1,3,115.0 -0.572958 ,442,85.0,70.00,16.45,7.0,10,o, 0.458366 ));
    //   dataArray.push_back(DataInfo(403,"beam_18215230158",1,3,120.0 -0.572958 ,442,85.0,70.00,16.45,7.0,10,o, 0.343775 ));
    //   dataArray.push_back(DataInfo(403,"beam_18215231920",1,3,125.0 -0.286479 ,442,85.0,70.00,16.45,7.0,10,o, 0.572958 ));
    //   dataArray.push_back(DataInfo(403,"beam_18215233358",1,3,130.0 -0.000000 ,442,85.0,70.00,16.45,7.0,10,o, 0.572958 ));
    //   dataArray.push_back(DataInfo(403,"beam_18215234920",1,3,135.0 -0.229183 ,442,85.0,70.00,16.45,7.0,10,o, 0.40107  ));
    //   dataArray.push_back(DataInfo(403,"beam_18216000538",1,3,140.0 -0.458366 ,442,85.0,70.00,16.45,7.0,10,o, 0.458366 ));
    // }

    study[403]="Theta scan @ phi=0 deg. Bar, 3LS lens, grease everywere";
    {
      // study id | run name | radiator | lens | angle | z pos | x pos | x step | y step | momentum | beam dim | sim offset
      dataArray.push_back(DataInfo(403,"beam_18215163732",1,3,20.0   ,442,85.0,70.00,16.45,7.0,10,o ));
      dataArray.push_back(DataInfo(403,"beam_18215164946",1,3,25.0   ,442,85.0,70.00,16.45,7.0,10,o ));
      dataArray.push_back(DataInfo(403,"beam_18215170338",1,3,30.0   ,442,85.0,70.00,16.45,7.0,10,o ));
      dataArray.push_back(DataInfo(403,"beam_18215172752",1,3,35.0   ,442,85.0,70.00,16.45,7.0,10,o ));
      dataArray.push_back(DataInfo(403,"beam_18215174559",1,3,40.0   ,442,85.0,70.00,16.45,7.0,10,o ));
      dataArray.push_back(DataInfo(403,"beam_18215180516",1,3,45.0   ,442,85.0,70.00,16.45,7.0,10,o ));
      dataArray.push_back(DataInfo(403,"beam_18215182758",1,3,50.0   ,442,85.0,70.00,16.45,7.0,10,o ));
      dataArray.push_back(DataInfo(403,"beam_18215184548",1,3,55.0   ,442,85.0,70.00,16.45,7.0,10,o ));
      dataArray.push_back(DataInfo(403,"beam_18215190401",1,3,60.0   ,442,85.0,70.00,16.45,7.0,10,o ));
      dataArray.push_back(DataInfo(403,"beam_18215192247",1,3,65.0   ,442,85.0,70.00,16.45,7.0,10,o ));
      dataArray.push_back(DataInfo(403,"beam_18215194145",1,3,70.0   ,442,85.0,70.00,16.45,7.0,10,o ));
      dataArray.push_back(DataInfo(403,"beam_18215200044",1,3,75.0   ,442,85.0,70.00,16.45,7.0,10,o ));
      dataArray.push_back(DataInfo(403,"beam_18215201936",1,3,80.0   ,442,85.0,70.00,16.45,7.0,10,o ));
      dataArray.push_back(DataInfo(403,"beam_18215203823",1,3,85.0   ,442,85.0,70.00,16.45,7.0,10,o ));
      dataArray.push_back(DataInfo(403,"beam_18215210112",1,3,90.0   ,442,85.0,70.00,16.45,7.0,10,o ));
      dataArray.push_back(DataInfo(403,"beam_18215212118",1,3,95.0   ,442,85.0,70.00,16.45,7.0,10,o ));
      dataArray.push_back(DataInfo(403,"beam_18215213958",1,3,100.0  ,442,85.0,70.00,16.45,7.0,10,o ));
      dataArray.push_back(DataInfo(403,"beam_18215220217",1,3,105.0  ,442,85.0,70.00,16.45,7.0,10,o ));
      dataArray.push_back(DataInfo(403,"beam_18215222213",1,3,110.0  ,442,85.0,70.00,16.45,7.0,10,o ));
      dataArray.push_back(DataInfo(403,"beam_18215224129",1,3,115.0  ,442,85.0,70.00,16.45,7.0,10,o ));
      dataArray.push_back(DataInfo(403,"beam_18215230158",1,3,120.0  ,442,85.0,70.00,16.45,7.0,10,o ));
      dataArray.push_back(DataInfo(403,"beam_18215231920",1,3,125.0  ,442,85.0,70.00,16.45,7.0,10,o ));
      dataArray.push_back(DataInfo(403,"beam_18215233358",1,3,130.0  ,442,85.0,70.00,16.45,7.0,10,o ));
      dataArray.push_back(DataInfo(403,"beam_18215234920",1,3,135.0  ,442,85.0,70.00,16.45,7.0,10,o ));
      dataArray.push_back(DataInfo(403,"beam_18216000538",1,3,140.0  ,442,85.0,70.00,16.45,7.0,10,o ));
    }

    study[404]="High stat run @ theta=20 @ phi=0 deg. Bar, 3LS lens, grease everywere";
    {
      // study id | run name | radiator | lens | angle | z pos | x pos | x step | y step | momentum | beam dim | sim offset
      dataArray.push_back(DataInfo(404,"beam_404_20",1,3,20.0,442,85.0,70.00,16.45,7.0,10,o));
      dataArray.push_back(DataInfo(404,"beam_404_25",1,3,25.0,442,85.0,70.00,16.45,7.0,10,o));
    }
    
    study[405]="Threshold scan @ 20 deg @ 7 GeV/c. Bar, 3LS lens, grease everywere";
    {
      // study id | run name | radiator | lens | angle | z pos | x pos | x step | y step | momentum | beam dim | sim offset | phi | test
      dataArray.push_back(DataInfo(405,"beam_18216095817",1,3,20.0,442,85.0,70.00,16.45,7.0,10,o,0,3000));
      dataArray.push_back(DataInfo(405,"beam_18216101618",1,3,20.0,442,85.0,70.00,16.45,7.0,10,o,0,2700));
      dataArray.push_back(DataInfo(405,"beam_18216103400",1,3,20.0,442,85.0,70.00,16.45,7.0,10,o,0,2400));
      dataArray.push_back(DataInfo(405,"beam_18216105109",1,3,20.0,442,85.0,70.00,16.45,7.0,10,o,0,2100));
      dataArray.push_back(DataInfo(405,"beam_18216110707",1,3,20.0,442,85.0,70.00,16.45,7.0,10,o,0,1800));
      dataArray.push_back(DataInfo(405,"beam_18216112220",1,3,20.0,442,85.0,70.00,16.45,7.0,10,o,0,1500));
      dataArray.push_back(DataInfo(405,"beam_18216113629",1,3,20.0,442,85.0,70.00,16.45,7.0,10,o,0,1200));
      dataArray.push_back(DataInfo(405,"beam_18216115110",1,3,20.0,442,85.0,70.00,16.45,7.0,10,o,0,900));   
      dataArray.push_back(DataInfo(405,"beam_18216120423",1,3,20.0,442,85.0,70.00,16.45,7.0,10,o,0,600));
      dataArray.push_back(DataInfo(405,"beam_18216121904",1,3,20.0,442,85.0,70.00,16.45,7.0,10,o,0,300));   
    }
    
    study[406]="Theta scan at phi=5 deg. Bar, 3LS lens, grease everywere";
    {
      // study id | run name | radiator | lens | angle | z pos | x pos | x step | y step | momentum | beam dim | sim offset | phi
      dataArray.push_back(DataInfo(406,"beam_18216135640",1,3,20.0,442,85.0,70.00,16.45,7.0,10,o,5));
      dataArray.push_back(DataInfo(406,"beam_18216141102",1,3,30.0,442,85.0,70.00,16.45,7.0,10,o,5));
      dataArray.push_back(DataInfo(406,"beam_18216142557",1,3,40.0,442,85.0,70.00,16.45,7.0,10,o,5));
      dataArray.push_back(DataInfo(406,"beam_18216144258",1,3,50.0,442,85.0,70.00,16.45,7.0,10,o,5));
      dataArray.push_back(DataInfo(406,"beam_18216151341",1,3,60.0,442,85.0,70.00,16.45,7.0,10,o,5));
      dataArray.push_back(DataInfo(406,"beam_18216153216",1,3,70.0,442,85.0,70.00,16.45,7.0,10,o,5));
      dataArray.push_back(DataInfo(406,"beam_18216155109",1,3,80.0,442,85.0,70.00,16.45,7.0,10,o,5));
      dataArray.push_back(DataInfo(406,"beam_18216161002",1,3,90.0,442,85.0,70.00,16.45,7.0,10,o,5));
      dataArray.push_back(DataInfo(406,"beam_18216162947",1,3,100.0,442,85.0,70.00,16.45,7.0,10,o,5));
      dataArray.push_back(DataInfo(406,"beam_18216165048",1,3,110.0,442,85.0,70.00,16.45,7.0,10,o,5));
      dataArray.push_back(DataInfo(406,"beam_18216171151",1,3,120.0,442,85.0,70.00,16.45,7.0,10,o,5));
      dataArray.push_back(DataInfo(406,"beam_18216173240",1,3,130.0,442,85.0,70.00,16.45,7.0,10,o,5));
      dataArray.push_back(DataInfo(406,"beam_18216175149",1,3,140.0,442,85.0,70.00,16.45,7.0,10,o,5));
    }

    study[407]="Theta scan at phi=10 deg. Bar, 3LS lens, grease everywere";
    {
      // study id | run name | radiator | lens | angle | z pos | x pos | x step | y step | momentum | beam dim | sim offset | phi
      dataArray.push_back(DataInfo(407,"beam_18217025337",1,3,20.0,442,85.0,70.00,16.45,7.0,10,o,10));
      dataArray.push_back(DataInfo(407,"beam_18217022618",1,3,30.0,442,85.0,70.00,16.45,7.0,10,o,10));
      dataArray.push_back(DataInfo(407,"beam_18217020804",1,3,40.0,442,85.0,70.00,16.45,7.0,10,o,10));
      dataArray.push_back(DataInfo(407,"beam_18217014909",1,3,50.0,442,85.0,70.00,16.45,7.0,10,o,10));
      dataArray.push_back(DataInfo(407,"beam_18217012945",1,3,60.0,442,85.0,70.00,16.45,7.0,10,o,10));
      dataArray.push_back(DataInfo(407,"beam_18217011020",1,3,70.0,442,85.0,70.00,16.45,7.0,10,o,10));
      dataArray.push_back(DataInfo(407,"beam_18217004709",1,3,80.0,442,85.0,70.00,16.45,7.0,10,o,10));
      dataArray.push_back(DataInfo(407,"beam_18217001626",1,3,90.0,442,85.0,70.00,16.45,7.0,10,o,10));
      dataArray.push_back(DataInfo(407,"beam_18216235611",1,3,100.0,442,85.0,70.00,16.45,7.0,10,o,10));
      dataArray.push_back(DataInfo(407,"beam_18216233650",1,3,110.0,442,85.0,70.00,16.45,7.0,10,o,10));
      dataArray.push_back(DataInfo(407,"beam_18216232116",1,3,120.0,442,85.0,70.00,16.45,7.0,10,o,10));
      dataArray.push_back(DataInfo(407,"beam_18216230450",1,3,130.0,442,85.0,70.00,16.45,7.0,10,o,10));
      dataArray.push_back(DataInfo(407,"beam_18216224323",1,3,140.0,442,85.0,70.00,16.45,7.0,10,o,10));
    }

    study[408]="Theta scan at phi=15 deg. Bar, 3LS lens, grease everywere";
    {
      Double_t o =  74.20;
      // study id | run name | radiator | lens | angle | z pos | x pos | x step | y step | momentum | beam dim | sim offset | phi
      dataArray.push_back(DataInfo(408,"beam_18217113055",1,3,20.0,442,85.0,70.00,16.45,7.0,10,o,15));
      dataArray.push_back(DataInfo(408,"beam_18217114922",1,3,30.0,442,85.0,70.00,16.45,7.0,10,o,15));
      dataArray.push_back(DataInfo(408,"beam_18217123110",1,3,40.0,442,85.0,70.00,16.45,7.0,10,o,15));
      dataArray.push_back(DataInfo(408,"beam_18217125724",1,3,50.0,442,85.0,70.00,16.45,7.0,10,o,15));
      dataArray.push_back(DataInfo(408,"beam_18217132900",1,3,60.0,442,85.0,70.00,16.45,7.0,10,o,15));
      dataArray.push_back(DataInfo(408,"beam_18217135557",1,3,70.0,442,85.0,70.00,16.45,7.0,10,o,15));
      dataArray.push_back(DataInfo(408,"beam_18217142313",1,3,80.0,442,85.0,70.00,16.45,7.0,10,o,15));
      dataArray.push_back(DataInfo(408,"beam_18217145344",1,3,90.0,442,85.0,70.00,16.45,7.0,10,o,15));
      dataArray.push_back(DataInfo(408,"beam_18217152331",1,3,100.0,442,85.0,70.00,16.45,7.0,10,o,15));
      dataArray.push_back(DataInfo(408,"beam_18217165144",1,3,110.0,442,85.0,70.00,16.45,7.0,10,o,15));
      dataArray.push_back(DataInfo(408,"beam_18217172018",1,3,120.0,442,85.0,70.00,16.45,7.0,10,o,15));
      dataArray.push_back(DataInfo(408,"beam_18217175308",1,3,130.0,442,85.0,70.00,16.45,7.0,10,o,15));
      dataArray.push_back(DataInfo(408,"beam_18217182228",1,3,140.0,442,85.0,70.00,16.45,7.0,10,o,15));
    }
    
    study[409]= "Theta scan at phi=2.5 deg. Bar, 3LS lens, grease everywere";
    {
      // study id | run name | radiator | lens | angle | z pos | x pos | x step | y step | momentum | beam dim | sim offset | phi
      dataArray.push_back(DataInfo(409,"beam_18218021228",1,3,20.0,442,85.0,70.00,16.45,7.0,10,o,2.5));
      dataArray.push_back(DataInfo(409,"beam_18218014740",1,3,30.0,442,85.0,70.00,16.45,7.0,10,o,2.5));
      dataArray.push_back(DataInfo(409,"beam_18218012045",1,3,40.0,442,85.0,70.00,16.45,7.0,10,o,2.5));
      dataArray.push_back(DataInfo(409,"beam_18218004711",1,3,50.0,442,85.0,70.00,16.45,7.0,10,o,2.5));
      dataArray.push_back(DataInfo(409,"beam_18218001701",1,3,60.0,442,85.0,70.00,16.45,7.0,10,o,2.5));
      dataArray.push_back(DataInfo(409,"beam_18217234746",1,3,70.0,442,85.0,70.00,16.45,7.0,10,o,2.5));
      dataArray.push_back(DataInfo(409,"beam_18217231747",1,3,80.0,442,85.0,70.00,16.45,7.0,10,o,2.5));
      dataArray.push_back(DataInfo(409,"beam_18217224253",1,3,90.0,442,85.0,70.00,16.45,7.0,10,o,2.5));
      dataArray.push_back(DataInfo(409,"beam_18217221819",1,3,100.0,442,85.0,70.00,16.45,7.0,10,o,2.5));
      dataArray.push_back(DataInfo(409,"beam_18217215134",1,3,110.0,442,85.0,70.00,16.45,7.0,10,o,2.5));
      dataArray.push_back(DataInfo(409,"beam_18217212916",1,3,120.0,442,85.0,70.00,16.45,7.0,10,o,2.5));
      dataArray.push_back(DataInfo(409,"beam_18217210456",1,3,130.0,442,85.0,70.00,16.45,7.0,10,o,2.5));
      dataArray.push_back(DataInfo(409,"beam_18217204226",1,3,140.0,442,85.0,70.00,16.45,7.0,10,o,2.5));
    }

    study[415]= "Momentum scan at theta=25 deg and phi=0 deg. Bar, 3LS lens, grease everywere";
    {
      // study id | run name | radiator | lens | angle | z pos | x pos | x step | y step | momentum | beam dim | sim offset | phi
      dataArray.push_back(DataInfo(415,"beam_415_10",1,3,25.0,424,85.0,70.00,16.45,10.0,10,o,0));
      dataArray.push_back(DataInfo(415,"beam_415_9", 1,3,25.0,424,85.0,70.00,16.45, 9.0,10,o,0));
      dataArray.push_back(DataInfo(415,"beam_415_8", 1,3,25.0,424,85.0,70.00,16.45, 8.0,10,o,0));
      dataArray.push_back(DataInfo(415,"beam_415_7", 1,3,25.0,424,85.0,70.00,16.45, 7.0,10,o,0));
      dataArray.push_back(DataInfo(415,"beam_415_6", 1,3,25.0,424,85.0,70.00,16.45, 6.0,10,o,0));
      dataArray.push_back(DataInfo(415,"beam_415_5", 1,3,25.0,424,85.0,70.00,16.45, 5.0,10,o,0));
      dataArray.push_back(DataInfo(415,"beam_415_4", 1,3,25.0,424,85.0,70.00,16.45, 4.0,10,o,0));
      dataArray.push_back(DataInfo(415,"beam_415_3", 1,3,25.0,424,85.0,70.00,16.45, 3.0,10,o,0));
      dataArray.push_back(DataInfo(415,"beam_415_2", 1,3,25.0,424,85.0,70.00,16.45, 2.0,10,o,0));
    }

    study[420]= "Z scan at theta=90 deg and phi=0 deg @ 7 GeV/c. Bar, 3LS lens, grease everywere";
    {
      // study id | run name | radiator | lens | angle | z pos | x pos | x step | y step | momentum | beam dim | sim offset | phi
      dataArray.push_back(DataInfo(420,"beam_18219154515",1,3,90.0,250,85.0,70.00,16.45,7.0,10,o,0));
      dataArray.push_back(DataInfo(420,"beam_18219142232",1,3,90.0,424,85.0,70.00,16.45,7.0,10,o,0));
      dataArray.push_back(DataInfo(420,"beam_18219161250",1,3,90.0,550,85.0,70.00,16.45,7.0,10,o,0));
      dataArray.push_back(DataInfo(420,"beam_18219163624",1,3,90.0,650,85.0,70.00,16.45,7.0,10,o,0));
      dataArray.push_back(DataInfo(420,"beam_18219165621",1,3,90.0,750,85.0,70.00,16.45,7.0,10,o,0));
      dataArray.push_back(DataInfo(420,"beam_18219171427",1,3,90.0,850,85.0,70.00,16.45,7.0,10,o,0));
    }
    study[421]= "With w/o EDD at theta=25 deg and phi=0 deg @ 7 GeV/c. Bar, 3LS lens, grease everywere";
    {
      // study id | run name | radiator | lens | angle | z pos | x pos | x step | y step | momentum | beam dim | sim offset | phi
      dataArray.push_back(DataInfo(421,"beam_421_wo",1,3,25.0,424,85.0,70.00,16.45,7.0,10,o,0));
      dataArray.push_back(DataInfo(421,"beam_421_wi",1,3,25.0,424,85.0,70.00,16.45,7.0,10,o,0));
    }
    study[422]= "Theta scan at phi=0 deg @ 4 GeV/c. Bar, 3LS lens, grease everywere";
    {
      // study id | run name | radiator | lens | angle | z pos | x pos | x step | y step | momentum | beam dim | sim offset | phi
      dataArray.push_back(DataInfo(422,"beam_422_20",1,3,20.0,424,85.0,70.00,16.45,4.0,10,o,0));
      dataArray.push_back(DataInfo(422,"beam_422_25",1,3,25.0,424,85.0,70.00,16.45,4.0,10,o,0));
      dataArray.push_back(DataInfo(422,"beam_422_30",1,3,30.0,424,85.0,70.00,16.45,4.0,10,o,0));
      dataArray.push_back(DataInfo(422,"beam_422_40",1,3,40.0,424,85.0,70.00,16.45,4.0,10,o,0));
      dataArray.push_back(DataInfo(422,"beam_422_45",1,3,45.0,424,85.0,70.00,16.45,4.0,10,o,0));
      dataArray.push_back(DataInfo(422,"beam_422_50",1,3,50.0,424,85.0,70.00,16.45,4.0,10,o,0));
      dataArray.push_back(DataInfo(422,"beam_422_55",1,3,55.0,424,85.0,70.00,16.45,4.0,10,o,0));
      dataArray.push_back(DataInfo(422,"beam_422_60",1,3,60.0,424,85.0,70.00,16.45,4.0,10,o,0));
      dataArray.push_back(DataInfo(422,"beam_422_70",1,3,70.0,424,85.0,70.00,16.45,4.0,10,o,0));
      dataArray.push_back(DataInfo(422,"beam_422_80",1,3,80.0,424,85.0,70.00,16.45,4.0,10,o,0));
      dataArray.push_back(DataInfo(422,"beam_422_90",1,3,90.0,424,85.0,70.00,16.45,4.0,10,o,0));
      dataArray.push_back(DataInfo(422,"beam_422_100",1,3,100.0,424,85.0,70.00,16.45,4.0,10,o,0));
      dataArray.push_back(DataInfo(422,"beam_422_110",1,3,110.0,424,85.0,70.00,16.45,4.0,10,o,0));
      dataArray.push_back(DataInfo(422,"beam_422_120",1,3,120.0,424,85.0,70.00,16.45,4.0,10,o,0));
      dataArray.push_back(DataInfo(422,"beam_422_130",1,3,130.0,424,85.0,70.00,16.45,4.0,10,o,0));
      dataArray.push_back(DataInfo(422,"beam_422_140",1,3,140.0,424,85.0,70.00,16.45,4.0,10,o,0));			
    }

    study[423]= "Fine theta scan at phi=0 deg @ 1.5 GeV/c. Bar, 3LS lens, grease everywere";
    {
      dataArray.push_back(DataInfo(423,"beam_18220132601",1,3,90.0,424,85.0,70.00,16.45,1.5,10,o,0));
      dataArray.push_back(DataInfo(423,"beam_18220132829",1,3,91.0,424,85.0,70.00,16.45,1.5,10,o,0));
      dataArray.push_back(DataInfo(423,"beam_18220133026",1,3,92.0,424,85.0,70.00,16.45,1.5,10,o,0));
      dataArray.push_back(DataInfo(423,"beam_18220133224",1,3,93.0,424,85.0,70.00,16.45,1.5,10,o,0));
      dataArray.push_back(DataInfo(423,"beam_18220133424",1,3,94.0,424,85.0,70.00,16.45,1.5,10,o,0));
      dataArray.push_back(DataInfo(423,"beam_18220133641",1,3,95.0,424,85.0,70.00,16.45,1.5,10,o,0));
      dataArray.push_back(DataInfo(423,"beam_18220133858",1,3,96.0,424,85.0,70.00,16.45,1.5,10,o,0));
      dataArray.push_back(DataInfo(423,"beam_18220134101",1,3,97.0,424,85.0,70.00,16.45,1.5,10,o,0));
      dataArray.push_back(DataInfo(423,"beam_18220134302",1,3,98.0,424,85.0,70.00,16.45,1.5,10,o,0));
      dataArray.push_back(DataInfo(423,"beam_18220134524",1,3,99.0,424,85.0,70.00,16.45,1.5,10,o,0));
      dataArray.push_back(DataInfo(423,"beam_18220134739",1,3,100.0,424,85.0,70.00,16.45,1.5,10,o,0));
      dataArray.push_back(DataInfo(423,"beam_18220134940",1,3,101.0,424,85.0,70.00,16.45,1.5,10,o,0));
      dataArray.push_back(DataInfo(423,"beam_18220135147",1,3,102.0,424,85.0,70.00,16.45,1.5,10,o,0));
      dataArray.push_back(DataInfo(423,"beam_18220135359",1,3,103.0,424,85.0,70.00,16.45,1.5,10,o,0));
      dataArray.push_back(DataInfo(423,"beam_18220135600",1,3,104.0,424,85.0,70.00,16.45,1.5,10,o,0));
      dataArray.push_back(DataInfo(423,"beam_18220135800",1,3,105.0,424,85.0,70.00,16.45,1.5,10,o,0));
      dataArray.push_back(DataInfo(423,"beam_18220140009",1,3,110.0,424,85.0,70.00,16.45,1.5,10,o,0));
      dataArray.push_back(DataInfo(423,"beam_18220140437",1,3,115.0,424,85.0,70.00,16.45,1.5,10,o,0));
      dataArray.push_back(DataInfo(423,"beam_18220140652",1,3,120.0,424,85.0,70.00,16.45,1.5,10,o,0));
      dataArray.push_back(DataInfo(423,"beam_18220140859",1,3,130.0,424,85.0,70.00,16.45,1.5,10,o,0));
      dataArray.push_back(DataInfo(423,"beam_18220141114",1,3,135.0,424,85.0,70.00,16.45,1.5,10,o,0));
      dataArray.push_back(DataInfo(423,"beam_18220141430",1,3,140.0,424,85.0,70.00,16.45,1.5,10,o,0));
      dataArray.push_back(DataInfo(423,"beam_18220141845",1,3,80.0,424,85.0,70.00,16.45,1.5,10,o,0));
      dataArray.push_back(DataInfo(423,"beam_18220142854",1,3,85.0,424,85.0,70.00,16.45,1.5,10,o,0));
    }
    
    study[424]= "Fine theta scan at phi=0 deg @ 2.0 GeV/c. Bar, 3LS lens, grease everywere";
    {
      dataArray.push_back(DataInfo(424,"beam_18220144642",1,3,20.0,424,85.0,70.00,16.45,2.0,10,o,0));
      dataArray.push_back(DataInfo(424,"beam_18220145129",1,3,30.0,424,85.0,70.00,16.45,2.0,10,o,0));
      dataArray.push_back(DataInfo(424,"beam_18220145842",1,3,40.0,424,85.0,70.00,16.45,2.0,10,o,0));
      dataArray.push_back(DataInfo(424,"beam_18220150324",1,3,50.0,424,85.0,70.00,16.45,2.0,10,o,0));
      dataArray.push_back(DataInfo(424,"beam_18220150922",1,3,60.0,424,85.0,70.00,16.45,2.0,10,o,0));
      dataArray.push_back(DataInfo(424,"beam_18220151434",1,3,70.0,424,85.0,70.00,16.45,2.0,10,o,0));
      dataArray.push_back(DataInfo(424,"beam_18220152011",1,3,75.0,424,85.0,70.00,16.45,2.0,10,o,0));
      dataArray.push_back(DataInfo(424,"beam_18220152458",1,3,80.0,424,85.0,70.00,16.45,2.0,10,o,0));
      dataArray.push_back(DataInfo(424,"beam_18220153017",1,3,81.0,424,85.0,70.00,16.45,2.0,10,o,0));
      dataArray.push_back(DataInfo(424,"beam_18220153503",1,3,82.0,424,85.0,70.00,16.45,2.0,10,o,0));
      dataArray.push_back(DataInfo(424,"beam_18220153941",1,3,83.0,424,85.0,70.00,16.45,2.0,10,o,0));
      dataArray.push_back(DataInfo(424,"beam_18220154406",1,3,84.0,424,85.0,70.00,16.45,2.0,10,o,0));
      dataArray.push_back(DataInfo(424,"beam_18220155023",1,3,85.0,424,85.0,70.00,16.45,2.0,10,o,0));
      dataArray.push_back(DataInfo(424,"beam_18220155459",1,3,86.0,424,85.0,70.00,16.45,2.0,10,o,0));
      dataArray.push_back(DataInfo(424,"beam_18220155943",1,3,87.0,424,85.0,70.00,16.45,2.0,10,o,0));
      dataArray.push_back(DataInfo(424,"beam_18220160420",1,3,88.0,424,85.0,70.00,16.45,2.0,10,o,0));
      dataArray.push_back(DataInfo(424,"beam_18220161102",1,3,89.0,424,85.0,70.00,16.45,2.0,10,o,0));
      dataArray.push_back(DataInfo(424,"beam_18220161752",1,3,90.0,424,85.0,70.00,16.45,2.0,10,o,0));
      dataArray.push_back(DataInfo(424,"beam_18220162300",1,3,92.0,424,85.0,70.00,16.45,2.0,10,o,0));
      dataArray.push_back(DataInfo(424,"beam_18220163058",1,3,95.0,424,85.0,70.00,16.45,2.0,10,o,0));
      dataArray.push_back(DataInfo(424,"beam_18220163736",1,3,100.0,424,85.0,70.00,16.45,2.0,10,o,0));
    }

    study[425]= "Theta scan at phi=0 deg @ 8 GeV/c. Bar, 3LS lens, grease everywere";
    {
      dataArray.push_back(DataInfo(425,"beam_425_20",1,3,20.0,424,85.0,70.00,16.45,8.0,10,o,0));
      dataArray.push_back(DataInfo(425,"beam_425_25",1,3,25.0,424,85.0,70.00,16.45,8.0,10,o,0));
      dataArray.push_back(DataInfo(425,"beam_425_30",1,3,30.0,424,85.0,70.00,16.45,8.0,10,o,0));
      dataArray.push_back(DataInfo(425,"beam_425_35",1,3,35.0,424,85.0,70.00,16.45,8.0,10,o,0));
      dataArray.push_back(DataInfo(425,"beam_425_40",1,3,40.0,424,85.0,70.00,16.45,8.0,10,o,0));
      dataArray.push_back(DataInfo(425,"beam_425_50",1,3,50.0,424,85.0,70.00,16.45,8.0,10,o,0));
      dataArray.push_back(DataInfo(425,"beam_425_60",1,3,60.0,424,85.0,70.00,16.45,8.0,10,o,0));
      dataArray.push_back(DataInfo(425,"beam_425_70",1,3,70.0,424,85.0,70.00,16.45,8.0,10,o,0));
      dataArray.push_back(DataInfo(425,"beam_425_80",1,3,80.0,424,85.0,70.00,16.45,8.0,10,o,0));
      dataArray.push_back(DataInfo(425,"beam_425_90",1,3,90.0,424,85.0,70.00,16.45,8.0,10,o,0));
      dataArray.push_back(DataInfo(425,"beam_425_100",1,3,100.0,424,85.0,70.00,16.45,8.0,10,o,0));
      dataArray.push_back(DataInfo(425,"beam_425_110",1,3,110.0,424,85.0,70.00,16.45,8.0,10,o,0));
      dataArray.push_back(DataInfo(425,"beam_425_120",1,3,120.0,424,85.0,70.00,16.45,8.0,10,o,0));
      dataArray.push_back(DataInfo(425,"beam_425_130",1,3,130.0,424,85.0,70.00,16.45,8.0,10,o,0));
      dataArray.push_back(DataInfo(425,"beam_425_140",1,3,140.0,424,85.0,70.00,16.45,8.0,10,o,0));
    }

    study[430]= "Theta scan at phi=0 deg @ 7 GeV/c. Bar, 3LS lens, cookies between bar lens and prism";
    {
      dataArray.push_back(DataInfo(430,"beam_430_20",1,3,20.0,424,85.0,70.00,16.45,7.0,10,o,0));
      dataArray.push_back(DataInfo(430,"beam_430_25",1,3,25.0,424,85.0,70.00,16.45,7.0,10,o,0));
      dataArray.push_back(DataInfo(430,"beam_430_80",1,3,80.0,424,85.0,70.00,16.45,7.0,10,o,0));     
    }
    int id;
    study[431]= "Theta scan at phi=0 deg @ 7 GeV/c. Bar, 3LS lens, all grease again";
    {
      dataArray.push_back(DataInfo(431,"beam_431_20",1,3,20.0,424,85.0,70.00,16.45,7.0,10,o,0));
      dataArray.push_back(DataInfo(431,"beam_18223103847",1,3,20.0,424,85.0,70.00,16.45,7.0,10,o,0));
    }  

    study[434]= "Test; X step; theta scan at phi=0 deg @ 7 GeV/c. Bar, 3LS lens, all grease again";
    {
      // study id | run name | radiator | lens | angle | z pos | x pos | x step | y step | momentum | beam dim | sim offset | phi 
      dataArray.push_back(DataInfo(434,"beam_18223212206",1,3,20.0,424,145.0,125.00,16.45,7.0,10,o,0));
      dataArray.push_back(DataInfo(434,"beam_18223214220",1,3,20.0,424,145.0,125.00,16.45,7.0,10,o,0));
      dataArray.push_back(DataInfo(434,"beam_18223222546",1,3,20.0,424,145.0,125.00,16.45,7.0,10,o,0));
      dataArray.push_back(DataInfo(434,"beam_18223221520",1,3,20.0,424,145.0,125.00,16.45,7.0,10,o,0));
    }

    study[435]= "X step; Theta scan; phi=0 deg @ 7 GeV/c. Bar, 3LS lens, all grease";
    {
      // study id | run name | radiator | lens | angle | z pos | x pos | x step | y step | momentum | beam dim | sim offset | phi 
      dataArray.push_back(DataInfo(435,"beam_435_20",1,3,20.0,424,145.0,125.00,16.45,7.0,10,o,0));
      dataArray.push_back(DataInfo(435,"beam_435_25",1,3,25.0,424,145.0,125.00,16.45,7.0,10,o,0));
      dataArray.push_back(DataInfo(435,"beam_435_30",1,3,30.0,424,145.0,125.00,16.45,7.0,10,o,0));
      dataArray.push_back(DataInfo(435,"beam_435_40",1,3,40.0,424,145.0,125.00,16.45,7.0,10,o,0));
      dataArray.push_back(DataInfo(435,"beam_435_50",1,3,50.0,424,145.0,125.00,16.45,7.0,10,o,0));
      dataArray.push_back(DataInfo(435,"beam_435_60",1,3,60.0,424,145.0,125.00,16.45,7.0,10,o,0));
      dataArray.push_back(DataInfo(435,"beam_435_70",1,3,70.0,424,145.0,125.00,16.45,7.0,10,o,0));
      dataArray.push_back(DataInfo(435,"beam_435_80",1,3,80.0,424,145.0,125.00,16.45,7.0,10,o,0));
      dataArray.push_back(DataInfo(435,"beam_435_90",1,3,90.0,424,145.0,125.00,16.45,7.0,10,o,0));
      dataArray.push_back(DataInfo(435,"beam_435_100",1,3,100.0,424,145.0,125.00,16.45,7.0,10,o,0));
      dataArray.push_back(DataInfo(435,"beam_435_110",1,3,110.0,424,145.0,125.00,16.45,7.0,10,o,0));
      dataArray.push_back(DataInfo(435,"beam_435_120",1,3,120.0,424,145.0,125.00,16.45,7.0,10,o,0));
      dataArray.push_back(DataInfo(435,"beam_435_130",1,3,130.0,424,145.0,125.00,16.45,7.0,10,o,0));
      dataArray.push_back(DataInfo(435,"beam_435_140",1,3,140.0,424,145.0,125.00,16.45,7.0,10,o,0));
    }
    
    study[436]= "X step; Momentum scan; theta=25 deg, phi=0 deg. Bar, 3LS lens, all grease";
    {
      // study id | run name | radiator | lens | angle | z pos | x pos | x step | y step | momentum | beam dim | sim offset | phi 
      dataArray.push_back(DataInfo(436,"beam_436_9",1,3,25.0,424,145.0,125.00,16.45,9.0,10,o,0));
      dataArray.push_back(DataInfo(436,"beam_436_8",1,3,25.0,424,145.0,125.00,16.45,8.0,10,o,0));
      dataArray.push_back(DataInfo(436,"beam_436_7",1,3,25.0,424,145.0,125.00,16.45,7.0,10,o,0));
      dataArray.push_back(DataInfo(436,"beam_436_7_0",1,3,25.0,424,145.0,125.00,16.45,7.01,10,o,0));
      dataArray.push_back(DataInfo(436,"beam_436_6",1,3,25.0,424,145.0,125.00,16.45,6.0,10,o,0));
      dataArray.push_back(DataInfo(436,"beam_436_5",1,3,25.0,424,145.0,125.00,16.45,5.0,10,o,0));
      dataArray.push_back(DataInfo(436,"beam_436_4",1,3,25.0,424,145.0,125.00,16.45,4.0,10,o,0));
      dataArray.push_back(DataInfo(436,"beam_436_3",1,3,25.0,424,145.0,125.00,16.45,3.0,10,o,0));				    
    }
      
    study[437]= "X step; Pdf study at theta=25 deg; phi=0 deg @ 7 GeV/c. Bar, 3LS lens, all grease";
    {
      // study id | run name | radiator | lens | angle | z pos | x pos | x step | y step | momentum | beam dim | sim offset | phi 
      dataArray.push_back(DataInfo(437,"beam_437_25",1,3,25.0,424,145.0,125.00,16.45,7.0,10,o,0)); // high stat
      dataArray.push_back(DataInfo(437,"beam_18224151405",1,3,25.1,424,145.0,125.00,16.45,7.0,10,o,0));
      dataArray.push_back(DataInfo(437,"beam_18224150822",1,3,25.2,424,145.0,125.00,16.45,7.0,10,o,0));
      dataArray.push_back(DataInfo(437,"beam_18224150252",1,3,25.5,424,145.0,125.00,16.45,7.0,10,o,0));
      dataArray.push_back(DataInfo(437,"beam_18224145747",1,3,26.0,424,145.0,125.00,16.45,7.0,10,o,0));
      dataArray.push_back(DataInfo(437,"beam_18224161452",1,3,27.0,424,145.0,125.00,16.45,7.0,10,o,0));
      dataArray.push_back(DataInfo(437,"beam_18224152210",1,3,24.9,424,145.0,125.00,16.45,7.0,10,o,0));
      dataArray.push_back(DataInfo(437,"beam_18224152845",1,3,24.8,424,145.0,125.00,16.45,7.0,10,o,0));
      dataArray.push_back(DataInfo(437,"beam_18224153447",1,3,24.5,424,145.0,125.00,16.45,7.0,10,o,0));
      dataArray.push_back(DataInfo(437,"beam_18224154326",1,3,24.0,424,145.0,125.00,16.45,7.0,10,o,0));
      dataArray.push_back(DataInfo(437,"beam_18224155333",1,3,23.0,424,145.0,125.00,16.45,7.0,10,o,0));
      dataArray.push_back(DataInfo(437,"beam_18224155910",1,3,22.0,424,145.0,125.00,16.45,7.0,10,o,0));
      //bar-prism misaligment
      dataArray.push_back(DataInfo(437,"beam_18224170908",1,3,25.0,424,145.0,125.00,16.45,7.0,10,o,0)); //0.5mrad 
      dataArray.push_back(DataInfo(437,"beam_18224173000",1,3,25.0,424,145.0,125.00,16.45,7.0,10,o,0)); //1.0mrad
      dataArray.push_back(DataInfo(437,"beam_18224174751",1,3,25.0,424,145.0,125.00,16.45,7.0,10,o,0));	//2.25mrad
    }

    study[438]= "Theta scan; phi=0 deg @ 7 GeV/c. Plate, no lens, all grease";
    {
      // study id | run name | radiator | lens | angle | z pos | x pos | x step | y step | momentum | beam dim | sim offset | phi 
      dataArray.push_back(DataInfo(438,"beam_438_20",2,0,20.0,424,145.0,125.00,16.45,7.0,10,o,0));
      dataArray.push_back(DataInfo(438,"beam_438_25",2,0,25.0,424,145.0,125.00,16.45,7.0,10,o,0));
      dataArray.push_back(DataInfo(438,"beam_438_30",2,0,30.0,424,145.0,125.00,16.45,7.0,10,o,0));
      dataArray.push_back(DataInfo(438,"beam_438_40",2,0,40.0,424,145.0,125.00,16.45,7.0,10,o,0));
      dataArray.push_back(DataInfo(438,"beam_438_50",2,0,50.0,424,145.0,125.00,16.45,7.0,10,o,0));
      dataArray.push_back(DataInfo(438,"beam_438_60",2,0,60.0,424,145.0,125.00,16.45,7.0,10,o,0));
      dataArray.push_back(DataInfo(438,"beam_438_70",2,0,70.0,424,145.0,125.00,16.45,7.0,10,o,0));
      dataArray.push_back(DataInfo(438,"beam_438_80",2,0,80.0,424,145.0,125.00,16.45,7.0,10,o,0));
      dataArray.push_back(DataInfo(438,"beam_438_90",2,0,90.0,424,145.0,125.00,16.45,7.0,10,o,0));
      dataArray.push_back(DataInfo(438,"beam_438_100",2,0,100.0,424,145.0,125.00,16.45,7.0,10,o,0));
      dataArray.push_back(DataInfo(438,"beam_438_110",2,0,110.0,424,145.0,125.00,16.45,7.0,10,o,0));
      dataArray.push_back(DataInfo(438,"beam_438_120",2,0,120.0,424,145.0,125.00,16.45,7.0,10,o,0));
      dataArray.push_back(DataInfo(438,"beam_438_130",2,0,130.0,424,145.0,125.00,16.45,7.0,10,o,0));
      dataArray.push_back(DataInfo(438,"beam_438_140",2,0,140.0,424,145.0,125.00,16.45,7.0,10,o,0));
    }
    
    study[439]= "Momentum scan; theta=25 deg, phi=0 deg. Plate, no lens, all grease";
    {
      // study id | run name | radiator | lens | angle | z pos | x pos | x step | y step | momentum | beam dim | sim offset | phi 
      dataArray.push_back(DataInfo(439,"beam_439_8",1,3,25.0,424,145.0,125.00,16.45,8.0,10,o,0));
      dataArray.push_back(DataInfo(439,"beam_439_7",1,3,25.0,424,145.0,125.00,16.45,7.0,10,o,0));
      dataArray.push_back(DataInfo(439,"beam_439_6",1,3,25.0,424,145.0,125.00,16.45,6.0,10,o,0));
      dataArray.push_back(DataInfo(439,"beam_439_5",1,3,25.0,424,145.0,125.00,16.45,5.0,10,o,0));
      dataArray.push_back(DataInfo(439,"beam_439_4",1,3,25.0,424,145.0,125.00,16.45,4.0,10,o,0));
    }

    study[440]= "Theta scan; theta=25 deg, phi=0 deg @ 6 GeV/c. Plate, no lens, all grease";
    {
      // study id | run name | radiator | lens | angle | z pos | x pos | x step | y step | momentum | beam dim | sim offset | phi 
      dataArray.push_back(DataInfo(440,"beam_440_20",1,3,20.0,424,145.0,125.00,16.45,6.0,10,o,0));
      dataArray.push_back(DataInfo(440,"beam_440_30",1,3,30.0,424,145.0,125.00,16.45,6.0,10,o,0));
      dataArray.push_back(DataInfo(440,"beam_440_40",1,3,40.0,424,145.0,125.00,16.45,6.0,10,o,0));
      dataArray.push_back(DataInfo(440,"beam_440_50",1,3,50.0,424,145.0,125.00,16.45,6.0,10,o,0));
      dataArray.push_back(DataInfo(440,"beam_440_60",1,3,60.0,424,145.0,125.00,16.45,6.0,10,o,0));
      dataArray.push_back(DataInfo(440,"beam_440_70",1,3,70.0,424,145.0,125.00,16.45,6.0,10,o,0));
      dataArray.push_back(DataInfo(440,"beam_440_80",1,3,80.0,424,145.0,125.00,16.45,6.0,10,o,0));
      dataArray.push_back(DataInfo(440,"beam_440_90",1,3,90.0,424,145.0,125.00,16.45,6.0,10,o,0));
      dataArray.push_back(DataInfo(440,"beam_440_100",1,3,100.0,424,145.0,125.00,16.45,6.0,10,o,0));
      dataArray.push_back(DataInfo(440,"beam_440_110",1,3,110.0,424,145.0,125.00,16.45,6.0,10,o,0));
      dataArray.push_back(DataInfo(440,"beam_440_120",1,3,120.0,424,145.0,125.00,16.45,6.0,10,o,0));
    }

    study[441]= "Theta scan; phi=0 deg @ 7 GeV/c. Bar, 3LS lens, all grease";
    {
      // study id | run name | radiator | lens | angle | z pos | x pos | x step | y step | momentum | beam dim | sim offset | phi 
      dataArray.push_back(DataInfo(441,"beam_441_20",1,3,20.0,442,85.0,70.00,16.45,7.0,10,o));
      dataArray.push_back(DataInfo(441,"beam_441_25",1,3,25.0,442,85.0,70.00,16.45,7.0,10,o));
      //X collimator changed from +-10mm to +-5mm
      dataArray.push_back(DataInfo(441,"beam_441_20_1",1,3,25.0,442,85.0,70.00,16.45,7.0,10,o,0,1));
      //Y collimator changed from +-10mm to +-5mm, X colli still at +-5mm
      dataArray.push_back(DataInfo(441,"beam_18226135843",1,3,25.0,442,85.0,70.00,16.45,7.0,10,o,0,2));
      //no EDD no splitter
      dataArray.push_back(DataInfo(441,"beam_18226142449",1,3,25.0,442,85.0,70.00,16.45,7.0,10,o,0,3));
    }

    study[442]= "Theta scan; phi=0 deg @ 7 GeV/c. Bar, 3LS lens, all grease";
    {
      // study id | run name | radiator | lens | angle | z pos | x pos | x step | y step | momentum | beam dim | sim offset | phi 
      dataArray.push_back(DataInfo(442,"beam_442_20",1,3,20.0,442,85.0,70.00,16.45,7.0,10,o));
      dataArray.push_back(DataInfo(442,"beam_442_25",1,3,25.0,442,85.0,70.00,16.45,7.0,10,o));
      dataArray.push_back(DataInfo(442,"beam_18226173017",1,3,30.0,442,85.0,70.00,16.45,7.0,10,o));
      dataArray.push_back(DataInfo(442,"beam_18226180719",1,3,40.0,442,85.0,70.00,16.45,7.0,10,o));
      dataArray.push_back(DataInfo(442,"beam_18226183631",1,3,50.0,442,85.0,70.00,16.45,7.0,10,o));
      dataArray.push_back(DataInfo(442,"beam_18226185930",1,3,60.0,442,85.0,70.00,16.45,7.0,10,o));
      dataArray.push_back(DataInfo(442,"beam_18226192230",1,3,70.0,442,85.0,70.00,16.45,7.0,10,o));
      dataArray.push_back(DataInfo(442,"beam_18226194717",1,3,80.0,442,85.0,70.00,16.45,7.0,10,o));
      dataArray.push_back(DataInfo(442,"beam_442_90",1,3,90.0,442,85.0,70.00,16.45,7.0,10,o));
      dataArray.push_back(DataInfo(442,"beam_18226205028",1,3,100.0,442,85.0,70.00,16.45,7.0,10,o));
    }
    
    study[450]="Theta scan; phi=0; bar + grease + 3LS; mcp3,4 have custom CS";
    for(auto i=20; i<150; i=i+5){
      dataArray.push_back(DataInfo(450,Form("beam_450_%d",i),1,3,i,447,85.0,66.80,14.8,7.0,3,0,0.23));
    }

    study[451]="Theta scan; phi=0; bar + grease + 3LS ";    
    for(int i=20; i<=140; i=i+5){
      dataArray.push_back(DataInfo(451,Form("beam_451_%d",i),1,3,i,442,85.0,66.80,14.8,7.0,3,0,0));
    }
    
    study[452]="Theta scan; phi=5 bar + grease + 3LS ";
    for(int i=20; i<=140; i=i+5){
      dataArray.push_back(DataInfo(452,Form("beam_452_%d",i),1,3,i,442,85.0,66.80,14.8,7.0,3,0,5));
    }
    
    study[453]="Theta scan; phi=10 bar + grease + 3LS ";
    for(int i=20; i<=140; i=i+5){
      dataArray.push_back(DataInfo(453,Form("beam_453_%d",i),1,3,i,442,85.0,66.80,14.8,7.0,3,0,10));
    }
    
    study[454]="Theta scan; phi=15 bar + grease + 3LS ";
    for(int i=20; i<=140; i=i+5){
      dataArray.push_back(DataInfo(454,Form("beam_454_%d",i),1,3,i,442,85.0,66.80,14.8,7.0,3,0,15));
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

    aliasArray.push_back(dataArray[i]);
	  
    
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

  //std::sort(aliasArray.begin(), aliasArray.end());

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
      std::cout<< newset[i].getOpt()<<std::endl;
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
      std::cout<<"\"/d/proc/jul18/"<<newset[i].getStudyId()<<"/"<<newset[i].getChildRunId(0)<<"C.root\","
	       << newset[i].getStudyId() <<","<<i<<","<<newset[i].getMomentum() << ","<<newset[i].getRadiatorId()
	       <<","<<newset[i].getLensId()
	       <<","<<newset[i].getAngle()<<","<<newset[i].getZ()
	       <<","<<newset[i].getX()<<","<<newset[i].getXstep()
	       <<","<<newset[i].getYstep()<<std::endl;
    }
  }
  
  if(format==7){ // reco sim
    for(UInt_t i = 0; i != newset.size(); i++) {
      std::cout<<"\"/d/proc/jul18/"<<newset[i].getStudyId()<<"/"<<newset[i].getChildRunId(0)<<"S.root\","
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
    if(sname == name) return aliasArray[i];
  }
  return DataInfo();
}

void datainfo(Int_t studyId=0, Int_t format = 0){
  datainfo_init();
  createAliases();

  std::vector<DataInfo> newset= getStudy(studyId);

  if(format==8){
    for(Int_t i=studyId; i<gg_nstudies; i++){
      std::vector<DataInfo> newset= getStudy(i);
      if(newset.size()>0){
	p_print(newset, format);
	std::cout<<Form("desc[%d]=\"",i)<<study[i]<<"\";" <<std::endl;
      }
    }
    for(Int_t i=studyId; i<gg_nstudies; i++){
      std::vector<DataInfo> newset= getStudy(i);
      if(newset.size()>0){
	std::cout<<Form("<option value=%d> %d ",i,i)<<study[i]<<"</option>" <<std::endl;
      }
    }
  }else if(format==10){ //cp all
    for(Int_t i=0; i<gg_nstudies; i++){
      if(studyId==0 || i==studyId){
	std::vector<DataInfo> newset= getStudy(i);
	if(newset.size()>0){
	  std::cout<<"mkdir /d/proc/jul18/"<<i <<std::endl;
	  std::cout<<"cp ";
	  p_print(newset, format);
	  std::cout<<"/d/proc/jul18/"<<i<<std::endl;
	}
      }
    }
  }else if(format==11){ //print all ID starting from studyId
    for(Int_t i=studyId; i<gg_nstudies; i++){
      std::vector<DataInfo> newset= getStudy(i);
      if(newset.size()>0){
	std::cout<<i<<std::endl;
      }
    }
  }else if(format==12) p_exportinfo(); // html
  else p_print(newset, format);

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
