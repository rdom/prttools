#include "PrtHit.h"

ClassImp(PrtHit)

// -----   Default constructor   -------------------------------------------
PrtHit::PrtHit(){ 
  fType=-1;
  fCherenkovMC = 0;
  
  // fLocalPos = TVector3(0,0,0);
  // fGlobalPos = TVector3(0,0,0);
  // fDigiPos = TVector3(0,0,0);
  // fMomentum = TVector3(0,0,0);

  fMcpId=-1;
  fPixelId=-1;
  fChannel= -1;
  fTdc = -1;
  fTrb = -1;
  fMultiplicity = 0;
  fLeadTime = -1; 
  fTotTime = -1; 
}
