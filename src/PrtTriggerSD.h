// -----------------------------------------
// PrtTriggerSD.h
//
// Author  : R.Dzhygadlo at gsi.de
// -----------------------------------------

#ifndef PrtTriggerSD_h
#define PrtTriggerSD_h 1

#include <vector>
#include "G4VSensitiveDetector.hh"

#include "PrtEvent.h"

class G4Step;
class G4HCofThisEvent;

class PrtTriggerSD : public G4VSensitiveDetector
{
public:
  PrtTriggerSD(const G4String& name, 
	     const G4String& hitsCollectionName, 
	     G4int nofCells);
  virtual ~PrtTriggerSD();
  
  // methods from base class
  virtual void   Initialize(G4HCofThisEvent* hitCollection);
  virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* history);
  virtual void   EndOfEvent(G4HCofThisEvent* hitCollection);

};


#endif
