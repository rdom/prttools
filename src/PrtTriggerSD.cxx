#include "PrtTriggerSD.h"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4RunManager.hh"
#include "G4TransportationManager.hh"
#include <TVector3.h>

#include "PrtEvent.h"

#include "PrtRunAction.h"
#include "PrtManager.h"

PrtTriggerSD::PrtTriggerSD(
                            const G4String& name, 
                            const G4String& hitsCollectionName,
                            G4int nofCells)
  : G4VSensitiveDetector(name)
{
  collectionName.insert(hitsCollectionName);
}

PrtTriggerSD::~PrtTriggerSD() 
{ 
}

void PrtTriggerSD::Initialize(G4HCofThisEvent* hce)
{ 

}

G4bool PrtTriggerSD::ProcessHits(G4Step* aStep, G4TouchableHistory* hist)
{   
  // energy deposit
  G4double edep = aStep->GetTotalEnergyDeposit();
  G4ThreeVector globalpos = aStep->GetPostStepPoint()->GetPosition();
  G4Track* track = aStep->GetTrack(); 
  G4ThreeVector g4pos = track->GetVertexPosition();
  
  PrtManager::Instance()->Event()->SetTrigger(1);
  //PrtManager::Instance()->Event()->SetPosition(TVector3(globalpos.x(),globalpos.y(),globalpos.z()));

  return true;
}

void PrtTriggerSD::EndOfEvent(G4HCofThisEvent*)
{ 
 
}

