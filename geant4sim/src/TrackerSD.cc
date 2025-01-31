//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id$
//
/// \file TrackerSD.cc
/// \brief Implementation of the TrackerSD class

#include "TrackerSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "HelixTrack.h"
#include "file.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

namespace B1 {

  TrackerSD::TrackerSD(const G4String& name,
			   const G4String& hitsCollectionName) 
    : G4VSensitiveDetector(name),
      fHitsCollection(NULL)
  {
    collectionName.insert(hitsCollectionName);
  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  
  TrackerSD::~TrackerSD() 
  {}

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  
  void TrackerSD::Initialize(G4HCofThisEvent* hce)
  {
    // Create hits collection

    fHitsCollection 
      = new TrackerHitsCollection(SensitiveDetectorName, collectionName[0]); 

    // Add this collection in hce

    G4int hcID 
      = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
    hce->AddHitsCollection( hcID, fHitsCollection ); 
  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  
  G4bool TrackerSD::ProcessHits(G4Step* aStep, 
				  G4TouchableHistory*)
  {  
    //G4cout << "Processing hits " << G4endl;
    // energy deposit

    G4double edep = aStep->GetTotalEnergyDeposit();
    
    if (edep==0.) return false;
    
    TrackerHit* newHit = new TrackerHit();
    
    newHit->SetTrackID  (aStep->GetTrack()->GetTrackID());
    newHit->SetplaneID(aStep->GetPreStepPoint()->GetTouchableHandle()
			 ->GetCopyNumber());
    newHit->SetEdep(edep);
    newHit->SetPos (aStep->GetPostStepPoint()->GetPosition());
    
    fHitsCollection->insert( newHit );
    
    // For muon generation
    /*
    const std::vector<const G4Track*>* secondary = aStep->GetSecondaryInCurrentStep();
    for (int i=0; i < secondary->size(); i++) {
      const G4Track* particle = secondary->at(i);
      G4cout << particle->GetDefinition()->GetParticleName() << G4endl;
      if (particle->GetDefinition()->GetParticleName() != "e+") continue;
      // Initialize correctly the positron the first time
      if (track.trueMomentum < 0.) {
	const G4ThreeVector vertexMom = particle->GetVertexMomentumDirection();
	const G4ThreeVector vertexPos = particle->GetVertexPosition();
	track.trueMomentum = vertexMom.mag();
	G4cout << "Found a positron with energy " << track.trueMomentum << G4endl;
	track.azimuthalAngle = vertexMom.theta();
	track.polarAngle = vertexMom.phi();
	double spinAngle = vertexMom.phi() + vertexPos.phi();
	track.spinAngle = spinAngle;
	if (std::sin(spinAngle) > 0) {
	  track.emissionAngle = (vertexMom.theta() < M_PI / 2.) ? M_PI / 2. - vertexMom.theta() : - (M_PI / 2. + vertexMom.theta());
	}
	else {
	  track.emissionAngle = (vertexMom.theta() < M_PI / 2.) ? M_PI / 2. + vertexMom.theta() : - (M_PI / 2. - vertexMom.theta());
	}
	track.addOrigin(TVector3(vertexPos.x(), vertexPos.y(), vertexPos.z()));
	
      }
      
      // Add hit informations
      G4ThreeVector position = aStep->GetPostStepPoint()->GetPosition();
      track.addHit((double)position.x()*0.1, (double)position.y()*0.1, (double)position.z()*0.1, aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber());
      
      // Skip other positrons
      break;
    }
    //newHit->Print();
    */

    // For positron generation

    
    // Add hits to the helix track event if they belong to the positron
    if (aStep->GetTrack()->GetTrackID() == 1) {
      //G4cout << "Adding a hit " << G4endl;
      G4ThreeVector position = aStep->GetPostStepPoint()->GetPosition();
      track.addHit((double)position.x()*0.1, (double)position.y()*0.1, (double)position.z()*0.1, 
		   aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber());
      // check the current number of hits: kill particles making too many turns
      if (track.hitsCoordinates.size() > 100) {
	G4Track* theTrack = aStep->GetTrack();
	theTrack->SetTrackStatus(fStopAndKill);
	return true;
      }
    }
    
    else { // Do not track other particles
      //G4cout << "Killing a particle " << G4endl;
      G4Track* theTrack = aStep->GetTrack();
      theTrack->SetTrackStatus(fStopAndKill);
      return true;
    }
    

    return true;
  }
  
  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  
  void TrackerSD::EndOfEvent(G4HCofThisEvent*)
  {
    if ( verboseLevel>1 ) { 
      G4int nofHits = fHitsCollection->entries();
      G4cout << "\n-------->Hits Collection: in this event they are " << nofHits 
	     << " hits in the tracker chambers: " << G4endl;
      for ( G4int i=0; i<nofHits; i++ ) (*fHitsCollection)[i]->Print();
    }
  }
  
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  
}
