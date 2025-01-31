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
//
/// \file B1/src/SteppingAction.cc
/// \brief Implementation of the B1::SteppingAction class

#include "SteppingAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"

#include "G4Track.hh"
#include "G4Step.hh"
#include "G4Event.hh"
#include "G4ThreeVector.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "TMath.h"
#include "file.h"

namespace B1
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  SteppingAction::SteppingAction(EventAction* event) : event(event) {}

  SteppingAction::~SteppingAction() {}

  void SteppingAction::UserSteppingAction(const G4Step* step) {
    G4cout << "In a stepping action" << G4endl;
    // check that we are following the positron
    G4Track* aTrack = step->GetTrack();
    if (aTrack->GetParticleDefinition()->GetParticleName() != "e+" ||
	aTrack->GetTrackID() != 1) {
      aTrack->SetTrackStatus(fStopAndKill);
      return;
    }
    // Check if the particle has entered a detector
    if (step->IsFirstStepInVolume()) {
    //if (step->GetPostStepPoint()->GetStepStatus() == fGeomBoundary) {
      // Get position of the particle at the current step
      G4ThreeVector position = step->GetPostStepPoint()->GetPosition();

      // calculate plane ID
      int id = 0;
      if (TMath::Hypot((double)position.x()*0.1, (double)position.y()*0.1) > 4.0 && TMath::Hypot((double)position.x()*0.1, (double)position.y()*0.1) < 10.0) {
	// id for petals
	id = (int)((TMath::ATan2((double)position.y(), (double)position.x()) + TMath::Pi())*(180./(TMath::Pi()))/36.);
      }
      else if (TMath::Hypot((double)position.x()*0.1, (double)position.y()*0.1) < 4.0) {
	// id for endplates
	if ((double)position.z() < 0) {
	  id = 36 + (int)((TMath::Abs((double)position.z()*0.1) - 5.)/0.5);
	}
	else {
	  id = 40 + (int)((TMath::Abs((double)position.z()*0.1) - 5.)/0.5);
	}
      }

      // Store position information into the helix track
      //G4cout << "Hit at position: " << position.x() << ", " << position.y() << ", " << position.z() << G4endl; 
      if (TMath::Hypot((double)position.x()*0.1, (double)position.y()*0.1) < 10.0) {
	track.addHit((double)position.x() * 0.1, (double)position.y() * 0.1, (double)position.z() * 0.1, id); // in cm
      }
    }
    G4ThreeVector position = step->GetPostStepPoint()->GetPosition();
    track.addTrackPoint((double)position.x() * 0.1, (double)position.y() * 0.1, (double)position.z() * 0.1); // in cm
  }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
