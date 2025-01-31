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
/// \file B1/src/RunAction.cc
/// \brief Implementation of the B1::RunAction class

#include "RunAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"
#include "file.h"

#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4AccumulableManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "HelixTrack.h"
#include "TString.h"

namespace B1

{
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  RunAction::RunAction() {}

  void RunAction::BeginOfRunAction(const G4Run*) {
    // Open ROOT file at the beginning of the run
    //InitializeTTree();
    tracks.clear();

  }

  void RunAction::EndOfRunAction(const G4Run* run) {
    // Close ROOT file at the end of the run
    //tree->SetEntries();
    //tree->Write();
    //file->Write();
    //file->Close();
    //delete rootFile;
    int id = run->GetRunID();
    G4cout << "Writing output file of run " << id << G4endl;
    writeToFile(Form("chet_sim_z20_FullGeo_7Cyl_Nopetals_%d.root", id + 3), tracks);
  }

  //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  void RunAction::writeToFile(const char* filename, std::vector<HelixTrack>& trkx) {
    TFile outputFile(filename, "RECREATE");
    if (!outputFile.IsOpen()) {
      std::cerr << "Error: Could not open ROOT file for writing: " << filename << std::endl;
      return;
    }

    TTree tree("HelixTrackTree", "Tree containing HelixTrack objects");

    // Define TBranches corresponding to each member of the HelixTrack class
    double trueMomentum;
    double polarAngle;
    double azimuthalAngle;
    double spinAngle;
    double emissionAngle;
    TVector3 origin;
    std::vector<std::vector<double>> hitsCoordinates;
    std::vector<std::vector<double>> trackCoordinates;
    std::vector<int> planeID;

    tree.Branch("trueMomentum", &trueMomentum);
    tree.Branch("polarAngle", &polarAngle);
    tree.Branch("azimuthalAngle", &azimuthalAngle);
    tree.Branch("spinAngle", &spinAngle);
    tree.Branch("emissionAngle", &emissionAngle);
    tree.Branch("origin", &origin);
    tree.Branch("hitsCoordinates", &hitsCoordinates);
    tree.Branch("trackCoordinates", &trackCoordinates);
    tree.Branch("planeID", &planeID);

    // Fill the TTree with data from each HelixTrack object
    for (auto & tracki : trkx) {
      trueMomentum = tracki.trueMomentum;
      polarAngle = tracki.polarAngle;
      azimuthalAngle = tracki.azimuthalAngle;
      spinAngle = tracki.spinAngle;
      emissionAngle = tracki.emissionAngle;
      origin = tracki.origin;
      hitsCoordinates = tracki.hitsCoordinates;
      trackCoordinates = tracki.trackCoordinates;
      planeID = tracki.planeID;
      tree.Fill();
    }

    // Write the TTree to the ROOT file and close the file
    tree.Write();
    outputFile.Close();
  }

}
