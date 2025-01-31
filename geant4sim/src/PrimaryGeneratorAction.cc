// PrimaryGeneratorAction.cc
#include "PrimaryGeneratorAction.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4DecayWithSpin.hh"
#include "G4ProcessManager.hh"
#include "Randomize.hh"
#include "TVector3.h"
#include "TMath.h"
#include "file.h"

namespace B1 {

  PrimaryGeneratorAction::PrimaryGeneratorAction()  {

    particleGun = new G4ParticleGun(1);
    
    // Set positron as the particle
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* particle = particleTable->FindParticle("e+");
    //G4ParticleDefinition* particle = G4MuonPlus::Definition();
    
    // Initialize muon decay
    //G4DecayWithSpin* muonDecayProcess = new G4DecayWithSpin();
    //G4ProcessManager* muonProcessManager = particle->GetProcessManager();
    //muonProcessManager->AddProcess(muonDecayProcess, 1, 1, 1);
    //muonProcessManager->SetProcessOrdering(muonDecayProcess, G4ProcessVectorDoItIndex::idxAtRest);
    //muonProcessManager->SetProcessOrdering(muonDecayProcess, G4ProcessVectorDoItIndex::idxAlongStep);
    //muonProcessManager->SetProcessOrdering(muonDecayProcess, G4ProcessVectorDoItIndex::idxPostStep);

    particleGun->SetParticleDefinition(particle);
  }
  
  PrimaryGeneratorAction::~PrimaryGeneratorAction() {
    delete particleGun;
  }
  
  void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event) {
    
    // Let's generate a rotating muon
    /*
    G4double radius = 3.0 * cm;
    G4double theta = G4UniformRand() * 2 * M_PI;
    G4double x0 = radius * std::cos(theta);
    G4double y0 = radius * std::sin(theta);
    G4double z0 = 0.0;
    G4double ux = std::cos(theta - M_PI / 2.);
    G4double uy = std::sin(theta - M_PI / 2.);
    G4double uz = 0.;
    G4double momentum = 0.0 * MeV; // monochromatic muons on 3.0 cm orbit

    particleGun->SetParticleMomentumDirection(G4ThreeVector(ux, uy, uz));
    // Set momentum magnitude
    particleGun->SetParticleMomentum(momentum * G4ThreeVector(ux, uy, uz));
    
    particleGun->SetParticlePosition(G4ThreeVector(x0, y0, z0));
    
    // Set polarization along muon direction
    particleGun->SetParticlePolarization(G4ThreeVector(ux, uy, uz));

    // Initialize dummy track
    track = HelixTrack(0., 0., 0.);
    track.trueMomentum = -999.;

    G4cout << "Generated a muon " << G4endl;
    
    particleGun->GeneratePrimaryVertex(event);
    */
    // Here we are generating positrons
    
    //G4cout << "Initializing a new muon " << G4endl;
    track = HelixTrack(0., 0., 0.);

    // Set momentum in the range [0; 68.9] MeV/c
    G4double momentum = G4UniformRand() * 64.7 * MeV;
    track.trueMomentum = momentum;
    
    // Set position randomly in a circle of radius 3 cm lying in the x-y plane with z = 0
    G4double radius = 3.1 * cm;
    G4double theta = G4UniformRand() * 2 * M_PI;
    G4double x0 = radius * std::cos(theta);
    G4double y0 = radius * std::sin(theta);
    G4double z0 = 0.0;
    track.addOrigin(TVector3(x0, y0, z0));
    
    // Set random emission angle on a sphere
    G4double costheta = 2 * G4UniformRand() - 1;
    G4double sintheta = std::sqrt(1 - costheta * costheta);
    G4double phi = G4UniformRand() * 2 * M_PI;
    G4double ux = sintheta * std::cos(phi);
    G4double uy = sintheta * std::sin(phi);
    G4double uz = costheta;

    double azimuthalAngle = std::acos(costheta);
    track.azimuthalAngle = azimuthalAngle;
    track.polarAngle = phi;
    //G4cout << "Theta = " << azimuthalAngle << " Phi = " << phi << G4endl;
    // Define emission angle theta_e = between 0 and 2pi, 0 being along tangent direction, pi/2 along z
    TVector3 radial_axis = {TMath::Cos(theta), TMath::Sin(theta), 0.};
    TVector3 emission_axis = {ux, uy, uz};
    track.emissionAngle = (uz > 0) ? radial_axis.Angle(emission_axis) : - radial_axis.Angle(emission_axis);

    track.spinAngle = phi + theta;

    particleGun->SetParticleMomentumDirection(G4ThreeVector(ux, uy, uz));
    
    // Set momentum magnitude
    particleGun->SetParticleMomentum(momentum * G4ThreeVector(ux, uy, uz));
    
    particleGun->SetParticlePosition(G4ThreeVector(x0, y0, z0));
  
    particleGun->GeneratePrimaryVertex(event);
    

  }    
  
}
