#include <iostream>
#include <vector>
#include <cmath>

#include "ConstField.h"
#include "Exception.h"
#include "FieldManager.h"
#include "KalmanFitterRefTrack.h"
#include "StateOnPlane.h"
#include "Track.h"
#include "TrackCand.h"
#include "TrackPoint.h"

#include "MeasurementProducer.h"
#include "MeasurementFactory.h"

//#include "mySpacepointDetectorHit.h"
//#include "mySpacepointMeasurement.h"

#include "SpacepointMeasurement.h"
#include "PlanarMeasurement.h"
#include "MaterialEffects.h"
#include "RKTrackRep.h"
#include "TGeoMaterialInterface.h"

#include <TTree.h>
#include <TChain.h>
#include <TFile.h>
#include <TProfile.h>
#include <TRandom.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TPad.h>
//#include "TRectangle.h"
#include <TEfficiency.h>
#include <TH2F.h>

#include <EventDisplay.h>

#include <TDatabasePDG.h>
#include <TEveManager.h>
#include <TGeoManager.h>
#include <TVector3.h>

#include "TDatabasePDG.h"
#include <TMath.h>

#include "NonUniformBField.h"

// Constants
const double PI = 3.14159265358979323846;

TMatrixDSym CovarianceToCylinder(TMatrixDSym cov, TVector3 mom) {
  // Transform a covariance matrix in x,y,z, momx, momy, momz
  // into a covariance matrix in x, y, z, mom, theta, phi
  TMatrixDSym covCyl(cov);

  TMatrixD Jac(6, 6);
  Jac.Zero();

  double p = mom.Mag();
  double pt = TMath::Hypot(mom.X(), mom.Y());

  if (p == 0 || pt == 0) return covCyl;

  // Calculate Jacobian
  Jac[0][0] = 1.;
  Jac[1][1] = 1.;
  Jac[2][2] = 1.;

  Jac[3][3] = mom.X() / p;
  Jac[3][4] = mom.Y() / p;
  Jac[3][5] = mom.Z() / p;

  Jac[4][3] = mom.X() * mom.Z() / p / p / pt;
  Jac[4][4] = mom.Y() * mom.Z() / p / p / pt;
  Jac[4][5] = - pt / p / p;

  Jac[5][3] = - mom.Y() / pt / pt;
  Jac[5][4] = mom.X() / pt / pt;
  Jac[5][5] = 0.;

  covCyl.Similarity(Jac);
  return covCyl;
  
}

std::vector<std::vector<double>> sort_vector_z(std::vector<std::vector<double>> vectors) {
  // sort vectors according to z(third) value

  std::vector<double> z_of_vectors;
  for (auto v : vectors) {
    z_of_vectors.push_back(v.at(2));
  }
  
  // initialize original index locations
  vector<size_t> idx(z_of_vectors.size());
  iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  // using std::stable_sort instead of std::sort
  // to avoid unnecessary index re-orderings
  // when v contains elements of equal values 
  std::stable_sort(idx.begin(), idx.end(),
       [&z_of_vectors](size_t i1, size_t i2) {return z_of_vectors[i1] < z_of_vectors[i2];});

  std::vector<std::vector<double>> copy_vector;
  for (auto i : idx) {
    copy_vector.push_back(vectors.at(i));
  }
  return copy_vector;
}

std::vector<int> sort_vector_id(std::vector<std::vector<double>> vectors, std::vector<int> vector_id) {
  // sort vectors according to z(third) value

  std::vector<double> z_of_vectors;
  for (auto v : vectors) {
    z_of_vectors.push_back(v.at(2));
  }
  
  // initialize original index locations
  vector<size_t> idx(z_of_vectors.size());
  iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  // using std::stable_sort instead of std::sort
  // to avoid unnecessary index re-orderings
  // when v contains elements of equal values 
  std::stable_sort(idx.begin(), idx.end(),
       [&z_of_vectors](size_t i1, size_t i2) {return z_of_vectors[i1] < z_of_vectors[i2];});

  std::vector<int> copy_vector;
  for (auto i : idx) {
    copy_vector.push_back(vector_id.at(i));
  }
  return copy_vector;
}

// Draw the hits in xy view
TGraph* drawXYView_hits(std::vector<std::vector<double>> hitsCoordinates, bool eff) {
  int nHits = hitsCoordinates.size();
  double* x = new double[nHits];
  double* y = new double[nHits];
  for (int i = 0; i < nHits; ++i) {
    x[i] = hitsCoordinates[i].at(0);
    y[i] = hitsCoordinates[i].at(1);
  }
  TGraph* gr = new TGraph(nHits, x, y);
  gr->SetMarkerStyle(20);
  if (eff) {
    gr->SetLineColor(kBlue);
  }
  else {
    gr->SetLineColor(kRed);
  }
  return gr;
}
  
  // Draw the hits in yz view
TGraph* drawYZView_hits(std::vector<std::vector<double>> hitsCoordinates, bool eff) {
  int nHits = hitsCoordinates.size();
  double* y = new double[nHits];
  double* z = new double[nHits];
  for (int i = 0; i < nHits; ++i) {
    y[i] = hitsCoordinates[i].at(1);
    z[i] = hitsCoordinates[i].at(2);
    }
  TGraph* gr = new TGraph(nHits, z, y);
  gr->SetMarkerStyle(20);
  if (eff) {
    gr->SetLineColor(kBlue);
  }
  else {
    gr->SetLineColor(kRed);
  }
  return gr;
}

TH2F *h_polarAngle_bias = new TH2F("h_polarAngle_bias", ";#phi [rad]; Bias #phi_{true} - #phi_{fit} [rad]", 20, 0., 2 * PI, 50, -0.25, 0.25);
TH2F *h_polarAngle_std = new TH2F("h_polarAngle_std", ";#phi [rad]; #sigma #phi_{fit} [rad]", 20, 0., 2 * PI, 50, 0., 0.05);
TH2F *h_azimuthalAngle_bias = new TH2F("h_azimuthalAngle_bias", ";#theta^{geo} [rad]; Bias #theta^{geo}_{true} - #theta^{geo}_{fit} [rad]", 10, 0., PI, 40, - 0.1, 0.1);
TH2F *h_azimuthalAngle_std = new TH2F("h_azimuthalAngle_std", ";#theta [rad]; #sigma #theta_{fit} [rad]", 20, -PI, PI, 20, 0., 0.1);
TH2F *h_mom_bias = new TH2F("h_mom_bias", "; Momentum [MeV]; Bias Mom_{true} - Mom_{fit} [MeV]", 10, 0., 68.9, 40, -10, 10);
TH2F *h_mom_std = new TH2F("h_mom_std", "; Momentum [MeV]; #sigma Mom_{fit} [MeV]", 10, 0., 68.9, 100, 0, 20);
TH2F *h_XOrbit_bias = new TH2F("h_XOrbit_bias", "; x_{0} (on orbit) [cm]; Bias x_{0, true} - x_{0, fit} [cm]", 60, -3.0, 3-0, 20, -0.5, 0.5);
TH2F *h_XOrbit_std = new TH2F("h_XOrbit_std", "; x_{0} (on orbit) [cm]; #sigma x_{0, fit} [cm]", 60, -3.0, 3-0, 20, 0., 0.5);
TH2F *h_YOrbit_bias = new TH2F("h_YOrbit_bias", "; y_{0} (on orbit) [cm]; Bias y_{0, true} - y_{0, fit} [cm]", 60, -3.0, 3-0, 20, -0.5, 0.5);
TH2F *h_YOrbit_std = new TH2F("h_YOrbit_std", "; y_{0} (on orbit) [cm]; #sigma y_{0, fit} [cm]", 60, -3.0, 3-0, 20, 0., 0.5);

TCanvas* canvas = new TCanvas("canvas", "Muon Decay Simulation", 1200, 600);

TEfficiency* detector_eff_phi = new TEfficiency("detEffPhi", "Tracking Efficiency; momentum [MeV]; Angle in orbit plane #phi [rad]", 10, 0, 68.9, 10, 0., 2 * PI);
TEfficiency* detector_eff_theta = new TEfficiency("detEffTheta", "Tracking Efficiency; momentum [MeV]; Emission Angle #theta [rad]", 10, 0, 68.9, 10, 0., PI);
TEfficiency* detector_acc_phi = new TEfficiency("detAccPhi", "Tracker Acceptance; momentum [MeV]; Angle in orbit plane #phi [rad]", 10, 0, 68.9, 10, 0., 2 * PI);
TEfficiency* detector_acc_theta = new TEfficiency("detAccTheta", "Tracker Acceptance; momentum [MeV]; Emission Angle #theta [rad]", 10, 0, 68.9, 10, 0., PI);

void chet_genfit() {

  canvas->Divide(2, 1);
  
  // Load events
  TChain *tracks = new TChain("HelixTrackTree");
  for (int i=0; i<4; i++) {
    //tracks->Add(Form("chet_sim_z8_with4CylinderEndCaps_%d.root", i));
    tracks->Add(Form("chet_sim_z20_FullGeo_7Cyl_Nopetals_%d.root", i));
  }
  std::cout << tracks->GetEntries() << std::endl;
  //TFile *file = TFile::Open(filename);
  //TTree *tracks = file->Get<TTree>("HelixTrackTree");
  
  // Define TBranches corresponding to each member of the HelixTrack class
  double trueMomentum;
  double polarAngle;
  double azimuthalAngle;
  double spinAngle;
  double emissionAngle;
  TVector3* origin = 0;
  std::vector<std::vector<double>>* hitsCoordinates = 0;
  std::vector<std::vector<double>>* trackCoordinates = 0;
  std::vector<int>* planeID = 0;

  TBranch* bmom = 0;
  TBranch* bpol = 0;
  TBranch* baz = 0;
  TBranch* bspin = 0;
  //TBranch* bemission = 0;
  TBranch* borigin = 0;
  TBranch* bhits = 0;
  TBranch* btrack = 0;
  TBranch *bid = 0;

  tracks->SetBranchAddress("trueMomentum", &trueMomentum, &bmom);
  tracks->SetBranchAddress("polarAngle", &polarAngle, &bpol);
  tracks->SetBranchAddress("azimuthalAngle", &azimuthalAngle, &baz);
  tracks->SetBranchAddress("spinAngle", &spinAngle, &bspin);
  //tracks->SetBranchAddress("emissionAngle", &emissionAngle, &bemission);
  tracks->SetBranchAddress("origin", &origin, &borigin);
  tracks->SetBranchAddress("hitsCoordinates", &hitsCoordinates, &bhits);
  tracks->SetBranchAddress("trackCoordinates", &trackCoordinates, &btrack); 
  tracks->SetBranchAddress("planeID", &planeID, &bid); 
    
  // Output TTree for storing tracks informations
  TTree *outTree = new TTree("fittedTracks", "Fitted Tracks State");
  TVectorD parVar(6);
  TVectorD parFit(6);
  TVectorD parTrue(6);
  double chisquare;
  int ngoodhits;

  outTree->Branch("parVar", "TVectorD", &parVar);
  outTree->Branch("parFit", "TVectorD", &parFit);
  outTree->Branch("parTrue", "TVectorD", &parTrue);
  outTree->Branch("chisquare", &chisquare);
  outTree->Branch("ngoodhits", &ngoodhits);
  
  // init geometry and mag. field
  new TGeoManager("DetectorGeometry", "CHET geometry");
  //TGeoManager::Import("detectorGeometry_with4CylinderEndCaps_design2.root");
  TGeoManager::Import("detectorGeometry_z20_fullGeo_7Cyl_Nopetals.root");
  genfit::MaterialEffects::getInstance()->init(new genfit::TGeoMaterialInterface());
  double B = 22.0; // kGaus // 2.2 T
  genfit::FieldManager::getInstance()->init(new genfit::ConstField(0. , 0., B));
  //std::shared_ptr<NonUniformBField> bField = std::make_shared<NonUniformBField>("fieldMap.txt");
  //genfit::FieldManager::getInstance()->init(bField.get()); // Non constant field
  // particle id in pdg
  const int pdg = -11; // positron
  
  // init event display
  //genfit::EventDisplay* display = genfit::EventDisplay::getInstance();


  // FIT WITH CANDIDATE TRACK
  /*
  TClonesArray myDetectorHitArray("genfit::mySpacepointDetectorHit");
  
  // init the factory
  int myDetId(3);
  genfit::MeasurementFactory<genfit::AbsMeasurement> factory;
  genfit::MeasurementProducer<genfit::mySpacepointDetectorHit, genfit::mySpacepointMeasurement> myProducer(&myDetectorHitArray);
  factory.addProducer(myDetId, &myProducer);
  */

  int drawn = 0;
  
  int nEvents = tracks->GetEntriesFast();
  std::cout << nEvents << std::endl;

  int fitted = 0;
  int in_acceptance = 0;
  
  for (int iev=0; iev<nEvents; iev++) {

    tracks->GetEntry(iev);

    int nhits = hitsCoordinates->size();

    
    // Emission angle with respect to z axis: it is pi/2 aligned to it. Goes from -pi -> pi, needs to know if trajectory is in-going or out-going
    // with respect to muon orbit.
    double decay_angle_orbit_plane = TMath::ATan2(origin->Y(), origin->X());
    double r = TMath::Hypot(origin->X(), origin->Y());
    TVector3 radial_axis = {TMath::Cos(decay_angle_orbit_plane), TMath::Sin(decay_angle_orbit_plane), 0.};
    TVector3 positron_vector = {TMath::Sin(azimuthalAngle)*TMath::Cos(polarAngle), TMath::Sin(azimuthalAngle)*TMath::Sin(polarAngle), TMath::Cos(azimuthalAngle)};
    int is_decay_inner = (TMath::Hypot(origin->X() / r - TMath::Sin(polarAngle), TMath::Cos(polarAngle) + origin->Y() / r) < TMath::Hypot(origin->X() / r - origin->Y() / r, origin->X() / r + origin->Y() / r)) ? 1 : 0;
    emissionAngle = 0;
    if (TMath::Cos(azimuthalAngle) >= 0) {
      if (is_decay_inner) {
	emissionAngle = azimuthalAngle + PI / 2.; 
      }
      else {
	emissionAngle = PI / 2. - azimuthalAngle;
      }
    }
    else {
      if (is_decay_inner) {
	emissionAngle = azimuthalAngle - PI / 2.;
      }
      else {
	emissionAngle = - (azimuthalAngle + PI / 2.);
      }
    }
    //std::cout << "Emission Angle = " << emissionAngle << std::endl;
    /*
    if (TMath::Sin(polarAngle) > 0) {
      emissionAngle = (azimuthalAngle < PI/2.) ? PI/2. - azimuthalAngle :  (PI/2. - azimuthalAngle);
    }
    else {
      emissionAngle = (azimuthalAngle < PI/2.) ? PI/2. + azimuthalAngle : (- 3 * PI/2. + azimuthalAngle);
    }
    */
    // skip tracking on events with not many hits
    if (nhits < 3) {
      //std::cout << "Skip event " << iev << std::endl;
      // redefine emissionAngle being taken with respect to spin direction
      detector_acc_theta->Fill(false, trueMomentum, emissionAngle);
      detector_acc_phi->Fill(false, trueMomentum, polarAngle);
      continue;
    }
    in_acceptance++;
    detector_acc_theta->Fill(true, trueMomentum, emissionAngle);
    detector_acc_phi->Fill(true, trueMomentum, polarAngle);
    //std::cout << "Track " << iev << " mom = " << trueMomentum << std::endl;
    //std::cout << "Origin at " << origin->X() * 0.1 << " " << origin->Y() * 0.1 << " " << origin->Z() * 0.1 << " with exit angles: phi = " << polarAngle << " and theta = " << azimuthalAngle << std::endl;
    // Fill efficiency histogram with good events
    

    // Sort hits
    
    *planeID = sort_vector_id(*hitsCoordinates, *planeID);
    *hitsCoordinates = sort_vector_z(*hitsCoordinates);
    
    drawn = 1;

    // Drawing
    /*
    // Draw xy view
    canvas->cd(1);
    TH2F* h2_xy = new TH2F("h2_xy", "XY View", 100, -10., 10., 100, -10., 10.);
    h2_xy->SetStats(0);
    h2_xy->GetXaxis()->SetTitle("X (cm)");
    h2_xy->GetYaxis()->SetTitle("Y (cm)");
    h2_xy->Draw("same");
    
    TGraph* gr_xy = drawXYView_hits(*hitsCoordinates, true);
    //gr_xy->SetTitle(Form("%d xy", i));
    gr_xy->Draw("CPsame");
    
    canvas->Update();
    
    // Draw yz view
    canvas->cd(2);
    TH2F* h2_yz = new TH2F("h2_yz", "YZ View", 100, -11., 11., 100, -10., 10.);
    h2_yz->SetStats(0);
    h2_yz->GetXaxis()->SetTitle("Z (cm)");
    h2_yz->GetYaxis()->SetTitle("Y (cm)");
    h2_yz->Draw("same");
    
    
    TGraph* gr_yz = drawYZView_hits(*hitsCoordinates, true);
    //gr_yz->SetTitle();
    gr_yz->Draw("CPsame");
    canvas->Update();
    */
    // Fit tracks

    
    // true start values
    TVector3 pos_cm = {origin->X(), origin->Y(), origin->Z()};
    TVector3 pos = pos_cm;
    TVector3 mom(1.,0,0);
    mom.SetPhi(polarAngle);
    mom.SetTheta(azimuthalAngle);
    mom.SetMag(trueMomentum * 1e-3);

    // Initial covariance
    TMatrixDSym cov(6);
    for (int i=0; i<6; i++) {
      for (int j=0; j<6; j++) {
	cov(i, j) = 0.1 * 0.1;
      }
    }
    
    // init fitter
    genfit::AbsKalmanFitter* fitter = new genfit::KalmanFitterRefTrack();
  
    // FIT WITH PLANES
    
    // trackrep
    genfit::AbsTrackRep* rep = new genfit::RKTrackRep(pdg);
    genfit::MeasuredStateOnPlane stateRef(rep);
    rep->setPosMomCov(stateRef, pos, mom, cov);

    TVectorD seedState(6);
    TMatrixDSym seedCov(6);
    rep->get6DStateCov(stateRef, seedState, seedCov);
    // create track
    //genfit::Track fitTrack(rep, pos, mom);
    genfit::Track fitTrack(rep, seedState, seedCov);
    
    const int detId(3); // detector ID
    int planeId(0); // detector plane ID
    int hitId(0); // hit ID

    double detectorResolution(0.1); // resolution of planar detectors
    TMatrixDSym hitCov(2);
    hitCov.UnitMatrix();
    hitCov *= detectorResolution*detectorResolution;
    
    // Bad resolution on virtual planes
    double virtualResolution(1.);
    TMatrixDSym virtualhitCov(3);
    virtualhitCov.UnitMatrix();
    virtualhitCov *= virtualResolution * virtualResolution;
    
    TVector3 previous_hit;
    int cylinderID = 100;
    // Loop over hit
    for (int i=0; i<nhits; i++) {
      TVectorD hitCoords(2);
      std::vector<double> AbsCoords = hitsCoordinates->at(i);
      double x_i = AbsCoords.at(0);
      double y_i = AbsCoords.at(1);
      double z_i = AbsCoords.at(2);

      
      TVector3 this_hit = {x_i, y_i, z_i};
      /*
      // Add virtual plane measurements on track when points are very far
      if (i != 0 &&
	  (this_hit.Phi() - previous_hit.Phi()) > 0.15)
	{
	  TVector3 mid(pos_cm.X(), pos_cm.Y(), (this_hit.Z() + previous_hit.Z())/2.);
	  TVector3 dir = (mid - previous_hit);
	  double len = dir.Mag();
	  dir = dir.Unit();
	  for (int j=0; j < 1; j++) {
	    TVectorD virtualhitCoords(3);
	    virtualhitCoords[0] = mid.X();//previous_hit.X() + dir.X() * (j + 1) * 3.;
	    virtualhitCoords[1] = mid.Y();//previous_hit.Y() + dir.Y() * (j + 1) * 3.;
	    virtualhitCoords[2] = mid.Z(); //previous_hit.Z() + dir.Z() * (j + 1) * 3.;
	    fitTrack.insertMeasurement(new genfit::SpacepointMeasurement(virtualhitCoords, virtualhitCov, 3, (int) (i * 1000 + j), nullptr));
	  }
	  
	  dir = (this_hit - mid);
	  len = dir.Mag();
	  dir = dir.Unit();
	  for (int j=0; j < (int)(len/3) - 1; j++) {
	    TVectorD virtualhitCoords(3);
	    virtualhitCoords[0] = mid.X() + dir.X() * (j + 1) * 3.;
	    virtualhitCoords[1] = mid.Y() + dir.Y() * (j + 1) * 3.;
	    virtualhitCoords[2] = mid.Z() + dir.Z() * (j + 1) * 3.;
	    fitTrack.insertMeasurement(new genfit::SpacepointMeasurement(virtualhitCoords, virtualhitCov, 3, (int) (i * 1000 + j), nullptr));
	  }
	  
	}
      
      // update previous hit
      previous_hit = this_hit;
      */
	
      // "Petals" planes
      // (u, v) measurements: u = radius (centered with respect to plane), v = z coordinate
      /*
      if (planeID->at(i) < 30) {
	//std::cout << "Hit on petal " << planeID->at(i) << " at " << x_i << " , " << y_i << " , " << z_i << std::endl;
	hitCoords[0] = + TMath::Hypot(x_i, y_i) - 3.0 - (8.5 - 4.5)/2 ; // cm
	hitCoords[1] = z_i; // cm
	// Plane center and orientation
	double Phi = TMath::ATan2(y_i, x_i);
	double X_C = (3.0 + (8.5 - 4.5)/2) * TMath::Cos(Phi);
	double Y_C = (3.0 + (8.5 - 4.5)/2) * TMath::Sin(Phi); 
	genfit::PlanarMeasurement* measurement = new genfit::PlanarMeasurement(hitCoords, hitCov, detId, i, nullptr);
	measurement->setPlane(genfit::SharedPlanePtr(new genfit::DetPlane(TVector3(X_C, Y_C, 0.), TVector3(TMath::Cos(Phi), TMath::Sin(Phi), 0), TVector3(0, 0, 1))), planeID->at(i));
	fitTrack.insertPoint(new genfit::TrackPoint(measurement, &fitTrack));
      }
      */
      
      // End Caps disks
      /*
      else {
	//std::cout << "Hit on endcap " << planeID->at(i) <<  " at " << x_i << " , " << y_i << " , " << z_i << std::endl;
	// u, v = x, y
	hitCoords[0] = x_i;
	hitCoords[1] = y_i;
	genfit::PlanarMeasurement* measurement = new genfit::PlanarMeasurement(hitCoords, hitCov, detId, i, nullptr);
	measurement->setPlane(genfit::SharedPlanePtr(new genfit::DetPlane(TVector3(0., 0., z_i), TVector3(1, 0, 0), TVector3(0, 1, 0))), planeID->at(i));
	fitTrack.insertPoint(new genfit::TrackPoint(measurement, &fitTrack));
      }
      */
      // End cap cylinders
      if (true) {
	// Plane center and orientation
	double Phi = TMath::ATan2(y_i, x_i);
	double X_C = x_i; // if not smearing is added, the origin of the virtual plane of the cylinder is set at the hit coordinate (fixed R and Phi)
	double Y_C = y_i;
	double Z_C = 0. ; // US and DS end cap
	// u, v = z, rphi is known since it is a cylinder
	hitCoords[0] = 0.0;
	hitCoords[1] = z_i - Z_C;
	 
	genfit::PlanarMeasurement* measurement = new genfit::PlanarMeasurement(hitCoords, hitCov, detId, i, nullptr);
	measurement->setPlane(genfit::SharedPlanePtr(new genfit::DetPlane(TVector3(X_C, Y_C, Z_C), TVector3(TMath::Sin(Phi), -TMath::Cos(Phi), 0), TVector3(0, 0, 1))), planeID->at(i));
	fitTrack.insertPoint(new genfit::TrackPoint(measurement, &fitTrack));
      }

    }

    // Add a virtual measurement at the origin of the track (can be provided by a first fit)
    /*
    TVectorD originCoords(3);
    originCoords[0] = pos_cm.X();
    originCoords[1] = pos_cm.Y();
    originCoords[2] = pos_cm.Z();
    fitTrack.insertMeasurement(new genfit::SpacepointMeasurement(originCoords, virtualhitCov, 3, nhits, nullptr));
    */
    
    // FIT WITH TRACK CAND
    /* 
    // track cand
    genfit::TrackCand myCand;

    myDetectorHitArray.Clear();
    
    // covariance
    double resolution = 0.001;
    TMatrixDSym cov(3);
    for (int i = 0; i < 3; ++i)
      cov(i,i) = resolution*resolution;

    // Loop over hits
    for (int i=0; i<nhits; i++) {
      //std::cout << "Creating hit " << i << endl;
      std::vector<double> hit = track.hitsCoordinates.at(i);
      TVector3 currentPos;
      currentPos.SetX(gRandom->Gaus(hit.at(0), resolution));
      currentPos.SetY(gRandom->Gaus(hit.at(1), resolution));
      currentPos.SetZ(gRandom->Gaus(hit.at(2), resolution));
      // Fill the TClonesArray and the TrackCand
      // In a real experiment, you detector code would deliver mySpacepointDetectorHits and fill the TClonesArray.
      // The patternRecognition would create the TrackCand.
      new(myDetectorHitArray[i]) genfit::mySpacepointDetectorHit(currentPos, cov);
      myCand.addHit(myDetId, i);
    }

    // Create track candidate (from Pattern Recognition)
    // smeared start values (would come from the pattern recognition)
    const bool smearPosMom = false;   // init the Reps with smeared pos and mom
    const double posSmear = 0.01;     // cm
    const double momSmear = 1. /180.*TMath::Pi();     // rad
    const double momMagSmear = 0.001;   // relative

    TVector3 posM(pos);
    TVector3 momM(mom);
    if (smearPosMom) {
      posM.SetX(gRandom->Gaus(posM.X(),posSmear));
      posM.SetY(gRandom->Gaus(posM.Y(),posSmear));
      posM.SetZ(gRandom->Gaus(posM.Z(),posSmear));

      momM.SetPhi(gRandom->Gaus(mom.Phi(),momSmear));
      momM.SetTheta(gRandom->Gaus(mom.Theta(),momSmear));
      momM.SetMag(gRandom->Gaus(mom.Mag(), momMagSmear*mom.Mag()));
    }

    // initial guess for cov
    TMatrixDSym covSeed(6);
    for (int i = 0; i < 3; ++i)
      covSeed(i,i) = resolution*resolution;
    for (int i = 3; i < 6; ++i)
      covSeed(i,i) = pow(resolution / nhits / sqrt(3), 2);

    // set start values and pdg to cand
    myCand.setPosMomSeedAndPdgCode(posM, momM, pdg);
    myCand.setCovSeed(covSeed);

    // create track
    genfit::Track fitTrack(myCand, factory, new genfit::RKTrackRep(pdg));
    */
    
    
    // do the fit
    
    try {
      fitter->processTrack(&fitTrack);
    }
    catch(genfit::Exception& e){
      detector_eff_theta->Fill(false, trueMomentum, emissionAngle);
      detector_eff_phi->Fill(false, trueMomentum, polarAngle);
      std::cerr << e.what();
      std::cerr << "Exception, next track" << std::endl;
      continue;
    }
    
    fitTrack.checkConsistency();
    
    // print fit result
    //fitTrack.getFittedState().Print();
    
    // add track to event display
    //display->addEvent(&fitTrack);
    // print fit result
    //fitTrack.getFittedState().Print();

    genfit::FitStatus *fitStatus = fitTrack.getFitStatus(rep);
    bool is_tracked = false;
    if (fitStatus->isFitConverged() && fitStatus->getNdf() > 0) {
      fitted++;
      // Fill efficiency histogram
      is_tracked = true;
      detector_eff_theta->Fill(true, trueMomentum, emissionAngle);
      detector_eff_phi->Fill(true, trueMomentum, polarAngle);
      // Project to orbit plane
      try {
	const genfit::MeasuredStateOnPlane &stFirst = fitTrack.getFittedState();
	//std::cout << "Fitted mom = " << rep->getMomMag(stFirst) * 1e3 << " Diff = " << rep->getMomMag(stFirst) * 1e3 - trueMomentum << std::endl;
	TVector3 pos_o;
	TVector3 mom_o;
	TMatrixDSym cov_o(6);
	stFirst.getPosMomCov(pos_o, mom_o, cov_o);
	genfit::MeasuredStateOnPlane state_orbit(rep);
	rep->setPosMomCov(state_orbit, pos_o, mom_o, cov_o);
	rep->extrapolateToPlane(state_orbit, genfit::SharedPlanePtr(new genfit::DetPlane(TVector3(0., 0., 0.), TVector3(1, 0, 0), TVector3(0, 1, 0))));
	state_orbit.getPosMomCov(pos_o, mom_o, cov_o);
	// Fill histograms for resolutions etc
	cov_o = CovarianceToCylinder(cov_o, mom_o);
	double mom_phi = (mom_o.Phi() > 0 ) ? mom_o.Phi() : mom_o.Phi() + 2. * PI;
	
	h_polarAngle_bias->Fill(polarAngle, polarAngle - mom_phi);
	h_azimuthalAngle_bias->Fill(azimuthalAngle, azimuthalAngle - mom_o.Theta());
	h_mom_bias->Fill(trueMomentum, trueMomentum - mom_o.Mag() * 1e3);
	h_XOrbit_bias->Fill(pos_cm.X(), pos_cm.X() - pos_o.X());
	h_YOrbit_bias->Fill(pos_cm.Y(), pos_cm.Y() - pos_o.Y());
	h_polarAngle_std->Fill(polarAngle, TMath::Sqrt(cov_o(3, 3)));
	h_azimuthalAngle_std->Fill(emissionAngle, TMath::Sqrt(cov_o(4, 4)));
	h_mom_std->Fill(trueMomentum, TMath::Sqrt(state_orbit.getMomVar())*1e3);
	h_XOrbit_std->Fill(pos_cm.X(), TMath::Sqrt(cov_o(0, 0)));
	h_YOrbit_std->Fill(pos_cm.Y(), TMath::Sqrt(cov_o(1, 1)));
	
	// Fill outTree
	for (int i=0; i<6; i++) {
	  parVar[i] = TMath::Sqrt(cov_o(i, i)); 
	}

	// Emission angle between radial and z axis
	double decay_angle_orbit_plane_fit = TMath::ATan2(pos_o.Y(), pos_o.X()); 
	TVector3 radial_axis_fit = {TMath::Cos(decay_angle_orbit_plane_fit), TMath::Sin(decay_angle_orbit_plane_fit), 0.};
	double r_fit = TMath::Hypot(pos_o.X(), pos_o.Y());
	int is_decay_fit_inner = (TMath::Hypot(pos_o.X() / r_fit - TMath::Sin(mom_o.Phi()), TMath::Cos(mom_o.Phi()) + pos_o.Y() / r_fit) < TMath::Hypot(pos_o.X() / r_fit - pos_o.Y() / r_fit, pos_o.X() / r_fit + pos_o.Y() / r_fit)) ? 1 : 0;
	TVector3 positron_vector_fit = {TMath::Sin(mom_o.Theta())*TMath::Cos(mom_phi), TMath::Sin(mom_o.Theta())*TMath::Sin(mom_phi), TMath::Cos(mom_o.Theta())};

	double emissionAngle_fit = 0.;
	
	if (TMath::Cos(mom_o.Theta()) >= 0) {
	  if (is_decay_fit_inner) {
	    emissionAngle_fit = mom_o.Theta() + PI / 2.; 
	  }
	  else {
	    emissionAngle_fit = PI / 2. - mom_o.Theta();
	  }
	}
	else {
	  if (is_decay_fit_inner) {
	    emissionAngle_fit = mom_o.Theta() - PI / 2.;
	  }
	  else {
	    emissionAngle_fit = - (mom_o.Theta() + PI / 2.);
	  }
	}
	

	parFit[0] = pos_o.X();
	parFit[1] = pos_o.Y();
	parFit[2] = 0.;
	parFit[3] = mom_o.Mag() * 1e3;
	parFit[4] = emissionAngle_fit;
	parFit[5] = mom_phi;
	parTrue[0] = pos_cm.X();
	parTrue[1] = pos_cm.Y();
	parTrue[2] = 0.;
	parTrue[3] = trueMomentum;
	parTrue[4] = emissionAngle;
	parTrue[5] = spinAngle;
	chisquare = fitStatus->getChi2();
	ngoodhits = (int)(nhits - fitStatus->getNFailedPoints());

	outTree->Fill();
	
	//std::cout << "Residual position in orbit plane " << (pos_o - pos_cm).X() << " cm, " << (pos_o - pos_cm).Y() << " cm, " << (pos_o - pos_cm).Z() << " cm " << std::endl;
      }
      catch(genfit::Exception& e){
	detector_eff_theta->Fill(false, trueMomentum, emissionAngle);
	detector_eff_phi->Fill(false, trueMomentum, polarAngle);
	std::cerr << e.what();
	std::cerr << "Exception, next track" << std::endl;
	continue;
      }
    }
    else {
      detector_eff_theta->Fill(false, trueMomentum, emissionAngle);
      detector_eff_phi->Fill(false, trueMomentum, polarAngle);
    }
    
    delete fitter;
    
  }

  
  // Draw muon trajectory
  /*
  TEllipse *muon_trajectory = new TEllipse(0., 0., 3., 3.);
  muon_trajectory->SetFillStyle(0);
  muon_trajectory->SetLineWidth(3);
  muon_trajectory->SetLineColor(kRed);
  canvas->cd(1);
  muon_trajectory->Draw("same");
  
  TLine *l = new TLine(0., -3., 0., 3.);
  l->SetLineColor(kRed);
  l->SetLineWidth(3);
  canvas->cd(2);
  l->Draw("same");
  
  canvas->Update();
  
  // Draw detectors in xy plane
  for (int i=0; i<36; i++) {
    double phi = 2 * i * PI / 36.;
    double xmin = 4.2 * TMath::Cos(phi);
    double xmax = 9.7 * TMath::Cos(phi);
    double ymin = 4.2 * TMath::Sin(phi);
    double ymax = 9.7 * TMath::Sin(phi);
    TLine* l = new TLine(xmin, ymin, xmax, ymax);
    l->SetLineColor(kGray);
    l->SetLineWidth(2);
    canvas->cd(1);
    l->Draw("same");
    
    TBox *b = new TBox(-5, ymin, 5, ymax);
    b->SetFillStyle(0);
    b->SetLineColorAlpha(kGray, 0.5);
    b->SetLineWidth(2);
    canvas->cd(2);
    b->Draw("same");
    
  }
    
  // Draw disks
  TEllipse* disk_out = new TEllipse(0., 0., 4., 4.);
  disk_out->SetLineWidth(2);
  disk_out->SetFillStyle(0);
  disk_out->SetLineColor(kOrange-7);
  
  TEllipse* disk_in = new TEllipse(0., 0., 2., 2.);
  disk_in->SetLineWidth(2);
  disk_in->SetFillStyle(0);
  disk_in->SetLineColor(kOrange-7);
  
  canvas->cd(1);
  disk_out->Draw("same");
  disk_in->Draw("same");
  
  // Look for hits on disks located at z < 0 and z > 0 between |zmin| and |zmax|
  double ndisks = 2.;
  double zmin = 6.; // cm
  double zstep = 1.0; // cm
  double disk_inner_radius = 2.0; // cm
  double disk_outer_radius = 4.0; // cm
  double zmax = zmin + (ndisks/2. - 1.) * zstep;

  
  for (int i=0; i<ndisks; i++) {
    // calculate hit position
    double zpos = (i < ndisks/2.) ? - zmax + zstep * i : zmin + zstep * (i - ndisks/2);
    TLine *l = new TLine(zpos, -4., zpos, 4.);
    l->SetLineWidth(2);
    l->SetLineColor(kOrange-7);
    canvas->cd(2);
    l->Draw("same");
  }
  
  
  canvas->Update();
  */
  
  //display->open();
  

  // Fill output file
  //TFile *outFile = TFile::Open("fittedTracksProperties_500umPrecision.root", "RECREATE");
  TFile *outFile = TFile::Open("fittedTrackPropertied_z20_fullGeo_7Cyl_Nopetals.root", "RECREATE");
  outTree->Write();
  outFile->Close();
  
  // Draw TEfficiency                                                                       
  
  std::cout << "Efficiency = " << fitted << " / " << in_acceptance << std::endl;
  
  TCanvas *cEff = new TCanvas();
  cEff->Divide(2, 2);
  
  gStyle->SetPalette(kTemperatureMap);
  
  TH2 *h_eff_th = detector_eff_theta->CreateHistogram();
  
  //h_eff_th->SaveAs("hEfficiencyTheta_3hits_chet_withEndCaps_30planes_500umPrecision.C");
  h_eff_th->SaveAs("hEfficiencyTheta_3hits_chet_z20_fullGeo_7Cyl_Nopetals.C");
  TH2 *h_eff_phi = detector_eff_phi->CreateHistogram();
  
  //h_eff_phi->SaveAs("hEfficiencyPhi_3hits_chet_withEndCaps_30planes_500umPrecision.C");
  h_eff_phi->SaveAs("hEfficiencyPhi_3hits_chet_z20_fullGeo_7Cyl_Nopetals.C");
  TH2 *h_acc_th = detector_acc_theta->CreateHistogram();
  
  //h_acc_th->SaveAs("hAcceptanceTheta_3hits_chet_withEndCaps_30planes_500umPrecision.C");
  h_acc_th->SaveAs("hAcceptanceTheta_3hits_chet_z20_fullGeo_7Cyl_Nopetals.C");
  TH2 *h_acc_phi = detector_acc_phi->CreateHistogram();
  
  //h_acc_phi->SaveAs("hAcceptancePhi_3hits_chet_withEndCaps_30planes_500umPrecision.C");
  h_acc_phi->SaveAs("hAcceptancePhi_3hits_chet_z20_fullGeo_7Cyl_Nopetals.C");

  cEff->cd(1);
  h_eff_th->Draw("COLZ0");
  cEff->cd(3);
  h_eff_phi->Draw("COLZ0");
  cEff->cd(2);
  h_acc_th->Draw("COLZ0");
  cEff->cd(4);
  h_acc_phi->Draw("COLZ0");

  cEff->Draw();
  
  
  std::vector<TH2F*> histos = {h_polarAngle_bias, h_azimuthalAngle_bias, h_mom_bias, h_XOrbit_bias, h_YOrbit_bias, h_polarAngle_std, h_azimuthalAngle_std, h_mom_std, h_XOrbit_std, h_YOrbit_std};

  for (auto h : histos) {
    TCanvas *c = new TCanvas();
    TProfile *p = h->ProfileX();
    p->SetMarkerStyle(20);
    p->SetMarkerSize(0.7);
    p->SetMarkerColor(kBlack);
    p->SetLineColor(kBlack);
    h->Draw("colz"); p->Draw("p0same");
    c->SaveAs(Form("%s_z20_fullGeo_7Cyl_Nopetals.C",h->GetName()));
  }
  
  TFile *resolution_file = TFile::Open("resolution_file_z20_fullGeo_7Cyl_Nopetals.root", "RECREATE");

  for (auto h : histos) {
    h->Write();
  }
  
  h_polarAngle_bias->Write();
  h_polarAngle_std->Write();
  h_azimuthalAngle_bias->Write();
  h_azimuthalAngle_std->Write();
  h_mom_bias->Write();
  h_mom_std->Write();
  h_XOrbit_bias->Write();
  h_XOrbit_std->Write();
  h_YOrbit_bias->Write();
  h_YOrbit_std->Write();

  resolution_file->Close();
  
  
}
