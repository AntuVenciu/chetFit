#ifndef FITTERALGORITHMS_HH
#define FITTERALGORITHMS_HH

#include <iostream>
#include <unistd.h>
#include <vector>
#include <cmath>
#include <numeric>
#include <random>

#include <TROOT.h>
#include <TRandom.h>
#include <TMath.h>
#include <TChain.h>
#include <TApplication.h>
#include <TEveManager.h>
#include <TGeoManager.h>
#include <TGeoMaterialInterface.h>
#include <TVector3.h>
#include <TDatabasePDG.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TGraph2D.h>
#include <TEllipse.h>
#include <TBox.h>
#include <TH1I.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TEfficiency.h>
#include <TMatrixD.h>
#include <TArc.h>
#include <Fit/Fitter.h>
#include <Math/Functor.h>
#include <TF1.h>
#include <TDecompChol.h>
#include <TFitResult.h>
#include <TMatrixDSymEigen.h>
#include <TDecompSVD.h>
#include <ROOT/RVec.hxx>

#include <ConstField.h>
#include <Exception.h>
#include <FieldManager.h>
#include <KalmanFitter.h>
#include <KalmanFitterRefTrack.h>
#include <DAF.h>
#include <StateOnPlane.h>
#include <Track.h>
#include <TrackPoint.h>
#include <TrackCand.h>
#include <MaterialEffects.h>
#include <RKTrackRep.h>
#include <EventDisplay.h>
#include <AbsMeasurement.h>
#include <PlanarMeasurement.h>
#include <SpacepointMeasurement.h>
#include <MeasurementProducer.h>
#include <MeasurementFactory.h>
#include <mySpacepointDetectorHit.h>
#include <mySpacepointMeasurement.h>
#include <KalmanFitterInfo.h>

#include "globalsettings.hh"
#include "options.hh"
#include "trackdatamanager.hh"
#include "patternalgorithms.hh"
#include "auxiliaryalgorithms.hh"


namespace FITALG
{
    void PlanarFitter(Options opts);
    void SpacepointFitter(Options opts);
    void HelixFitter(Options opts);

    std::array<Double_t, 6> HelixPrefitter(const std::vector<std::vector<Double_t>>& hitsCoordinates, const std::vector<Int_t>& cylinders, Options opts);
};


#endif  // FITTERALGORITHMS_HH