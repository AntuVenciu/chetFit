#ifndef FITTERALGORITHMS_HH
#define FITTERALGORITHMS_HH

#include <iostream>
#include <unistd.h>
#include <vector>
#include <cmath>
#include <numeric>
#include <random>

#include <TRandom.h>
#include <TMath.h>
#include <TChain.h>
#include <TEveManager.h>
#include <TGeoManager.h>
#include <TGeoMaterialInterface.h>
#include <TVector3.h>
#include <TDatabasePDG.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TEllipse.h>
#include <TBox.h>
#include <TH1I.h>
#include <TEfficiency.h>
#include <TMatrixD.h>

#include <ConstField.h>
#include <Exception.h>
#include <FieldManager.h>
#include <KalmanFitterRefTrack.h>
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


#include "auxiliaryalgorithms.hh"
#include "options.hh"


namespace FITALG
{
    void PlanarFitter(Options opts);
    void SpacepointFitter(Options opts);
};



struct CHeTResolutions
{
    CHeTResolutions(Double_t correlationPhiZ) : corrPhiZ(correlationPhiZ) {}

    // Resolutions (cm)
    Double_t sigmaR = 0.1 / sqrt(12);
    inline Double_t sigmaPhi(Int_t cylID = 2) const { return 0.1 / Radii[cylID]; };
    Double_t sigmaZ = 0.15;
    Double_t covRPhi = 0.;
    Double_t covRZ = 0.;
    inline Double_t covPhiZ(Int_t cylID = 2) const { return corrPhiZ * sigmaPhi(cylID) * sigmaZ; }
    
    // Correlation
    Double_t corrPhiZ;

    // Radii and transformation matrix (will be moved in a globals namespace)
    const Float_t Radii[7] = {2.1, 2.4, 3.7, 4.5, 6.5, 7.5, 8.5};

    // Matrix
    TMatrixDSym GetMatrixCylindrical(Int_t cylID) const
    {
        TMatrixDSym C_RphiZ(3);

        C_RphiZ(0,0) = sigmaR*sigmaR;
        C_RphiZ(0,1) = covRPhi;
        C_RphiZ(0,2) = covRZ;
        C_RphiZ(1,0) = covRPhi;
        C_RphiZ(1,1) = sigmaPhi(cylID)*sigmaPhi(cylID);
        C_RphiZ(1,2) = covPhiZ(cylID);
        C_RphiZ(2,0) = covRZ;
        C_RphiZ(2,1) = covPhiZ(cylID);
        C_RphiZ(2,2) = sigmaZ*sigmaZ;
        
        //C_RphiZ.Print();
        return C_RphiZ;
    };


    TMatrixDSym GetMatrixCartesian(Int_t cylID, Double_t phi) const
    {
        const Double_t R = Radii[cylID];
        const TMatrixDSym C_RphiZ = GetMatrixCylindrical(cylID);
        
        // Define the Jacobian matrix J
        TMatrixD J(3,3);
        
        J(0,0) = cos(phi);   J(0,1) = -R * sin(phi); J(0,2) = 0;
        J(1,0) = sin(phi);   J(1,1) = R * cos(phi);  J(1,2) = 0;
        J(2,0) = 0;          J(2,1) = 0;             J(2,2) = 1;

        // Compute transformed covariance: C_xyz = J * C_RphiZ * J^T
        TMatrixDSym C_xyz(3);
        TMatrixD temp = J * C_RphiZ;
        temp = temp * J.T();

        // Copy the symmetric part into TMatrixDSym
        for(Int_t i = 0; i < 3; i++)
        {
            for(Int_t j = 0; j < 3; j++)
            { // Fill only lower triangle
                C_xyz(i,j) = temp(i,j);
            }
        }

        //C_xyz.Print();
        return C_xyz;
    }
};


#endif  // FITTERALGORITHMS_HH