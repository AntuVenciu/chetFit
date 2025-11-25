#ifndef GLOBALSETTINGS_HH
#define GLOBALSETTINGS_HH

#include <cmath>

#include <TROOT.h>
#include <TMath.h>
#include <TVector3.h>
#include <TMatrixD.h>
#include <TMatrixDSym.h>


namespace muEDM
{
    extern Double_t B;  // kGaus // 2.2 T
};



namespace CHeT
{
  // CHeT detector global values
    // Structure [cm]
    constexpr Int_t nCylinders = 6;
    constexpr Float_t Radii[nCylinders] = {1.70, 2.10, 3.70, 6.50, 6.90, 7.30};
    constexpr Float_t Length = 30;
        // NFibers ?
        // Stereo angles [rad]
    constexpr Double_t inStereoAngle[nCylinders] = {0.3365, 0.4091, 0.6553, 0.7854, 0.7854, 0.7854};
    constexpr Double_t outStereoAngle[nCylinders] = {0.3476, 0.4196, 0.6632, 0.7854, 0.7854, 0.7854};
    inline Double_t stereoAngle(Int_t cylID) { return (inStereoAngle[cylID] + outStereoAngle[cylID]) / 2.; };

    constexpr Double_t fiberWidth = 0.05;
    constexpr Int_t nFibersPerSiPM = 4;

    // Resolutions [cm]
    struct Resolutions
    {
        Double_t sigmaPitch = (nFibersPerSiPM*fiberWidth) / sqrt(12);
        Double_t sigmaR = (2.*fiberWidth) / sqrt(12);

        // Fitting tricks
        Double_t scaleCov = 1.;

        // Matrix
        // (I'm using r, Rphi, z)
        inline TMatrixDSym GetMatrixCylindrical(Int_t cylID) const
        {
            TMatrixDSym C_RphiZ(3);
            const Double_t k = 0.5*sigmaPitch*sigmaPitch;

            C_RphiZ(0,0) = sigmaR*sigmaR;
            C_RphiZ(0,1) = 0.;
            C_RphiZ(0,2) = 0.;
            C_RphiZ(1,0) = 0.;
            C_RphiZ(1,1) = k * pow(cos(stereoAngle(cylID)), -2.);
            C_RphiZ(1,2) = 0.;
            C_RphiZ(2,0) = 0.;
            C_RphiZ(2,1) = 0.;
            C_RphiZ(2,2) = k * pow(sin(stereoAngle(cylID)), -2.);

            //C_RphiZ.Print();
            return scaleCov*scaleCov*C_RphiZ;
        };

        inline TMatrixDSym GetMatrixCartesian(Int_t cylID, Double_t phi) const
        {
            TMatrixDSym C_xyz = GetMatrixCylindrical(cylID);

            //std::cout << ">>> C_RphiZ = " << std::endl;
            //C_xyz.Print();

            // Define the Jacobian matrix J
            TMatrixD J(3,3);

            J(0,0) = cos(phi);   J(0,1) = -sin(phi); J(0,2) = 0;
            J(1,0) = sin(phi);   J(1,1) = cos(phi);  J(1,2) = 0;
            J(2,0) = 0;          J(2,1) = 0;             J(2,2) = 1;

            // Compute transformed covariance: C_xyz = J * C_RphiZ * J^T
            C_xyz.Similarity(J);

            //std::cout << ">>> C_xyz = " << std::endl;
            //C_xyz.Print();
            return C_xyz;
        };
    };
};



namespace ANS
{
  // Analysis tools global values
    // Detector response
    constexpr Int_t smearNFibers = 4;

    // Pattern recognition response
    constexpr Double_t sigmaSeedPos = 0.1;    // cm
    constexpr Double_t sigmaSeedMomMag = 0.1;    // relative
    constexpr Double_t sigmaSeedMomDir = 3. /180.*TMath::Pi();    // rad

    TMatrixDSym CovFromCardinalToCylindricalMom(TMatrixDSym cov, TVector3 mom);
    TMatrixD ComputeHelixJacobian(
        const TVector3& pivot,
        Double_t xC, Double_t yC, Double_t R,
        Double_t z0, Double_t tanLambda,
        Double_t B = 2.2
    );

};


#endif  // GLOBALSETTINGS_HH 