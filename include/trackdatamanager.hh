#ifndef TRACKDATAMANAGER_HH
#define TRACKDATAMANAGER_HH

#include <iostream>
#include <vector>

#include <TChain.h>
#include <TEfficiency.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TVector3.h>
#include <TMath.h>


class TrackDataManager
{
  public:
    TrackDataManager(Int_t maxEvents = -1);
    ~TrackDataManager();
    
    inline TChain* GetChain() const { return tracksChain; };
    inline Int_t GetNEvents() const { return nEvents; };

    // Data
    Double_t trueMomentum;
    Double_t polarAngle;
    Double_t azimuthalAngle;
    Double_t theThetaAngle;
    Double_t spinAngle;
    Double_t emissionAngle;
    TVector3* fOrigin;
    std::vector<std::vector<Double_t>>* hitsCoordinates;
    std::vector<std::vector<Double_t>>* trackCoordinates;
    std::vector<Int_t>* cylinderID;

    // Histos
    TEfficiency *accPhi, *accTheta, *effPhi, *effTheta;
    TH1I *histTurns, *histCylinders, *histFakeHits;
    TEfficiency *effTurns, *effCylinders;
    TProfile *histCylVMom, *histTurnsVMom;
    TH2D *graphMom, *graphTheta, *graphPhi;
    TH1D *histDiffX, *histDiffY, *histDiffZ;
    TH1D *histDiffMom, *histDiffTheta, *histDiffPhi;
    TH2D *hist2MomRes, *hist2ThetaRes, *hist2PhiRes;
    TProfile *profMomRes, *profThetaRes, *profPhiRes;
    
  private:
    TChain* tracksChain;
    Int_t nEvents;
};
    
#endif  // TRACKDATAMANAGER_HH
