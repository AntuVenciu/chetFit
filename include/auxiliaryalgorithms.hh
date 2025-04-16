#ifndef AUXILIARYALGORITHMS_HH
#define AUXILIARYALGORITHMS_HH

#include <iostream>
#include <unistd.h>
#include <vector>
#include <cmath>
#include <numeric>
#include <random>

#include <TMath.h>
#include <TVector3.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TEllipse.h>
#include <TBox.h>
#include <TMatrixD.h>
#include <TH1F.h>
#include <TRandom.h>

#include <Track.h>
#include <RKTrackRep.h>


namespace AUXALG
{
    std::pair<std::vector<std::vector<Double_t>>, std::vector<Int_t>> SortVectorZ(std::vector<std::vector<Double_t>> vectors, std::vector<Int_t> cylinders);

    std::pair<std::vector<std::vector<Double_t>>, std::vector<Int_t>> ShuffleVectorZ(std::vector<std::vector<Double_t>> vectors, std::vector<Int_t> cylinders);

    void DrawXYView_hits(TVector3* origin, std::vector<std::vector<Double_t>> hitsCoordinates, TCanvas *canvas);

    void DrawYZView_hits(TVector3* origin, std::vector<std::vector<Double_t>> hitsCoordinates, TCanvas *canvas);

    void DrawXYZView_hits(std::vector<std::vector<Double_t>> hitsCoordinates, TCanvas *canvas);

    Int_t CountTurns(const std::vector<std::vector<Double_t>> hitsCoordinates);

    std::pair<std::vector<std::vector<Double_t>>, std::vector<Int_t>> SelectTurn(Float_t turnID, const std::vector<std::vector<Double_t>>& hitsCoordinates, const std::vector<Int_t>& cylinders);

    std::vector<Int_t> SplitTurns(const std::vector<std::vector<Double_t>>& hitsCoordinates);

    Int_t CountCylinders(const std::vector<Int_t>& cylinders);

    TMatrixDSym CovFromCardinalToCylindricalMom(TMatrixDSym cov, TVector3 mom);

    std::tuple<std::vector<Double_t>, std::vector<Double_t>, std::vector<Double_t>> GetResults(genfit::Track *fitTrack, genfit::AbsTrackRep *rep, Double_t trueMom, TVector3 truePos, Double_t trueTheta, Double_t truePhi);

    std::vector<Double_t> SmearMeasurement(Int_t cylID, std::vector<Double_t> hitCoords);
};


#endif  // AUXILIARYALGORITHMS_HH