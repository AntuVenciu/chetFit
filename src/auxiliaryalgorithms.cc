#include "auxiliaryalgorithms.hh"

using namespace std;



pair<vector<vector<Double_t>>, vector<Int_t>> AUXALG::SortVectorZ(vector<vector<Double_t>> vectors, vector<Int_t> cylinders)
{
    // Sort vectors according to z(third) value

    vector<Double_t> z_of_vectors;
    for(auto v : vectors)
        z_of_vectors.push_back(abs(v.at(2)));

    // initialize original index locations
    vector<size_t> idx(z_of_vectors.size());
    iota(idx.begin(), idx.end(), 0);

    // sort indexes based on comparing values in v
    // using stable_sort instead of sort
    // to avoid unnecessary index re-orderings
    // when v contains elements of equal values 
    stable_sort(idx.begin(), idx.end(),
        [&z_of_vectors](size_t i1, size_t i2) {return z_of_vectors[i1] < z_of_vectors[i2];});

    vector<vector<Double_t>> copy_vector;
    for(auto i : idx)
        copy_vector.push_back(vectors.at(i));

    vector<Int_t> copy_cylinders;
    for(auto i : idx)
        copy_cylinders.push_back(cylinders.at(i));

    return {copy_vector, copy_cylinders};
}



pair<vector<vector<Double_t>>, vector<Int_t>> AUXALG::ShuffleVectorZ(vector<vector<Double_t>> vectors, vector<Int_t> cylinders)
{
    // Shuffle the vectors and cylinders randomly
    
    // Create a random engine and a distribution
    random_device rd;
    mt19937 g(rd());

    // Combine vectors and cylinders into a single vector of pairs
    vector<pair<vector<Double_t>, Int_t>> combined;
    for (size_t i = 0; i < vectors.size(); ++i) {
        combined.push_back({vectors[i], cylinders[i]});
    }

    // Shuffle the combined vector randomly
    shuffle(combined.begin(), combined.end(), g);

    // Extract the shuffled vectors and cylinders back
    vector<vector<Double_t>> shuffled_vectors;
    vector<Int_t> shuffled_cylinders;
    for (const auto& p : combined) {
        shuffled_vectors.push_back(p.first);
        shuffled_cylinders.push_back(p.second);
    }

    return {shuffled_vectors, shuffled_cylinders};
}



void AUXALG::DrawXYView_hits(TVector3* origin, vector<vector<Double_t>> hitsCoordinates, TCanvas *canvas)
{
    canvas->cd();
    auto frame = canvas->DrawFrame(-9, -9, 9, 9);
    auto hframe = (TH1F*) gPad->GetPrimitive("hframe");
    hframe->SetLineWidth(0);
    frame->Draw();
    frame->SetTitle("X-Y View;X [cm];Y [cm]");

    Int_t nHits = hitsCoordinates.size();
    Double_t* x = new Double_t[nHits];
    Double_t* y = new Double_t[nHits];
    for(Int_t i = 0; i < nHits; i++)
    {
        x[i] = hitsCoordinates[i].at(0);
        y[i] = hitsCoordinates[i].at(1);
    }
    TGraph *gr = new TGraph(nHits, x, y);
    gr->SetMarkerStyle(20);
    gr->SetLineColor(kBlue);

    Double_t x0 = origin->X()*1E-1;
    Double_t y0 = origin->Y()*1E-1;
    TGraph *ogr = new TGraph(1, &x0, &y0);
    ogr->SetMarkerStyle(20);
    ogr->SetMarkerColor(kRed);

    TEllipse *ell[7];
    Float_t R[7] = {8.55, 7.55, 6.55, 3.9, 3.7, 2.1, 1.7};
    for(auto i = 0; i < 7; i++)
    {
        ell[i] = new TEllipse(0, 0, R[i]);
        ell[i]->Draw("same");
    }

    gr->Draw("PC same");
    ogr->Draw("P same");
    canvas->Update();
}



void AUXALG::DrawYZView_hits(TVector3* origin, vector<vector<Double_t>> hitsCoordinates, TCanvas *canvas)
{
    canvas->cd();
    auto frame = canvas->DrawFrame(-40, -9, 40, 9);
    auto hframe = (TH1F*) gPad->GetPrimitive("hframe");
    hframe->SetLineWidth(0);
    frame->Draw();
    frame->SetTitle("Z-Y View;Z [cm]; Y [cm]");

    Int_t nHits = hitsCoordinates.size();
    Double_t* y = new Double_t[nHits];
    Double_t* z = new Double_t[nHits];
    for(Int_t i = 0; i < nHits; i++)
    {
        y[i] = hitsCoordinates[i].at(1);
        z[i] = hitsCoordinates[i].at(2);
    }
    TGraph* gr = new TGraph(nHits, z, y);
    gr->SetMarkerStyle(20);
    gr->SetLineColor(kBlue);

    Double_t z0 = origin->Z()*1E-1;
    Double_t y0 = origin->Y()*1E-1;
    TGraph *ogr = new TGraph(1, &z0, &y0);
    ogr->SetMarkerStyle(20);
    ogr->SetMarkerColor(kRed);

    TBox *box[7];
    Float_t L = 30;
    Float_t R[7] = {8.55, 7.55, 6.55, 3.9, 3.7, 2.1, 1.7};
    for(auto i = 0; i < 7; i++)
    {
        box[i] = new TBox(-L, -R[i], L, R[i]);
        box[i]->SetFillStyle(0);
        box[i]->SetLineColor(kBlack);
        box[i]->Draw("same");
    }

    gr->Draw("PC same");
    ogr->Draw("P same");

    canvas->Update();
}



void AUXALG::DrawXYZView_hits(vector<vector<Double_t>> hitsCoordinates, TCanvas *canvas)
{
    canvas->cd();

    Int_t nHits = hitsCoordinates.size();
    Double_t* x = new Double_t[nHits];
    Double_t* y = new Double_t[nHits];
    Double_t* z = new Double_t[nHits];

    for(Int_t i = 0; i < nHits; ++i)
    {
        x[i] = hitsCoordinates[i].at(0);
        y[i] = hitsCoordinates[i].at(1);
        z[i] = hitsCoordinates[i].at(2);
    }

    TGraph2D *gr = new TGraph2D(nHits, z, x, y);
    gr->SetTitle("3D View; Z [cm]; X [cm]; Y [cm]");
    gr->SetMarkerStyle(20);
    gr->SetLineColor(kBlue);

    gr->Draw("P LINE");

    canvas->Update();
}



Int_t AUXALG::CountTurns(const vector<vector<Double_t>> hitsCoordinates) 
{
    if(hitsCoordinates.size() < 3) 
        return 0; // Servono almeno 3 punti per trovare un massimo o minimo

    Int_t nTurns = 0;
    Bool_t foundMax = false, foundMin = false;

    for(size_t i = 1; i < hitsCoordinates.size() - 1; i++)
    {
        Double_t yPrev = hitsCoordinates[i - 1][1];
        Double_t yCurr = hitsCoordinates[i][1];
        Double_t yNext = hitsCoordinates[i + 1][1];

        // Controlliamo se è un massimo locale
        if(yCurr > yPrev && yCurr > yNext)
        {
            foundMax = true;
        }
        // Controlliamo se è un minimo locale
        else if(yCurr < yPrev && yCurr < yNext)
        {
            foundMin = true;
        }

        // Se abbiamo sia un massimo che un minimo -> un giro completato
        if(foundMax && foundMin)
        {
            nTurns++;
            foundMax = false;
            foundMin = false;
        }
    }

    return nTurns;
}



pair<vector<vector<Double_t>>, vector<Int_t>> AUXALG::SelectTurn(Float_t turnID, const vector<vector<Double_t>>& hitsCoordinates, const vector<Int_t>& cylinders)
{
    vector<vector<Double_t>> turnHits;
    vector<Int_t> turnCylinders;

    if (hitsCoordinates.size() < 3)
        return {turnHits, turnCylinders}; // Troppi pochi punti per definire un giro

    Int_t totalTurns = CountTurns(hitsCoordinates);
    if (turnID > totalTurns)
        return {hitsCoordinates, cylinders}; // Se voglio più giri di quelli presenti, prendo tutta la traccia

    Double_t nHalfTurns = 0.0;
    turnHits.push_back(hitsCoordinates.front()); // Includi il primo punto
    turnCylinders.push_back(cylinders.front());

    for (size_t i = 1; i < hitsCoordinates.size() - 1; i++)
    {
        Double_t yPrev = hitsCoordinates[i - 1][1];
        Double_t yCurr = hitsCoordinates[i][1];
        Double_t yNext = hitsCoordinates[i + 1][1];

        // Identificazione di un massimo o minimo locale
        if ((yCurr > yPrev && yCurr > yNext) || (yCurr < yPrev && yCurr < yNext))
        {
            nHalfTurns += 1.0; // Ora conto direttamente i mezzi giri
        }

        // Se il numero di **giri completi** supera `turnID`, interrompo
        if (nHalfTurns / 2.0 >= turnID)
        {
            break;
        }

        turnHits.push_back(hitsCoordinates[i]);
        turnCylinders.push_back(cylinders[i]);
    }

    // Includi sempre l'ultimo punto del semigiro
    turnHits.push_back(hitsCoordinates[turnHits.size()]);
    turnCylinders.push_back(cylinders[turnCylinders.size()]);

    return {turnHits, turnCylinders};
}



vector<Int_t> AUXALG::SplitTurns(const vector<vector<Double_t>>& hitsCoordinates) 
{
    vector<Int_t> turnIndices;
    if(hitsCoordinates.size() < 3) 
        return turnIndices;
    
    turnIndices.push_back(0); // Il primo indice è sempre 0
    
    Bool_t foundMax = false, foundMin = false;
    
    for(size_t i = 1; i < hitsCoordinates.size() - 1; i++)
    {
        Double_t yPrev = hitsCoordinates[i - 1][1];
        Double_t yCurr = hitsCoordinates[i][1];
        Double_t yNext = hitsCoordinates[i + 1][1];

        // Controlliamo se è un massimo locale
        if(yCurr > yPrev && yCurr > yNext)
        {
            foundMax = true;
        }
        // Controlliamo se è un minimo locale
        else if(yCurr < yPrev && yCurr < yNext)
        {
            foundMin = true;
        }

        // Se troviamo un massimo e poi un minimo (o viceversa), aggiungiamo l'indice
        if(foundMax && foundMin)
        {
            turnIndices.push_back(i);
            foundMax = false;
            foundMin = false;
        }
    }
    
    return turnIndices;
}



Int_t AUXALG::CountCylinders(const vector<Int_t>& cylinders)
{
    set<Int_t> uniqueValues(cylinders.begin(), cylinders.end());
    return uniqueValues.size();
}



TMatrixDSym AUXALG::CovFromCardinalToCylindricalMom(TMatrixDSym cov, TVector3 mom)
{
    // Transform a covariance matrix in x,y,z, momx, momy, momz
    // into a covariance matrix in x, y, z, mom, theta, phi
    TMatrixDSym covCyl(cov);

    TMatrixD Jac(6, 6);
    Jac.Zero();

    Double_t p = mom.Mag();
    Double_t pt = TMath::Hypot(mom.X(), mom.Y());

    if(p == 0 || pt == 0)
        return covCyl;

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



tuple<vector<Double_t>, vector<Double_t>, vector<Double_t>> AUXALG::GetResults(genfit::Track *fitTrack, genfit::AbsTrackRep *rep, Double_t trueMom, TVector3 truePos, Double_t trueTheta, Double_t truePhi)
{   
    try
    {
        // Extrapolate to orbit
        const genfit::MeasuredStateOnPlane &stFirst = fitTrack->getFittedState();
        TVector3 posProj;
        TVector3 momProj;
        TMatrixDSym covProj;

        stFirst.getPosMomCov(posProj, momProj, covProj);
        genfit::MeasuredStateOnPlane stateOrbit(rep);
        rep->setPosMomCov(stateOrbit, posProj, momProj, covProj);
        rep->extrapolateToPlane(stateOrbit, genfit::SharedPlanePtr(new genfit::DetPlane(TVector3(0., 0., 0.), TVector3(1, 0, 0), TVector3(0, 1, 0))));
        stateOrbit.getPosMomCov(posProj, momProj, covProj);
        covProj = AUXALG::CovFromCardinalToCylindricalMom(covProj, momProj);
        
        // Compute angles
        Double_t x = posProj.X();
        Double_t y = posProj.Y();
        Double_t theta = momProj.Theta();
        Double_t momProjPhi = (momProj.Phi() > 0 ) ? momProj.Phi() : momProj.Phi() + TMath::TwoPi();

        Double_t pz = TMath::Cos(theta);
        Double_t pr = TMath::Sin(theta) * (x * TMath::Cos(momProjPhi) + y * TMath::Sin(momProjPhi)) / TMath::Sqrt(x*x + y*y);

        Double_t momProjTheta = TMath::ATan2(pz, pr); // angle in plane (e_r, z) in radiants, in (-pi, pi)

        // Pulls: 
            // Position
        Double_t dX = (posProj.X() - truePos.X()) / sqrt(covProj(0,0));
        Double_t dY = (posProj.Y() - truePos.Y()) / sqrt(covProj(1,1));
        Double_t dZ = (posProj.Z() - truePos.Z()) / sqrt(covProj(2,2));
            // Momentum
        Double_t dMom = (momProj.Mag()*1E3 - trueMom) / (sqrt(covProj(3,3))*1E3);
        Double_t dTheta = (momProjTheta - trueTheta) / sqrt(covProj(4,4));
        Double_t dPhi = TMath::ATan2(sin(momProjPhi - truePhi), cos(momProjPhi - truePhi)) / sqrt(covProj(5,5));
        //Double_t dMom = (momProj.Mag()*1E3 - trueMom);
        //Double_t dTheta = (momProjTheta - trueTheta);
        //Double_t dPhi = TMath::ATan2(sin(momProjPhi - truePhi), cos(momProjPhi - truePhi));

        return make_tuple(
            vector<Double_t>{posProj.X(), posProj.Y(), posProj.Z(), momProj.Mag()*1E3, momProjTheta, momProjPhi},
            vector<Double_t>{sqrt(covProj(0,0)), sqrt(covProj(1,1)), sqrt(covProj(2,2)), sqrt(covProj(3,3))*1E3, sqrt(covProj(4,4)), sqrt(covProj(5,5))},
            vector<Double_t>{dX, dY, dZ, dMom, dTheta, dPhi}
        );
    }
    catch(genfit::Exception& e)
    {
        cerr << "Exception, next track" << endl;
        cerr << e.what();
    }

    return make_tuple(vector<Double_t>(), vector<Double_t>(), vector<Double_t>());
}



vector<Double_t> AUXALG::SmearMeasurement(Int_t cylID, vector<Double_t> hitCoords)
{
    const Float_t Radii[7] = {1.7, 2.1, 3.7, 3.9, 6.55, 7.55, 8.55};

    if(cylID < 0 || cylID >= 7)
        throw out_of_range("Invalid cylinder ID");

    Double_t x = hitCoords.at(0);
    Double_t y = hitCoords.at(1);
    Double_t z = hitCoords.at(2);

    Double_t r_nominal = Radii[cylID];
    Double_t phi_nominal = TMath::ATan2(y, x);

    // --- Radial smearing ---
    Double_t r_min = r_nominal - 0.05;
    Double_t r_max = r_nominal + 0.05;
    Double_t r = sqrt(r_min*r_min + (r_max*r_max - r_min*r_min)*gRandom->Rndm());

    // --- Angular smearing ---
    Double_t dphi = 0.05 / r;
    Double_t phi = gRandom->Uniform(phi_nominal - dphi, phi_nominal + dphi);

    // --- Longitudinal smearing ---
    z += gRandom->Uniform(-0.05, 0.05);

    // --- Back to cartesian ---
    x = r*cos(phi);
    y = r*sin(phi);

    return {x, y, z};
}
