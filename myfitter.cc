#include <iostream>
#include <unistd.h>
#include <vector>
#include <cmath>
#include <numeric>

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
#include <PlanarMeasurement.h>

using namespace std;

constexpr Int_t DEBUG_LVL = 0;



vector<vector<Double_t>> SortVectorZ(vector<vector<Double_t>> vectors)
{
    // Sort vectors according to z(third) value

    vector<Double_t> z_of_vectors;
    for(auto v : vectors)
        z_of_vectors.push_back(v.at(2));
  
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

    return copy_vector;
}



void DrawXYView_hits(TVector3* origin, vector<vector<Double_t>> hitsCoordinates, TCanvas *canvas)
{
    canvas->cd();
    canvas->DrawFrame(-9, -9, 9, 9);

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
    Float_t R[7] = {8.5, 7.5, 6.5, 4.5, 3.7, 2.4, 2.1};
    for(auto i = 0; i < 7; i++)
    {
        ell[i] = new TEllipse(0, 0, R[i]);
        ell[i]->Draw("same");
    }

    gr->Draw("PC same");
    ogr->Draw("P same");
    canvas->Update();
}



void DrawYZView_hits(vector<vector<Double_t>> hitsCoordinates, TCanvas *canvas)
{
    canvas->cd();
    canvas->DrawFrame(-40, -9, 40, 9);

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

    TBox *box[7];
    Float_t L = 30;
    Float_t R[7] = {8.5, 7.5, 6.5, 4.5, 3.7, 2.4, 2.1};
    for(auto i = 0; i < 7; i++)
    {
        box[i] = new TBox(-L, -R[i], L, R[i]);
        box[i]->SetFillStyle(0);
        box[i]->SetLineColor(kBlack);
        box[i]->Draw("same");
    }

    gr->Draw("PC same");

    canvas->Update();
}



TGraph2D* DrawXYZView_hits(vector<vector<Double_t>> hitsCoordinates)
{
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
    gr->SetMarkerStyle(20);
    gr->SetLineColor(kBlue);
    return gr;
}





Int_t main(Int_t argc, char** argv)
{
    // Set event number
    Bool_t processAll = true;
    Bool_t saveMode = false;
    Int_t event = -1;
    Int_t eventMax = -1;
    Int_t opt;

    while((opt = getopt(argc, argv, "e:M:s:")) != -1)
    {
        switch(opt)
        {
            case 'e':
                processAll = false;
                event = stoi(optarg);
                break;
            case 'M':
                eventMax = stoi(optarg);
                break;
            case 's':
                saveMode = true;
                break;
            case '?':
                cerr << "Options:\n \t-e <_eventID_>\n \t-M <_MaxEvent_>\n";
                return 1;
        }
    }

    // Load events
    TChain *tracksChain = new TChain("HelixTrackTree");
    for(Int_t i=0; i<4; i++)
        tracksChain->Add(Form("../chet_sim_z20_FullGeo_7Cyl_Nopetals_%d.root", i));
    Int_t nEvents = tracksChain->GetEntries(); 
    cout << "\n>>> There are " << nEvents << " events\n" << endl;

    if(eventMax != -1)
        nEvents = eventMax;

    Double_t trueMomentum;
    Double_t polarAngle;
    Double_t azimuthalAngle;
    Double_t spinAngle;
    Double_t emissionAngle;
    TVector3* fOrigin = 0;
    vector<vector<Double_t>>* hitsCoordinates = 0;
    vector<vector<Double_t>>* trackCoordinates = 0;
    vector<Int_t>* planeID = 0;

    tracksChain->SetBranchAddress("trueMomentum", &trueMomentum);
    tracksChain->SetBranchAddress("polarAngle", &polarAngle);
    tracksChain->SetBranchAddress("azimuthalAngle", &azimuthalAngle);
    tracksChain->SetBranchAddress("spinAngle", &spinAngle);
    tracksChain->SetBranchAddress("origin", &fOrigin);
    tracksChain->SetBranchAddress("hitsCoordinates", &hitsCoordinates);
    tracksChain->SetBranchAddress("trackCoordinates", &trackCoordinates);
    tracksChain->SetBranchAddress("planeID", &planeID);
    



    // Init geometry and magnetic field
    new TGeoManager("DetectorGeometry", "CHET geometry");
    TGeoManager::Import("../detectorGeometry_z20_fullGeo_7Cyl_Nopetals.root");
    genfit::MaterialEffects::getInstance()->init(new genfit::TGeoMaterialInterface());
    Double_t B = 22.0; // kGaus // 2.2 T
    genfit::FieldManager::getInstance()->init(new genfit::ConstField(0., 0., B));

    // PID for positron
    const Int_t pdg = -11;

    // Init event display
    genfit::EventDisplay* display = genfit::EventDisplay::getInstance();

    // Init fitter (maxIterations, deltaPVal) (Possible values = 20, 1.E-3)
    genfit::AbsKalmanFitter* fitter = new genfit::KalmanFitterRefTrack();
    fitter->setDebugLvl(DEBUG_LVL);

    // Create Track
    genfit::Track* fitTrack = nullptr;
    
    // Get a specific event for single view mode
    gRandom->SetSeed(0);
    if(event == -1)
        event = gRandom->Integer(nEvents);

    // Event loop
    Int_t inAcceptance = 0; Int_t inEfficiency = 0;
    for(Int_t ev = 0; ev < nEvents; ev++)
    {
        if(!processAll)
            if(ev != event)
                continue;

        // Clean up
        delete fitTrack; fitTrack = nullptr;

        // Get event
        tracksChain->GetEntry(ev);
        Int_t nHits = hitsCoordinates->size();

        if(!processAll)
        {
            cout << "\n>>> Event number = " << ev << endl;
            cout << ">>> NHits = " << nHits << endl;
        }

        // Check acceptance
        if(nHits < 3)
        {
            if(!processAll)
                cout << ">>> Track is not in acceptance!" << endl;

            continue;
        }
        inAcceptance++;

        // Draw hits
        if(!processAll)
        {
            *hitsCoordinates = SortVectorZ(*hitsCoordinates);
    
            TCanvas *canvHitsXY = new TCanvas("canvHitsXY", "canvHitsXY", 700, 700);
            TCanvas *canvHitsYZ = new TCanvas("canvHitsYZ");
            TCanvas *canvHitsXYZ = new TCanvas("canvHitsXYZ");
    
            DrawXYView_hits(fOrigin, *hitsCoordinates, canvHitsXY);
            DrawYZView_hits(*hitsCoordinates, canvHitsYZ);
    
            canvHitsXYZ->cd();
            auto grXYZ = DrawXYZView_hits(*hitsCoordinates);
            grXYZ->Draw("P LINE");
            canvHitsXYZ->Update();
        }

        // Start values for the fit
        TVector3 pos = {fOrigin->X()*1E-1, fOrigin->Y()*1E-1, fOrigin->Z()*1E-1};
        TVector3 mom = {1, 0, 0};
        mom.SetPhi(polarAngle);
        mom.SetTheta(azimuthalAngle);
        mom.SetMag(trueMomentum * 1E-3);

        if(!processAll)
        {
            cout << Form(">>> Vertex position = (%f, %f, %f) cm", pos[0], pos[1], pos[2]) << endl;
            cout << Form(">>> Vertex momentum = (%f, %f, %f) MeV", mom[0]*1E3, mom[1]*1E3, mom[2]*1E3) << endl;
        }

        // Track Rep
        genfit::RKTrackRep *rep = new genfit::RKTrackRep(pdg);

        // Create Track
        fitTrack = new genfit::Track(rep, pos, mom);
    
        // IDs for detector, planes and hits
        const Int_t detId = 0;
        Int_t planeId = 0;
        Int_t hitId = 0;

        // Resolution of planar detectors
        const Double_t detectorResolution(0.1);
        TMatrixDSym hitCov(2);
        hitCov.UnitMatrix();
        hitCov *= detectorResolution*detectorResolution;

        // Add hits to track with coordinates
        for(Int_t i=0; i<nHits; i++)
        {
            TVectorD hitCoords(2);

            vector<Double_t> AbsCoords = hitsCoordinates->at(i);
            Double_t x_i = AbsCoords.at(0);
            Double_t y_i = AbsCoords.at(1);
            Double_t z_i = AbsCoords.at(2);

            Double_t phi = TMath::ATan2(y_i, x_i);
            Double_t X_C = x_i;
            Double_t Y_C = y_i;
            Double_t Z_C = 0.;

            // u, v = z, rphi is known since it is a cylinder
            hitCoords[0] = 0.;
	        hitCoords[1] = z_i - Z_C;
        
            genfit::PlanarMeasurement* measurement = new genfit::PlanarMeasurement(hitCoords, hitCov, detId, ++hitId, nullptr);
            measurement->setPlane(genfit::SharedPlanePtr(new genfit::DetPlane(TVector3(X_C, Y_C, Z_C), TVector3(-TMath::Sin(phi), TMath::Cos(phi), 0), TVector3(0, 0, 1))), ++planeId);

            fitTrack->insertPoint(new genfit::TrackPoint(measurement, fitTrack));
        }
        
        // Check
        fitTrack->checkConsistency();
    
        // Do the fit
        try
        {
            fitter->processTrack(fitTrack);
        }
        catch(genfit::Exception& e)
        {
            cerr << e.what();
            cerr << "Exception, next track" << endl;
            continue;
        }

        // Fit result
        if(fitTrack->getFitStatus(rep)->isFitConverged())
            inEfficiency++;

        if(!processAll)
        {
            cout << "\n\n>>> Did FIT converge? " << (fitTrack->getFitStatus(rep)->isFitConverged() ? "Yes" : "No") << "\n\n" << endl;
            //fitTrack->getFittedState().Print();
        }
    
        // Check
        fitTrack->checkConsistency();

        // Last steps
        display->addEvent(fitTrack);

        cout << "\r>>> Processed event number " << ev << flush;
    }
    cout << endl;

    // Print results
    if(processAll)
    {
        cout << ">>> Acceptance = " << (Float_t) (inAcceptance * 100) / nEvents << endl;
        cout << ">>> Efficiency = " << (Float_t) (inEfficiency * 100) / inAcceptance << endl;
    }

    // Delete fitter
    delete fitter;
    
    // Open event display
    display->setOptions("ABDEFGHMPT");
    display->open();

    // Finally
    return 0;
}