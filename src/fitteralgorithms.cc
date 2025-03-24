#include "fitteralgorithms.hh"

using namespace std;
// Units are in cm

constexpr Int_t DEBUG_LVL = 0;
constexpr Double_t CORR_PHIZ = -0.01;


void FITALG::PlanarFitter(Options opts)
{
    // Load events
    TChain *tracksChain = new TChain("HelixTrackTree");
    for(Int_t i=0; i<4; i++)
        tracksChain->Add(Form("../chet_sim_z20_FullGeo_7Cyl_Nopetals_%d.root", i));
    Int_t nEvents = tracksChain->GetEntries(); 
    cout << "\n>>> There are " << nEvents << " events\n" << endl;

    if(opts.eventMax != -1)
        nEvents = opts.eventMax;

    Double_t trueMomentum;
    Double_t polarAngle;
    Double_t azimuthalAngle;
    Double_t spinAngle;
    Double_t emissionAngle;
    TVector3* fOrigin = 0;
    vector<vector<Double_t>>* hitsCoordinates = 0;
    vector<vector<Double_t>>* trackCoordinates = 0;
    vector<Int_t>* cylinderID = 0;

    tracksChain->SetBranchAddress("trueMomentum", &trueMomentum);
    tracksChain->SetBranchAddress("polarAngle", &polarAngle);
    tracksChain->SetBranchAddress("azimuthalAngle", &azimuthalAngle);
    tracksChain->SetBranchAddress("spinAngle", &spinAngle);
    tracksChain->SetBranchAddress("origin", &fOrigin);
    tracksChain->SetBranchAddress("hitsCoordinates", &hitsCoordinates);
    tracksChain->SetBranchAddress("trackCoordinates", &trackCoordinates);
    tracksChain->SetBranchAddress("planeID", &cylinderID);

    // Useful objects
    TH1I *histEvents = new TH1I("histEvents", "Efficiency per-Event;eventID;Efficiency", nEvents, 0, nEvents);
    TH1I *histTurns = new TH1I("histTurns", "Number of Turns;nTurns;Counts", 20, 0, 20);
    TEfficiency *effTurns = new TEfficiency("effTurns","nTurns Efficiency;nTurns;Efficiency", 10, 0, 10);



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
    display->reset();

    // Init fitter (maxIterations, deltaPVal) (Possible values = 20, 1.E-3)
    genfit::AbsKalmanFitter* fitter = new genfit::KalmanFitterRefTrack(20, 1.E-3);
    fitter->setDebugLvl(DEBUG_LVL);

    // Create Track
    genfit::Track* fitTrack = nullptr;
    
    // Get a specific event for single view mode
    gRandom->SetSeed(0);
    if(opts.event == -1)
        opts.event = gRandom->Integer(nEvents);

    // Event loop
    Int_t inAcceptance = 0; Int_t inEfficiency = 0;
    Int_t nTurns;
    for(Int_t ev = 0; ev < nEvents; ev++)
    {
        if(!opts.processAll)
            if(ev != opts.event)
                continue;

        // Clean up
        delete fitTrack; fitTrack = nullptr;

        // Get event
        tracksChain->GetEntry(ev);
        Int_t nHits = hitsCoordinates->size();

        if(!opts.processAll)
        {
            cout << "\n>>> Event number = " << ev << endl;
            cout << ">>> NHits = " << nHits << endl;
        }

        // Check acceptance
        if(nHits < 3)
        {
            if(!opts.processAll)
                cout << ">>> Track is not in acceptance!" << endl;

            continue;
        }
        inAcceptance++;

        // Sort and if in single event mode draw hits
            // will need a revision when origin point is not in Z = 0 anymore
        auto sortedHits = AUXALG::SortVectorZ(*hitsCoordinates, *cylinderID);
        *hitsCoordinates = sortedHits.first;
        *cylinderID = sortedHits.second;
        nTurns = AUXALG::CountTurns(*hitsCoordinates);

        if(opts.turnMode)
        {            
            if(!opts.processAll)
                cout << ">>> nTurns = " << nTurns << endl;

            auto turnHits = AUXALG::SelectTurn(opts.turnID, *hitsCoordinates, *cylinderID);
            *hitsCoordinates = turnHits.first;
            *cylinderID = turnHits.second;

            nHits = hitsCoordinates->size();

            // Check acceptance again
            if(nHits < 3)
            {
                if(!opts.processAll)
                    cout << ">>> Track is not in acceptance anymore!" << endl;

                inAcceptance--;
                continue;
            }
        }

        if(!opts.processAll)
        {
            TCanvas *canvHitsXY = new TCanvas("canvHitsXY", "canvHitsXY", 700, 700);
            TCanvas *canvHitsYZ = new TCanvas("canvHitsYZ");
            TCanvas *canvHitsXYZ = new TCanvas("canvHitsXYZ");

            AUXALG::DrawXYView_hits(fOrigin, *hitsCoordinates, canvHitsXY);
            AUXALG::DrawYZView_hits(fOrigin, *hitsCoordinates, canvHitsYZ);

            canvHitsXYZ->cd();
            auto grXYZ = AUXALG::DrawXYZView_hits(*hitsCoordinates);
            grXYZ->Draw("P LINE");
            canvHitsXYZ->Update();
        }

        // Start values for the fit
        TVector3 pos = {fOrigin->X()*1E-1, fOrigin->Y()*1E-1, fOrigin->Z()*1E-1};
        TVector3 mom = {1, 0, 0};
        mom.SetPhi(polarAngle);
        mom.SetTheta(azimuthalAngle);
        mom.SetMag(trueMomentum * 1E-3);

        if(!opts.processAll)
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
        const Double_t detectorResolution = 0.1;
        TMatrixDSym hitCov(2);
        hitCov.UnitMatrix();
        hitCov *= detectorResolution*detectorResolution;

        // Add hits to track with coordinates
        for(Int_t i = 0; i < nHits; i++)
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

            genfit::PlanarMeasurement* measurement = new genfit::PlanarMeasurement(hitCoords, hitCov, (*cylinderID)[i], ++hitId, nullptr);
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
        //fitTrack->Print();
        Bool_t isFitConverged = fitTrack->getFitStatus(rep)->isFitConverged();
        
        if(isFitConverged)
        {
            inEfficiency++;
            
            histEvents->Fill(ev);
        }

        if(!opts.processAll)
        {
            cout << "\n\n>>> Did FIT converge? " << (isFitConverged ? "Yes" : "No") << "\n\n" << endl;
            //fitTrack->getFittedState().Print();
        }

        // Store turns data
        histTurns->Fill(nTurns);
        effTurns->Fill(isFitConverged, nTurns);

        // Check
        fitTrack->checkConsistency();

        // Last steps
        display->addEvent(fitTrack);

        cout << "\r>>> Processed event number " << ev << flush;
    }
    cout << endl;

    // Print results
    if(opts.processAll)
    {
        cout << "\n---------------------------------------------------------" << endl;
        cout << ">>> Acceptance = " << (Float_t) (inAcceptance * 100) / nEvents << endl;
        cout << ">>> Efficiency = " << (Float_t) (inEfficiency * 100) / inAcceptance << endl;
        cout << "---------------------------------------------------------\n" << endl;
    }

    // Draw Graphs
    if(opts.processAll)
    {
        TCanvas *canvEvents = new TCanvas("canvEvents");
        canvEvents->cd();
        histEvents->SetLineColor(0);
        histEvents->SetFillColor(kBlack);
        histEvents->Draw();
    
        TCanvas *canvTurns = new TCanvas("canvTurns");
        canvTurns->Divide(2);
        canvTurns->cd(1);
        histTurns->Draw();
        canvTurns->cd(2);
        effTurns->Draw("AP");
    }

    // Delete fitter
    delete fitter;

    // Open event display
    display->setOptions("ABDEFGHMPT");
    display->open();

    // Finally
    return;
}



void FITALG::SpacepointFitter(Options opts)
{
    // Load events
    TChain *tracksChain = new TChain("HelixTrackTree");
    for(Int_t i=0; i<4; i++)
        tracksChain->Add(Form("../chet_sim_z20_FullGeo_7Cyl_Nopetals_%d.root", i));
    Int_t nEvents = tracksChain->GetEntries(); 
    cout << "\n>>> There are " << nEvents << " events\n" << endl;

    if(opts.eventMax != -1)
        nEvents = opts.eventMax;

    Double_t trueMomentum;
    Double_t polarAngle;
    Double_t azimuthalAngle;
    Double_t spinAngle;
    Double_t emissionAngle;
    TVector3* fOrigin = 0;
    vector<vector<Double_t>>* hitsCoordinates = 0;
    vector<vector<Double_t>>* trackCoordinates = 0;
    vector<Int_t>* cylinderID = 0;

    tracksChain->SetBranchAddress("trueMomentum", &trueMomentum);
    tracksChain->SetBranchAddress("polarAngle", &polarAngle);
    tracksChain->SetBranchAddress("azimuthalAngle", &azimuthalAngle);
    tracksChain->SetBranchAddress("spinAngle", &spinAngle);
    tracksChain->SetBranchAddress("origin", &fOrigin);
    tracksChain->SetBranchAddress("hitsCoordinates", &hitsCoordinates);
    tracksChain->SetBranchAddress("trackCoordinates", &trackCoordinates);
    tracksChain->SetBranchAddress("planeID", &cylinderID);

    // Useful objects
    TH1I *histEvents = new TH1I("histEvents", "Efficiency per-Event;eventID;Efficiency", nEvents, 0, nEvents);
    TH1I *histTurns = new TH1I("histTurns", "Number of Turns;nTurns;Counts", 20, 0, 20);
    TEfficiency *effTurns = new TEfficiency("effTurns","nTurns Efficiency;nTurns;Efficiency", 10, 0, 10);



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
    display->reset();

    // Init fitter (maxIterations, deltaPVal) (Possible values = 20, 1.E-3)
    genfit::AbsKalmanFitter* fitter = new genfit::KalmanFitterRefTrack(20, 1.E-3);
    fitter->setDebugLvl(DEBUG_LVL);

    // Create array of hits
    TClonesArray myDetectorHitArray("genfit::mySpacepointDetectorHit");

    // Init the factory
    Int_t detId = 0;
    genfit::MeasurementFactory<genfit::AbsMeasurement> factory;
    genfit::MeasurementProducer<genfit::mySpacepointDetectorHit, genfit::mySpacepointMeasurement> myProducer(&myDetectorHitArray);
    factory.addProducer(detId, &myProducer);

    // Create Track
    genfit::Track* fitTrack = nullptr;
    
    // Get a specific event for single view mode
    gRandom->SetSeed(0);
    if(opts.event == -1)
        opts.event = gRandom->Integer(nEvents);

    // Event loop
    Int_t inAcceptance = 0; Int_t inEfficiency = 0;
    Int_t nTurns;
    for(Int_t ev = 0; ev < nEvents; ev++)
    {
        if(!opts.processAll)
            if(ev != opts.event)
                continue;

        // Clean up
        delete fitTrack; fitTrack = nullptr;
        myDetectorHitArray.Clear();

        // Get event
        tracksChain->GetEntry(ev);
        Int_t nHits = hitsCoordinates->size();

        if(!opts.processAll)
        {
            cout << "\n>>> Event number = " << ev << endl;
            cout << ">>> NHits = " << nHits << endl;
        }

        // Check acceptance
        if(nHits < 3)
        {
            if(!opts.processAll)
                cout << ">>> Track is not in acceptance!" << endl;

            continue;
        }
        inAcceptance++;

        // Sort and if in single event mode draw hits
            // will need a revision when origin point is not in Z = 0 anymore
        auto sortedHits = AUXALG::SortVectorZ(*hitsCoordinates, *cylinderID);
        *hitsCoordinates = sortedHits.first;
        *cylinderID = sortedHits.second;
        nTurns = AUXALG::CountTurns(*hitsCoordinates);
        
        if(opts.turnMode)
        {            
            if(!opts.processAll)
                cout << ">>> nTurns = " << nTurns << endl;

            auto turnHits = AUXALG::SelectTurn(opts.turnID, *hitsCoordinates, *cylinderID);
            *hitsCoordinates = turnHits.first;
            *cylinderID = turnHits.second;

            nHits = hitsCoordinates->size();

            // Check acceptance again
            if(nHits < 3)
            {
                if(!opts.processAll)
                    cout << ">>> Track is not in acceptance anymore!" << endl;

                inAcceptance--;
                continue;
            }
        }

        if(!opts.processAll)
        {
            TCanvas *canvHitsXY = new TCanvas("canvHitsXY", "canvHitsXY", 700, 700);
            TCanvas *canvHitsYZ = new TCanvas("canvHitsYZ");
            TCanvas *canvHitsXYZ = new TCanvas("canvHitsXYZ");

            AUXALG::DrawXYView_hits(fOrigin, *hitsCoordinates, canvHitsXY);
            AUXALG::DrawYZView_hits(fOrigin, *hitsCoordinates, canvHitsYZ);

            canvHitsXYZ->cd();
            auto grXYZ = AUXALG::DrawXYZView_hits(*hitsCoordinates);
            grXYZ->Draw("P LINE");
            canvHitsXYZ->Update();
        }

        // Start values for the fit
        TVector3 pos = {fOrigin->X()*1E-1, fOrigin->Y()*1E-1, fOrigin->Z()*1E-1};
        TVector3 mom = {1, 0, 0};
        mom.SetPhi(polarAngle);
        mom.SetTheta(azimuthalAngle);
        mom.SetMag(trueMomentum * 1E-3);

        if(!opts.processAll)
        {
            cout << Form(">>> Vertex position = (%f, %f, %f) cm", pos[0], pos[1], pos[2]) << endl;
            cout << Form(">>> Vertex momentum = (%f, %f, %f) MeV", mom[0]*1E3, mom[1]*1E3, mom[2]*1E3) << endl;
        }

        // Track candidate
        genfit::TrackCand trackCand;

        // Resolution of detectors
        const Double_t detectorResolution = 0.1;
        const CHeTResolutions hitCov(CORR_PHIZ);
    
        //TMatrixDSym hitCov(3);
        //hitCov.UnitMatrix();
        //hitCov *= detectorResolution*detectorResolution;

        // Fill the candidate
        for(Int_t i = 0; i < nHits; i++)
        {
            TVector3 hitCoords;

            vector<Double_t> AbsCoords = hitsCoordinates->at(i);
            hitCoords[0] = AbsCoords.at(0);
            hitCoords[1] = AbsCoords.at(1);
            hitCoords[2] = AbsCoords.at(2);
        
            new(myDetectorHitArray[i]) genfit::mySpacepointDetectorHit(hitCoords, hitCov.GetMatrixCartesian((*cylinderID)[i], TMath::ATan2(hitCoords[1], hitCoords[0])));
            //new(myDetectorHitArray[i]) genfit::mySpacepointDetectorHit(hitCoords, hitCov);
            trackCand.addHit(detId, i);
        }

        // Smearing?

        // Initial guess for cov
        TMatrixDSym covSeed(6);
        for(Int_t i = 0; i < 3; i++)
            covSeed(i,i) = detectorResolution*detectorResolution;
        for(Int_t i = 3; i < 6; i++)
            covSeed(i,i) = pow(detectorResolution / nHits / sqrt(3), 2);

        // Set start values
        trackCand.setPosMomSeedAndPdgCode(pos, mom, pdg);
        trackCand.setCovSeed(covSeed);

        // Track rep
        genfit::AbsTrackRep *rep = new genfit::RKTrackRep(pdg);

        // Create track
        fitTrack = new genfit::Track(trackCand, factory, rep);
    
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
        //fitTrack->Print();
        Bool_t isFitConverged = fitTrack->getFitStatus(rep)->isFitConverged();
        
        if(isFitConverged)
        {
            inEfficiency++;
            
            histEvents->Fill(ev);
        }

        if(!opts.processAll)
        {
            cout << "\n\n>>> Did FIT converge? " << (isFitConverged ? "Yes" : "No") << "\n\n" << endl;
            //fitTrack->getFittedState().Print();
        }

        // Store turns data
        histTurns->Fill(nTurns);
        effTurns->Fill(isFitConverged, nTurns);

        // Check
        fitTrack->checkConsistency();

        // Last steps
        display->addEvent(fitTrack);

        cout << "\r>>> Processed event number " << ev << flush;
    }
    cout << endl;

    // Print results
    if(opts.processAll)
    {
        cout << "\n---------------------------------------------------------" << endl;
        cout << ">>> Acceptance = " << (Float_t) (inAcceptance * 100) / nEvents << endl;
        cout << ">>> Efficiency = " << (Float_t) (inEfficiency * 100) / inAcceptance << endl;
        cout << "---------------------------------------------------------\n" << endl;
    }

    // Draw Graphs
    if(opts.processAll)
    {
        TCanvas *canvEvents = new TCanvas("canvEvents");
        canvEvents->cd();
        histEvents->SetLineColor(0);
        histEvents->SetFillColor(kBlack);
        histEvents->Draw();
    
        TCanvas *canvTurns = new TCanvas("canvTurns");
        canvTurns->Divide(2);
        canvTurns->cd(1);
        histTurns->Draw();
        canvTurns->cd(2);
        effTurns->Draw("AP");
    }

    // Delete fitter
    delete fitter;

    // Open event display
    display->setOptions("ABDEFGHMPT");
    display->open();

    // Finally
    return;
}