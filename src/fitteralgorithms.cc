#include "fitteralgorithms.hh"

using namespace std;
using namespace ROOT;
// Units are in cm

constexpr Int_t DEBUG_LVL = 0;
constexpr Bool_t NORMALIZED_PULLS = true;


void FITALG::PlanarFitter(Options opts)
{
    // Load events
    TChain *tracksChain = new TChain("HelixTrackTree");
    for(Int_t i=0; i<2; i++)
        tracksChain->Add(Form("../chet_sim_dataset_%d.root", i));
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
    tracksChain->SetBranchAddress("polarAngle", &azimuthalAngle); // Note that now the convention on polar/azimuth is corrected!
    tracksChain->SetBranchAddress("azimuthalAngle", &polarAngle);
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
    TGeoManager::Import("../chet_sim_geometry.gdml");
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
        auto sortedHits = PTTALG::SortVectorZ(*hitsCoordinates, *cylinderID);
        *hitsCoordinates = sortedHits.first;
        *cylinderID = sortedHits.second;
        nTurns = PTTALG::CountTurns(*hitsCoordinates);

        if(opts.turnMode)
        {            
            if(!opts.processAll)
                cout << ">>> nTurns = " << nTurns << endl;

            auto turnHits = PTTALG::SelectTurn(opts.turnID, *hitsCoordinates, *cylinderID);
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
            AUXALG::DrawXYZView_hits(*hitsCoordinates, canvHitsXYZ);
        }

        // Start values for the fit
        TVector3 pos = {fOrigin->X()*1E-1, fOrigin->Y()*1E-1, fOrigin->Z()*1E-1};
        TVector3 mom = {1, 0, 0};
        mom.SetPhi(azimuthalAngle);
        mom.SetTheta(polarAngle);
        mom.SetMag(trueMomentum * 1E-3);

        if(!opts.processAll)
        {
            cout << Form(">>> Vertex position = (%f, %f, %f) cm", pos[0], pos[1], pos[2]) << endl;
            cout << Form(">>> Vertex momentum = (%f, %f, %f) MeV/c", mom[0]*1E3, mom[1]*1E3, mom[2]*1E3) << endl;
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
            measurement->setPlane(genfit::SharedPlanePtr(new genfit::DetPlane(TVector3(X_C, Y_C, Z_C), TVector3(-sin(phi), TMath::Cos(phi), 0), TVector3(0, 0, 1))), ++planeId);

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
    exit(0);
}



void FITALG::SpacepointFitter(Options opts)
{
    // Load events
    TrackDataManager data(opts.eventMax);

    Int_t nEvents = data.GetNEvents();
    Int_t processedEvents = nEvents;

    // Init geometry and magnetic field
    new TGeoManager("DetectorGeometry", "CHET geometry");
    TGeoManager::Import("../chet_sim_geometry.gdml");
    genfit::MaterialEffects::getInstance()->init(new genfit::TGeoMaterialInterface());
    Double_t B = 22.0; // kGaus // 2.2 T
    genfit::FieldManager::getInstance()->init(new genfit::ConstField(0., 0., B));

    // PID for positron
    const Int_t pdg = -11;

    // Init event display
    genfit::EventDisplay* display = genfit::EventDisplay::getInstance();
    display->reset();
    cout << endl;

    // Init fitter (maxIterations, deltaPVal) (Possible values = 20, 1.E-3)
    genfit::AbsKalmanFitter* fitter = new genfit::KalmanFitterRefTrack(20, 1.E-3);
    //genfit::AbsKalmanFitter* fitter = new genfit::KalmanFitter(20, 1.E-3);
    
    //genfit::AbsKalmanFitter* fitter = new genfit::DAF(true);

    // Set the annealing scheme
    //static_cast<genfit::DAF*>(fitter)->setAnnealingScheme(1000., 1., 10);
    //static_cast<genfit::DAF*>(fitter)->setConvergenceDeltaWeight(0.00001);
    //fitter->setMaxIterations(20);

    fitter->setDebugLvl(DEBUG_LVL);

    // Create array of hits
    TClonesArray chetHitArray("genfit::mySpacepointDetectorHit");

    // Init the factory
    Int_t detId = 0;
    genfit::MeasurementFactory<genfit::AbsMeasurement> chetFactory;
    genfit::MeasurementProducer<genfit::mySpacepointDetectorHit, genfit::mySpacepointMeasurement> cylProducer(&chetHitArray);
    chetFactory.addProducer(detId, &cylProducer);

    // Create Track
    genfit::Track* fitTrack = nullptr;

    // Get a specific event for single view mode
    gRandom->SetSeed(0);
    if(opts.event == -1)
        opts.event = gRandom->Integer(nEvents);

    // Event loop
    Int_t inAcceptance = 0; Int_t inEfficiency = 0;
    Int_t nTurns, nCylinders;
    for(Int_t ev = 0; ev < nEvents; ev++)
    {
        if(!opts.processAll)
            if(ev != opts.event)
                continue;

        // Clean up
        delete fitTrack; fitTrack = nullptr;
        chetHitArray.Clear();

        // Get event
        data.GetChain()->GetEntry(ev);
        Int_t nHits = data.hitsCoordinates->size();

        if(!opts.processAll)
        {
            cout << "\n>>> Event number = " << ev << endl;
            cout << ">>> NHits = " << nHits << endl;
        }

        // Get True MC data
        TVector3 pos = {data.fOrigin->X()*1E-1, data.fOrigin->Y()*1E-1, data.fOrigin->Z()*1E-1};
        TVector3 mom = {1, 0, 0};
        mom.SetPhi(data.azimuthalAngle);
        mom.SetTheta(data.polarAngle);
        mom.SetMag(data.trueMomentum*1E-3);

        // Check momentum: RKTrackRep can handle if > 4 MeV
        if(mom.Mag()*1E3 < 5)
        {
            processedEvents--;
            continue;
        }

        // Compute the super emission angle theta for future analysis
        Double_t x = data.fOrigin->X()*1E-1;
        Double_t y = data.fOrigin->Y()*1E-1;
        Double_t pz = cos(data.polarAngle);
        Double_t pr = sin(data.polarAngle) * (x * cos(data.azimuthalAngle) + y * sin(data.azimuthalAngle)) / sqrt(x*x + y*y);

        Double_t theThetaAngle = TMath::ATan2(pz, pr); // angle in plane (e_r, z) in radiants, in (-pi, pi)

        if(!opts.processAll)
        {
            cout << Form(">>> Vertex position = (%f, %f, %f) cm", pos[0], pos[1], pos[2]) << endl;
            cout << Form(">>> Vertex momentum = (%f, %f, %f) MeV/c", mom[0]*1E3, mom[1]*1E3, mom[2]*1E3) << endl;
            cout << Form(">>> (p, theta, phi) = (%f MeV/c, %f pi rad, %f pi rad)", mom.Mag()*1E3, theThetaAngle / TMath::Pi(), data.azimuthalAngle / TMath::Pi()) << endl;
        }

        // Check acceptance
        if(nHits < 3)
        {
            if(!opts.processAll)
                cout << ">>> Track is not in acceptance!" << endl;

            data.accTheta->Fill(false, data.trueMomentum, theThetaAngle);
            data.accPhi->Fill(false, data.trueMomentum, data.azimuthalAngle);

            continue;
        }
        inAcceptance++;

        // Sort and if in single event mode draw hits
            // will need a revision when origin point is not in Z = 0 anymore
        auto sortedHits = PTTALG::SortVectorZ(*(data.hitsCoordinates), *(data.cylinderID));
        *(data.hitsCoordinates) = sortedHits.first;
        *(data.cylinderID) = sortedHits.second;
        nTurns = PTTALG::CountTurns(*(data.hitsCoordinates));
        nCylinders = PTTALG::CountCylinders(*(data.cylinderID));

        // Apply turn analysis
        if(opts.turnMode)
        {
            if(!opts.processAll)
            {
                cout << ">>> nTurns = " << nTurns << endl;
                cout << ">>> nCylinders = " << nCylinders << endl;
            }

            auto turnHits = PTTALG::SelectTurn(opts.turnID, *(data.hitsCoordinates), *(data.cylinderID));
            *(data.hitsCoordinates) = turnHits.first;
            *(data.cylinderID) = turnHits.second;

            nHits = data.hitsCoordinates->size();

            // Check acceptance again
            if(nHits < 3)
            {
                if(!opts.processAll)
                    cout << ">>> Track is not in acceptance anymore!" << endl;

                inAcceptance--;

                data.accTheta->Fill(false, data.trueMomentum, theThetaAngle);
                data.accPhi->Fill(false, data.trueMomentum, data.azimuthalAngle);

                continue;
            }
        }

        // Track is in acceptance!
        data.accTheta->Fill(true, data.trueMomentum, theThetaAngle);
        data.accPhi->Fill(true, data.trueMomentum, data.azimuthalAngle);

        
        // Fast detector simulation
        if(opts.useSmearing)
        {
            for(Int_t i = 0; i < nHits; i++)
                data.hitsCoordinates->at(i) = AUXALG::SmearMeasurement((*(data.cylinderID))[i], data.hitsCoordinates->at(i));
        }


        // Prefitter
        auto helixPars = HelixPrefitter(*(data.hitsCoordinates), *(data.cylinderID), opts);
        Double_t xC = helixPars[0],
                 yC = helixPars[1],
                 R = helixPars[2],
                 phi0 = helixPars[3],
                 z0 = helixPars[4],
                 tanLambda = helixPars[5];

        // Cumulative arc length computation
        vector<Double_t> s_cumulative;
        Double_t previous_phi = phi0;
        Double_t previous_s = 0;

        for(const auto& point : *(data.hitsCoordinates))
        {
            Double_t dx = point[0] - xC;
            Double_t dy = point[1] - yC;
            Double_t phi = atan2(dy, dx);
            Double_t dphi = phi - previous_phi;
            if(dphi > M_PI) dphi -= 2 * M_PI;
            if(dphi < -M_PI) dphi += 2 * M_PI;
            dphi *= -1;
            Double_t s = previous_s + R * dphi;
            s_cumulative.push_back(s);
            previous_phi = phi;
            previous_s = s;
        }

        // Track candidate
        genfit::TrackCand trackCand;

        // Resolution of detectors
        const CHeT::Resolutions hitCov;
        
        // Fill the candidate
        vector<TVector3> measuredCoordinates;
        Int_t hitID = 0;
        Int_t nAddedHits = 0;
        for(Int_t i = 0; i < nHits; i++)
        {   
            TVector3 hitCoords;
            vector<Double_t> measuredCoords = data.hitsCoordinates->at(i);

            hitCoords[0] = measuredCoords.at(0);
            hitCoords[1] = measuredCoords.at(1);
            hitCoords[2] = measuredCoords.at(2);

            measuredCoordinates.push_back(hitCoords);

            new(chetHitArray[hitID]) genfit::mySpacepointDetectorHit(
                hitCoords,
                hitCov.GetMatrixCartesian((*(data.cylinderID))[i], TMath::ATan2(hitCoords[1], hitCoords[0]))
            );
            
            trackCand.addHit(detId, hitID, -1, i);
            hitID++;

            if(((i + 1) < nHits) && ((s_cumulative[i+1] - s_cumulative[i])/ R > TMath::Pi()/4))
            {
                const Double_t s_middle = 0.5 * (s_cumulative[i] + s_cumulative[i + 1]);
                AUXALG::AddFakeHitFromHelix(trackCand, hitID, i + 0.5,
                                            s_middle, xC, yC, R, z0, phi0, tanLambda,
                                            chetHitArray, 1.);
                nAddedHits++;
                hitID++;
            }
        }
        
        // Plotting
        if(!opts.processAll)
        {
            vector<vector<Double_t>> plottedCoordsVec;
            for(const auto& vec : measuredCoordinates)
                plottedCoordsVec.push_back({vec.X(), vec.Y(), vec.Z()});

            TCanvas *canvHitsXY = new TCanvas("canvHitsXY", "canvHitsXY", 700, 700);
            TCanvas *canvHitsYZ = new TCanvas("canvHitsYZ", "canvHitsYZ", 900, 500);
            TCanvas *canvHitsXYZ = new TCanvas("canvHitsXYZ");

            AUXALG::DrawXYView_hits(data.fOrigin, plottedCoordsVec, canvHitsXY);
            AUXALG::DrawYZView_hits(data.fOrigin, plottedCoordsVec, canvHitsYZ);
            AUXALG::DrawXYZView_hits(plottedCoordsVec, canvHitsXYZ);
        }


        // Start values for fit (would come from pattern recognition)
        // Initial guess for cov
        const Double_t seedResolution = 0.1; // ?
        TMatrixDSym covSeed(6);
        for(Int_t i = 0; i < 3; i++)
            covSeed(i,i) = seedResolution*seedResolution;
        for(Int_t i = 3; i < 6; i++)
            covSeed(i,i) = pow(seedResolution / nHits / sqrt(3), 2);

        // Set start values
        TVector3 posSeed(pos);
        TVector3 momSeed(mom);
        if(opts.pttrecMode)
            tie(posSeed, momSeed) = PTTALG::SmearSeed(pos, mom);

        if(!opts.processAll)
        {
            cout << endl;
            cout << Form(">>> Seed position = (%f, %f, %f) cm", posSeed[0], posSeed[1], posSeed[2]) << endl;
            cout << Form(">>> Seed momentum = (%f, %f, %f) MeV", momSeed[0]*1E3, momSeed[1]*1E3, momSeed[2]*1E3) << endl;
        }
        
        trackCand.setPosMomSeedAndPdgCode(posSeed, momSeed, pdg);
        trackCand.setCovSeed(covSeed);


        // Track rep
        genfit::AbsTrackRep *rep = new genfit::RKTrackRep(pdg);

        // Create track
        fitTrack = new genfit::Track(trackCand, chetFactory, rep);

        // Check
        fitTrack->checkConsistency();
        //fitTrack->Print();

        // Do the fit
        try
        {
            fitter->processTrack(fitTrack, true);
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
            
            auto [fRes, fSigma, fPulls] = AUXALG::GetResults(fitTrack, rep, (*data.fOrigin)*1E-1, data.trueMomentum, theThetaAngle, data.azimuthalAngle, NORMALIZED_PULLS);
            if(fRes.size() == 0) continue;
            data.histDiffX->Fill(fPulls[0]);
            data.histDiffY->Fill(fPulls[1]);
            data.histDiffZ->Fill(fPulls[2]);
            data.histDiffMom->Fill(fPulls[3]);
            data.histDiffTheta->Fill(fPulls[4]);
            data.histDiffPhi->Fill(fPulls[5]);

            data.graphMom->Fill(data.trueMomentum, fRes[3]);
            data.graphTheta->Fill(theThetaAngle, fRes[4]);
            data.graphPhi->Fill(data.azimuthalAngle, fRes[5]);

            data.hist2MomRes->Fill(data.trueMomentum, fSigma[3]);
            data.hist2ThetaRes->Fill(theThetaAngle, fSigma[4]);
            data.hist2PhiRes->Fill(data.azimuthalAngle, fSigma[5]);

            data.profMomRes->Fill(data.trueMomentum, fSigma[3]);
            data.profThetaRes->Fill(theThetaAngle, fSigma[4]);
            data.profPhiRes->Fill(data.azimuthalAngle, fSigma[5]);

            if(!opts.processAll)
                cout << Form(">>> Fitted (p, theta, phi) = (%f MeV/c, %f pi rad, %f pi rad)", mom.Mag()*1E3, fRes[4] / TMath::Pi(), fRes[5] / TMath::Pi()) << endl;
        }


        if(!opts.processAll)
        {
            cout << "\n\n>>> Did FIT converge? " << (isFitConverged ? "Yes" : "No") << "\n\n" << endl;
        }

        // Track is in efficiency?
        data.effTheta->Fill(isFitConverged, data.trueMomentum, theThetaAngle);
        data.effPhi->Fill(isFitConverged, data.trueMomentum, data.azimuthalAngle);

        // Store turns and cylinders data
        data.histTurns->Fill(nTurns);
        data.effTurns->Fill(isFitConverged, nTurns);
        data.histCylinders->Fill(nCylinders);
        data.effCylinders->Fill(isFitConverged, nCylinders);
        data.histTurnsVMom->Fill(data.trueMomentum, nTurns);
        data.histCylVMom->Fill(data.trueMomentum, nCylinders);

        // Check
        fitTrack->checkConsistency();

        // Last steps
        display->addEvent(fitTrack);

        cout << "\r>>> Processed event number " << ev << flush;
    }
    cout << endl;

    // Print results and draw graphs
    if(opts.processAll)
    {
        // Recap
        cout << "\n---------------------------------------------------------" << endl;
        cout << ">>> Acceptance = " << (Float_t) (inAcceptance * 100) / processedEvents << endl;
        cout << ">>> Efficiency = " << (Float_t) (inEfficiency * 100) / inAcceptance << endl;
        cout << "---------------------------------------------------------\n" << endl;

        // Graphs
        if(opts.quietMode)
            gROOT->SetBatch(true);

        TCanvas *canvEfficiency = new TCanvas("canvEfficiency");
        canvEfficiency->Divide(2,2);
        canvEfficiency->cd(1);
        data.accTheta->Draw("COLZ TEXT");
        canvEfficiency->cd(2);
        data.accPhi->Draw("COLZ TEXT");
        canvEfficiency->cd(3);
        data.effTheta->Draw("COLZ TEXT");
        canvEfficiency->cd(4);
        data.effPhi->Draw("COLZ TEXT");
    
        TCanvas *canvTurns = new TCanvas("canvTurns");
        canvTurns->Divide(2);
        canvTurns->cd(1);
        data.histTurns->Draw();
        canvTurns->cd(2);
        data.effTurns->Draw("AP");

        TCanvas *canvCylinders = new TCanvas("canvCylinders");
        canvCylinders->Divide(2, 2);
        canvCylinders->cd(1);
        data.histCylinders->Draw();
        canvCylinders->cd(2);
        data.effCylinders->Draw("AP");
        canvCylinders->cd(3);
        data.histTurnsVMom->Draw();
        canvCylinders->cd(4);
        data.histCylVMom->Draw();
        
        TCanvas *canvLinearity = new TCanvas("canvLinearity");
        canvLinearity->Divide(3);
        canvLinearity->cd(1);
        data.graphMom->SetMarkerStyle(20);
        data.graphMom->Draw("SCAT");
        canvLinearity->cd(2);
        data.graphTheta->SetMarkerStyle(20);
        data.graphTheta->Draw("SCAT");
        canvLinearity->cd(3);
        data.graphPhi->SetMarkerStyle(20);
        data.graphPhi->Draw("SCAT");

        TCanvas *canvfPullsPos = new TCanvas("Pulls Position");
        canvfPullsPos->Divide(3);
        
        canvfPullsPos->cd(1);
        data.histDiffX->Draw();
        
        canvfPullsPos->cd(2);
        data.histDiffY->Draw();
        
        canvfPullsPos->cd(3);
        data.histDiffZ->Draw();
        
        TCanvas *canvfPullsMom = new TCanvas("Pulls Momentum");   
        canvfPullsMom->Divide(3);
        
        canvfPullsMom->cd(1);
        data.histDiffMom->Draw();

        canvfPullsMom->cd(2);
        data.histDiffTheta->Draw();

        canvfPullsMom->cd(3);
        data.histDiffPhi->Draw();

        TCanvas *canvProfRes = new TCanvas("ProfResolutions");
        canvProfRes->Divide(3);
        canvProfRes->cd(1);
        data.profMomRes->Draw();
        canvProfRes->cd(2);
        data.profThetaRes->Draw();
        canvProfRes->cd(3);
        data.profPhiRes->Draw();
    
        TCanvas *canvResMom = new TCanvas("canvResMom");
        canvResMom->cd();
        data.hist2MomRes->Draw();

        TCanvas *canvResTheta = new TCanvas("canvResTheta");
        canvResTheta->cd();
        data.hist2ThetaRes->Draw();

        TCanvas *canvResPhi = new TCanvas("canvResPhi");
        canvResPhi->cd();
        data.hist2PhiRes->Draw();

        if(opts.quietMode)
        {
            canvEfficiency->SaveAs("canvEfficiency.pdf");
            canvResMom->SaveAs("canvResMom.pdf");
            canvResTheta->SaveAs("canvResTheta.pdf");
            canvResPhi->SaveAs("canvResPhi.pdf");
        }
    }

    // Delete fitter
    delete fitter;

    // Open event display
    display->setOptions("ABDEFGHMPT");
    display->open();

    // Finally
    exit(0);
}



void FITALG::HelixFitter(Options opts)
{
    // Load events
    TrackDataManager data(opts.eventMax);

    Int_t nEvents = data.GetNEvents();

    genfit::EventDisplay* display = genfit::EventDisplay::getInstance();

    // Get a specific event for single view mode
    gRandom->SetSeed(0);
    if(opts.event == -1)
        opts.event = gRandom->Integer(nEvents);

    // Event loop
    Int_t inAcceptance = 0; Int_t inEfficiency = 0;
    Int_t nTurns, nCylinders;
    for(Int_t ev = 0; ev < nEvents; ev++)
    {
        if(!opts.processAll)
            if(ev != opts.event)
                continue;

        // Get event
        data.GetChain()->GetEntry(ev);
        Int_t nHits = data.hitsCoordinates->size();

        if(!opts.processAll)
        {
            cout << "\n>>> Event number = " << ev << endl;
            cout << ">>> NHits = " << nHits << endl;
        }

        // Get True MC data
        TVector3 pos = {data.fOrigin->X()*1E-1, data.fOrigin->Y()*1E-1, data.fOrigin->Z()*1E-1};
        TVector3 mom = {1, 0, 0};
        mom.SetPhi(data.azimuthalAngle);
        mom.SetTheta(data.polarAngle);
        mom.SetMag(data.trueMomentum*1E-3);

        // Compute the super emission angle theta for future analysis
        Double_t x = data.fOrigin->X()*1E-1;
        Double_t y = data.fOrigin->Y()*1E-1;
        Double_t pz = cos(data.polarAngle);
        Double_t pr = sin(data.polarAngle) * (x * cos(data.azimuthalAngle) + y * sin(data.azimuthalAngle)) / sqrt(x*x + y*y);

        Double_t theThetaAngle = TMath::ATan2(pz, pr); // angle in plane (e_r, z) in radiants, in (-pi, pi)

        if(!opts.processAll)
        {
            cout << Form(">>> Vertex position = (%f, %f, %f) cm", pos[0], pos[1], pos[2]) << endl;
            cout << Form(">>> Vertex momentum = (%f, %f, %f) MeV/c", mom[0]*1E3, mom[1]*1E3, mom[2]*1E3) << endl;
            cout << Form(">>> (p, theta, phi) = (%f MeV/c, %f pi rad, %f pi rad)", mom.Mag()*1E3, theThetaAngle / TMath::Pi(), data.azimuthalAngle / TMath::Pi()) << endl;
        }

        // Check acceptance
        if(nHits < 3)
        {
            if(!opts.processAll)
                cout << ">>> Track is not in acceptance!" << endl;

            data.accTheta->Fill(false, data.trueMomentum, theThetaAngle);
            data.accPhi->Fill(false, data.trueMomentum, data.azimuthalAngle);

            continue;
        }
        inAcceptance++;

        // Sort and if in single event mode draw hits
            // will need a revision when origin point is not in Z = 0 anymore
        auto sortedHits = PTTALG::SortVectorZ(*(data.hitsCoordinates), *(data.cylinderID));
        *(data.hitsCoordinates) = sortedHits.first;
        *(data.cylinderID) = sortedHits.second;
        nTurns = PTTALG::CountTurns(*(data.hitsCoordinates));
        nCylinders = PTTALG::CountCylinders(*(data.cylinderID));

        // Apply turn analysis
        if(opts.turnMode)
        {
            if(!opts.processAll)
            {
                cout << ">>> nTurns = " << nTurns << endl;
                cout << ">>> nCylinders = " << nCylinders << endl;
            }

            auto turnHits = PTTALG::SelectTurn(opts.turnID, *(data.hitsCoordinates), *(data.cylinderID));
            *(data.hitsCoordinates) = turnHits.first;
            *(data.cylinderID) = turnHits.second;

            nHits = data.hitsCoordinates->size();

            // Check acceptance again
            if(nHits < 3)
            {
                if(!opts.processAll)
                    cout << ">>> Track is not in acceptance anymore!" << endl;

                inAcceptance--;

                data.accTheta->Fill(false, data.trueMomentum, theThetaAngle);
                data.accPhi->Fill(false, data.trueMomentum, data.azimuthalAngle);

                continue;
            }
        }

        // Track is in acceptance!
        data.accTheta->Fill(true, data.trueMomentum, theThetaAngle);
        data.accPhi->Fill(true, data.trueMomentum, data.azimuthalAngle);


        // Resolution of detectors
        const CHeT::Resolutions hitCov;

        // Fill the candidate
        vector<TVector3> measuredCoordinates;

        RVecD r_1, r_2;
        RVecD w;
        RVec<RVecD> V(2*nHits, RVecD(2*nHits));
        
        for(Int_t i = 0; i < nHits; i++)
        {
            TVector3 hitCoords;

            // Measurements
            vector<Double_t> measuredCoords;
            if(opts.useSmearing)
                measuredCoords = AUXALG::SmearMeasurement((*(data.cylinderID))[i], data.hitsCoordinates->at(i));
            else
                measuredCoords = data.hitsCoordinates->at(i);

            hitCoords[0] = measuredCoords.at(0);
            hitCoords[1] = measuredCoords.at(1);
            hitCoords[2] = measuredCoords.at(2);

            measuredCoordinates.push_back(hitCoords);

            Double_t uu = measuredCoords.at(0);
            Double_t vv =  measuredCoords.at(1);
            Double_t phi_global = atan2(vv, uu);

            TMatrixDSym cov_xy = hitCov.GetMatrixCartesian((*(data.cylinderID))[i], phi_global).GetSub(0,1,0,1);
            
            TMatrixDSym cov_rphiz = hitCov.GetMatrixCylindrical((*(data.cylinderID))[i]); 
            
            // Fill m, V and w
            r_1.push_back(uu);
            r_2.push_back(vv);
            w.push_back(1./cov_rphiz(1,1));

            for(auto j = 0; j < 2; j++)
                for(auto k = 0; k < 2; k++)
                    V[nHits*j + i][nHits*k + i] = cov_xy(j,k);
        }
        
        // ... Centering ...
        //u -= Mean(u);
        //v -= Mean(v);
        RVecD m_c = Concatenate(r_1,r_2);
        
        // ... Scaling ...
        //Double_t b = 0.5;
        //auto q = Dot(m_c, m_c);
        //auto Q = sqrt(q/nHits);
        //auto m_cs = m_c * b/Q;
        
        TMatrixD V_11(nHits, nHits),
                 V_12(nHits, nHits),
                 V_21(nHits, nHits),
                 V_22(nHits, nHits);
        
        for(auto i = 0; i < nHits; i++)
            for(auto j = 0; j < nHits; j++)
            {
                V_11(i,j) = V[i][j];
                V_12(i,j) = V[i][nHits + j];
                V_21(i,j) = V[nHits + i][j];
                V_22(i,j) = V[nHits + i][nHits + j];
            }
        
        // ... Mapping ...
        RVecD r_3 = r_1*r_1 + r_2*r_2;
        
        RVecD r = Concatenate(m_c, r_3);
        
        // C Matrix
        map<pair<Int_t,Int_t>, TMatrixD> C;
        for(auto i = 1; i <= 3; ++i)
            for(auto j = 1; j <= 3; ++j)
                C.insert({{i, j}, TMatrixD(nHits, nHits)});
        
        C[{1,1}] = V_11;
        C[{1,2}] = V_12;
        C[{2,1}] = V_21;
        C[{2,2}] = V_22;
        
        // C_13, C_23
        TMatrixD C13(nHits, nHits), C23(nHits, nHits);
        for(auto i = 0; i < nHits; ++i)
            for(auto j = 0; j < nHits; ++j)
            {
                C13(i,j) = 2*V_11(i,j)*r_1[j] + 2*V_12(i,j)*r_2[j];
                C23(i,j) = 2*V_21(i,j)*r_1[j] + 2*V_22(i,j)*r_2[j];
            }

        C[{1,3}] = C13;
        C[{2,3}] = C23;
        C[{3,1}] = TMatrixD(TMatrixD::kTransposed, C13);
        C[{3,2}] = TMatrixD(TMatrixD::kTransposed, C23);
        
        // C_33
        TMatrixD C33(nHits, nHits);
        for(auto i = 0; i < 2; ++i)
            for(auto j = 0; j < 2; ++j)
            {
                const TMatrixD &Vii = (i == 0 ? V_11 : V_22);
                const TMatrixD &Vij = (i == 0 && j == 0) ? V_11 :
                                      (i == 0 && j == 1) ? V_12 :
                                      (i == 1 && j == 0) ? V_21 : V_22;

                const RVecD &ri = (i == 0 ? r_1 : r_2);
                const RVecD &rj = (j == 0 ? r_1 : r_2);

                for(auto m = 0; m < nHits; ++m)
                    for(auto n = 0; n < nHits; ++n)
                        C33(m,n) += 2*Vii(m,n)*Vij(m,n) + 4*Vij(m,n)*ri[m]*rj[n];
            }

        C[{3,3}] = C33;

        
        // ... Center of gravity ...
        w /= Sum(w);
        TVectorD w_vec(w.size());
        for(auto i = 0; i < w.size(); ++i)
            w_vec(i) = w[i];

        TMatrixD r_mat(nHits, 3);
        for(auto j = 0; j < 3; ++j)
            for(auto i = 0; i < nHits; ++i)    
                r_mat(i,j) = r[j*nHits + i];
        
        TMatrixD r_mat_tr(TMatrixD::kTransposed, r_mat);
        TVectorD r_0 = r_mat_tr * w_vec;
        
        // Var(r_0)
        TMatrixD C_0(3,3);

        for(auto i = 1; i <= 3; ++i)
        {
            for(auto j = 1; j <= 3; ++j)
            {
                const TMatrixD &Cij = C[{i,j}];
                C_0(i-1,j-1) = Cij.Similarity(w_vec);
            }
        }

        // ... Substract ...
        TMatrixD H(nHits, nHits);
        for(auto i = 0; i < nHits; ++i)
            for(auto j = 0; j < nHits; ++j)
                H(i,j) = (i == j ? 1 : 0) - w_vec(j); 

        // s and D Matrix
        TMatrixD s = H*r_mat;
        
        map<Int_t, TVectorD> s_map;
        for(auto i = 1; i <= 3; ++i)
            s_map.insert({i, TVectorD(nHits)});

        for(auto i = 0; i < nHits; ++i)
        {
            s_map[1](i) = s(i,0);
            s_map[2](i) = s(i,1);
            s_map[3](i) = s(i,2);
        }
        
        
        map<pair<Int_t,Int_t>, TMatrixD> D;
        for(auto i = 1; i <= 3; ++i)
            for(auto j = 1; j <= 3; ++j)
                D.insert({{i, j}, TMatrixD(nHits, nHits)});

        for(auto i = 1; i <= 3; ++i)
        {
            for(auto j = 1; j <= 3; ++j)
            {
                TMatrixD &Dij = D[{i,j}];
                const TMatrixD &Cij = C[{i,j}];
                TMatrixD H_tr(TMatrixD::kTransposed, H);
                
                Dij = H * Cij * H_tr;
            }
        }
        
        // ... Computation of weighted sample covariance matrix 𝑨 ...
        map<Int_t, pair<Int_t,Int_t>> nu = {
            {1, {1,1}},
            {2, {1,2}},
            {3, {1,3}},
            {4, {2,2}},
            {5, {2,3}},
            {6, {3,3}}
        };
        
        map<Int_t, Double_t> A_alpha;
        for(auto alpha = 1; alpha <= 6; ++alpha)
        {
            auto [i, j] = nu[alpha];
            TVectorD w_sj(nHits);
            for(auto k = 0; k < w.size(); ++k)
                w_sj[k] = w_vec(k) * s_map[j](k); // w ⊙ s_j
    
            A_alpha[alpha] = s_map[i] * w_sj; // s_iᵀ ⋅ (w ⊙ s_j)
        }
        
        
        // Compute covariance matrix of A E_{alpha,beta} for alpha, beta = 1..6
        TMatrixDSym E(6);
        
        // Precompute W2 = w * w^T (outer product)
        TMatrixD W2(nHits, nHits);
        for(auto i = 0; i < nHits; ++i)
            for(auto j = 0; j < nHits; ++j)
                W2(i,j) = w_vec(i) * w_vec(j);
        
        for(auto alpha = 1; alpha <= 6; ++alpha)
        {
            auto [i, j] = nu[alpha];
            const auto& si = s_map[i];
            const auto& sj = s_map[j];
            
            for(auto beta = alpha; beta <= 6; ++beta)
            {
                auto [k, l] = nu[beta];
                const auto& sk = s_map[k];
                const auto& sl = s_map[l];

                const TMatrixD& Dik = D[{i,k}];
                const TMatrixD& Dil = D[{i,l}];
                const TMatrixD& Djk = D[{j,k}];
                const TMatrixD& Djl = D[{j,l}];
            
                // Hadamard products
                Double_t S_term = 0.;

                for(auto qi = 0; qi < nHits; ++qi)
                    for(auto qj = 0; qj < nHits; ++qj)
                        S_term += (Dik(qi,qj) * W2(qi,qj) * Djl(qi, qj) + Dil(qi,qj) * W2(qi,qj) * Djk(qi, qj));

                TMatrixD temp_jl2(nHits, nHits), temp_jk2(nHits, nHits), temp_il2(nHits, nHits), temp_ik2(nHits, nHits);
                for(auto qu = 0; qu < nHits; ++qu)
                    for(auto qv = 0; qv < nHits; ++qv)
                    {
                        temp_jl2(qu,qv) = Djl(qu,qv) * W2(qu,qv);
                        temp_jk2(qu,qv) = Djk(qu,qv) * W2(qu,qv);
                        temp_il2(qu,qv) = Dil(qu,qv) * W2(qu,qv);
                        temp_ik2(qu,qv) = Dik(qu,qv) * W2(qu,qv);
                    }
                
                Double_t scalar = si * (temp_jl2 * sk) + si * (temp_jk2 * sl) + sj * (temp_il2 * sk) + sj * (temp_ik2 * sl);
                
                // Finally
                E(alpha - 1, beta - 1) = S_term + scalar;
                
                if(alpha != beta)
                    E(beta - 1, alpha - 1) = E(alpha - 1, beta - 1); // symmetry
            }
        }

        //E.Print();
        
        // ... Computation of n and c
        // A matrix
        TMatrixDSym A(3); A.Zero();
        for(auto alpha = 1; alpha <= 6; ++alpha)
        {
            auto [i, j] = nu[alpha];
            A(i-1,j-1) = A_alpha[alpha];
            A(j-1,i-1) = A_alpha[alpha];
        }
        
        // Diagonalizzazione
        TVectorD eigenVals(3);
        TMatrixD eigenVecs = A.EigenVectors(eigenVals);

        // Normale = autovettore con autovalore minimo
        Int_t minIdx = (eigenVals(0) < eigenVals(1)) ?
                        ((eigenVals(0) < eigenVals(2)) ? 0 : 2) :
                        ((eigenVals(1) < eigenVals(2)) ? 1 : 2);

        TVectorD n(3);
        for(Int_t i = 0; i < 3; ++i)
            n[i] = eigenVecs(i, minIdx);
        
        Double_t c = -(n*r_0);
        
        // Compute the Jacobian
        TMatrixD J2(3, 6);
        const Double_t epsilon = 1e-2;
        
        for(auto alpha = 1; alpha <= 6; ++alpha)
        {
            // Copia A_alpha e perturba il solo alpha-esimo parametro
            auto A_plus = A_alpha;
            auto A_minus = A_alpha;

            A_plus[alpha] += epsilon;
            A_minus[alpha] -= epsilon;

            auto calc_n = [&](const map<Int_t, Double_t>& A_mod) -> TVectorD
            {
                TMatrixDSym A_mat(3); A_mat.Zero();
                for(auto a = 1; a <= 6; ++a)
                {
                    auto [i, j] = nu[a];
                    A_mat(i-1,j-1) = A_mod.at(a);
                    A_mat(j-1,i-1) = A_mod.at(a);
                }

                TVectorD evals(3);
                TMatrixD evecs = A_mat.EigenVectors(evals);

                Int_t minIdx = (evals(0) < evals(1)) ? ((evals(0) < evals(2)) ? 0 : 2)
                                           : ((evals(1) < evals(2)) ? 1 : 2);

                TVectorD n(3);
                
                for(auto i = 0; i < 3; ++i)
                    n[i] = evecs(i, minIdx);

                return n;
            };

            TVectorD n_plus  = calc_n(A_plus);
            TVectorD n_minus = calc_n(A_minus);

            for(auto i = 0; i < 3; ++i)
                J2(i, alpha - 1) = (n_plus[i] - n_minus[i]) / (2. * epsilon);
        }

        //TMatrixDSym C_n(3, 3);
        auto C_n = E.Similarity(J2);      // J2 * E * J2ᵀ

        
        // Joint covariance matrix of n and c
        // Parte alta-sinistra: Cn
        TMatrixD Cnc(4, 4);
        for(auto i = 0; i < 3; ++i)
            for(auto j = 0; j < 3; ++j)
                Cnc(i, j) = C_n(i, j);

        // Parte in alto a destra: -Cn * r0
        TVectorD Cn_r0(3);
        Cn_r0 = C_n * r_0;

        for(auto i = 0; i < 3; ++i)
        {
            Cnc(i, 3) = -Cn_r0[i];       // colonna finale
            Cnc(3, i) = -Cn_r0[i];       // riga finale (simmetrico)
        }

        // Calcolo var[c]
        Double_t ncnrcr = C_0.Similarity(n) + C_n.Similarity(r_0);
        Double_t S_trace = 0.;
        for(auto i = 0; i < 3; ++i)
            for(auto j = 0; j < 3; ++j)
                S_trace += C_n(i,j) * C_0(i,j); // Hadamard product

        Double_t var_c = ncnrcr + S_trace;
        Cnc(3, 3) = var_c;

        
        // ... Circle parameters ...
        const Double_t xC = -n(0) / (2 * n(2));
        const Double_t yC = -n(1) / (2 * n(2));
        Double_t r2 = (1 - n(2)*n(2) - 4 * c * n(2)) / (4 * n(2)*n(2));
        const Double_t R  = sqrt(r2);
        
        // Jacobian
        Double_t h = sqrt(1 - n(2)*n(2) - 4*c*n(2));
        TMatrixD J3(3,4);   J3.Zero();
        J3(0,0) = -1./(2*n(2));     J3(0,2) = n(0)/(2*n(2)*n(2));
        J3(1,1) = -1./(2*n(2));     J3(1,2) = n(1)/(2*n(2)*n(2));
        J3(2,2) = -h/(2*n(2)*n(2)) - (4*c + 2*n(2)) / (4*h*n(2));
        J3(2,3) = -1./h;
        
        
        // Covariance matrix
        TMatrixD J3_tr = TMatrixD(TMatrixD::kTransposed, J3);
        TMatrixD covCircle(3,3);
        covCircle = J3 * Cnc * J3_tr;
        
        if(!opts.processAll)
            cout << Form("\nxC = %f +/- %f\nyC = %f +/- %f\nR = %f +/- %f\n\n", xC, sqrt(covCircle(0,0)), yC, sqrt(covCircle(1,1)), R, sqrt(covCircle(2,2)));
        

        
        // --- Fit helix in Z vs arc length s ---
        vector<Double_t> s_values;
        vector<Double_t> z_values;

        // Choose the pivot
        const TVector3 pivot = measuredCoordinates[0];

        // Compute phi0 and dr
        const Double_t phi0 = atan2(pivot.Y() - yC, pivot.X() - xC);
        const Double_t cos_phi0 = cos(phi0);
        const Double_t sin_phi0 = sin(phi0);
        const Double_t dr = sqrt(pow(pivot.X() - xC, 2) + pow(pivot.Y() - yC, 2)) - R;

        // Compute arc lengths s
        Double_t previous_phi = phi0;
        Double_t previous_s = 0;
        for(const auto& point : measuredCoordinates)
        {
            Double_t dx = point.X() - xC;
            Double_t dy = point.Y() - yC;
            Double_t phi = atan2(dy, dx);

            // Angle unwrapping
            Double_t dphi = phi - previous_phi;
            if(dphi > TMath::Pi()) dphi -= 2 * TMath::Pi();
            if(dphi < -TMath::Pi()) dphi += 2 * TMath::Pi();
            dphi *= -1;

            Double_t s = previous_s + R * dphi;

            s_values.push_back(s);
            z_values.push_back(point.Z());

            previous_phi = phi;
            previous_s = s;
        }

        // --- Fit z vs s ---
        auto *graphZvsS = new TGraphErrors(s_values.size());
        for(size_t i = 0; i < s_values.size(); ++i)
        {
            Double_t s = s_values[i];
            Double_t z = z_values[i];

            // Get the hit            
            const auto& hit = measuredCoordinates[i];
            const TMatrixDSym& matCov = hitCov.GetMatrixCartesian((*(data.cylinderID))[i], atan2(hit.Y(), hit.X()));
            
            // Uncertainties on Z
            Double_t sigmaZ = sqrt(matCov(2, 2));

            graphZvsS->SetPoint(i, s, z);
            graphZvsS->SetPointError(i, 0, sigmaZ);
        }

        TF1 *fitZvsS = new TF1("fitZvsS", "[0] + x*[1]", -10, 10);
        fitZvsS->SetParNames("z0", "tan(lambda)");
        fitZvsS->SetParameters(0, 1);
        TFitResultPtr fitlinePtr = nullptr;
        if(!opts.processAll) 
            fitlinePtr= graphZvsS->Fit(fitZvsS, "S");
        else
            fitlinePtr= graphZvsS->Fit(fitZvsS, "SQ");

        // Get results
        Double_t z0 = fitZvsS->GetParameter(0);
        Double_t tanLambda = fitZvsS->GetParameter(1);
        
        TMatrixDSym covLine = fitlinePtr->GetCovarianceMatrix();
        
        
        // Helix Fit result
        Bool_t isFitConverged = fitlinePtr->IsValid();
        if(isFitConverged)
        {
            inEfficiency++;
            
            // Construct the full params vector
            TVectorD parHelix(5);
            parHelix(0) = xC; parHelix(1) = yC; parHelix(2) = R;
            parHelix(3) = z0; parHelix(4) = tanLambda;

            // Construct the full covariance matrix
            TMatrixDSym covHelix(5);
            for(Int_t i = 0; i < 3; ++i)
                for(Int_t j = 0; j < 3; ++j)
                    covHelix(i, j) = covCircle(i, j);
            for(Int_t i = 0; i < 2; ++i)
                for(Int_t j = 0; j < 2; ++j)
                    covHelix(i + 3, j + 3) = covLine(i, j);


            // Fit info and plots
            if(!opts.processAll)
            {
                // --- Output parameters ---
                cout << "\n\nHelix parameters:";
                parHelix.Print();
                cout << "Phi0 = " << phi0 << endl;

                // --- Covariance matrix ---
                covHelix.Print();
            
            
                // Convert hits
                vector<vector<Double_t>> plottedCoordsVec;
                for(const auto& vec : measuredCoordinates)
                    plottedCoordsVec.push_back({vec.X(), vec.Y(), vec.Z()});

                // Canvas
                TCanvas *canvHitsXY = new TCanvas("canvHitsXY", "XY View", 700, 700);
                TCanvas *canvHitsYZ = new TCanvas("canvHitsYZ", "YZ View", 900, 500);
                TCanvas *canvHitsXYZ = new TCanvas("canvHitsXYZ", "3D Helix Fit", 800, 600);
                TCanvas *canvZvsS = new TCanvas("canvZvsS", "z vs s", 700, 500);
            
                AUXALG::DrawXYView_hits(data.fOrigin, plottedCoordsVec, canvHitsXY);
                AUXALG::DrawYZView_hits(data.fOrigin, plottedCoordsVec, canvHitsYZ);
                AUXALG::DrawXYZView_hits(plottedCoordsVec, canvHitsXYZ);
            
                AUXALG::DrawXYView_arc(xC, yC, R, plottedCoordsVec, canvHitsXY, nTurns, opts.turnID);
                AUXALG::DrawZvsSFit(graphZvsS, fitZvsS, canvZvsS);
                AUXALG::DrawXYZView_helixFromHits(xC, yC, R, z0, phi0, tanLambda, plottedCoordsVec, canvHitsXYZ);
            }
        
            // Extrapolate status at decay vertex
            // Find s at z = 0
            Double_t s_at_z0 = -z0 / tanLambda;
            Double_t phi_at_z0 = s_at_z0 / R;
        
            // Vertex decay
            Double_t x_at_z0 = xC + R * cos(phi0 - phi_at_z0);
            Double_t y_at_z0 = yC + R * sin(phi0 - phi_at_z0);
                //Double_t x_at_z0 = pivot.X() - dr * cos_phi0 + R * (cos(phi0 - phi_at_z0) - cos_phi0);
                //Double_t y_at_z0 = pivot.Y() - dr * sin_phi0 + R * (sin(phi0 - phi_at_z0) - sin_phi0);
            // Vertex momentum
            const Double_t k = 2.99792458; // MeV/c * T * cm      
            TVector3 pFitted(sin(phi0 - phi_at_z0), -cos(phi0 - phi_at_z0), tanLambda);
            pFitted *= k * 2.2 * R;
        
            // Uncertainties
            TMatrixD Jac = ANS::ComputeHelixJacobian(pivot, xC, yC, R, z0, tanLambda, 2.2);
            TMatrixDSym covFittedState = covHelix.Similarity(Jac);
            covFittedState = ANS::CovFromCardinalToCylindricalMom(covFittedState, pFitted);
            
            // Compute angles with muEDM convention
            Double_t pFittedPhi = (pFitted.Phi() > 0 ) ? pFitted.Phi() : pFitted.Phi() + TMath::TwoPi();
            Double_t pFitted_z = cos(pFitted.Theta());
            Double_t pFitted_r = sin(pFitted.Theta()) * (x_at_z0 * cos(pFittedPhi) + y_at_z0 * sin(pFittedPhi)) / sqrt(x_at_z0*x_at_z0 + y_at_z0*y_at_z0);
            Double_t pFittedTheTheta = atan2(pFitted_z, pFitted_r);
            
            // Fitted final extrapolated state
                // x, y, z, p, TheTheta, phi
            TVectorD fittedState(6); 
            fittedState(0) = x_at_z0;
            fittedState(1) = y_at_z0;
            fittedState(2) = 0.;
            fittedState(3) = pFitted.Mag();
            fittedState(4) = pFittedTheTheta;
            fittedState(5) = pFittedPhi;

            // Print results
            if(!opts.processAll)
            {
                cout << "\nFitted vertex position: (" << fittedState(0) << ", " << fittedState(1) << ", 0) cm" << endl;
                cout << Form("Fidded p = (%.2f, %.2f, %.2f) MeV/c", pFitted.X(), pFitted.Y(), pFitted.Z()) << endl;
                cout << Form("(p, theta, phi) = (%.2f MeV/c, %.2f pi rad, %.2f pi rad)", fittedState(3), fittedState(4)/TMath::Pi(), fittedState(5)/TMath::Pi());
                cout << endl;

                covFittedState.Print();
            }

            // Fill histos
            auto [fRes, fSigma, fPulls] = AUXALG::GetResults(fittedState, covFittedState, (*data.fOrigin)*1E-1, data.trueMomentum, theThetaAngle, data.azimuthalAngle, NORMALIZED_PULLS);
            if(fRes.size() == 0)
                continue;

            data.histDiffX->Fill(fPulls[0]);
            data.histDiffY->Fill(fPulls[1]);
            data.histDiffZ->Fill(fPulls[2]);
            data.histDiffMom->Fill(fPulls[3]);
            data.histDiffTheta->Fill(fPulls[4]);
            data.histDiffPhi->Fill(fPulls[5]);

            data.graphMom->Fill(data.trueMomentum, fRes[3]);
            data.graphTheta->Fill(theThetaAngle, fRes[4]);
            data.graphPhi->Fill(data.azimuthalAngle, fRes[5]);

            data.hist2MomRes->Fill(data.trueMomentum, fSigma[3]);
            data.hist2ThetaRes->Fill(theThetaAngle, fSigma[4]);
            data.hist2PhiRes->Fill(data.azimuthalAngle, fSigma[5]);

            data.profMomRes->Fill(data.trueMomentum, fSigma[3]);
            data.profThetaRes->Fill(theThetaAngle, fSigma[4]);
            data.profPhiRes->Fill(data.azimuthalAngle, fSigma[5]);
        }
        
        if(!opts.processAll)
        {
            cout << "\n\n>>> Did FIT converge? " << (isFitConverged ? "Yes" : "No") << "\n\n" << endl;
        }

        // Track is in efficiency?
        data.effTheta->Fill(isFitConverged, data.trueMomentum, theThetaAngle);
        data.effPhi->Fill(isFitConverged, data.trueMomentum, data.azimuthalAngle);

        // Store turns and cylinders data
        data.histTurns->Fill(nTurns);
        data.effTurns->Fill(isFitConverged, nTurns);
        data.histCylinders->Fill(nCylinders);
        data.effCylinders->Fill(isFitConverged, nCylinders);
        data.histTurnsVMom->Fill(data.trueMomentum, nTurns);
        data.histCylVMom->Fill(data.trueMomentum, nCylinders);

        
        cout << "\r>>> Processed event number " << ev << flush;
        
    }
    
    cout << endl;

    // Print results and draw graphs
    if(opts.processAll)
    {
        // Recap
        cout << "\n---------------------------------------------------------" << endl;
        cout << ">>> Acceptance = " << (Float_t) (inAcceptance * 100) / nEvents << endl;
        cout << ">>> Efficiency = " << (Float_t) (inEfficiency * 100) / inAcceptance << endl;
        cout << "---------------------------------------------------------\n" << endl;

        // Graphs
        if(opts.quietMode)
            gROOT->SetBatch(true);

        TCanvas *canvEfficiency = new TCanvas("canvEfficiency");
        canvEfficiency->Divide(2,2);
        canvEfficiency->cd(1);
        data.accTheta->Draw("COLZ TEXT");
        canvEfficiency->cd(2);
        data.accPhi->Draw("COLZ TEXT");
        canvEfficiency->cd(3);
        data.effTheta->Draw("COLZ TEXT");
        canvEfficiency->cd(4);
        data.effPhi->Draw("COLZ TEXT");
    
        TCanvas *canvTurns = new TCanvas("canvTurns");
        canvTurns->Divide(2);
        canvTurns->cd(1);
        data.histTurns->Draw();
        canvTurns->cd(2);
        data.effTurns->Draw("AP");

        TCanvas *canvCylinders = new TCanvas("canvCylinders");
        canvCylinders->Divide(2, 2);
        canvCylinders->cd(1);
        data.histCylinders->Draw();
        canvCylinders->cd(2);
        data.effCylinders->Draw("AP");
        canvCylinders->cd(3);
        data.histTurnsVMom->Draw();
        canvCylinders->cd(4);
        data.histCylVMom->Draw();
        
        TCanvas *canvLinearity = new TCanvas("canvLinearity");
        canvLinearity->Divide(3);
        canvLinearity->cd(1);
        data.graphMom->SetMarkerStyle(20);
        data.graphMom->Draw("SCAT");
        canvLinearity->cd(2);
        data.graphTheta->SetMarkerStyle(20);
        data.graphTheta->Draw("SCAT");
        canvLinearity->cd(3);
        data.graphPhi->SetMarkerStyle(20);
        data.graphPhi->Draw("SCAT");

        TCanvas *canvfPullsPos = new TCanvas("Pulls Position");
        canvfPullsPos->Divide(3);
        
        canvfPullsPos->cd(1);
        data.histDiffX->Draw();
        
        canvfPullsPos->cd(2);
        data.histDiffY->Draw();
        
        canvfPullsPos->cd(3);
        data.histDiffZ->Draw();
        
        TCanvas *canvfPullsMom = new TCanvas("Pulls Momentum");   
        canvfPullsMom->Divide(3);
        
        canvfPullsMom->cd(1);
        data.histDiffMom->Draw();

        canvfPullsMom->cd(2);
        data.histDiffTheta->Draw();

        canvfPullsMom->cd(3);
        data.histDiffPhi->Draw();

        TCanvas *canvProfRes = new TCanvas("ProfResolutions");
        canvProfRes->Divide(3);
        canvProfRes->cd(1);
        data.profMomRes->Draw();
        canvProfRes->cd(2);
        data.profThetaRes->Draw();
        canvProfRes->cd(3);
        data.profPhiRes->Draw();

        TCanvas *canvResMom = new TCanvas("canvResMom");
        canvResMom->cd();
        data.hist2MomRes->Draw();

        TCanvas *canvResTheta = new TCanvas("canvResTheta");
        canvResTheta->cd();
        data.hist2ThetaRes->Draw();

        TCanvas *canvResPhi = new TCanvas("canvResPhi");
        canvResPhi->cd();
        data.hist2PhiRes->Draw();

        if(opts.quietMode)
        {
            canvEfficiency->SaveAs("canvEfficiency.pdf");
            canvResMom->SaveAs("canvResMom.pdf");
            canvResTheta->SaveAs("canvResTheta.pdf");
            canvResPhi->SaveAs("canvResPhi.pdf");
        }
    }

    // Trick
    display->open();
    
    // Finally
    exit(0);

}



array<Double_t, 6> FITALG::HelixPrefitter(const vector<vector<Double_t>>& hitsCoordinates, const vector<Int_t>& cylinders, Options opts)
{
    // Return
    array<Double_t, 6> fParameters = {0, 0, 0, 0, 0, 0};

    // Resolution of detectors
    const CHeT::Resolutions hitCov;

    // Fill the candidate
    vector<TVector3> measuredCoordinates;
    Int_t nHits = hitsCoordinates.size();

    RVecD r_1, r_2;
    RVecD w;
    RVec<RVecD> V(2*nHits, RVecD(2*nHits));
    
    for(Int_t i = 0; i < nHits; i++)
    {
        TVector3 hitCoords;

        // Measurements
        vector<Double_t> measuredCoords = hitsCoordinates[i];

        hitCoords[0] = measuredCoords.at(0);
        hitCoords[1] = measuredCoords.at(1);
        hitCoords[2] = measuredCoords.at(2);

        measuredCoordinates.push_back(hitCoords);

        Double_t uu = measuredCoords.at(0);
        Double_t vv =  measuredCoords.at(1);
        Double_t phi_global = atan2(vv, uu);

        TMatrixDSym cov_xy = hitCov.GetMatrixCartesian(cylinders[i], phi_global).GetSub(0,1,0,1);
        
        TMatrixDSym cov_rphiz = hitCov.GetMatrixCylindrical(cylinders[i]); 
        
        // Fill m, V and w
        r_1.push_back(uu);
        r_2.push_back(vv);
        w.push_back(1./cov_rphiz(1,1));

        for(auto j = 0; j < 2; j++)
            for(auto k = 0; k < 2; k++)
                V[nHits*j + i][nHits*k + i] = cov_xy(j,k);
    }

    RVecD m_c = Concatenate(r_1,r_2);
    
    TMatrixD V_11(nHits, nHits),
             V_12(nHits, nHits),
             V_21(nHits, nHits),
             V_22(nHits, nHits);
    
    for(auto i = 0; i < nHits; i++)
        for(auto j = 0; j < nHits; j++)
        {
            V_11(i,j) = V[i][j];
            V_12(i,j) = V[i][nHits + j];
            V_21(i,j) = V[nHits + i][j];
            V_22(i,j) = V[nHits + i][nHits + j];
        }
    
    // ... Mapping ...
    RVecD r_3 = r_1*r_1 + r_2*r_2;
    
    RVecD r = Concatenate(m_c, r_3);
    
    // C Matrix
    map<pair<Int_t,Int_t>, TMatrixD> C;
    for(auto i = 1; i <= 3; ++i)
        for(auto j = 1; j <= 3; ++j)
            C.insert({{i, j}, TMatrixD(nHits, nHits)});
    
    C[{1,1}] = V_11;
    C[{1,2}] = V_12;
    C[{2,1}] = V_21;
    C[{2,2}] = V_22;
    
    // C_13, C_23
    TMatrixD C13(nHits, nHits), C23(nHits, nHits);
    for(auto i = 0; i < nHits; ++i)
        for(auto j = 0; j < nHits; ++j)
        {
            C13(i,j) = 2*V_11(i,j)*r_1[j] + 2*V_12(i,j)*r_2[j];
            C23(i,j) = 2*V_21(i,j)*r_1[j] + 2*V_22(i,j)*r_2[j];
        }

    C[{1,3}] = C13;
    C[{2,3}] = C23;
    C[{3,1}] = TMatrixD(TMatrixD::kTransposed, C13);
    C[{3,2}] = TMatrixD(TMatrixD::kTransposed, C23);

    // C_33
    TMatrixD C33(nHits, nHits);
    for(auto i = 0; i < 2; ++i)
        for(auto j = 0; j < 2; ++j)
        {
            const TMatrixD &Vii = (i == 0 ? V_11 : V_22);
            const TMatrixD &Vij = (i == 0 && j == 0) ? V_11 :
                                  (i == 0 && j == 1) ? V_12 :
                                  (i == 1 && j == 0) ? V_21 : V_22;

            const RVecD &ri = (i == 0 ? r_1 : r_2);
            const RVecD &rj = (j == 0 ? r_1 : r_2);

            for(auto m = 0; m < nHits; ++m)
                for(auto n = 0; n < nHits; ++n)
                    C33(m,n) += 2*Vii(m,n)*Vij(m,n) + 4*Vij(m,n)*ri[m]*rj[n];
        }

    C[{3,3}] = C33;

    
    // ... Center of gravity ...
    w /= Sum(w);
    TVectorD w_vec(w.size());
    for(auto i = 0; i < w.size(); ++i)
        w_vec(i) = w[i];

    TMatrixD r_mat(nHits, 3);
    for(auto j = 0; j < 3; ++j)
        for(auto i = 0; i < nHits; ++i)    
            r_mat(i,j) = r[j*nHits + i];
    
    TMatrixD r_mat_tr(TMatrixD::kTransposed, r_mat);
    TVectorD r_0 = r_mat_tr * w_vec;
    
    // Var(r_0)
    TMatrixD C_0(3,3);

    for(auto i = 1; i <= 3; ++i)
    {
        for(auto j = 1; j <= 3; ++j)
        {
            const TMatrixD &Cij = C[{i,j}];
            C_0(i-1,j-1) = Cij.Similarity(w_vec);
        }
    }

    // ... Substract ...
    TMatrixD H(nHits, nHits);
    for(auto i = 0; i < nHits; ++i)
        for(auto j = 0; j < nHits; ++j)
            H(i,j) = (i == j ? 1 : 0) - w_vec(j); 

    // s and D Matrix
    TMatrixD s = H*r_mat;
    
    map<Int_t, TVectorD> s_map;
    for(auto i = 1; i <= 3; ++i)
        s_map.insert({i, TVectorD(nHits)});

    for(auto i = 0; i < nHits; ++i)
    {
        s_map[1](i) = s(i,0);
        s_map[2](i) = s(i,1);
        s_map[3](i) = s(i,2);
    }
    
    
    map<pair<Int_t,Int_t>, TMatrixD> D;
    for(auto i = 1; i <= 3; ++i)
        for(auto j = 1; j <= 3; ++j)
            D.insert({{i, j}, TMatrixD(nHits, nHits)});

    for(auto i = 1; i <= 3; ++i)
    {
        for(auto j = 1; j <= 3; ++j)
        {
            TMatrixD &Dij = D[{i,j}];
            const TMatrixD &Cij = C[{i,j}];
            TMatrixD H_tr(TMatrixD::kTransposed, H);
            
            Dij = H * Cij * H_tr;
        }
    }
    
    // ... Computation of weighted sample covariance matrix 𝑨 ...
    map<Int_t, pair<Int_t,Int_t>> nu = {
        {1, {1,1}},
        {2, {1,2}},
        {3, {1,3}},
        {4, {2,2}},
        {5, {2,3}},
        {6, {3,3}}
    };

    map<Int_t, Double_t> A_alpha;
    for(auto alpha = 1; alpha <= 6; ++alpha)
    {
        auto [i, j] = nu[alpha];
        TVectorD w_sj(nHits);
        for(auto k = 0; k < w.size(); ++k)
            w_sj[k] = w_vec(k) * s_map[j](k); // w ⊙ s_j

        A_alpha[alpha] = s_map[i] * w_sj; // s_iᵀ ⋅ (w ⊙ s_j)
    }
    
    
    // Compute covariance matrix of A E_{alpha,beta} for alpha, beta = 1..6
    TMatrixDSym E(6);

    // Precompute W2 = w * w^T (outer product)
    TMatrixD W2(nHits, nHits);
    for(auto i = 0; i < nHits; ++i)
        for(auto j = 0; j < nHits; ++j)
            W2(i,j) = w_vec(i) * w_vec(j);

    for(auto alpha = 1; alpha <= 6; ++alpha)
    {
        auto [i, j] = nu[alpha];
        const auto& si = s_map[i];
        const auto& sj = s_map[j];
        
        for(auto beta = alpha; beta <= 6; ++beta)
        {
            auto [k, l] = nu[beta];
            const auto& sk = s_map[k];
            const auto& sl = s_map[l];

            const TMatrixD& Dik = D[{i,k}];
            const TMatrixD& Dil = D[{i,l}];
            const TMatrixD& Djk = D[{j,k}];
            const TMatrixD& Djl = D[{j,l}];
        
            // Hadamard products
            Double_t S_term = 0.;

            for(auto qi = 0; qi < nHits; ++qi)
                for(auto qj = 0; qj < nHits; ++qj)
                    S_term += (Dik(qi,qj) * W2(qi,qj) * Djl(qi, qj) + Dil(qi,qj) * W2(qi,qj) * Djk(qi, qj));

            TMatrixD temp_jl2(nHits, nHits), temp_jk2(nHits, nHits), temp_il2(nHits, nHits), temp_ik2(nHits, nHits);
            for(auto qu = 0; qu < nHits; ++qu)
                for(auto qv = 0; qv < nHits; ++qv)
                {
                    temp_jl2(qu,qv) = Djl(qu,qv) * W2(qu,qv);
                    temp_jk2(qu,qv) = Djk(qu,qv) * W2(qu,qv);
                    temp_il2(qu,qv) = Dil(qu,qv) * W2(qu,qv);
                    temp_ik2(qu,qv) = Dik(qu,qv) * W2(qu,qv);
                }
            
            Double_t scalar = si * (temp_jl2 * sk) + si * (temp_jk2 * sl) + sj * (temp_il2 * sk) + sj * (temp_ik2 * sl);
            
            // Finally
            E(alpha - 1, beta - 1) = S_term + scalar;
            
            if(alpha != beta)
                E(beta - 1, alpha - 1) = E(alpha - 1, beta - 1); // symmetry
        }
    }
    

    // ... Computation of n and c
    // A matrix
    TMatrixDSym A(3); A.Zero();
    for(auto alpha = 1; alpha <= 6; ++alpha)
    {
        auto [i, j] = nu[alpha];
        A(i-1,j-1) = A_alpha[alpha];
        A(j-1,i-1) = A_alpha[alpha];
    }

    // Diagonalizzazione
    TVectorD eigenVals(3);
    TMatrixD eigenVecs = A.EigenVectors(eigenVals);

    // Normale = autovettore con autovalore minimo
    Int_t minIdx = (eigenVals(0) < eigenVals(1)) ?
                    ((eigenVals(0) < eigenVals(2)) ? 0 : 2) :
                    ((eigenVals(1) < eigenVals(2)) ? 1 : 2);

    TVectorD n(3);
    for(Int_t i = 0; i < 3; ++i)
        n[i] = eigenVecs(i, minIdx);
    
    Double_t c = -(n*r_0);
    
    // Compute the Jacobian
    TMatrixD J2(3, 6);
    const Double_t epsilon = 1e-2;

    for(auto alpha = 1; alpha <= 6; ++alpha)
    {
        // Copia A_alpha e perturba il solo alpha-esimo parametro
        auto A_plus = A_alpha;
        auto A_minus = A_alpha;

        A_plus[alpha] += epsilon;
        A_minus[alpha] -= epsilon;

        auto calc_n = [&](const map<Int_t, Double_t>& A_mod) -> TVectorD
        {
            TMatrixDSym A_mat(3); A_mat.Zero();
            for(auto a = 1; a <= 6; ++a)
            {
                auto [i, j] = nu[a];
                A_mat(i-1,j-1) = A_mod.at(a);
                A_mat(j-1,i-1) = A_mod.at(a);
            }

            TVectorD evals(3);
            TMatrixD evecs = A_mat.EigenVectors(evals);

            Int_t minIdx = (evals(0) < evals(1)) ? ((evals(0) < evals(2)) ? 0 : 2)
                                       : ((evals(1) < evals(2)) ? 1 : 2);

            TVectorD n(3);
            
            for(auto i = 0; i < 3; ++i)
                n[i] = evecs(i, minIdx);

            return n;
        };

        TVectorD n_plus  = calc_n(A_plus);
        TVectorD n_minus = calc_n(A_minus);

        for(auto i = 0; i < 3; ++i)
            J2(i, alpha - 1) = (n_plus[i] - n_minus[i]) / (2. * epsilon);
    }

    //TMatrixDSym C_n(3, 3);
    auto C_n = E.Similarity(J2);      // J2 * E * J2ᵀ


    // Joint covariance matrix of n and c
    // Parte alta-sinistra: Cn
    TMatrixD Cnc(4, 4);
    for(auto i = 0; i < 3; ++i)
        for(auto j = 0; j < 3; ++j)
            Cnc(i, j) = C_n(i, j);

    // Parte in alto a destra: -Cn * r0
    TVectorD Cn_r0(3);
    Cn_r0 = C_n * r_0;

    for(auto i = 0; i < 3; ++i)
    {
        Cnc(i, 3) = -Cn_r0[i];       // colonna finale
        Cnc(3, i) = -Cn_r0[i];       // riga finale (simmetrico)
    }

    // Calcolo var[c]
    Double_t ncnrcr = C_0.Similarity(n) + C_n.Similarity(r_0);
    Double_t S_trace = 0.;
    for(auto i = 0; i < 3; ++i)
        for(auto j = 0; j < 3; ++j)
            S_trace += C_n(i,j) * C_0(i,j); // Hadamard product

    Double_t var_c = ncnrcr + S_trace;
    Cnc(3, 3) = var_c;


    // ... Circle parameters ...
    const Double_t xC = -n(0) / (2 * n(2));
    const Double_t yC = -n(1) / (2 * n(2));
    Double_t r2 = (1 - n(2)*n(2) - 4 * c * n(2)) / (4 * n(2)*n(2));
    const Double_t R  = sqrt(r2);
    
    // Jacobian
    Double_t h = sqrt(1 - n(2)*n(2) - 4*c*n(2));
    TMatrixD J3(3,4);   J3.Zero();
    J3(0,0) = -1./(2*n(2));     J3(0,2) = n(0)/(2*n(2)*n(2));
    J3(1,1) = -1./(2*n(2));     J3(1,2) = n(1)/(2*n(2)*n(2));
    J3(2,2) = -h/(2*n(2)*n(2)) - (4*c + 2*n(2)) / (4*h*n(2));
    J3(2,3) = -1./h;
    

    // Covariance matrix
    TMatrixD J3_tr = TMatrixD(TMatrixD::kTransposed, J3);
    TMatrixD covCircle(3,3);
    covCircle = J3 * Cnc * J3_tr;
    
    if(!opts.processAll)
        cout << Form("\n>>> Prefitter =\nxC = %f +/- %f\nyC = %f +/- %f\nR = %f +/- %f\n\n", xC, sqrt(covCircle(0,0)), yC, sqrt(covCircle(1,1)), R, sqrt(covCircle(2,2)));
    

    
    // --- Fit helix in Z vs arc length s ---
    vector<Double_t> s_values;
    vector<Double_t> z_values;

    // Choose the pivot
    const TVector3 pivot = measuredCoordinates[0];

    // Compute phi0 and dr
    const Double_t phi0 = atan2(pivot.Y() - yC, pivot.X() - xC);
    const Double_t cos_phi0 = cos(phi0);
    const Double_t sin_phi0 = sin(phi0);
    const Double_t dr = sqrt(pow(pivot.X() - xC, 2) + pow(pivot.Y() - yC, 2)) - R;

    // Compute arc lengths s
    Double_t previous_phi = phi0;
    Double_t previous_s = 0;
    for(const auto& point : measuredCoordinates)
    {
        Double_t dx = point.X() - xC;
        Double_t dy = point.Y() - yC;
        Double_t phi = atan2(dy, dx);

        // Angle unwrapping
        Double_t dphi = phi - previous_phi;
        if(dphi > TMath::Pi()) dphi -= 2 * TMath::Pi();
        if(dphi < -TMath::Pi()) dphi += 2 * TMath::Pi();
        dphi *= -1;

        Double_t s = previous_s + R * dphi;

        s_values.push_back(s);
        z_values.push_back(point.Z());

        previous_phi = phi;
        previous_s = s;
    }

    // --- Fit z vs s ---
    auto *graphZvsS = new TGraphErrors(s_values.size());
    for(size_t i = 0; i < s_values.size(); ++i)
    {
        Double_t s = s_values[i];
        Double_t z = z_values[i];

        // Get the hit            
        const auto& hit = measuredCoordinates[i];
        const TMatrixDSym& matCov = hitCov.GetMatrixCartesian(cylinders[i], atan2(hit.Y(), hit.X()));
        
        // Uncertainties on Z
        Double_t sigmaZ = sqrt(matCov(2, 2));

        graphZvsS->SetPoint(i, s, z);
        graphZvsS->SetPointError(i, 0, sigmaZ);
    }

    TF1 *fitZvsS = new TF1("fitZvsS", "[0] + x*[1]", -10, 10);
    fitZvsS->SetParNames("z0", "tan(lambda)");
    fitZvsS->SetParameters(0, 1);
    TFitResultPtr fitlinePtr = nullptr;
    if(!opts.processAll)
        fitlinePtr= graphZvsS->Fit(fitZvsS, "S");
    else
        fitlinePtr= graphZvsS->Fit(fitZvsS, "SQ");

    // Get results
    Double_t z0 = fitZvsS->GetParameter(0);
    Double_t tanLambda = fitZvsS->GetParameter(1);
    
    TMatrixDSym covLine = fitlinePtr->GetCovarianceMatrix();


    // Helix Fit result
    // Construct the full params vector
    fParameters = {xC, yC, R, phi0, z0, tanLambda};

    // Construct the full covariance matrix
    TMatrixDSym covHelix(5);
    for(Int_t i = 0; i < 3; ++i)
        for(Int_t j = 0; j < 3; ++j)
            covHelix(i, j) = covCircle(i, j);
    for(Int_t i = 0; i < 2; ++i)
        for(Int_t j = 0; j < 2; ++j)
            covHelix(i + 3, j + 3) = covLine(i, j);

    return fParameters;
}





//fitTrack->addTrackRep(repHelix);

        //for(auto i = 0; i < nHits; ++i)
        //{
        //    auto tP = fitTrack->getPoint(i);
        //    auto kFI = new genfit::KalmanFitterInfo(tP, rep);
        //    //kFI->setRefenceState();
        //    tP->setFitterInfo(kFI);
        //    cout << tP->getKalmanFitterInfo() << endl;
        //}