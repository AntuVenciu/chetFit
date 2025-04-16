#include "fitteralgorithms.hh"

using namespace std;
// Units are in cm

constexpr Int_t DEBUG_LVL = 0;
constexpr Double_t CORR_PHIZ = 0.;
constexpr Double_t SCALE_COV = 1.;


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
    exit(0);
}



void FITALG::SpacepointFitter(Options opts)
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
    Double_t theThetaAngle;
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
        // Acceptance and efficiency
    TEfficiency *accPhi = new TEfficiency("accPhi", "Acceptance: Phi vs Momentum; Momentum [MeV/c];#phi [rad]", 10, 0, 68.9, 10, 0., TMath::TwoPi());
    TEfficiency *accTheta = new TEfficiency("accTheta", "Acceptance: Theta vs Momentum; Momentum [MeV/c]; #theta [rad]", 10, 0, 68.9, 20, -TMath::Pi(), TMath::Pi());
    TEfficiency *effPhi = new TEfficiency("effPhi", "Efficiency: Phi vs Momentum; Momentum [MeV/c];#phi [rad]", 10, 0, 68.9, 10, 0., TMath::TwoPi());
    TEfficiency *effTheta = new TEfficiency("effTheta", "Efficiency: Theta vs Momentum; Momentum [MeV/c]; #theta [rad]", 10, 0, 68.9, 20, -TMath::Pi(), TMath::Pi());

        // Turns info
    TH1I *histTurns = new TH1I("histTurns", "Number of Turns;nTurns;Counts", 20, 0, 20);
    TEfficiency *effTurns = new TEfficiency("effTurns","nTurns Efficiency;nTurns;Efficiency", 10, 0, 10);
    TH1I *histCylinders = new TH1I("histCylinders", "Number of Cylinders;nCylinders;Counts", 7, 0, 7);
    TEfficiency *effCylinders = new TEfficiency("effCylinders","nCylinders Efficiency;nCylinders;Efficiency", 7, 0, 7);
    TProfile *histCylVMom = new TProfile("histCylVMom", "Number of Cylinders vs Momentum;Momentum [MeV/c];nCylinders", 20, 0, 68.9, 0, 7);
    TProfile *histTurnsVMom = new TProfile("histTurnsVMom", "Number of Turns vs Momentum;Momentum [MeV/c];nTurns", 20, 0, 68.9, 0, 20);

        // Linearity
    auto *graphMom = new TH2D("graphMom", "Linearity: Momentum; Momentum_{MC} [MeV/c]; Momentum_{fit} [MeV/c]", 100, 0., 0., 100, 0., 100.);
    auto *graphTheta = new TH2D("graphTheta", "Linearity: Theta; #theta_{MC} [rad]; #theta_{fit} [rad]", 100, 0., 0., 100, 0., 0.);
    auto *graphPhi = new TH2D("graphPhi", "Linearity: Phi; #phi_{MC} [rad]; #phi_{fit} [rad]", 100, 0., 0., 100, 0., 0.);

        // Pulls
    TH1D *histDiffX = new TH1D("histDiffX", "Pull Plot: decay X position;Pulls;Counts", 50, -10, 10);
    TH1D *histDiffY = new TH1D("histDiffY", "Pull Plot: decay Y position;Pulls;Counts", 50, -10, 10);
    TH1D *histDiffZ = new TH1D("histDiffZ", "Pull Plot: decay Z position;Pulls;Counts", 50, -10, 10);
    TH1D *histDiffMom = new TH1D("histDiffMom", "Pull Plot: Momentum;Pulls;Counts", 50, -10, 10);
    TH1D *histDiffTheta = new TH1D("histDiffTheta", "Pull Plot: #theta;Pulls;Counts", 50, -10, 10);
    TH1D *histDiffPhi = new TH1D("histDiffPhi", "Pull Plot: #phi;Pulls;Counts", 50, -10, 10);

        // Resolutions
    TH2D *hist2MomRes = new TH2D("hist2MomRes", "Histo 2D: Momentum resolution;Momentum [MeV/c];#sigma_{p} [MeV/c]",  20, 0., 68.9, 50, 0, 20);
    TH2D *hist2ThetaRes = new TH2D("hist2ThetaRes", "Histo 2D: Theta resolution;#theta [rad];#sigma_{#theta} [rad]", 40, -TMath::Pi(), TMath::Pi(), 50, 0., 0.5);
    TH2D *hist2PhiRes = new TH2D("hist2PhiRes", "Histo 2D: Phi resolution;#phi [rad];#sigma_{#phi} [rad]", 20, 0., TMath::TwoPi(), 50, 0., 0.5);

    TProfile *profMomRes = new TProfile("profMomRes", "Profile plot: Momentum resolution;Momentum [MeV/c];#sigma_{p} [MeV/c]", 20, 0., 68.9, 0, 40.);
    TProfile *profThetaRes = new TProfile("profThetaRes", "Profile plot: Theta resolution;#theta [rad];#sigma_{#theta} [rad]", 40, -TMath::Pi(), TMath::Pi(), 0., 0.5);
    TProfile *profPhiRes = new TProfile("profPhiRes", "Profile plot: Phi resolution;#phi [rad];#sigma_{#phi} [rad]", 20, 0., TMath::TwoPi(), 0., 0.5);


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
        tracksChain->GetEntry(ev);
        Int_t nHits = hitsCoordinates->size();

        if(!opts.processAll)
        {
            cout << "\n>>> Event number = " << ev << endl;
            cout << ">>> NHits = " << nHits << endl;
        }

        // Compute the super emission angle theta for future analysis
        Double_t x = fOrigin->X()*1E-1;
        Double_t y = fOrigin->Y()*1E-1;
        Double_t pz = cos(polarAngle);
        Double_t pr = sin(polarAngle) * (x * cos(azimuthalAngle) + y * sin(azimuthalAngle)) / sqrt(x*x + y*y);

        theThetaAngle = TMath::ATan2(pz, pr); // angle in plane (e_r, z) in radiants, in (-pi, pi)

        // Check acceptance
        if(nHits < 3)
        {
            if(!opts.processAll)
                cout << ">>> Track is not in acceptance!" << endl;

            accTheta->Fill(false, trueMomentum, theThetaAngle);
            accPhi->Fill(false, trueMomentum, azimuthalAngle);

            continue;
        }
        inAcceptance++;

        // Sort and if in single event mode draw hits
            // will need a revision when origin point is not in Z = 0 anymore
        auto sortedHits = AUXALG::SortVectorZ(*hitsCoordinates, *cylinderID);
        *hitsCoordinates = sortedHits.first;
        *cylinderID = sortedHits.second;
        nTurns = AUXALG::CountTurns(*hitsCoordinates);
        nCylinders = AUXALG::CountCylinders(*cylinderID);

        if(opts.turnMode)
        {            
            if(!opts.processAll)
            {
                cout << ">>> nTurns = " << nTurns << endl;
                cout << ">>> nCylinders = " << nCylinders << endl;
            }

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

                accTheta->Fill(false, trueMomentum, theThetaAngle);
                accPhi->Fill(false, trueMomentum, azimuthalAngle);

                continue;
            }
        }

        // Track is in acceptance!
        accTheta->Fill(true, trueMomentum, theThetaAngle);
        accPhi->Fill(true, trueMomentum, azimuthalAngle);

        if(!opts.processAll)
        {
            TCanvas *canvHitsXY = new TCanvas("canvHitsXY", "canvHitsXY", 700, 700);
            TCanvas *canvHitsYZ = new TCanvas("canvHitsYZ", "canvHitsYZ", 900, 500);
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
        mom.SetMag(trueMomentum*1E-3);

        if(!opts.processAll)
        {
            cout << Form(">>> Vertex position = (%f, %f, %f) cm", pos[0], pos[1], pos[2]) << endl;
            cout << Form(">>> Vertex momentum = (%f, %f, %f) MeV", mom[0]*1E3, mom[1]*1E3, mom[2]*1E3) << endl;
            cout << Form(">>> Angles (theta, phi) = (%f pi, %f pi) rad", theThetaAngle / TMath::Pi(), azimuthalAngle / TMath::Pi()) << endl;
        }

        // Track candidate
        genfit::TrackCand trackCand;

        // Resolution of detectors
        const Double_t detectorResolution = 0.1;
        const CHeTResolutions hitCov(CORR_PHIZ, SCALE_COV);
    
        //TMatrixDSym hitCov(3);
        //hitCov.UnitMatrix();
        //hitCov *= detectorResolution*detectorResolution;

        // Fill the candidate
        for(Int_t i = 0; i < nHits; i++)
        {
            TVector3 hitCoords;

            vector<Double_t> AbsCoords = hitsCoordinates->at(i);
            //vector<Double_t> measuredCoords = AUXALG::SmearMeasurement((*cylinderID)[i], hitsCoordinates->at(i));
            
            hitCoords[0] = AbsCoords.at(0);
            hitCoords[1] = AbsCoords.at(1);
            hitCoords[2] = AbsCoords.at(2);

            // Smearing?
            //hitCoords[0] = measuredCoords.at(0);
            //hitCoords[1] = measuredCoords.at(1);
            //hitCoords[2] = measuredCoords.at(2);

            new(chetHitArray[i]) genfit::mySpacepointDetectorHit(hitCoords, hitCov.GetMatrixCartesian((*cylinderID)[i], TMath::ATan2(hitCoords[1], hitCoords[0])));
            //new(chetHitArray[i]) genfit::mySpacepointDetectorHit(hitCoords, hitCov);
            trackCand.addHit(detId, i);
        }

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
        fitTrack = new genfit::Track(trackCand, chetFactory, rep);
    
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
            
            TVector3 truePos = (*fOrigin) * 1E-1; // cm
            auto [fRes, fSigma, fPulls] = AUXALG::GetResults(fitTrack, rep, trueMomentum, truePos, theThetaAngle, azimuthalAngle);
            if(fRes.size() == 0) continue;
            histDiffX->Fill(fPulls[0]);
            histDiffY->Fill(fPulls[1]);
            histDiffZ->Fill(fPulls[2]);
            histDiffMom->Fill(fPulls[3]);
            histDiffTheta->Fill(fPulls[4]);
            histDiffPhi->Fill(fPulls[5]);

            graphMom->Fill(trueMomentum, fRes[3]);
            graphTheta->Fill(theThetaAngle, fRes[4]);
            graphPhi->Fill(azimuthalAngle, fRes[5]);

            hist2MomRes->Fill(trueMomentum, fSigma[3]);
            hist2ThetaRes->Fill(theThetaAngle, fSigma[4]);
            hist2PhiRes->Fill(azimuthalAngle, fSigma[5]);

            profMomRes->Fill(trueMomentum, fSigma[3]);
            profThetaRes->Fill(theThetaAngle, fSigma[4]);
            profPhiRes->Fill(azimuthalAngle, fSigma[5]);

            if(!opts.processAll)
                cout << Form(">>> Fitted angles (theta, phi) = (%f pi, %f pi) rad", fRes[4] / TMath::Pi(), fRes[5] / TMath::Pi()) << endl;
        }


        if(!opts.processAll)
        {
            cout << "\n\n>>> Did FIT converge? " << (isFitConverged ? "Yes" : "No") << "\n\n" << endl;
        }

        // Track is in efficiency?
        effTheta->Fill(isFitConverged, trueMomentum, theThetaAngle);
        effPhi->Fill(isFitConverged, trueMomentum, azimuthalAngle);

        // Store turns and cylinders data
        histTurns->Fill(nTurns);
        effTurns->Fill(isFitConverged, nTurns);
        histCylinders->Fill(nCylinders);
        effCylinders->Fill(isFitConverged, nCylinders);
        histTurnsVMom->Fill(trueMomentum, nTurns);
        histCylVMom->Fill(trueMomentum, nCylinders);

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
        cout << ">>> Acceptance = " << (Float_t) (inAcceptance * 100) / nEvents << endl;
        cout << ">>> Efficiency = " << (Float_t) (inEfficiency * 100) / inAcceptance << endl;
        cout << "---------------------------------------------------------\n" << endl;

        // Graphs
        if(opts.saveMode)
            gROOT->SetBatch(true);

        TCanvas *canvEfficiency = new TCanvas("canvEfficiency");
        canvEfficiency->Divide(2,2);
        canvEfficiency->cd(1);
        accTheta->Draw();
        canvEfficiency->cd(2);
        accPhi->Draw();
        canvEfficiency->cd(3);
        effTheta->Draw();
        canvEfficiency->cd(4);
        effPhi->Draw();
    
        TCanvas *canvTurns = new TCanvas("canvTurns");
        canvTurns->Divide(2);
        canvTurns->cd(1);
        histTurns->Draw();
        canvTurns->cd(2);
        effTurns->Draw("AP");

        TCanvas *canvCylinders = new TCanvas("canvCylinders");
        canvCylinders->Divide(2, 2);
        canvCylinders->cd(1);
        histCylinders->Draw();
        canvCylinders->cd(2);
        effCylinders->Draw("AP");
        canvCylinders->cd(3);
        histTurnsVMom->Draw();
        canvCylinders->cd(4);
        histCylVMom->Draw();
        
        TCanvas *canvLinearity = new TCanvas("canvLinearity");
        canvLinearity->Divide(3);
        canvLinearity->cd(1);
        graphMom->SetMarkerStyle(20);
        graphMom->Draw("SCAT");
        canvLinearity->cd(2);
        graphTheta->SetMarkerStyle(20);
        graphTheta->Draw("SCAT");
        canvLinearity->cd(3);
        graphPhi->SetMarkerStyle(20);
        graphPhi->Draw("SCAT");

        TCanvas *canvfPullsPos = new TCanvas("Pulls Position");
        canvfPullsPos->Divide(3);
        
        canvfPullsPos->cd(1);
        histDiffX->Draw();
        
        canvfPullsPos->cd(2);
        histDiffY->Draw();
        
        canvfPullsPos->cd(3);
        histDiffZ->Draw();
        
        TCanvas *canvfPullsMom = new TCanvas("Pulls Momentum");   
        canvfPullsMom->Divide(3);
        
        canvfPullsMom->cd(1);
        histDiffMom->Draw();

        canvfPullsMom->cd(2);
        histDiffTheta->Draw();

        canvfPullsMom->cd(3);
        histDiffPhi->Draw();

        TCanvas *canvProfRes = new TCanvas("ProfResolutions");
        canvProfRes->Divide(3);
        canvProfRes->cd(1);
        profMomRes->Draw();
        canvProfRes->cd(2);
        profThetaRes->Draw();
        canvProfRes->cd(3);
        profPhiRes->Draw();
    
        TCanvas *canvResMom = new TCanvas("canvResMom");
        canvResMom->cd();
        hist2MomRes->Draw();

        TCanvas *canvResTheta = new TCanvas("canvResTheta");
        canvResTheta->cd();
        hist2ThetaRes->Draw();

        TCanvas *canvResPhi = new TCanvas("canvResPhi");
        canvResPhi->cd();
        hist2PhiRes->Draw();

        if(opts.saveMode)
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