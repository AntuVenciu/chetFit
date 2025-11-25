#include "trackdatamanager.hh"

using namespace std;



TrackDataManager::TrackDataManager(Int_t maxEvents)
{
    tracksChain = new TChain("HelixTrackTree");
    for(Int_t i = 0; i < 1; i++)
        tracksChain->Add(Form("../chet_sim_dataset_%d.root", i));

    nEvents = tracksChain->GetEntries();
    cout << "\n>>> There are " << nEvents << " events\n" << endl;

    if(maxEvents != -1)
        nEvents = maxEvents;

    // Init branches
    fOrigin = nullptr;
    hitsCoordinates = nullptr;
    trackCoordinates = nullptr;
    cylinderID = nullptr;

    tracksChain->SetBranchAddress("trueMomentum", &trueMomentum);
    tracksChain->SetBranchAddress("polarAngle", &azimuthalAngle);
    tracksChain->SetBranchAddress("azimuthalAngle", &polarAngle);
    tracksChain->SetBranchAddress("spinAngle", &spinAngle);
    tracksChain->SetBranchAddress("origin", &fOrigin);
    tracksChain->SetBranchAddress("hitsCoordinates", &hitsCoordinates);
    tracksChain->SetBranchAddress("trackCoordinates", &trackCoordinates);
    tracksChain->SetBranchAddress("planeID", &cylinderID);

    // Init histograms
    accPhi = new TEfficiency("accPhi", "Acceptance: Phi vs Momentum; Momentum [MeV/c];#phi [rad]", 10, 0, 68.9, 10, 0., TMath::TwoPi());
    accTheta = new TEfficiency("accTheta", "Acceptance: Theta vs Momentum; Momentum [MeV/c]; #theta [rad]", 10, 0, 68.9, 20, -TMath::Pi(), TMath::Pi());
    effPhi = new TEfficiency("effPhi", "Efficiency: Phi vs Momentum; Momentum [MeV/c];#phi [rad]", 10, 0, 68.9, 10, 0., TMath::TwoPi());
    effTheta = new TEfficiency("effTheta", "Efficiency: Theta vs Momentum; Momentum [MeV/c]; #theta [rad]", 10, 0, 68.9, 20, -TMath::Pi(), TMath::Pi());

    histTurns = new TH1I("histTurns", "Number of Turns;nTurns;Counts", 20, 0, 20);
    effTurns = new TEfficiency("effTurns", "nTurns Efficiency;nTurns;Efficiency", 10, 0, 10);
    histCylinders = new TH1I("histCylinders", "Number of Cylinders;nCylinders;Counts", 7, 0, 7);
    effCylinders = new TEfficiency("effCylinders", "nCylinders Efficiency;nCylinders;Efficiency", 7, 0, 7);
    histCylVMom = new TProfile("histCylVMom", "Number of Cylinders vs Momentum;Momentum [MeV/c];nCylinders", 20, 0, 68.9, 0, 7);
    histTurnsVMom = new TProfile("histTurnsVMom", "Number of Turns vs Momentum;Momentum [MeV/c];nTurns", 20, 0, 68.9, 0, 20);

    graphMom = new TH2D("graphMom", "Linearity: Momentum; Momentum_{MC} [MeV/c]; Momentum_{fit} [MeV/c]", 100, 0., 70., 100, 0., 100.);
    graphTheta = new TH2D("graphTheta", "Linearity: Theta; #theta_{MC} [rad]; #theta_{fit} [rad]", 100, 0., 0., 100, 0., 0.);
    graphPhi = new TH2D("graphPhi", "Linearity: Phi; #phi_{MC} [rad]; #phi_{fit} [rad]", 100, 0., 0., 100, 0., 0.);

    histDiffX = new TH1D("histDiffX", "Pull Plot: decay X position;Pulls;Counts", 50, -10, 10);
    histDiffY = new TH1D("histDiffY", "Pull Plot: decay Y position;Pulls;Counts", 50, -10, 10);
    histDiffZ = new TH1D("histDiffZ", "Pull Plot: decay Z position;Pulls;Counts", 50, -10, 10);
    histDiffMom = new TH1D("histDiffMom", "Pull Plot: Momentum;Pulls;Counts", 50, -10, 10);
    histDiffTheta = new TH1D("histDiffTheta", "Pull Plot: #theta;Pulls;Counts", 50, -10, 10);
    histDiffPhi = new TH1D("histDiffPhi", "Pull Plot: #phi;Pulls;Counts", 50, -10, 10);

    hist2MomRes = new TH2D("hist2MomRes", "Histo 2D: Momentum resolution;Momentum [MeV/c];#sigma_{p} [MeV/c]", 20, 0., 68.9, 50, 0, 20);
    hist2ThetaRes = new TH2D("hist2ThetaRes", "Histo 2D: Theta resolution;#theta [rad];#sigma_{#theta} [rad]", 40, -TMath::Pi(), TMath::Pi(), 50, 0., 0.5);
    hist2PhiRes = new TH2D("hist2PhiRes", "Histo 2D: Phi resolution;#phi [rad];#sigma_{#phi} [rad]", 20, 0., TMath::TwoPi(), 50, 0., 0.5);

    profMomRes = new TProfile("profMomRes", "Profile plot: Momentum resolution;Momentum [MeV/c];#sigma_{p} [MeV/c]", 20, 0., 68.9, 0, 40.);
    profThetaRes = new TProfile("profThetaRes", "Profile plot: Theta resolution;#theta [rad];#sigma_{#theta} [rad]", 40, -TMath::Pi(), TMath::Pi(), 0., 0.5);
    profPhiRes = new TProfile("profPhiRes", "Profile plot: Phi resolution;#phi [rad];#sigma_{#phi} [rad]", 20, 0., TMath::TwoPi(), 0., 0.5);

    histFakeHits = new TH1I("histFakeHits", "Number of pre-fitted hits added", 10, 0, 0);
}



TrackDataManager::~TrackDataManager()
{
    delete tracksChain;
}
