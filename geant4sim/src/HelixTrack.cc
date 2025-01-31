#include "HelixTrack.h"

void HelixTrack::addHit(double x, double y, double z, int i) {
  std::vector<double> hitCoordinates = {x, y, z};
  hitsCoordinates.push_back(hitCoordinates);
  planeID.push_back(i);
}

void HelixTrack::addTrackPoint(double x, double y, double z) {
  std::vector<double> hitCoordinates = {x, y, z};
  trackCoordinates.push_back(hitCoordinates);
}

void HelixTrack::addOrigin(TVector3 track0) {
  origin = track0;
}

TGraph* HelixTrack::drawXYView_hits() {
  int nHits = hitsCoordinates.size();
  double* x = new double[nHits];
  double* y = new double[nHits];
  for (int i = 0; i < nHits; ++i) {
    x[i] = hitsCoordinates[i].at(0);
    y[i] = hitsCoordinates[i].at(1);
  }
  TGraph* gr = new TGraph(nHits, x, y);
  gr->SetMarkerStyle(20);
  return gr;
}

TGraph* HelixTrack::drawYZView_hits() {
  int nHits = hitsCoordinates.size();
  double* y = new double[nHits];
  double* z = new double[nHits];
  for (int i = 0; i < nHits; ++i) {
    y[i] = hitsCoordinates[i].at(1);
    z[i] = hitsCoordinates[i].at(2);
  }
  TGraph* gr = new TGraph(nHits, z, y);
  gr->SetMarkerStyle(20);
  return gr;
}

TGraph* HelixTrack::drawXYView_track() {
  int nHits = trackCoordinates.size();
  double* x = new double[nHits];
  double* y = new double[nHits];
  for (int i = 0; i < nHits; ++i) {
    x[i] = trackCoordinates[i].at(0);
    y[i] = trackCoordinates[i].at(1);
  }
  TGraph* gr = new TGraph(nHits, x, y);
  gr->SetLineColor(kBlue);
  return gr;
}

TGraph* HelixTrack::drawYZView_track() {
  int nHits = trackCoordinates.size();
  double* y = new double[nHits];
  double* z = new double[nHits];
  for (int i = 0; i < nHits; ++i) {
    y[i] = trackCoordinates[i].at(1);
    z[i] = trackCoordinates[i].at(2);
  }
  TGraph* gr = new TGraph(nHits, z, y);
  gr->SetLineColor(kBlue);
  return gr;
}




