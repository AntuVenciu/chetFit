// HelixTrack.h
#ifndef HELIXTRACK_H
#define HELIXTRACK_H

#include <TObject.h>
#include <TGraph.h>
#include <TVector3.h>
#include <vector>

class HelixTrack {
 public:
  double trueMomentum;
  double polarAngle; // angle in x-y frame with respect to x
  double azimuthalAngle; // angle in y-z frame with respect to z
  double spinAngle; // angle with respect to radial axis (ortogonal to spin direction)
  double emissionAngle; // angle with respect to BField: 0 perpendicular to B Field, Pi/2 along it
  TVector3 origin; // origin of the track
  std::vector<std::vector<double>> hitsCoordinates; // coordinate of the hits on detector
  std::vector<std::vector<double>> trackCoordinates; // coordinates for the track (equal spacing)
  std::vector<int> planeID;

  HelixTrack(double momentum, double polar, double azimuthal)
    : trueMomentum(momentum), polarAngle(polar), azimuthalAngle(azimuthal) {
    hitsCoordinates.clear();
    trackCoordinates.clear();
    origin = TVector3(0., 0., 0.);
  }

  void addHit(double x, double y, double z, int i);

  void addTrackPoint(double x, double y, double z);

  void addOrigin(TVector3 track0);

  // Draw the hits in xy view
  TGraph* drawXYView_hits();

  // Draw the hits in yz view
  TGraph* drawYZView_hits();

  // Draw the track in xy view
  TGraph* drawXYView_track();

  // Draw the track in yz view
  TGraph* drawYZView_track();
};


#endif
