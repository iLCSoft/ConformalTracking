#ifndef KDCLUSTER_H
#define KDCLUSTER_H 1

#include <EVENT/TrackerHitPlane.h>
#include <IMPL/TrackImpl.h>
#include <IMPL/TrackerHitPlaneImpl.h>
#include <math.h>
#include <boost/array.hpp>
#include <boost/multi_array.hpp>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iosfwd>
#include <iostream>
#include <sstream>
#include <vector>

#include <UTIL/BitSet32.h>
#include <UTIL/CellIDEncoder.h>
#include <UTIL/ILDConf.h>
#include <UTIL/LCRelationNavigator.h>

// ------------------------------------------------------------------------------------
// The KDCluster class is a simple hit class used in the Cellular Automaton tracking.
// They are designed to be a lightweight object containing all of the information
// needed for tracking in conformal space, to be ordered in a binary tree (KDtree,
// hence the name). The hits contain their co-ordinates in conformal space, both in
// cartesian (u-v) and polar (r-theta) notation. Cartesian positions allow fast nearest
// neighbour searches, while the polar positions allow both nearest neighbour searches
// in theta (avoiding sector definitions) and directed tracking to flow inside-out or
// outside-in. They additionally contain minimal detector information (id, layer and side)
// ------------------------------------------------------------------------------------

class KDCluster {
public:
  // Constructors, main initialisation is with tracker hit
  KDCluster() {}
  KDCluster(TrackerHitPlane* hit) {
    // Calculate conformal position in cartesian co-ordinates
    double radius2 = (hit->getPosition()[0] * hit->getPosition()[0] + hit->getPosition()[1] * hit->getPosition()[1]);
    m_u            = hit->getPosition()[0] / radius2;
    m_v            = hit->getPosition()[1] / radius2;
    // Note the position in polar co-ordinates
    m_r     = 1. / (sqrt(radius2));
    m_theta = atan2(m_v, m_u) + M_PI;
    // Store the (unaltered) z position
    m_z = hit->getPosition()[2];
  }

  // Destructor
  virtual ~KDCluster() {}

  // Calls to get co-ordinates
  double getU() { return m_u; }
  double getV() { return m_v; }
  double getR() { return m_r; }
  double getTheta() { return m_theta; }
  double getZ() { return m_z; }

  // Manually set co-ordinates
  void setU(double u) { m_u = u; }
  void setV(double v) { m_v = v; }
  void setR(double r) { m_r = r; }
  void setTheta(double theta) { m_theta = theta; }
  void setZ(double z) { m_z = z; }

  // Subdetector information
  void setDetectorInfo(int subdet, int side, int layer) {
    m_subdet = subdet;
    m_side   = side;
    m_layer  = layer;
  }
  int getSubdetector() { return m_subdet; }
  int getSide() { return m_side; }
  int getLayer() { return m_layer; }

  // Check if another hit is on the same detecting layer
  bool sameLayer(KDCluster* kdhit) {
    if (kdhit->getSubdetector() == m_subdet && kdhit->getSide() == m_side && kdhit->getLayer() == m_layer)
      return true;
    return false;
  }

private:
  // Each hit contains the conformal co-ordinates in cartesian
  // and polar notation, plus the subdetector information
  double m_u;
  double m_v;
  double m_r;
  double m_z;
  double m_theta;
  int    m_subdet;
  int    m_side;
  int    m_layer;
};

#endif
