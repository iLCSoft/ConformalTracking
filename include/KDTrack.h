#ifndef KDTRACK_H
#define KDTRACK_H 1

#include "KDCluster.h"
#include "Parameters.h"

#include <cmath>

class TH2F;

// ------------------------------------------------------------------------------------
// The KDTrack class is a simple track class designed to allow fast linear fitting in
// conformal space. The class holds a vector of KDCluster objects (the hits attached to
// the track) and a psuedo-kalman filter object.
// ------------------------------------------------------------------------------------

class KalmanTrack;

class KDTrack {
public:
  //--- Constructor and destructor
  KDTrack(Parameters const& par);
  ~KDTrack();

  KDTrack(const KDTrack&) = default;
  KDTrack& operator=(const KDTrack&) = delete;

  //--- Functions to add and remove clusters
  void add(SKDCluster cluster) {
    m_clusters.push_back(cluster);
    m_nPoints++;
  }
  void insert(SKDCluster cluster) {
    m_clusters.insert(m_clusters.begin(), cluster);
    m_nPoints++;
  }
  void remove(int clusterN) {
    m_clusters.erase(m_clusters.begin() + clusterN);
    m_nPoints--;
  }

  //--- Fit functions
  double calculateChi2();
  double calculateChi2SZ(TH2F* histo = NULL, bool debug = false);
  void   linearRegression(bool highPTfit = false);
  void   linearRegressionConformal(bool debug = false);
  double sinc(double);
  void   FillDistribution(TH2F*);

  //--- Functions to set and return member variables

  // UV fit parameters
  double intercept() { return m_intercept; }
  double gradient() { return m_gradient; }
  double quadratic() { return m_quadratic; }
  double chi2() { return m_chi2; }
  double chi2ndof() { return m_chi2ndof; }
  bool   rotated() { return m_rotated; }

  // SZ fit parameters
  double interceptZS() { return m_interceptZS; }
  double gradientZS() { return m_gradientZS; }
  double chi2ZS() { return m_chi2ZS; }
  double chi2ndofZS() { return m_chi2ndofZS; }

  // Clusters and kalman track pointer
  int              nPoints() { return m_nPoints; }
  SharedKDClusters clusters() { return m_clusters; }
  KalmanTrack*     kalmanTrack() { return m_kalmanTrack; }
  void             setKalmanTrack(KalmanTrack* track) { m_kalmanTrack = track; }

  double pt() { return m_pT; }

  //--- Each KDTrack contains parameters for the two separate fits,
  //--- along with the errors and list of clusters used

  // UV fit parameters
  double m_gradient       = 0.0;
  double m_gradientError  = 0.0;
  double m_intercept      = 0.0;
  double m_interceptError = 0.0;
  double m_quadratic      = 0.0;
  double m_chi2           = 0.0;
  double m_chi2ndof       = 0.0;
  bool   m_rotated        = false;

  // SZ fit parameters
  double m_gradientZS  = 0.0;
  double m_interceptZS = 0.0;
  double m_chi2ZS      = 0.0;
  double m_chi2ndofZS  = 0.0;
  bool   m_rotatedSZ   = 0.0;
  bool   fillFit       = false;

  // Clusters and kalman track pointer
  int              m_nPoints = 0;
  double           m_pT      = 0.0;
  SharedKDClusters m_clusters{};
  KalmanTrack*     m_kalmanTrack      = nullptr;
  bool             m_kalmanFitForward = true;

private:
};

typedef std::vector<std::unique_ptr<KDTrack>> UniqueKDTracks;
typedef std::unique_ptr<KDTrack>              UKDTrack;

#endif
