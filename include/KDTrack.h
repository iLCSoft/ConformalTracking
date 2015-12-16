#ifndef KDTRACK_H
#define KDTRACK_H 1

#include <math.h>
#include "KDCluster.h"

// ------------------------------------------------------------------------------------
//
// ------------------------------------------------------------------------------------

class KDTrack {
public:
  // Constructor
  KDTrack();

  // Destructor
  virtual ~KDTrack();

  // Functions to add clusters
  void add(KDCluster* cluster) {
    m_clusters.push_back(cluster);
    m_nPoints++;
  }
  void insert(KDCluster* cluster) {
    m_clusters.insert(m_clusters.begin(), cluster);
    m_nPoints++;
  }
  void remove(int clusterN) {
    m_clusters.erase(m_clusters.begin() + clusterN);
    m_nPoints--;
  }

  // Fit functions
  void         fit(double, double);
  const double calculateChi2();

  // Minuit interface functions
  double operator()(const double* x);

  // Functions to set member variables
  void setGradient(double gradient) { m_gradient = gradient; }
  void setIntercept(double intercept) { m_intercept = intercept; }

  // Functions to return member variables
  double                  chi2ndof() { return m_chi2ndof; }
  double                  chi2() { return m_chi2; }
  std::vector<KDCluster*> clusters() { return m_clusters; }
  double                  gradient() { return m_gradient; }
  double                  intercept() { return m_intercept; }
  int                     nPoints() { return m_nPoints; }

  std::vector<KDCluster*> m_clusters;

private:
  // Each KDTrack contains a list of clusters, gradient and intercept
  double m_gradient;
  double m_intercept;
  double m_chi2;
  double m_ndof;
  double m_chi2ndof;
  int    m_nPoints;
};

#endif
