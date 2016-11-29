#ifndef KDTRACK_H
#define KDTRACK_H 1

#include "KDCluster.h"
#include "TH2F.h"
#include <math.h>

// ------------------------------------------------------------------------------------
// The KDTrack class is a simple track class designed to allow fast linear fitting in
// conformal space. As the errors in the u (global x) direction cannot be ignored in the
// linear track fit, a minimiser function is needed that projects these errors onto the
// v axis. The class holds a vector of KDCluster objects, the hits attached to this track
// ------------------------------------------------------------------------------------

class KalmanTrack;

class KDTrack{
public:
  
  // Constructor
  KDTrack();
  
  // Destructor
  virtual ~KDTrack();
  
  // Functions to add clusters
  void add(KDCluster* cluster){m_clusters.push_back(cluster); m_nPoints++;}
  void insert(KDCluster* cluster){m_clusters.insert(m_clusters.begin(),cluster); m_nPoints++;}
  void remove(int clusterN){m_clusters.erase(m_clusters.begin()+clusterN); m_nPoints--;}
  
		// Fit functions
  void fit();
  double calculateChi2();
  double calculateChi2SZ(TH2F* histo = NULL);
  void linearRegression();
  void linearRegressionConformal();
  
  // Minuit interface functions
  double operator() (const double *x);
  
  // Functions to set member variables
  void setGradient(double gradient){m_gradient = gradient;}
  void setIntercept(double intercept){m_intercept = intercept;}
  void setGradientError(double gradientError){m_gradientError = gradientError;}
  void setInterceptError(double interceptError){m_interceptError = interceptError;}
  void setGradientZS(double gradientZS){m_gradientZS = gradientZS;}
  void setInterceptZS(double interceptZS){m_interceptZS = interceptZS;}
  void setConformalFit(bool fit){m_conformalFit = fit;}
  void setKalmanTrack(KalmanTrack* track){m_kalmanTrack = track;}
  
  // Functions to return member variables
  double chi2ndof(){return m_chi2ndof;}
  double chi2(){return m_chi2;}
  double chi2ndofZS(){return m_chi2ndofZS;}
  double chi2ZS(){return m_chi2ZS;}
  std::vector<KDCluster*> clusters(){return m_clusters;}
  double gradient(){return m_gradient;}
  double intercept(){return m_intercept;}
  int nPoints(){return m_nPoints;}
  KalmanTrack* kalmanTrack(){return m_kalmanTrack;}

  void FillDistribution(TH2F*);
  
  std::vector<KDCluster*> m_clusters;
  bool m_conformalFit;
  KalmanTrack* m_kalmanTrack;

private:
  
  // Each KDTrack contains a list of clusters, gradient and intercept
  double m_gradient;
  double m_intercept;
  double m_quadratic;
  double m_chi2;
  double m_chi2ZS;
  double m_ndof;
  double m_chi2ndof;
  double m_chi2ndofZS;
  int m_nPoints;
  double m_gradientZS;
  double m_interceptZS;
  double m_gradientError;
  double m_interceptError;
//  bool m_conformalFit;
  bool fillFit;
  
};

#endif
