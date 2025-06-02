#ifndef KALMANTRACK_H
#define KALMANTRACK_H 1

#include "KDCluster.h"
#include "KDTrack.h"

#include <memory>

// ------------------------------------------------------------------------------------
// The Kalman track is a simple extension for the KDTrack, which applies a kalman
// filter to the given track + extension, in order to calculate the delta chi2.
// ------------------------------------------------------------------------------------

// Additional node class to store the states at each point in the filter
class KalmanNode {
public:
  // Constructor
  KalmanNode() = default;
  KalmanNode(SKDCluster measurement, bool rotated) {
    if (!rotated) {
      m_uMeasured     = measurement->getU();
      m_vMeasured     = measurement->getV();
      m_errorMeasured = measurement->getErrorV();
    } else {
      m_uMeasured     = measurement->getV();
      m_vMeasured     = (-1.) * measurement->getU();
      m_errorMeasured = measurement->getErrorU();
    }
  };

  // Destructor
  ~KalmanNode() = default;

  // Functions

  // Member variables. The gradient at this node after filtering, and
  // the predicted, measured and filtered positions
  double m_gradient       = 0.0;
  double m_quadratic      = 0.0;
  double m_uMeasured      = 0.0;
  double m_vMeasured      = 0.0;
  double m_errorMeasured  = 0.0;
  double m_uPredicted     = 0.0;
  double m_vPredicted     = 0.0;
  double m_errorPredicted = 0.0;
  double m_uFiltered      = 0.0;
  double m_vFiltered      = 0.0;

  double m_gradientZS       = 0.0;
  double m_sMeasured        = 0.0;
  double m_zMeasured        = 0.0;
  double m_errorMeasuredZS  = 0.0;
  double m_sPredicted       = 0.0;
  double m_zPredicted       = 0.0;
  double m_errorPredictedZS = 0.0;
  double m_sFiltered        = 0.0;
  double m_zFiltered        = 0.0;

  double m_deltaChi2 = 0.0;

private:
};

// The Kalman track class
class KalmanTrack {
public:
  // Constructors
  KalmanTrack();
  KalmanTrack(KDTrack*);

  KalmanTrack(const KalmanTrack&)            = default;
  KalmanTrack& operator=(const KalmanTrack&) = default;

  // Destructor
  ~KalmanTrack() = default;

  // Functions
  double addCluster(SKDCluster);  // returns delta chi2 for adding this cluster

  // Member variables
  KDTrack*         m_conformalTrack = nullptr;  // The parent conformal track
  SharedKDClusters m_kalmanClusters{};          // The additional clusters to which the kalman filter is applied
  std::vector<std::shared_ptr<KalmanNode>> m_nodes{};
  double                                   m_moliere = 0.0;  // Multiple scattering angle
  double                                   m_theta   = 0.0;  // Polar angle of the parent KDTrack
  bool                                     m_rotated = false;

private:
};

#endif
