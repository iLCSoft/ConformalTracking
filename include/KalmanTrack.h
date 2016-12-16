#ifndef KALMANTRACK_H
#define KALMANTRACK_H 1

#include "KDCluster.h"
#include "KDTrack.h"

// ------------------------------------------------------------------------------------
// The Kalman track is a simple extension for the KDTrack, which applies a kalman
// filter to the given track + extension, in order to calculate the delta chi2.
// ------------------------------------------------------------------------------------

// Additional node class to store the states at each point in the filter
class KalmanNode{
public:
  
  // Constructor
  KalmanNode(){};
  KalmanNode(KDCluster* measurement, bool rotated){

    if(!rotated){
      m_uMeasured=measurement->getU();
      m_vMeasured=measurement->getV();
      m_errorMeasured=measurement->getErrorV();
    }else{
      m_uMeasured=measurement->getV();
      m_vMeasured=(-1.)*measurement->getU();
      m_errorMeasured=measurement->getErrorU();
    }
    
  };
  
  // Destructor
  ~KalmanNode(){};
  
  // Functions
  
  // Member variables. The gradient at this node after filtering, and
  // the predicted, measured and filtered positions
  double m_gradient;
  double m_quadratic;
  double m_uMeasured;
  double m_vMeasured;
  double m_errorMeasured;
  double m_uPredicted;
  double m_vPredicted;
  double m_errorPredicted;
  double m_uFiltered;
  double m_vFiltered;
  double m_deltaChi2;
  
private:
  
};

// The Kalman track class
class KalmanTrack{
public:
  
  // Constructors
  KalmanTrack();
  KalmanTrack(KDTrack*);
  
  // Destructor
  ~KalmanTrack();
  
  // Functions
  double addCluster(KDCluster*); // returns delta chi2 for adding this cluster
  
  // Member variables
  KDTrack* m_conformalTrack; // The parent conformal track
  std::vector<KDCluster*> m_kalmanClusters; // The additional clusters to which the kalman filter is applied
  std::vector<KalmanNode*> m_nodes;
  double m_moliere; // Multiple scattering angle
  double m_theta; // Polar angle of the parent KDTrack
  bool m_rotated;

private:
  
};


#endif
