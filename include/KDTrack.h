#ifndef KDTRACK_H
#define KDTRACK_H 1

#include "KDCluster.h"
#include "TH2F.h"
#include <math.h>

// ------------------------------------------------------------------------------------
// The KDTrack class is a simple track class designed to allow fast linear fitting in
// conformal space. The class holds a vector of KDCluster objects (the hits attached to
// the track) and a psuedo-kalman filter object.
// ------------------------------------------------------------------------------------

class KalmanTrack;

class KDTrack{
public:
    
    //--- Constructor and destructor
    KDTrack();
    virtual ~KDTrack();
    
    //--- Functions to add and remove clusters
    void add(KDCluster* cluster){m_clusters.push_back(cluster); m_nPoints++;}
    void insert(KDCluster* cluster){m_clusters.insert(m_clusters.begin(),cluster); m_nPoints++;}
    void remove(int clusterN){m_clusters.erase(m_clusters.begin()+clusterN); m_nPoints--;}
    
    //--- Fit functions
    double calculateChi2();
    double calculateChi2SZ(TH2F* histo = NULL, bool debug = false);
    void linearRegression();
    void linearRegressionConformal(bool debug = false);
    double sinc(double);
    void FillDistribution(TH2F*);
    
    //--- Functions to set and return member variables
    
    // UV fit parameters
    double intercept(){return m_intercept;}
    double gradient(){return m_gradient;}
    double quadratic(){return m_quadratic;}
    double chi2(){return m_chi2;}
    double chi2ndof(){return m_chi2ndof;}
    bool rotated(){return m_rotated;}

    // SZ fit parameters
    double interceptZS(){return m_interceptZS;}
    double gradientZS(){return m_gradientZS;}
	double chi2ZS(){return m_chi2ZS;}
    double chi2ndofZS(){return m_chi2ndofZS;}

    // Clusters and kalman track pointer
    int nPoints(){return m_nPoints;}
    std::vector<KDCluster*> clusters(){return m_clusters;}
    KalmanTrack* kalmanTrack(){return m_kalmanTrack;}
    void setKalmanTrack(KalmanTrack* track){m_kalmanTrack = track;}
    
    //--- Each KDTrack contains parameters for the two separate fits,
    //--- along with the errors and list of clusters used
    
    // UV fit parameters
    double m_gradient;
    double m_gradientError;
    double m_intercept;
    double m_interceptError;
    double m_quadratic;
    double m_chi2;
    double m_chi2ndof;
    bool m_rotated;
    
    // SZ fit parameters
    double m_gradientZS;
    double m_interceptZS;
    double m_chi2ZS;
    double m_chi2ndofZS;
    bool m_rotatedSZ;
    bool fillFit;

    // Clusters and kalman track pointer
    int m_nPoints;
    std::vector<KDCluster*> m_clusters;
    KalmanTrack* m_kalmanTrack;
    
private:
    
};

#endif
