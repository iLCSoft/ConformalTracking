#include "KDTrack.h"

// Constructor
KDTrack::KDTrack() {
  m_gradient  = 0.;
  m_intercept = 0.;
  m_nPoints   = 0;
}

// Destructor
KDTrack::~KDTrack() {}

// Once the Minuit minimisation has finished, need to assign the best fit parameters
// and calculate the chi2.
void KDTrack::fit(double gradient, double intercept) {
  // Set the gradient and intercept
  m_gradient  = gradient;
  m_intercept = intercept;

  // Calculate the chi2
  m_chi2 = this->calculateChi2();

  // Set the chi2/ndof
  m_chi2ndof = m_chi2 / (m_nPoints - 2);
}

// Function to calculate the chi2
const double KDTrack::calculateChi2() {
  // Value to return
  double chi2 = 0.;

  // Loop over all hits on the track and calculate the residual.
  // The y error includes a projection of the x error onto the y axis
  for (int hit = 0; hit < m_nPoints; hit++) {
    double residualY = (m_gradient * m_clusters[hit]->getU() + m_intercept) - m_clusters[hit]->getV();
    //    double du = m_clusters[hit]->getError()*m_clusters[hit]->getR()*m_clusters[hit]->getR()*sin(m_clusters[hit]->getTheta());
    //    double dv = m_clusters[hit]->getError()*m_clusters[hit]->getR()*m_clusters[hit]->getR()*cos(m_clusters[hit]->getTheta());
    //    double dx = du;
    //    double dx = m_clusters[hit]->getErrorU();

    double dx  = m_clusters[hit]->getErrorU();
    double dv  = m_clusters[hit]->getErrorV();
    double dy2 = (dv * dv) + (m_gradient * m_gradient * dx * dx);
    chi2 += (residualY * residualY) / (dy2);
  }

  return chi2;
}

// Minimisation operator used by Minuit. Minuit passes the current iteration of
// the parameters and checks if the chi2 is better or worse
double KDTrack::operator()(const double* parameters) {
  //  std::cout<<" --setting gradient "<<parameters[0]<<" and intercept "<<parameters[1]<<std::endl;
  // Update the track gradient and intercept
  this->setGradient(parameters[0]);
  this->setIntercept(parameters[1]);

  // Calculate the chi2
  const double chi2 = this->calculateChi2();

  // Return this to minuit
  return chi2;
}