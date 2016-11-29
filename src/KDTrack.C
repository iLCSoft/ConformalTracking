#include "KDTrack.h"

// Constructor
KDTrack::KDTrack() {
  m_gradient     = 0.;
  m_intercept    = 0.;
  m_nPoints      = 0;
  m_conformalFit = true;
  fillFit        = false;
  m_kalmanTrack  = NULL;
}

// Destructor
KDTrack::~KDTrack() {}

// Once the Minuit minimisation has finished, need to assign the best fit parameters
// and calculate the chi2.
void KDTrack::fit() {
  // Calculate the chi2
  m_chi2 = this->calculateChi2();

  //  if(m_chi2 < 10.){
  //    std::cout<<"- track chi2 is "<<m_chi2<<std::endl;
  //    std::cout<<"-- gradientZS is "<<m_gradientZS<<std::endl;
  //    std::cout<<"-- interceptZS is "<<m_interceptZS<<std::endl;
  //    std::cout<<"-- conformal gradient is "<<m_gradient<<std::endl;
  //    std::cout<<"-- conformal intercept is "<<m_intercept<<std::endl;
  //    std::cout<<"-- error on conformal gradient is "<<m_gradientError<<std::endl;
  //    std::cout<<"-- error on conformal intercept is "<<m_interceptError<<std::endl;
  //    std::cout<<"-- track chi2ZS is "<<this->calculateChi2SZ()<<std::endl;
  //  }

  // Set the chi2/ndof
  m_chi2ndof = m_chi2 / (m_nPoints - 2);
}

// Function to calculate the chi2
double KDTrack::calculateChi2() {
  // Value to return
  double chi2 = 0.;

  // Loop over all hits on the track and calculate the residual.
  // The y error includes a projection of the x error onto the y axis
  for (int hit = 0; hit < m_nPoints; hit++) {
    double residualY = (m_gradient * m_clusters[hit]->getU() +
                        m_quadratic * m_clusters[hit]->getU() * m_clusters[hit]->getU() + m_intercept) -
                       m_clusters[hit]->getV();
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
  double chi2 = 0.;

  // Update the track gradient and intercept
  if (m_conformalFit) {
    this->setGradient(parameters[0]);
    this->setIntercept(parameters[1]);
    chi2 += this->calculateChi2();
  } else {
    this->setGradientZS(parameters[0]);
    this->setInterceptZS(parameters[1]);
    chi2 += this->calculateChi2SZ();
  }

  // Return this to minuit
  return chi2;
}

// Function to calculate the chi2 of the S-Z fit of the helix
double KDTrack::calculateChi2SZ(TH2F* histo) {
  // Value to return
  double chi2 = 0.;

  // Calculate a and b from the conformal fit
  double b = 1. / (2. * m_intercept);
  double a = -1. * b * m_gradient;

  double db = 0;  //m_interceptError / (2.*m_intercept*m_intercept);
  double da = 0;  //sqrt( db*db*m_gradient*m_gradient + m_gradientError*m_gradientError*b*b );

  // Loop over all hits on the track and calculate the residual.
  // The y error includes a projection of the x error onto the y axis
  //  std::cout<<"calculateChi2SZ running "<<std::endl;
  if (fillFit)
    std::cout << "== note that a is " << a << " +/- " << da << " and b is " << b << " +/- " << db << std::endl;

  for (int hit = 0; hit < m_nPoints; hit++) {
    double xMa = m_clusters[hit]->getX() - a;
    double yMb = m_clusters[hit]->getY() - b;

    double dx = m_clusters[hit]->getErrorX();
    double dy = m_clusters[hit]->getErrorY();

    double s = atan2(yMb, xMa);
    double errorS =
        (1. / (1 + (yMb / xMa) * (yMb / xMa))) * sqrt((1. / xMa) * (1. / xMa) * (dy * dy + db * db) +
                                                      ((yMb / (xMa * xMa)) * (yMb / (xMa * xMa))) * (dx * dx + da * da));

    //    std::cout<<" - s = "<<s<<", predicted s is "<<(m_gradientZS*m_clusters[hit]->getZ() + m_interceptZS)<<std::endl;
    double residualS = (m_gradientZS * m_clusters[hit]->getZ() + m_interceptZS) - s;
    double ds        = errorS;
    double dz        = m_clusters[hit]->getErrorZ();

    double ds2 = (ds * ds) + (m_gradientZS * m_gradientZS * dz * dz);
    //    double ds2 = (m_gradientZS*m_gradientZS * dz*dz);

    //    std::cout<<"- residual is "<<residualS<<std::endl;
    chi2 += (residualS * residualS) / (ds2);

    if (fillFit) {
      histo->Fill(m_clusters[hit]->getZ(), s);
      std::cout << "== hit has s = " << s << ", z = " << m_clusters[hit]->getZ() << std::endl;
      std::cout << "== ds = " << errorS << ", error on predicted value is " << m_gradientZS * dz
                << ", total error (squared) is " << ds2 << std::endl;
    }
  }

  //  std::cout<<"calculateChi2SZ returning "<<chi2<<std::endl;
  return chi2;
}

void KDTrack::FillDistribution(TH2F* histo) {
  fillFit = true;
  calculateChi2SZ(histo);
  fillFit = false;
}

// Fit the track (linear regression)

void KDTrack::linearRegression() {
  //  double vecx[2] = {0., 0.};
  double vecx[3] = {0., 0., 0.};                                   // quadratic expression
                                                                   //  double matx[2][2] = {{0., 0.}, {0., 0.}};
  double matx[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};  // quadratic expression

  // Loop over all hits and fill the matrices
  for (int hit = 0; hit < m_nPoints; hit++) {
    // Get the global point details
    double x   = m_clusters[hit]->getU();
    double y   = m_clusters[hit]->getV();
    double er2 = 1.;  // errors are assumed to be the same for all points here
    // Fill the matrices
    vecx[0] += y / er2;
    vecx[1] += x * y / er2;
    vecx[2] += x * x * y / er2;  // quadratic expression
    matx[0][0] += 1. / er2;
    matx[1][0] += x / er2;
    matx[1][1] += x * x / er2;
    matx[2][0] += x * x / er2;          // quadratic expression
    matx[2][1] += x * x * x / er2;      // quadratic expression
    matx[2][2] += x * x * x * x / er2;  // quadratic expression
  }

  // Invert the matrix

  // Get the determinant
  double detx = matx[0][0] * matx[1][1] - matx[1][0] * matx[1][0];
  //  std::cout<<"- determinant is "<<detx<<std::endl;
  double detxq = matx[0][0] * (matx[1][1] * matx[2][2] - matx[2][1] * matx[2][1]) -
                 matx[1][0] * (matx[1][0] * matx[2][2] - matx[2][1] * matx[2][0]) +
                 matx[2][0] * (matx[1][0] * matx[2][1] - matx[1][1] * matx[2][0]);  // quadratic expression

  // Check for singularities.
  if (detx == 0.)
    return;

  // Now make the adjoint (for 3*3 matrix inverse)
  double adjx[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};

  adjx[0][0] = (matx[1][1] * matx[2][2] - matx[2][1] * matx[2][1]);
  adjx[1][0] = (matx[1][0] * matx[2][2] - matx[2][0] * matx[2][1]) * (-1.);
  adjx[1][1] = (matx[0][0] * matx[2][2] - matx[2][0] * matx[2][0]);
  adjx[2][0] = (matx[1][0] * matx[2][1] - matx[2][0] * matx[1][1]);
  adjx[2][1] = (matx[0][0] * matx[2][1] - matx[2][0] * matx[1][0]) * (-1.);
  adjx[2][2] = (matx[0][0] * matx[1][1] - matx[1][0] * matx[1][0]);

  // Get the track parameters
  double slopex     = (vecx[1] * matx[0][0] - vecx[0] * matx[1][0]) / detx;
  double interceptx = (vecx[0] * matx[1][1] - vecx[1] * matx[1][0]) / detx;

  double interceptxq = (vecx[0] * adjx[0][0] + vecx[1] * adjx[1][0] + vecx[2] * adjx[2][0]) / detxq;
  double slopexq     = (vecx[0] * adjx[1][0] + vecx[1] * adjx[1][1] + vecx[2] * adjx[2][1]) / detxq;
  double quadratic   = (vecx[0] * adjx[2][0] + vecx[1] * adjx[2][1] + vecx[2] * adjx[2][2]) / detxq;

  //  std::cout<<"- slope from linear relation is: "<<interceptx<<" with intercept "<<slopex<<std::endl;
  //  std::cout<<"- slope from quadratic relation is: "<<interceptxq<<" with intercept "<<slopexq<<std::endl;

  // Set the track parameters
  m_gradient  = slopexq;
  m_intercept = interceptxq;
  m_quadratic = quadratic;

  // Calculate the chi2
  m_chi2     = this->calculateChi2();
  m_chi2ndof = m_chi2 / (m_nPoints - 2);
  //  std::cout<<"- chi2 comes out as "<<m_chi2<<std::endl;
}

void KDTrack::linearRegressionConformal() {
  double vecx[2]    = {0., 0.};
  double matx[2][2] = {{0., 0.}, {0., 0.}};

  // Calculate a and b from the conformal fit
  double b = 1. / (2. * m_intercept);
  double a = -1. * b * m_gradient;

  // Loop over all hits and fill the matrices
  for (int hit = 0; hit < m_nPoints; hit++) {
    // Get the global point details
    double xMa = m_clusters[hit]->getX() - a;
    double yMb = m_clusters[hit]->getY() - b;
    double s   = atan2(yMb, xMa);

    double x   = m_clusters[hit]->getZ();
    double y   = s;
    double er2 = 1.;  // errors are assumed to be the same for all points here

    // Fill the matrices
    vecx[0] += y / er2;
    vecx[1] += x * y / er2;
    matx[0][0] += 1. / er2;
    matx[1][0] += x / er2;
    matx[1][1] += x * x / er2;
  }

  // Invert the matrices
  double detx = matx[0][0] * matx[1][1] - matx[1][0] * matx[1][0];

  // Check for singularities.
  if (detx == 0.)
    return;

  // Get the track parameters
  double slopex     = (vecx[1] * matx[0][0] - vecx[0] * matx[1][0]) / detx;
  double interceptx = (vecx[0] * matx[1][1] - vecx[1] * matx[1][0]) / detx;

  // Set the track parameters
  m_gradientZS  = slopex;
  m_interceptZS = interceptx;

  // Calculate the chi2
  m_chi2ZS     = this->calculateChi2SZ();
  m_chi2ndofZS = m_chi2ZS / (m_nPoints - 2);
}
