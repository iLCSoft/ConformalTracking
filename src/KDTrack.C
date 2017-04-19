#include "KDTrack.h"

// Constructor
KDTrack::KDTrack() {
  m_gradient    = 0.;
  m_intercept   = 0.;
  m_nPoints     = 0;
  fillFit       = false;
  m_rotated     = false;
  m_kalmanTrack = NULL;
}

// Destructor
KDTrack::~KDTrack() {}

// Function to calculate the chi2
double KDTrack::calculateChi2() {
  // Value to return
  double chi2 = 0.;

  // Loop over all hits on the track and calculate their residuals.
  // The y error includes a projection of the x error onto the y axis
  double residual2 = 0;
  for (int hit = 0; hit < m_nPoints; hit++) {
    // Get the u and v for this track. If the track has been
    // rotated for the fit, perform the rotation
    double xMeasured = m_clusters[hit]->getU();
    double yMeasured = m_clusters[hit]->getV();
    double dx        = m_clusters[hit]->getErrorU();
    double dv        = m_clusters[hit]->getErrorV();
    if (m_rotated) {
      double newxMeasured = yMeasured;
      double newyMeasured = -1. * xMeasured;
      double newdx        = dv;
      double newdv        = dx;
      xMeasured           = newxMeasured;
      yMeasured           = newyMeasured;
      dx                  = newdx;
      dv                  = newdv;
    }

    // Get the residual between the track fit and the hit
    double residualY = (m_gradient * xMeasured + m_quadratic * xMeasured * xMeasured + m_intercept) - yMeasured;

    // Get the error on the hit position
    double term = m_gradient + 2 * m_quadratic * xMeasured;
    double dy2  = (dv * dv) + (term * term * dx * dx);

    // Increment the chi2
    chi2 += (residualY * residualY) / (dy2);
    residual2 += (residualY * residualY);
  }

  // Set errors on the gradient and intercept
  m_interceptError *= (residual2 / (m_nPoints - 3));
  m_gradientError *= (residual2 / (m_nPoints - 3));
  m_interceptError = sqrt(m_interceptError);
  m_gradientError  = sqrt(m_gradientError);

  // Return the final chi2
  return chi2;
}

// Function to calculate the chi2 of the S-Z fit of the helix
double KDTrack::calculateChi2SZ(TH2F* histo, bool debug) {
  // Value to return
  double chi2 = 0.;

  // Calculate a and b from the conformal fit
  double b = 1. / (2. * m_intercept);
  double a = -1. * b * m_gradient;

  // Get the errors on a and b
  double db = m_interceptError / (2. * m_intercept * m_intercept);
  double da = sqrt(db * db * m_gradient * m_gradient + m_gradientError * m_gradientError * b * b);

  // Calculate the initial phi0 and its error, used for the calculation of s
  double x0      = -1. * m_clusters[0]->getX();
  double y0      = m_clusters[0]->getY();
  double errorx0 = m_clusters[0]->getErrorX();
  double errory0 = m_clusters[0]->getErrorY();
  ;
  if (m_rotated) {
    x0      = m_clusters[0]->getY();
    y0      = m_clusters[0]->getX();
    errorx0 = m_clusters[0]->getErrorY();
    errory0 = m_clusters[0]->getErrorX();
  }
  double phi0  = atan2(y0 - b, x0 - a);
  double cPhi0 = cos(phi0);
  double sPhi0 = sin(phi0);
  double ratio = (y0 - b) / (x0 - a);
  double errorPhi0 =
      (1. / ((y0 - b) * (1 + (ratio * ratio)))) * sqrt((((errorx0 * errorx0 + da * da) * ratio * ratio * ratio * ratio) +
                                                        (errory0 * errory0 + db * db) * ratio * ratio));

  // Loop over all hits on the track and calculate the residual.
  // The y error includes a projection of the x error onto the y axis
  if (debug)
    std::cout << "== Calculating track chi2 in sz" << std::endl;
  if (fillFit)
    std::cout << "== note that a is " << a << " +/- " << da << " and b is " << b << " +/- " << db << std::endl;
  for (int hit = 1; hit < m_nPoints; hit++) {
    // Get the global point details
    double xi     = -1. * m_clusters[hit]->getX();
    double yi     = m_clusters[hit]->getY();
    double errorx = m_clusters[hit]->getErrorX();
    double errory = m_clusters[hit]->getErrorY();
    if (m_rotated) {
      xi     = m_clusters[hit]->getY();
      yi     = m_clusters[hit]->getX();
      errorx = m_clusters[hit]->getErrorY();
      errory = m_clusters[hit]->getErrorX();
    }

    // Calculate the delta phi w.r.t the first hit, and then s
    double deltaPhi = atan2(yi - y0, xi - x0) - phi0offset;
    double s        = ((xi - x0) * cPhi0 + (yi - y0) * sPhi0) / (sinc(deltaPhi));

    // Calculate the errors on everything
    double dsdx        = deltaPhi * cPhi0 / sin(deltaPhi);
    double dsdy        = deltaPhi * sPhi0 / sin(deltaPhi);
    double dsdx0       = deltaPhi * cPhi0 / sin(deltaPhi);
    double dsdy0       = deltaPhi * sPhi0 / sin(deltaPhi);
    double dsdPhi0     = ((yi - y0) * cPhi0 - (xi - x0) * sPhi0) / (sinc(deltaPhi));
    double dsdDeltaPhi = ((xi - x0) * cPhi0 + (yi - y0) * sPhi0) * (sin(deltaPhi) - deltaPhi * cos(deltaPhi)) /
                         (sin(deltaPhi) * sin(deltaPhi));
    double newRatio      = (yi - y0) / (xi - x0);
    double errorDeltaPhi = (1. / ((yi - y0) * (1 + (newRatio * newRatio)))) *
                           sqrt((((errorx * errorx + errorx0 * errorx0) * newRatio * newRatio * newRatio * newRatio) +
                                 (errory * errory + errory0 * errory0) * newRatio * newRatio));
    double errorS2 = (dsdx * dsdx * errorx * errorx) + (dsdy * dsdy * errory * errory) +
                     (dsdx0 * dsdx0 * errorx0 * errorx0) + (dsdy0 * dsdy0 * errory0 * errory0) +
                     (dsdPhi0 * dsdPhi0 * errorPhi0 * errorPhi0) +
                     (dsdDeltaPhi * dsdDeltaPhi * errorDeltaPhi * errorDeltaPhi);

    // Now calculate the residual at this point
    double z         = m_clusters[hit]->getZ();
    double residualS = (m_gradientZS * z + m_interceptZS) - s;

    // Calculate the corresponding hit error
    double errorZ = m_clusters[hit]->getErrorZ();
    double ds2    = (errorS2 + m_gradientZS * m_gradientZS * errorZ * errorZ);

    // Increment the chi2
    chi2 += (residualS * residualS) / (ds2);

    // Debug info
    if (debug)
      std::cout << "- hit " << hit << " has residualS = " << residualS << ", with error dz = " << errorZ
                << " and error ds = " << sqrt(errorS2) << std::endl;
    if (fillFit) {
      histo->Fill(m_clusters[hit]->getZ(), s);
      std::cout << "== hit has s = " << s << ", z = " << m_clusters[hit]->getZ() << std::endl;
    }
  }

  // Return the final summed chi2
  return chi2;
}

// Debug function to fill sz histogram
void KDTrack::FillDistribution(TH2F* histo) {
  fillFit = true;
  calculateChi2SZ(histo);
  fillFit = false;
}

// Fit the track in uv (multiple linear regression)
void KDTrack::linearRegression() {
  // Decide if this track should be rotated for the fit (fits fail where
  // the track points along the y-axis. Check the theta of the first hit,
  // those close to the y-axis should be rotated.

  // If track has not yet been fitted
  if (m_gradient == 0.) {
    // Check for hit within 90 degrees of +/- y-axis
    double theta = m_clusters[0]->getTheta();
    if ((theta > M_PI / 4. && theta < 3. * M_PI / 4.) || (theta > 5. * M_PI / 4. && theta < 7. * M_PI / 4.))
      m_rotated = true;
  }

  // Initialise variables for the linear regression
  double vecx[3]    = {0., 0., 0.};
  double matx[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};

  // Loop over all hits and fill the matrices
  for (int hit = 0; hit < m_nPoints; hit++) {
    // Get the global point details
    double x   = m_clusters[hit]->getU();
    double y   = m_clusters[hit]->getV();
    double er2 = m_clusters[hit]->getErrorV() * m_clusters[hit]->getErrorV();

    // If rotated then perform the rotation
    if (m_rotated) {
      double newx   = y;
      double newy   = -1. * x;
      double newer2 = m_clusters[hit]->getErrorU() * m_clusters[hit]->getErrorU();
      x             = newx;
      y             = newy;
      er2           = newer2;
    }

    // Fill the matrices
    const double inverseEr2 = 1. / er2;
    vecx[0] += y * inverseEr2;
    vecx[1] += x * y * inverseEr2;
    vecx[2] += x * x * y * inverseEr2;  // quadratic expression
    matx[0][0] += inverseEr2;
    matx[1][0] += x * inverseEr2;
    matx[1][1] += x * x * inverseEr2;
    matx[2][0] += x * x * inverseEr2;          // quadratic expression
    matx[2][1] += x * x * x * inverseEr2;      // quadratic expression
    matx[2][2] += x * x * x * x * inverseEr2;  // quadratic expression
  }

  // Invert the matrix

  // Get the determinant
  double detx = matx[0][0] * (matx[1][1] * matx[2][2] - matx[2][1] * matx[2][1]) -
                matx[1][0] * (matx[1][0] * matx[2][2] - matx[2][1] * matx[2][0]) +
                matx[2][0] * (matx[1][0] * matx[2][1] - matx[1][1] * matx[2][0]);

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
  double intercept = (vecx[0] * adjx[0][0] + vecx[1] * adjx[1][0] + vecx[2] * adjx[2][0]) / detx;
  double gradient  = (vecx[0] * adjx[1][0] + vecx[1] * adjx[1][1] + vecx[2] * adjx[2][1]) / detx;
  double quadratic = (vecx[0] * adjx[2][0] + vecx[1] * adjx[2][1] + vecx[2] * adjx[2][2]) / detx;

  // Set the track parameters
  m_intercept = intercept;
  m_gradient  = gradient;
  m_quadratic = quadratic;

  // Set the corresponding errors
  m_interceptError = adjx[0][0] / detx;  // to be multipled by sigma^2 in chi2 calculation
  m_gradientError  = adjx[1][1] / detx;  // to be multipled by sigma^2 in chi2 calculation

  // Calculate the chi2
  m_chi2     = this->calculateChi2();
  m_chi2ndof = m_chi2 / (m_nPoints - 3);
}

// Fit the track in sz (linear regression)
void KDTrack::linearRegressionConformal(bool debug) {
  // Initialise variables for the linear regression
  double vecx[2]    = {0., 0.};
  double matx[2][2] = {{0., 0.}, {0., 0.}};

  // Calculate a and b from the conformal fit
  double b = 1. / (2. * m_intercept);
  double a = -1. * b * m_gradient;

  // Get the errors on a and b
  double db = m_interceptError / (2. * m_intercept * m_intercept);
  double da = sqrt(db * db * m_gradient * m_gradient + m_gradientError * m_gradientError * b * b);

  // Calculate the initial phi0, used for the calculation of s
  double x0 = -1. * m_clusters[0]->getX();
  double y0 = m_clusters[0]->getY();
  if (m_rotated) {
    x0 = m_clusters[0]->getY();
    y0 = m_clusters[0]->getX();
  }
  double phi0  = atan2(y0 - b, x0 - a);
  double cPhi0 = cos(phi0);
  double sPhi0 = sin(phi0);

  // Loop over all hits and fill the matrices
  if (debug)
    std::cout << "== Fitting track in sz" << std::endl;
  for (int hit = 1; hit < m_nPoints; hit++) {
    // Get the global point details
    double xi = -1. * m_clusters[hit]->getX();
    double yi = m_clusters[hit]->getY();
    if (m_rotated) {
      xi = m_clusters[hit]->getY();
      yi = m_clusters[hit]->getX();
    }

    // Calculate the delta phi w.r.t the first hit, and then s
    double deltaPhi = atan2(yi - y0, xi - x0);
    double s        = ((xi - x0) * cPhi0 + (yi - y0) * sPhi0) / (sinc(deltaPhi));
    if (debug)
      std::cout << "- for hit " << hit << " s = " << s << std::endl;

    // Now set the values for the fit
    double       y          = s;
    double       x          = m_clusters[hit]->getZ();
    double       er2        = 1;  // errors are assumed to be the same for all points here
    const double inverseEr2 = 1. / er2;

    // Fill the matrices
    vecx[0] += y * inverseEr2;
    vecx[1] += x * y * inverseEr2;
    matx[0][0] += inverseEr2;
    matx[1][0] += x * inverseEr2;
    matx[1][1] += x * x * inverseEr2;
  }

  // Invert the matrices
  double detx = matx[0][0] * matx[1][1] - matx[1][0] * matx[1][0];

  // Check for singularities.
  if (detx == 0.)
    return;

  // Get the track parameters
  double slope     = (vecx[1] * matx[0][0] - vecx[0] * matx[1][0]) / detx;
  double intercept = (vecx[0] * matx[1][1] - vecx[1] * matx[1][0]) / detx;

  // Set the track parameters
  m_gradientZS  = slope;
  m_interceptZS = intercept;

  // Calculate the chi2
  m_chi2ZS     = this->calculateChi2SZ();
  m_chi2ndofZS = m_chi2ZS / (m_nPoints - 2);
}

// Mathematical sinc function
double KDTrack::sinc(double phi) { return (sin(phi) / phi); }
