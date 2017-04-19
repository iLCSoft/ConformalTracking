#include "KDTrack.h"

// Constructor
KDTrack::KDTrack(){
  m_gradient = 0.;
  m_intercept = 0.;
  m_nPoints = 0;
  fillFit = false;
  m_rotated = false;
  m_kalmanTrack = NULL;
}

// Destructor
KDTrack::~KDTrack(){
}

// Function to calculate the chi2
double KDTrack::calculateChi2(bool setErrors){
  
  // Value to return
  double chi2 = 0.;
  
  // Loop over all hits on the track and calculate the residual.
  // The y error includes a projection of the x error onto the y axis
  double residual2 = 0;
  for(int hit=0;hit<m_nPoints;hit++){

    double xMeasured = m_clusters[hit]->getU();
    double yMeasured = m_clusters[hit]->getV();
    double dx = m_clusters[hit]->getErrorU();
    double dv = m_clusters[hit]->getErrorV();
    if(m_rotated){
      double newxMeasured = yMeasured;
      double newyMeasured = -1.*xMeasured;
      double newdx = dv;
      double newdv = dx;
      xMeasured=newxMeasured;
      yMeasured=newyMeasured;
      dx=newdx;
      dv=newdv;
    }
    
    double residualY = (m_gradient*xMeasured + m_quadratic*xMeasured*xMeasured + m_intercept) - yMeasured;
    double term = m_gradient+2*m_quadratic*xMeasured;
    double dy2 = (dv*dv) + ( term*term * dx*dx);
    chi2 += (residualY*residualY)/(dy2);
    residual2+=(residualY*residualY);
  }
  
  if(setErrors){
    m_interceptError *= (residual2/(m_nPoints-3));
    m_gradientError *= (residual2/(m_nPoints-3));
    
    m_interceptError = sqrt(m_interceptError);
    m_gradientError = sqrt(m_gradientError);
    
//    std::cout<<"-- Gradient is "<<m_gradient<<" +/- "<<m_gradientError<<std::endl;
//    std::cout<<"-- Intercept is "<<m_intercept<<" +/- "<<m_interceptError<<std::endl;
  }
  
  return chi2;
  
}

// Function to calculate the chi2 of the S-Z fit of the helix
double KDTrack::calculateChi2SZ(TH2F* histo, bool debug){
  
  // Value to return
  double chi2 = 0.;
  
  // Calculate a and b from the conformal fit
  double b = 1./(2.*m_intercept);
  double a = -1.*b*m_gradient;
  
//  double R = sqrt(a*a+b*b);
//  double pt = 0.3 * 4. * R;
  
//  double db = 0;//m_interceptError / (2.*m_intercept*m_intercept);
//  double da = 0;//sqrt( db*db*m_gradient*m_gradient + m_gradientError*m_gradientError*b*b );
  double db = m_interceptError / (2.*m_intercept*m_intercept);
  double da = sqrt( db*db*m_gradient*m_gradient + m_gradientError*m_gradientError*b*b );
  
  // Loop over all hits on the track and calculate the residual.
  // The y error includes a projection of the x error onto the y axis
//  std::cout<<"calculateChi2SZ running "<<std::endl;
  if(fillFit) std::cout<<"== note that a is "<<a<<" +/- "<<da<<" and b is "<<b<<" +/- "<<db<<std::endl;

  double x0 = -1.*m_clusters[0]->getX();
  double y0 = m_clusters[0]->getY();
  double errorx0 = m_clusters[0]->getErrorX();
  double errory0 = m_clusters[0]->getErrorY();;
  if(m_rotated){
    x0 = m_clusters[0]->getY();
    y0 = m_clusters[0]->getX();
    errorx0 = m_clusters[0]->getErrorY();
    errory0 = m_clusters[0]->getErrorX();
  }
  double phi0 = atan2(y0-b,x0-a);
  double phi0offset = 0.;
  
//  if( (phi0 > -M_PI/8. && M_PI < M_PI/8.) ||
//      (phi0 > 3.*M_PI/8. && M_PI < 5.*M_PI/8.) ||
//      (phi0 > -5.*M_PI/8. && M_PI < -3.*M_PI/8.) ||
//     	(phi0 > 7.*M_PI/8. || M_PI < -7.*M_PI/8.)){
//    m_rotatedSZ = true;
//  }
  
  if(m_rotatedSZ){
//    double newx0 = (x0-y0)/sqrt(2.);
//    double newy0 = (y0+x0)/sqrt(2.);
//    double newerrorx0 = (errorx0-errory0)/sqrt(2.);
//    double newerrory0 = (errory0+errorx0)/sqrt(2.);
//    std::cout<<"---- phi0 was: "<<phi0<<", x= "<<x0<<" and y = "<<y0<<". track being rotated"<<std::endl;
//    x0 = newx0;
//    y0 = newy0;
//    errorx0 = newerrorx0;
//    errory0 = newerrory0;
//    phi0 = atan2(y0-b,x0-a);
//    std::cout<<"---- new phi0 is: "<<phi0<<", x = "<<x0<<" = "<<newx0<<" and y = "<<y0<<std::endl;
//    double newa = (a-b)/sqrt(2.);
//    double newb = (b+a)/sqrt(2.);
//    double newda = (da-db)/sqrt(2.);
//    double newdb = (db+da)/sqrt(2.);
//    a = newa; b = newb; da = newda; db = newdb;
    phi0offset = 0.35;
  }
  
  phi0+=phi0offset;
  
  double cPhi0 = cos(phi0);
  double sPhi0 = sin(phi0);
  double ratio = (y0-b)/(x0-a);
  double errorPhi0 = ( 1./((y0-b)*(1+(ratio*ratio))) ) * sqrt( ( ((errorx0*errorx0+da*da)*ratio*ratio*ratio*ratio) + (errory0*errory0+db*db)*ratio*ratio) );

  if(debug) std::cout<<"== Calculating track chi2 in sz"<<std::endl;
  
  for(int hit=1;hit<m_nPoints;hit++){
    
//    double xMa = m_clusters[hit]->getX() - a;
//    double yMb = m_clusters[hit]->getY() - b;
//  
//    double dx = m_clusters[hit]->getErrorX();
//    double dy = m_clusters[hit]->getErrorY();
//  
//    if(m_rotated){
//      xMa = m_clusters[hit]->getY() - a;
//      yMb = (-1.*m_clusters[hit]->getX()) - b;
//      dx = m_clusters[hit]->getErrorY();
//      dy = m_clusters[hit]->getErrorX();
//    }
//    
//    double s = atan2(yMb,xMa);
    
    double xi = -1.*m_clusters[hit]->getX();
    double yi = m_clusters[hit]->getY();
    double errorx = m_clusters[hit]->getErrorX();
    double errory = m_clusters[hit]->getErrorY();
    if(m_rotated){
      xi = m_clusters[hit]->getY();
      yi = m_clusters[hit]->getX();
      errorx = m_clusters[hit]->getErrorY();
      errory = m_clusters[hit]->getErrorX();
    }
//    if(m_rotatedSZ){
//      double newxi = (xi-yi)/sqrt(2.);
//      double newyi = (yi+xi)/sqrt(2.);
//      double newerrorx = (errorx-errory)/sqrt(2.);
//      double newerrory = (errory+errorx)/sqrt(2.);
//      xi = newxi; yi = newyi; errorx = newerrorx; errory = newerrory;
//    }
    
    double deltaPhi = atan2(yi-y0,xi-x0) - phi0offset;
    double s = ((xi-x0)*cPhi0 + (yi-y0)*sPhi0)/(sinc(deltaPhi));
    if(hit == 0) s = 0.;
    
    double dsdx = deltaPhi*cPhi0/sin(deltaPhi);
    double dsdy = deltaPhi*sPhi0/sin(deltaPhi);
    double dsdx0 = deltaPhi*cPhi0/sin(deltaPhi);
    double dsdy0 = deltaPhi*sPhi0/sin(deltaPhi);
    double dsdPhi0 = ((yi-y0)*cPhi0 - (xi-x0)*sPhi0)/(sinc(deltaPhi));
    double dsdDeltaPhi = ((xi-x0)*cPhi0 + (yi-y0)*sPhi0)*(sin(deltaPhi)-deltaPhi*cos(deltaPhi))/(sin(deltaPhi)*sin(deltaPhi));
    
    double newRatio = (yi-y0)/(xi-x0);
    double errorDeltaPhi = ( 1./((yi-y0)*(1+(newRatio*newRatio))) ) * sqrt( ( ((errorx*errorx+errorx0*errorx0)*newRatio*newRatio*newRatio*newRatio) + (errory*errory+errory0*errory0)*newRatio*newRatio) );
    
    double errorS2 = (dsdx*dsdx*errorx*errorx) + (dsdy*dsdy*errory*errory) + (dsdx0*dsdx0*errorx0*errorx0) + (dsdy0*dsdy0*errory0*errory0) + (dsdPhi0*dsdPhi0*errorPhi0*errorPhi0) + (dsdDeltaPhi*dsdDeltaPhi*errorDeltaPhi*errorDeltaPhi);
    if(hit == 0) errorS2 = 0.;
    
    double z = m_clusters[hit]->getZ();
//    double ratio = yMb/xMa;
//    double errorS = ( 1./(yMb*(1+(ratio*ratio))) ) * sqrt( ( (dx*dx+da*da) + (dy*dy+db*db)*ratio*ratio) );
    double errorZ = m_clusters[hit]->getErrorZ();
    
//    if(m_rotatedSZ){
//      z = s;
//      s = m_clusters[hit]->getZ();
//      errorZ = errorS;
//      errorS = m_clusters[hit]->getErrorZ();
//    }

    
//    std::cout<<" - s = "<<s<<", predicted s is "<<(m_gradientZS*m_clusters[hit]->getZ() + m_interceptZS)<<std::endl;
    double residualS = (m_gradientZS*z + m_interceptZS) - s;
//    double residualZ = (m_gradientZS*s + m_interceptZS) - z;
//    double ds = errorS;
//    double dz = errorZ;
    
//    double ds2 = (dz*dz);// + (m_gradientZS*m_gradientZS * ds*ds);
//    double ds2 = (dz*dz + m_gradientZS*m_gradientZS*errorS2);
      double ds2 = (errorS2 + m_gradientZS*m_gradientZS*errorZ*errorZ);
//    double ds2 = (m_gradientZS*m_gradientZS * dz*dz);
    
//    std::cout<<"- residual is "<<residualS<<std::endl;
    chi2 += (residualS*residualS)/(ds2);
//    chi2 += (residualZ*residualZ)/(ds2);
    
    if(debug) std::cout<<"- hit "<<hit<<" has residualS = "<<residualS<<", with error dz = "<<errorZ<<" and error ds = "<<sqrt(errorS2)<<std::endl;
//    if(debug) std::cout<<"- hit "<<hit<<" has residualZ = "<<residualZ<<", with error dz = "<<dz<<" and error ds = "<<sqrt(errorS2)<<std::endl;
    if(fillFit){
      histo->Fill(m_clusters[hit]->getZ(),s);
      std::cout<<"== hit has s = "<<s<<", z = "<<m_clusters[hit]->getZ()<<std::endl;
//      std::cout<<"== ds = "<<errorS<<", error on predicted value is "<<m_gradientZS * dz<<", total error (squared) is "<<ds2<<std::endl;
    }
  }
  
//  std::cout<<"calculateChi2SZ returning "<<chi2<<std::endl;
  return chi2;

  
}

void KDTrack::FillDistribution(TH2F* histo){
  
  fillFit = true;
  calculateChi2SZ(histo);
  fillFit = false;
  
}

// Fit the track (linear regression)

void KDTrack::linearRegression(){
  
  // Decide if this track should be rotated for the fit
  // Check theta of first hit, those close to the y-axis should
  // be rotated.
  if(m_gradient == 0.){
    double theta = m_clusters[0]->getTheta();
    if( (theta > M_PI/4. && theta < 3.*M_PI/4.) || (theta > 5.*M_PI/4. && theta < 7.*M_PI/4.) ) m_rotated = true;
  }
  
//  double vecx[2] = {0., 0.};
  double vecx[3] = {0., 0., 0.}; // quadratic expression
//  double matx[2][2] = {{0., 0.}, {0., 0.}};
  double matx[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}}; // quadratic expression
  
  // Loop over all hits and fill the matrices
  for(int hit=0;hit<m_nPoints;hit++){
    // Get the global point details
    double x = m_clusters[hit]->getU();
    double y =  m_clusters[hit]->getV();
//    double er2 = 1.; // errors are assumed to be the same for all points here
    double er2 = m_clusters[hit]->getErrorV()*m_clusters[hit]->getErrorV(); //FIXME: should this not be errorV??
    
    if(m_rotated){
      double newx = y;
      double newy = -1.*x;
      double newer2 = m_clusters[hit]->getErrorU()*m_clusters[hit]->getErrorU();
      x=newx;
      y=newy;
      er2=newer2;
    }
    
    const double andreEr2 = 1./er2;
    // Fill the matrices
    vecx[0] += y * andreEr2;
    vecx[1] += x * y * andreEr2;
    vecx[2] += x * x * y * andreEr2; // quadratic expression
    matx[0][0] += andreEr2;
    matx[1][0] += x * andreEr2;
    matx[1][1] += x * x * andreEr2;
    matx[2][0] += x * x * andreEr2; // quadratic expression
    matx[2][1] += x * x * x * andreEr2; // quadratic expression
    matx[2][2] += x * x * x * x * andreEr2; // quadratic expression
  }
  
  // Invert the matrix
  
  // Get the determinant
  double detx = matx[0][0] * matx[1][1] - matx[1][0] * matx[1][0];
//  std::cout<<"- determinant is "<<detx<<std::endl;
  double detxq = matx[0][0]*(matx[1][1]*matx[2][2] - matx[2][1]*matx[2][1]) - matx[1][0]*(matx[1][0]*matx[2][2] - matx[2][1]*matx[2][0]) + matx[2][0]*(matx[1][0]*matx[2][1] - matx[1][1]*matx[2][0]) ; // quadratic expression

  // Check for singularities.
  if (detx == 0.) return;

  // Now make the adjoint (for 3*3 matrix inverse)
  double adjx[3][3] = {{0., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};

  adjx[0][0] = (matx[1][1]*matx[2][2] - matx[2][1]*matx[2][1]);
  adjx[1][0] = (matx[1][0]*matx[2][2] - matx[2][0]*matx[2][1])*(-1.);
  adjx[1][1] = (matx[0][0]*matx[2][2] - matx[2][0]*matx[2][0]);
  adjx[2][0] = (matx[1][0]*matx[2][1] - matx[2][0]*matx[1][1]);
  adjx[2][1] = (matx[0][0]*matx[2][1] - matx[2][0]*matx[1][0])*(-1.);
  adjx[2][2] = (matx[0][0]*matx[1][1] - matx[1][0]*matx[1][0]);
  
  // Get the track parameters
  double slopex = (vecx[1] * matx[0][0] - vecx[0] * matx[1][0]) / detx;
  double interceptx = (vecx[0] * matx[1][1] - vecx[1] * matx[1][0]) / detx;
  
  double interceptxq = (vecx[0]*adjx[0][0] + vecx[1]*adjx[1][0] + vecx[2]*adjx[2][0]) / detxq;
  double slopexq =     (vecx[0]*adjx[1][0] + vecx[1]*adjx[1][1] + vecx[2]*adjx[2][1]) / detxq;
  double quadratic =   (vecx[0]*adjx[2][0] + vecx[1]*adjx[2][1] + vecx[2]*adjx[2][2]) / detxq;

//  std::cout<<"- slope from linear relation is: "<<interceptx<<" with intercept "<<slopex<<std::endl;
//  std::cout<<"- slope from quadratic relation is: "<<interceptxq<<" with intercept "<<slopexq<<std::endl;

  // Set the track parameters
  m_gradient = slopexq;
  m_intercept = interceptxq;
  m_quadratic = quadratic;
  
  m_interceptError = adjx[0][0]/detxq; // to be multipled by sigma^2 in chi2 calculation
  m_gradientError = adjx[1][1]/detxq; // to be multipled by sigma^2 in chi2 calculation
  
  // Calculate the chi2
  m_chi2 = this->calculateChi2(true);
  m_chi2ndof = m_chi2/(m_nPoints-3);
//  std::cout<<"- chi2 comes out as "<<m_chi2<<std::endl;
  
}

void KDTrack::linearRegressionConformal(bool debug){
  
  double vecx[2] = {0., 0.};
  double matx[2][2] = {{0., 0.}, {0., 0.}};
  
  // Calculate a and b from the conformal fit
  double b = 1./(2.*m_intercept);
  double a = -1.*b*m_gradient;
  
  double db = m_interceptError / (2.*m_intercept*m_intercept);
  double da = sqrt( db*db*m_gradient*m_gradient + m_gradientError*m_gradientError*b*b );

  double x0 = -1.*m_clusters[0]->getX();
  double y0 = m_clusters[0]->getY();
  if(m_rotated){
    x0 = m_clusters[0]->getY();
    y0 = m_clusters[0]->getX();
  }
  double phi0 = atan2(y0-b,x0-a);
  double cPhi0 = cos(phi0);
  double sPhi0 = sin(phi0);
  if(debug)std::cout<<"== Fitting track in sz"<<std::endl;
  // Loop over all hits and fill the matrices
  for(int hit=1;hit<m_nPoints;hit++){

    // Get the global point details
//    double xMa = m_clusters[hit]->getX() - a;
//    double yMb = m_clusters[hit]->getY() - b;
//    double dx = m_clusters[hit]->getErrorX();
//    double dy = m_clusters[hit]->getErrorY();
//
//    if(m_rotated){
//      xMa = m_clusters[hit]->getY() - a;
//      yMb = (-1.*m_clusters[hit]->getX()) - b;
//      dx = m_clusters[hit]->getErrorY();
//      dy = m_clusters[hit]->getErrorX();
//    }
//
//    double s = atan2(yMb,xMa);

    double xi = -1.*m_clusters[hit]->getX();
    double yi = m_clusters[hit]->getY();
    if(m_rotated){
      xi = m_clusters[hit]->getY();
      yi = m_clusters[hit]->getX();
    }
    
    double deltaPhi = atan2(yi-y0,xi-x0);
    double s = ((xi-x0)*cPhi0 + (yi-y0)*sPhi0)/(sinc(deltaPhi));
    if(hit == 0) s = 0.;
    
    if(debug) std::cout<<"- for hit "<<hit<<" s = "<<s<<std::endl;
    
    double y =  s;
    double x = m_clusters[hit]->getZ();
    double er2 = 1;//m_clusters[hit]->getErrorZ()*m_clusters[hit]->getErrorZ(); // errors are assumed to be the same for all points here
    
//    if(m_rotatedSZ){
//      x = m_clusters[hit]->getZ();
//      y = s;
//      double ratio = yMb/xMa;
//      er2 = ( 1./(yMb*(1+(ratio*ratio))) ) * sqrt( ( (dx*dx+da*da) + (dy*dy+db*db)*ratio*ratio) );
//      er2*=er2;
//    }
    
    const double andreEr2 = 1./er2;

    // Fill the matrices
    vecx[0] += y * andreEr2;
    vecx[1] += x * y * andreEr2;
    matx[0][0] += andreEr2;
    matx[1][0] += x * andreEr2;
    matx[1][1] += x * x * andreEr2;
    
  }
  
  // Invert the matrices
  double detx = matx[0][0] * matx[1][1] - matx[1][0] * matx[1][0];
  
  // Check for singularities.
  if (detx == 0.) return;
  
  // Get the track parameters
  double slopex = (vecx[1] * matx[0][0] - vecx[0] * matx[1][0]) / detx;
  double interceptx = (vecx[0] * matx[1][1] - vecx[1] * matx[1][0]) / detx;
  
  // Set the track parameters
  m_gradientZS = slopex;
  m_interceptZS = interceptx;
  
  // Calculate the chi2
  m_chi2ZS = this->calculateChi2SZ();
  m_chi2ndofZS = m_chi2ZS/(m_nPoints-2);

}

double KDTrack::sinc(double phi){
  return (sin(phi)/phi);
}


