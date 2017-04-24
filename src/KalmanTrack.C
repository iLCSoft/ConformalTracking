#include "KalmanTrack.h"

// Constructors
KalmanTrack::KalmanTrack() {}

KalmanTrack::KalmanTrack(KDTrack* seedTrack) {
  // Save the link to this track
  m_conformalTrack = seedTrack;
  m_rotated        = seedTrack->rotated();
  // Save a seed node to start the kalman filter from
  KDCluster*                  seed     = seedTrack->m_clusters[0];
  std::shared_ptr<KalmanNode> seedNode = std::make_shared<KalmanNode>(seed, m_rotated);
  seedNode->m_uFiltered                = seedNode->m_uMeasured;
  seedNode->m_vFiltered                = seedNode->m_uMeasured * seedTrack->gradient() + seedTrack->intercept();
  seedNode->m_gradient                 = seedTrack->gradient();
  seedNode->m_quadratic                = seedTrack->quadratic();
  m_nodes.push_back(seedNode);
  m_theta = seed->getTheta();  // temporary theta?
                               //  m_moliere = 0.; // Scattering angle has to be set explicitly

  // Approximate scattering angle for transverse momentum pT GeV particle with negligible mass (2% X0)..
  m_moliere = 13.6 * sqrt(0.02) * (1. + 0.038 * log(0.02)) / (1000. * seedTrack->m_pT);
}

// Destructor
//KalmanTrack::~KalmanTrack(){
//  for(int n=0;n<m_nodes.size();n++) delete m_nodes[n];
//}

// Function to add a cluster to the track. Returns the delta chi2
double KalmanTrack::addCluster(KDCluster* cluster) {
  //  std::cout<<"adding cluster to kalman filter"<<std::endl;
  // Store the cluster and make a new node
  m_kalmanClusters.push_back(cluster);
  std::shared_ptr<KalmanNode> node = std::make_shared<KalmanNode>(cluster, m_rotated);
  //  std::cout<<"- new node has U,V("<<node->m_uMeasured<<","<<node->m_vMeasured<<")"<<std::endl;
  //  std::cout<<"- previous node has U,V("<<m_nodes.back()->m_uFiltered<<","<<m_nodes.back()->m_vFiltered<<")"<<std::endl;

  // Get the predicted state at the node
  node->m_uPredicted = node->m_uMeasured;

  // Extrapolate the previous node along its gradient
  node->m_vPredicted = m_nodes.back()->m_vFiltered -
                       m_nodes.back()->m_gradient * fabs(node->m_uPredicted - m_nodes.back()->m_uFiltered) +
                       m_nodes.back()->m_quadratic * pow((node->m_uPredicted - m_nodes.back()->m_uFiltered), 2);

  //  node->m_zPredicted = m_nodes.back()->m_zFiltered + m_nodes.back()->m_gradientZS*(node->m_sPredicted-m_nodes.back()->m_sFiltered);

  // Set the error on the predicted position
  node->m_errorPredicted = (node->m_uPredicted - m_nodes.back()->m_uFiltered) * tan(m_moliere);

  //  std::cout<<"- predicted node has U,V("<<node->m_uPredicted<<","<<node->m_vPredicted<<")"<<std::endl;

  // Calculate the delta chi2 for the measured hit
  double deltaChi2  = 0.;
  double du         = node->m_errorPredicted * sin(m_theta);
  double dv         = node->m_errorPredicted * cos(m_theta);
  double dtotal2    = dv * dv + (m_nodes.back()->m_gradient * m_nodes.back()->m_gradient * du * du);
  double deltaV     = (node->m_vPredicted - node->m_vMeasured);
  deltaChi2         = deltaV * deltaV / dtotal2;
  node->m_deltaChi2 = deltaChi2;

  // Make the filtered node which will be propagated in the next step
  // set m_gradient, m_uFiltered and m_vFiltered
  node->m_uFiltered = node->m_uMeasured;
  node->m_vFiltered = (node->m_vMeasured * node->m_errorPredicted + node->m_vPredicted * node->m_errorMeasured) /
                      (node->m_errorMeasured + node->m_errorPredicted);
  node->m_gradient  = (node->m_vFiltered - m_nodes.back()->m_vFiltered) / (node->m_uFiltered - m_nodes.back()->m_uFiltered);
  node->m_quadratic = m_nodes.back()->m_quadratic;
  // Save this node for subsequent steps
  m_nodes.push_back(node);

  // Return the increase in chi2

  // Calculate a and b from the conformal fit
  double intercept = m_conformalTrack->intercept();
  double gradient  = m_conformalTrack->gradient();
  double b         = 1. / (2. * intercept);
  double a         = -1. * b * gradient;

  // Get the errors on a and b
  double db = m_conformalTrack->m_gradientError / (2. * intercept * intercept);
  double da =
      sqrt(db * db * gradient * gradient + m_conformalTrack->m_gradientError * m_conformalTrack->m_gradientError * b * b);

  // Calculate the initial phi0 and its error, used for the calculation of s
  double x0      = -1. * m_conformalTrack->m_clusters[0]->getX();
  double y0      = m_conformalTrack->m_clusters[0]->getY();
  double errorx0 = m_conformalTrack->m_clusters[0]->getErrorX();
  double errory0 = m_conformalTrack->m_clusters[0]->getErrorY();
  ;
  if (m_rotated) {
    x0      = m_conformalTrack->m_clusters[0]->getY();
    y0      = m_conformalTrack->m_clusters[0]->getX();
    errorx0 = m_conformalTrack->m_clusters[0]->getErrorY();
    errory0 = m_conformalTrack->m_clusters[0]->getErrorX();
  }
  double phi0  = atan2(y0 - b, x0 - a);
  double cPhi0 = cos(phi0);
  double sPhi0 = sin(phi0);
  double ratio = (y0 - b) / (x0 - a);
  double errorPhi0 =
      (1. / ((y0 - b) * (1 + (ratio * ratio)))) * sqrt((((errorx0 * errorx0 + da * da) * ratio * ratio * ratio * ratio) +
                                                        (errory0 * errory0 + db * db) * ratio * ratio));

  // Get the global point details
  double xi     = -1. * cluster->getX();
  double yi     = cluster->getY();
  double errorx = cluster->getErrorX();
  double errory = cluster->getErrorY();
  if (m_rotated) {
    xi     = cluster->getY();
    yi     = cluster->getX();
    errorx = cluster->getErrorY();
    errory = cluster->getErrorX();
  }

  // Calculate the delta phi w.r.t the first hit, and then s
  double deltaPhi = atan2(yi - y0, xi - x0);
  double s        = ((xi - x0) * cPhi0 + (yi - y0) * sPhi0) / (m_conformalTrack->sinc(deltaPhi));

  // Calculate the errors on everything
  double dsdx    = deltaPhi * cPhi0 / sin(deltaPhi);
  double dsdy    = deltaPhi * sPhi0 / sin(deltaPhi);
  double dsdx0   = deltaPhi * cPhi0 / sin(deltaPhi);
  double dsdy0   = deltaPhi * sPhi0 / sin(deltaPhi);
  double dsdPhi0 = ((yi - y0) * cPhi0 - (xi - x0) * sPhi0) / (m_conformalTrack->sinc(deltaPhi));
  double dsdDeltaPhi =
      ((xi - x0) * cPhi0 + (yi - y0) * sPhi0) * (sin(deltaPhi) - deltaPhi * cos(deltaPhi)) / (sin(deltaPhi) * sin(deltaPhi));
  double newRatio      = (yi - y0) / (xi - x0);
  double errorDeltaPhi = (1. / ((yi - y0) * (1 + (newRatio * newRatio)))) *
                         sqrt((((errorx * errorx + errorx0 * errorx0) * newRatio * newRatio * newRatio * newRatio) +
                               (errory * errory + errory0 * errory0) * newRatio * newRatio));
  double errorS2 = (dsdx * dsdx * errorx * errorx) + (dsdy * dsdy * errory * errory) + (dsdx0 * dsdx0 * errorx0 * errorx0) +
                   (dsdy0 * dsdy0 * errory0 * errory0) + (dsdPhi0 * dsdPhi0 * errorPhi0 * errorPhi0) +
                   (dsdDeltaPhi * dsdDeltaPhi * errorDeltaPhi * errorDeltaPhi);

  // Now calculate the residual at this point
  double z         = cluster->getZ();
  double residualS = (m_conformalTrack->m_gradientZS * z + m_conformalTrack->m_interceptZS) - s;

  // Calculate the corresponding hit error
  double errorZ = cluster->getErrorZ();
  double ds2    = (errorS2 + m_conformalTrack->m_gradientZS * m_conformalTrack->m_gradientZS * errorZ * errorZ);

  // Increment the chi2
  deltaChi2 += (residualS * residualS) / (ds2);

  //    std::cout<<"Increase in chi2: "<<deltaChi2<<std::endl;

  /*
  double b = 1./(2.*m_conformalTrack->intercept());
  double a = -1.*b*m_conformalTrack->gradient();
 
  double db = 0;
  double da = 0;

  double xMa = cluster->getX() - a;
  double yMb = cluster->getY() - b;
  double dx = cluster->getErrorX();
  double dy = cluster->getErrorY();

  if(m_rotated){
    xMa = cluster->getY() - a;
    yMb = (-1.*cluster->getX()) - b;
    dx = cluster->getErrorY();
    dy = cluster->getErrorX();
  }
  
  double s = atan2(yMb,xMa);
  double ratio = yMb/xMa;

  double errorS = ( 1./(yMb*(1+(ratio*ratio))) ) * sqrt( ( (dx*dx+da*da) + (dy*dy+db*db)*ratio*ratio) );
  double residualZ = (m_conformalTrack->gradientZS()*s + m_conformalTrack->interceptZS()) - cluster->getZ();
  double ds = errorS;
  double dz = cluster->getErrorZ();
  
  double ds2 = (dz*dz) + (m_conformalTrack->gradientZS()*m_conformalTrack->gradientZS() * ds*ds);
  //    double ds2 = (m_gradientZS*m_gradientZS * dz*dz);
  
  //    std::cout<<"- residual is "<<residualS<<std::endl;
  //    chi2 += (residualS*residualS)/(ds2);
  deltaChi2 += (residualZ*residualZ)/(ds2);
//*/

  return deltaChi2;
}
