#include "KalmanTrack.h"

// Constructors
KalmanTrack::KalmanTrack(){
}

KalmanTrack::KalmanTrack(KDTrack* seedTrack){
  // Save the link to this track
  m_conformalTrack = seedTrack;
  // Save a seed node to start the kalman filter from
  KDCluster* seed = seedTrack->m_clusters[0];
  KalmanNode* seedNode = new KalmanNode(seed);
  seedNode->m_uFiltered = seedNode->m_uMeasured;
  seedNode->m_vFiltered = seedNode->m_uMeasured * seedTrack->gradient() + seedTrack->intercept();
  seedNode->m_gradient = seedTrack->gradient();
  m_nodes.push_back(seedNode);
  m_theta = seed->getTheta(); // temporary theta?
//  m_moliere = 0.; // Scattering angle has to be set explicitly
  
  // Approximate scattering angle for 5000 MeV particle with negligible mass..
  m_moliere = 13.6*sqrt(0.01)*(1.+0.038*log(0.01))/(5000.);
  
}

// Destructor
KalmanTrack::~KalmanTrack(){
  for(int n=0;n<m_nodes.size();n++) delete m_nodes[n];
}

// Function to add a cluster to the track. Returns the delta chi2
double KalmanTrack::addCluster(KDCluster* cluster){
  
//  std::cout<<"adding cluster to kalman filter"<<std::endl;
  // Store the cluster and make a new node
  m_kalmanClusters.push_back(cluster);
  KalmanNode* node = new KalmanNode(cluster);
//  std::cout<<"- new node has U,V("<<node->m_uMeasured<<","<<node->m_vMeasured<<")"<<std::endl;
//  std::cout<<"- previous node has U,V("<<m_nodes.back()->m_uFiltered<<","<<m_nodes.back()->m_vFiltered<<")"<<std::endl;
  
  // Get the predicted state at the node
  node->m_uPredicted = node->m_uMeasured;
  
  // Extrapolate the previous node along its gradient
  node->m_vPredicted = m_nodes.back()->m_vFiltered + m_nodes.back()->m_gradient*(node->m_uPredicted-m_nodes.back()->m_uFiltered);
  
  // Set the error on the predicted position
  node->m_errorPredicted = (node->m_uPredicted-m_nodes.back()->m_uFiltered)*tan(m_moliere);
  
//  std::cout<<"- predicted node has U,V("<<node->m_uPredicted<<","<<node->m_vPredicted<<")"<<std::endl;

  // Calculate the delta chi2 for the measured hit
  double deltaChi2=0.;
  double du = node->m_errorPredicted * sin(m_theta);
  double dv = node->m_errorPredicted * cos(m_theta);
  double dtotal2 = dv*dv + (m_nodes.back()->m_gradient*m_nodes.back()->m_gradient * du*du);
  double deltaV = fabs(node->m_vPredicted - node->m_vMeasured);
  deltaChi2 = deltaV*deltaV/dtotal2;
  node->m_deltaChi2 = deltaChi2;
  
  // Make the filtered node which will be propagated in the next step
  // set m_gradient, m_uFiltered and m_vFiltered
  node->m_uFiltered = node->m_uMeasured;
  node->m_vFiltered = (node->m_vMeasured*node->m_errorPredicted + node->m_vPredicted*node->m_errorMeasured)/(node->m_errorMeasured + node->m_errorPredicted);
  node->m_gradient = (node->m_vFiltered - m_nodes.back()->m_vFiltered) / (node->m_uFiltered - m_nodes.back()->m_uFiltered);
  
  // Save this node for subsequent steps
  m_nodes.push_back(node);
  
  // Return the increase in chi2
  return deltaChi2;
}
