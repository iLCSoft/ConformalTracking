#ifndef KDCLUSTER_H
#define KDCLUSTER_H 1

#include <iosfwd>
#include <vector>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <math.h>
#include <boost/multi_array.hpp>
#include <boost/array.hpp>
#include <EVENT/TrackerHitPlane.h>
#include <IMPL/TrackerHitPlaneImpl.h>
#include <IMPL/TrackImpl.h>

#include <UTIL/CellIDEncoder.h>
#include <UTIL/ILDConf.h>
#include <UTIL/BitSet32.h>
#include <UTIL/LCRelationNavigator.h>

// ------------------------------------------------------------------------------------
// The KDCluster class is a simple hit class used in the Cellular Automaton tracking.
// They are designed to be a lightweight object containing all of the information
// needed for tracking in conformal space, to be ordered in a binary tree (KDtree,
// hence the name). The hits contain their co-ordinates in conformal space, both in
// cartesian (u-v) and polar (r-theta) notation. Cartesian positions allow fast nearest
// neighbour searches, while the polar positions allow both nearest neighbour searches
// in theta (avoiding sector definitions) and directed tracking to flow inside-out or
// outside-in. They additionally contain minimal detector information (id, layer and side)
// ------------------------------------------------------------------------------------

class KDCluster
{
	public:
  
  // Constructors, main initialisation is with tracker hit
  KDCluster() :
    m_x(0.0),
    m_y(0.0),
    m_u(0.0),
    m_v(0.0),
    m_r(0.0),
    m_z(0.0),
    m_s(0.0),
    m_error(0.0),
    m_errorX(0.0),
    m_errorY(0.0),
    m_errorU(0.0),
    m_errorV(0.0),
    m_errorZ(0.0),
    m_errorS(0.0),
    m_theta(0.0),
    m_subdet(0),
    m_side(0),
    m_layer(0),
    m_removed(false),
    m_endcap(false)
  {}
  KDCluster(TrackerHitPlane* hit, bool endcap):
    m_x(hit->getPosition()[0]),
    m_y(hit->getPosition()[1]),
    m_u(0.0),
    m_v(0.0),
    m_r(0.0),
    m_z(hit->getPosition()[2]), // Store the (unaltered) z position
    m_s(0.0),
    m_error(0.0),
    m_errorX(0.0),
    m_errorY(0.0),
    m_errorU(0.0),
    m_errorV(0.0),
    m_errorZ(0.0),
    m_errorS(0.0),
    m_theta(0.0),
    m_subdet(0),
    m_side(0),
    m_layer(0),
    m_removed(false),
    m_endcap(endcap)
  {
      // Calculate conformal position in cartesian co-ordinates
    const double radius2 = (m_x*m_x + m_y*m_y);
      m_u = m_x / radius2 ;
      m_v = m_y / radius2 ;
      // Note the position in polar co-ordinates
			m_r = 1./(sqrt(radius2));
			m_theta = atan2( m_v, m_u ) + M_PI;

      // Get the error in the conformal (uv) plane
      // This is the xy error projected. Unfortunately, the
      // dU is not always aligned with the xy plane, it might
      // be dV. Check and take the smallest
      m_error = hit->getdU(); m_errorZ = hit->getdV();
      if(hit->getdV() < m_error){
        m_error = hit->getdV();
        m_errorZ = hit->getdU();
      }
      
      if(endcap){
        m_errorU = (m_error*fabs(sin(m_theta))+m_errorZ*fabs(cos(m_theta)))*(m_r*m_r);
        m_errorV = (m_error*fabs(cos(m_theta))+m_errorZ*fabs(sin(m_theta)))*(m_r*m_r);
        m_errorX = m_error*sin(m_theta);
        m_errorY = m_error*cos(m_theta);
      }else{
        m_errorU = m_error*m_r*m_r*sin(m_theta);
        m_errorV = m_error*m_r*m_r*cos(m_theta);
        m_errorX = m_error*sin(m_theta);
        m_errorY = m_error*cos(m_theta);
      }
    }
  
  	// Destructor
		virtual ~KDCluster(){}
  
  	// Calls to get co-ordinates
		double getX(){return m_x;}
		double getY(){return m_y;}
		double getU(){return m_u;}
		double getV(){return m_v;}
		double getR(){return m_r;}
	  double getTheta(){return m_theta;}
		double getZ(){return m_z;}
		double getS(){return m_s;}
		double getError(){return m_error;}
		double getErrorX(){return m_errorX;}
		double getErrorY(){return m_errorY;}
		double getErrorU(){return m_errorU;}
		double getErrorV(){return m_errorV;}
		double getErrorZ(){return m_errorZ;}
		double getErrorS(){return m_errorS;}
  	bool removed(){return m_removed;}

  	// Manually set co-ordinates
  	void setU(double u){m_u=u;}
		void setV(double v){m_v=v;}
		void setR(double r){m_r=r;}
		void setTheta(double theta){m_theta=theta;}
		void setZ(double z){m_z=z;}
		void setError(double error){m_error=error;}
		void setErrorS(double errorS){m_errorS=errorS;}
  	void remove(){m_removed = true;}

  	// Subdetector information
		void setDetectorInfo(int subdet,int side, int layer){
			m_subdet=subdet;
			m_side=side;
			m_layer=layer;
		}
		int getSubdetector(){return m_subdet;}
		int getSide(){return m_side;}
		int getLayer(){return m_layer;}
  
  	// Check if another hit is on the same detecting layer
		bool sameLayer(KDCluster* kdhit){
			if(kdhit->getSubdetector() == m_subdet &&
				 kdhit->getSide() == m_side &&
				 kdhit->getLayer() == m_layer) return true;
			return false;
		}
	
	private:

  	// Each hit contains the conformal co-ordinates in cartesian
	  // and polar notation, plus the subdetector information
		double m_x;
		double m_y;
		double m_u;
		double m_v;
		double m_r;
		double m_z;
		double m_s;
		double m_error;
		double m_errorX;
		double m_errorY;
		double m_errorU;
		double m_errorV;
		double m_errorZ;
		double m_errorS;
		double m_theta;
		int m_subdet;
		int m_side;
		int m_layer;
  	bool m_removed;
  	bool m_endcap;

};

// Vector of kd clusters
//typedef std::vector<KDCluster*> KDTrack;


#endif
