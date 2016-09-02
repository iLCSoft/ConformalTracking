#include <algorithm>
#include "KDTree.h"

//-----------------------------------------------------------------------------
// Implementation file for class : KDTree.C
//
// 2011-07-21 : Matthew Reid
//
//-----------------------------------------------------------------------------

// could set the number of neighbous here to protect user changing it, 
// will leave option in for now in nearestNeighbours(..., N,...);
const int KDTree::k(1);

KDTree::KDTree(const VecCluster& pts, double overlapTheta)
	: array(boost::extents[pts.size()][2]), arrayTheta(boost::extents[pts.size()][2])
{
	// Fill multi_array
	det = pts;
	VecCluster::const_iterator iter = pts.begin();
	const VecCluster::const_iterator end = pts.end();
  unsigned long idx(0);
  unsigned long idtheta(0);

	for(; iter != end; ++iter)
	{
		array[idx][0] = (*iter)->getU();
		array[idx][1] = (*iter)->getV();
    ++idx;

    double theta = (*iter)->getTheta();
		arrayTheta[idtheta][0] = theta;
		arrayTheta[idtheta][1] = 0;
    idtheta++;
    thetaLookup[theta] = (*iter);
    
    if(theta > (2.*M_PI-overlapTheta)){
      arrayTheta.resize(boost::extents[arrayTheta.shape()[0]+1][2]);
      arrayTheta[idtheta][0] = theta - 2.*M_PI;
      arrayTheta[idtheta][1] = 0;
      thetaLookup[theta - 2.*M_PI] = (*iter);
      idtheta++;
    }else if(theta < overlapTheta){
      arrayTheta.resize(boost::extents[arrayTheta.shape()[0]+1][2]);
      arrayTheta[idtheta][0] = theta + 2.*M_PI;
      arrayTheta[idtheta][1] = 0;
      thetaLookup[theta + 2.*M_PI] = (*iter);
      idtheta++;
    }
    
	}
  
	// Build kdtree
	tree = new kdtree2::KDTree(array);
	treeTheta = new kdtree2::KDTree(arrayTheta);

}

// d'tor
KDTree::~KDTree()
{
	delete tree;
	delete treeTheta;
	tree = 0;
	treeTheta = 0;
}


bool distComparator(const kdtree2::KDTreeResult& a, const kdtree2::KDTreeResult& b)
{
	return (a.dis < b.dis);
}

void KDTree::nearestNeighbours(KDCluster* pt, int N, VecCluster& result)
{
	// Search kdtree for N points around query pt
	KDTreeResultVector vec;
	std::vector<double> qv(2);
	qv[0] = pt->getU(); //xyplane x position;
	qv[1] = pt->getV(); //xyplane y position;
	tree->n_nearest(qv, N, vec);

	// Sort and transform results
	this->transformResults(vec, result);
}

void KDTree::allNeighboursInRadius(KDCluster* pt, const double radius, VecCluster& result)
{
	// Search kdtree for all points in radius
	KDTreeResultVector vec;
	std::vector<double> qv(2);// should be 2 if using (x,y)-plane
	qv[0] = pt->getU();
	qv[1] = pt->getV();
	//qv[2] = pt.z; // again remove if (x,y) plane only
	tree->r_nearest(qv, radius*radius, vec);

	// Sort and transform results
	this->transformResults(vec, result);
}

void KDTree::allNeighboursInTheta(KDCluster* pt, const double thetaRange, VecCluster& result)
{
	// Search kdtree for all points in radius
	KDTreeResultVector vec;
	std::vector<double> qv(2);// should be 2 if using (x,y)-plane
	qv[0] = pt->getTheta();
	qv[1] = 0.;
	//qv[2] = pt.z; // again remove if (x,y) plane only
	treeTheta->r_nearest(qv, thetaRange*thetaRange, vec);
	
	// Sort and transform results
	this->transformThetaResults(vec, result);
}

void KDTree::allNeighboursInTheta(double theta, const double thetaRange, VecCluster& result)
{
	// Search kdtree for all points in radius
	KDTreeResultVector vec;
	std::vector<double> qv(2);// should be 2 if using (x,y)-plane
	qv[0] = theta;
	qv[1] = 0.;
	//qv[2] = pt.z; // again remove if (x,y) plane only
	treeTheta->r_nearest(qv, thetaRange*thetaRange, vec);
	
	// Sort and transform results
	this->transformThetaResults(vec, result);
}

void KDTree::transformResults(KDTreeResultVector& vec, VecCluster& result)
{
	// Transform results to our VecCluster format
	if(vec.size() > 1) std::sort(vec.begin(), vec.end(), distComparator);
	result.clear();
	result.reserve(vec.size());
	
	KDTreeResultVector::const_iterator iter = vec.begin();
	const KDTreeResultVector::const_iterator end = vec.end();
	const VecCluster::const_iterator end2 = det.end();

	KDCluster* res;

	// Assign back the z value to the NN cluster
	for(; iter != end; ++iter)
	{
		int idx = iter->idx;
		//res.first->globalX(array[idx][0]);
		//res.first->globalY(array[idx][1]);

		//int check(0);
		/*
		for(VecCluster::const_iterator it = det.begin(); it != end2; ++it)
		{
			if (array[idx][0] == (*it)->getU() && array[idx][1] == (*it)->getV() )
			{
				//double z = it->globalZ();
				res = (*it);
				// once assigned no need to continue looping fo 1NN
				if(k==1)break;
			}

		}
		//*/
    res = det[idx];
		result.push_back(res);
	}
}

void KDTree::transformThetaResults(KDTreeResultVector& vec, VecCluster& result)
{
  // Transform results to our VecCluster format
  if(vec.size() > 1) std::sort(vec.begin(), vec.end(), distComparator);
  result.clear();
  result.reserve(vec.size());
  
  KDTreeResultVector::const_iterator iter = vec.begin();
  const KDTreeResultVector::const_iterator end = vec.end();
  const VecCluster::const_iterator end2 = det.end();
  
  KDCluster* res;
  
  // Assign back the z value to the NN cluster
  for(; iter != end; ++iter)
  {
    int idx = iter->idx;
    //res.first->globalX(array[idx][0]);
    //res.first->globalY(array[idx][1]);
    
    //int check(0);
    res = thetaLookup[arrayTheta[idx][0]];
    
//    for(VecCluster::const_iterator it = det.begin(); it != end2; ++it)
//    {
//      if (arrayTheta[idx][0] == (*it)->getTheta())
//      {
//        //double z = it->globalZ();
//        res = (*it);
//        // once assigned no need to continue looping fo 1NN
//        if(k==1)break;
//      }
//      
//    }
    
    result.push_back(res);
  }
}
