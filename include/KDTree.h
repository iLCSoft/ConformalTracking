#ifndef KDTREE_WRAPPER_HPP
#define KDTREE_WRAPPER_HPP

#include <boost/array.hpp>
#include <boost/multi_array.hpp>
#include <fstream>
#include <iomanip>
#include <iosfwd>
#include <iostream>
#include <sstream>
#include <vector>

#include "KDCluster.h"
#include "kdtree2.h"
//#include "Point.h"

//-----------------------------------------------------------------------------
// Class : KDTree
//
// Info : Wrapper interface to the kdtree2 class, to be passed a vector of
// clusters on a given plane. Will return the kNN to a given seed point.
//
// 2011-07-21 : Matthew Reid
//
//-----------------------------------------------------------------------------

// How to use:-
// initialise KDTree nn(<VecCluster> plane_i);
// VecCluster result(N);
// nn2.nearestNeighbours(<TestBeamCluster> seed, N, result);
// Where N is the number of nearest neighbours to be returned

// Define some typedefs for use later
typedef std::vector<KDCluster*> VecCluster;

// build the class
class KDTree {
public:
  typedef kdtree2::KDTreeResultVector KDTreeResultVector;

  explicit KDTree(const VecCluster& pts, double overlapTheta);
  ~KDTree();

  void nearestNeighbours(KDCluster* pt, int N, VecCluster& result);
  void allNeighboursInRadius(KDCluster* pt, const double radius, VecCluster& result);
  void allNeighboursInTheta(KDCluster* pt, const double thetaRange, VecCluster& result);
  void allNeighboursInTheta(double theta, const double thetaRange, VecCluster& result);

private:
  void transformResults(KDTreeResultVector& vec, VecCluster& result);
  void transformThetaResults(KDTreeResultVector& vec, VecCluster& result);

  static const int k;
  boost::multi_array<double, 2> array;
  boost::multi_array<double, 2> arrayTheta;
  kdtree2::KDTree* tree;
  kdtree2::KDTree* treeTheta;
  VecCluster       det;
  std::map<double, KDCluster*> thetaLookup;
};

bool distComparator(const kdtree2::KDTreeResult& a, const kdtree2::KDTreeResult& b);

#endif
