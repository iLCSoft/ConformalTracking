#ifndef KDTREE_WRAPPER_HPP
#define KDTREE_WRAPPER_HPP

#include <boost/array.hpp>
#include <boost/multi_array.hpp>
#include <fstream>
#include <functional>
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
typedef std::vector<SKDCluster> VecCluster;

// build the class
class KDTree {
public:
  typedef kdtree2::KDTreeResultVector KDTreeResultVector;

  explicit KDTree(const SharedKDClusters& pts, double overlapTheta, bool sort);
  ~KDTree();

  KDTree(const KDTree&) = delete;
  KDTree& operator=(const KDTree&) = delete;

  void nearestNeighbours(SKDCluster const& pt, int N, SharedKDClusters& result,
                         std::function<bool(SKDCluster const&)> const& filter = [](SKDCluster const&) { return false; });
  void allNeighboursInRadius(SKDCluster const& pt, const double radius, SharedKDClusters& result,
                             std::function<bool(SKDCluster const&)> const& filter = [](SKDCluster const&) { return false; });
  void allNeighboursInTheta(SKDCluster const& pt, const double thetaRange, SharedKDClusters& result,
                            std::function<bool(SKDCluster const&)> const& filter = [](SKDCluster const&) { return false; });
  void allNeighboursInTheta(double theta, const double thetaRange, SharedKDClusters& result,
                            std::function<bool(SKDCluster const&)> const& filter = [](SKDCluster const&) { return false; });

private:
  void transformResults(KDTreeResultVector& vec, SharedKDClusters& result,
                        std::function<bool(SKDCluster const&)> const& filter = [](SKDCluster const&) { return false; });
  void transformThetaResults(KDTreeResultVector& vec, SharedKDClusters& result,
                             std::function<bool(SKDCluster const&)> const& filter = [](SKDCluster const&) { return false; });

  static const int k;
  boost::multi_array<double, 2> array{};
  boost::multi_array<double, 2> arrayTheta{};
  kdtree2::KDTree* tree      = nullptr;
  kdtree2::KDTree* treeTheta = nullptr;
  SharedKDClusters det{};
  std::map<double, SKDCluster> thetaLookup{};
  bool sortTreeResults = true;
};

bool distComparator(const kdtree2::KDTreeResult& a, const kdtree2::KDTreeResult& b);

using UKDTree = std::unique_ptr<KDTree>;

#endif
