#ifndef __KDTREE2_HPP
#define __KDTREE2_HPP

//
// (c) Matthew B. Kennel, Institute for Nonlinear Science, UCSD (2004)
//
// Licensed under the Academic Free License version 1.1 found in file LICENSE
// with additional provisions in that same file.

//
// Implement a kd tree for fast searching of points in a fixed data base
// in k-dimensional Euclidean space.
//
//

#include <stdio.h>  // bug noticed when compiling on Fedora 14, could not find printf
#include <algorithm>
#include <vector>

#include <boost/array.hpp>
#include <boost/multi_array.hpp>
namespace kdtree2 {
  typedef boost::multi_array<double, 2>           KDTreeArray;
  typedef boost::const_multi_array_ref<double, 2> KDTreeROArray;  // read only ref

  typedef struct {
    double lower, upper;
  } interval;

  // let the compiler know that this is a names of classes.
  class KDTreeNode;
  class SearchRecord;

  //
  // struct KDTreeResult
  // class  KDTreeResult
  //

  struct KDTreeResult {
    //
    // the search routines return a (wrapped) vector
    // of these.
    //
  public:
    double dis;  // its square Euclidean distance
    int    idx;  // which neighbor was found
    KDTreeResult(double d, int i) : dis(d), idx(i) {}
  };

  class KDTreeResultVector : public std::vector<KDTreeResult> {
    // inherit a std::vector<KDTreeResult>
    // but, optionally maintain it in heap form as a priority
    // queue.
  public:
    //
    // add one new element to the list of results, and
    // keep it in heap order.  To keep it in ordinary, as inserted,
    // order, then simply use push_back() as inherited
    // via std::vector<>

    void   push_element_and_heapify(KDTreeResult&);
    double replace_maxpri_elt_return_new_maxpri(KDTreeResult&);

    double max_value();
    // return the distance which has the maximum value of all on list,
    // assuming that ALL insertions were made by
    // push_element_and_heapify()
  };

  //
  // class KDTree
  //
  // The main data structure, one for each k-d tree, pointing
  // to a tree of an indeterminate number of "KDTreeNode"s.
  //

  class KDTree {
  public:
    const KDTreeArray& the_data;
    // "the_data" is a reference to the underlying multi_array of the
    // data to be included in the tree.
    //
    // NOTE: this structure does *NOT* own the storage underlying this.
    // Hence, it would be a very bad idea to change the underlying data
    // during use of the search facilities of this tree.
    // Also, the user must deallocate the memory underlying it.

    const int  N;             // number of data points
    int        dim;           //
    bool       sort_results;  // USERS set to 'true'.
    const bool rearrange;     // are we rearranging?

  public:
    //
    // constructor, passing in a multi_array<double,2> , aka
    // KDTreeArray.
    //
    // constructor, has optional 'dim_in' feature, to use only
    // first 'dim_in' components for definition of nearest neighbors.
    //

    KDTree(KDTreeArray& data_in, bool rearrange_in = true, int dim_in = -1);

    KDTree(const KDTree&)            = delete;
    KDTree& operator=(const KDTree&) = delete;

    // destructor
    ~KDTree();

  public:
    // search routines

    void n_nearest_brute_force(std::vector<double>& qv, KDTreeResultVector& result);
    // search for n nearest to a given query vector 'qv' usin
    // exhaustive slow search.  For debugging, usually. first element is the NN

    void n_nearest(std::vector<double>& qv, int nn, KDTreeResultVector& result);
    // search for n nearest to a given query vector 'qv'.

    void n_nearest_around_point(int idxin, int correltime, int nn, KDTreeResultVector& result);
    // search for 'nn' nearest to point [idxin] of the input data, excluding
    // neighbors within correltime

    void r_nearest(std::vector<double>& qv, double r2, KDTreeResultVector& result);
    // search for all neighbors in ball of size (square Euclidean distance)
    // r2.   Return number of neighbors in 'result.size()',

    void r_nearest_around_point(int idxin, int correltime, double r2, KDTreeResultVector& result);
    // like 'r_nearest', but around existing point, with decorrelation
    // interval.

    int r_count(std::vector<double>& qv, double r2);
    // count number of neighbors within square distance r2.
    int r_count_around_point(int idxin, int correltime, double r2);
    // like r_count, c

    friend class KDTreeNode;
    friend class SearchRecord;

  private:
    // private data members

    KDTreeNode* root;  // the root pointer

    const KDTreeArray* data;
    // pointing either to the_data or an internal
    // rearranged data as necessary

    std::vector<int> ind;
    // the index for the tree leaves.  Data in a leaf with bounds [l,u] are
    // in  'the_data[ind[l],*] to the_data[ind[u],*]

    KDTreeArray rearranged_data{};
    // if rearrange is true then this is the rearranged data storage.

    static const int bucketsize = 12;  // global constant.

  private:
    void        set_data(KDTreeArray& din);
    void        build_tree();  // builds the tree.  Used upon construction.
    KDTreeNode* build_tree_for_range(int l, int u, KDTreeNode* parent);
    void        select_on_coordinate(int c, int k, int l, int u);
    int         select_on_coordinate_value(int c, double alpha, int l, int u);
    void        spread_in_coordinate(int c, int l, int u, interval& interv);
  };

  //
  // class KDTREENODE
  //
  // a node in the tree.  Many are created per tree dynamically..
  //

  class KDTreeNode {
  public:
    // constructor
    KDTreeNode(int dim);
    //, int cut_dim_in,
    // 	       double cut_val_in, double cut_val_left_in,
    //	       double cut_val_right_in);
    // destructor
    ~KDTreeNode();

    KDTreeNode(const KDTreeNode&)            = delete;
    KDTreeNode& operator=(const KDTreeNode&) = delete;

  private:
    // visible to self and KDTree.
    friend class KDTree;  // allow kdtree2 to access private

    int    cut_dim = 0;                                             // dimension to cut;
    double cut_val = 0.0, cut_val_left = 0.0, cut_val_right = 0.0;  //cut value
    int    l = 0, u = 0;                                            // extents in index array for searching

    std::vector<interval> box;  // [min,max] of the box enclosing all points

    KDTreeNode *left = nullptr, *right = nullptr;  // pointers to left and right nodes.

    void search(SearchRecord& sr);
    // recursive innermost core routine for searching..

    bool box_in_search_range(SearchRecord& sr);
    // return true if the bounding box for this node is within the
    // search range given by the searchvector and maximum ballsize in 'sr'.

    void check_query_in_bound(SearchRecord& sr);  // debugging only

    // for processing final buckets.
    void process_terminal_node(SearchRecord& sr);
    void process_terminal_node_fixedball(SearchRecord& sr);
  };

}  // namespace kdtree2

#endif
