/*
 * This file belongs to the Galois project, a C++ library for exploiting
 * parallelism. The code is being released under the terms of the 3-Clause BSD
 * License (a copy is located in LICENSE.txt at the top-level directory).
 *
 * Copyright (C) 2018, The University of Texas at Austin. All rights reserved.
 * UNIVERSITY EXPRESSLY DISCLAIMS ANY AND ALL WARRANTIES CONCERNING THIS
 * SOFTWARE AND DOCUMENTATION, INCLUDING ANY WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR ANY PARTICULAR PURPOSE, NON-INFRINGEMENT AND WARRANTIES OF
 * PERFORMANCE, AND ANY WARRANTY THAT MIGHT OTHERWISE ARISE FROM COURSE OF
 * DEALING OR USAGE OF TRADE.  NO WARRANTY IS EITHER EXPRESS OR IMPLIED WITH
 * RESPECT TO THE USE OF THE SOFTWARE OR DOCUMENTATION. Under no circumstances
 * shall University be liable for incidental, special, indirect, direct or
 * consequential damages or loss of profits, interruption of business, or
 * related expenses which may arise from use of Software or Documentation,
 * including but not limited to those resulting from defects in Software and/or
 * Documentation, or loss or inaccuracy of data of any kind.
 */

/**
 * @file LC_CSR_CSC_Graph.h
 *
 * Contains the implementation of a bidirectional LC_CS_Graph.
 */
#ifndef GALOIS_GRAPHS_LC_CSR_CSC_GRAPH_H
#define GALOIS_GRAPHS_LC_CSR_CSC_GRAPH_H

#include "galois/config.h"

#include "galois/graphs/LC_CSR_Graph.h"

namespace galois {
namespace graphs {

/**
 * An bidirectional LC_CSR_Graph that allows the construction of in-edges from
 * its outedges.
 *
 * @tparam NodeTy type of the node data
 * @tparam EdgeTy type of the edge data
 * @tparam EdgeDataByValue If set to true, the in-edges will have their own
 * copy of the edge data. Otherwise, the in-edge edge data will be shared with
 * its corresponding out-edge.
 * @tparam HasNoLockable If set to true, then node accesses will cannot acquire
 * an abstract lock. Otherwise, accessing nodes can get a lock.
 * @tparam UseNumaAlloc If set to true, allocate data in a possibly more NUMA
 * friendly way.
 * @tparam HasOutOfLineLockable
 * @tparam FileEdgeTy
 */
template <typename NodeTy, typename EdgeTy, bool EdgeDataByValue = false,
          bool HasNoLockable = false, bool UseNumaAlloc = false,
          bool HasOutOfLineLockable = false, typename FileEdgeTy = EdgeTy,
          typename NodeIndexTy = uint32_t, typename EdgeIndexTy = uint64_t>
class LC_CSR_CSC_Graph
    : public LC_CSR_Graph<NodeTy, EdgeTy, HasNoLockable, UseNumaAlloc,
                          HasOutOfLineLockable, FileEdgeTy, NodeIndexTy,
                          EdgeIndexTy> {
  // typedef to make it easier to read
  //! Typedef referring to base LC_CSR_Graph
  using BaseGraph =
      LC_CSR_Graph<NodeTy, EdgeTy, HasNoLockable, UseNumaAlloc,
                   HasOutOfLineLockable, FileEdgeTy, NodeIndexTy, EdgeIndexTy>;
  //! Typedef referring to this class itself
  using ThisGraph =
      LC_CSR_CSC_Graph<NodeTy, EdgeTy, EdgeDataByValue, HasNoLockable,
                       UseNumaAlloc, HasOutOfLineLockable, FileEdgeTy,
                       NodeIndexTy, EdgeIndexTy>;

public:
  //! Graph node typedef
  using GraphNode = NodeIndexTy;

protected:
  // redefinitions of base class typedefs
  //! large array for edge data
  using EdgeData = LargeArray<EdgeTy>;
  //! large array for edge destinations
  using EdgeDst = LargeArray<NodeIndexTy>;
  //! large array for edge index data
  using EdgeIndData = LargeArray<EdgeIndexTy>;

public:
  //! iterator for edges
  using edge_iterator =
      boost::counting_iterator<typename EdgeIndData::value_type>;
  //! reference to edge data
  using edge_data_reference = typename EdgeData::reference;

protected:
  //! edge index data for the reverse edges
  EdgeIndData in_edge_ind_data_;
  //! edge destination data for the reverse edges
  EdgeDst in_edge_dst_;
  //! Edge data of inedges can be a value copy of the outedges (i.e. in and
  //! out edges have separate edge values) or inedges can refer to the same
  //! data as its corresponding outedge; this is what this typedef is for
  using EdgeDataRep =
      typename std::conditional<EdgeDataByValue, EdgeData, EdgeIndData>::type;
  //! The data for the reverse edges
  EdgeDataRep inEdgeData;
  EdgeIndData in_edge_to_out_edge_;

  //! redefinition of the edge sort iterator in LC_CSR_Graph
  using edge_sort_iterator =
      internal::EdgeSortIterator<GraphNode, typename EdgeIndData::value_type,
                                 EdgeDst, EdgeDataRep>;

  //! beginning iterator to an edge sorter for in-edges
  edge_sort_iterator in_edge_sort_begin(GraphNode N) {
    return edge_sort_iterator(*in_raw_begin(N), &in_edge_dst_, &inEdgeData);
  }

  //! ending iterator to an edge sorter for in-edges
  edge_sort_iterator in_edge_sort_end(GraphNode N) {
    return edge_sort_iterator(*in_raw_end(N), &in_edge_dst_, &inEdgeData);
  }

  /**
   * Copy the data of outedge by value to inedge.
   *
   * @param e_new position of out-edge to copy as an in-edge
   * @param e position of in-edge
   */
  template <bool A                            = EdgeDataByValue,
            typename std::enable_if<A>::type* = nullptr>
  void createEdgeData(const uint64_t e_new, const uint64_t e) {
    BaseGraph::edgeDataCopy(inEdgeData, BaseGraph::edgeData, e_new, e);
  }

  //! Do nothing; getting edge data will be done via pointer
  template <bool A                             = EdgeDataByValue,
            typename std::enable_if<!A>::type* = nullptr>
  void createEdgeData(const uint64_t, const uint64_t) {
    // do nothing
  }

  /**
   * Determine the in-edge indices for every node by accumulating how many
   * in-edges each node has, getting a prefix sum, and saving it to the
   * in edge index data array.
   *
   * @param dataBuffer temporary buffer that is used to accumulate in-edge
   * counts; at the end of this function, it will contain a prefix sum of
   * in-edges
   */
  void determineInEdgeIndices(EdgeIndData& dataBuffer) {
    // counting outgoing edges in the tranpose graph by
    // counting incoming edges in the original graph
    galois::do_all(galois::iterate(UINT64_C(0), BaseGraph::numEdges),
                   [&](uint64_t e) {
                     auto dst = BaseGraph::edgeDst[e];
                     __sync_add_and_fetch(&(dataBuffer[dst]), 1);
                   });

    // prefix sum calculation of the edge index array
    for (uint32_t n = 1; n < BaseGraph::numNodes; ++n) {
      dataBuffer[n] += dataBuffer[n - 1];
    }

    // copy over the new tranposed edge index data
    in_edge_ind_data_.allocateInterleaved(BaseGraph::numNodes);
    galois::do_all(galois::iterate(UINT64_C(0), BaseGraph::numNodes),
                   [&](uint64_t n) { in_edge_ind_data_[n] = dataBuffer[n]; });
  }

  /**
   * Determine the destination of each in-edge and copy the data associated
   * with an edge (or point to it).
   *
   * @param dataBuffer A prefix sum of in-edges
   */
  void determineInEdgeDestAndData(EdgeIndData& dataBuffer) {
    // after this block dataBuffer[i] will now hold number of edges that all
    // nodes before the ith node have; used to determine where to start
    // saving an edge for a node
    if (BaseGraph::numNodes >= 1) {
      dataBuffer[0] = 0;
      galois::do_all(
          galois::iterate(UINT64_C(1), BaseGraph::numNodes),
          [&](uint64_t n) { dataBuffer[n] = in_edge_ind_data_[n - 1]; });
    }

    // allocate edge dests and data
    in_edge_dst_.allocateInterleaved(BaseGraph::numEdges);

    if (!std::is_void<EdgeTy>::value) {
      inEdgeData.allocateInterleaved(BaseGraph::numEdges);
    }
    in_edge_to_out_edge_.allocateInterleaved(BaseGraph::numEdges);

    galois::do_all(
        galois::iterate(UINT64_C(0), BaseGraph::numNodes), [&](uint64_t src) {
          // e = start index into edge array for a particular node
          uint64_t e = (src == 0) ? 0 : BaseGraph::edgeIndData[src - 1];

          // get all outgoing edges of a particular node in the non-transpose
          // and convert to incoming
          while (e < BaseGraph::edgeIndData[src]) {
            // destination nodde
            auto dst = BaseGraph::edgeDst[e];
            // location to save edge
            auto e_new = __sync_fetch_and_add(&(dataBuffer[dst]), 1);
            // save src as destination
            in_edge_dst_[e_new] = src;
            // edge data to "new" array
            createEdgeData(e_new, e);
            in_edge_to_out_edge_[e_new] = e;
            e++;
          }
        });
  }

public:
  /////////////////////////////////////////////////////////////////////////////
  // Manual construction functions
  /////////////////////////////////////////////////////////////////////////////

  // no edge data support at the moment for these functions because not required
  // for the current use case

  //! Reallocate memory for the CSC part of the graph
  void CSCAllocate() {
    // assumes nodes and edges set from CSR version of this call
    in_edge_dst_.deallocate();
    in_edge_ind_data_.deallocate();

    if (UseNumaAlloc) {
      in_edge_ind_data_.allocateBlocked(BaseGraph::numNodes);
      in_edge_dst_.allocateBlocked(BaseGraph::numEdges);
    } else {
      in_edge_ind_data_.allocateInterleaved(BaseGraph::numNodes);
      in_edge_dst_.allocateInterleaved(BaseGraph::numEdges);
    }
  }
  //! Construct the in edge for some edge index by setting the destination
  void ConstructInEdge(EdgeIndexTy e, NodeIndexTy dst) {
    in_edge_dst_[e] = dst;
  }
  //! In-edge index setting
  void FixEndInEdge(NodeIndexTy n, EdgeIndexTy e) { in_edge_ind_data_[n] = e; }

  /////////////////////////////////////////////////////////////////////////////
  // Construction functions
  /////////////////////////////////////////////////////////////////////////////

  /**
   * Call only after the LC_CSR_Graph part of this class is fully constructed.
   * Creates the in edge data by reading from the out edge data.
   */
  void constructIncomingEdges() {
    galois::StatTimer incomingEdgeConstructTimer("IncomingEdgeConstruct");
    incomingEdgeConstructTimer.start();

    // initialize the temp array
    EdgeIndData dataBuffer;
    dataBuffer.allocateInterleaved(BaseGraph::numNodes);
    galois::do_all(galois::iterate(UINT64_C(0), BaseGraph::numNodes),
                   [&](uint64_t n) { dataBuffer[n] = 0; });

    determineInEdgeIndices(dataBuffer);
    determineInEdgeDestAndData(dataBuffer);

    incomingEdgeConstructTimer.stop();
  }

  /////////////////////////////////////////////////////////////////////////////
  // Access functions
  /////////////////////////////////////////////////////////////////////////////

  /**
   * Grabs in edge beginning without lock/safety.
   *
   * @param N node to get edge beginning of
   * @returns Iterator to first in edge of node N
   */
  edge_iterator in_raw_begin(GraphNode N) const {
    return edge_iterator((N == 0) ? 0 : in_edge_ind_data_[N - 1]);
  }

  /**
   * Grabs in edge end without lock/safety.
   *
   * @param N node to get edge end of
   * @returns Iterator to end of in edges of node N (i.e. first edge of
   * node N+1)
   */
  edge_iterator in_raw_end(GraphNode N) const {
    return edge_iterator(in_edge_ind_data_[N]);
  }

  /**
   * Wrapper to get the in edge end of a node; lock if necessary.
   *
   * @param N node to get edge beginning of
   * @param mflag how safe the acquire should be
   * @returns Iterator to first in edge of node N
   */
  edge_iterator in_edge_begin(GraphNode N,
                              MethodFlag mflag = MethodFlag::WRITE) {
    BaseGraph::acquireNode(N, mflag);
    if (!HasNoLockable && galois::runtime::shouldLock(mflag)) {
      for (edge_iterator ii = in_raw_begin(N), ee = in_raw_end(N); ii != ee;
           ++ii) {
        BaseGraph::acquireNode(in_edge_dst_[*ii], mflag);
      }
    }
    return in_raw_begin(N);
  }

  /**
   * Wrapper to get the in edge end of a node; lock if necessary.
   *
   * @param N node to get in edge end of
   * @param mflag how safe the acquire should be
   * @returns Iterator to end of in edges of node N (i.e. first in edge of N+1)
   */
  edge_iterator in_edge_end(GraphNode N, MethodFlag mflag = MethodFlag::WRITE) {
    BaseGraph::acquireNode(N, mflag);
    return in_raw_end(N);
  }

  uint64_t getInDegree(GraphNode N) const {
    return (in_raw_end(N) - in_raw_begin(N));
  }

  /**
   * @param N node to get in edges for
   * @param mflag how safe the acquire should be
   * @returns Range to in edges of node N
   */
  runtime::iterable<NoDerefIterator<edge_iterator>>
  in_edges(GraphNode N, MethodFlag mflag = MethodFlag::WRITE) {
    return internal::make_no_deref_range(in_edge_begin(N, mflag),
                                         in_edge_end(N, mflag));
  }

  /**
   * Given an edge id for in edges, get the destination of the edge
   *
   * @param ni edge id
   * @returns destination for that in edge
   */
  GraphNode getInEdgeDst(edge_iterator ni) const { return in_edge_dst_[*ni]; }

  /**
   * Given an edge id for in edge, get the data associated with that edge.
   * Returns a constant reference.
   *
   * In-edge has own copy of edge-data version.
   *
   * @param ni in-edge id
   * @returns data of the edge
   */
  template <bool A                            = EdgeDataByValue,
            typename std::enable_if<A>::type* = nullptr>
  edge_data_reference
  getInEdgeData(edge_iterator ni, MethodFlag = MethodFlag::UNPROTECTED) const {
    return inEdgeData[*ni];
  }

  /**
   * Given an edge id for in edge, get the data associated with that edge.
   * Returns a non-constant reference.
   *
   * In-edge has own copy of edge-data version.
   *
   * @param ni in-edge id
   * @returns data of the edge
   */
  template <bool A                            = EdgeDataByValue,
            typename std::enable_if<A>::type* = nullptr>
  edge_data_reference getInEdgeData(edge_iterator ni,
                                    MethodFlag = MethodFlag::UNPROTECTED) {
    return inEdgeData[*ni];
  }

  /**
   * Given an edge id for in edge, get the data associated with that edge.
   * Returns a constant reference.
   *
   * In-edge and out-edge share edge data version.
   *
   * @param ni in-edge id
   * @returns data of the edge
   */
  template <bool A                             = EdgeDataByValue,
            typename std::enable_if<!A>::type* = nullptr>
  edge_data_reference
  getInEdgeData(edge_iterator ni, MethodFlag = MethodFlag::UNPROTECTED) const {
    return BaseGraph::edgeData[in_edge_to_out_edge_[*ni]];
  }

  /**
   * Given an edge id for in edge, get the data associated with that edge.
   * Returns a non-constant reference.
   *
   * In-edge and out-edge share edge data version.
   *
   * @param ni in-edge id
   * @returns data of the edge
   */
  template <bool A                             = EdgeDataByValue,
            typename std::enable_if<!A>::type* = nullptr>
  edge_data_reference getInEdgeData(edge_iterator ni,
                                    MethodFlag = MethodFlag::UNPROTECTED) {
    return BaseGraph::edgeData[in_edge_to_out_edge_[*ni]];
  }

  //! Returns corresponding index for the out-edge corresponding to
  //! an in-edge.
  size_t InEdgeToOutEdge(edge_iterator ni) const {
    return in_edge_to_out_edge_[*ni];
  }

  /**
   * @returns the prefix sum of in-edges
   */
  const EdgeIndData& getInEdgePrefixSum() const { return in_edge_ind_data_; }

  /////////////////////////////////////////////////////////////////////////////
  // Utility
  /////////////////////////////////////////////////////////////////////////////

  /**
   * Sorts outgoing edges of a node. Comparison is over getEdgeDst(e).
   */
  void sortInEdgesByDst(GraphNode N, MethodFlag mflag = MethodFlag::WRITE) {
    BaseGraph::acquireNode(N, mflag);
    // depending on value/ref the type of EdgeSortValue changes
    using EdgeSortVal = EdgeSortValue<
        GraphNode,
        typename std::conditional<EdgeDataByValue, EdgeTy, uint64_t>::type>;

    std::sort(in_edge_sort_begin(N), in_edge_sort_end(N),
              [=](const EdgeSortVal& e1, const EdgeSortVal& e2) {
                return e1.dst < e2.dst;
              });
  }

  /**
   * Sorts all incoming edges of all nodes in parallel. Comparison is over
   * getEdgeDst(e).
   */
  void sortAllInEdgesByDst(MethodFlag mflag = MethodFlag::WRITE) {
    galois::do_all(
        galois::iterate((size_t)0, this->size()),
        [=](GraphNode N) { this->sortInEdgesByDst(N, mflag); },
        galois::no_stats(), galois::steal());
  }

  /**
   * Directly reads the GR file to construct CSR graph
   * and then constructs reverse edges based on that.
   */
  void readAndConstructBiGraphFromGRFile(const std::string& filename) {
    this->readGraphFromGRFile(filename);
    constructIncomingEdges();
  }
};

} // namespace graphs
} // namespace galois
#endif
