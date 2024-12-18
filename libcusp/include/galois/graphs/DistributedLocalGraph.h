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
 * @file DistributedLocalGraph.h
 *
 * Contains the implementation for DistLocalGraph. Command line argument
 * definitions are found in DistributedGraph.cpp.
 */

#ifndef _GALOIS_DISTRIBUTED_LOCAL_GRAPH_H
#define _GALOIS_DISTRIBUTED_LOCAL_GRAPH_H

#include <unordered_map>
#include <fstream>

#include "galois/graphs/DistributedGraph.h"
#include "galois/graphs/LS_LC_CSR_Graph.h"
#include "galois/graphs/BufferedGraph.h"
#include "galois/runtime/DistStats.h"
#include "galois/graphs/OfflineGraph.h"
#include "galois/DynamicBitset.h"

/*
 * Headers for boost serialization
 */

namespace galois {
namespace graphs {

/**
 * Base DistLocalGraph class that all distributed graphs extend from.
 *
 * @tparam NodeTy type of node data for the graph
 * @tparam EdgeTy type of edge data for the graph
 */
template <typename NodeTy, typename EdgeTy>
class DistLocalGraph {
private:
  //! Graph name used for printing things
  constexpr static const char* const GRNAME = "dGraph";

  using GraphTy = galois::graphs::LS_LC_CSR_Graph<NodeTy, EdgeTy>;

  // vector for determining range objects for master nodes + nodes
  // with edges (which includes masters)
  //! represents split of all nodes among threads to balance edges
  std::vector<uint32_t> allNodesRanges;
  //! represents split of master nodes among threads to balance edges
  std::vector<uint32_t> masterRanges;
  //! represents split of nodes with edges (includes masters) among threads to
  //! balance edges
  std::vector<uint32_t> withEdgeRanges;
  //! represents split of all nodes among threads to balance in-edges
  std::vector<uint32_t> allNodesRangesIn;
  //! represents split of master nodes among threads to balance in-edges
  std::vector<uint32_t> masterRangesIn;

  using NodeRangeType =
      galois::runtime::SpecificRange<boost::counting_iterator<size_t>>;

  //! Vector of ranges that stores the 3 different range objects that a user is
  //! able to access
  std::vector<NodeRangeType> specificRanges;
  //! Like specificRanges, but for in edges
  std::vector<NodeRangeType> specificRangesIn;

protected:
  //! The internal graph used by DistLocalGraph to represent the graph
  GraphTy* graph;

  //! Marks if the graph is transposed or not.
  bool transposed;

  // global graph variables
  uint64_t numGlobalNodes; //!< Total nodes in the global unpartitioned graph.
  uint64_t numGlobalEdges; //!< Total edges in the global unpartitioned graph.
  uint32_t numNodes;       //!< Num nodes in this graph in total
  uint64_t numEdges;       //!< Num edges in this graph in total

  const unsigned id;       //!< ID of the machine.
  const uint32_t numHosts; //!< Total number of machines

  // local graph
  // size() = Number of nodes created on this host (masters + mirrors)
  uint32_t numOwned;     //!< Number of nodes owned (masters) by this host.
                         //!< size() - numOwned = mirrors on this host
  uint32_t numOwnedInit; //!< Number of nodes owned (masters) by this host that
                         //!< was loaded initially (static graph)
  uint32_t beginMaster;  //!< Local id of the beginning of master nodes.
                         //!< beginMaster + numOwned = local id of the end of
                         //!< master nodes
  uint32_t numNodesWithEdges; //!< Number of nodes (masters + mirrors) that have
                              //!< outgoing edges

  std::vector<uint32_t>
      ownedNodesIndices; //!< Indices of owned nodes that are added dynamically
  //! Information that converts host to range of nodes that host reads
  std::vector<std::pair<uint64_t, uint64_t>> gid2host;
  //! Mirror nodes from different hosts. For reduce
  std::vector<std::vector<size_t>> mirrorNodes;

  //! GID = localToGlobalVector[LID]
  std::vector<uint64_t> localToGlobalVector;
  //! LID = globalToLocalMap[GID]
  std::unordered_map<uint64_t, uint32_t> globalToLocalMap;

  //! Increments evilPhase, a phase counter used by communication.
  void inline increment_evilPhase() {
    ++galois::runtime::evilPhase;
    if (galois::runtime::evilPhase >=
        static_cast<uint32_t>(
            std::numeric_limits<int16_t>::max())) { // limit defined by MPI or
                                                    // LCI
      galois::runtime::evilPhase = 1;
    }
  }

  //! Returns evilPhase + 1, handling loop around as necessary
  unsigned inline evilPhasePlus1() {
    unsigned result = galois::runtime::evilPhase + 1;

    // limit defined by MPI or LCI
    if (result >= uint32_t{std::numeric_limits<int16_t>::max()}) {
      return 1;
    }
    return result;
  }

  //! used to sort edges in the sort edges function
  template <typename GraphNode, typename ET>
  struct IdLess {
    bool
    operator()(const galois::graphs::EdgeSortValue<GraphNode, ET>& e1,
               const galois::graphs::EdgeSortValue<GraphNode, ET>& e2) const {
      return e1.dst < e2.dst;
    }
  };

private:
  /**
   * Given an OfflineGraph, compute the masters for each node by
   * evenly (or unevenly as specified by scale factor)
   * blocking the nodes off to assign to each host. Considers
   * ONLY nodes and not edges.
   *
   * @param g The offline graph which has loaded the graph you want
   * to get the masters for
   * @param scalefactor A vector that specifies if a particular host
   * should have more or less than other hosts
   * @param DecomposeFactor Specifies how decomposed the blocking
   * of nodes should be. For example, a factor of 2 will make 2 blocks
   * out of 1 block had the decompose factor been set to 1.
   */
  void computeMastersBlockedNodes(galois::graphs::OfflineGraph& g,
                                  const std::vector<unsigned>& scalefactor,
                                  unsigned DecomposeFactor = 1) {
    uint64_t numNodes_to_divide = g.size();
    if (scalefactor.empty() || (numHosts * DecomposeFactor == 1)) {
      for (unsigned i = 0; i < numHosts * DecomposeFactor; ++i)
        gid2host.push_back(galois::block_range(uint64_t{0}, numNodes_to_divide,
                                               i, numHosts * DecomposeFactor));
      return;
    }

    // TODO: not compatible with DecomposeFactor.
    assert(scalefactor.size() == numHosts);

    unsigned numBlocks = 0;

    for (unsigned i = 0; i < numHosts; ++i) {
      numBlocks += scalefactor[i];
    }

    std::vector<std::pair<uint64_t, uint64_t>> blocks;
    for (unsigned i = 0; i < numBlocks; ++i) {
      blocks.push_back(
          galois::block_range(uint64_t{0}, numNodes_to_divide, i, numBlocks));
    }

    std::vector<unsigned> prefixSums;
    prefixSums.push_back(0);

    for (unsigned i = 1; i < numHosts; ++i) {
      prefixSums.push_back(prefixSums[i - 1] + scalefactor[i - 1]);
    }

    for (unsigned i = 0; i < numHosts; ++i) {
      unsigned firstBlock = prefixSums[i];
      unsigned lastBlock  = prefixSums[i] + scalefactor[i] - 1;
      gid2host.push_back(
          std::make_pair(blocks[firstBlock].first, blocks[lastBlock].second));
    }
  }

  /**
   * Given an OfflineGraph, compute the masters for each node by
   * evenly (or unevenly as specified by scale factor)
   * blocking the nodes off to assign to each host while taking
   * into consideration the only edges of the node to get
   * even blocks.
   *
   * @param g The offline graph which has loaded the graph you want
   * to get the masters for
   * @param scalefactor A vector that specifies if a particular host
   * should have more or less than other hosts
   * @param DecomposeFactor Specifies how decomposed the blocking
   * of nodes should be. For example, a factor of 2 will make 2 blocks
   * out of 1 block had the decompose factor been set to 1.
   */
  void computeMastersBalancedEdges(galois::graphs::OfflineGraph& g,
                                   const std::vector<unsigned>& scalefactor,
                                   uint32_t edgeWeight,
                                   unsigned DecomposeFactor = 1) {
    if (edgeWeight == 0) {
      edgeWeight = 1;
    }

    auto& net = galois::runtime::getSystemNetworkInterface();

    gid2host.resize(numHosts * DecomposeFactor);
    for (unsigned d = 0; d < DecomposeFactor; ++d) {
      auto r = g.divideByNode(0, edgeWeight, (id + d * numHosts),
                              numHosts * DecomposeFactor, scalefactor);
      gid2host[id + d * numHosts].first  = *(r.first.first);
      gid2host[id + d * numHosts].second = *(r.first.second);
    }

    for (unsigned h = 0; h < numHosts; ++h) {
      if (h == id) {
        continue;
      }
      galois::runtime::SendBuffer b;
      for (unsigned d = 0; d < DecomposeFactor; ++d) {
        galois::runtime::gSerialize(b, gid2host[id + d * numHosts]);
      }
      net.sendTagged(h, galois::runtime::evilPhase, std::move(b));
    }
    net.flush();
    unsigned received = 1;
    while (received < numHosts) {
      decltype(net.recieveTagged(galois::runtime::evilPhase)) p;
      do {
        p = net.recieveTagged(galois::runtime::evilPhase);
      } while (!p);
      assert(p->first != id);
      auto& b = p->second;
      for (unsigned d = 0; d < DecomposeFactor; ++d) {
        galois::runtime::gDeserialize(b, gid2host[p->first + d * numHosts]);
      }
      ++received;
    }
    increment_evilPhase();

#ifndef NDEBUG
    for (unsigned h = 0; h < numHosts; h++) {
      if (h == 0) {
        assert(gid2host[h].first == 0);
      } else if (h == numHosts - 1) {
        assert(gid2host[h].first == gid2host[h - 1].second);
        assert(gid2host[h].second == g.size());
      } else {
        assert(gid2host[h].first == gid2host[h - 1].second);
        assert(gid2host[h].second == gid2host[h + 1].first);
      }
    }
#endif
  }

  /**
   * Given an OfflineGraph, compute the masters for each node by
   * evenly (or unevenly as specified by scale factor)
   * blocking the nodes off to assign to each host while taking
   * into consideration the edges of the node AND the node itself.
   *
   * @param g The offline graph which has loaded the graph you want
   * to get the masters for
   * @param scalefactor A vector that specifies if a particular host
   * should have more or less than other hosts
   * @param DecomposeFactor Specifies how decomposed the blocking
   * of nodes should be. For example, a factor of 2 will make 2 blocks
   * out of 1 block had the decompose factor been set to 1. Ignored
   * in this function currently.
   *
   * @todo make this function work with decompose factor
   */
  void computeMastersBalancedNodesAndEdges(
      galois::graphs::OfflineGraph& g, const std::vector<unsigned>& scalefactor,
      uint32_t nodeWeight, uint32_t edgeWeight, unsigned) {
    if (nodeWeight == 0) {
      nodeWeight = g.sizeEdges() / g.size(); // average degree
    }
    if (edgeWeight == 0) {
      edgeWeight = 1;
    }

    auto& net = galois::runtime::getSystemNetworkInterface();
    gid2host.resize(numHosts);
    auto r = g.divideByNode(nodeWeight, edgeWeight, id, numHosts, scalefactor);
    gid2host[id].first  = *r.first.first;
    gid2host[id].second = *r.first.second;
    for (unsigned h = 0; h < numHosts; ++h) {
      if (h == id)
        continue;
      galois::runtime::SendBuffer b;
      galois::runtime::gSerialize(b, gid2host[id]);
      net.sendTagged(h, galois::runtime::evilPhase, std::move(b));
    }
    net.flush();
    unsigned received = 1;
    while (received < numHosts) {
      decltype(net.recieveTagged(galois::runtime::evilPhase)) p;
      do {
        p = net.recieveTagged(galois::runtime::evilPhase);
      } while (!p);
      assert(p->first != id);
      auto& b = p->second;
      galois::runtime::gDeserialize(b, gid2host[p->first]);
      ++received;
    }
    increment_evilPhase();
  }

protected:
  /**
   * Wrapper call that will call into more specific compute masters
   * functions that compute masters based on nodes, edges, or both.
   *
   * @param masters_distribution method of masters distribution to use
   * @param g The offline graph which has loaded the graph you want
   * to get the masters for
   * @param scalefactor A vector that specifies if a particular host
   * should have more or less than other hosts
   * @param nodeWeight weight to give nodes when computing balance
   * @param edgeWeight weight to give edges when computing balance
   * @param DecomposeFactor Specifies how decomposed the blocking
   * of nodes should be. For example, a factor of 2 will make 2 blocks
   * out of 1 block had the decompose factor been set to 1.
   */
  uint64_t computeMasters(MASTERS_DISTRIBUTION masters_distribution,
                          galois::graphs::OfflineGraph& g,
                          const std::vector<unsigned>& scalefactor,
                          uint32_t nodeWeight = 0, uint32_t edgeWeight = 0,
                          unsigned DecomposeFactor = 1) {
    galois::Timer timer;
    timer.start();
    g.reset_seek_counters();

    uint64_t numNodes_to_divide = g.size();

    // compute masters for all nodes
    switch (masters_distribution) {
    case BALANCED_MASTERS:
      computeMastersBlockedNodes(g, scalefactor, DecomposeFactor);
      break;
    case BALANCED_MASTERS_AND_EDGES:
      computeMastersBalancedNodesAndEdges(g, scalefactor, nodeWeight,
                                          edgeWeight, DecomposeFactor);
      break;
    case BALANCED_EDGES_OF_MASTERS:
    default:
      computeMastersBalancedEdges(g, scalefactor, edgeWeight, DecomposeFactor);
      break;
    }

    timer.stop();

    galois::runtime::reportStatCond_Tmax<MORE_DIST_STATS>(
        GRNAME, "MasterDistTime", timer.get());

    galois::gPrint(
        "[", id, "] Master distribution time : ", timer.get_usec() / 1000000.0f,
        " seconds to read ", g.num_bytes_read(), " bytes in ", g.num_seeks(),
        " seeks (", g.num_bytes_read() / (float)timer.get_usec(), " MBPS)\n");
    return numNodes_to_divide;
  }

  //! reader assignment from a file
  //! corresponds to master assignment if using an edge cut
  void readersFromFile(galois::graphs::OfflineGraph& g, std::string filename) {
    // read file lines
    std::ifstream mappings(filename);
    std::string curLine;

    unsigned timesToRead = id + 1;

    for (unsigned i = 0; i < timesToRead; i++) {
      std::getline(mappings, curLine);
    }

    std::vector<char> modifyLine(curLine.begin(), curLine.end());
    char* tokenizedString = modifyLine.data();
    char* token;
    token = strtok(tokenizedString, " ");

    // loop 6 more times
    for (unsigned i = 0; i < 6; i++) {
      token = strtok(NULL, " ");
    }
    std::string left(token);

    // 3 more times for right
    for (unsigned i = 0; i < 3; i++) {
      token = strtok(NULL, " ");
    }
    std::string right(token);

    gid2host.resize(numHosts);
    gid2host[id].first  = std::stoul(left);
    gid2host[id].second = std::stoul(right) + 1;
    galois::gPrint("[", id, "] Left: ", gid2host[id].first,
                   ", Right: ", gid2host[id].second, "\n");

    /////////////////////////
    // send/recv from other hosts
    /////////////////////////
    auto& net = galois::runtime::getSystemNetworkInterface();

    for (unsigned h = 0; h < numHosts; ++h) {
      if (h == id)
        continue;
      galois::runtime::SendBuffer b;
      galois::runtime::gSerialize(b, gid2host[id]);
      net.sendTagged(h, galois::runtime::evilPhase, std::move(b));
    }
    net.flush();
    unsigned received = 1;
    while (received < numHosts) {
      decltype(net.recieveTagged(galois::runtime::evilPhase)) p;
      do {
        p = net.recieveTagged(galois::runtime::evilPhase);
      } while (!p);
      assert(p->first != id);
      auto& b = p->second;
      galois::runtime::gDeserialize(b, gid2host[p->first]);
      ++received;
    }
    increment_evilPhase();

    // sanity checking assignment
    for (unsigned h = 0; h < numHosts; h++) {
      if (h == 0) {
        GALOIS_ASSERT(gid2host[h].first == 0);
      } else if (h == numHosts - 1) {
        GALOIS_ASSERT(gid2host[h].first == gid2host[h - 1].second,
                      gid2host[h].first, " ", gid2host[h - 1].second);
        GALOIS_ASSERT(gid2host[h].second == g.size(), gid2host[h].second, " ",
                      g.size());
      } else {
        GALOIS_ASSERT(gid2host[h].first == gid2host[h - 1].second,
                      gid2host[h].first, " ", gid2host[h - 1].second);
        GALOIS_ASSERT(gid2host[h].second == gid2host[h + 1].first,
                      gid2host[h].second, " ", gid2host[h + 1].first);
      }
    }
  }

  uint32_t G2L(uint64_t gid) const {
    assert(isLocal(gid));
    return globalToLocalMap.at(gid);
  }

  uint64_t L2G(uint32_t lid) const { return localToGlobalVector[lid]; }

public:
  //! Type representing a node in this graph
  using GraphNode = typename GraphTy::VertexTopologyID;
  //! Type representing an edge data in this graph
  using EdgeType = EdgeTy;
  //! iterator type over edges
  using edge_iterator = typename GraphTy::EdgeIterator;

  /**
   * Constructor for DistLocalGraph. Initializes metadata fields.
   *
   * @param host host number that this graph resides on
   * @param numHosts total number of hosts in the currently executing program
   */
  DistLocalGraph(unsigned host, unsigned numHosts)
      : transposed(false), id(host), numHosts(numHosts) {
    mirrorNodes.resize(numHosts);
    numGlobalNodes = 0;
    numGlobalEdges = 0;
  }

  /**
   * Return a vector of pairs denoting mirror node ranges.
   *
   * Assumes all mirror nodes occur after the masters: this invariant should be
   * held by CuSP.
   */
  std::vector<std::pair<uint32_t, uint32_t>> getMirrorRanges() const {
    std::vector<std::pair<uint32_t, uint32_t>> mirrorRangesVector;
    // order of nodes locally is masters, outgoing mirrors, incoming mirrors,
    // so just get from numOwned to end
    if (numOwned != numNodes) {
      assert(numOwned < numNodes);
      mirrorRangesVector.push_back(std::make_pair(numOwned, numNodes));
    }
    return mirrorRangesVector;
  }

  std::vector<std::vector<size_t>>& getMirrorNodes() { return mirrorNodes; }

private:
  virtual unsigned getHostIDImpl(uint64_t) const = 0;
  virtual bool isOwnedImpl(uint64_t) const       = 0;
  virtual bool isLocalImpl(uint64_t) const       = 0;
  virtual bool isVertexCutImpl() const           = 0;
  virtual std::pair<unsigned, unsigned> cartesianGridImpl() const {
    return std::make_pair(0u, 0u);
  }

public:
  virtual ~DistLocalGraph() {}
  void initGraph(uint64_t numNodes) { graph = new GraphTy(numNodes); }
  //! Determines which host has the master for a particular node
  //! @returns Host id of node in question
  inline unsigned getHostID(uint64_t gid) const { return getHostIDImpl(gid); }
  //! Determine if a node has a master on this host.
  //! @returns True if passed in global id has a master on this host
  inline bool isOwned(uint64_t gid) const { return isOwnedImpl(gid); }
  //! Determine if a node has a proxy on this host
  //! @returns True if passed in global id has a proxy on this host
  inline bool isLocal(uint64_t gid) const { return isLocalImpl(gid); }
  /**
   * Returns true if current partition is a vertex cut
   * @returns true if partition being stored in this graph is a vertex cut
   */
  inline bool is_vertex_cut() const { return isVertexCutImpl(); }
  /**
   * Returns Cartesian split (if it exists, else returns pair of 0s
   */
  inline std::pair<unsigned, unsigned> cartesianGrid() const {
    return cartesianGridImpl();
  }

  bool isTransposed() { return transposed; }

  /**
   * Converts a local node id into a global node id
   *
   * @param nodeID local node id
   * @returns global node id corresponding to the local one
   */
  inline uint64_t getGID(const uint32_t nodeID) const { return L2G(nodeID); }

  /**
   * Converts a global node id into a local node id
   *
   * @param nodeID global node id
   * @returns local node id corresponding to the global one
   */
  inline uint32_t getLID(const uint64_t nodeID) const { return G2L(nodeID); }

  /**
   * Get data of a node.
   *
   * @param N node to get the data of
   * @param mflag access flag for node data
   * @returns A node data object
   */
  inline NodeTy& getData(GraphNode N) {
    auto& r = graph->getData(N);
    return r;
  }

  /**
   * Get the edge data for a particular edge in the graph.
   *
   * @param ni edge to get the data of
   * @param mflag access flag for edge data
   * @returns The edge data for the requested edge
   */
  inline EdgeTy& getEdgeData(GraphNode src, edge_iterator ni) {
    GraphNode dst = getEdgeDst(ni);
    auto& r       = graph->getEdgeData(std::make_pair(src, dst));
    return r;
  }

  inline EdgeTy& getEdgeData(edge_iterator ni) {
    auto& r = graph->getEdgeData(*ni);
    return r;
  }

  /**
   * Gets edge destination of edge ni.
   *
   * @param ni edge id to get destination of
   * @returns Local ID of destination of edge ni
   */
  GraphNode getEdgeDst(edge_iterator ni) { return graph->getEdgeDst(*ni); }

  /**
   * Gets the first edge of some node.
   *
   * @param N node to get the edge of
   * @returns iterator to first edge of N
   */
  inline edge_iterator edge_begin(GraphNode N) {
    return graph->edges(N).begin();
  }

  /**
   * Gets the end edge boundary of some node.
   *
   * @param N node to get the edge of
   * @returns iterator to the end of the edges of node N, i.e. the first edge
   * of the next node (or an "end" iterator if there is no next node)
   */
  inline edge_iterator edge_end(GraphNode N) { return graph->edges(N).end(); }

  /**
   * Return the degree of the edge in the local graph
   **/
  inline uint64_t localDegree(GraphNode N) { return graph->getDegree(N); }

  /**
   * Returns an iterable object over the edges of a particular node in the
   * graph.
   *
   * @param N node to get edges iterator over
   */
  inline galois::runtime::iterable<galois::NoDerefIterator<edge_iterator>>
  edges(GraphNode N) {
    return galois::graphs::internal::make_no_deref_range(edge_begin(N),
                                                         edge_end(N));
  }

  /**
   * Gets number of nodes on this (local) graph.
   *
   * @returns number of nodes present in this (local) graph
   */
  inline size_t size() const { return graph->size(); }

  /**
   * Gets number of edges on this (local) graph.
   *
   * @returns number of edges present in this (local) graph
   */
  inline size_t sizeEdges() { return graph->sizeEdges(); }

  /**
   * Gets number of nodes on this (local) graph.
   *
   * @returns number of nodes present in this (local) graph
   */
  inline size_t numMasters() const { return numOwned; }

  /**
   * Gets number of nodes with edges (may include nodes without edges)
   * on this (local) graph.
   *
   * @returns number of nodes with edges (may include nodes without edges
   * as it measures a contiguous range)
   */
  inline size_t getNumNodesWithEdges() const { return numNodesWithEdges; }

  /**
   * Gets number of nodes on the global unpartitioned graph.
   *
   * @returns number of nodes present in the global unpartitioned graph
   */
  inline size_t globalSize() const { return numGlobalNodes; }

  /**
   * Gets number of edges on the global unpartitioned graph.
   *
   * @returns number of edges present in the global unpartitioned graph
   */
  inline size_t globalSizeEdges() const { return numGlobalEdges; }

  /**
   * Returns a range object that encapsulates all nodes of the graph.
   *
   * @returns A range object that contains all the nodes in this graph
   */
  inline const NodeRangeType& allNodesRange() const {
    assert(specificRanges.size() == 3);
    return specificRanges[0];
  }

  size_t getValue(uint32_t lid) { return ownedNodesIndices[lid]; }

  /**
   * Returns a range object that encapsulates only master nodes in this
   * graph.
   *
   * @returns A range object that contains the master nodes in this graph
   */
  inline const NodeRangeType& masterNodesRange() const {
    assert(specificRanges.size() == 3);
    return specificRanges[1];
  }

  /**
   * Returns a range object that encapsulates master nodes and nodes
   * with edges in this graph.
   *
   * @returns A range object that contains the master nodes and the nodes
   * with outgoing edges in this graph
   */
  inline const NodeRangeType& allNodesWithEdgesRange() const {
    assert(specificRanges.size() == 3);
    return specificRanges[2];
  }

  void setNumOwnedInit(uint32_t num) {
    numOwnedInit = num;
    for (uint32_t i = 0; i < num; i++) {
      ownedNodesIndices.push_back(i);
    }
  }

  /**
   * Returns a vector object that contains the global IDs (in order) of
   * the master nodes in this graph.
   *
   * @returns A vector object that contains the global IDs (in order) of
   * the master nodes in this graph
   */
  std::vector<uint64_t> getMasterGlobalIDs() {
    std::vector<uint64_t> IDs;

    IDs.reserve(numMasters());
    for (auto node : masterNodesRange()) {
      IDs.push_back(getGID(node));
    }

    return IDs;
  }

protected:
  /**
   * Uses a pre-computed prefix sum to determine division of nodes among
   * threads.
   *
   * The call uses binary search to determine the ranges.
   */
  inline void determineThreadRanges() {
    allNodesRanges = galois::graphs::determineUnitRangesFromPrefixSum(
        galois::runtime::activeThreads, graph->getEdgePrefixSum());
  }

  /**
   * Determines the thread ranges for master nodes only and saves them to
   * the object.
   *
   * Only call after graph is constructed + only call once
   */
  inline void determineThreadRangesMaster() {
    // make sure this hasn't been called before
    if (masterRanges.size() != 0)
      masterRanges.clear();

    // first check if we even need to do any work; if already calculated,
    // use already calculated vector
    if (beginMaster == 0 && (beginMaster + numOwned) == size()) {
      masterRanges = allNodesRanges;
    } else if (beginMaster == 0 &&
               (beginMaster + numOwned) == numNodesWithEdges &&
               withEdgeRanges.size() != 0) {
      masterRanges = withEdgeRanges;
    } else {
      galois::gDebug("Manually det. master thread ranges");
      masterRanges = galois::graphs::determineUnitRangesFromGraph(
          *graph, galois::runtime::activeThreads, beginMaster,
          beginMaster + numOwned, 0, true);
    }
  }

  /**
   * Determines the thread ranges for nodes with edges only and saves them to
   * the object.
   *
   * Only call after graph is constructed + only call once
   */
  inline void determineThreadRangesWithEdges() {
    // make sure not called before
    if (withEdgeRanges.size() != 0)
      withEdgeRanges.clear();

    // first check if we even need to do any work; if already calculated,
    // use already calculated vector
    if (numNodesWithEdges == size()) {
      withEdgeRanges = allNodesRanges;
    } else if (beginMaster == 0 &&
               (beginMaster + numOwned) == numNodesWithEdges &&
               masterRanges.size() != 0) {
      withEdgeRanges = masterRanges;
    } else {
      galois::gDebug("Manually det. with edges thread ranges");
      withEdgeRanges = galois::graphs::determineUnitRangesFromGraph(
          *graph, galois::runtime::activeThreads, 0, numNodesWithEdges, 0);
    }
  }

  /**
   * Initializes the 3 range objects that a user can access to iterate
   * over the graph in different ways.
   */
  void initializeSpecificRanges() {
    if (specificRanges.size() != 0)
      specificRanges.clear();

    // TODO/FIXME assertion likely not safe if a host gets no nodes
    // make sure the thread ranges have already been calculated
    // for the 3 ranges
    assert(allNodesRanges.size() != 0);
    assert(masterRanges.size() != 0);
    assert(withEdgeRanges.size() != 0);

    // 0 is all nodes
    specificRanges.push_back(galois::runtime::makeSpecificRange(
        boost::counting_iterator<size_t>(0),
        boost::counting_iterator<size_t>(size()), allNodesRanges.data()));

    // 1 is master nodes
    specificRanges.push_back(galois::runtime::makeSpecificRange(
        boost::counting_iterator<size_t>(beginMaster),
        boost::counting_iterator<size_t>(beginMaster + numOwned),
        masterRanges.data()));

    // 2 is with edge nodes
    specificRanges.push_back(galois::runtime::makeSpecificRange(
        boost::counting_iterator<size_t>(0),
        boost::counting_iterator<size_t>(numNodesWithEdges),
        withEdgeRanges.data()));

    assert(specificRanges.size() == 3);
  }

  /**
   * Specific range editor: makes the range for edges equivalent to the range
   * for masters.
   */
  void edgesEqualMasters() { specificRanges[2] = specificRanges[1]; }

  void recalculateG2LMap() {
    for (uint64_t i = 0; i < localToGlobalVector.size(); i++) {
      globalToLocalMap[localToGlobalVector[i]] = i;
    }
  }

public:
  /**
   * Write the local LC_CSR graph to the file on a disk.
   *
   * @todo revive this
   */
  void save_local_graph_to_file(std::string) { GALOIS_DIE("not implemented"); }

  /**
   * Read the local LC_CSR graph from the file on a disk.
   *
   * @todo revive this
   */
  void read_local_graph_from_file(std::string) {
    GALOIS_DIE("not implemented");
  }

  /**
   * Deallocates underlying LC CSR Graph
   */
  void deallocate() {
    galois::gDebug("Deallocating CSR in DistLocalGraph");
    graph->deallocate();
  }

  /**
   * Sort the underlying LC_CSR_Graph by ID (destinations)
   * It sorts edges of the nodes by destination.
   */
  void sortEdgesByDestination() {
    galois::do_all(
        galois::iterate(graph->vertices().begin(), graph->vertices().end()),
        [&](GraphNode n) { graph->sortEdges(n); }, galois::no_stats(),
        galois::loopname("CSREdgeSort"), galois::steal());
  }

  //! Used by substrate to determine if some stats are to be reported
  bool is_a_graph() const { return true; }
  inline NodeTy& getTopologyID(uint64_t nodeID) {
    return graph.getData(getLID(nodeID));
  }

  inline NodeTy& getTopologyIDFromIndex(uint64_t index) {
    return graph.getData(index);
  }

  uint64_t getTokenID(NodeTy& vertex) {
    return getGID(&vertex - &graph.getData(0));
  }

  uint32_t getVertexIndex(NodeTy& vertex) {
    return (&vertex - &graph.getData(0));
  }

  uint64_t getLocalityVertex(NodeTy& vertex) {
    uint64_t gid = getTopologyID(vertex);
    return getHostIDImpl(gid);
  }

  /** Edge Manipulation **/
  edge_iterator mintEdgeHandle(NodeTy& src, std::uint64_t off) {
    return edge_begin(src) + off;
  }

  // template <typename T = NodeTy>
  // typename std::enable_if<!std::is_void<T>::value>::type
  // setData(typename GraphTy::node_data_reference vertex, T data) {
  //   graph.setData(vertex, data);
  // }

  ///** Data Manipulations **/

  // typename GraphTy::node_data_reference
  // getData(typename GraphTy::node_data_reference vertex) {
  //   return graph.getData(getTokenID(vertex));
  // }

  template <typename T = NodeTy>
  typename std::enable_if<!std::is_void<T>::value>::type
  setEdgeData(edge_iterator eh, T data) {
    graph.setEdgeData(eh, data);
  }

  template <typename T = NodeTy>
  typename std::enable_if<!std::is_void<T>::value, EdgeTy&>::type
  getEdgeData(edge_iterator eh) {
    return graph.getEdgeData(eh);
  }

  enum Task {
    ADD_VERTEX,
    ADD_VERTEX_TOPOLOGY_ONLY,
    ADD_EDGES,
    ADD_EDGES_TOPOLOGY_ONLY,
    DELETE_VERTEX,
    DELETE_EDGES
  };

  template <typename... Args>
  void sendModifyRequest(uint32_t host, Args... args) {
    galois::runtime::SendBuffer b;
    galois::runtime::gSerialize(b, args...);
    galois::runtime::getSystemNetworkInterface().sendTagged(
        host, galois::runtime::evilPhase, std::move(b));
  }

  void updateRanges() {
    determineThreadRanges();
    determineThreadRangesMaster();
    determineThreadRangesWithEdges();
    initializeSpecificRanges();
  }

  // Assumptions:
  //  1. A vertex is added before any edges are added to it
  //  2. No support for deleting edges/vertices yet
  //  3. Only works for OEC
  void
  updateVariables(bool isVertex, uint64_t src,
                  std::optional<std::vector<uint64_t>> dsts  = std::nullopt,
                  std::optional<std::vector<NodeTy>> dstData = std::nullopt) {

    if (isVertex) {
      assert(globalToLocalMap.find(src) == globalToLocalMap.end());
      localToGlobalVector.push_back(src);
      globalToLocalMap[src] = localToGlobalVector.size() - 1;
      numNodes++;
      ownedNodesIndices.push_back(numNodes - 1);
      numOwned++;
    } else {
      assert(globalToLocalMap.find(src) != globalToLocalMap.end());
      uint64_t srcLID = globalToLocalMap[src];
      if (edge_begin(srcLID) == edge_end(srcLID)) {
        numNodesWithEdges++;
      }
      uint32_t i = 0;
      for (auto token : dsts.value()) {
        if (globalToLocalMap.find(token) == globalToLocalMap.end()) {
          localToGlobalVector.push_back(token);
          globalToLocalMap[token] = localToGlobalVector.size() - 1;
          numNodes++;
          numNodesWithEdges++;
          std::vector<NodeTy> data;
          data.push_back(dstData.value()[i]);
          graph->addVertices(data);
          mirrorNodes[getHostID(token)].push_back(token);
        }
        i++;
        if (isOwned(token) &&
            (edge_begin(getLID(token)) == edge_end(getLID(token)))) {
          numNodesWithEdges++;
        }
      }
      numEdges += dsts.value().size();
    }
  }

  /** Topology Modifications **/
  void addVertexTopologyOnly(uint32_t token) {
    updateVariables(true, token);
    graph->addVertexTopologyOnly();
  }

  template <typename T>
  void addVertex(T data) {
    uint64_t token = data.id;
    std::vector<T> dataVec;
    dataVec.push_back(data);
    updateVariables(true, token);
    graph->addVertices(dataVec);
  }

  void addEdgesTopologyOnly(uint64_t src, std::vector<uint64_t> dsts) {
    updateVariables(false, src, dsts);
    std::vector<uint64_t> lids;
    for (uint32_t i = 0; i < dsts.size(); i++) {
      lids.push_back(getLID(dsts[i]));
    }
    graph->addEdgesTopologyOnly(getLID(src), lids);
  }

  void addEdges(uint64_t src, std::vector<uint64_t> dsts,
                std::vector<EdgeTy> data, std::vector<NodeTy> dstData) {
    updateVariables(false, src, dsts, dstData);
    std::vector<uint64_t> lids;
    for (uint32_t i = 0; i < dsts.size(); i++) {
      lids.push_back(getLID(dsts[i]));
    }
    graph->addEdges(getLID(src), lids, data);
  }

  void deleteVertex(uint64_t src) {
    // TODO(Divija): Uncomment when we have the graph API
    // graph.deleteVertex(getLID(src));
  }

  void deleteEdges(uint64_t src, std::vector<edge_iterator> edges) {
    // TODO:Remove dst tokens from local map?
    // TODO(Divija): Uncomment when we have the graph API
    // return graph.deleteEdges(getLID(src), edges);
  }
};

template <typename NodeTy, typename EdgeTy>
constexpr const char* const
    galois::graphs::DistLocalGraph<NodeTy, EdgeTy>::GRNAME;
} // end namespace graphs
} // end namespace galois

#endif //_GALOIS_DISTRIBUTED_LOCAL_GRAPH_H
