/*
 * This file belongs to the Galois project, a C++ library for exploiting
 * parallelism. The code is being released under the terms of the 3-Clause BSD
 * License (a copy is located in LICENSE.txt at the top-level directory).
 *
 * Copyright (C) 2024, The University of Texas at Austin. All rights reserved.
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

#ifndef GALOIS_GRAPHS_LC_CSR_GRAPH_H
#define GALOIS_GRAPHS_LC_CSR_GRAPH_H

#include <unordered_set>

#include "boost/range/counting_range.hpp"

namespace galois::graphs {

/**
 * Local computation graph.
 */
class LS_LC_CSR_Graph : private boost::noncopyable {

public:
  using VertexTopologyID = uint64_t;
  using EdgeHandle       = uint64_t;
  using VertexRange =
      boost::iterator_range<boost::counting_iterator<VertexTopologyID>>;
  using EdgeRange = boost::iterator_range<boost::counting_iterator<EdgeHandle>>;

protected:
  struct VertexMetadata;
  struct EdgeMetadata;
  std::vector<VertexMetadata> Vertices;
  std::vector<EdgeMetadata> Edges;

public:
  LS_LC_CSR_Graph(uint64_t num_vertices)
      : Vertices(num_vertices, VertexMetadata{0, 0}) {}

  VertexRange vertices() {
    return VertexRange(static_cast<VertexTopologyID>(0),
                       static_cast<VertexTopologyID>(Vertices.size()));
  }

  EdgeRange edges(VertexTopologyID node) {
    const auto& node_meta = Vertices[node];
    return EdgeRange(node_meta.begin, node_meta.end);
  }

  int addEdgesTopologyOnly(VertexTopologyID src,
                           const std::vector<VertexTopologyID> dsts) {
    // todo: synchronization?

    auto& node_meta = Vertices[src];

    const auto new_begin = static_cast<EdgeHandle>(Edges.size());
    for (auto ii = dsts.begin(); ii != dsts.end(); ++ii) {
      Edges.emplace_back(EdgeMetadata{
          .flags = 0,
          .dst   = static_cast<VertexTopologyID>(*ii),
      });
    }
    for (auto i = node_meta.begin; i < node_meta.end; ++i) {
      const EdgeMetadata& edge = Edges[i];
      if (!edge.is_tomb())
        Edges.push_back(edge);
    }

    node_meta.begin = new_begin;
    node_meta.end   = static_cast<EdgeHandle>(Edges.size());

    return 0;
  }

  int deleteEdges(VertexTopologyID src,
                  const std::vector<VertexTopologyID> edges) {
    std::unordered_set<VertexTopologyID> edges_set(edges.begin(), edges.end());

    auto& vertex_meta = Vertices[src];
    for (auto i = vertex_meta.begin; i < vertex_meta.end; ++i) {
      EdgeMetadata& edge = Edges[i];
      if (!edge.is_tomb() && edges_set.find(edge.dst) != edges_set.end()) {
        edge.tomb();
        // todo(meyer): optimization of changing vertex begin/end if tomb at
        // front or back
      }
    }

    return 0;
  }

  VertexTopologyID getEdgeDst(EdgeHandle edge) { return Edges[edge].dst; }

protected:
  struct VertexMetadata {
    EdgeHandle begin : 48; // inclusive
    EdgeHandle end : 48;   // exclusive
  };

  struct EdgeMetadata {
    enum Flags : uint16_t { TOMB = 0x0 };

    uint16_t flags : 16;
    VertexTopologyID dst : 48;

    bool is_tomb() const noexcept { return flags & TOMB; }
    void tomb() { flags |= TOMB; }
  } __attribute__((packed));

  static_assert(sizeof(EdgeMetadata) <= sizeof(uint64_t));
};
}; // namespace galois::graphs

#endif
