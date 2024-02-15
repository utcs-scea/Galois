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
#include <iterator>
#include <cstddef>

#include "boost/range/counting_range.hpp"

#include "galois/LargeVector.h"

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

protected:
  // forward-declarations
  struct VertexMetadata;
  struct EdgeMetadata;

  class EdgeIterator;
  class EdgeRange;

  std::vector<VertexMetadata> m_vertices;
  LargeVector<EdgeMetadata> m_edges;

public:
  LS_LC_CSR_Graph(uint64_t num_vertices)
      : m_vertices(num_vertices, VertexMetadata{0, 0}), m_edges() {}

  VertexRange vertices() {
    return VertexRange(static_cast<VertexTopologyID>(0),
                       static_cast<VertexTopologyID>(m_vertices.size()));
  }

  EdgeRange edges(VertexTopologyID node) {
    const auto& node_meta = m_vertices[node];
    return EdgeRange(this, node_meta.begin, node_meta.end);
  }

  int addEdgesTopologyOnly(VertexTopologyID src,
                           const std::vector<VertexTopologyID> dsts) {
    // todo: synchronization?

    auto& node_meta = m_vertices[src];

    const auto new_begin = static_cast<EdgeHandle>(m_edges.size());
    for (auto ii = dsts.begin(); ii != dsts.end(); ++ii) {
      m_edges.emplace_back(EdgeMetadata{
          .flags = 0,
          .dst   = static_cast<VertexTopologyID>(*ii),
      });
    }
    for (auto i = node_meta.begin; i < node_meta.end; ++i) {
      const EdgeMetadata& edge = m_edges[i];
      if (!edge.is_tomb())
        m_edges.push_back(edge);
    }

    node_meta.begin = new_begin;
    node_meta.end   = static_cast<EdgeHandle>(m_edges.size());

    return 0;
  }

  int deleteEdges(VertexTopologyID src,
                  const std::vector<VertexTopologyID>& edges) {
    std::unordered_set<VertexTopologyID> edges_set(edges.begin(), edges.end());
    auto& vertex_meta = m_vertices[src];
    for (auto i = vertex_meta.begin; i < vertex_meta.end; ++i) {
      EdgeMetadata& edge = m_edges[i];
      if (!edge.is_tomb() && edges_set.find(edge.dst) != edges_set.end()) {
        edge.tomb();

        // remove tombstoned edges from the start of the edge list
        if (i == vertex_meta.begin)
          ++vertex_meta.begin;
      }
    }

    // remove tombstoned edges from the end of the edge list
    for (auto i = vertex_meta.end; i > vertex_meta.begin; --i) {
      if (m_edges[i - 1].is_tomb()) {
        --vertex_meta.end;
      } else {
        break;
      }
    }

    return 0;
  }

  VertexTopologyID getEdgeDst(EdgeHandle edge) { return m_edges[edge].dst; }

protected:
  struct VertexMetadata {
    EdgeHandle begin : 48; // inclusive
    EdgeHandle end : 48;   // exclusive
  };

  struct EdgeMetadata {
    enum Flags : uint16_t { TOMB = 0x1 };

    uint16_t flags : 16;
    VertexTopologyID dst : 48;

    bool is_tomb() const noexcept { return (flags & TOMB) > 0; }
    void tomb() { flags |= TOMB; }
  } __attribute__((packed));

  static_assert(sizeof(EdgeMetadata) <= sizeof(uint64_t));

  class EdgeIterator {
    using iterator_category = std::input_iterator_tag;
    using difference_type   = std::ptrdiff_t;
    using value_type        = EdgeHandle;
    using pointer           = void*;

  private:
    const LS_LC_CSR_Graph* m_graph;
    EdgeHandle m_handle;
    const EdgeHandle m_end;

  protected:
    friend class EdgeRange;

    EdgeIterator(const LS_LC_CSR_Graph* graph, EdgeHandle handle,
                 EdgeHandle end)
        : m_graph(graph), m_handle(handle), m_end(end) {}

  public:
    bool operator==(const EdgeIterator& other) const {
      return m_graph == other.m_graph && m_handle == other.m_handle;
    }

    bool operator!=(const EdgeIterator& other) const {
      return !(*this == other);
    }

    EdgeHandle operator*() const {
      assert(!m_graph->m_edges[m_handle].is_tomb());
      return m_handle;
    }

    EdgeIterator& operator++() {
      ++m_handle;

      while (m_handle < m_end && m_graph->m_edges[m_handle].is_tomb())
        ++m_handle;

      return *this;
    }
  };

  class EdgeRange {
  protected:
    friend class LS_LC_CSR_Graph;

    const LS_LC_CSR_Graph* m_graph;
    const EdgeHandle m_begin;
    const EdgeHandle m_end;

    EdgeRange(LS_LC_CSR_Graph* graph, EdgeHandle begin, EdgeHandle end)
        : m_graph(graph), m_begin(begin), m_end(end) {}

  public:
    EdgeIterator begin() { return EdgeIterator(m_graph, m_begin, m_end); }
    EdgeIterator end() { return EdgeIterator(m_graph, m_end, m_end); }
  };
};
}; // namespace galois::graphs

#endif
