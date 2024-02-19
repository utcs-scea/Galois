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

#include <boost/range/iterator_range_core.hpp>
#include <boost/range/counting_range.hpp>
#include <boost/iterator/iterator_facade.hpp>

#include "galois/LargeVector.h"

namespace galois::graphs {

/**
 * Local computation graph.
 */
class LS_LC_CSR_Graph : private boost::noncopyable {
public:
  using VertexTopologyID = uint64_t;
  using VertexRange =
      boost::iterator_range<boost::counting_iterator<VertexTopologyID>>;

  struct EdgeHandle {
  private:
    uint8_t buffer : 1;
    uint64_t index : 48;

    EdgeHandle(uint8_t buffer, uint64_t index) : buffer(buffer), index(index) {}

    EdgeHandle(uint64_t const& v) : buffer(v >> 63), index(v) {}

    friend class LS_LC_CSR_Graph;

  public:
    EdgeHandle(EdgeHandle const&) = default;
    EdgeHandle(EdgeHandle&&)      = default;

  } __attribute__((packed));

private:
  // forward-declarations
  struct VertexMetadata;
  struct EdgeMetadata;

  class EdgeIterator;
  using EdgeRange = boost::iterator_range<EdgeIterator>;

  std::vector<VertexMetadata> m_vertices;
  LargeVector<EdgeMetadata> m_edges[2];

  // returns a reference to the metadata for the pointed-to edge
  inline EdgeMetadata& getEdgeMetadata(EdgeHandle const& handle) {
    return getEdgeMetadata(handle.buffer, handle.index);
  }

  inline EdgeMetadata& getEdgeMetadata(uint8_t buffer, uint64_t index) const {
    return m_edges[buffer][index];
  }

public:
  LS_LC_CSR_Graph(uint64_t num_vertices)
      : m_vertices(num_vertices, VertexMetadata()) {}

  VertexRange vertices() {
    return VertexRange(static_cast<VertexTopologyID>(0),
                       static_cast<VertexTopologyID>(m_vertices.size()));
  }

  EdgeRange edges(VertexTopologyID node) {
    auto& vertex_meta = m_vertices[node];
    auto const* ii    = &getEdgeMetadata(vertex_meta.buffer, vertex_meta.begin);
    auto const* ee    = &getEdgeMetadata(vertex_meta.buffer, vertex_meta.end);

    return EdgeRange(EdgeIterator(ii, ee), EdgeIterator(ee, ee));
  }

  int addEdgesTopologyOnly(VertexTopologyID src,
                           const std::vector<VertexTopologyID> dsts) {
    auto& vertex_meta = m_vertices[src];

    // Copies the edge list to the end of m_edges[1], prepending
    // the new edges.

    // todo: acquire lock on vertex

    // todo: acquire lock on m_edges[1]

    const uint64_t new_degree = vertex_meta.degree + dsts.size();
    const uint64_t new_begin  = m_edges[1].size();
    m_edges[1].resize(new_begin + new_degree);
    const uint64_t new_end = m_edges[1].size();

    // todo: release lock on m_edges[1]

    // insert new edges
    std::transform(dsts.begin(), dsts.end(), &getEdgeMetadata(1, new_begin),
                   [](VertexTopologyID dst) {
                     return EdgeMetadata{.flags = 0, .dst = dst};
                   });

    // copy old, non-tombstoned edges
    std::copy_if(&getEdgeMetadata(vertex_meta.buffer, vertex_meta.begin),
                 &getEdgeMetadata(vertex_meta.buffer, vertex_meta.end),
                 &getEdgeMetadata(1, new_begin + dsts.size()),
                 [](EdgeMetadata& edge) { return !edge.is_tomb(); });

    // update vertex metadata
    vertex_meta.buffer = 1;
    vertex_meta.begin  = new_begin;
    vertex_meta.end    = new_end;
    vertex_meta.degree += dsts.size();

    // todo: release vertex lock

    return 0;
  }

  int deleteEdges(VertexTopologyID src,
                  const std::vector<VertexTopologyID>& edges) {
    std::unordered_set<VertexTopologyID> edges_set(edges.begin(), edges.end());

    auto& vertex_meta = m_vertices[src];
    // todo: acquire vertex lock

    for (auto i = vertex_meta.begin; i < vertex_meta.end; ++i) {
      EdgeMetadata& edge_meta =
          getEdgeMetadata(EdgeHandle(vertex_meta.buffer, i));
      if (!edge_meta.is_tomb() &&
          edges_set.find(edge_meta.dst) != edges_set.end()) {
        edge_meta.tomb();
        --vertex_meta.degree;
        // remove tombstoned edges from the start of the edge list
        if (i == vertex_meta.begin)
          ++vertex_meta.begin;
      }
    }

    // remove tombstoned edges from the end of the edge list
    for (auto i = vertex_meta.end; i > vertex_meta.begin; --i) {
      if (getEdgeMetadata(EdgeHandle(vertex_meta.buffer, i - 1)).is_tomb()) {
        --vertex_meta.end;
        --vertex_meta.degree;
      } else {
        break;
      }
    }

    // todo: release vertex lock
    return 0;
  }

  VertexTopologyID getEdgeDst(EdgeHandle edge) {
    return getEdgeMetadata(edge).dst;
  }

  // Performs the compaction algorithm by copying any vertices left in buffer 0
  // to buffer 1, then swapping the buffers.
  //
  // Should not be called from within a Galois parallel kernel.
  void compact() {
    using std::swap;

    // move from buffer 0 to buffer 1
    galois::do_all(galois::iterate(vertices().begin(), vertices().end()),
                   [&](VertexTopologyID vertex_id) {
                     // todo: acquire vertex lock
                     VertexMetadata& vertex_meta = m_vertices[vertex_id];
                     if (vertex_meta.buffer == 0) {
                       // todo: acquire m_edges[1] lock
                       const uint64_t new_begin = m_edges[1].size();
                       m_edges[1].resize(new_begin + vertex_meta.degree);
                       // todo: release m_edges[1] lock
                       const uint64_t new_end = new_begin + vertex_meta.degree;

                       std::copy_if(&getEdgeMetadata(0, vertex_meta.begin),
                                    &getEdgeMetadata(0, vertex_meta.end),
                                    &getEdgeMetadata(1, new_begin),
                                    [](EdgeMetadata& edge_meta) {
                                      return !edge_meta.is_tomb();
                                    });

                       vertex_meta.begin = new_begin;
                       vertex_meta.end   = new_end;
                     }

                     // we are about to swap the buffers
                     vertex_meta.buffer = 0;

                     // don't release the vertex lock until after the edge
                     // arrays are swapped
                   });

    // at this point, we are guaranteed no live edges in m_edges[0], and we hold
    // the lock on all vertices, so we can swap the buffers safely.

    // no need for a lock on m_edges[0] since this is the only place it is
    // updated

    m_edges[0].resize(0);
    swap(m_edges[0], m_edges[1]);

    galois::do_all(galois::iterate(vertices().begin(), vertices().end()),
                   [&](VertexTopologyID vertex_id) {
                     // todo: release vertex lock
                   });
  }

private:
  struct VertexMetadata {
    uint8_t buffer : 1;
    uint64_t begin : 48; // inclusive
    uint64_t end : 48;   // exclusive
    uint64_t degree;

    VertexMetadata() : buffer(0), begin(0), end(0), degree(0) {}
  };

  struct EdgeMetadata {
    enum Flags : uint16_t { TOMB = 0x1 };

    uint16_t flags : 16;
    VertexTopologyID dst : 48;

    bool is_tomb() const noexcept { return (flags & TOMB) > 0; }
    void tomb() { flags |= TOMB; }
  } __attribute__((packed));

  static_assert(sizeof(EdgeMetadata) <= sizeof(uint64_t));

  class EdgeIterator
      : public boost::iterator_facade<EdgeIterator, EdgeMetadata const,
                                      boost::forward_traversal_tag,
                                      VertexTopologyID> {
  private:
    EdgeMetadata const* m_ptr;
    EdgeMetadata const* const m_end;

    explicit EdgeIterator(EdgeMetadata const* ptr, EdgeMetadata const* end)
        : m_ptr(ptr), m_end(end) {}

    void increment() {
      while (++m_ptr < m_end && m_ptr->is_tomb())
        ;
    };

    // note: equality fails across generations
    bool equal(EdgeIterator const& other) const { return m_ptr == other.m_ptr; }

    VertexTopologyID dereference() const { return m_ptr->dst; }

    friend class LS_LC_CSR_Graph;
    friend class boost::iterator_core_access;
  };
};

}; // namespace galois::graphs

#endif
