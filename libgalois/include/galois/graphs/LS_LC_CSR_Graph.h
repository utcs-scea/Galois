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
#include <atomic>
#include <new>
#include <type_traits>

#include <boost/range/iterator_range_core.hpp>
#include <boost/range/counting_range.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/functional/hash.hpp>

#include <parallel_hashmap/phmap.h>

#include "galois/config.h"
#include "galois/LargeVector.h"
#include "galois/LargeArray.h"

#ifdef __cpp_lib_hardware_interference_size
using std::hardware_destructive_interference_size;
#else
constexpr std::size_t hardware_destructive_interference_size = 64;
#endif

namespace galois::graphs {

/**
 * Local computation graph.
 */
template <typename VertexData, typename EdgeData, bool concurrent = true>
class LS_LC_CSR_Graph : private boost::noncopyable {
public:
  using VertexTopologyID = uint64_t;
  using VertexRange =
      boost::iterator_range<boost::counting_iterator<VertexTopologyID>>;

  using EdgeHandle = std::pair<VertexTopologyID, VertexTopologyID>;

private:
  using SpinLock = galois::substrate::PaddedLock<concurrent>;
  static constexpr bool HasVertexData = !std::is_same_v<VertexData, void>;
  static constexpr bool HasEdgeData   = !std::is_same_v<EdgeData, void>;

  using VertexDataStore =
      std::conditional_t<HasVertexData, typename std::vector<VertexData>,
                         typename std::tuple<>>;
  using EdgeDataStore = std::conditional_t<
      HasEdgeData,
      phmap::flat_hash_map<
          std::pair<VertexTopologyID, VertexTopologyID>, EdgeData,
          boost::hash<std::pair<VertexTopologyID, VertexTopologyID>>>,
      std::tuple<>>;

  // forward-declarations of internal structs
  struct VertexMetadata;
  struct EdgeMetadata;

  class EdgeIterator;
  using EdgeRange = boost::iterator_range<EdgeIterator>;

  VertexDataStore m_vertex_data;
  std::vector<VertexMetadata> m_vertices;

  // m_edges[0] is the CSR with gaps, m_edges[1] is the update log.
  LargeVector<EdgeMetadata> m_edges[2];
  SpinLock m_edges_lock; // guards resizing of edges vectors
  EdgeDataStore m_edge_data;

  alignas(hardware_destructive_interference_size) std::atomic_uint64_t
      m_edges_tail = ATOMIC_VAR_INIT(0);

  // m_holes is the number of holes in the log (m_edges[1])
  alignas(hardware_destructive_interference_size) std::atomic_uint64_t m_holes =
      ATOMIC_VAR_INIT(0);

  // returns a reference to the metadata for the pointed-to edge
  inline EdgeMetadata& getEdgeMetadata(uint8_t buffer, uint64_t index) const {
    return m_edges[buffer][index];
  }

public:
  LS_LC_CSR_Graph(uint64_t num_vertices)
      : m_vertices(num_vertices, VertexMetadata()) {
    if constexpr (HasVertexData) {
      m_vertex_data.resize(num_vertices);
    }
  }

  inline uint64_t size() const noexcept { return m_vertices.size(); }

  /** Data Manipulations **/

  template <typename V = VertexData, typename = std::enable_if<HasVertexData>>
  inline void setData(VertexTopologyID vertex, V data) {
    m_vertex_data[vertex] = data;
  }

  // return data associated with a vertex
  template <typename V = VertexData, typename = std::enable_if<HasVertexData>>
  inline V& getData(VertexTopologyID vertex) {
    return m_vertex_data[vertex];
  }

  template <typename E = EdgeData, typename = std::enable_if<HasEdgeData>>
  void setEdgeData(EdgeHandle handle, E data) {
    m_edge_data[handle] = data;
  }

  inline VertexTopologyID begin() const noexcept {
    return static_cast<VertexTopologyID>(0);
  }

  inline VertexTopologyID end() const noexcept {
    return static_cast<VertexTopologyID>(m_vertices.size());
  }

  VertexRange vertices() { return VertexRange(begin(), end()); }

  VertexTopologyID addVertexTopologyOnly() {
    m_vertices.emplace_back();
    if constexpr (HasVertexData) {
      m_vertex_data.resize(m_vertices.size());
    }
    return m_vertices.size() - 1;
  }

  // Adds multiple vertices to the graph. The new vertices will be assigned
  // consecutive topology IDs, and the lowest new ID is returned.
  template <typename V = VertexData, typename = std::enable_if<HasVertexData>>
  VertexTopologyID addVertices(std::vector<V> data) {
    VertexTopologyID const start = m_vertices.size();
    m_vertices.resize(m_vertices.size() + data.size());
    m_vertex_data.resize(m_vertices.size());

    galois::do_all(
        galois::iterate(0ul, data.size()),
        [&](VertexTopologyID const& off) { setData(start + off, data[off]); });
    return start;
  }

  inline size_t getDegree(VertexTopologyID id) { return m_vertices[id].degree; }

  inline VertexTopologyID getEdgeDst(EdgeHandle eh) { return eh.second; }

  template <typename E = EdgeData, typename = std::enable_if<HasEdgeData>>
  inline E& getEdgeData(EdgeHandle handle) {
    return m_edge_data[handle];
  }

  EdgeRange edges(VertexTopologyID node) {
    auto& vertex_meta = m_vertices[node];

    EdgeMetadata const* const start =
        &getEdgeMetadata(vertex_meta.buffer, vertex_meta.begin);

    EdgeMetadata const* const end =
        &getEdgeMetadata(vertex_meta.buffer, vertex_meta.end);

    return EdgeRange(EdgeIterator(node, start, end),
                     EdgeIterator(node, end, end));
  }

  void addEdges(VertexTopologyID src, const std::vector<VertexTopologyID> dsts,
                std::vector<EdgeData> data) {
    GALOIS_ASSERT(data.size() == dsts.size());
    this->addEdgesTopologyOnly(src, dsts);
    for (size_t i = 0; i < dsts.size(); ++i) {
      m_edge_data[std::make_pair(src, dsts[i])] = data[i];
    }
    // todo: save edge data
  }

  void addEdgesTopologyOnly(VertexTopologyID src,
                            const std::vector<VertexTopologyID> dsts) {

    // Copies the edge list to the end of m_edges[1], prepending
    // the new edges.

    auto& vertex_meta = m_vertices[src];

    uint64_t const new_degree = vertex_meta.degree + dsts.size();
    uint64_t const new_begin =
        m_edges_tail.fetch_add(new_degree, std::memory_order_relaxed);
    uint64_t const new_end = new_begin + new_degree;

    if (m_edges[1].size() < new_end) {
      m_edges_lock.lock();
      {
        if (m_edges[1].size() < new_end)
          m_edges[1].resize(std::max(m_edges[1].size() * 2, new_end));
      }
      m_edges_lock.unlock();
    }

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

    m_holes.fetch_add(vertex_meta.degree, std::memory_order_relaxed);
    vertex_meta.degree += dsts.size();
  }

  void deleteEdges(VertexTopologyID src,
                   const std::vector<VertexTopologyID>& edges) {
    std::unordered_set<VertexTopologyID> edges_set(edges.begin(), edges.end());

    auto& vertex_meta    = m_vertices[src];
    uint64_t holes_added = 0;
    for (auto i = vertex_meta.begin; i < vertex_meta.end; ++i) {
      EdgeMetadata& edge_meta = getEdgeMetadata(vertex_meta.buffer, i);
      if (!edge_meta.is_tomb() &&
          edges_set.find(edge_meta.dst) != edges_set.end()) {
        edge_meta.set_tomb();
        --vertex_meta.degree;
        ++holes_added;
        // remove tombstoned edges from the start of the edge list
        if (i == vertex_meta.begin)
          ++vertex_meta.begin;
      }
    }

    // remove tombstoned edges from the end of the edge list
    for (auto i = vertex_meta.end; i > vertex_meta.begin; --i) {
      if (getEdgeMetadata(vertex_meta.buffer, i - 1).is_tomb()) {
        --vertex_meta.end;
        --vertex_meta.degree;
      } else {
        break;
      }
    }

    m_holes.fetch_add(holes_added, std::memory_order_relaxed);
  }

  // Performs the compaction algorithm by copying any vertices left in buffer 0
  // to buffer 1, then swapping the buffers.
  //
  // Not safe to call in parallel with insertions/deletions.
  void compact() {
    using std::swap;

    // move from buffer 0 to buffer 1
    galois::do_all(galois::iterate(vertices().begin(), vertices().end()),
                   [&](VertexTopologyID vertex_id) {
                     VertexMetadata& vertex_meta = m_vertices[vertex_id];

                     if (vertex_meta.buffer == 0) {
                       this->addEdgesTopologyOnly(vertex_id, {});
                     }

                     // we are about to swap the buffers, so all vertices will
                     // be in buffer 0
                     vertex_meta.buffer = 0;
                   });

    // At this point, there are no more live edges in buffer 0.
    m_edges_lock.lock();
    {
      m_edges[0].resize(0);
      swap(m_edges[0], m_edges[1]);
      // relaxed is fine because of locks held:
      m_edges_tail.store(0, std::memory_order_relaxed);
      m_holes.store(0, std::memory_order_relaxed);
    }
    m_edges_lock.unlock();
  }

  /*
    Compaction policy utilities.
  */

  // Returns an estimated memory usage in bytes for the entire data structure.
  inline size_t getMemoryUsageBytes() {
    size_t estimate = m_vertices.size() * sizeof(VertexMetadata);
    if constexpr (HasVertexData) {
      estimate += m_vertices.size() * sizeof(VertexData);
    }
    m_edges_lock.lock();
    {
      estimate +=
          (m_edges[0].size() + m_edges_tail.load(std::memory_order_relaxed)) *
          sizeof(EdgeMetadata);
    }
    m_edges_lock.unlock();
    if constexpr (HasEdgeData) {
      estimate += m_edge_data.size() *
                  (sizeof(EdgeData) +
                   sizeof(std::pair<VertexTopologyID, VertexTopologyID>));
    }
    return estimate;
  }

  // Returns the number of bytes used for holes in the log.
  inline size_t getLogHolesMemoryUsageBytes() {
    return m_holes.load(std::memory_order_relaxed) * sizeof(EdgeMetadata);
  }

private:
  struct VertexMetadata {
    uint8_t buffer : 1;
    uint64_t begin : 48; // inclusive
    uint64_t end : 48;   // exclusive
    uint64_t degree;

    VertexMetadata() : buffer(0), begin(0), end(0), degree(0) {}

    VertexMetadata(VertexMetadata const& other)
        : buffer(other.buffer), begin(other.begin), end(other.end),
          degree(other.degree) {}

    VertexMetadata(VertexMetadata&& other)
        : buffer(std::move(other.buffer)), begin(std::move(other.begin)),
          end(std::move(other.end)), degree(std::move(other.degree)) {}
  };

  struct EdgeMetadata {
    enum Flags : uint16_t { TOMB = 0x1 };

    uint16_t flags : 16;
    VertexTopologyID dst : 48;

    bool is_tomb() const noexcept { return (flags & TOMB) > 0; }
    void set_tomb() { flags |= TOMB; }
  } __attribute__((packed));

  static_assert(sizeof(EdgeMetadata) <= sizeof(uint64_t));

  class EdgeIterator
      : public boost::iterator_facade<EdgeIterator, EdgeHandle,
                                      boost::forward_traversal_tag,
                                      EdgeHandle const> {
  private:
    VertexTopologyID const src;
    EdgeMetadata const* curr;
    EdgeMetadata const* const end;

    EdgeIterator(VertexTopologyID src, EdgeMetadata const* start,
                 EdgeMetadata const* end)
        : src(src), curr(start), end(end) {}

    void increment() {
      while (++curr < end && curr->is_tomb())
        ;
    }

    // updates to the graph will invalidate iterators
    bool equal(EdgeIterator const& other) const { return curr == other.curr; }

    EdgeHandle dereference() const { return EdgeHandle(src, curr->dst); }

    friend class LS_LC_CSR_Graph;
    friend class boost::iterator_core_access;
  };
};

}; // namespace galois::graphs

#endif
