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
#include "galois/PrefixSum.h"

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

  class EdgeIterator;
  using EdgeRange = boost::iterator_range<EdgeIterator>;

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
  using EdgeMetadata = VertexTopologyID;
  struct VertexMetadata {
    EdgeMetadata* begin = nullptr; // inclusive
    EdgeMetadata* end   = nullptr; // exclusive

    inline uint64_t degree() const { return end - begin; }
  };

  VertexDataStore m_vertex_data;
  std::vector<VertexMetadata> m_vertices;

  // m_edges is the CSR with gaps, m_edit_log is the update log.
  EdgeDataStore m_edge_data;

  SpinLock m_edges_lock; // guards resizing of both edge vectors:
  LargeVector<EdgeMetadata> m_edges;
  LargeVector<EdgeMetadata> m_edit_log;
  alignas(hardware_destructive_interference_size) std::atomic_uint64_t
      m_edges_tail = ATOMIC_VAR_INIT(0);

  // m_holes is the number of holes in the log (m_edit_log)
  alignas(hardware_destructive_interference_size) std::atomic_uint64_t m_holes =
      ATOMIC_VAR_INIT(0);

  /*
   * Prefix Sum utilities
   */
  std::vector<uint64_t> m_pfx_sum_cache;
  static uint64_t transmute(const VertexMetadata& vertex_meta) {
    return vertex_meta.degree();
  }
  static uint64_t scan_op(const VertexMetadata& p, const uint64_t& l) {
    return p.degree() + l;
  }
  static uint64_t combiner(const uint64_t& f, const uint64_t& s) {
    return f + s;
  }
  PrefixSum<VertexMetadata, uint64_t, transmute, scan_op, combiner,
            CacheLinePaddedArr>
      m_pfx{&m_vertices[0], &m_pfx_sum_cache[0]};

  alignas(hardware_destructive_interference_size)
      std::atomic<bool> m_prefix_valid = ATOMIC_VAR_INIT(false);

  void resetPrefixSum() {
    m_pfx_sum_cache.resize(m_vertices.size());
    m_pfx.src = &m_vertices[0];
    m_pfx.dst = &m_pfx_sum_cache[0];
  }

  // Compute the prefix sum using the two level method
  void computePrefixSum() {
    m_pfx.computePrefixSum(m_vertices.size());
    m_prefix_valid.store(true, std::memory_order_release);
  }

public:
  LS_LC_CSR_Graph(uint64_t num_vertices)
      : m_vertices(num_vertices, VertexMetadata()) {
    if constexpr (HasVertexData) {
      m_vertex_data.resize(num_vertices);
    }
    resetPrefixSum();
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
    resetPrefixSum();
    return m_vertices.size() - 1;
  }

  // Adds multiple vertices to the graph. The new vertices will be assigned
  // consecutive topology IDs, and the lowest new ID is returned.
  template <typename V = VertexData, typename = std::enable_if<HasVertexData>>
  VertexTopologyID addVertices(std::vector<V> data) {
    VertexTopologyID const start = m_vertices.size();
    m_vertices.resize(m_vertices.size() + data.size());
    m_vertex_data.resize(m_vertices.size());
    resetPrefixSum();
    galois::do_all(
        galois::iterate(0ul, data.size()),
        [&](VertexTopologyID const& off) { setData(start + off, data[off]); });
    return start;
  }

  inline size_t getDegree(VertexTopologyID id) {
    return m_vertices[id].degree();
  }

  inline VertexTopologyID getEdgeDst(EdgeHandle eh) { return eh.second; }

  template <typename E = EdgeData, typename = std::enable_if<HasEdgeData>>
  inline E& getEdgeData(EdgeHandle handle) {
    return m_edge_data[handle];
  }

  /*
   * Count the total number of edges in parallel.
   */
  uint64_t sizeEdges() {
    galois::GAccumulator<uint64_t> num_edges;
    num_edges.reset();
    galois::do_all(galois::iterate(begin(), end()),
                   [&](VertexTopologyID const& vertex) {
                     num_edges += getDegree(vertex);
                   });
    return num_edges.reduce();
  }

  inline EdgeIterator edge_begin(VertexTopologyID vertex) {
    auto const& vertex_meta = m_vertices[vertex];
    return EdgeIterator(vertex, vertex_meta.begin);
  }

  inline EdgeIterator edge_end(VertexTopologyID vertex) {
    VertexMetadata const& vertex_meta = m_vertices[vertex];
    return EdgeIterator(vertex, vertex_meta.end);
  }

  inline EdgeRange edges(VertexTopologyID node) {
    return EdgeRange(edge_begin(node), edge_end(node));
  }

  /*
   * Sort the outgoing edges for the given vertex.
   */
  void sortEdges(VertexTopologyID node) {
    auto& vertex_meta = m_vertices[node];
    std::sort(vertex_meta.begin, vertex_meta.end);
  }

  /*
   * Returns whether the given edge exists.
   *
   * Assumes sortEdges was already called for the vertex!
   */
  bool findEdgeSorted(VertexTopologyID src, VertexTopologyID dst) {
    auto const& vertex_meta = m_vertices[src];
    return std::binary_search(vertex_meta.begin, vertex_meta.end,
                              EdgeMetadata(dst));
  }

  template <bool sorted = false>
  void addEdges(VertexTopologyID src, const std::vector<VertexTopologyID> dsts,
                std::vector<EdgeData> data) {
    GALOIS_ASSERT(data.size() == dsts.size());
    this->addEdgesTopologyOnly<sorted>(src, dsts);
    for (size_t i = 0; i < dsts.size(); ++i) {
      m_edge_data[std::make_pair(src, dsts[i])] = data[i];
    }
  }

  /*
   * Adds outgoing edges from the given src to all dsts. If `sorted`, assume
   * both `dsts` and the existing edge array is sorted ascending, and maintain
   * sorted order.
   */
  template <bool sorted = false>
  void addEdgesTopologyOnly(VertexTopologyID src,
                            const std::vector<VertexTopologyID> dsts) {
    // Copies the edge list to the end of m_edit_log together with the new
    // edges.

    auto& vertex_meta = m_vertices[src];
    m_holes.fetch_add(vertex_meta.degree(), std::memory_order_relaxed);

    uint64_t const new_degree = vertex_meta.degree() + dsts.size();
    uint64_t const new_begin_offset =
        m_edges_tail.fetch_add(new_degree, std::memory_order_relaxed);
    EdgeMetadata* const new_begin = m_edit_log.begin() + new_begin_offset;
    EdgeMetadata* const new_end   = new_begin + new_degree;

    while (m_edit_log.size() < new_begin_offset + new_degree) {
      m_edges_lock.lock();
      {
        if (m_edit_log.size() < new_begin_offset + new_degree) {
          m_edit_log.resize(
              std::max(m_edit_log.size() * 2, new_begin_offset + new_degree));
        }
      }
      m_edges_lock.unlock();
    }

    EdgeMetadata* log_dst = new_begin;
    if constexpr (sorted) {
      std::merge(dsts.begin(), dsts.end(), vertex_meta.begin, vertex_meta.end,
                 log_dst);
    } else {
      // copy old edges
      log_dst = std::copy(vertex_meta.begin, vertex_meta.end, log_dst);

      // insert new edges
      std::copy(dsts.begin(), dsts.end(), log_dst);
    }

    // update vertex metadata
    vertex_meta.begin = new_begin;
    vertex_meta.end   = new_end;

    m_prefix_valid.store(false, std::memory_order_release);
  }

  // Performs the compaction algorithm by copying any vertices left in buffer 0
  // to buffer 1, then swapping the buffers.
  //
  // Not safe to call in parallel with insertions/deletions.
  void compact() {
    using std::swap;

    EdgeMetadata* const log_head = m_edit_log.begin();
    EdgeMetadata* const log_tail =
        log_head + m_edges_tail.load(std::memory_order_acquire);

    // move from buffer 0 to buffer 1
    galois::do_all(galois::iterate(vertices().begin(), vertices().end()),
                   [&](VertexTopologyID vertex_id) {
                     VertexMetadata& vertex_meta = m_vertices[vertex_id];

                     if (vertex_meta.begin < log_head ||
                         vertex_meta.begin >= log_tail) {
                       this->addEdgesTopologyOnly<false>(vertex_id, {});
                     }
                   });

    // At this point, there are no more live edges in buffer 0.
    m_edges_lock.lock();
    {
      m_edges.resize(0);
      swap(m_edges, m_edit_log);
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
          (m_edges.size() + m_edges_tail.load(std::memory_order_relaxed)) *
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

  /**
   * DO NOT USE WHILE MODIFYING THE GRAPH!
   * ONLY USE IF GRAPH HAS BEEN LOADED
   *
   * @param n Index into edge prefix sum
   * @returns The value that would be located at index n in an edge prefix sum
   * array
   */
  uint64_t operator[](uint64_t n) {
    if (!m_prefix_valid.load(std::memory_order_acquire))
      computePrefixSum();
    return m_pfx_sum_cache[n];
  }

  std::vector<uint64_t> const& getEdgePrefixSum() {
    if (!m_prefix_valid.load(std::memory_order_acquire))
      computePrefixSum();
    return m_pfx_sum_cache;
  }

public:
  class EdgeIterator
      : public boost::iterator_facade<EdgeIterator, EdgeHandle,
                                      boost::random_access_traversal_tag,
                                      EdgeHandle const> {
  private:
    VertexTopologyID const src;
    EdgeMetadata const* ptr;

    EdgeIterator(VertexTopologyID src, EdgeMetadata const* ptr)
        : src(src), ptr(ptr) {}

    void advance(std::ptrdiff_t n) { ptr += n; }

    std::ptrdiff_t distance_to(EdgeIterator const& y) const {
      return y.ptr - ptr;
    }

    void increment() { ++ptr; }

    void decrement() { --ptr; }

    bool equal(EdgeIterator const& other) const { return ptr == other.ptr; }

    EdgeHandle dereference() const { return EdgeHandle(src, *ptr); }

    friend class LS_LC_CSR_Graph;
    friend class boost::iterator_core_access;
  };
};

}; // namespace galois::graphs

#endif
