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

#ifndef GALOIS_GRAPHS_LS_LC_CSR_GRAPH_H
#define GALOIS_GRAPHS_LS_LC_CSR_GRAPH_H

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

  // todo: should we use a galois spinlock here instead of a mutex?
  using EdgeDataMap = phmap::parallel_flat_hash_map_m<EdgeHandle, EdgeData>;

  using EdgeDataStore =
      std::conditional_t<HasEdgeData, EdgeDataMap, std::tuple<>>;

  // forward-declarations of internal structs
  struct VertexMetadata;
  using EdgeMetadata = VertexTopologyID;

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
  // todo(meyer): there is currently a memory leak in prefix sum :/
  PrefixSum<VertexMetadata, uint64_t, transmute, scan_op, combiner,
            CacheLinePaddedArr>
      m_pfx{&m_vertices[0], &m_pfx_sum_cache[0]};

  bool m_prefix_valid;

  void resetPrefixSum() {
    m_pfx_sum_cache.resize(m_vertices.size());
    m_pfx.src      = &m_vertices[0];
    m_pfx.dst      = &m_pfx_sum_cache[0];
    m_prefix_valid = false;
  }

  // Compute the prefix sum using the two level method
  void computePrefixSum() {
    constexpr uint64_t PARALLEL_PREFIX_SUM_VERTEX_THRESHOLD =
        static_cast<uint64_t>(1) << 30;
    if (m_vertices.size() > PARALLEL_PREFIX_SUM_VERTEX_THRESHOLD) {
      m_pfx.computePrefixSumSerially(m_vertices.size());
    } else {
      m_pfx.computePrefixSum(m_vertices.size());
    }
    m_prefix_valid = true;
  }

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
  inline void setEdgeData(EdgeHandle handle, E data) {
    m_edge_data.lazy_emplace_l(
        handle, [&](auto& v) { v.second = data; },
        [&](auto const& cons) { cons(handle, data); });
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

  /**
   * Adds multiple vertices to the graph. The new vertices will be assigned
   * consecutive topology IDs, and the lowest new ID is returned.
   */
  VertexTopologyID addVerticesTopologyOnly(size_t count) {
    VertexTopologyID const start = m_vertices.size();
    m_vertices.resize(start + count);
    return start;
  }

  // Adds multiple vertices to the graph. The new vertices will be assigned
  // consecutive topology IDs, and the lowest new ID is returned.
  template <typename V = VertexData, typename = std::enable_if<HasVertexData>>
  VertexTopologyID addVertices(std::vector<V> data) {
    auto const start = addVerticesTopologyOnly(data.size());
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

  template <typename E = EdgeData, typename = std::enable_if<HasEdgeData>>
  inline E& getEdgeData(EdgeIterator const& it) {
    return m_edge_data[*it];
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
    return EdgeIterator(
        vertex, &getEdgeMetadata(vertex_meta.buffer, vertex_meta.begin));
  }

  inline EdgeIterator edge_end(VertexTopologyID vertex) {
    auto const& vertex_meta = m_vertices[vertex];
    return EdgeIterator(vertex,
                        &getEdgeMetadata(vertex_meta.buffer, vertex_meta.end));
  }

  inline EdgeRange edges(VertexTopologyID node) {
    return EdgeRange(edge_begin(node), edge_end(node));
  }

  /*
   * Iterates over the outgoing edges, calling the callback with the
   * VertexTopologyID of each edge.
   */
  template <typename Callback>
  void for_each_edge(VertexTopologyID vertex, Callback const& callback) {
    auto const& vertex_meta = m_vertices[vertex];
    EdgeMetadata const* begin =
        &getEdgeMetadata(vertex_meta.buffer, vertex_meta.begin);
    for (uint64_t i = 0; i < vertex_meta.degree(); ++i) {
      callback(static_cast<VertexTopologyID>(*begin++));
    }
  }

  /*
   * Sort the outgoing edges for the given vertex.
   */
  void sortEdges(VertexTopologyID node) {
    auto& vertex_meta = m_vertices[node];

    EdgeMetadata* start =
        &getEdgeMetadata(vertex_meta.buffer, vertex_meta.begin);

    EdgeMetadata* end = &getEdgeMetadata(vertex_meta.buffer, vertex_meta.end);

    std::sort(start, end);
  }

  /*
   * Returns whether the given edge exists.
   *
   * Assumes sortEdges was already called for the vertex!
   */
  bool findEdgeSorted(VertexTopologyID src, VertexTopologyID dst) {
    auto const& vertex_meta = m_vertices[src];
    EdgeMetadata* start =
        &getEdgeMetadata(vertex_meta.buffer, vertex_meta.begin);
    EdgeMetadata* end = &getEdgeMetadata(vertex_meta.buffer, vertex_meta.end);

    return std::binary_search(start, end, EdgeMetadata(dst));
  }

  template <bool sorted = false>
  void addEdges(VertexTopologyID src, const std::vector<VertexTopologyID> dsts,
                std::vector<EdgeData> data) {
    GALOIS_ASSERT(data.size() == dsts.size());
    this->addEdgesTopologyOnly<sorted>(src, dsts);
    for (size_t i = 0; i < dsts.size(); ++i) {
      auto key = std::make_pair(src, dsts[i]);
      setEdgeData(key, data[i]);
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
    // Copies the edge list to the end of m_edges[1] together with the new
    // edges.

    auto& vertex_meta = m_vertices[src];
    if (vertex_meta.buffer == 1)
      m_holes.fetch_add(vertex_meta.degree(), std::memory_order_relaxed);

    uint64_t const new_degree = vertex_meta.degree() + dsts.size();
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

    EdgeMetadata* log_dst = &getEdgeMetadata(1, new_begin);
    if constexpr (sorted) {
      std::merge(dsts.begin(), dsts.end(),
                 &getEdgeMetadata(vertex_meta.buffer, vertex_meta.begin),
                 &getEdgeMetadata(vertex_meta.buffer, vertex_meta.end),
                 log_dst);
    } else {
      // copy old edges
      log_dst = std::copy(
          &getEdgeMetadata(vertex_meta.buffer, vertex_meta.begin),
          &getEdgeMetadata(vertex_meta.buffer, vertex_meta.end), log_dst);

      // insert new edges
      std::copy(dsts.begin(), dsts.end(), log_dst);
    }

    // update vertex metadata
    vertex_meta.buffer = 1;
    vertex_meta.begin  = new_begin;
    vertex_meta.end    = new_end;

    m_prefix_valid = false;
  }

  // Performs the compaction algorithm by copying any vertices left in buffer 0
  // to buffer 1, then swapping the buffers.
  //
  // Not safe to call in parallel with insertions/deletions.
  void compact() {
    using std::swap;

    if (m_vertices.empty())
      return;

    auto const& prefix_sum = getEdgePrefixSum();
    std::vector<VertexMetadata> new_vertices(m_vertices.size());
    new_vertices[0].buffer = 0;
    new_vertices[0].begin  = 0;
    new_vertices[0].end    = prefix_sum[0];
    galois::do_all(galois::iterate(1ul, m_vertices.size()),
                   [&](VertexTopologyID vertex_id) {
                     new_vertices[vertex_id].buffer = 0;
                     new_vertices[vertex_id].begin  = prefix_sum[vertex_id - 1];
                     new_vertices[vertex_id].end    = prefix_sum[vertex_id];
                   });

    LargeVector<EdgeMetadata> new_edges(prefix_sum.back());

    galois::do_all(
        galois::iterate(0ul, m_vertices.size()),
        [&](VertexTopologyID vertex_id) {
          auto const& m_vertex_meta   = m_vertices[vertex_id];
          auto const& new_vertex_meta = new_vertices[vertex_id];
          std::copy(&getEdgeMetadata(m_vertex_meta.buffer, m_vertex_meta.begin),
                    &getEdgeMetadata(m_vertex_meta.buffer, m_vertex_meta.end),
                    &new_edges[new_vertex_meta.begin]);
        },
        galois::steal());

    swap(m_vertices, new_vertices);
    swap(m_edges[0], new_edges);
    m_edges[1].resize(0);

    m_edges_tail   = 0;
    m_holes        = 0;
    m_prefix_valid = false;
  }

  /*
    Compaction policy utilities.
  */

  inline size_t getCSRMemoryUsageBytes() {
    return m_edges[0].size() * sizeof(EdgeMetadata);
  }

  // Returns the number of bytes used for the log.
  inline size_t getLogMemoryUsageBytes() {
    return m_edges_tail.load(std::memory_order_relaxed) * sizeof(EdgeMetadata);
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
    if (!m_prefix_valid)
      computePrefixSum();
    return m_pfx_sum_cache[n];
  }

  std::vector<uint64_t> const& getEdgePrefixSum() {
    if (!m_prefix_valid)
      computePrefixSum();
    return m_pfx_sum_cache;
  }

private:
  struct VertexMetadata {
    uint64_t begin; // inclusive
    uint64_t end;   // exclusive
    uint8_t buffer;

    VertexMetadata() : begin(0), end(0), buffer(0) {}

    VertexMetadata(VertexMetadata const& other)
        : begin(other.begin), end(other.end), buffer(other.buffer) {}

    VertexMetadata(VertexMetadata&& other)
        : begin(std::move(other.begin)), end(std::move(other.end)),
          buffer(std::move(other.buffer)) {}

    inline uint64_t degree() const { return end - begin; }
  };

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
