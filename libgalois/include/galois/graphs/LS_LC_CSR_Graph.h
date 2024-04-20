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
#include <vector>
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
  LargeVector<VertexMetadata> m_vertices;

  // m_edges[0] is the CSR with gaps, m_edges[1] is the update log.
  LargeVector<EdgeMetadata> m_edges[2];
  SpinLock m_edges_lock; // guards resizing of edges vectors
  EdgeDataStore m_edge_data;
  alignas(hardware_destructive_interference_size) std::atomic_uint64_t
      m_edges_tail = ATOMIC_VAR_INIT(0);

  /*
   * Prefix Sum utilities
   */
  static constexpr uint64_t PARALLEL_PREFIX_SUM_VERTEX_THRESHOLD = 1ul << 25;
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

  static uint64_t transmute_b0(VertexMetadata const& vertex_meta) {
    return (vertex_meta.buffer) ? 0 : vertex_meta.degree();
  }
  static uint64_t scan_op_b0(VertexMetadata const& p, const uint64_t& l) {
    return (p.buffer) ? l : p.degree() + l;
  }
  PrefixSum<VertexMetadata, uint64_t, transmute_b0, scan_op_b0, combiner,
            CacheLinePaddedArr>
      m_pfx_b0{&m_vertices[0], &m_pfx_sum_cache[0]};

  // todo(meyer): there is currently a memory leak in prefix sum :/

  bool m_prefix_valid;

  void resetPrefixSum() {
    m_pfx_sum_cache.resize(m_vertices.size());
    m_pfx.src      = &m_vertices[0];
    m_pfx.dst      = &m_pfx_sum_cache[0];
    m_pfx_b0.src   = &m_vertices[0];
    m_pfx_b0.dst   = &m_pfx_sum_cache[0];
    m_prefix_valid = false;
  }

  // Compute the prefix sum using the two level method
  void computePrefixSum() {
    resetPrefixSum();
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
  LS_LC_CSR_Graph(uint64_t num_vertices) : m_vertices(num_vertices) {
    galois::do_all(galois::iterate(0ul, num_vertices),
                   [&](VertexTopologyID const& vertex) {
                     m_vertices[vertex].buffer = 0;
                     m_vertices[vertex].begin  = 0;
                     m_vertices[vertex].end    = 0;
                   });
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
    galois::do_all(galois::iterate(start, start + count),
                   [&](VertexTopologyID const& vertex) {
                     m_vertices[vertex].buffer = 0;
                     m_vertices[vertex].begin  = 0;
                     m_vertices[vertex].end    = 0;
                   });
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
    GALOIS_ASSERT(data.size() == dsts.size(), "Data size mismatch");
    this->addEdgesTopologyOnly<sorted>(src, dsts);
    for (size_t i = 0; i < dsts.size(); ++i) {
      auto key = std::make_pair(src, dsts[i]);
      setEdgeData(key, data[i]);
    }
  }

  template <bool sorted = false>
  void addBatchTopologyOnly(
      std::vector<std::pair<VertexTopologyID, std::vector<VertexTopologyID>>>
          edges) {
    if (edges.empty())
      return;

    static constexpr size_t PARALLEL_VERTEX_SORT_THRESHOLD = (1ul << 16);
    if (edges.size() >= PARALLEL_VERTEX_SORT_THRESHOLD) {
      galois::ParallelSTL::sort(
          edges.begin(), edges.end(),
          [](auto const& a, auto const& b) { return a.first < b.first; });
    } else {
      std::sort(edges.begin(), edges.end(),
                [](auto const& a, auto const& b) { return a.first < b.first; });
    }

    std::vector<uint64_t> pfx_sum(edges.size());
    galois::GReduceMax<uint64_t> max_vertex_id;
    max_vertex_id.reset();
    galois::GAccumulator<uint64_t> old_degree_total;
    old_degree_total.reset();
    galois::do_all(
        galois::iterate(0ul, edges.size()),
        [&](size_t idx) {
          auto const vertex_id = edges[idx].first;
          max_vertex_id.update(vertex_id);
          auto const old_degree =
              vertex_id < m_vertices.size() ? getDegree(vertex_id) : 0;
          old_degree_total += old_degree;
          pfx_sum[idx] = old_degree + edges[idx].second.size();
        },
        galois::loopname("ComputeVertexDegrees"));

    uint64_t const prev_num_vertices = m_vertices.size();
    m_vertices.resize(std::max(max_vertex_id.reduce() + 1, m_vertices.size()));
    galois::do_all(galois::iterate(m_vertices.begin() + prev_num_vertices,
                                   m_vertices.end()),
                   [&](VertexMetadata& vertex_meta) {
                     vertex_meta.buffer = 0;
                     vertex_meta.begin  = 0;
                     vertex_meta.end    = 0;
                   });

    for (size_t i = 1; i < pfx_sum.size(); ++i)
      pfx_sum[i] += pfx_sum[i - 1];
    auto const num_new_edges = pfx_sum.back();
    std::cout << "old degrees total: " << old_degree_total.reduce()
              << ", new degrees total: " << num_new_edges << std::endl;

    auto const start =
        m_edges_tail.fetch_add(num_new_edges, std::memory_order_relaxed);
    if (m_edges[1].size() < start + num_new_edges)
      m_edges[1].resize(start + num_new_edges);

    galois::do_all(
        galois::iterate(0ul, edges.size()),
        [&](size_t idx) {
          auto const& [src, dsts] = edges[idx];
          auto const new_begin    = (idx) ? (start + pfx_sum[idx - 1]) : start;
          auto const new_end      = start + pfx_sum[idx];
          auto& vertex_meta       = m_vertices[src];
          EdgeMetadata* log_dst   = &getEdgeMetadata(1, new_begin);
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
        },
        galois::steal(), galois::loopname("CopyEdgesToLog"));

    m_prefix_valid = false;
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

    auto& vertex_meta         = m_vertices[src];
    uint64_t const new_degree = vertex_meta.degree() + dsts.size();
    uint64_t const new_begin =
        m_edges_tail.fetch_add(new_degree, std::memory_order_relaxed);
    uint64_t const new_end = new_begin + new_degree;

    if (m_edges[1].size() < new_end) {
      m_edges_lock.lock();
      {
        if (m_edges[1].size() < new_end)
          m_edges[1].resize(new_end);
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

public:
  // Performs the compaction algorithm by copying any vertices left in buffer 0
  // to buffer 1, then swapping the buffers.
  //
  // Not safe to call in parallel with insertions/deletions.
  void compact() {
    using std::swap;

    if (m_vertices.empty())
      return;
    auto const num_vertices = m_vertices.size();

    // step 1: copy from CSR to log
    {
      resetPrefixSum();
      if (m_vertices.size() > PARALLEL_PREFIX_SUM_VERTEX_THRESHOLD) {
        m_pfx_b0.computePrefixSum(m_vertices.size());
      } else {
        m_pfx_b0.computePrefixSumSerially(m_vertices.size());
      }
      auto const log_num_edges = m_pfx_sum_cache.back();
      auto const start =
          m_edges_tail.fetch_add(log_num_edges, std::memory_order_relaxed);
      if (m_edges[1].size() < start + log_num_edges) {
        m_edges[1].resize(start + log_num_edges);
      }

      galois::do_all(
          galois::iterate(0ul, num_vertices),
          [&](VertexTopologyID ii) {
            auto& vertex_meta = m_vertices[ii];
            if (vertex_meta.buffer)
              return; // already on the log
            auto const new_begin = start + ((ii) ? m_pfx_sum_cache[ii - 1] : 0);
            auto const new_end   = start + m_pfx_sum_cache[ii];
            std::copy(&getEdgeMetadata(0, vertex_meta.begin),
                      &getEdgeMetadata(0, vertex_meta.end),
                      &getEdgeMetadata(1, new_begin));
            // vertex_meta.buffer = 1;
            vertex_meta.buffer = 0; // not accurate, but it soon will be...
            vertex_meta.begin  = new_begin;
            vertex_meta.end    = new_end;
          },
          galois::steal(), galois::loopname("CopyCSRToLog"));
    }
    // At this point, all edges are on the log (and the prefix sum is invalid).
    // We can now compact into the CSR.
    {
      // compute the actual prefix sum
      auto const& prefix_sum = getEdgePrefixSum();
      m_edges[0].resize(prefix_sum.back());

      galois::do_all(
          galois::iterate(0ul, m_vertices.size()),
          [&](size_t idx) {
            auto& vertex_meta    = m_vertices[idx];
            auto const new_begin = (idx) ? prefix_sum[idx - 1] : 0;
            auto const new_end   = prefix_sum[idx];
            std::copy(&getEdgeMetadata(1, vertex_meta.begin),
                      &getEdgeMetadata(1, vertex_meta.end),
                      &getEdgeMetadata(0, new_begin));
            // vertex_meta.buffer = 0; // already done above
            vertex_meta.begin = new_begin;
            vertex_meta.end   = new_end;
          },
          galois::steal(), galois::loopname("CompactLogToCSR"));

      swap(m_edges[0], m_edges[1]);
      // m_edges[1].resize(0);
      m_edges_tail   = 0;
      m_prefix_valid = false;
    }
  }

  /*
    Compaction policy utilities.
  */

  inline size_t getCSRMemoryUsageBytes() {
    return m_edges[0].size() * sizeof(EdgeMetadata);
  }

  // Returns the number of bytes used for the log.
  inline size_t getLogMemoryUsageBytes() {
    return m_edges_tail * sizeof(EdgeMetadata);
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

    VertexMetadata() = delete;

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

} // namespace galois::graphs

#endif
