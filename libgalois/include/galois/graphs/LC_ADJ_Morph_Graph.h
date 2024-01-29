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
 * @file LC_ADJ_Morph_Graph.h
 *
 * Contains LC_ADJ_Morph_Graph and associated helpers.
 */

#pragma once

#include <algorithm>
#include <boost/iterator/counting_iterator.hpp>
#include <stdint.h>
#include <utility>
#include <vector>
#include <mutex>

namespace galois {
namespace graphs {

/**
 * A graph that can have new edges added to it.
 */
class LC_ADJ_Morph_Graph {
public:
  using iterator  = boost::counting_iterator<uint64_t>;
  using edge      = std::pair<uint64_t, uint64_t>;
  using edge_list = std::vector<edge>;

  using GraphNode = uint64_t;

private:
public:
  LC_ADJ_Morph_Graph(uint64_t num_nodes, uint64_t,
                     std::function<uint64_t(uint64_t)>,
                     std::function<uint64_t(uint64_t, uint64_t)>,
                     std::function<uint64_t(uint64_t, uint64_t)>)
      : edge_lists(num_nodes, edge_list()), lock() {}

  LC_ADJ_Morph_Graph(const LC_ADJ_Morph_Graph&)            = delete;
  LC_ADJ_Morph_Graph& operator=(const LC_ADJ_Morph_Graph&) = delete;

  void addEdge(uint64_t src, uint64_t dst) {
    auto guard = std::lock_guard(lock);
    edge_lists[src].emplace_back(src, dst);
  }

  void sortAllEdgesByDst(void) {
    auto guard = std::lock_guard(lock);
    for (auto node = begin(); node != end(); ++node) {
      auto& adj = edge_lists[*node];
      std::sort(adj.begin(), adj.end(),
                [](auto& a, auto& b) { return a.second < b.second; });
    }
  }

  size_t size(void) { return edge_lists.size(); };

  uint64_t getDegree(uint64_t node) { return edge_lists[node].size(); }

  iterator begin() { return iterator(0); }

  iterator end() { return iterator(this->size()); }

  const edge_list& out_edges(uint64_t node) { return edge_lists[node]; }

  constexpr inline uint64_t getEdgeDst(const edge& e) { return e.second; }

private:
  std::vector<edge_list> edge_lists;
  std::mutex lock;
};

} // namespace graphs
} // namespace galois
