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

#include <iostream>
#include <queue>

#include "galois/Galois.h"
#include "galois/graphs/LS_LC_CSR_Graph.h"

template <typename GraphTy>
void check() {
  GraphTy g(4);

  g.addEdgesTopologyOnly(0, {1, 2});
  g.addEdgesTopologyOnly(1, {2, 3});
  g.addEdgesTopologyOnly(2, {3});

  auto print_graph = [&g](std::string_view msg) {
    std::cout << "- " << msg << " -" << std::endl;
    for (auto src : g.vertices()) {
      for (auto edge : g.edges(src)) {
        std::cout << src << "->" << g.getEdgeDst(edge) << std::endl;
      }
    }
  };

  print_graph("initial graph");

  g.addEdgesTopologyOnly(2, {1});

  print_graph("added 2->1");

  g.compact();

  print_graph("compacted");

  std::cout << std::endl << std::endl;
}

int main() {
  galois::SharedMemSys Galois_runtime;
  check<galois::graphs::LS_LC_CSR_Graph<void, void>>();
  check<galois::graphs::LS_LC_CSR_Graph<float, float>>();

  // check that we can access data on nodes/edges
  galois::graphs::LS_LC_CSR_Graph<uint32_t, uint32_t> g(4);

  g.setData(0, 0);
  GALOIS_ASSERT(g.getData(0) == 0);
  g.setData(1, 1);
  GALOIS_ASSERT(g.getData(1) == 1);
  g.setData(2, 2);
  GALOIS_ASSERT(g.getData(2) == 2);
  g.setData(3, 3);
  GALOIS_ASSERT(g.getData(3) == 3);

  uint64_t four = g.addVertices({4, 5, 6, 7});

  for (size_t ii = 0; ii < 4; ++ii) {
    // make sure previous data survived the resize
    GALOIS_ASSERT(g.getData(ii) == ii);
    // check the new vertex data
    GALOIS_ASSERT(g.getData(four + ii) == 4 + ii);
  }

  g.addEdges(0, {1, 2, 3}, {1, 2, 3});
  for (auto const& handle : g.edges(0)) {
    GALOIS_ASSERT(g.getEdgeDst(handle) == g.getEdgeData(handle));
  }

  uint64_t eight = g.addVertexTopologyOnly();
  GALOIS_ASSERT(eight == 8);

  g.addEdgesTopologyOnly(eight, {3, 2, 1, 0});
  GALOIS_ASSERT(g.getEdgeDst(*g.edge_begin(eight)) == 3);
  GALOIS_ASSERT(g.getEdgeDst(*(++g.edge_begin(eight))) == 2);
  GALOIS_ASSERT(g.getEdgeDst(*(++(++g.edge_begin(eight)))) == 1);
  GALOIS_ASSERT(g.getEdgeDst(*(++(++(++g.edge_begin(eight))))) == 0);
  GALOIS_ASSERT(g.getEdgeDst(*(--g.edge_end(eight))) == 0);

  g.sortEdges(eight);

  GALOIS_ASSERT(g.getEdgeDst(*g.edge_begin(eight)) == 0);
  GALOIS_ASSERT(g.getEdgeDst(*(++g.edge_begin(eight))) == 1);
  GALOIS_ASSERT(g.getEdgeDst(*(++(++g.edge_begin(eight)))) == 2);
  GALOIS_ASSERT(g.getEdgeDst(*(g.edge_begin(eight) + 2)) == 2);
  GALOIS_ASSERT(g.getEdgeDst(*(--g.edge_end(eight))) == 3);

  // check searching
  GALOIS_ASSERT(g.findEdgeSorted(eight, 0));
  GALOIS_ASSERT(g.findEdgeSorted(eight, 1));
  GALOIS_ASSERT(g.findEdgeSorted(eight, 2));
  GALOIS_ASSERT(g.findEdgeSorted(eight, 3));
  GALOIS_ASSERT(!g.findEdgeSorted(eight, 4));

  // check prefix sum
  GALOIS_ASSERT(g[0] == 3);
  GALOIS_ASSERT(g[1] == 3);
  GALOIS_ASSERT(g[2] == 3);
  // ...
  GALOIS_ASSERT(g[8] == 7);

  uint64_t num_vertices = (1 << 22) + 67;
  galois::graphs::LS_LC_CSR_Graph<void, void> big(num_vertices);
  for (uint64_t i = 0; i < num_vertices; ++i) {
    big.addEdgesTopologyOnly(i, {(i + 1) % num_vertices});
  }
  for (uint64_t i = 0; i < num_vertices; ++i) {
    GALOIS_ASSERT(big[i] == i + 1);
  }

  return 0;
}
