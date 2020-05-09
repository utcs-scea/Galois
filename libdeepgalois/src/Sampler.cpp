#include "deepgalois/utils.h"
#include "deepgalois/Sampler.h"
#include "galois/Galois.h"
#include <time.h>
#include <vector>
#define PARALLEL_GEN

namespace deepgalois {
inline unsigned getDegree(Graph* g, index_t v) {
  // return g->get_degree(v);
  // return std::distance(g->edge_begin(v), g->edge_end(v));
  return g->edge_end(v) - g->edge_begin(v);
}

void Sampler::initializeMaskedGraph(size_t count, mask_t* masks, Graph* g) {
  this->count_ = count;
  this->masks_ = masks;
  // save original graph
  Sampler::graph = g;
  // allocate the object for the new masked graph
  Sampler::masked_graph = new Graph();

  std::vector<uint32_t> degrees(g->size(), 0);
  // get degrees of nodes that will be in new graph
  this->get_masked_degrees(g->size(), masks, g, degrees);
  auto offsets = deepgalois::parallel_prefix_sum(degrees);
  size_t ne    = offsets[g->size()];

  // save ids (on original graph) of training nodes to vector
  for (size_t i = 0; i < g->size(); i++) {
    if (masks[i] == 1)
      Sampler::node_train.push_back(i);
  }

  Sampler::masked_graph->allocateFrom(g->size(), ne);
  Sampler::masked_graph->constructNodes();
  // same as original graph, except keep only edges involved in masks
  galois::do_all(
      galois::iterate((size_t)0, g->size()),
      [&](const auto src) {
        Sampler::masked_graph->fixEndEdge(src, offsets[src + 1]);
        if (masks[src] == 1) {
          auto idx = offsets[src];
          for (auto e = g->edge_begin(src); e != g->edge_end(src); e++) {
            const auto dst = g->getEdgeDst(e);
            if (masks[dst] == 1)
              Sampler::masked_graph->constructEdge(idx++, dst, 0);
          }
        }
      },
      galois::loopname("gen_subgraph"));

  Sampler::masked_graph->degree_counting();
  Sampler::avg_deg  = masked_graph->sizeEdges() / masked_graph->size();
  Sampler::subg_deg = (avg_deg > SAMPLE_CLIP) ? SAMPLE_CLIP : avg_deg;

  // size_t idx = 0;
  // vertices_.resize(count);
  // for (size_t i = begin; i < end; i++) {
  //  if (masks_[i] == 1)
  //    vertices_[idx++] = i;
  //}
}

//! determine degree of each vertex in a masked graph (given by masks and g)
void Sampler::get_masked_degrees(size_t n, mask_t* masks, Graph* g,
                                 std::vector<uint32_t>& degrees) {
  assert(degrees.size() == n);
#ifdef PARALLEL_GEN
  galois::do_all(
      galois::iterate(size_t(0), n),
      [&](const auto src) {
#else
  for (size_t src = 0; src < n; src++) {
#endif
        if (masks[src] == 1) {
          for (auto e = g->edge_begin(src); e != g->edge_end(src); e++) {
            const auto dst = g->getEdgeDst(e);
            if (masks[dst] == 1)
              degrees[src]++;
          }
        }
      }
#ifdef PARALLEL_GEN
      ,
      galois::loopname("update_degrees"));
#endif
}

//! returns a graph in the variable sub: it is g with the mask applied
void Sampler::getMaskedGraph(size_t n, mask_t* masks, Graph* g, Graph& sub) {
  std::vector<uint32_t> degrees(n, 0);
  this->get_masked_degrees(n, masks, g, degrees);
  // auto offsets = deepgalois::parallel_prefix_sum(degrees);
  auto offsets = deepgalois::prefix_sum(degrees);
  size_t ne    = offsets[n];
  // galois::gPrint("Generate masked graph: num_vertices=", n, ", num_edges=",
  // ne, "\n");
  //

  // note this constructs the full graph's nodes; just trims edges
  sub.allocateFrom(n, ne);
  sub.constructNodes();

#ifdef PARALLEL_GEN
  galois::do_all(
      galois::iterate((size_t)0, n),
      [&](const auto src) {
#else
  for (size_t src = 0; src < n; src++) {
#endif
        sub.fixEndEdge(src, offsets[src + 1]);
        if (masks[src] == 1) {
          auto idx = offsets[src];
          for (auto e = g->edge_begin(src); e != g->edge_end(src); e++) {
            const auto dst = g->getEdgeDst(e);
            if (masks[dst] == 1)
              sub.constructEdge(idx++, dst, 0);
          }
        }
      }
#ifdef PARALLEL_GEN
      ,
      galois::loopname("gen_subgraph"));
#endif
}

// helper function for graph saint implementation below
void Sampler::checkGSDB(std::vector<db_t>& DB0, std::vector<db_t>& DB1,
                        std::vector<db_t>& DB2, size_t size) {
  if (DB0.capacity() < size) {
    DB0.reserve(DB0.capacity() * 2);
    DB1.reserve(DB1.capacity() * 2);
    DB2.reserve(DB2.capacity() * 2);
  }
  DB0.resize(size);
  DB1.resize(size);
  DB2.resize(size);
}

//! debug function: prints out sets of vertices
void print_vertex_set(VertexSet vertex_set) {
  unsigned counter = 0;
  unsigned n       = vertex_set.size();
  galois::gPrint("( ");
  for (int i : vertex_set) {
    counter++;
    if (counter > 16 && counter < n - 16)
      continue;
    galois::gPrint(i, " ");
  }
  galois::gPrint(")\n");
}

// implementation from GraphSAINT
// https://github.com/GraphSAINT/GraphSAINT/blob/master/ipdps19_cpp/sample.cpp
void Sampler::select_vertices(size_t n, int m, VertexSet& st, unsigned seed) {
  unsigned myseed = seed;

  // unsigned myseed = tid;
  // DBx: Dashboard line x, IAx: Index array line x
  std::vector<db_t> DB0, DB1, DB2, IA0, IA1, IA2, IA3, IA4, nDB0, nDB1, nDB2;
  DB0.reserve(subg_deg * m * ETA);
  DB1.reserve(subg_deg * m * ETA);
  DB2.reserve(subg_deg * m * ETA);
  IA0.reserve(n);
  IA1.reserve(n);
  IA2.reserve(n);
  IA3.reserve(n);
  IA4.reserve(n);
  IA0.resize(m);
  IA1.resize(m);
  IA2.resize(m);
  IA3.resize(m);

  // galois::gPrint("seed ", myseed, " m ", m, "\n");
  // galois::gPrint("node_train size: ", node_train.size(), "\n");
  // printf("( ");
  // for (size_t i = 0; i < 10; i++) std::cout << node_train[i] << " ";
  // printf(")\n");

  for (int i = 0; i < m; i++) {
    auto rand_idx = rand_r(&myseed) % Sampler::node_train.size();
    db_t v = IA3[i] = Sampler::node_train[rand_idx];
    st.insert(v);
    IA0[i] = getDegree(Sampler::masked_graph, v);
    IA0[i] = (IA0[i] > SAMPLE_CLIP) ? SAMPLE_CLIP : IA0[i];
    IA1[i] = 1;
    IA2[i] = 0;
  }
  // calculate prefix sum for IA0 and store in IA2 to compute the address for
  // each frontier in DB
  IA2[0] = IA0[0];
  for (int i = 1; i < m; i++)
    IA2[i] = IA2[i - 1] + IA0[i];
  // now fill DB accordingly
  checkGSDB(DB0, DB1, DB2, IA2[m - 1]);
  for (int i = 0; i < m; i++) {
    db_t DB_start = (i == 0) ? 0 : IA2[i - 1];
    db_t DB_end   = IA2[i];
    for (auto j = DB_start; j < DB_end; j++) {
      DB0[j] = IA3[i];
      DB1[j] = (j == DB_start) ? (j - DB_end) : (j - DB_start);
      DB2[j] = i + 1;
    }
  }

  db_t choose, neigh_v, newsize, tmp;
  for (size_t itr = 0; itr < n - m; itr++) {
    choose = db_t(-1);
    while (choose == db_t(-1)) {
      tmp = rand_r(&myseed) % DB0.size();
      if (size_t(tmp) < DB0.size())
        if (DB0[tmp] != db_t(-1))
          choose = tmp;
    }
    choose      = (DB1[choose] < 0) ? choose : (choose - DB1[choose]);
    db_t v      = DB0[choose];
    auto degree = getDegree(Sampler::masked_graph, v);
    neigh_v     = (degree != 0) ? rand_r(&myseed) % degree : db_t(-1);
    if (neigh_v != db_t(-1)) {
      neigh_v = Sampler::masked_graph->getEdgeDst(
          Sampler::masked_graph->edge_begin(v) + neigh_v);
      st.insert(neigh_v);
      IA1[DB2[choose] - 1] = 0;
      IA0[DB2[choose] - 1] = 0;
      for (auto i = choose; i < choose - DB1[choose]; i++)
        DB0[i] = db_t(-1);
      newsize = getDegree(Sampler::masked_graph, neigh_v);
      newsize = (newsize > SAMPLE_CLIP) ? SAMPLE_CLIP : newsize;
    } else
      newsize = 0;
    // shrink DB to remove sampled nodes, also shrink IA accordingly
    bool cond = DB0.size() + newsize > DB0.capacity();
    if (cond) {
      // compute prefix sum for the location in shrinked DB
      IA4.resize(IA0.size());
      IA4[0] = IA0[0];
      for (size_t i = 1; i < IA0.size(); i++)
        IA4[i] = IA4[i - 1] + IA0[i];
      nDB0.resize(IA4.back());
      nDB1.resize(IA4.back());
      nDB2.resize(IA4.back());
      IA2.assign(IA4.begin(), IA4.end());
      for (size_t i = 0; i < IA0.size(); i++) {
        if (IA1[i] == 0)
          continue;
        db_t DB_start = (i == 0) ? 0 : IA4[i - 1];
        db_t DB_end   = IA4[i];
        for (auto j = DB_start; j < DB_end; j++) {
          nDB0[j] = IA3[i];
          nDB1[j] = (j == DB_start) ? (j - DB_end) : (j - DB_start);
          nDB2[j] = i + 1;
        }
      }
      // remap the index in DB2 by compute prefix of IA1 (new idx in IA)
      IA4.resize(IA1.size());
      IA4[0] = IA1[0];
      for (size_t i = 1; i < IA1.size(); i++)
        IA4[i] = IA4[i - 1] + IA1[i];
      DB0.assign(nDB0.begin(), nDB0.end());
      DB1.assign(nDB1.begin(), nDB1.end());
      DB2.assign(nDB2.begin(), nDB2.end());
      for (auto i = DB2.begin(); i < DB2.end(); i++)
        *i = IA4[*i - 1];
      db_t curr = 0;
      for (size_t i = 0; i < IA0.size(); i++) {
        if (IA0[i] != 0) {
          IA0[curr] = IA0[i];
          IA1[curr] = IA1[i];
          IA2[curr] = IA2[i];
          IA3[curr] = IA3[i];
          curr++;
        }
      }
      IA0.resize(curr);
      IA1.resize(curr);
      IA2.resize(curr);
      IA3.resize(curr);
    }
    checkGSDB(DB0, DB1, DB2, newsize + DB0.size());
    IA0.push_back(newsize);
    IA1.push_back(1);
    IA2.push_back(IA2.back() + IA0.back());
    IA3.push_back(neigh_v);
    db_t DB_start = (*(IA2.end() - 2));
    db_t DB_end   = IA2.back();
    for (auto j = DB_start; j < DB_end; j++) {
      DB0[j] = IA3.back();
      DB1[j] = (j == DB_start) ? (j - DB_end) : (j - DB_start);
      DB2[j] = IA3.size();
    }
  }
  // galois::gPrint("Done selection, vertex_set size: ", st.size(), ", set: ");
  // print_vertex_set(st);
}

// API function for user-defined selection strategy
// Select n vertices from vertices and put them in vertex_set.
// nv: number of vertices in the original graph;
// n: number of vertices in the subgraph;
// m: number of vertices in the frontier.
// our implementation of GraphSAINT sampling
void Sampler::select_vertices(size_t nv, size_t n, int m, Graph* g,
                              VertexList vertices, VertexSet& vertex_set) {
  // galois::gPrint("Select a vertex set of size ", n, " from ", nv, " vertices,
  // graph size: ", g->size(), "\n");
  assert(nv == vertices.size());
  auto frontier_indices = deepgalois::select_k_items(
      m, 0, (int)nv); // randomly select m vertices from vertices as frontier
  VertexList frontier(m);
  for (int i = 0; i < m; i++)
    frontier[i] = vertices[frontier_indices[i]];
  vertex_set.insert(frontier.begin(), frontier.end());
  // galois::gPrint("vertex_set size: ", vertex_set.size(), "\n");
  int* degrees = new int[m];
  for (int i = 0; i < m; i++) {
    // galois::do_all(galois::iterate(size_t(0), size_t(m)), [&](const auto i) {
    degrees[i] = (int)getDegree(g, frontier[i]);
  } //, galois::loopname("compute_degrees"));
  for (size_t i = 0; i < n - m; i++) {
    auto pos    = select_one_item((int)m, degrees);
    auto u      = frontier[pos];
    auto degree = degrees[pos];
    int j       = 0;
    for (; j < degree; j++) {
      auto neighbor_id = rand() % degree; // randomly select a neighbor
      auto dst         = g->getEdgeDst(g->edge_begin(u) + neighbor_id);
      if (vertex_set.find(dst) == vertex_set.end()) {
        frontier[pos] = dst;
        degrees[pos]  = getDegree(g, frontier[pos]);
        vertex_set.insert(dst);
        break;
      }
    }
    if (j == degree)
      galois::gPrint("Not found from ", degree, " neighbors\n");
  }
  /*
  assert(n == vertex_set.size()); // size of vertex_set could be slightly
  smaller than n galois::gPrint("Done selection, vertex_set size: ",
  vertex_set.size(), ", set: "); print_vertex_set(vertex_set);
  */
}

void Sampler::getMasks(size_t n, VertexSet vertices, mask_t* masks) {
  // galois::gPrint("Updating masks, size = ", vertices.size(), "\n");
  std::fill(masks, masks + n, 0);
  for (auto v : vertices)
    masks[v] = 1;
}

inline VertexList Sampler::reindexVertices(size_t n, VertexSet vertex_set) {
  VertexList new_ids(n, 0);
  int vid = 0;
  for (auto v : vertex_set) {
    new_ids[v] = vid++; // reindex
  }
  return new_ids;
}

// Given a subset of vertices and a graph g, generate a subgraph sg from the
// graph g
void Sampler::reindexSubgraph(VertexSet& keptVertices, Graph& origGraph,
                              Graph& reindexGraph) {
  // auto n = origGraph.size(); // old graph size
  auto nv            = keptVertices.size(); // new graph (subgraph) size
  VertexList new_ids = this->reindexVertices(graph->size(), keptVertices);
  std::vector<uint32_t> degrees(nv, 0); // degrees of vertices in the subgraph
  for (auto v : keptVertices) {
    degrees[new_ids[v]] = getDegree(&origGraph, v);
  }
  // auto offsets = deepgalois::parallel_prefix_sum(degrees);
  auto offsets = deepgalois::prefix_sum(degrees);
  auto ne      = offsets[nv];
  // galois::gPrint("Generate subgraph: num_vertices=", nv, ", num_edges=", ne,
  // "\n");
  reindexGraph.allocateFrom(nv, ne);
  reindexGraph.constructNodes();
  VertexList old_ids(keptVertices.begin(),
                     keptVertices.end()); // vertex ID mapping
#ifdef PARALLEL_GEN
  galois::do_all(
      galois::iterate((size_t)0, nv),
      [&](const auto i) {
#else
  for (size_t i = 0; i < nv; i++) {
#endif
        reindexGraph.fixEndEdge(i, offsets[i + 1]);
        unsigned j  = 0;
        auto old_id = old_ids[i];
        for (auto e = origGraph.edge_begin(old_id);
             e != origGraph.edge_end(old_id); e++) {
          auto dst = new_ids[origGraph.getEdgeDst(e)];
          assert(dst < nv);
          reindexGraph.constructEdge(offsets[i] + j, dst, 0);
          j++;
        }
      }
#ifdef PARALLEL_GEN
      ,
      galois::loopname("construct_graph"));
#endif
}

void Sampler::subgraph_sample(size_t n, Graph& sg, mask_t* masks,
                              unsigned tid) {
  VertexSet sampledSet;
  // n = 9000 by default
  // this->select_vertices(count_, n, m_, masked_graph, vertices_, sampledSet);

  // do the sampling of vertices from training set + using masked graph
  this->select_vertices(n, m_, sampledSet, tid); // m = 1000 by default

  // create the masks on the masked_graph
  getMasks(Sampler::graph->size(), sampledSet, masks);

  Graph masked_sg;
  this->getMaskedGraph(
      Sampler::graph->size(), masks, Sampler::masked_graph,
      masked_sg); // remove edges whose destination is not masked
  this->reindexSubgraph(sampledSet, masked_sg, sg);
}

} // namespace deepgalois