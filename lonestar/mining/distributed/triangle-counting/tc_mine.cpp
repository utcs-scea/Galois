// #include "MiningBench/Start.h"
#include "pangolin/BfsMining/vertex_miner.h"
#define TRIANGLE

const char* name = "TC";
const char* desc =
    "Counts the triangles in a graph (inputs do NOT need to be symmetrized)";
const char* url = nullptr;

// #include "pangolin/BfsMining/vertex_miner_api.h"
// class MyAPI : public VertexMinerAPI<BaseEmbedding> {
// public:
//   // toExtend (only extend the last vertex in the embedding)
//   static bool toExtend(unsigned n, const BaseEmbedding&, unsigned pos) {
//     return pos == n - 1;
//   }
//   // toAdd (only add vertex connected to all the vertices in the embedding)
//   static bool toAdd(unsigned, PangolinGraph&, const BaseEmbedding&, unsigned,
//                     VertexId) {
//     return true;
//   }
// };

// class AppMiner : public VertexMiner<SimpleElement, BaseEmbedding, MyAPI, true> {
// public:
//   AppMiner(unsigned ms, int nt)
//       : VertexMiner<SimpleElement, BaseEmbedding, MyAPI, true>(ms, nt, 0) {
//     if (ms != 3) {
//       printf("ERROR: command line argument k must be 3\n");
//       exit(1);
//     }
//     set_num_patterns(1);
//   }
//   ~AppMiner() {}
//   void print_output() {
//     std::cout << "\n\ttotal_num_triangles = " << get_total_count() << "\n";
//   }

//   void tc_vertex_solver() {  // vertex parallel
//         galois::do_all(
//             galois::iterate(this->graph.masterNodesRange()),
//             [&](const GNode& src) {
//                 for (auto e : this->graph.edges(src)) {
//                     auto dst = this->graph.getEdgeDst(e);
//                     accumulators[0] += this->intersect(src, dst);
//                 }
//             },
//             galois::chunk_size<CHUNK_SIZE>(), galois::steal(), galois::loopname("TC")
//         );
//     }

//     inline unsigned intersect_merge(unsigned src, unsigned dst) { // src->dst [0->2]
//         unsigned count = 0;
//         auto host_masters = graph.masterNodesRange();
//         for(auto e1 : graph.edges(src)){ // src->alt [0->1]
//             GNode to = graph.getEdgeDst(e1); // 1

//             // Looking for directed edge from min->max
//             GNode min_vertexID = std::min(to, dst); // 1
//             GNode max_vertexID = std::max(to, dst); // 2

//             // If min_vertexID in master, local triangles!
//             if(std::find(host_masters.begin(), host_masters.end(), min_vertexID) != host_masters.end()){
//                 for(auto e : graph.edges(min_vertexID)){
//                     GNode min_vertexID_dst = graph.getEdgeDst(e);
//                     if (min_vertexID_dst == max_vertexID) {
//                         count += 1;
//                         break;
//                     }
//                     // if (max_vertexID > min_vertexID_dst) break; // QUESTION: Guaranteed 2b in order?
//                 }
//             }
//             // Oh no! Dst on a different host!: Gotta request edges: min_vertexID->max_vertexID!
//             else graph.getData(min_vertexID).requested_edges.push_back(max_vertexID);  
//         }
//         return count;
//     }
// };
#include "dist_engine.h"
