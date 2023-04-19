// // ****************************************************
// //                   dist_engine.h
// // ****************************************************
// #pragma once
// #include "DistBench/Input.h"
// #include "DistBench/Start.h"
// #include "galois/AtomicHelpers.h"
// #include "galois/DReducible.h"
// #include "galois/DTerminationDetector.h"
// #include "galois/DistGalois.h"
// #include "galois/Galois.h"
// #include "galois/Version.h"
// #include "galois/graphs/GenericPartitioners.h"
// #include "galois/graphs/GluonSubstrate.h"
// #include "galois/graphs/MiningPartitioner.h"
// #include "galois/gstl.h"
// #include "galois/runtime/Tracer.h"
// #include "galois/substrate/SimpleLock.h"
// #include "llvm/Support/CommandLine.h"
// #include "pangolin/BfsMining/embedding_list.h"
// #include "pangolin/res_man.h"

// #define DEBUG 1


// // ##################################################################
// //                      GENERIC HELPER FUNCS
// // ##################################################################
// template <typename T>
// size_t size_of_set_intersection(std::set<T> org_set, std::vector<T> vec) {
//     size_t count = 0;
//     for (auto i : vec) {
//         if (std::find(org_set.begin(), org_set.end(), i) != org_set.end()) count++;
//     }
//     return count;
// }


// template<typename T>
// void remove_duplicates(std::vector<T>& vec){
//     std::sort(vec.begin(), vec.end()); 
//     auto exclusive_end = std::unique(vec.begin(), vec.end());
//     vec.erase(exclusive_end, vec.end());
// }




// // ##################################################################
// //                      GRAPH STRUCTS + SEMANTICS
// // ##################################################################
// struct NodeData {
//     std::set<uint64_t> dests;
//     std::vector<uint64_t> requested_edges;
// };

// // Command Line
// namespace cll = llvm::cl;

// typedef galois::graphs::DistGraph<NodeData, void> Graph; //typedef galois::graphs::MiningGraph<NodeData, void, MiningPolicyDegrees> Graph;
// typedef typename Graph::GraphNode PGNode;
// using DGAccumulatorTy = galois::DGAccumulator<uint64_t>;

// // Template Semantics: GraphPtr and SubstratePtr
// template <typename NodeData, typename EdgeData>
// using GraphPtr = std::unique_ptr<Graph>;

// template <typename NodeData, typename EdgeData>
// using SubstratePtr = std::unique_ptr<galois::graphs::GluonSubstrate<Graph>>;


// // ##################################################################
// //                      GLOBAL (aka host) VARS
// // ##################################################################
// galois::DynamicBitSet bitset_dests;
// galois::DynamicBitSet bitset_requested_edges;
// std::unique_ptr<galois::graphs::GluonSubstrate<Graph>> syncSubstrate; 
// auto& net = galois::runtime::getSystemNetworkInterface();


// // ##################################################################
// //                      READ + INIT GRAPH
// // ##################################################################
// // template <typename NodeData, typename EdgeData, bool iterateOutEdges = true>
// // static GraphPtr<NodeData, EdgeData> vertexLoadDGraph(bool loadProxyEdges) {
// //     galois::StatTimer dGraphTimer("GraphConstructTime", "DistBench");

// //     dGraphTimer.start();
// //     const auto& net = galois::runtime::getSystemNetworkInterface();
// //     GraphPtr<NodeData, EdgeData> loadedGraph = std::make_unique<Graph>(inputFile, net.ID, net.Num);//std::make_unique<Graph>(inputFile, net.ID, net.Num, loadProxyEdges, loadProxyEdges);
// //     assert(loadedGraph != nullptr);
// //     dGraphTimer.stop();

// //     return loadedGraph;
// // }

// // template <typename NodeData, typename EdgeData, bool iterateOutEdges = true>
// // std::pair<GraphPtr<NodeData, EdgeData>, SubstratePtr<NodeData, EdgeData>>
// // vertexBasedDistGraphInitialization(bool loadProxyEdges = true) {
// //     galois::StatTimer initTimer("DistGraphInitialization", "DistMiningBench");
// //     std::vector<unsigned> scaleFactor;

// //     // Create graph + substrate
// //     GraphPtr<NodeData, EdgeData> g;
// //     SubstratePtr<NodeData, EdgeData> s;

// //     // Load graph
// //     g = vertexLoadDGraph<NodeData, EdgeData, iterateOutEdges>(loadProxyEdges);
    
// //     // Load substrate
// //     const auto& net = galois::runtime::getSystemNetworkInterface();
// //     s = std::make_unique<galois::graphs::GluonSubstrate<Graph>>(*g, net.ID, net.Num, g->isTransposed(), g->cartesianGrid(), partitionAgnostic, commMetadata);
// //     return std::make_pair(std::move(g), std::move(s));
// // }


// // ##################################################################
// //                   PRINT GRAPH (for debugging)
// // ##################################################################
// std::shared_ptr<galois::substrate::SimpleLock> lock_ptr = std::make_shared<galois::substrate::SimpleLock>();
// std::unique_ptr<Graph> print_graph(std::unique_ptr<Graph> hg) {
//     galois::do_all(
//         galois::iterate(hg->allNodesWithEdgesRange()),  // hg = local subgraph
//         [&](const PGNode& node) {                        // Plug current guy in here
//             std::stringstream sout;
//             sout << "\tSrc: " << node << std::endl;
//             sout << "\tDsts: " << std::endl;
//             for (auto i = hg->edge_begin(node); i != hg->edge_end(node); i++) {
//                 sout << "\t\t" << hg->getEdgeDst(i) << std::endl;
//             }

//             // Atomized print!
//             lock_ptr->lock();
//             printf("Host %d / %d:\n", net.ID, net.Num);
//             std::cout << sout.str() << std::endl;
//             lock_ptr->unlock();

//         },
//         galois::steal());
//     MPI_Barrier(MPI_COMM_WORLD);
//     return hg;
// }


// // ##################################################################
// //             SYNCHRONIZATION [from SyncStructures.h]
// // * tell HOW to communicate ... what node data
// // ##################################################################
// struct SyncPhase1 {
//     typedef std::vector<uint64_t> ValTy;

//     // Extract Vector! Tell what to communicate:
//     // For a node, tell what part of NodeData to Communicate
//     static ValTy extract(uint32_t, const struct NodeData& node) {
//         ValTy vec_ver(node.dests.begin(), node.dests.end());
//         return vec_ver;
//         // return node.dests;  
//     }

//     // Reduce() -- what do after master receive
//     // Reduce will send the data to the host
//     // Masters = receivers, so will update master dests
//     static bool reduce(uint32_t, struct NodeData& node, ValTy y) {
//         for (ValTy::iterator it = y.begin(); it != y.end(); ++it) {
//             lock_ptr->lock();
//             node.dests.insert(*it);  // insertBag = per thread-vector:
//             lock_ptr->unlock();
//         }
//         return true;
//     }

//     static void reset(uint32_t, struct NodeData& node) {
//         std::set<uint64_t> empty_set;
//         galois::set(node.dests, empty_set);
//     }

//     static void setVal(uint32_t, struct NodeData& node, ValTy y) {
//         std::set<uint64_t> set_ver(y.begin(), y.end());
//         node.dests = set_ver;
//         // node.dests = y;  // deep copy
//     }

//     static bool extract_batch(unsigned, uint8_t*, size_t*, DataCommMode*) {return false;}
//     static bool extract_batch(unsigned, uint8_t*) { return false; }
//     static bool extract_reset_batch(unsigned, uint8_t*, size_t*, DataCommMode*) {return false;}
//     static bool extract_reset_batch(unsigned, uint8_t*) { return false; }
//     static bool reset_batch(size_t, size_t) { return false; }
//     static bool reduce_batch(unsigned, uint8_t*, DataCommMode) {return false;}
//     static bool reduce_mirror_batch(unsigned, uint8_t*, DataCommMode) {return false;}
//     static bool setVal_batch(unsigned, uint8_t*, DataCommMode) {return false;}  
// };


// struct SyncPhase2 : public SyncPhase1 {
//     typedef std::vector<uint64_t> ValTy;

//     static ValTy extract(uint32_t, const struct NodeData& node) {
//         return node.requested_edges;
//     }

//     static bool reduce(uint32_t, struct NodeData& node, ValTy y) {
//         for (ValTy::iterator it = y.begin(); it != y.end(); ++it) {
//             lock_ptr->lock();
//             node.requested_edges.push_back(*it);
//             lock_ptr->unlock();
//         }
//         return true;
//     }

    
//     static void reset(uint32_t, struct NodeData& node) {
//         galois::set(node.requested_edges, (ValTy)0);
//     }

//     static void setVal(uint32_t, struct NodeData& node, ValTy y) {
//         node.requested_edges = y;  // deep copy
//     }
// };


// // ##################################################################
// //                BITSETs [also from SyncStructures.h]
// // * tell which nodes to communicate
// // ##################################################################
// struct BitsetPhase1 {
//     static constexpr bool is_vector_bitset() { return false; }
//     static constexpr bool is_valid() { return true; }
//     static galois::DynamicBitSet& get() { return bitset_dests; }
//     static void reset_range(size_t begin, size_t end) {
//         bitset_dests.reset(begin, end);
//     }
// };

// struct BitsetPhase2 {
//     static constexpr bool is_vector_bitset() { return false; }
//     static constexpr bool is_valid() { return true; }
//     static galois::DynamicBitSet& get() { return bitset_requested_edges; }
//     static void reset_range(size_t begin, size_t end) {
//         bitset_requested_edges.reset(begin, end);
//     }
// };



// // ##################################################################
// //         TRIANGLE COUNTING HELPERS [from Pangolin Vertex Miner]
// // Should be using dests
// // ##################################################################
// template <typename NodeData, typename EdgeData>
// std::pair<GraphPtr<NodeData, EdgeData>, size_t> intersect_merge(GraphPtr<NodeData, EdgeData> hg, unsigned src, unsigned global_dst) {  // src->dst [0->2]
//     if(DEBUG) printf("\n******** intersect_merge (src = %d, global_dst = %d) ********\n", src, global_dst);
//     size_t count = 0;
//     auto host_masters = hg->masterNodesRange();
//     PGNode local_dstID = hg->getLID(global_dst);
    
//     if (std::find(host_masters.begin(), host_masters.end(), local_dstID) != host_masters.end()){
//         // Local Triangle!
//         for (auto global_e1_dst : hg->getData(src).dests){
//             if(DEBUG) printf("Have: %ld->%d, %ld->%ld ... Missing e = %d -> %ld\n", hg->getGID(src), global_dst, hg->getGID(src), global_e1_dst, global_dst, global_e1_dst);
//             for (auto dst_dst_globalID : hg->getData(local_dstID).dests) {
//                 if(DEBUG) printf("\t dst_dst_globalID = %ld\n", dst_dst_globalID);
//                 if (dst_dst_globalID == global_e1_dst) {
//                     if(DEBUG) printf("*** Local Triangle! %ld - %d - %ld\n", hg->getGID(src), global_dst, dst_dst_globalID);
//                     count += 1;
//                     break;
//                 }
//                 // ENSURE you dont double-count triangles!
//                 if (global_e1_dst > dst_dst_globalID){
//                     if(DEBUG) printf("*** BREAK! %ld > %ld\n", global_e1_dst, dst_dst_globalID);
//                     break;
//                 } 
//             }
//         }
//     }
//     else{
//         // Request missing Edge!
//         for (auto global_e1_dst : hg->getData(src).dests){
//             if(DEBUG) printf("Have: %ld->%d, %ld->%ld ... Missing e = %d -> %ld\n", hg->getGID(src), global_dst, hg->getGID(src), global_e1_dst, global_dst, global_e1_dst);
//             if (global_dst < global_e1_dst) hg->getData(local_dstID).requested_edges.push_back(global_e1_dst);
//         }

//     }

//     return std::make_pair(std::move(hg), count);
// }





// // ##################################################################
// //                          main()
// // ##################################################################
// int main(int argc, char** argv) {
//     // Initialize Galois Runtime
//     galois::DistMemSys G;
//     DistBenchStart(argc, argv, name, desc, url);
//     DGAccumulatorTy numTriangles;
//     numTriangles.reset();
//     auto myrank = net.ID;
    
//     // Initialize Graph
//     std::unique_ptr<Graph> hg;
//     std::tie(hg, syncSubstrate) = distGraphInitialization<NodeData, void>();  // QUESTION
//     bitset_dests.resize(hg->size());
//     bitset_requested_edges.resize(hg->size());

//     if(DEBUG) hg = print_graph(std::move(hg));

//     // ****************************************************
//     //               1. FIND + SEND PAIRS
//     // ****************************************************
//     /*
//     1. Find Pairs ... Where would get stored on host?
//         * Iter in jaccard ... send the edge (bc its iden by the src, dest)
//         * Adjust the graph struct to accept tuple
//     2. for each pair, id the src .. QUESTION: How get src + dest ids? -- Do i use extract?
//         * Gotta convert local ID to global ID in BFSS Analytics
//         * L2G: libcusp DistributedGraph.h
//     */
//     if(DEBUG) printf("Step 1 -- Creating msgs\n");
//     galois::do_all(
//         galois::iterate(hg->allNodesWithEdgesRange()),  // hg = local subgraph
//         [&](const PGNode& node) {                        // Plug current guy in here
//             if(DEBUG) printf("* Node = %d\n", node);
            
//             for (auto e : hg->edges(node)){
//                 auto edge_dest = hg->getEdgeDst(e);
//                 if(DEBUG) printf("\t* Dst = %ld\n", edge_dest);
//                 auto src_globalID = hg->getGID(node);  // node.L2G();
//                 auto dest_globalID = hg->getGID(edge_dest);
//                 if(DEBUG) printf("\t* GlobalEdge = %ld -> %ld\n", src_globalID, dest_globalID);

//                 // Assume: Bi-Directional Edges: So if this is false, dw the dest will handle it!
//                 // Use globalID for trianlge stuff
//                 auto min_globalID = std::min(src_globalID, dest_globalID);
//                 auto max_globalID = std::max(src_globalID, dest_globalID);
//                 hg->getData(hg->getLID(min_globalID)).dests.insert(max_globalID);
//                 if(DEBUG) if (src_globalID < dest_globalID) hg->getData(node).dests.insert(dest_globalID);
//             }
//         },
//         galois::steal());

//     // if i am a mirror node, then i do the bitset.set(key=localID) (similar to insert[key])
//     // Bitset: Identify WHICH nodes gotta send messages!: Mirror Nodes
//     if(DEBUG) printf("***** Creating bitset to send mirror nodes:\n");
//     auto mirrorNodes = hg->getMirrorNodes()[myrank];
//     galois::do_all(
//         galois::iterate(mirrorNodes.begin(), mirrorNodes.end()),
//         [&](const PGNode& node) {
//             // GOTTA BE LOCAL ID, since per host communication
//             bitset_dests.set(node);
//         },
//         galois::steal());

//     // Sync -- Broadcase from mirrors -> masters
//     if(DEBUG) printf("\n **************** PHASE 1: B4 ****************\n");

//     if(DEBUG) galois::do_all(
//         galois::iterate(hg->allNodesWithEdgesRange()),  // hg = local subgraph
//         [&](const PGNode& node) {                        // Plug current guy in here
//             printf("NODE: %d -- \n", node);
//             for (auto dst : hg->getData(node).dests) {
//                 printf("\t dst = %ld\n", dst);
//             }
//         });

    
//     syncSubstrate->sync<writeDestination, readSource, SyncPhase1, BitsetPhase1, false>("TC");
//     bitset_dests.reset();

//     if(DEBUG) printf("\n **************** PHASE 1: AFTER ****************\n");

//     if(DEBUG) galois::do_all(
//         galois::iterate(hg->allNodesWithEdgesRange()),  // hg = local subgraph
//         [&](const PGNode& node) {                        // Plug current guy in here
//             printf("NODE: %d -- \n", node);
//             for (auto dst : hg->getData(node).dests) {
//                 printf("\t dst = %ld\n", dst);
//             }
//         });





//     if(DEBUG) printf("Step 2\n");
//     // ****************************************************
//     //         2. CALC MISSING EDGES + ASK
//     // ****************************************************
//     /* // Iterate over master Nodes:
//     - Look thru master dests
//     - Master Sender: Fire and forget for each pair
//     - Master Receiver: Counts or doesnt
//     */
//     // Counts triangles on each host and place requested edges on each mirror
//     galois::do_all(
//         galois::iterate(hg->masterNodesRange()),
//         [&](const PGNode& src) {
//             for (auto global_dst : hg->getData(src).dests) {
//                 size_t local_num_triangles = 0;
//                 std::tie(hg, local_num_triangles) = intersect_merge<NodeData, void>(std::move(hg), src, global_dst);
//                 numTriangles += local_num_triangles;
//                 if(DEBUG) printf("Local Num Triangles = %ld", local_num_triangles);
//             }
//         },
//         galois::chunk_size<CHUNK_SIZE>(), galois::steal(), galois::loopname("TC"));

//     // Send requests mirrors -> masters
//     galois::do_all(
//         galois::iterate(mirrorNodes.begin(), mirrorNodes.end()),
//         [&](const PGNode& node) {
//             // GOTTA BE LOCAL ID, since per host communication
//             bitset_requested_edges.set(node);
//         },
//         galois::steal());


//     if(DEBUG) printf("\n **************** REQUESTS: BEFORE ****************\n");
//      if(DEBUG) galois::do_all(
//         galois::iterate(hg->allNodesWithEdgesRange()),  // hg = local subgraph
//         [&](const PGNode& node) {                        // Plug current guy in here
//             printf("NODE: %d -- \n", node);
//             for (auto dst : hg->getData(node).requested_edges) {
//                 printf("\t dst = %ld\n", dst);
//             }
//         });


//     // Actual sending/braodcasting
//     syncSubstrate->sync<writeDestination, readSource, SyncPhase2, BitsetPhase2, false>("TC");  // QUESTION

//     // Sync -- Broadcase from mirrors -> masters
//     if(DEBUG) printf("\n **************** REQUESTS: AFTER ****************\n");

   
//     if(DEBUG) galois::do_all(
//         galois::iterate(hg->allNodesWithEdgesRange()),  // hg = local subgraph
//         [&](const PGNode& node) {                        // Plug current guy in here
//             printf("NODE: %d -- \n", node);
//             for (auto dst : hg->getData(node).requested_edges) {
//                 printf("\t dst = %ld\n", dst);
//             }
//         });




//     if(DEBUG) printf("Step 3\n");
//     // ****************************************************
//     // 3. Each Master check requests + accumulate numTriangles
//     // ****************************************************
//     galois::do_all(
//         galois::iterate(hg->masterNodesRange()),
//         [&](const PGNode& node) {
//             numTriangles += size_of_set_intersection(hg->getData(node).dests, hg->getData(node).requested_edges);
//         },
//         galois::chunk_size<CHUNK_SIZE>(), galois::steal(), galois::loopname("TC"));


//     if(DEBUG) printf("\n************** Step 4 ****************\n");
//     // ****************************************************
//     //    4. GLOBAL REDUCE of accumulator (num triangles)
//     // ****************************************************
//     uint64_t total_triangles = numTriangles.reduce();
//     if (galois::runtime::getSystemNetworkInterface().ID == 0) {
//         galois::gPrint("Total number of triangles ", total_triangles, "\n");
//     }

//     return 0;
// }
