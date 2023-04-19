// ****************************************************
//                   dist_engine.h
// ****************************************************
#pragma once
#include "DistBench/Input.h"
#include "DistBench/Start.h"
#include "galois/AtomicHelpers.h"
#include "galois/DReducible.h"
#include "galois/DTerminationDetector.h"
#include "galois/DistGalois.h"
#include "galois/Galois.h"
#include "galois/Version.h"
#include "galois/graphs/GenericPartitioners.h"
#include "galois/graphs/GluonSubstrate.h"
#include "galois/graphs/MiningPartitioner.h"
#include "galois/gstl.h"
#include "galois/runtime/Tracer.h"
#include "galois/substrate/SimpleLock.h"
#include "llvm/Support/CommandLine.h"
#include "pangolin/BfsMining/embedding_list.h"
#include "pangolin/res_man.h"


// ##################################################################
//                      GENERIC HELPER FUNCS
// ##################################################################
template <typename T>
size_t size_of_set_intersection(std::vector<T> vec1, std::vector<T> vec2) {
    size_t count = 0;
    if (vec1.size() <= vec2.size()) {
        std::set<T> vec1_set(vec1.begin(), vec1.end());
        for (auto i : vec2) {
            if (std::find(vec1.begin(), vec1.end(), i) != vec1.end()) count++;
        }
    } else {
        std::set<T> vec2_set(vec2.begin(), vec2.end());
        for (auto i : vec1) {
            if (std::find(vec2.begin(), vec2.end(), i) != vec2.end()) count++;
        }
    }
    return count;
}


// ##################################################################
//                      GRAPH STRUCTS + SEMANTICS
// ##################################################################
struct NodeData {
    std::vector<uint64_t> dests;
    std::vector<uint64_t> requested_edges;
};

// Command Line
namespace cll = llvm::cl;

typedef galois::graphs::MiningGraph<NodeData, void, MiningPolicyDegrees> Graph;
typedef typename Graph::GraphNode GNode;
using DGAccumulatorTy = galois::DGAccumulator<uint64_t>;

// Template Semantics: GraphPtr and SubstratePtr
template <typename NodeData, typename EdgeData>
using MiningGraphPtr = std::unique_ptr<Graph>;

template <typename NodeData, typename EdgeData>
using MiningSubstratePtr = std::unique_ptr<galois::graphs::GluonSubstrate<Graph>>;


// ##################################################################
//                      GLOBAL (aka host) VARS
// ##################################################################
galois::DynamicBitSet bitset_dests;
std::unique_ptr<galois::graphs::GluonSubstrate<Graph>> syncSubstrate;
uint64_t numTriangles; // DGAccumulatorTy numTriangles;
auto& net = galois::runtime::getSystemNetworkInterface();


// ##################################################################
//                      READ + INIT GRAPH
// ##################################################################
template <typename NodeData, typename EdgeData, bool iterateOutEdges = true>
static MiningGraphPtr<NodeData, EdgeData> vertexLoadDGraph(bool loadProxyEdges) {
    galois::StatTimer dGraphTimer("GraphConstructTime", "DistBench");

    dGraphTimer.start();
    const auto& net = galois::runtime::getSystemNetworkInterface();
    MiningGraphPtr<NodeData, EdgeData> loadedGraph = std::make_unique<Graph>(inputFile, net.ID, net.Num, loadProxyEdges, loadProxyEdges);
    assert(loadedGraph != nullptr);
    dGraphTimer.stop();

    return loadedGraph;
}

template <typename NodeData, typename EdgeData, bool iterateOutEdges = true>
std::pair<MiningGraphPtr<NodeData, EdgeData>, MiningSubstratePtr<NodeData, EdgeData>>
vertexBasedDistGraphInitialization(bool loadProxyEdges = true) {
    galois::StatTimer initTimer("DistGraphInitialization", "DistMiningBench");
    std::vector<unsigned> scaleFactor;

    // Create graph + substrate
    MiningGraphPtr<NodeData, EdgeData> g;
    MiningSubstratePtr<NodeData, EdgeData> s;

    // Load graph
    g = vertexLoadDGraph<NodeData, EdgeData, iterateOutEdges>(loadProxyEdges);
    
    // Load substrate
    const auto& net = galois::runtime::getSystemNetworkInterface();
    s = std::make_unique<galois::graphs::GluonSubstrate<Graph>>(*g, net.ID, net.Num, g->isTransposed(), g->cartesianGrid(), partitionAgnostic, commMetadata);
    return std::make_pair(std::move(g), std::move(s));
}


// ##################################################################
//                   PRINT GRAPH (for debugging)
// ##################################################################
std::shared_ptr<galois::substrate::SimpleLock> lock_ptr = std::make_shared<galois::substrate::SimpleLock>();
std::unique_ptr<Graph> print_graph(std::unique_ptr<Graph> hg) {
    galois::do_all(
        galois::iterate(hg->allNodesWithEdgesRange()),  // hg = local subgraph
        [&](const GNode& node) {                        // Plug current guy in here
            std::stringstream sout;
            sout << "\tSrc: " << node << std::endl;
            sout << "\tDsts: " << std::endl;
            for (auto i = hg->edge_begin(node); i != hg->edge_end(node); i++) {
                sout << "\t\t" << hg->getEdgeDst(i) << std::endl;
            }

            // Atomized print!
            lock_ptr->lock();
            printf("Host %d / %d:\n", net.ID, net.Num);
            std::cout << sout.str() << std::endl;
            lock_ptr->unlock();

        },
        galois::steal());
    MPI_Barrier(MPI_COMM_WORLD);
    return hg;
}


// ##################################################################
//             SYNCHRONIZATION [from SyncStructures.h]
// * tell HOW to communicate ... what node data
// ##################################################################
struct SyncPhase1 {
    typedef std::vector<uint64_t> ValTy;

    // Extract Vector! Tell what to communicate:
    // For a node, tell what part of NodeData to Communicate
    static ValTy extract(uint32_t, const struct NodeData& node) {
        return node.dests;  
    }

    // Reduce() -- what do after master receive
    // Reduce will send the data to the host
    // Masters = receivers, so will update master dests
    static bool reduce(uint32_t, struct NodeData& node, ValTy y) {
        for (ValTy::iterator it = y.begin(); it != y.end(); ++it) {
            lock_ptr->lock();
            node.dests.push_back(*it);  // insertBag = per thread-vector:
            lock_ptr->unlock();
        }

        return true;
    }

    static void reset(uint32_t, struct NodeData& node) {
        galois::set(node.dests, (ValTy)0);
    }

    static void setVal(uint32_t, struct NodeData& node, ValTy y) {
        node.dests = y;  // deep copy
    }

    static bool extract_batch(unsigned, uint8_t*, size_t*, DataCommMode*) {return false;}
    static bool extract_batch(unsigned, uint8_t*) { return false; }
    static bool extract_reset_batch(unsigned, uint8_t*, size_t*, DataCommMode*) {return false;}
    static bool extract_reset_batch(unsigned, uint8_t*) { return false; }
    static bool reset_batch(size_t, size_t) { return false; }
    static bool reduce_batch(unsigned, uint8_t*, DataCommMode) {return false;}
    static bool reduce_mirror_batch(unsigned, uint8_t*, DataCommMode) {return false;}
    static bool setVal_batch(unsigned, uint8_t*, DataCommMode) {return false;}  
};


struct SyncPhase2 : public SyncPhase1 {
    typedef std::vector<uint64_t> ValTy;

    static ValTy extract(uint32_t, const struct NodeData& node) {
        return node.requested_edges;
    }

    static bool reduce(uint32_t, struct NodeData& node, ValTy y) {
        for (ValTy::iterator it = y.begin(); it != y.end(); ++it) {
            lock_ptr->lock();
            node.requested_edges.push_back(*it);
            lock_ptr->unlock();
        }

        return true;
    }

    static void reset(uint32_t, struct NodeData& node) {
        galois::set(node.requested_edges, (ValTy)0);
    }

    static void setVal(uint32_t, struct NodeData& node, ValTy y) {
        node.requested_edges = y;  // deep copy
    }
};


// ##################################################################
//                BITSETs [also from SyncStructures.h]
// * tell which nodes to communicate
// ##################################################################
struct BitsetPhase1 {
    static constexpr bool is_vector_bitset() { return false; }
    static constexpr bool is_valid() { return true; }
    static galois::DynamicBitSet& get() { return bitset_dests; }
    static void reset_range(size_t begin, size_t end) {
        bitset_dests.reset(begin, end);
    }
};



// ##################################################################
//         TRIANGLE COUNTING HELPERS [from Pangolin Vertex Miner]
// ##################################################################
template <typename NodeData, typename EdgeData>
std::pair<MiningGraphPtr<NodeData, EdgeData>, size_t> intersect_merge(MiningGraphPtr<NodeData, EdgeData> hg, unsigned src, unsigned dst) {  // src->dst [0->2]
    printf("\n******** intersect_merge (src = %d, dst = %d) ********\n", src, dst);
    size_t count = 0;
    auto host_masters = hg->masterNodesRange();
    for (auto e1 : hg->edges(src)) {    // src->alt [0->1]
        // Looking for directed edge from min->max
        GNode bigger_dstID = hg->getEdgeDst(e1);  // 1
        printf("bigger_dstID = %d\n", bigger_dstID);

        

        // If min_vertexID in master, local triangles!
        if (std::find(host_masters.begin(), host_masters.end(), dst) != host_masters.end()) {
            printf("Dst = %d EDGES:\n", dst);
            for (auto e : hg->edges(dst)) {
                GNode dest_dstID = hg->getEdgeDst(e);
                printf("min_vertexID_dst = %d\n", dest_dstID);
                if (dest_dstID == bigger_dstID) {
                    count += 1;
                    break;
                }
                // ENSURE you dont double-count triangles!
                if (bigger_dstID > dest_dstID) break; 
            }
        }
        // Oh no! Dst on a different host!: Gotta request edges: min_vertexID->max_vertexID!
        else{
            auto min_vertexID = std::min(bigger_dstID, dst);
            auto max_vertexID = std::max(bigger_dstID, dst);
            hg->getData(min_vertexID).requested_edges.push_back(max_vertexID);
        } 
    }

    return std::make_pair(std::move(hg), count);
}

template <typename NodeData, typename EdgeData>
std::unique_ptr<Graph> local_tc_vertex_solver(std::unique_ptr<Graph> hg) {  // vertex parallel
    galois::do_all(
        galois::iterate(hg->masterNodesRange()),
        [&](const GNode& src) {
            for (auto e : hg->edges(src)) {
                auto dst = hg->getEdgeDst(e);
                size_t local_num_triangles = 0;
                std::tie(hg, local_num_triangles) = intersect_merge<NodeData, EdgeData>(std::move(hg), src, dst);
                numTriangles += local_num_triangles;
                printf("Local Num Triangles = %ld", local_num_triangles);
            }
        },
        galois::chunk_size<CHUNK_SIZE>(), galois::steal(), galois::loopname("TC"));
    return hg;
}



// ##################################################################
//                          main()
// ##################################################################
int main(int argc, char** argv) {
    // Initialize Galois Runtime
    galois::DistMemSys G;
    DistBenchStart(argc, argv, name, desc, url);
    auto myrank = net.ID;
    
    // Initialize Graph
    std::unique_ptr<Graph> hg;
    std::tie(hg, syncSubstrate) = vertexBasedDistGraphInitialization<NodeData, void>(false);  // QUESTION
    bitset_dests.resize(hg->size());

    hg = print_graph(std::move(hg));

    // ****************************************************
    //               1. FIND + SEND PAIRS
    // ****************************************************
    /*
    1. Find Pairs ... Where would get stored on host?
        * Iter in jaccard ... send the edge (bc its iden by the src, dest)
        * Adjust the graph struct to accept tuple
    2. for each pair, id the src .. QUESTION: How get src + dest ids? -- Do i use extract?
        * Gotta convert local ID to global ID in BFSS Analytics
        * L2G: libcusp DistributedGraph.h
    */
    printf("Step 1\n");
    std::vector<std::pair<uint64_t, uint64_t>> edge_list;

    printf("\t Creating msgs\n");
    galois::do_all(
        galois::iterate(hg->allNodesWithEdgesRange()),  // hg = local subgraph
        [&](const GNode& node) {                        // Plug current guy in here
            printf("* Node = %d\n", node);
            
            for (auto e : hg->edges(node)){
                auto edge_dest = hg->getEdgeDst(e);
                printf("\t* Dst = %d\n", edge_dest);
                auto src_globalID = hg->getGID(node);  // node.L2G();
                auto dest_globalID = hg->getGID(edge_dest);
                printf("\t* GlobalEdge = %ld -> %ld\n", src_globalID, dest_globalID);

                // Assume: Bi-Directional Edges: So if this is false, dw the dest will handle it!
                // Use globalID for trianlge stuff
                if (src_globalID < dest_globalID) hg->getData(node).dests.push_back(dest_globalID);
            }
        },
        galois::steal());

    // if i am a mirror node, then i do the bitset.set(key=localID) (similar to insert[key])
    // Bitset: Identify WHICH nodes gotta send messages!: Mirror Nodes
    printf("\t Creating bitset\n");
    auto mirrorNodes = hg->getMirrorNodes()[myrank];
    galois::do_all(
        galois::iterate(mirrorNodes.begin(), mirrorNodes.end()),
        [&](const GNode& node) {
            // GOTTA BE LOCAL ID, since per host communication
            bitset_dests.set(node);
        },
        galois::steal());

    // Sync -- Broadcase from mirrors -> masters
    printf("\t PHASE 1: Sync\n");
    syncSubstrate->sync<writeDestination, readSource, SyncPhase1, BitsetPhase1, false>("TC");
    bitset_dests.reset();

    printf("Step 2\n");
    // ****************************************************
    //         2. CALC MISSING EDGES + ASK
    // ****************************************************
    /* // Iterate over master Nodes:
    - Look thru master dests
    - Master Sender: Fire and forget for each pair
    - Master Receiver: Counts or doesnt
    */
    // Counts triangles on each host and place requested edges on each mirror
    hg = local_tc_vertex_solver<NodeData, void>(std::move(hg));

    // Send requests mirrors -> masters
    galois::do_all(
        galois::iterate(mirrorNodes.begin(), mirrorNodes.end()),
        [&](const GNode& node) {
            // GOTTA BE LOCAL ID, since per host communication
            bitset_dests.set(node);
        },
        galois::steal());

    // Actual sending/braodcasting
    syncSubstrate->sync<writeDestination, readSource, SyncPhase2, BitsetPhase1, false>("TC");  // QUESTION

    printf("Step 3\n");
    // ****************************************************
    // 3. Each Master check requests + accumulate numTriangles
    // ****************************************************
    galois::do_all(
        galois::iterate(hg->masterNodesRange()),
        [&](const GNode& node) {
            size_of_set_intersection(hg->getData(node).dests, hg->getData(node).requested_edges);
            numTriangles += size_of_set_intersection(hg->getData(node).dests, hg->getData(node).requested_edges);
        },
        galois::chunk_size<CHUNK_SIZE>(), galois::steal(), galois::loopname("TC"));


    printf("Step 4\n");
    // ****************************************************
    //    4. GLOBAL REDUCE of accumulator (num triangles)
    // ****************************************************
    uint64_t total_triangles = numTriangles;//numTriangles.reduce();
    if (galois::runtime::getSystemNetworkInterface().ID == 0) {
        galois::gPrint("Total number of triangles ", total_triangles, "\n");
    }

    return 0;
}
