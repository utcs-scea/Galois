// ****************************************************
//                   dist_engine.h
// ****************************************************
#pragma once
#include "DistBench/MiningStart.h"
#include "galois/DReducible.h"
#include "galois/DTerminationDetector.h"
#include "galois/DistGalois.h"
#include "galois/Galois.h"
#include "galois/graphs/GenericPartitioners.h"
#include "galois/graphs/GluonSubstrate.h"
#include "galois/graphs/MiningPartitioner.h"
#include "galois/gstl.h"
#include "galois/runtime/Tracer.h"
#include "galois/substrate/SimpleLock.h"
#include "pangolin/BfsMining/embedding_list.h"
#include "pangolin/res_man.h"

// ./gen-rand-tricount 1 1 2 2 57 19 19

/*
distGraphInitialization() {
#endif
  using Graph     = galois::graphs::DistGraph<NodeData, EdgeData>;
  using Substrate = galois::graphs::GluonSubstrate<Graph>;
  std::vector<unsigned> scaleFactor;
  DistGraphPtr<NodeData, EdgeData> g;
  DistSubstratePtr<NodeData, EdgeData> s;

#ifdef GALOIS_ENABLE_GPU
  internal::heteroSetup(scaleFactor);
#endif
  g = loadDistGraph<NodeData, EdgeData, iterateOutEdges>(scaleFactor);
  // load substrate
  const auto& net = galois::runtime::getSystemNetworkInterface();
  s = std::make_unique<Substrate>(*g, net.ID, net.Num, g->isTransposed(),
                                  g->cartesianGrid(), partitionAgnostic,
                                  commMetadata);

// marshal graph to GPU as necessary
#ifdef GALOIS_ENABLE_GPU
  marshalGPUGraph(s, cuda_ctx);
#endif

  return std::make_pair(std::move(g), std::move(s));
}

*/

struct NodeData {
    std::vector<uint64_t> dests;
};

galois::DynamicBitSet bitset_dests;

namespace cll = llvm::cl;
typedef galois::graphs::MiningGraph<NodeData, void, MiningPolicyDegrees> Graph;
std::unique_ptr<galois::graphs::GluonEdgeSubstrate<Graph>> syncSubstrate;
typedef typename Graph::GraphNode GNode;

auto& net = galois::runtime::getSystemNetworkInterface();

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

// SyncStructures.h
struct Phase1Sync {
    typedef std::vector<uint64_t> ValTy;

    // for a node: heres what i communicate -- Nodedata
    // Extract Vector! Tell what to comm
    static ValTy extract(uint32_t node_id, const struct NodeData& node) {
        return node.dests;  // should return vec?
    }

    static bool extract_batch(unsigned, uint8_t*, size_t*, DataCommMode*) {
        return false;
    }

    static bool extract_batch(unsigned, uint8_t*) { return false; }

    static bool extract_reset_batch(unsigned, uint8_t*, size_t*,
                                    DataCommMode*) {
        return false;
    }

    static bool extract_reset_batch(unsigned, uint8_t*) { return false; }

    static bool reset_batch(size_t, size_t) { return false; }

    // Reduce() -- what do after master receive: ALSO modeify p_map with msg (embed_list) .. if mirror dont add it lmao
    // Per-Host: p_map 4 global IDs; Per-Vertex: num_msgs (missing edges)
    // Reduce will send the data to the host --> Update pangolin structure
    static bool reduce(uint32_t node_id, struct NodeData& node, ValTy y) {
        /*
        1. receveir msg
        2. Add to host p_map (global var)
        */
        // Add to master

        for (ValTy::iterator it = y.begin(); it != y.end(); ++it) {
            lock_ptr->lock();
            node.dests.push_back(*it);
            lock_ptr->unlock();
        }

        return true;
    }

    static bool reduce_batch(unsigned, uint8_t*, DataCommMode) {
        return false;
    }

    static bool reduce_mirror_batch(unsigned, uint8_t*, DataCommMode) {
        return false;
    }

    static void reset(uint32_t, struct NodeData& node) {
        galois::set(node.dests, (ValTy)0);
    }

    static void setVal(uint32_t, struct NodeData& node, ValTy y) {
        node.dests = y;  // deep copy
    }

    static bool setVal_batch(unsigned, uint8_t*, DataCommMode) {
        return false;
    }
};

struct My_Bitset {
    static constexpr bool is_vector_bitset() { return false; }

    static constexpr bool is_valid() { return true; }

    static galois::DynamicBitSet& get() { return bitset_dests; }

    static void reset_range(size_t begin, size_t end) {
        bitset_dests.reset(begin, end);
    }
};

// filter which nodes u actually communicate
//

// Make a pmap -- EXAMPLE in engine.h

int main(int argc, char** argv) {
    // Initialize Galois Runtime
    galois::DistMemSys G;

    DistBenchStart(argc, argv, name, desc, url);

    auto myrank = net.ID;
    // auto nprocs = net.Num;

    // Initialize Graph
    std::unique_ptr<Graph> hg;
    std::tie(hg, syncSubstrate) = distGraphInitialization<NodeData, void>(); // QUESTION
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
    std::vector<std::pair<uint64_t, uint64_t>> edge_list;

    galois::do_all(
        galois::iterate(hg->allNodesWithEdgesRange()),  // hg = local subgraph
        [&](const GNode& node) {                        // Plug current guy in here
            auto edges_iter = hg->edge_begin(node);
            auto edges_end = hg->edge_end(node);

            while (edges_iter != edges_end) {
                auto edge_dest = hg->getEdgeDst(*edges_iter);
                auto src_globalID = hg->getGID(node);  // node.L2G();
                auto dest_globalID = hg->getGID(edge_dest);

                // Assume: Bi-Directional Edges: So if this is false, dw the dest will handle it!
                // Use globalID for trianlge stuff
                if (src_globalID < dest_globalID) hg->getData(node).dests.push_back(dest_globalID);
            }
        },
        galois::steal());

    // if i am a mirror node, then i do the bitset.set(key=localID) (similar to insert[key])
    // Bitset: Identify WHICH nodes gotta send messages!: Mirror Nodes
    auto mirrorNodes = hg->getMirrorNodes()[myrank];
    galois::do_all(
        galois::iterate(mirrorNodes.begin(), mirrorNodes.end()),
        [&](const GNode& node) {
            // GOTTA BE LOCAL ID, since per host communication
            bitset_dests.set(node);
        },
        galois::steal());

    // Sync -- Broadcast + get/sync Responses ... syncs to the master and communicates res to mirrors
    // Bitset def in SyncStructures.h
    syncSubstrate->sync<writeDestination, readSource, Phase1Sync, My_Bitset, false>("TC"); // QUESTION
    bitset_dests.reset();

    // Print line! -- in reduce, add a print

    // ****************************************************
    //         2. CALC MISSING EDGES + COMM BROADCASE
    // ****************************************************
    // Iterate over master Nodes:
    /*
    - Look thru master dests
    - Master Sender: Fire and forget for each pair
    - Master Receiver: Counts or doesnt
    */

    // ****************************************************
    //    3. GLOBAL REDUCE -- Get final triangle count
    // ****************************************************

    /*
     * The master can request missing edge? --> if get response u inc numTriangles at dest
     * syncSubtrate: WriteAtDest, ReadAtSrc, XX, YY
     * sendEdges <-- QUESTION: Can i j send the edge? will it know where to send it to?
     * SEND DOUBLE SIDED EDGE
     * QUESTION: OR need to create a sendParitialPattern?
     */

    // ****************************************************
    //         3. RUN NORMAL TRIANGLE COUNTING
    // ****************************************************

    // using DGAccumulatorTy = galois::DGAccumulator<uint64_t>;
    // DGAccumulatorTy& numTriangles;

    // Make an edge iterator
    // Cant nest a do_all
    // Each thread pinned to core automatically -- hwtoql?

    // galois use a for loop using hg.begin, hg.end
    // use lock around cout ... + build local string <-- before, gotta convert to global ID

    // use sync substrate

    /*
    *********** PREPROCESSING ***********
    // What is miner.read_graph? miner.initizlize? miner.tc_solver? exactly doing?
        * read_graph:
            - reads .gr into struct
    //     * init: (Vertex_Miner)
    //         - galois::on_each -- Inits counters
    //         - init_emb_list -- this->emb_list.init(this->graph, this->max_size, enable_dag);
    //         - checks pattern type
    // // QUESTION: Need to change read_graph to be DistGraph?
    //     * Make a dist-graph pangolin Graph on new branch
    //     * Ensure normal DistGraph works, then check if proting pangolin works
    // // QUESTION: Difference between LonestarMineStart and DistBenchStart?
    //     * Make own Command Line Options
    //     * Can inport a basic file?

    // //
    // // 3. Send embedding/edgeList to master .. Ref MiningPartitioner.h: 1041
    //


    // *********** For each host, call engine.h stuff and aggregate counts ***********
    // // 4. For each master, count the triangles using the engine.h semantics:: galois forall loop
    // QUESTION: Would smth like this work?

    // unsigned _num_iterations = 0;
    // DGAccumulatorTy numTriangles;
    // syncSubstrate->set_num_round(_num_iterations);
    // numTriangles.reset();

    // galois::do_all(
    //         galois::iterate(allMasterNodes), TC(&_graph, numTriangles),
    //         galois::steal(),
    //         galois::loopname(syncSubstrate->get_run_identifier("TC").c_str()));

    // uint64_t total_triangles = numTriangles.reduce();
    // if (galois::runtime::getSystemNetworkInterface().ID == 0) {
    //   galois::gPrint("Total number of triangles ", total_triangles, "\n");
    // }
    // */

    // // ---------------------- From engine.h ----------------------
    // LonestarMineStart(argc, argv, name, desc, url);

    // galois::StatTimer totalTime("TimerTotal");
    // totalTime.start();

    // AppMiner miner(k, numThreads);
    // galois::StatTimer Tinitial("GraphReadingTime");
    // Tinitial.start();
    // miner.read_graph(filetype, inputFile);
    // Tinitial.stop();
    // ResourceManager rm;
    // for (unsigned nt = 0; nt < num_trials; nt++) {
    //     galois::StatTimer Tinitemb("EmbInitTime");
    //     Tinitemb.start();
    //     miner.initialize(pattern_filename); // doesnt actually use the pattern_filenae
    //     Tinitemb.stop();

    //     galois::StatTimer execTime("Timer_0");
    //     execTime.start();

    //     // Use triangle-specific solver (has DAG Optimization)
    //     #ifdef TRIANGLE
    //     miner.tc_solver();
    //     #else
    //     miner.solver();
    //     #endif

    //     execTime.stop();
    //     miner.print_output();
    //     miner.clean();
    // }
    // std::cout << "\n\t" << rm.get_peak_memory() << "\n\n";

    // totalTime.stop();

    return 0;
}
