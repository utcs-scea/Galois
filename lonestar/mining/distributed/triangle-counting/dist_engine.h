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

#include <iostream>

#define BENCH 1
#define SORTED_EDGES 1

static cll::opt<uint64_t> query_size("q",
                            cll::desc("Maximum queries per round (default 1)"),
                            cll::init(1));

// ##################################################################
//                      GRAPH STRUCTS + SEMANTICS
// ##################################################################
struct NodeData {
    std::set<uint64_t> queries;
    boost::iterators::counting_iterator<uint64_t> edge_iterator;
};

// Command Line
namespace cll = llvm::cl;

static cll::opt<std::string> graphName("graphName", cll::desc("Name of the input graph"), cll::init("temp"));

typedef galois::graphs::DistGraph<NodeData, void> Graph;
typedef typename Graph::GraphNode PGNode;
using DGAccumulatorTy = galois::DGAccumulator<uint64_t>;

// ##################################################################
//                      GLOBAL (aka host) VARS
// ##################################################################
galois::DynamicBitSet bitset_queries;
std::unique_ptr<galois::graphs::GluonSubstrate<Graph>> syncSubstrate;

// ##################################################################
//             SYNCHRONIZATION [from SyncStructures.h]
// * tell HOW to communicate ... what node data
// ##################################################################
struct SyncQueries {
    typedef std::vector<uint64_t> ValTy;

    // Extract Vector of dests! Tell what to communicate:
    // For a node, tell what part of NodeData to Communicate
    static ValTy extract(uint32_t, const struct NodeData& node) {
        ValTy vec_ver(node.queries.begin(), node.queries.end());
        return vec_ver;
    }

    // Reduce() -- what do after master receive
    // Reduce will send the data to the host
    // Masters = receivers, so will update master dests
    static bool reduce(uint32_t, struct NodeData& node, ValTy y) {
        for (ValTy::iterator it = y.begin(); it != y.end(); ++it) {
            node.queries.insert(*it);  // insertBag = per thread-vector:
        }
        return true;
    }

    static void reset(uint32_t, struct NodeData& node) {
        node.queries.clear();
    }

    static void setVal(uint32_t, struct NodeData& node, ValTy y) {
        std::set<uint64_t> set_ver(y.begin(), y.end());
        node.queries = set_ver; // deep copy
    }

    static bool extract_batch(unsigned, uint8_t*, size_t*, DataCommMode*) { return false; }
    static bool extract_batch(unsigned, uint8_t*) { return false; }
    static bool extract_reset_batch(unsigned, uint8_t*, size_t*, DataCommMode*) { return false; }
    static bool extract_reset_batch(unsigned, uint8_t*) { return false; }
    static bool reset_batch(size_t, size_t) { return false; }
    static bool reduce_batch(unsigned, uint8_t*, DataCommMode) { return false; }
    static bool reduce_mirror_batch(unsigned, uint8_t*, DataCommMode) { return false; }
    static bool setVal_batch(unsigned, uint8_t*, DataCommMode) { return false; }
};

// ##################################################################
//                BITSETs [also from SyncStructures.h]
// * tell which nodes to communicate
// ##################################################################
GALOIS_SYNC_STRUCTURE_BITSET(queries);

// ##################################################################
//                          main()
// ##################################################################
int main(int argc, char** argv) {
    // galois::StatTimer e2e_timer("TimeE2E", "TC");
    // galois::StatTimer algo_timer("TimeAlgo", "TC");
    // e2e_timer.start();
    double e2e_start = 0;
    double e2e_end = 0;
    double mk_graph_start = 0;
    double mk_graph_end = 0;
    double algo_start = 0;
    double algo_end = 0;

    if(BENCH) e2e_start = MPI_Wtime();
    
    // Initialize Galois Runtime
    galois::DistMemSys G;
    DistBenchStart(argc, argv, name, desc, url);
    DGAccumulatorTy num_triangles;
    num_triangles.reset();

    // Initialize Graph
    if(BENCH){
        galois::runtime::getHostBarrier().wait();
        mk_graph_start = MPI_Wtime();
    }
    std::unique_ptr<Graph> hg;
    std::tie(hg, syncSubstrate) = distGraphInitialization<NodeData, void>();
    if(BENCH){
        galois::runtime::getHostBarrier().wait();
        mk_graph_end = MPI_Wtime();
        std::cout << "Time_Graph_Creation, " << mk_graph_end - mk_graph_start << "\n";
    }

    Graph& hg_ref = *hg;
    bitset_queries.resize(hg_ref.size());

    auto& net = galois::runtime::getSystemNetworkInterface();
	uint32_t num_hosts = hg->getNumHosts();
	uint64_t host_id = net.ID;
    if (host_id == 0) {
        galois::gPrint("#Nodes: ", hg->size(), "\n#Edges: ", hg->sizeEdges(), "\n");
    }

	galois::DGAccumulator<uint64_t> local_read_stream;
	galois::DGAccumulator<uint64_t> local_read_random;
	galois::DGAccumulator<uint64_t> local_write_random;
	galois::DGAccumulator<uint64_t> remote_read_stream;
	galois::DGAccumulator<uint64_t> remote_read_random;
	galois::DGAccumulator<uint64_t> remote_write_stream;
	galois::DGAccumulator<uint64_t> remote_write_random;
	galois::DGAccumulator<uint64_t> remote_read_to_host[num_hosts];
	galois::DGAccumulator<uint64_t> remote_write_to_host[num_hosts];

	std::ofstream file;
	file.open(graphName + "_" + std::to_string(num_hosts) + "procs_id" + std::to_string(host_id));
	file << "#####   Stat   #####" << std::endl;
	file << "host " << host_id << " total edges: " << hg->sizeEdges() << std::endl;

    const auto& allNodes = hg_ref.allNodesRange();
    const auto& nodesWithEdges = hg_ref.allNodesWithEdgesRange();

    galois::do_all(
            galois::iterate(allNodes),
            [&](auto node) {
                auto& data = hg_ref.getData(node);
                data.edge_iterator = hg_ref.edge_begin(node);
                // local_read_stream += 1;
                data.queries.clear();
                },
            galois::no_stats(), galois::steal(),
            galois::loopname(syncSubstrate->get_run_identifier("Initialization").c_str()));

    // ALGORITHM TIME
    if(BENCH){
        galois::runtime::getHostBarrier().wait();
        algo_start = MPI_Wtime();
    }

    uint64_t _num_iterations = 0;
    // using DGTerminatorDetector =
    // typename std::conditional<async, galois::DGTerminator<unsigned int>,
    //                         galois::DGAccumulator<unsigned int>>::type;
    galois::DGAccumulator<uint64_t> dga;

    do {
      syncSubstrate->set_num_round(_num_iterations);
      dga.reset();

      local_read_stream.reset();
      local_read_random.reset();
      local_write_random.reset();
      remote_read_stream.reset();
      remote_read_random.reset();
      remote_write_stream.reset();
      remote_write_random.reset();

      for (uint32_t i=0; i<num_hosts; i++) {
        remote_read_to_host[i].reset();
        remote_write_to_host[i].reset();
      }

      galois::do_all(
            galois::iterate(allNodes),
            [&](Graph::GraphNode A) {
                auto& data = hg_ref.getData(A);
                data.queries.clear();
                auto A_gid = hg_ref.getGID(A);
                local_read_stream += 1;
                for (;data.queries.size() < query_size && data.edge_iterator != hg_ref.edge_end(A); data.edge_iterator++) {
                    local_read_stream += 1;
                    auto C = hg_ref.getEdgeDst(data.edge_iterator);
                    auto C_gid = hg_ref.getGID(C);
                    if (hg_ref.isOwned(C_gid)) {
                        local_read_random += 1;
                    }
                    else {
                        remote_read_random += 1;
                        
                        unsigned to_host = hg_ref.getHostID(C_gid);
                        remote_read_to_host[to_host] += 1;
                    }
                    if (A_gid > C_gid) {  // A > C; only query lower neighbors
                        data.queries.insert(C_gid);
                        if (hg_ref.isOwned(A_gid)) {
                            local_write_random += 1;
                        } else {
                            remote_write_random += 1;
                            unsigned to_host = hg_ref.getHostID(A_gid);
                            remote_write_to_host[to_host] += 1;
                        }
                        bitset_queries.set(A);
                        dga += 1;
                    }
#if SORTED_EDGES
                    else {
                        break;
                    }
#endif
                }
            },
            galois::no_stats(), galois::steal(),
            galois::loopname(syncSubstrate->get_run_identifier("Initialization").c_str()));

      syncSubstrate->sync<writeSource, readDestination, SyncQueries,
                          Bitset_queries, /*async*/ false>("TC");

        galois::do_all(
            galois::iterate(nodesWithEdges),
            [&](Graph::GraphNode B) {
                
#if SORTED_EDGES
                auto& data = hg_ref.getData(B);
                auto e = data.edge_iterator;
#else
                auto e = hg_ref.edge_begin(B);
#endif
                auto B_gid = hg_ref.getGID(B);
                local_read_stream += 1;
                for (; e != hg_ref.edge_end(B); e++) {
                    local_read_stream += 1;
                    auto A = hg_ref.getEdgeDst(e);
                    auto A_gid = hg_ref.getGID(A);
                    if (hg_ref.isOwned(A_gid)) {
                        local_read_random += 1;
                    }
                    else {
                        remote_read_random += 1;
                        
                        unsigned to_host = hg_ref.getHostID(A_gid);
                        remote_read_to_host[to_host] += 1;
                    }
                    // A > B; only pull queries from upper neghbors
                    if (B_gid >= A_gid) {
                        continue;
                    }
                    auto& A_data = hg_ref.getData(A);
                    for (auto&& C_gid: A_data.queries) {
                        if (C_gid >= B_gid) { // B > C; only query lower neighbors
#if SORTED_EDGES
                            break;
#else
                            continue;
#endif
                        }
                        for (auto&& e_: hg_ref.edges(B)) {
                            local_read_stream += 1;
                            auto C_ = hg_ref.getEdgeDst(e_);
                            auto C_gid_ = hg_ref.getGID(C_);
                            if (hg_ref.isOwned(C_gid_)) {
                                local_read_random += 1;
                            }
                            else {
                                remote_read_random += 1;
                                
                                unsigned to_host = hg_ref.getHostID(C_gid_);
                                remote_read_to_host[to_host] += 1;
                            }
                            if (C_gid == C_gid_) {
                                num_triangles += 1;
                                continue;
                            }
#if SORTED_EDGES
                            else if (C_gid_ > C_gid) {
                                break;
                            }
#endif
                        }
                    }
                }
            },
            galois::no_stats(), galois::steal(),
            galois::loopname(syncSubstrate->get_run_identifier("MakeQueries").c_str()));

    auto num_works = dga.reduce(syncSubstrate->get_run_identifier());
    if (num_works < hg_ref.size()) query_size = (uint64_t)query_size << 1;
    // galois::gPrint(_num_iterations, " ", num_works, "\n");

    file << "#####   Round " << _num_iterations << "   #####" << std::endl;
	file << "host " << host_id << " round local read stream: " << local_read_stream.read_local() << std::endl;
	file << "host " << host_id << " round local read random: " << local_read_random.read_local() << std::endl;
	file << "host " << host_id << " round local write random: " << local_write_random.read_local() << std::endl;
	file << "host " << host_id << " round remote read stream: " << remote_read_stream.read_local() << std::endl;
	file << "host " << host_id << " round remote read random: " << remote_read_random.read_local() << std::endl;
	file << "host " << host_id << " round remote write stream: " << remote_write_stream.read_local() << std::endl;
	file << "host " << host_id << " round remote write random: " << remote_write_random.read_local() << std::endl;

	for (uint32_t i=0; i<num_hosts; i++) {
		file << "host " << host_id << " remote read stream to host " << i << ": 0" << std::endl;
		file << "host " << host_id << " remote read random to host " << i << ": " << remote_read_to_host[i].read_local() << std::endl;
		file << "host " << host_id << " remote write stream to host " << i << ": 0" << std::endl;
		file << "host " << host_id << " remote write random to host " << i << ": " << remote_write_to_host[i].read_local() << std::endl;
	}

    //   galois::runtime::reportStat_Tsum(
    //       REGION_NAME, syncSubstrate->get_run_identifier("NumWorkItems"),
    //       (unsigned long)dga.read_local());
      ++_num_iterations;
    } while (dga.reduce(syncSubstrate->get_run_identifier()));

    auto total_triangles = num_triangles.reduce();
    if (net.ID == 0) galois::gPrint("Total number of triangles ", total_triangles, "\n");

    if(BENCH){
        galois::runtime::getHostBarrier().wait();
        algo_end = MPI_Wtime();
        e2e_end = MPI_Wtime();
        std::cout << "Time_TC_Algo, " << algo_end - algo_start << "\n";
        std::cout << "Time_E2E, " << e2e_end - e2e_start << "\n";
    }

    return 0;
}
