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

#define DEBUG 0
#define BENCH 1

// ##################################################################
//                      GENERIC HELPER FUNCS
// ##################################################################
template <typename T>
size_t size_of_set_intersection(std::set<T> org_set, std::vector<T> vec) {
    size_t count = 0;
    for (auto i : vec) {
        if (std::find(org_set.begin(), org_set.end(), i) != org_set.end()) count++;
    }
    return count;
}

// ##################################################################
//                      GRAPH STRUCTS + SEMANTICS
// ##################################################################
struct NodeData {
    std::set<uint64_t> dests;
    std::vector<uint64_t> requested_edges;
    std::shared_ptr<galois::substrate::SimpleLock> my_lock = std::make_shared<galois::substrate::SimpleLock>();
    std::atomic<uint64_t> num_intersections = 0;
};

// Command Line
namespace cll = llvm::cl;

static cll::opt<std::string> graphName("graphName", cll::desc("Name of the input graph"), cll::init("temp"));

typedef galois::graphs::DistGraph<NodeData, void> Graph;
typedef typename Graph::GraphNode PGNode;
using DGAccumulatorTy = galois::DGAccumulator<uint64_t>;
using GAccumulatorTy = galois::GAccumulator<uint64_t>;
DGAccumulatorTy* global_num_tri_ptr;
GAccumulatorTy* thread_num_tri_ptr;

// Template Semantics: GraphPtr and SubstratePtr
template <typename NodeData, typename EdgeData>
using GraphPtr = std::unique_ptr<Graph>;

template <typename NodeData, typename EdgeData>
using SubstratePtr = std::unique_ptr<galois::graphs::GluonSubstrate<Graph>>;

// ##################################################################
//                      GLOBAL (aka host) VARS
// ##################################################################
galois::DynamicBitSet bitset_dests;
galois::DynamicBitSet bitset_requested_edges;
std::unique_ptr<galois::graphs::GluonSubstrate<Graph>> syncSubstrate;
auto& net = galois::runtime::getSystemNetworkInterface();

std::shared_ptr<galois::substrate::SimpleLock> lock_ptr = std::make_shared<galois::substrate::SimpleLock>();
std::shared_ptr<galois::substrate::SimpleLock> print_lock_ptr = std::make_shared<galois::substrate::SimpleLock>();

// ##################################################################
//                   PRINT GRAPH (for debugging)
// ##################################################################
void print_graph(Graph& hg_ref) {
    print_lock_ptr->lock();
    printf("Host %d / %d:\n", net.ID, net.Num);
    galois::do_all(
        galois::iterate(hg_ref.allNodesWithEdgesRange()),  // hg = local subgraph
        [&](const PGNode& node) {                       // Plug current guy in here
            std::stringstream sout;
            sout << "\tSrc: " << node << "; GLOBAL = " << hg_ref.getGID(node) << std::endl;
            sout << "\tDsts: " << std::endl;
            for (auto i = hg_ref.edge_begin(node); i != hg_ref.edge_end(node); i++) {
                auto ldst = hg_ref.getEdgeDst(i);
                sout << "\t\tLocal = " << ldst << ", global = " << hg_ref.getGID(ldst) << std::endl;
            }
            std::cout << sout.str() << std::endl;
        },
        galois::steal());
    print_lock_ptr->unlock();
    galois::runtime::getHostBarrier().wait(); 
}

// ##################################################################
//             SYNCHRONIZATION [from SyncStructures.h]
// * tell HOW to communicate ... what node data
// ##################################################################
struct SyncPhase1 {
    typedef std::vector<uint64_t> ValTy;

    // Extract Vector of dests! Tell what to communicate:
    // For a node, tell what part of NodeData to Communicate
    static ValTy extract(uint32_t, const struct NodeData& node) {
        ValTy vec_ver(node.dests.begin(), node.dests.end());
        return vec_ver;
    }

    // Reduce() -- what do after master receive
    // Reduce will send the data to the host
    // Masters = receivers, so will update master dests
    static bool reduce(uint32_t, struct NodeData& node, ValTy y) {
        for (ValTy::iterator it = y.begin(); it != y.end(); ++it) {
            node.my_lock->lock();
            node.dests.insert(*it);  // insertBag = per thread-vector:
            node.my_lock->unlock();
        }
        return true;
    }

    static void reset(uint32_t, struct NodeData& node) {
        std::set<uint64_t> empty_set;
        galois::set(node.dests, empty_set);
    }

    static void setVal(uint32_t, struct NodeData& node, ValTy y) {
        std::set<uint64_t> set_ver(y.begin(), y.end());
        node.dests = set_ver; // deep copy
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

struct SyncPhase2 {
    typedef std::vector<uint64_t> ValTy;

    static ValTy extract(uint32_t, const struct NodeData& node) {
        return node.requested_edges;
    }

    // Count number of requests found in the dests
    static bool reduce(uint32_t, struct NodeData& node, ValTy y) {
        for (ValTy::iterator it = y.begin(); it != y.end(); ++it) {
            // if (std::find(node.dests.begin(), node.dests.end(), *it) != node.dests.end()) *thread_num_tri_ptr += 1;
			node.my_lock->lock();
			node.requested_edges.push_back(*it);
			node.my_lock->unlock();
        }
        return true;
    }

    static void reset(uint32_t, struct NodeData& node) {
        galois::set(node.requested_edges, (ValTy)0);
    }

    static void setVal(uint32_t, struct NodeData& node, ValTy y) {
        node.requested_edges = y;  // deep copy
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
struct BitsetPhase1 {
    static constexpr bool is_vector_bitset() { return false; }
    static constexpr bool is_valid() { return true; }
    static galois::DynamicBitSet& get() { return bitset_dests; }
    static void reset_range(size_t begin, size_t end) {
        bitset_dests.reset(begin, end);
    }
};

struct BitsetPhase2 {
    static constexpr bool is_vector_bitset() { return false; }
    static constexpr bool is_valid() { return true; }
    static galois::DynamicBitSet& get() { return bitset_requested_edges; }
    static void reset_range(size_t begin, size_t end) {
        bitset_requested_edges.reset(begin, end);
    }
};

// ##################################################################
//         TRIANGLE COUNTING HELPERS [from Pangolin Vertex Miner]
// Should be using dests
// ##################################################################
template <typename NodeData, typename EdgeData>
size_t intersect_merge(Graph& hg, 
						unsigned src, 
						unsigned global_dst,
						galois::DGAccumulator<uint64_t>& local_read_stream,
						galois::DGAccumulator<uint64_t>& local_read_random,
						galois::DGAccumulator<uint64_t>& remote_write_stream,
						galois::DGAccumulator<uint64_t>* remote_write_to_host) {
	size_t count = 0;
    PGNode local_dstID = hg.getLID(global_dst);

    if (hg.isOwned(global_dst)) {
        // Local Triangle!
        for (auto global_e1_dst : hg.getData(src).dests) {
			local_read_stream += 1;

            if (global_dst < global_e1_dst) {  // Ignore self-loops (ie 0-1, 0-1, missing 1-1)
                // Find Missing Edge from global_dst -> global_e1_dst on this host
                auto& dstData = hg.getData(local_dstID);
				local_read_random += 1;
				
				for (auto dst_dst_globalID : dstData.dests) {
					local_read_stream += 1;

                    if (dst_dst_globalID == global_e1_dst) {
                        count += 1;
                        break;
                    }
                }
            }
        }
    } else {
        // Request missing Edge!
        for (auto global_e1_dst : hg.getData(src).dests) {
			local_read_stream += 1;
			
            // Find missing edge on other hosts
            if (global_dst < global_e1_dst) {
                auto& myLocalData = hg.getData(local_dstID);
                myLocalData.my_lock->lock();
                myLocalData.requested_edges.push_back(global_e1_dst);
                myLocalData.my_lock->unlock();

				remote_write_stream += 1;
				
				unsigned to_host = hg.getHostID(global_dst);
				remote_write_to_host[to_host] += 1;

                // Already know local_dstID is a mirror node since we are in the else block
                bitset_requested_edges.set(local_dstID);
            }
        }
    }

    return count;
}

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

    double ph1_gen_msg_start = 0;
    double ph1_gen_msg_end = 0;
    double ph1_comm_start = 0;
    double ph1_comm_end = 0;
    double ph2_gen_msg_start = 0;
    double ph2_gen_msg_end = 0;
    double ph2_comm_start = 0;
    double ph2_comm_end = 0;
	double ph3_comp_start = 0;
	double ph3_comp_end = 0;
    if(BENCH) e2e_start = MPI_Wtime();
    
    // Initialize Galois Runtime
    galois::DistMemSys G;
    DistBenchStart(argc, argv, name, desc, url);
    global_num_tri_ptr = new DGAccumulatorTy();
    thread_num_tri_ptr = new GAccumulatorTy();
    global_num_tri_ptr->reset();
    thread_num_tri_ptr->reset();

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
    bitset_dests.resize(hg_ref.size());
    bitset_requested_edges.resize(hg_ref.size());

	uint32_t num_hosts = hg->getNumHosts();
	uint64_t host_id = galois::runtime::getSystemNetworkInterface().ID;

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

    // ALGORITHM TIME
    if(BENCH){
        galois::runtime::getHostBarrier().wait();
        algo_start = MPI_Wtime();
        ph1_gen_msg_start = MPI_Wtime();
    }

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

    // ****************************************************
    //               1. FIND + SEND PAIRS
    // ****************************************************
    // Iterate over all host nodes and store directed edges from smallerGlobalID -> largerGlobalID
    galois::do_all(
        galois::iterate(hg_ref.allNodesWithEdgesRange()),
        [&](const PGNode& node) {
            auto src_globalID = hg_ref.getGID(node);
			local_read_stream += 1;

            // Iterate over all edges in a node
            for (auto e : hg_ref.edges(node)) {
				local_read_stream += 1;

                auto edge_dest = hg_ref.getEdgeDst(e);
                auto dest_globalID = hg_ref.getGID(edge_dest);
				
				if (hg_ref.isOwned(dest_globalID)) {
					local_read_random += 1;
				}
				else {
					remote_read_random += 1;
					
					unsigned to_host = hg_ref.getHostID(dest_globalID);
					remote_read_to_host[to_host] += 1;
				}

                // Assume: Bi-Directional Edges: So if this is false, dw the dest will handle it!
                // Based on the current edge, add smallerGlobalID -> largerGlobalID to smallerGlobalID's NodeData
                auto min_globalID = std::min(src_globalID, dest_globalID);
                auto max_globalID = std::max(src_globalID, dest_globalID);
				auto min_node = hg_ref.getLID(min_globalID);

                // Ensure you dont have self-loops
                if(min_globalID < max_globalID){
                    auto& minData = hg_ref.getData(hg_ref.getLID(min_globalID));

                    // fine-grain locking for the nodes
                    minData.my_lock->lock();
                    minData.dests.insert(max_globalID);  // ok since its a set .. no duplicates
                    minData.my_lock->unlock();
					
					if (hg_ref.isOwned(min_globalID)) {
						local_write_random += 1;
					}
					else {
						remote_write_random += 1;

						unsigned to_host = hg_ref.getHostID(min_globalID);
						remote_write_to_host[to_host] += 1;

                    	// Add Mirrors to the bitset, bc want only mirrors to communicate to masters
                    	bitset_dests.set(min_node);  // GOTTA BE LOCAL ID, since per host communication

					}
                } 
            }
        },
        galois::steal());

	file << "#####   Round 0   #####" << std::endl;
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

    if(BENCH){
        galois::runtime::getHostBarrier().wait();
        ph1_gen_msg_end = MPI_Wtime();
        std::cout << "Time_Phase1_GenMsg, " << ph1_gen_msg_end - ph1_gen_msg_start << "\n";
        ph1_comm_start = MPI_Wtime();
    }

    // PHASE 1 BROADCAST [Works!]: Send dests from mirrors to masters!
    syncSubstrate->sync<writeDestination, readSource, SyncPhase1, BitsetPhase1, false>("Phase1Communication");
    bitset_dests.reset();

    if(BENCH){
        ph1_comm_end = MPI_Wtime();
        std::cout << "Time_Phase1_Comm, " << ph1_comm_end - ph1_comm_start << "\n";
        ph2_gen_msg_start = MPI_Wtime();
    }
	
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

    // ****************************************************
    //         2. CALC MISSING EDGES + ASK
    // ****************************************************
    // Counts triangles on each host and place requested edges on each mirror
    galois::do_all(
        galois::iterate(hg_ref.masterNodesRange()),
        [&](const PGNode& src) {
			auto& srcData = hg_ref.getData(src);
			local_read_stream += 1;

			for (auto global_dst : srcData.dests) {
				local_read_stream += 1;

				// assert(hg_ref.getGID(src) < global_dst);
                *thread_num_tri_ptr += intersect_merge<NodeData, void>(hg_ref, src, global_dst, local_read_stream, local_read_random, remote_write_stream, remote_write_to_host);
            }
        },
        galois::chunk_size<CHUNK_SIZE>(), galois::steal(), galois::loopname("LocalTCCount"));
	
	file << "#####   Round 1   #####" << std::endl;
	file << "host " << host_id << " round local read stream: " << local_read_stream.read_local() << std::endl;
	file << "host " << host_id << " round local read random: " << local_read_random.read_local() << std::endl;
	file << "host " << host_id << " round local write random: " << local_write_random.read_local() << std::endl;
	file << "host " << host_id << " round remote read stream: " << remote_read_stream.read_local() << std::endl;
	file << "host " << host_id << " round remote read random: " << remote_read_random.read_local() << std::endl;
	file << "host " << host_id << " round remote write stream: " << remote_write_stream.read_local() << std::endl;
	file << "host " << host_id << " round remote write random: " << remote_write_random.read_local() << std::endl;

	for (uint32_t i=0; i<num_hosts; i++) {
		file << "host " << host_id << " remote read stream to host " << i << ": " << remote_read_to_host[i].read_local() << std::endl;
		file << "host " << host_id << " remote read random to host " << i << ": 0" << std::endl;
		file << "host " << host_id << " remote write stream to host " << i << ": " << remote_write_to_host[i].read_local() << std::endl;
		file << "host " << host_id << " remote write random to host " << i << ": 0" << std::endl;
	}

    if(BENCH){
        galois::runtime::getHostBarrier().wait();
        ph2_gen_msg_end = MPI_Wtime();
        std::cout << "Time_Phase2_GenMsg, " << ph2_gen_msg_end - ph2_gen_msg_start << "\n";
        ph2_comm_start = MPI_Wtime();
    }

    // SEND REQUESTS from mirrors to masters!
    syncSubstrate->sync<writeDestination, readSource, SyncPhase2, BitsetPhase2, false>("Phase2Communication");
    bitset_requested_edges.reset();

    if(BENCH){
        ph2_comm_end = MPI_Wtime();
        std::cout << "Time_Phase2_Comm, " << ph2_comm_end - ph2_comm_start << "\n";
		ph3_comp_start = MPI_Wtime();
    }
	
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

    // ****************************************************
    //         3. HANDLE REQUESTS
    // ****************************************************
    // Handle requested edges on each master
    galois::do_all(
        galois::iterate(hg_ref.masterNodesRange()),
        [&](const PGNode& src) {
			auto& srcData = hg_ref.getData(src);
			local_read_stream += 1;

			for (auto global_requested_dst : srcData.requested_edges) {
				local_read_stream += 1;

				for (auto global_dst : srcData.dests) {
					local_read_stream += 1;

					if (global_requested_dst == global_dst) {
						*thread_num_tri_ptr += 1;
						break;
					}
				}
            }
        },
        galois::chunk_size<CHUNK_SIZE>(), galois::steal(), galois::loopname("LocalTCCount"));
	
	file << "#####   Round 2   #####" << std::endl;
	file << "host " << host_id << " round local read stream: " << local_read_stream.read_local() << std::endl;
	file << "host " << host_id << " round local read random: " << local_read_random.read_local() << std::endl;
	file << "host " << host_id << " round local write random: " << local_write_random.read_local() << std::endl;
	file << "host " << host_id << " round remote read stream: " << remote_read_stream.read_local() << std::endl;
	file << "host " << host_id << " round remote read random: " << remote_read_random.read_local() << std::endl;
	file << "host " << host_id << " round remote write stream: " << remote_write_stream.read_local() << std::endl;
	file << "host " << host_id << " round remote write random: " << remote_write_random.read_local() << std::endl;

	for (uint32_t i=0; i<num_hosts; i++) {
		file << "host " << host_id << " remote read stream to host " << i << ": 0" << std::endl;
		file << "host " << host_id << " remote read random to host " << i << ": 0" << std::endl;
		file << "host " << host_id << " remote write stream to host " << i << ": 0" << std::endl;
		file << "host " << host_id << " remote write random to host " << i << ": 0" << std::endl;
	}

    if(BENCH){
        galois::runtime::getHostBarrier().wait();
        ph3_comp_end = MPI_Wtime();
        std::cout << "Time_Phase3_Comp, " << ph3_comp_end - ph3_comp_start << "\n";
    }

    // ****************************************************
    //    4. GLOBAL REDUCE of accumulator (num triangles)
    // ****************************************************
    // Ensure all triangles/requests been counted
    // galois::runtime::getHostBarrier().wait(); // necessary? No bc sync adds to count and sync is non-blocking 

    // Count All the local triangle counts!
    *global_num_tri_ptr += thread_num_tri_ptr->reduce();
    uint64_t total_triangles = global_num_tri_ptr->reduce();
    if (net.ID == 0) galois::gPrint("Total number of triangles ", total_triangles, "\n");

    if(BENCH){
        galois::runtime::getHostBarrier().wait();
        algo_end = MPI_Wtime();
        e2e_end = MPI_Wtime();
        std::cout << "Time_TC_Algo, " << algo_end - algo_start << "\n";
        std::cout << "Time_E2E, " << e2e_end - e2e_start << "\n";
    }

    delete(global_num_tri_ptr);
    delete(thread_num_tri_ptr);


    return 0;
}
