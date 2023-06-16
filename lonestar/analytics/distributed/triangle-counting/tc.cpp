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

/* This is an implementation of Distributed multi-GPU triangle counting code.
 * The single GPU code which is executed on GPU is generated using the IrGL
 * compiler. Currently, it does not support distributed multi-CPU code.
 *
 * TODO implement CPU kernel
 */

#include "DistBench/MiningStart.h"
#include "galois/DistGalois.h"
#include "galois/DReducible.h"
#include "galois/DTerminationDetector.h"
#include "galois/gstl.h"
#include "galois/graphs/GenericPartitioners.h"
#include "galois/graphs/MiningPartitioner.h"
#include "galois/runtime/Tracer.h"

#include <iostream>
#include <limits>

#ifdef GALOIS_ENABLE_GPU
#include "tc_cuda.h"
struct CUDA_Context* cuda_ctx;
#else
enum { CPU, GPU_CUDA };
int personality = CPU;
#endif

namespace cll = llvm::cl;
static cll::opt<std::string> graphName("graphName", cll::desc("Name of the input graph"), cll::init("temp"));

constexpr static const char* const REGION_NAME = "TC";

/*******************************************************************************
 * Graph structure declarations + other initialization
 ******************************************************************************/

typedef galois::graphs::MiningGraph<void, void, MiningPolicyDegrees> Graph;
typedef typename Graph::GraphNode GNode;

std::unique_ptr<galois::graphs::GluonEdgeSubstrate<Graph>> syncSubstrate;

template <bool async>
struct TC {
  Graph* graph;
  using DGAccumulatorTy = galois::DGAccumulator<uint64_t>;
  DGAccumulatorTy& numTriangles;
    
  galois::DGAccumulator<uint64_t>& local_read_stream_phase1;
  galois::DGAccumulator<uint64_t>& local_read_stream_phase2;
  galois::DGAccumulator<uint64_t>& local_read_stream_phase3;
  galois::DGAccumulator<uint64_t>& local_read_random_phase1;
  galois::DGAccumulator<uint64_t>& local_read_random_phase2;
  galois::DGAccumulator<uint64_t>& local_write_random_phase1;
  galois::DGAccumulator<uint64_t>& remote_read_random_phase1;
  galois::DGAccumulator<uint64_t>& remote_write_stream_phase2;
  galois::DGAccumulator<uint64_t>& remote_write_random_phase1;
  galois::DGAccumulator<uint64_t>* remote_read_random_to_host;
  galois::DGAccumulator<uint64_t>* remote_write_stream_to_host;
  galois::DGAccumulator<uint64_t>* remote_write_random_to_host;
  galois::DGAccumulator<uint64_t>* local_read_stream_to_host;

  std::ofstream& file;

  TC(Graph* _graph, 
     DGAccumulatorTy& _numTriangles, 
     galois::DGAccumulator<uint64_t>& _local_read_stream_phase1,
     galois::DGAccumulator<uint64_t>& _local_read_stream_phase2,
     galois::DGAccumulator<uint64_t>& _local_read_stream_phase3,
     galois::DGAccumulator<uint64_t>& _local_read_random_phase1,
     galois::DGAccumulator<uint64_t>& _local_read_random_phase2,
     galois::DGAccumulator<uint64_t>& _local_write_random_phase1,
     galois::DGAccumulator<uint64_t>& _remote_read_random_phase1,
     galois::DGAccumulator<uint64_t>& _remote_write_stream_phase2,
     galois::DGAccumulator<uint64_t>& _remote_write_random_phase1,
     galois::DGAccumulator<uint64_t>* _remote_read_random_to_host,
     galois::DGAccumulator<uint64_t>* _remote_write_stream_to_host,
     galois::DGAccumulator<uint64_t>* _remote_write_random_to_host,
     galois::DGAccumulator<uint64_t>* _local_read_stream_to_host,
     std::ofstream& _file)
      : graph(_graph), 
        numTriangles(_numTriangles), 
        local_read_stream_phase1(_local_read_stream_phase1),
        local_read_stream_phase2(_local_read_stream_phase2),
        local_read_stream_phase3(_local_read_stream_phase3),
        local_read_random_phase1(_local_read_random_phase1),
        local_read_random_phase2(_local_read_random_phase2),
        local_write_random_phase1(_local_write_random_phase1),
        remote_read_random_phase1(_remote_read_random_phase1),
        remote_write_stream_phase2(_remote_write_stream_phase2),
        remote_write_random_phase1(_remote_write_random_phase1),
        remote_read_random_to_host(_remote_read_random_to_host),
        remote_write_stream_to_host(_remote_write_stream_to_host),
        remote_write_random_to_host(_remote_write_random_to_host),
        local_read_stream_to_host(_local_read_stream_to_host),
        file(_file) {}

  // use the below line once CPU code is added
  void static go(Graph& _graph, 
                 galois::DGAccumulator<uint64_t>& local_read_stream_phase1,
                 galois::DGAccumulator<uint64_t>& local_read_stream_phase2,
                 galois::DGAccumulator<uint64_t>& local_read_stream_phase3,
                 galois::DGAccumulator<uint64_t>& local_read_random_phase1,
                 galois::DGAccumulator<uint64_t>& local_read_random_phase2,
                 galois::DGAccumulator<uint64_t>& local_write_random_phase1,
                 galois::DGAccumulator<uint64_t>& remote_read_random_phase1,
                 galois::DGAccumulator<uint64_t>& remote_write_stream_phase2,
                 galois::DGAccumulator<uint64_t>& remote_write_random_phase1,
                 galois::DGAccumulator<uint64_t>* remote_read_random_to_host,
                 galois::DGAccumulator<uint64_t>* remote_write_stream_to_host,
                 galois::DGAccumulator<uint64_t>* remote_write_random_to_host,
                 galois::DGAccumulator<uint64_t>* local_read_stream_to_host,
                 std::ofstream& file) {
    unsigned _num_iterations = 0;
    DGAccumulatorTy numTriangles;
    syncSubstrate->set_num_round(_num_iterations);
    numTriangles.reset();
    const auto& allMasterNodes = _graph.masterNodesRange();

    uint32_t num_hosts = _graph.getNumHosts();
    uint64_t host_id = galois::runtime::getSystemNetworkInterface().ID;
	
    local_read_stream_phase1.reset();
    local_read_stream_phase2.reset();
    local_read_stream_phase3.reset();
    local_read_random_phase1.reset();
    local_read_random_phase2.reset();
    local_write_random_phase1.reset();
    remote_read_random_phase1.reset();
    remote_write_stream_phase2.reset();
    remote_write_random_phase1.reset();

    for (uint32_t i=0; i<num_hosts; i++) {
      remote_read_random_to_host[i].reset();
      remote_write_stream_to_host[i].reset();
      remote_write_random_to_host[i].reset();
      local_read_stream_to_host[i].reset();
    }

#ifdef GALOIS_ENABLE_GPU
    if (personality == GPU_CUDA) { ///< GPU TC.
      std::string impl_str(syncSubstrate->get_run_identifier("TC"));
      galois::StatTimer StatTimer_cuda(impl_str.c_str(), REGION_NAME);
      StatTimer_cuda.start();
      uint64_t num_local_triangles = 0;
      TC_masterNodes_cuda(num_local_triangles, cuda_ctx);
      numTriangles += num_local_triangles;
      StatTimer_cuda.stop();
    } else { ///< CPU TC.
#endif
      galois::do_all(
          galois::iterate(allMasterNodes), 
          TC(&_graph, 
             numTriangles,
             local_read_stream_phase1,
             local_read_stream_phase2,
             local_read_stream_phase3,
             local_read_random_phase1,
             local_read_random_phase2,
             local_write_random_phase1,
             remote_read_random_phase1,
             remote_write_stream_phase2,
             remote_write_random_phase1,
             remote_read_random_to_host,
             remote_write_stream_to_host,
             remote_write_random_to_host,
             local_read_stream_to_host,
             file),
          galois::steal(),
          galois::loopname(syncSubstrate->get_run_identifier("TC").c_str()));
#ifdef GALOIS_ENABLE_GPU
    }
#endif
    
    file << "#####   Round 0   #####" << std::endl;
	file << "host " << host_id << " round local read stream: " << local_read_stream_phase1.read_local() << std::endl;
	file << "host " << host_id << " round local read random: " << local_read_random_phase1.read_local() << std::endl;
	file << "host " << host_id << " round local write random: " << local_write_random_phase1.read_local() << std::endl;
	file << "host " << host_id << " round remote read stream: 0" << std::endl;
	file << "host " << host_id << " round remote read random: " << remote_read_random_phase1.read_local() << std::endl;
	file << "host " << host_id << " round remote write stream: 0" << std::endl;
	file << "host " << host_id << " round remote write random: " << remote_write_random_phase1.read_local() << std::endl;

	for (uint32_t i=0; i<num_hosts; i++) {
		file << "host " << host_id << " remote read stream to host " << i << ": 0" << std::endl;
		file << "host " << host_id << " remote read random to host " << i << ": " << remote_read_random_to_host[i].read_local() << std::endl;
		file << "host " << host_id << " remote write stream to host " << i << ": 0" << std::endl;
		file << "host " << host_id << " remote write random to host " << i << ": " << remote_write_random_to_host[i].read_local() << std::endl;
	}
    
    file << "#####   Round 1   #####" << std::endl;
	file << "host " << host_id << " round local read stream: " << local_read_stream_phase2.read_local() << std::endl;
	file << "host " << host_id << " round local read random: " << local_read_random_phase2.read_local() << std::endl;
	file << "host " << host_id << " round local write random: 0" << std::endl;
	file << "host " << host_id << " round remote read stream: 0" << std::endl;
	file << "host " << host_id << " round remote read random: 0" << std::endl;
	file << "host " << host_id << " round remote write stream: " << remote_write_stream_phase2.read_local() << std::endl;
	file << "host " << host_id << " round remote write random: 0" << std::endl;

	for (uint32_t i=0; i<num_hosts; i++) {
		file << "host " << host_id << " remote read stream to host " << i << ": 0" << std::endl;
		file << "host " << host_id << " remote read random to host " << i << ": 0" << std::endl;
		file << "host " << host_id << " remote write stream to host " << i << ": " << remote_write_stream_to_host[i].read_local() << std::endl;
		file << "host " << host_id << " remote write random to host " << i << ": 0" << std::endl;
	}

    for (uint32_t i=0; i<num_hosts; i++) {
      remote_write_stream_to_host[i].reduce();
      local_read_stream_to_host[i].reduce();
    }
    
    file << "#####   Round 2   #####" << std::endl;
	file << "host " << host_id << " round local read stream: " << local_read_stream_phase3.read_local() + remote_write_stream_to_host[host_id].read_local() + local_read_stream_to_host[host_id].read_local() << std::endl;
	file << "host " << host_id << " round local read random: 0" << std::endl;
	file << "host " << host_id << " round local write random: 0" << std::endl;
	file << "host " << host_id << " round remote read stream: 0" << std::endl;
	file << "host " << host_id << " round remote read random: 0" << std::endl;
	file << "host " << host_id << " round remote write stream: 0" << std::endl;
	file << "host " << host_id << " round remote write random: 0" << std::endl;

	for (uint32_t i=0; i<num_hosts; i++) {
		file << "host " << host_id << " remote read stream to host " << i << ": 0" << std::endl;
		file << "host " << host_id << " remote read random to host " << i << ": 0" << std::endl;
		file << "host " << host_id << " remote write stream to host " << i << ": 0" << std::endl;
		file << "host " << host_id << " remote write random to host " << i << ": 0" << std::endl;
	}

    uint64_t total_triangles = numTriangles.reduce();
    if (galois::runtime::getSystemNetworkInterface().ID == 0) {
      galois::gPrint("Total number of triangles ", total_triangles, "\n");
    }
  }

  void operator()(GNode v) const {
    // uint64_t v_GID = graph->getGID(v);
      
    local_read_stream_phase1 += 1;
    local_read_stream_phase2 += 1;
    local_read_stream_phase3 += 1;

    size_t numTriangles_local = 0;
      
    for (auto vIter : graph->edges(v)) {
      auto w = graph->getEdgeDst(vIter);
      uint64_t w_GID = graph->getGID(w);
      bool w_owned = graph->isOwned(w_GID);
      unsigned to_host = graph->getHostID(w_GID);

      if (w_owned) {
          local_read_stream_phase1 += 2;
          local_read_random_phase1 += 2;

          local_write_random_phase1 = 1;
      }
      else {
          local_read_stream_phase1 += 1;
          remote_read_random_phase1 += 1;
          remote_read_random_to_host[to_host] += 1;
          
          remote_write_random_phase1 += 1;
          remote_write_random_to_host[to_host] += 1;
      }

      local_read_stream_phase2 += 1;

      for (auto vIter2 : graph->edges(v)) {
        auto x = graph->getEdgeDst(vIter2);      
        uint64_t x_GID = graph->getGID(x);
        // bool x_owned = graph->isOwned(x_GID);

        local_read_stream_phase2 += 1;

        if (w_GID < x_GID) {
            if (w_owned) {
                local_read_random_phase2 += 1;
            }
            else {
                remote_write_stream_phase2 += 1;        
                remote_write_stream_to_host[to_host] += 1;
            }

            for (auto wIter: graph->edges(w)) {
                auto m = graph->getEdgeDst(wIter);             
                // uint64_t m_GID = graph->getGID(m);
                // bool m_owned = graph->isOwned(m_GID);

                if (w_owned) {
                    local_read_stream_phase2 += 1;
                }
                else {
                    local_read_stream_to_host[to_host] += 1;
                }
                
                if (m == x) {
                    ++numTriangles_local;
                }
            }
        }
      }
    } ///< Finding triangles is done.
    numTriangles += numTriangles_local;
  } ///< CPU operator is done.
};

/*******************************************************************************
 * Main
 ******************************************************************************/

constexpr static const char* const name =
    "TC - Distributed Multi-GPU Triangle Counting ";
constexpr static const char* const desc = "TC on Distributed GPU (D-IrGL).";
constexpr static const char* const url  = nullptr;

int main(int argc, char** argv) {
  galois::DistMemSys G;
  DistBenchStart(argc, argv, name, desc, url);

  if (!symmetricGraph) {
    GALOIS_DIE("This application requires a symmetric graph input;"
               " please use the -symmetricGraph flag "
               " to indicate the input is a symmetric graph.");
  }

  const auto& net = galois::runtime::getSystemNetworkInterface();

  galois::StatTimer StatTimer_total("TimerTotal", REGION_NAME);

  StatTimer_total.start();
  std::unique_ptr<Graph> hg;
#ifdef GALOIS_ENABLE_GPU
  std::tie(hg, syncSubstrate) =
      distGraphInitialization<void, void>(&cuda_ctx, false);
#else
  std::tie(hg, syncSubstrate) = distGraphInitialization<void, void>(false);
#endif

  if (personality == GPU_CUDA) {
#ifdef GALOIS_ENABLE_GPU
    std::string timer_str("SortEdgesGPU");
    galois::StatTimer edgeSortTime("SortEdgesGPU", REGION_NAME);
    edgeSortTime.start();
    sortEdgesByDestination_cuda(cuda_ctx);
    edgeSortTime.stop();
#else
    abort();
#endif
  } else if (personality == CPU) {
    galois::StatTimer edgeSortTime("SortEdgesCPU", REGION_NAME);
    edgeSortTime.start();
    hg->sortEdgesByDestination();
    edgeSortTime.stop();
  }

  uint32_t num_hosts = hg->getNumHosts();
  uint64_t host_id = net.ID;

  ///! accumulators for use in operators
  galois::DGAccumulator<uint64_t> DGAccumulator_numTriangles;
  galois::DGAccumulator<uint64_t> local_read_stream_phase1;
  galois::DGAccumulator<uint64_t> local_read_stream_phase2;
  galois::DGAccumulator<uint64_t> local_read_stream_phase3;
  galois::DGAccumulator<uint64_t> local_read_random_phase1;
  galois::DGAccumulator<uint64_t> local_read_random_phase2;
  galois::DGAccumulator<uint64_t> local_write_random_phase1;
  galois::DGAccumulator<uint64_t> remote_read_random_phase1;
  galois::DGAccumulator<uint64_t> remote_write_stream_phase2;
  galois::DGAccumulator<uint64_t> remote_write_random_phase1;
  galois::DGAccumulator<uint64_t> remote_read_random_to_host[num_hosts];
  galois::DGAccumulator<uint64_t> remote_write_stream_to_host[num_hosts];
  galois::DGAccumulator<uint64_t> remote_write_random_to_host[num_hosts];
  galois::DGAccumulator<uint64_t> local_read_stream_to_host[num_hosts];

  std::ofstream file;
  file.open(graphName + "_" + std::to_string(num_hosts) + "procs_id" + std::to_string(host_id));
  file << "#####   Stat   #####" << std::endl;
  file << "host " << host_id << " total edges: " << hg->sizeEdges() << std::endl;

  for (auto run = 0; run < numRuns; ++run) {
    galois::gPrint("[", net.ID, "] TC::go run ", run, " called\n");
    std::string timer_str("Timer_" + std::to_string(run));
    galois::StatTimer StatTimer_main(timer_str.c_str(), REGION_NAME);

    StatTimer_main.start();
    TC<false>::go(*hg,
                  local_read_stream_phase1,
                  local_read_stream_phase2,
                  local_read_stream_phase3,
                  local_read_random_phase1,
                  local_read_random_phase2,
                  local_write_random_phase1,
                  remote_read_random_phase1,
                  remote_write_stream_phase2,
                  remote_write_random_phase1,
                  remote_read_random_to_host,
                  remote_write_stream_to_host,
                  remote_write_random_to_host,
                  local_read_stream_to_host,
                  file);
    StatTimer_main.stop();

    syncSubstrate->set_num_run(run + 1);
  }
  StatTimer_total.stop();

  if (output) {
    galois::gError("output requested but this application doesn't support it");
    return 1;
  }

  return 0;
}
