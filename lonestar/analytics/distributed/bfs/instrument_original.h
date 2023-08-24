#include <map>

static cll::opt<std::string> graphName("graphName", cll::desc("Name of the input graph"), cll::init("temp"));

template <typename Graph>
struct Instrument {
  Graph* graph;
  uint64_t hostID;
  uint64_t numHosts;

  std::unique_ptr<galois::DGAccumulator<uint64_t>> local_read_stream;
  std::unique_ptr<galois::DGAccumulator<uint64_t>> master_read;
  std::unique_ptr<galois::DGAccumulator<uint64_t>> master_write;
  std::unique_ptr<galois::DGAccumulator<uint64_t>> mirror_read;
  std::unique_ptr<galois::DGAccumulator<uint64_t>> mirror_write;
  std::unique_ptr<galois::DGAccumulator<uint64_t>[]> remote_comm_to_host;
  std::ofstream file;

  void init(uint64_t hid, uint64_t numH, std::unique_ptr<Graph>& g) {
    graph    = g.get();
    hostID   = hid;
    numHosts = numH;

    local_read_stream   = std::make_unique<galois::DGAccumulator<uint64_t>>();
    master_read         = std::make_unique<galois::DGAccumulator<uint64_t>>();
    master_write        = std::make_unique<galois::DGAccumulator<uint64_t>>();
    mirror_read         = std::make_unique<galois::DGAccumulator<uint64_t>>();
    mirror_write        = std::make_unique<galois::DGAccumulator<uint64_t>>();
    remote_comm_to_host = std::make_unique<galois::DGAccumulator<uint64_t>[]>(numH);
    clear();

    // start instrumentation
    file.open(graphName + "_" + std::to_string(numH) + "procs_id" + std::to_string(hid));
    file << "#####   Stat   #####" << std::endl;
    file << "host " << hid << " total edges: " << graph->sizeEdges() << std::endl;
  }

  void clear() {
    local_read_stream->reset();
    master_read->reset();
    master_write->reset();
    mirror_read->reset();
    mirror_write->reset();
    
    for (auto i = 0ul; i < numHosts; i++) {
        remote_comm_to_host[i].reset();
    }
  }

  void record_local_read_stream() {
      *local_read_stream += 1;
  }

  void record_read_random(typename Graph::GraphNode node) {
    auto gid = graph->getGID(node);

    if (graph->isOwned(gid)) { // master
      *master_read += 1;
    }
    else { // mirror
      *mirror_read += 1;
    }
  }

  void record_write_random(typename Graph::GraphNode node, bool comm) {
    auto gid = graph->getGID(node);
    
    if (graph->isOwned(gid)) { // master
      *master_write += 1;
    }
    else { // mirror
      *mirror_write += 1;

      if (comm) {
        remote_comm_to_host[graph->getHostID(gid)] += 1;
      }
    }
  }

  void log_run(uint64_t run) {
    file << "#####   Run " << run << "   #####" << std::endl;
  }

  void log_round(uint64_t num_iterations) {
    auto host_id   = hostID;
    auto num_hosts = numHosts;
    file << "#####   Round " << num_iterations << "   #####" << std::endl;
    file << "host " << host_id << " local read (stream): " << local_read_stream->read_local() << std::endl;
    file << "host " << host_id << " master reads: " << master_read->read_local() << std::endl;
    file << "host " << host_id << " master writes: " << master_write->read_local() << std::endl;
    file << "host " << host_id << " mirror reads: " << mirror_read->read_local() << std::endl;
    file << "host " << host_id << " mirror writes: " << mirror_write->read_local() << std::endl;

    for (uint32_t i = 0; i < num_hosts; i++) {
        file << "host " << host_id << " remote communication for host " << i << ": " << remote_comm_to_host[i].read_local() << std::endl;
    }
  }

};
