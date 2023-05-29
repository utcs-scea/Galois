#include <map>

static cll::opt<std::string>
    logID("logid", cll::desc("Name of the input graph"), cll::init("ego"));

constexpr unsigned CACHE_LEVELS = 10;

template <typename Graph>
struct Instrument {
  Graph* graph;
  uint64_t hostID;
  uint64_t numHosts;

  // size of total nodes; cache level for each mirror (lower level higher
  // in-degree)
  std::vector<unsigned> node_cache_level;

  std::unique_ptr<galois::DGAccumulator<uint64_t>> local_read_stream;
  std::unique_ptr<galois::DGAccumulator<uint64_t>> master_read;
  std::unique_ptr<galois::DGAccumulator<uint64_t>> master_write;
  std::unique_ptr<galois::DGAccumulator<uint64_t>[]> mirror_read;
  std::unique_ptr<galois::DGAccumulator<uint64_t>[]> mirror_write;
  std::unique_ptr<galois::DGAccumulator<uint64_t>[]> remote_read;
  std::unique_ptr<galois::DGAccumulator<uint64_t>[]> remote_write;
  std::unique_ptr<std::unique_ptr<galois::DGAccumulator<uint64_t>[]>[]>
      remote_read_to_host;
  std::unique_ptr<std::unique_ptr<galois::DGAccumulator<uint64_t>[]>[]>
      remote_write_to_host;
  std::unique_ptr<std::unique_ptr<galois::DGAccumulator<uint64_t>[]>[]>
      remote_comm_to_host;
  std::ofstream file;

  void init(uint64_t hid, uint64_t numH, std::unique_ptr<Graph>& g) {
    /**
     * Counts cache level for each mirror node.
     */
    graph    = g.get();
    hostID   = hid;
    numHosts = numH;

    // set up final result container
    node_cache_level.resize(graph->size());

    local_read_stream = std::make_unique<galois::DGAccumulator<uint64_t>>();
    master_read       = std::make_unique<galois::DGAccumulator<uint64_t>>();
    master_write      = std::make_unique<galois::DGAccumulator<uint64_t>>();
    mirror_read =
        std::make_unique<galois::DGAccumulator<uint64_t>[]>(CACHE_LEVELS + 1);
    mirror_write =
        std::make_unique<galois::DGAccumulator<uint64_t>[]>(CACHE_LEVELS + 1);
    remote_read =
        std::make_unique<galois::DGAccumulator<uint64_t>[]>(CACHE_LEVELS + 1);
    remote_write =
        std::make_unique<galois::DGAccumulator<uint64_t>[]>(CACHE_LEVELS + 1);
    remote_read_to_host =
        std::make_unique<std::unique_ptr<galois::DGAccumulator<uint64_t>[]>[]>(
            numH);
    remote_write_to_host =
        std::make_unique<std::unique_ptr<galois::DGAccumulator<uint64_t>[]>[]>(
            numH);
    remote_comm_to_host =
        std::make_unique<std::unique_ptr<galois::DGAccumulator<uint64_t>[]>[]>(
            numH);
    for (uint32_t i = 0; i < numH; i++) {
      remote_read_to_host[i] =
          std::make_unique<galois::DGAccumulator<uint64_t>[]>(CACHE_LEVELS + 1);
      remote_write_to_host[i] =
          std::make_unique<galois::DGAccumulator<uint64_t>[]>(CACHE_LEVELS + 1);
      remote_comm_to_host[i] =
          std::make_unique<galois::DGAccumulator<uint64_t>[]>(CACHE_LEVELS + 1);
    }
    clear();

    // load transposed graph (to count incoming degrees)
    std::vector<unsigned> _;
    auto tgr =
        loadDistGraph<typename Graph::GraphNode, typename Graph::EdgeType,
                      /*iterateOutEdges*/ false>(_);

    // in-degree counting (using transposed graph)
    galois::InsertBag<std::pair<uint64_t, typename Graph::GraphNode>>
        indeg_nodes;
    const auto& allNodes = tgr->allNodesRange();
    galois::do_all(
        galois::iterate(allNodes),
        [&](auto node) {
          // ignore master nodes
          if (tgr->isOwned(tgr->getGID(node))) {
            return;
          }
          indeg_nodes.emplace(std::make_pair(
              std::distance(tgr->edge_begin(node), tgr->edge_end(node)), node));
        },
        galois::steal(), galois::no_stats());

    /**
     * NOTE: alternative way to count incoming degree in original graph
     *
    constexpr uint64_t TEST_NODE = 45355;
    std::atomic<uint64_t> indeg{0};
    const auto& allNodes = graph->allNodesRange();
    galois::do_all(
        galois::iterate(allNodes),
        [&](auto src) {
          for (auto e : graph->edges(src)) {
            auto node = graph->getEdgeDst(e);
            if (graph->isOwned(graph->getGID(node))) {
              continue;
            }
            if (node == TEST_NODE) {
              galois::atomicAdd(indeg, 1ul);
            }
          }
        },
        galois::steal(), galois::no_stats());
    assert(indeg.load(std::memory_order_relaxed) == 27);
     */

    // descending sort
    std::multimap<uint64_t, typename Graph::GraphNode, std::greater<int>>
        sorted_indeg_nodes(indeg_nodes.begin(), indeg_nodes.end());
    // cut into levels
    auto [level_size, surplus] =
        std::div((int64_t)sorted_indeg_nodes.size(), (int64_t)CACHE_LEVELS);
    auto it = sorted_indeg_nodes.begin(), end = sorted_indeg_nodes.end();
    for (unsigned cache_level = 1; cache_level <= CACHE_LEVELS; cache_level++) {
      if (surplus) {
        end = std::next(it, level_size + 1);
        surplus--;
      } else {
        end = std::next(it, level_size);
      }
      galois::do_all(
          galois::iterate(it, end),
          [&](auto& deg_node) {
            auto& [deg, node]      = deg_node;
            node_cache_level[node] = cache_level;
          },
          galois::no_stats());
      it = end;
    }

    // start instrumentation
    file.open(logID + "_" + std::to_string(numH) + "procs_id" +
              std::to_string(hid));
    file << "#####   Stat   #####" << std::endl;
    file << "host " << hid << " total edges: " << graph->sizeEdges()
         << std::endl;
  }

  void clear() {
    local_read_stream->reset();
    master_read->reset();
    master_write->reset();
    for (auto i = 0ul; i < CACHE_LEVELS + 1; i++) {
      mirror_read[i].reset();
      mirror_write[i].reset();
      remote_read[i].reset();
      remote_write[i].reset();
    }
    for (auto i = 0ul; i < numHosts; i++) {
      for (auto j = 0ul; j < CACHE_LEVELS + 1; j++) {
        remote_read_to_host[i][j].reset();
        remote_write_to_host[i][j].reset();
        remote_comm_to_host[i][j].reset();
      }
    }
  }

  void record_local_read_stream() { *local_read_stream += 1; }

  void record_read_random(typename Graph::GraphNode node) {
    auto gid = graph->getGID(node);
    if (graph->isOwned(gid)) { // master
      *master_read += 1;
      return;
    }
    // mirror
    int cache_level = node_cache_level[node];

    for (int i = 0; i < cache_level; i++) { // different configs
      remote_read[i] += 1;
      remote_read_to_host[graph->getHostID(gid)][i] += 1;
    }

    for (int i = cache_level; i < 11; i++) {
      mirror_read[i] += 1;
    }
  }

  void record_write_random(typename Graph::GraphNode node, bool comm) {
    auto gid = graph->getGID(node);
    if (graph->isOwned(gid)) { // master
      *master_write += 1;
      return;
    }
    // mirror
    int cache_level = node_cache_level[node];

    for (int i = 0; i < cache_level; i++) {
      remote_write[i] += 1;

      remote_write_to_host[graph->getHostID(gid)][i] += 1;
    }

    for (int i = cache_level; i < 11; i++) {
      mirror_write[i] += 1;

      if (comm) {
        remote_comm_to_host[graph->getHostID(gid)][i] += 1;
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
    file << "host " << host_id
         << " local read (stream): " << local_read_stream->read_local()
         << std::endl;
    file << "host " << host_id << " master reads: " << master_read->read_local()
         << std::endl;
    file << "host " << host_id
         << " master writes: " << master_write->read_local() << std::endl;

    for (int i = 0; i < 11; i++) {
      file << "host " << host_id << " cache " << i
           << " mirror reads: " << mirror_read[i].read_local() << std::endl;
      file << "host " << host_id << " cache " << i
           << " mirror writes: " << mirror_write[i].read_local() << std::endl;
      file << "host " << host_id << " cache " << i
           << " remote reads: " << remote_read[i].read_local() << std::endl;
      file << "host " << host_id << " cache " << i
           << " remote writes: " << remote_write[i].read_local() << std::endl;

      for (uint32_t j = 0; j < num_hosts; j++) {
        file << "host " << host_id << " cache " << i << " remote read to host "
             << j << ": " << remote_read_to_host[j][i].read_local()
             << std::endl;
        file << "host " << host_id << " cache " << i << " remote write to host "
             << j << ": " << remote_write_to_host[j][i].read_local()
             << std::endl;
        file << "host " << host_id << " cache " << i
             << " dirty mirrors for host " << j << ": "
             << remote_comm_to_host[j][i].read_local() << std::endl;
      }
    }
  }
};