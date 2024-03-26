#include <iostream>
#include <thread>
#include <chrono>
#include <galois/Timer.h>

template <typename NodeData, typename EdgeData>
class graphUpdateManager {
public:
  graphUpdateManager() = default;
  graphUpdateManager(std::string inputFile, int period) {
    periodForCheck = period;
    graphFile      = inputFile;
  }
  graphUpdateManager(int period) { periodForCheck = period; }
  // disable copy constructor
  graphUpdateManager(const graphUpdateManager&)            = delete;
  graphUpdateManager& operator=(const graphUpdateManager&) = delete;

  // diable move constructor
  graphUpdateManager(graphUpdateManager&&)            = delete;
  graphUpdateManager& operator=(graphUpdateManager&&) = delete;

  using T = galois::graphs::DistGraph<NodeData, EdgeData>;

  void start(T& graph) {
    // start the dynamic changes
    startIngest();
    checkThread =
        std::thread(&graphUpdateManager::checkForMessages, this, graph);
  }

  void stop() { checkThread.join(); }

private:
  std::thread checkThread;
  uint64_t periodForCheck;
  std::string graphFile;

  void startIngest() {
    // TODO: start the ingestion thread in batches
  }

  template <typename N = NodeData, typename E = EdgeData>
  void checkForMessages(T& graph) {
    // check for messages
    auto& net = galois::runtime::getSystemNetworkInterface();
    while (true) {
      std::this_thread::sleep_for(std::chrono::milliseconds(periodForCheck));
      while (net.anyPendingReceives()) {
        auto m = net.recieveTagged(galois::runtime::evilPhase);
        typename T::Task m_task;
        uint64_t src_node;
        galois::runtime::gDeserialize(m->second, m_task);
        galois::runtime::gDeserialize(m->second, src_node);
        switch (m_task) {
        case T::Task::ADD_VERTEX_TOPOLOGY_ONLY:
          graph.addVertexTopologyOnly(src_node);
          break;
        case T::Task::ADD_VERTEX: {
          NodeData data;
          galois::runtime::gDeserialize(m->second, data);
          graph.addVertex(src_node, data);
          break;
        }
        case T::Task::ADD_EDGES_TOPOLOGY_ONLY: {
          std::vector<uint32_t> dsts;
          galois::runtime::gDeserialize(m->second, dsts);
          graph.addEdgesTopologyOnly(src_node, dsts);
          break;
        }
        case T::Task::ADD_EDGES: {
          if constexpr (!std::is_void<E>::value) {
            std::vector<E> edge_data;
            galois::runtime::gDeserialize(m->second, edge_data);
            std::vector<uint32_t> edge_dsts;
            galois::runtime::gDeserialize(m->second, edge_dsts);
            graph.addEdges(src_node, edge_dsts, edge_data);
          } else {
            std::vector<uint32_t> edge_dsts;
            galois::runtime::gDeserialize(m->second, edge_dsts);
            graph.addEdgesTopologyOnly(src_node, edge_dsts);
          }
          break;
        }
        case T::Task::DELETE_VERTEX: {
          graph.deleteVertex(src_node);
          break;
        }
        case T::Task::DELETE_EDGES: {
          std::vector<typename T::edge_iterator> edges;
          galois::runtime::gDeserialize(m->second, edges);
          graph.deleteEdges(src_node, edges);
          break;
        }
        default:
          std::cerr << "Unknown task\n";
          break;
        }
      }
    }
  }
};
