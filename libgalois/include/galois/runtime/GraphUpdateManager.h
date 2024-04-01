#include <iostream>
#include <thread>
#include <chrono>
#include <galois/Timer.h>

template <typename NodeData, typename EdgeData>
class graphUpdateManager {
public:
  using T              = galois::graphs::DistLocalGraph<NodeData, EdgeData>;
  graphUpdateManager() = default;
  graphUpdateManager(std::string inputFile, int period, T* distGraphPtr) {
    periodForCheck = period;
    graphFile      = inputFile;
    graph          = distGraphPtr;
  }
  graphUpdateManager(int period) { periodForCheck = period; }
  // disable copy constructor
  graphUpdateManager(const graphUpdateManager&)            = delete;
  graphUpdateManager& operator=(const graphUpdateManager&) = delete;

  // diable move constructor
  graphUpdateManager(graphUpdateManager&&)            = delete;
  graphUpdateManager& operator=(graphUpdateManager&&) = delete;

  void start() {
    // start the dynamic changes
    startIngest = std::thread(&graphUpdateManager::ingest, this);
    checkThread = std::thread(&graphUpdateManager::checkForMessages, this);
  }

  void setBatchSize(uint64_t size) { batchSize = size; }

  uint64_t getBatchSize() { return batchSize; }

  void setPeriod(uint64_t period) { periodForCheck = period; }

  uint64_t getPeriod() { return periodForCheck; }

  void stop() { checkThread.join(); }

private:
  std::thread checkThread;
  std::thread startIngest;
  uint64_t periodForCheck;
  std::string graphFile;
  T* graph;
  uint64_t batchSize = 1000;

  template <typename N = NodeData, typename E = EdgeData>
  void processBatch(std::vector<std::string>& lines) {
    // process the batch
    for (auto& line : lines) {
      std::istringstream iss(line);
      int opcode;
      iss >> opcode;
      switch (opcode) {
      case 0: { // add vertex
        uint64_t src;
        iss >> src;
        NodeData data;
        // graph->addVertex(src, data);
        break;
      }
      case 1: { // add edge
        uint64_t src;
        std::vector<uint64_t> dsts;
        iss >> src;
        if (iss.peek() == EOF) {
          // graph->addEdgesTopologyOnly(src, dsts);
          break;
        } else if constexpr (!std::is_void<E>::value) {
          std::vector<E> data;
          // graph->addEdges(src, dsts, data);
          break;
        }
      }
      case 2: { // delete vertex
        uint64_t src;
        iss >> src;
        // graph->deleteVertex(src);
        break;
      }
      case 3: { // delete edge
        uint64_t src;
        std::vector<uint64_t> dsts;
        iss >> src;
        while (iss.peek() != EOF) {
          uint64_t dst;
          iss >> dst;
          dsts.push_back(dst);
        }
        // graph->deleteEdges(src, dsts);
        break;
      }
      }
    }
  }

  template <typename N = NodeData, typename E = EdgeData>
  void ingest() {
    // read the file in batches
    std::ifstream inputFile(graphFile);
    if (!inputFile.is_open()) {
      std::cerr << "Error opening file: " << graphFile << "\n";
      return;
    }

    std::string line;
    std::vector<std::string> lines(batchSize);
    while (std::getline(inputFile, line)) {
      lines.push_back(line);
      if (lines.size() == batchSize) {
        processBatch(lines);
        lines.clear();
        std::this_thread::sleep_for(std::chrono::milliseconds(periodForCheck));
      }
    }
  }

  template <typename N = NodeData, typename E = EdgeData>
  void checkForMessages() {
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
          // graph->addVertexTopologyOnly(src_node);
          break;
        case T::Task::ADD_VERTEX: {
          NodeData data;
          galois::runtime::gDeserialize(m->second, data);
          graph->addVertex(src_node, data);
          break;
        }
        case T::Task::ADD_EDGES_TOPOLOGY_ONLY: {
          std::vector<uint64_t> dsts;
          galois::runtime::gDeserialize(m->second, dsts);
          graph->addEdgesTopologyOnly(src_node, dsts);
          break;
        }
        case T::Task::ADD_EDGES: {
          if constexpr (!std::is_void<E>::value) {
            std::vector<E> edge_data;
            galois::runtime::gDeserialize(m->second, edge_data);
            std::vector<uint64_t> edge_dsts;
            galois::runtime::gDeserialize(m->second, edge_dsts);
            graph->addEdges(src_node, edge_dsts, edge_data);
          } else {
            std::vector<uint64_t> edge_dsts;
            galois::runtime::gDeserialize(m->second, edge_dsts);
            graph->addEdgesTopologyOnly(src_node, edge_dsts);
          }
          break;
        }
        case T::Task::DELETE_VERTEX: {
          graph->deleteVertex(src_node);
          break;
        }
        case T::Task::DELETE_EDGES: {
          // std::vector<typename T::edge_iterator> edges;
          // galois::runtime::gDeserialize(m->second, edges);
          // graph->deleteEdges(src_node, edges);
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
