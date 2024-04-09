#include <iostream>
#include <thread>
#include <mutex>
#include <galois/Timer.h>
#include "galois/wmd/graphTypes.h"


using namespace agile::workflow1;


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

  // disable move constructor
  graphUpdateManager(graphUpdateManager&&)            = delete;
  graphUpdateManager& operator=(graphUpdateManager&&) = delete;

  void start() {
    // start the dynamic changes
//    galois::runtime::evilPhase = galois::runtime::evilPhase + 1;
    startIngest = std::thread(&graphUpdateManager::ingest, this);
    checkThread = std::thread(&graphUpdateManager::checkForMessages, this);
    std::cout << "Started the dynamic changes\n";
  }

  void setBatchSize(uint64_t size) { batchSize = size; }

  uint64_t getBatchSize() { return batchSize; }

  void setPeriod(uint64_t period) { periodForCheck = period; }

  uint64_t getPeriod() { return periodForCheck; }

  bool stop() { 
    if(stopFlag) {
      startIngest.join();
      std::this_thread::sleep_for(std::chrono::milliseconds(2*periodForCheck));
      stopFlag2 = true;
      checkThread.join();
      std::cout << "Already stopped stopped " << graph->id << std::endl;
    }
    return stopFlag;
  }

private:
  std::thread checkThread;
  std::thread startIngest;
  uint64_t periodForCheck;
  std::string graphFile;
  T* graph;
  uint64_t batchSize = 1;
  bool stopFlag = false;
  bool stopFlag2 = false;

  template <typename N = NodeData, typename E = EdgeData>
  void processBatch(std::vector<std::string>& lines) {
    // process the batch
    for (auto& line : lines) {
      std::istringstream iss(line);
      //int opcode;
      //iss >> opcode;
      //switch (opcode) {
      //case 0: { // add vertex
      //  uint64_t src;
      //  iss >> src;
      //  NodeData data;
      //  //graph->addVertexTopologyOnly(src);
      //  break;
      //}
      //case 1: { // add edge
        std::string token;
        std::vector<std::string> tokens;
        while(std::getline(iss, token, ',')) {
          tokens.push_back(token);
        }
        Edge edge(tokens);
        // Inverse edge
        agile::workflow1::TYPES inverseEdgeType = agile::workflow1::TYPES::NONE;
        if (tokens[0] == "Sale") {
          inverseEdgeType = agile::workflow1::TYPES::PURCHASE;
        } else if (tokens[0] == "Author") {
          inverseEdgeType = agile::workflow1::TYPES::WRITTENBY;
        } else if (tokens[0] == "Includes") {
          inverseEdgeType = agile::workflow1::TYPES::INCLUDEDIN;
        } else if (tokens[0] == "HasTopic") {
          inverseEdgeType = agile::workflow1::TYPES::TOPICIN;
        } else if (tokens[0] == "HasOrg") {
          inverseEdgeType = agile::workflow1::TYPES::ORGIN;
        } else {
          assert(false);
        }
        agile::workflow1::Edge inverseEdge    = edge;
        inverseEdge.type = inverseEdgeType;
        std::swap(inverseEdge.src, inverseEdge.dst);
        std::swap(inverseEdge.src_type, inverseEdge.dst_type);
        std::vector<uint64_t> dsts;
      //  //if (iss.peek() == EOF) {
      //  //  // graph->addEdgesTopologyOnly(src, dsts);
      //  //  break;
      //  //} else if constexpr (!std::is_void<E>::value) {
        //while (iss.peek() != EOF) {
        //  uint64_t dst;
        //  iss >> dst;
        //  dsts.push_back(dst);
        //}
        dsts.push_back(edge.dst);
        std::vector<E> data;
        data.push_back(edge);
        //std::cout << "EDGE src " << edge.src << " dst " << edge.dst << " dsts sz " << dsts.size() << " data sz " << data.size() << std::endl;
        graph->addEdges(edge.src, dsts, data);
        dsts.clear();
        dsts.push_back(inverseEdge.dst);
        data.clear();
        data.push_back(inverseEdge);
        //std::cout << "EDGE inverse src " << edge.src << " dst " << edge.dst << " dsts sz " << dsts.size() << " data sz " << data.size() << std::endl;
        graph->addEdges(inverseEdge.src, dsts, data);
      //  break;
      //  //}
      //}
      //case 2: { // delete vertex
      //  uint64_t src;
      //  iss >> src;
      //  // graph->deleteVertex(src);
      //  break;
      //}
      //case 3: { // delete edge
      //  uint64_t src;
      //  std::vector<uint64_t> dsts;
      //  iss >> src;
      //  while (iss.peek() != EOF) {
      //    uint64_t dst;
      //    iss >> dst;
      //    dsts.push_back(dst);
      //  }
      //  // graph->deleteEdges(src, dsts);
      //  break;
      //}
      //}
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
    std::vector<std::string> lines;
    while ((std::getline(inputFile, line))) {
      lines.push_back(line);
      if (lines.size() == batchSize) {
        processBatch(lines);
        lines.clear();
        std::this_thread::sleep_for(std::chrono::milliseconds(periodForCheck));
      }
    }
    if(lines.size() > 0) {
      processBatch(lines);
    }
    inputFile.close();
    std::this_thread::sleep_for(std::chrono::milliseconds(periodForCheck));
    std::cout << "stopping flag " << graph->id << std::endl;  
    stopFlag = true;
  }

  template <typename N = NodeData, typename E = EdgeData>
  void checkForMessages() {
    // check for messages
    auto& net = galois::runtime::getSystemNetworkInterface();
    while (true && (!stopFlag2)) {
        auto m = net.recieveTagged(galois::runtime::evilPhase);
        if(m.has_value()) {
        typename T::Task m_task;
        uint64_t src_node;
        //galois::runtime::gDeserialize(m->second, m_task);
      //  std::cout << "Checking for messages dbg4\n";
        galois::runtime::gDeserialize(m->second, src_node);
      //  std::cout << "Checking for messages dbg5\n";
        //switch (m_task) {
        //case T::Task::ADD_VERTEX_TOPOLOGY_ONLY:
        //  std::cout << "Adding vertex topology only\n";
        //  // graph->addVertexTopologyOnly(src_node);
        //  break;
        //case T::Task::ADD_VERTEX: {
        //  std::cout << "Adding vertex\n";
        //  NodeData data;
        //  galois::runtime::gDeserialize(m->second, data);
        //  graph->addVertex(src_node, data);
        //  break;
        //}
        //case T::Task::ADD_EDGES_TOPOLOGY_ONLY: {
        //  std::cout << "Adding edges topology only\n";
        //  std::vector<uint64_t> dsts;
        //  galois::runtime::gDeserialize(m->second, dsts);
        //  graph->addEdgesTopologyOnly(src_node, dsts);
        //  break;
        //}
        //case T::Task::ADD_EDGES: {
          //if constexpr (!std::is_void<E>::value) {
          //  std::vector<E> edge_data;
          //  galois::runtime::gDeserialize(m->second, edge_data);
          //  std::vector<uint64_t> edge_dsts;
          //  galois::runtime::gDeserialize(m->second, edge_dsts);
          //  graph->addEdges(src_node, edge_dsts, edge_data);
          //} else {
            std::vector<uint64_t> edge_dsts;
            galois::runtime::gDeserialize(m->second, edge_dsts);
            std::vector<E> edge_data;
            galois::runtime::gDeserialize(m->second, edge_data);
            graph->addEdges(src_node, edge_dsts, edge_data);
          //}
         // break;
        //}
      //  case T::Task::DELETE_VERTEX: {
      //    std::cout << "Deleting vertex\n";
      //    graph->deleteVertex(src_node);
      //    break;
      //  }
      //  case T::Task::DELETE_EDGES: {
      //    std::cout << "Deleting edges\n";
      //    // std::vector<typename T::edge_iterator> edges;
      //    // galois::runtime::gDeserialize(m->second, edges);
      //    // graph->deleteEdges(src_node, edges);
      //    break;
      //  }
      //  default:
      //    std::cout << "Unknown task\n";
      //    std::cerr << "Unknown task\n";
      //    break;
      }
      std::this_thread::sleep_for(std::chrono::milliseconds(periodForCheck));
    }
  }
};
