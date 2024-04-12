#include <iostream>
#include <thread>
#include <mutex>
#include <galois/Timer.h>
#include "galois/wmd/graphTypes.h"


//Usage: call start() to start the ingestion of the file
//       call stop() to stop the ingestion of the file
//       call setBatchSize() to set the batch size
//   Refer to wmd-graph-build for an example of how to use this class

using namespace agile::workflow1;


template <typename NodeData, typename EdgeData>
class graphUpdateManager {
public:
  using T              = galois::graphs::DistLocalGraph<NodeData, EdgeData>;
  graphUpdateManager() = default;
  graphUpdateManager(std::unique_ptr<galois::graphs::FileParser<NodeData, EdgeData>> parser, std::string inputFile, int period, T* distGraphPtr) {
    periodForCheck = period;
    graphFile      = inputFile;
    graph          = distGraphPtr;
    fileParser     = std::move(parser);
  }
  // disable copy constructor
  graphUpdateManager(const graphUpdateManager&)            = delete;
  graphUpdateManager& operator=(const graphUpdateManager&) = delete;

  // disable move constructor
  graphUpdateManager(graphUpdateManager&&)            = delete;
  graphUpdateManager& operator=(graphUpdateManager&&) = delete;

  void start() {
    // start the dynamic changes
    startIngest = std::thread(&graphUpdateManager::ingestFile, this);
    checkThread = std::thread(&graphUpdateManager::checkForMessages, this);
  }

  void setBatchSize(uint64_t size) { batchSize = size; }

  uint64_t getBatchSize() { return batchSize; }

  void setPeriod(uint64_t period) { periodForCheck = period; }

  uint64_t getPeriod() { return periodForCheck; }

  bool stop() { 
    if(stopIngest) {
      startIngest.join();
      std::this_thread::sleep_for(std::chrono::milliseconds(2*periodForCheck));
      stopCheck = true;
      checkThread.join();
    }
    return stopIngest;
  }

private:
  std::thread checkThread;
  std::thread startIngest;
  uint64_t periodForCheck;
  std::string graphFile;
  T* graph;
  uint64_t batchSize = 1;
  bool stopIngest = false;
  bool stopCheck = false;
  std::unique_ptr<galois::graphs::FileParser<NodeData, EdgeData>> fileParser;

  template <typename N = NodeData, typename E = EdgeData>
  void processLine(const char* line, size_t len) {
    galois::graphs::ParsedGraphStructure<N, E> value = fileParser->ParseLine(const_cast<char*>(line), len);
    for (auto& edge : value.edges) {
      std::vector<uint64_t> dsts;
      dsts.push_back(edge.dst);
      std::vector<E> data;
      data.push_back(edge);
      graph->addEdges(edge.src, dsts, data);
    }
  }

  template <typename N = NodeData, typename E = EdgeData>
  void ingestFile() {
    std::ifstream inputFile(graphFile);
    if (!inputFile.is_open()) {
      std::cerr << "Error opening file: " << graphFile << "\n";
      return;
    }

    // Read each line from the stringstream
    std::string line;
    uint64_t lineNumber = 0;
    while ((std::getline(inputFile, line))) {
      processLine(line.c_str(), line.size());
      lineNumber++;
      if(lineNumber == batchSize) {
        std::this_thread::sleep_for(std::chrono::milliseconds(periodForCheck));
        lineNumber = 0;
      }
    }
    inputFile.close();
    stopIngest = true;
  }

  template <typename N = NodeData, typename E = EdgeData>
  void checkForMessages() {
    // check for messages
    auto& net = galois::runtime::getSystemNetworkInterface();
    while (!stopCheck) {
        auto m = net.recieveTagged(galois::runtime::evilPhase);
        if(m.has_value()) {
        uint64_t src_node;
        galois::runtime::gDeserialize(m->second, src_node);
            std::vector<uint64_t> edge_dsts;
            galois::runtime::gDeserialize(m->second, edge_dsts);
            std::vector<E> edge_data;
            galois::runtime::gDeserialize(m->second, edge_data);
            graph->addEdges(src_node, edge_dsts, edge_data);
      }
      std::this_thread::sleep_for(std::chrono::milliseconds(periodForCheck));
    }
  }
};
