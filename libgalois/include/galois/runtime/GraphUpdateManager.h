#include <iostream>
#include <thread>
#include <atomic>
#include <galois/Timer.h>
#include "galois/wmd/graphTypes.h"

// Usage: call start() to start the ingestion of the file
//        call stop() to stop the ingestion of the file
//        call setBatchSize() to set the batch size
//    Refer to wmd-graph-build for an example of how to use this class

using namespace agile::workflow1;

template <typename NodeData, typename EdgeData, typename NodeTy,
          typename EdgeTy, typename OECPolicy>
class graphUpdateManager {
public:
  using T =
      galois::graphs::WMDGraph<NodeData, EdgeData, NodeTy, EdgeTy, OECPolicy>;
  graphUpdateManager() = default;
  graphUpdateManager(
      std::unique_ptr<galois::graphs::FileParser<NodeData, EdgeData>> parser,
      int period, T* distGraphPtr) {
    periodForCheck = period;
    graph          = distGraphPtr;
    fileParser     = std::move(parser);
  }
  // disable copy constructor
  graphUpdateManager(const graphUpdateManager&)            = delete;
  graphUpdateManager& operator=(const graphUpdateManager&) = delete;

  // disable move constructor
  graphUpdateManager(graphUpdateManager&&)            = delete;
  graphUpdateManager& operator=(graphUpdateManager&&) = delete;

  void update() {
    ingestFile();
    graph->updateRanges();
  }

  void setBatchSize(uint64_t size) { batchSize = size; }

  uint64_t getBatchSize() { return batchSize; }

  void setPeriod(uint64_t period) { periodForCheck = period; }

  uint64_t getPeriod() { return periodForCheck; }

private:
  std::thread checkThread;
  std::thread startIngest;
  uint64_t periodForCheck;
  std::string graphFile;
  T* graph;
  uint64_t batchSize = 10;
  bool stopIngest    = false;
  bool stopCheck     = false;
  std::unique_ptr<galois::graphs::FileParser<NodeData, EdgeData>> fileParser;

  void processNodes(std::vector<NodeData>& nodes) {
    for(auto& node : nodes) {
      graph->addVertexTopologyOnly(node.id);
    }
  }

  void processEdges(std::vector<EdgeData>& edges) {
    for(auto& edge : edges) {
      std::vector<uint64_t> dsts;
      dsts.push_back(edge.dst);
      graph->addEdgesTopologyOnly(edge.src, dsts);
    }
  }

  template <typename N = NodeData, typename E = EdgeData>
  void processUpdates(std::vector<galois::graphs::ParsedGraphStructure<N, E>>& updateVector) {
    std::vector<std::vector<EdgeData>> updateEdges;
    std::vector<std::vector<NodeData>> updateNodes;
    auto& net = galois::runtime::getSystemNetworkInterface();
    updateNodes.resize(net.Num);
    updateEdges.resize(net.Num);
    uint64_t numEdges = 0;
    for(auto& update : updateVector) {
        if(update.isNode) {
          if(graph->isOwned(update.node.id)) {
            graph->addVertexTopologyOnly(update.node.id);
          } else {
            updateNodes[graph->getHostID(update.node.id)].push_back(update.node);
          }
        } else {
          for(auto& edge : update.edges) {
            if(graph->isOwned(edge.src)) {
              std::vector<E> edges;
              edges.push_back(edge);
              processEdges(edges);
            } else {
              updateEdges[graph->getHostID(edge.src)].push_back(edge);
            }
          }
        }
    }

    //send updates to other hosts
    for(uint32_t i = 0; i < net.Num; i++) {
      if(i == net.ID) {
        continue;
      }
      galois::runtime::SendBuffer b;
      galois::runtime::gSerialize(b, updateNodes[i]);
      net.sendTagged(i, galois::runtime::evilPhase, std::move(b));
    }
    for(uint32_t i=0; i<net.Num-1; i++) {
      decltype(net.recieveTagged(galois::runtime::evilPhase)) p;
      do {
        p = net.recieveTagged(galois::runtime::evilPhase);
      } while (!p);
      std::vector<N> recvNodes;
      galois::runtime::gDeserialize(p->second, recvNodes);
      processNodes(recvNodes); 
    }
    galois::runtime::evilPhase++;

    for(uint32_t i = 0; i < net.Num; i++) {
      if(i == net.ID) {
        continue;
      }
      galois::runtime::SendBuffer b;
      galois::runtime::gSerialize(b, updateEdges[i]);
      net.sendTagged(i, galois::runtime::evilPhase, std::move(b));
    }
    for(uint32_t i=0; i<net.Num-1; i++) {
      decltype(net.recieveTagged(galois::runtime::evilPhase)) p;
      do {
        p = net.recieveTagged(galois::runtime::evilPhase);
      } while (!p);
      std::vector<E> recvEdges;
      galois::runtime::gDeserialize(p->second, recvEdges);
      processEdges(recvEdges); 
    }
    galois::runtime::evilPhase++;
  }

  template <typename N = NodeData, typename E = EdgeData>
  void ingestFile() {
    std::vector<std::string> files = fileParser->GetFiles();
    for (auto& file : files) {
      std::ifstream inputFile(file);
      if (!inputFile.is_open()) {
        std::cerr << "Error opening file: " << graphFile << "\n";
        return;
      }

      // Read each line from the stringstream
      std::string line;
      std::vector<galois::graphs::ParsedGraphStructure<N, E>> parsedData;
      while ((std::getline(inputFile, line))) {
        parsedData.push_back(fileParser->ParseLine(const_cast<char*>(line.c_str()), line.size()));
      }
      processUpdates(parsedData);
      inputFile.close();
    }
    auto& net = galois::runtime::getSystemNetworkInterface();
    net.flush();
  }

};
