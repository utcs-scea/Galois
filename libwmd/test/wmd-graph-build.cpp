/*
 * Run this script in Debug view and compare its result to result of
 * MiningPartitioner.h
 *
 * A testing result sheet could be found at
 * https://docs.google.com/spreadsheets/u/1/d/1D0dAab29uazRKVroBZdEvAsEUT92aNHsiIYuHK6dlDA
 *
 * TODO: include a script to gen dataset and compare result
 *
 */
#include "galois/wmd/graph.h"
#include "galois/wmd/WMDPartitioner.h"
#include "galois/shad/DataTypes.h"
#include "galois/wmd/graphTypes.h"

#include "galois/DistGalois.h"
#include "galois/graphs/GenericPartitioners.h"
#include "galois/runtime/GraphUpdateManager.h"

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <boost/program_options.hpp>

namespace po = boost::program_options;

using namespace agile::workflow1;

typedef galois::graphs::WMDGraph<Vertex, Edge, OECPolicy> Graph;

void insertEdge(
    Edge edge,
    std::unordered_map<std::uint64_t, std::pair<TYPES, std::vector<Edge>>>&
        vertices) {
  if (vertices.find(edge.src) != vertices.end()) {
    vertices[edge.src].second.push_back(edge);
  } else {
    assert(false);
  }
}

void parser(std::string line,
            std::unordered_map<std::uint64_t,
                               std::pair<TYPES, std::vector<Edge>>>& vertices) {
  if (line.find("//") != std::string::npos ||
      line.find("#") != std::string::npos) {
    return;
  } else if (line.find("/*") != std::string::npos ||
             line.find("*/") != std::string::npos) {
    return;
  } else {
    const char* ptr = line.c_str();
    std::istringstream ss(ptr);
    std::string token;
    std::vector<std::string> tokens;
    while (std::getline(ss, token, ',')) {
      tokens.push_back(token);
    }
    if (tokens.size() == 9)
      tokens.push_back("");
    if (tokens.size() == 0)
      return;
    assert(tokens.size() == 10);
    bool isNode = tokens[0] == "Person" || tokens[0] == "ForumEvent" ||
                  tokens[0] == "Forum" || tokens[0] == "Publication" ||
                  tokens[0] == "Topic";
    if (isNode) {
      uint64_t id                        = 0;
      agile::workflow1::TYPES vertexType = agile::workflow1::TYPES::NONE;
      if (tokens[0] == "Person") {
        vertexType = agile::workflow1::TYPES::PERSON;
        id         = std::stoull(tokens[1]);
      } else if (tokens[0] == "ForumEvent") {
        vertexType = agile::workflow1::TYPES::FORUMEVENT;
        id         = std::stoull(tokens[4]);
      } else if (tokens[0] == "Forum") {
        vertexType = agile::workflow1::TYPES::FORUM;
        id         = std::stoull(tokens[3]);
      } else if (tokens[0] == "Publication") {
        vertexType = agile::workflow1::TYPES::PUBLICATION;
        id         = std::stoull(tokens[5]);
      } else if (tokens[0] == "Topic") {
        vertexType = agile::workflow1::TYPES::TOPIC;
        id         = std::stoull(tokens[6]);
      } else {
        assert(false);
      }
      vertices[id] =
          std::pair<TYPES, std::vector<Edge>>(vertexType, std::vector<Edge>());
    } else {
      Edge edge(tokens);
      insertEdge(edge, vertices);
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
      agile::workflow1::Edge inverseEdge = edge;
      inverseEdge.type                   = inverseEdgeType;
      std::swap(inverseEdge.src, inverseEdge.dst);
      std::swap(inverseEdge.src_type, inverseEdge.dst_type);
      insertEdge(inverseEdge, vertices);
    }
  }
}

void getDataFromGraph(
    std::string& filename,
    std::unordered_map<std::uint64_t, std::pair<TYPES, std::vector<Edge>>>&
        vertices) {
  // read file line by line
  std::string line;
  std::ifstream myfile(filename);
  if (myfile.is_open()) {
    while (getline(myfile, line)) {
      parser(line, vertices);
    }
    myfile.close();
  } else {
    std::cout << "Unable to open file";
  }
}

int main(int argc, char* argv[]) {
  std::string dataFile;
  std::string dynFile;
  int threads;
  po::options_description desc("Allowed options");
  desc.add_options()("help", "print help info")(
      "staticFile", po::value<std::string>(&dataFile)->required(),
      "Input file for initial static graph")(
      "dynFile", po::value<std::string>(&dynFile)->default_value(""),
      "Input file for dynamic graph")(
      "threads", po::value<int>(&threads)->default_value(1),
      "Number of threads");

  po::variables_map vm;
  try {
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  }

  if (vm.count("help")) {
    std::cout << desc << "\n";
    return 1;
  }

  galois::DistMemSys G; // init galois memory
  auto& net = galois::runtime::getSystemNetworkInterface();
  galois::setActiveThreads(threads);

  if (net.ID == 0) {
    galois::gPrint("Testing building WMD graph from file.\n");
    galois::gPrint("Num Hosts: ", net.Num,
                   ", Active Threads Per Hosts: ", galois::getActiveThreads(),
                   "\n");
  }

  std::string file = dataFile;
  std::vector<std::string> filenames;
  filenames.emplace_back(dataFile);
  std::vector<std::unique_ptr<galois::graphs::FileParser<
      agile::workflow1::Vertex, agile::workflow1::Edge>>>
      parsers;
  parsers.emplace_back(
      std::make_unique<galois::graphs::WMDParser<agile::workflow1::Vertex,
                                                 agile::workflow1::Edge>>(
          10, filenames));
  Graph* graph = new Graph(parsers, net.ID, net.Num, true, false,
                           galois::graphs::BALANCED_EDGES_OF_MASTERS);
  assert(graph != nullptr);

  std::unordered_map<std::uint64_t, std::pair<TYPES, std::vector<Edge>>>
      vertices;

  if (dynFile != "") {
    std::string dynamicFile = dynFile + std::to_string(net.ID) + ".txt";

    graphUpdateManager<agile::workflow1::Vertex, agile::workflow1::Edge> GUM(
        std::make_unique<galois::graphs::WMDParser<agile::workflow1::Vertex,
                                                   agile::workflow1::Edge>>(
            10, filenames),
        dynamicFile, 100, graph);
    GUM.start();
    // wait for GUM to finish
    while (!GUM.stop()) {
      std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }
    galois::runtime::getHostBarrier().wait();
    GUM.stop2();
  }

  // generate a file with sorted token of all nodes and its outgoing edge dst
  // compare it with other implementation to verify the correctness
  std::vector<std::pair<uint64_t, std::vector<uint64_t>>> tokenAndEdges;
  tokenAndEdges.resize(graph->numMasters());
  galois::do_all(
      galois::iterate(graph->masterNodesRange()),
      [&](size_t lid) {
        auto token = graph->getData(lid).id;
        std::vector<uint64_t> edgeDst;
        auto end = graph->edge_end(lid);
        auto itr = graph->edge_begin(lid);
        for (; itr != end; itr++) {
          edgeDst.push_back(graph->getEdgeDst(itr));
        }
        std::vector<uint64_t> edgeDstDbg;
        for (auto& e : graph->edges(lid)) {
          edgeDstDbg.push_back(graph->getEdgeDst(e));
        }
        assert(edgeDst == edgeDstDbg);
        std::sort(edgeDst.begin(), edgeDst.end());
        tokenAndEdges[lid] = std::make_pair(token, std::move(edgeDst));
      },
      galois::steal());
  // gather node info from other hosts
  if (net.ID != 0) { // send token and degree pairs to host 0
    galois::runtime::SendBuffer sendBuffer;
    galois::runtime::gSerialize(sendBuffer, tokenAndEdges);
    net.sendTagged(0, galois::runtime::evilPhase, std::move(sendBuffer));
  } else { // recv node range from other hosts
    for (size_t i = 0; i < net.Num - 1; i++) {
      decltype(net.recieveTagged(galois::runtime::evilPhase)) p;
      do {
        p = net.recieveTagged(galois::runtime::evilPhase);
      } while (!p);
      std::vector<std::pair<uint64_t, std::vector<uint64_t>>>
          incomingtokenAndEdges;
      galois::runtime::gDeserialize(p->second, incomingtokenAndEdges);
      // combine data
      std::move(incomingtokenAndEdges.begin(), incomingtokenAndEdges.end(),
                std::back_inserter(tokenAndEdges));
    }
  }
  if (net.ID == 0) {
    getDataFromGraph(file, vertices);
    if (dynFile != "") {
      for (uint32_t i = 0; i < net.Num; i++) {
        std::string dynamicFile = dynFile + std::to_string(i) + ".txt";
        getDataFromGraph(dynamicFile, vertices);
      }
    }
    // compare with vertices
    assert(tokenAndEdges.size() == vertices.size());
    for (size_t i = 0; i < tokenAndEdges.size(); i++) {
      auto& tokenAndEdge = tokenAndEdges[i];
      auto& vertex       = vertices[tokenAndEdge.first];
      assert(vertex.second.size() == tokenAndEdge.second.size());
      std::sort(vertex.second.begin(), vertex.second.end(),
                [](const agile::workflow1::Edge& a,
                   const agile::workflow1::Edge& b) { return a.dst < b.dst; });
      for (size_t j = 0; j < vertex.second.size(); j++) {
        assert(vertex.second[j].dst == tokenAndEdge.second[j]);
      }
    }
  }
  return 0;
}
