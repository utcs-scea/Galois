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
#include "galois/wmd/WMDPartitioner.h"
#include "galois/shad/DataTypes.h"
#include "galois/wmd/graphTypes.h"

#include "galois/DistGalois.h"
#include "galois/graphs/GenericPartitioners.h"

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <unordered_map>

struct NodeData {
  uint32_t dist_current;
};

typedef galois::graphs::WMDGraph<
    galois::graphs::ELVertex, galois::graphs::ELEdge, NodeData, void, OECPolicy>
    Graph;

const char* elGetOne(const char* line, std::uint64_t& val) {
  bool found = false;
  val        = 0;
  char c;
  while ((c = *line++) != '\0' && isspace(c)) {
  }
  do {
    if (isdigit(c)) {
      found = true;
      val *= 10;
      val += (c - '0');
    } else if (c == '_') {
      continue;
    } else {
      break;
    }
  } while ((c = *line++) != '\0' && !isspace(c));
  if (!found)
    val = UINT64_MAX;
  return line;
}

void parser(
    const char* line,
    std::unordered_map<std::uint64_t, std::vector<uint64_t>>& vertices) {
  uint64_t src, dst;
  line = elGetOne(line, src);
  line = elGetOne(line, dst);
  if (vertices.find(src) == vertices.end()) {
    std::vector<uint64_t> edges;
    vertices.insert(std::make_pair(src, std::vector<uint64_t>{dst}));
  } else {
    vertices[src].push_back(dst);
  }
}

void getELDataFromGraph(
    std::string& filename,
    std::unordered_map<std::uint64_t, std::vector<uint64_t>>& vertices) {
  // read file line by line
  std::string line;
  std::ifstream myfile(filename);
  if (myfile.is_open()) {
    std::string line;
    while (getline(myfile, line)) {
      parser(line.c_str(), vertices);
    }
    myfile.close();
  } else {
    std::cout << "Unable to open file";
  }
}

int main(int argc, char* argv[]) {
  galois::DistMemSys G; // init galois memory
  auto& net = galois::runtime::getSystemNetworkInterface();

  if (argc == 5)
    galois::setActiveThreads(atoi(argv[4]));

  if (std::strcmp(argv[2], "--numVertices") != 0) {
    std::cerr << "Usage: " << argv[2] << " --numVertices <value>" << std::endl;
    return 1;
  }

  uint64_t numVertices = std::stoull(argv[3]);

  if (net.ID == 0) {
    galois::gPrint("Testing building WMD graph from file.\n");
    galois::gPrint("Num Hosts: ", net.Num,
                   ", Active Threads Per Hosts: ", galois::getActiveThreads(),
                   "\n");
  }

  std::string dataFile = argv[1];
  std::string file     = dataFile;
  std::vector<std::string> filenames;
  filenames.emplace_back(dataFile);
  std::unordered_map<std::uint64_t, std::vector<uint64_t>> vertices;
  getELDataFromGraph(file, vertices);

  std::vector<std::unique_ptr<galois::graphs::FileParser<
      galois::graphs::ELVertex, galois::graphs::ELEdge>>>
      parsers;
  parsers.emplace_back(
      std::make_unique<galois::graphs::ELParser<galois::graphs::ELVertex,
                                                galois::graphs::ELEdge>>(
          filenames));
  Graph* graph = new Graph(parsers, net.ID, net.Num, true, false, numVertices,
                           galois::graphs::BALANCED_EDGES_OF_MASTERS);
  assert(graph != nullptr);

  // generate a file with sorted token of all nodes and its outgoing edge dst
  // compare it with other implementation to verify the correctness
  std::vector<std::pair<uint64_t, std::vector<uint64_t>>> tokenAndEdges;
  tokenAndEdges.resize(graph->numMasters());
  galois::do_all(
      galois::iterate(graph->masterNodesRange()),
      [&](size_t lid) {
        auto token = graph->getGID(lid);
        std::vector<uint64_t> edgeDst;
        auto end = graph->edge_end(lid);
        auto itr = graph->edge_begin(lid);
        for (; itr != end; itr++) {
          edgeDst.push_back(graph->getGID(graph->getEdgeDst(itr)));
        }
        std::vector<uint64_t> edgeDstDbg;
        for (auto& e : graph->edges(lid)) {
          edgeDstDbg.push_back(graph->getGID(graph->getEdgeDst(e)));
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
  // sort the node info by token order
  // serilize it to file
  if (net.ID == 0) {
    // compare with vertices
    for (size_t i = 0; i < tokenAndEdges.size(); i++) {
      auto& tokenAndEdge = tokenAndEdges[i];
      auto& vertex       = vertices[tokenAndEdge.first];
      assert(vertices.find(tokenAndEdge.first) != vertices.end());
      assert(vertex.size() == tokenAndEdge.second.size());
      std::sort(vertex.begin(), vertex.end());
      for (size_t j = 0; j < vertex.size(); j++) {
        std::cout << "src " << tokenAndEdge.first << " " << vertex[j]
                  << " tokenAndEdge: " << tokenAndEdge.second[j] << std::endl;
        assert(vertex[j] == tokenAndEdge.second[j]);
      }
    }
  }
  return 0;
}
