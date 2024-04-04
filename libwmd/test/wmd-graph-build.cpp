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

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <unordered_map>

using namespace agile::workflow1;

typedef galois::graphs::WMDGraph<Vertex,
                                 Edge, OECPolicy>
    Graph;


struct Vrtx {
  uint64_t id;
  agile::workflow1::TYPES type;
  std::vector<agile::workflow1::Edge> edges;

  Vrtx(uint64_t id, agile::workflow1::TYPES type) : id(id), type(type) {}
  Vrtx() : id(0), type(agile::workflow1::TYPES::NONE) {}
  
};

void parser(std::string line, std::unordered_map<std::uint64_t, Vrtx>& vertices) {
  if (line.find("//") != std::string::npos || line.find("#") != std::string::npos) {
      return;
    } else if (line.find("/*") != std::string::npos || line.find("*/") != std::string::npos) {
      return;
    } else {
      const char* ptr = line.c_str();
      std::istringstream ss(ptr);
      std::string token;
      std::vector<std::string> tokens;
      while(std::getline(ss, token, ',')) {
        tokens.push_back(token);
      }
      bool isNode = tokens[0] == "Person" ||
                    tokens[0] == "ForumEvent" ||
                    tokens[0] == "Forum" ||
                    tokens[0] == "Publication" ||
                    tokens[0] == "Topic";
      if (isNode) {
        uint64_t id;
        agile::workflow1::TYPES vertexType;
        if (tokens[0] == "Person") {
          vertexType = agile::workflow1::TYPES::PERSON;
          id = std::stoull(tokens[1]);
        } else if (tokens[0] == "ForumEvent") {
          vertexType = agile::workflow1::TYPES::FORUMEVENT;
          id = std::stoull(tokens[4]);
        } else if (tokens[0] == "Forum") {
          vertexType = agile::workflow1::TYPES::FORUM;
          id = std::stoull(tokens[3]);
        } else if (tokens[0] == "Publication") {
          vertexType = agile::workflow1::TYPES::PUBLICATION;
          id = std::stoull(tokens[5]);
        } else if (tokens[0] == "Topic") {
          vertexType = agile::workflow1::TYPES::TOPIC;
          id = std::stoull(tokens[6]);
        } else {
          assert(false);
        }
        Vrtx vertex(id, vertexType);
        vertices.insert({vertex.id, vertex});
      } else {
        Edge edge(tokens);
        Vrtx srcVertex;
        if(vertices.find(edge.src) != vertices.end()) {
          srcVertex = vertices[edge.src];
          vertices[edge.src].edges.push_back(edge);
        } else {
          assert(false);
        }
        if(vertices.find(edge.src) != vertices.end()) {
          vertices.insert({edge.src, srcVertex});
        } else {
          assert(false);
        }
        // Inverse edge
        agile::workflow1::TYPES inverseEdgeType;
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
        Vrtx dstVertex;
        if(vertices.find(edge.dst) != vertices.end()) {
          dstVertex = vertices[edge.dst];
          vertices[edge.dst].edges.push_back(inverseEdge);
        } else {
          assert(false);
        }
        if(vertices.find(edge.dst) != vertices.end()) {
          vertices.insert({edge.dst, dstVertex});
        } else {
          assert(false);
        }
      }
  }

}

void getDataFromGraph(std::string& filename, std::unordered_map<std::uint64_t, Vrtx>& vertices) {
  //read file line by line
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
  galois::DistMemSys G; // init galois memory
  auto& net = galois::runtime::getSystemNetworkInterface();

  if (argc == 3)
    galois::setActiveThreads(atoi(argv[2]));

  if (net.ID == 0) {
    galois::gPrint("Testing building WMD graph from file.\n");
    galois::gPrint("Num Hosts: ", net.Num,
                   ", Active Threads Per Hosts: ", galois::getActiveThreads(),
                   "\n");
  }

  std::string dataFile = argv[1];
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
  I_INIT("tmp/WMDGraph", net.ID, net.Num, 0);
  Graph* graph = new Graph(parsers, net.ID, net.Num, true, false,
                           galois::graphs::BALANCED_EDGES_OF_MASTERS);
  I_ROUND(0);
  I_CLEAR();
  assert(graph != nullptr);

  std::unordered_map<std::uint64_t, Vrtx> vertices;
  getDataFromGraph(file, vertices);

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
        for(auto& e : graph->edges(lid)) {
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

  // sort the node info by token order
  // serilize it to file
  if (net.ID == 0) {
    //compare with vertices
    assert(tokenAndEdges.size() == vertices.size());
    for (size_t i = 0; i < tokenAndEdges.size(); i++) {
      auto& tokenAndEdge = tokenAndEdges[i];
      auto& vertex = vertices[tokenAndEdge.first];
      assert(vertex.id == tokenAndEdge.first);
      assert(vertex.edges.size() == tokenAndEdge.second.size());
      std::sort(vertex.edges.begin(), vertex.edges.end(),
           [](const agile::workflow1::Edge& a, const agile::workflow1::Edge& b) {
             return a.dst < b.dst;
           });
      for (size_t j = 0; j < vertex.edges.size(); j++) {
        assert(vertex.edges[j].dst == tokenAndEdge.second[j]);
      }
    }
  }
  I_DEINIT();
  return 0;
}
