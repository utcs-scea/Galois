#include "galois/Logging.h"
#include "GNNBench/Input.h"

namespace cll = llvm::cl;

// Self documented via the desc argument

llvm::cl::opt<std::string> input_directory(
    "inputDirectory",
    cll::desc("Base directory to find all files required for doing GNN "
              "training (features, graph topology, masks, etc.)"),
    cll::init(galois::default_gnn_dataset_path));

llvm::cl::opt<std::string> input_name(
    cll::Positional,
    cll::desc("Base name of graph: used to find csgr, features, etc."),
    cll::Required);

llvm::cl::opt<galois::graphs::GNNPartitionScheme> partition_scheme(
    "partition", cll::desc("Type of partitioning."),
    cll::values(clEnumValN(galois::graphs::GNNPartitionScheme::kOEC, "oec",
                           "Outgoing Edge-Cut (default)"),
                clEnumValN(galois::graphs::GNNPartitionScheme::kCVC, "cvc",
                           "Cartesian Vertex-Cut"),
                clEnumValN(galois::graphs::GNNPartitionScheme::kOCVC, "ocvc",
                           "Original Cartesian Vertex-Cut")),
    cll::init(galois::graphs::GNNPartitionScheme::kOEC));

cll::opt<bool> useShad("useShad", cll::desc("true if the input graph is"
                                            " SHAD WMD graph format."
                                            " Otheriwse, set false."),
                       cll::init(false));

llvm::cl::opt<unsigned> num_layers(
    "numLayers",
    cll::desc(
        "Number of intermediate layers in the neural network (default 2))"),
    cll::init(2));

// llvm::cl::list<unsigned> layer_sizes(
//    "layerSizes",
//    cll::desc(
//        "Comma separated list of numbers specifying "
//        "intermediate layer sizes (does not include output); default sizes are
//        " "16 until last layer which is the size of the # of labels"),
//    cll::CommaSeparated);

llvm::cl::opt<unsigned> layer_size(
    "layerSize",
    cll::desc(
        "Number specifying "
        "intermediate layer sizes (does not include output); default sizes are "
        "16 until last layer which is the size of the # of labels"),
    cll::init(16));

llvm::cl::opt<galois::GNNLayerType> cl_layer_type(
    "layerType",
    cll::desc("Layer type specifying "
              "intermediate layers (does not include output); default SAGE"),
    cll::values(
        clEnumValN(galois::GNNLayerType::kGraphConvolutional, "gcn",
                   "Graph Convolutional Layer (default)"),
        clEnumValN(galois::GNNLayerType::kSAGE, "sage",
                   "SAGE layer (GCN with concat + mean)"),
        clEnumValN(galois::GNNLayerType::kL2Norm, "l2norm", "L2 norm layer"),
        clEnumValN(galois::GNNLayerType::kDense, "dense", "Dense layer")),
    cll::init(galois::GNNLayerType::kSAGE));

llvm::cl::list<unsigned> cl_fan_out_vector(
    "samplingFanOut",
    cll::desc(
        "Comma separated list of layer fanout if sampling/batching is used"),
    cll::CommaSeparated);

llvm::cl::opt<bool>
    disable_dropout("disableDropout",
                    cll::desc("If true (off by default), disables dropout of "
                              "layer weights during training"),
                    cll::init(false));

llvm::cl::opt<float> dropout_rate(
    "dropoutRate",
    cll::desc("Specifies probability that any one weight is DROPPED (e.g., if "
              "0.1, then 10 percent chance of dropping) (default 0.5)"),
    cll::init(0.5));

llvm::cl::opt<bool> disable_activation(
    "disableActivation",
    cll::desc("If true (off by default), disable activation at the "
              "end of an intermediate layers"),
    cll::init(false));

llvm::cl::opt<bool> disable_normalization(
    "disableNormalization",
    cll::desc("If true (off by default), disable normalizing vertex "
              "features based on their degree"),
    cll::init(false));

llvm::cl::opt<galois::GNNOutputLayerType> output_layer_type(
    "outputLayer", cll::desc("Type of output layer"),
    cll::values(clEnumValN(galois::GNNOutputLayerType::kSoftmax, "softmax",
                           "Softmax (default)"),
                clEnumValN(galois::GNNOutputLayerType::kSigmoid, "sigmoid",
                           "Sigmoid")),
    cll::init(galois::GNNOutputLayerType::kSoftmax));

llvm::cl::opt<bool>
    multiclass_labels("multiclassLabels",
                      cll::desc("If true (off by default), use multi-class "
                                "ground truth; required for some inputs"),
                      cll::init(false));

llvm::cl::opt<bool> disable_agg_after_update(
    "disableAggregationAfterUpdate",
    cll::desc("If true (off by default), disables aggregate "
              "after update optimization"),
    cll::init(false));

llvm::cl::opt<bool> disable_self_aggregate(
    "disableSelfAggregation",
    cll::desc("If true (off by default), disables aggregate of self feature"),
    cll::init(false));

llvm::cl::opt<bool>
    do_graph_sampling("doGraphSampling",
                      cll::desc("If true (off by default), sample nodes for "
                                "use every epoch at a 50\% drop rate"),
                      cll::init(false));

llvm::cl::opt<bool> use_train_subgraph(
    "useTrainingSubgraph",
    cll::desc(
        "If true (off by default), during training "
        "only compute minimum required for training nodes in training phase"),
    cll::init(false));

llvm::cl::opt<bool> inductive_subgraph(
    "inductiveSubgraph",
    cll::desc("If true (off by default), only sample training/other nodes when "
              "constructing subgraph"),
    cll::init(false));

llvm::cl::opt<unsigned>
    train_minibatch_size("trainMinibatchSize",
                         cll::desc("Size of training minibatch (default 0)"),
                         cll::init(0));

llvm::cl::opt<unsigned>
    test_minibatch_size("testMinibatchSize",
                        cll::desc("Size of test minibatch (default 0)"),
                        cll::init(0));

llvm::cl::opt<unsigned> minibatch_test_interval(
    "minibatchTestInterval",
    cll::desc("Size of test intervals for minibatch (default 0)"),
    cll::init(0));

llvm::cl::opt<unsigned>
    val_interval("valInterval",
                 cll::desc("# of epochs to test validation set (default 0)"),
                 cll::init(0));

llvm::cl::opt<unsigned>
    test_interval("testInterval",
                  cll::desc("# of epochs to test test set (default 0)"),
                  cll::init(0));

llvm::cl::opt<float>
    learning_rate("learningRate",
                  cll::desc("Adam optimizer learning rate (default 0.01)"),
                  cll::init(0.01));

const char* GNNPartitionToString(galois::graphs::GNNPartitionScheme s) {
  switch (s) {
  case galois::graphs::GNNPartitionScheme::kOEC:
    return "oec";
  case galois::graphs::GNNPartitionScheme::kCVC:
    return "cvc";
  case galois::graphs::GNNPartitionScheme::kOCVC:
    return "ocvc";
  default:
    GALOIS_LOG_FATAL("Invalid partitioning scheme");
    return "";
  }
}

//! Initializes the vector of layer sizes from command line args + graph
std::vector<galois::GNNLayerType> CreateLayerTypesVector() {
  std::vector<galois::GNNLayerType> layer_types;
  for (size_t i = 0; i < num_layers; i++) {
    layer_types.emplace_back(cl_layer_type);
  }
  // if (!cl_layer_types.size()) {
  //  // default is all GCN layers
  //  for (size_t i = 0; i < num_layers; i++) {
  //    layer_types.emplace_back(galois::GNNLayerType::kGraphConvolutional);
  //  }
  //} else {
  //  GALOIS_LOG_VASSERT(cl_layer_types.size() == num_layers,
  //                     "Number layer types should be {} not {}", num_layers,
  //                     cl_layer_types.size());
  //  for (size_t i = 0; i < num_layers; i++) {
  //    layer_types.emplace_back(cl_layer_types[i]);
  //  }
  //}
  return layer_types;
}

//! Initializes the vector of layer sizes from command line args + graph
std::vector<size_t>
CreateLayerSizesVector(const galois::graphs::GNNGraph* gnn_graph) {
  // set layer sizes for intermdiate and output layers
  std::vector<size_t> layer_sizes_vector;

  // if (layer_sizes.size()) {
  //  GALOIS_LOG_ASSERT(layer_sizes.size() == num_layers);
  //  for (size_t i = 0; i < num_layers; i++) {
  //    layer_sizes_vector.emplace_back(layer_sizes[i]);
  //  }
  //  // verify user satisfies last intermediate layer needing to have same size
  //  // as # label classes
  //  if (layer_sizes_vector.back() != gnn_graph->GetNumLabelClasses()) {
  //    galois::gWarn(
  //        "Size of last layer (", layer_sizes_vector.back(),
  //        ") is not equal to # label classes: forcefully changing it to ",
  //        gnn_graph->GetNumLabelClasses());
  //    layer_sizes_vector.back()   = gnn_graph->GetNumLabelClasses();
  //    layer_sizes[num_layers - 1] = gnn_graph->GetNumLabelClasses();
  //  }

  //  GALOIS_LOG_ASSERT(layer_sizes_vector.back() ==
  //                    gnn_graph->GetNumLabelClasses());
  //} else {
  //  // default 16 for everything until last 2
  //  for (size_t i = 0; i < num_layers - 1; i++) {
  //    layer_sizes_vector.emplace_back(16);
  //  }
  //  // last 2 sizes must be equivalent to # label classes; this is the last
  //  // intermediate layer
  //  layer_sizes_vector.emplace_back(gnn_graph->GetNumLabelClasses());
  //}

  for (size_t i = 0; i < num_layers - 1; i++) {
    layer_sizes_vector.emplace_back(layer_size);
  }
  // last 2 sizes must be equivalent to # label classes; this is the last
  // intermediate layer
  layer_sizes_vector.emplace_back(gnn_graph->GetNumLabelClasses());
  // TODO
  // for now only softmax layer which dictates the output size of the last
  // intermediate layer + size of the output layer
  // output layer at the moment required to be same as # label classes
  layer_sizes_vector.emplace_back(gnn_graph->GetNumLabelClasses());

  return layer_sizes_vector;
}

//! Setup layer config struct based on cli args
galois::GNNLayerConfig CreateLayerConfig() {
  galois::GNNLayerConfig layer_config;
  layer_config.disable_dropout                = disable_dropout;
  layer_config.dropout_rate                   = dropout_rate;
  layer_config.disable_activation             = disable_activation;
  layer_config.disable_normalization          = disable_normalization;
  layer_config.disable_aggregate_after_update = disable_agg_after_update;
  layer_config.disable_self_aggregate         = disable_self_aggregate;
  return layer_config;
}

std::unique_ptr<galois::BaseOptimizer>
CreateOptimizer(const galois::graphs::GNNGraph* gnn_graph) {
  std::vector<size_t> opt_sizes;

  // optimizer sizes are based on intermediate layer sizes, input feats, and
  // # label classes
  // if (layer_sizes.size()) {
  //  GALOIS_LOG_ASSERT(layer_sizes.size() == num_layers);
  //  opt_sizes.emplace_back(gnn_graph->node_feature_length() * layer_sizes[0]);
  //  // assumption here is that if it reached this point then layer sizes were
  //  // already sanity checked previously (esp. last layer)
  //  for (size_t i = 1; i < num_layers; i++) {
  //    opt_sizes.emplace_back(layer_sizes[i] * layer_sizes[i - 1]);
  //  }
  //} else {
  //  // everything is size 16 until last
  //  if (num_layers == 1) {
  //    // single layer requires a bit of special handling
  //    opt_sizes.emplace_back(gnn_graph->node_feature_length() *
  //                           gnn_graph->GetNumLabelClasses());
  //  } else {
  //    // first
  //    opt_sizes.emplace_back(gnn_graph->node_feature_length() * 16);
  //    for (size_t i = 1; i < num_layers - 1; i++) {
  //      opt_sizes.emplace_back(16 * 16);
  //    }
  //    // last
  //    opt_sizes.emplace_back(16 * gnn_graph->GetNumLabelClasses());
  //  }
  //}

  // everything is size 16 until last
  if (num_layers == 1) {
    // single layer requires a bit of special handling
    opt_sizes.emplace_back(gnn_graph->node_feature_length() *
                           gnn_graph->GetNumLabelClasses());
  } else {
    // first
    opt_sizes.emplace_back(gnn_graph->node_feature_length() * layer_size);
    for (size_t i = 1; i < num_layers - 1; i++) {
      opt_sizes.emplace_back(layer_size * layer_size);
    }
    // last
    opt_sizes.emplace_back(layer_size * gnn_graph->GetNumLabelClasses());
  }
  GALOIS_LOG_ASSERT(opt_sizes.size() == num_layers);

  galois::AdamOptimizer::AdamConfiguration adam_config;
  adam_config.alpha = learning_rate;

  // TODO only adam works right now, add the others later
  return std::make_unique<galois::AdamOptimizer>(adam_config, opt_sizes,
                                                 num_layers);
}

std::vector<unsigned> CreateFanOutVector() {
  std::vector<unsigned> fan_out;
  // fan out only matters if graph sampling is enabled
  if (do_graph_sampling) {
    // assert fan out size is the same
    if (cl_fan_out_vector.size() == num_layers) {
      for (unsigned i = 0; i < num_layers; i++) {
        fan_out.emplace_back(cl_fan_out_vector[i]);
      }
    } else {
      galois::gWarn("Fan out specification does not equal number of layers: "
                    "using default 10 followed by 25s");
      fan_out.emplace_back(10);
      for (unsigned i = 1; i < num_layers; i++) {
        fan_out.emplace_back(25);
      }
    }
  }
  return fan_out;
}

std::unique_ptr<galois::GraphNeuralNetwork> InitializeGraphNeuralNetwork() {
  // partition/load graph
  auto gnn_graph = std::make_unique<galois::graphs::GNNGraph>(
      input_directory, input_name, partition_scheme, !multiclass_labels,
      useShad);

  // create layer types vector
  std::vector<galois::GNNLayerType> layer_types = CreateLayerTypesVector();
  // sizes
  std::vector<size_t> layer_sizes_vector =
      CreateLayerSizesVector(gnn_graph.get());
  // layer config object
  galois::GNNLayerConfig layer_config = CreateLayerConfig();
  // GNN config object
  galois::GraphNeuralNetworkConfig gnn_config(
      num_layers, layer_types, layer_sizes_vector, output_layer_type,
      do_graph_sampling, layer_config);
  gnn_config.use_train_subgraph_      = use_train_subgraph;
  gnn_config.validation_interval_     = val_interval;
  gnn_config.test_interval_           = test_interval;
  gnn_config.train_minibatch_size_    = train_minibatch_size;
  gnn_config.test_minibatch_size_     = test_minibatch_size;
  gnn_config.minibatch_test_interval_ = minibatch_test_interval;
  gnn_config.inductive_subgraph_      = inductive_subgraph;
  gnn_config.fan_out_vector_          = CreateFanOutVector();

  // optimizer
  std::unique_ptr<galois::BaseOptimizer> opt = CreateOptimizer(gnn_graph.get());

  // create the gnn
  return std::make_unique<galois::GraphNeuralNetwork>(
      std::move(gnn_graph), std::move(opt), std::move(gnn_config));
}
