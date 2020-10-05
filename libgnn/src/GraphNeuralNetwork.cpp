#include "galois/GraphNeuralNetwork.h"
#include "galois/layers/GraphConvolutionalLayer.h"
#include "galois/layers/SoftmaxLayer.h"

galois::GraphNeuralNetwork::GraphNeuralNetwork(
    std::unique_ptr<galois::graphs::GNNGraph> graph,
    galois::GraphNeuralNetworkConfig&& config)
    : graph_(std::move(graph)), config_(std::move(config)) {
  // max number of rows that can be passed as inputs; allocate space for it as
  // this will be the # of rows for each layer
  size_t max_rows = graph_->size();

  // create the intermediate layers
  for (size_t i = 0; i < config_.num_intermediate_layers(); i++) {
    GNNLayerType layer_type = config_.intermediate_layer_type(i);
    size_t prev_layer_columns;

    if (i != 0) {
      // grab previous layer's size
      prev_layer_columns = config_.intermediate_layer_size(i - 1);
    } else {
      // first layer means the input columns are # features in graph
      prev_layer_columns = graph_->node_feature_length();
    }

    GNNLayerDimensions layer_dims = {.input_rows    = max_rows,
                                     .input_columns = prev_layer_columns,
                                     .output_columns =
                                         config_.intermediate_layer_size(i)};

    switch (layer_type) {
    case GNNLayerType::kGraphConvolutional:
      gnn_layers_.push_back(std::move(std::make_unique<GraphConvolutionalLayer>(
          i, *graph_, layer_dims, config_.default_layer_config())));
      break;
    default:
      GALOIS_LOG_FATAL("Invalid layer type during network construction");
    }
  }

  // create the output layer
  GNNLayerDimensions output_dims = {
      .input_rows = max_rows,
      // get last intermediate layer column size
      .input_columns = config_.intermediate_layer_size(
          config_.num_intermediate_layers() - 1),
      .output_columns = config_.output_layer_size()};

  switch (config_.output_layer_type()) {
  case (GNNOutputLayerType::kSoftmax):
    gnn_layers_.push_back(std::move(std::make_unique<SoftmaxLayer>(
        config_.num_intermediate_layers(), *graph_, output_dims)));
    break;
  default:
    GALOIS_LOG_FATAL("Invalid layer type during network construction");
  }
}

const std::vector<galois::GNNFloat>* galois::GraphNeuralNetwork::DoInference() {
  // start with graph features and pass it through all layers of the network
  const std::vector<GNNFloat>* layer_input = &(graph_->GetLocalFeatures());
  for (std::unique_ptr<galois::GNNLayer>& ptr : gnn_layers_) {
    layer_input = &(ptr->ForwardPhase(*layer_input));
  }
  return layer_input;
}

void galois::GraphNeuralNetwork::GradientPropagation() {
  // from output layer get initial gradients
  std::vector<galois::GNNFloat> dummy;
  std::unique_ptr<galois::GNNLayer>& output_layer = gnn_layers_.back();
  std::vector<galois::GNNFloat>* current_gradients =
      output_layer->BackwardPhase(dummy, nullptr);

  // loops through intermediate layers in a backward fashion
  // -1 to ignore output layer which was handled above
  for (size_t i = 0; i < gnn_layers_.size() - 1; i++) {
    // note this assumes you have at least 2 layers
    size_t layer_index = gnn_layers_.size() - 2 - i;

    // get the input to the layer before this one
    const std::vector<galois::GNNFloat>* prev_layer_input;
    if (layer_index != 0) {
      prev_layer_input = &(gnn_layers_[layer_index - 1]->GetForwardOutput());
    } else {
      prev_layer_input = &(graph_->GetLocalFeatures());
    }

    // backward prop and get a new set of gradients
    current_gradients = gnn_layers_[layer_index]->BackwardPhase(
        *prev_layer_input, current_gradients);
    // at this point in the layer the gradients exist; use the gradients to
    // update the weights of the layer
    // XXX need optimizers
  }
}