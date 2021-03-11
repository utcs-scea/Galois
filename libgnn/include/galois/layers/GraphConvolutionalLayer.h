#pragma once
#include "galois/layers/GNNLayer.h"

#ifdef GALOIS_ENABLE_GPU
#include "galois/layers/GraphConvolutionalLayer.cuh"
#endif

namespace galois {

extern galois::DynamicBitSet graphs::bitset_graph_aggregate;

class GraphConvolutionalLayer : public GNNLayer {
public:
  //! Initializes the variables of the base class and also allocates additional
  //! memory for temporary matrices. Also initializes sync substrate for the
  //! weight matrix
  GraphConvolutionalLayer(size_t layer_num,
                          const galois::graphs::GNNGraph& graph,
                          PointerWithSize<GNNFloat>* backward_output_matrix,
                          const GNNLayerDimensions& dimensions,
                          const GNNLayerConfig& config);

  GraphConvolutionalLayer(size_t layer_num,
                          const galois::graphs::GNNGraph& graph,
                          PointerWithSize<GNNFloat>* backward_output_matrix,
                          const GNNLayerDimensions& dimensions)
      : GraphConvolutionalLayer(layer_num, graph, backward_output_matrix,
                                dimensions, GNNLayerConfig()) {}

  // Parent functions
  const PointerWithSize<galois::GNNFloat>
  ForwardPhase(const PointerWithSize<galois::GNNFloat> input_embeddings) final;

  PointerWithSize<galois::GNNFloat>
  BackwardPhase(PointerWithSize<galois::GNNFloat> prev_layer_input,
                PointerWithSize<galois::GNNFloat>* input_gradient) final;

private:
  static const constexpr char* kRegionName = "GCNLayer";
  // 2 temporaries the size of the forward input; used for dropout and
  // aggregation (if either are required)
  std::vector<GNNFloat> in_temp_1_;
  std::vector<GNNFloat> in_temp_2_;
  // Temporary matrix the size of the output of the forward pass; used if
  // an intermediate op occurs before writing to the final output matrix
  std::vector<GNNFloat> out_temp_;

  // Pointer with size versions
  PointerWithSize<GNNFloat> p_in_temp_1_;
  PointerWithSize<GNNFloat> p_in_temp_2_;
  PointerWithSize<GNNFloat> p_out_temp_;

  // Each thread has a vector of size # input columns or # output columns for
  // storing intermediate results during aggregation.
  // The one used depeneds on if aggregation occurs before or after the mxm.
  galois::substrate::PerThreadStorage<std::vector<GNNFloat>>
      input_column_intermediates_;
  galois::substrate::PerThreadStorage<std::vector<GNNFloat>>
      output_column_intermediates_;

  //! CPU aggregation
  void AggregateAllCPU(
      size_t column_length, const GNNFloat* node_embeddings,
      GNNFloat* aggregate_output,
      galois::substrate::PerThreadStorage<std::vector<GNNFloat>>* pts);

  //! Performs aggregation for all nodes of the graph given the length of the
  //! vector to aggregate, the features themselves, an output array, and per
  //! thread storage for the intermediate scaling via norm factor
  void
  AggregateAll(size_t column_length, const GNNFloat* node_embeddings,
               GNNFloat* aggregate_output,
               galois::substrate::PerThreadStorage<std::vector<GNNFloat>>* pts);
  void
  AggregateAll(size_t column_length, const GNNFloat* node_embeddings,
               GNNFloat* aggregate_output,
               galois::substrate::PerThreadStorage<std::vector<GNNFloat>>* pts,
               bool is_backward);

  //! Do embedding update via mxm with this layer's weights (forward)
  void UpdateEmbeddings(const GNNFloat* node_embeddings, GNNFloat* output);
  //! Calculate graident via mxm with last layer's gradients (backward)
  void UpdateEmbeddingsDerivative(const GNNFloat* gradients, GNNFloat* output);
#ifdef GALOIS_ENABLE_GPU
  GCNGPUAllocations gpu_object_;
#endif
};

} // namespace galois
