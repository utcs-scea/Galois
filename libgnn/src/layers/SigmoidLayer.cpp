#include "galois/layers/SigmoidLayer.h"
#include "galois/GNNMath.h"
#include <math.h>

// TODO(loc) GPU support

const galois::PointerWithSize<galois::GNNFloat>
galois::SigmoidLayer::ForwardPhaseCPU(
    const galois::PointerWithSize<galois::GNNFloat> input_embeddings) {
  // loss is ignored for now anyways
  // input_loss_.assign(input_loss_.size(), 0.0);
  forward_output_matrix_.assign(forward_output_matrix_.size(), 0.0);
  const size_t feature_length = layer_dimensions_.input_columns;

  galois::do_all(
      galois::iterate(graph_.begin_owned(), graph_.end_owned()),
      [&](const unsigned local_node) {
        if (graph_.IsValidForPhase(local_node, layer_phase_)) {
          size_t node_offset = feature_length * local_node;
          // sigmoid the values for this node
          for (unsigned index = 0; index < feature_length; index++) {
            forward_output_matrix_[node_offset + index] =
                1.0 / (1.0 + expf(-input_embeddings[node_offset + index]));
          }
          // TODO(loc) calculate loss (it's not even being used/not required
          // for correctness so I'm ignoring it for now)
        }
      },
      galois::steal(), galois::loopname("SigmoidForward"));

  return forward_output_matrix_;
}

const galois::PointerWithSize<galois::GNNFloat>
galois::SigmoidLayer::ForwardPhase(
    const galois::PointerWithSize<galois::GNNFloat> input_embeddings) {
#ifdef GALOIS_ENABLE_GPU
  // TODO(loc) when GPU needs it
  return 0;
#else
  return ForwardPhaseCPU(input_embeddings);
#endif
}

galois::PointerWithSize<galois::GNNFloat>
galois::SigmoidLayer::BackwardPhaseCPU() {
  const size_t feature_length = layer_dimensions_.input_columns;
  backward_output_matrix_.assign(backward_output_matrix_.size(), 0);

  galois::do_all(
      galois::iterate(graph_.begin_owned(), graph_.end_owned()),
      [&](const unsigned local_node) {
        if (graph_.IsValidForPhase(local_node, layer_phase_)) {
          // derivative cross entropy into norm grad
          const GNNLabel* ground_truth = graph_.GetMultiClassLabel(local_node);
          size_t node_offset           = feature_length * local_node;
          std::vector<GNNFloat>* norm_gradient =
              norm_gradient_vectors_.getLocal();
          GNNCrossEntropyDerivative(feature_length, ground_truth,
                                    &(forward_output_matrix_[node_offset]),
                                    norm_gradient->data());

          // sigmoid derivative
          for (unsigned index = 0; index < feature_length; index++) {
            backward_output_matrix_[node_offset + index] =
                (*norm_gradient)[index] *
                forward_output_matrix_[node_offset + index] *
                (1.0 - forward_output_matrix_[node_offset + index]);
          }
        }
      },
      galois::steal(), galois::loopname("SigmoidBackward"));

  return backward_output_matrix_;
}

galois::PointerWithSize<galois::GNNFloat>
galois::SigmoidLayer::BackwardPhase(const PointerWithSize<galois::GNNFloat>,
                                    PointerWithSize<galois::GNNFloat>*) {
#ifdef GALOIS_ENABLE_GPU
  // TODO(loc) when GPU needs it
  return 0;
#else
  return BackwardPhaseCPU();
#endif
}
