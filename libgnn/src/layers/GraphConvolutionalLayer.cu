#include "gg.h"
#include "ggcuda.h"
#include "galois/GNNMath.cuh"
#include "galois/layers/GraphConvolutionalLayer.cuh"
#include "galois/cuda/DynamicBitset.h"
#include "sharedptr.h"

// TODO(lhc) better way for this declaration is to declare it
//           inside of the cuda context, but this messed linking to Gluon
extern Shared<DynamicBitset> cuda_bitset_graph_aggregate;

galois::GCNGPUAllocations::~GCNGPUAllocations() {
  GALOIS_LOG_VERBOSE("Freeing GCN layer allocations");
  CUDA_FREE(in_temp_1_);
  CUDA_FREE(in_temp_2_);
  CUDA_FREE(out_temp_);
}

void galois::GCNGPUAllocations::AllocateInTemp1(const size_t size) {
  CUDA_CHECK(cudaMalloc((void**)(&in_temp_1_), size * sizeof(GNNFloat)));
}

void galois::GCNGPUAllocations::AllocateInTemp2(const size_t size) {
  CUDA_CHECK(cudaMalloc((void**)(&in_temp_2_), size * sizeof(GNNFloat)));
}

void galois::GCNGPUAllocations::AllocateOutTemp(const size_t size) {
  CUDA_CHECK(cudaMalloc((void**)(&out_temp_), size * sizeof(GNNFloat)));
}

namespace {
// GPU side aggregation call: no matrix multiply, just regular dst accesses
__global__ void AggregateAllKernel(
    unsigned num_nodes, size_t column_length, const int* edge_index,
    const int* edge_destination, const uint32_t* global_degrees,
    const galois::GNNFloat* node_embeddings, galois::GNNFloat* aggregate_output,
    bool disable_self_aggregate, size_t last_master,
    DynamicBitset* cuda_bitset_graph_aggregate) {
  const unsigned thread_id =
      BLOCK_SIZE * blockIdx.x + threadIdx.x; // global thread index
  const unsigned thread_lane =
      threadIdx.x & (WARP_SIZE - 1); // thread index within the warp
  const unsigned warp_id = thread_id / WARP_SIZE; // global warp index
  const unsigned warp_lane =
      threadIdx.x / WARP_SIZE; // warp index within the CTA
  const unsigned num_warps =
      (BLOCK_SIZE / WARP_SIZE) * gridDim.x; // total number of active warps

  // each warp gets a source: this var holds the first/last edge worked on by
  // that warp
  __shared__ int edge_begin_end[BLOCK_SIZE / WARP_SIZE][2];

  // each warp works on a source: threads in warp split the feature
  for (int src = warp_id; src < static_cast<int>(num_nodes); src += num_warps) {
    galois::GNNFloat src_norm    = 0.0;
    galois::GNNFloat dst_norm    = 0.0;
    galois::GNNFloat norm_to_use = 1.0;

    if (global_degrees != nullptr) {
      src_norm = (global_degrees[src])
                     ? (1.0 / sqrt(static_cast<float>(global_degrees[src] + 1)))
                     : 0.0;
    }

    if (thread_lane < 2) {
      edge_begin_end[warp_lane][thread_lane] = edge_index[src + thread_lane];
    }
    // essentially what this is doing is making 2 of the threads set edge
    // begin/end; all threads wait for sync
    __syncthreads();

    const int row_begin     = edge_begin_end[warp_lane][0];
    const int row_end       = edge_begin_end[warp_lane][1];
    unsigned base_src_index = src * column_length;

    if (!disable_self_aggregate) {
      cuda_bitset_graph_aggregate->set(src);
      if (src < last_master) {
        norm_to_use = src_norm * src_norm;
        for (int i = 0; i < column_length; i += WARP_SIZE) {
          if (thread_lane + i < column_length) {
            aggregate_output[base_src_index + thread_lane + i] =
                node_embeddings[base_src_index + thread_lane + i] * norm_to_use;
          }
        }
      }
    }

    for (int offset = row_begin; offset < row_end; offset++) {
      int dst                 = edge_destination[offset];
      unsigned base_dst_index = dst * column_length;
      cuda_bitset_graph_aggregate->set(src);

      if (global_degrees != nullptr) {
        dst_norm =
            (global_degrees[dst])
                ? (1.0 / sqrt(static_cast<float>(global_degrees[dst] + 1)))
                : 0.0;
        // note that otherwise it's 1.0, so a no-op when it comes to multiply
        norm_to_use = src_norm * dst_norm;
      }

      // NOTE: this is where warp diverges
      // the feature aggregation is split among thread in a warp
      for (int i = 0; i < column_length; i += WARP_SIZE) {
        if ((thread_lane + i) < column_length) {
          if (global_degrees != nullptr) {
            aggregate_output[base_src_index + thread_lane + i] +=
                node_embeddings[base_dst_index + thread_lane + i] * norm_to_use;
          } else {
            aggregate_output[base_src_index + thread_lane + i] +=
                node_embeddings[base_dst_index + thread_lane + i];
          }
        }
      }
    }
  }
}

} // namespace

void galois::GCNGPUAllocations::AggregateAllGPU(
    const graphs::GNNGraphGPUAllocations& gpu_graph, size_t num_nodes,
    size_t column_length, const GNNFloat* node_embeddings,
    GNNFloat* aggregate_output, bool use_norm, bool disable_self_aggregate,
    size_t last_master) {
  // num_nodes should be greater than 0 to avoid negative number of thread
  if (num_nodes == 0) {
    return;
  }

  CUDA_CHECK(cudaMemset(aggregate_output, 0,
                        num_nodes * column_length * sizeof(GNNFloat)));
  if (use_norm) {
    AggregateAllKernel<<<(num_nodes - 1) / WARPS_PER_BLOCK + 1, BLOCK_SIZE>>>(
        num_nodes, column_length, gpu_graph.edge_index(),
        gpu_graph.edge_destinations(), gpu_graph.get_global_degrees(),
        node_embeddings, aggregate_output, disable_self_aggregate, last_master,
        cuda_bitset_graph_aggregate.gpu_wr_ptr());
  } else {
    AggregateAllKernel<<<(num_nodes - 1) / WARPS_PER_BLOCK + 1, BLOCK_SIZE>>>(
        num_nodes, column_length, gpu_graph.edge_index(),
        gpu_graph.edge_destinations(), nullptr, node_embeddings,
        aggregate_output, disable_self_aggregate, last_master,
        cuda_bitset_graph_aggregate.gpu_wr_ptr());
  }
  CUDA_TEST("GPU aggregate all failure");
}

void galois::GCNGPUAllocations::UpdateEmbeddingsGPU(
    size_t num_nodes, size_t input_columns, size_t output_columns,
    const GNNFloat* node_embeddings, const GNNFloat* layer_weights,
    GNNFloat* output) {
  CBlasSGEMMGPU(CUBLAS_OP_N, CUBLAS_OP_N, num_nodes, input_columns,
                output_columns, node_embeddings, layer_weights, output);
}

void galois::GCNGPUAllocations::UpdateEmbeddingsDerivativeGPU(
    size_t num_nodes, size_t input_columns, size_t output_columns,
    const GNNFloat* gradients, const GNNFloat* layer_weights,
    GNNFloat* output) {
  // note output clumns/input columns are flipped due to transpose of the
  // layer weights
  CBlasSGEMMGPU(CUBLAS_OP_N, CUBLAS_OP_T, num_nodes, output_columns,
                input_columns, gradients, layer_weights, output);
}

void galois::GCNGPUAllocations::GetWeightGradientsGPU(
    size_t num_nodes, size_t input_columns, size_t output_columns,
    const GNNFloat* prev_input, const GNNFloat* gradients, GNNFloat* output) {
  CBlasSGEMMGPU(CUBLAS_OP_T, CUBLAS_OP_N, input_columns, num_nodes,
                output_columns, prev_input, gradients, output);
}
