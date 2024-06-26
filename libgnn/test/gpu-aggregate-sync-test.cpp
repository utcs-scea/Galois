//! @file gpu-aggregate-sync-test.cpp
//! GPU sync test to make sure it's sane
#include "galois/Logging.h"
#include "galois/GraphNeuralNetwork.h"
#include "galois/layers/GraphConvolutionalLayer.h"
#include "galois/CUDAUtilHostDecls.h"

int main() {
  galois::DistMemSys G;

  if (galois::runtime::getSystemNetworkInterface().Num == 1) {
    GALOIS_LOG_ERROR("This test should be run with multiple hosts/processes");
    exit(1);
  }
  device_personality = DevicePersonality::GPU_CUDA;
  gpudevice          = galois::runtime::getSystemNetworkInterface().ID;
  SetCUDADeviceId(gpudevice);

  auto test_graph = std::make_unique<galois::graphs::GNNGraph<char, void>>(
      "tester", galois::graphs::GNNPartitionScheme::kOEC, true, false);

  // create same layer from convlayer-test and make sure result is the same even
  // in multi-host environment
  galois::GNNLayerDimensions dimension_0;
  dimension_0.input_rows     = test_graph->size();
  dimension_0.input_columns  = 3;
  dimension_0.output_columns = 2;
  galois::GNNLayerConfig l_config;
  l_config.disable_aggregate_after_update = true;

  unsigned num_layers = 2;
  test_graph->ResizeGPULayerVector(num_layers);
  test_graph->InitLayerVectorMetaObjects(
      0, galois::runtime::getSystemNetworkInterface().Num,
      dimension_0.input_columns, dimension_0.output_columns);
  test_graph->InitLayerVectorMetaObjects(
      1, galois::runtime::getSystemNetworkInterface().Num,
      dimension_0.input_columns, dimension_0.output_columns);

  galois::PointerWithSize<galois::GNNFloat> p_null(nullptr, 0);
  std::vector<galois::GNNFloat> back_matrix(21);
  galois::PointerWithSize<galois::GNNFloat> p_back(back_matrix);

  // create the layer, no norm factor
  std::unique_ptr<galois::GraphConvolutionalLayer<char, void>> layer_0 =
      std::make_unique<galois::GraphConvolutionalLayer<char, void>>(
          0, *(test_graph.get()), &p_null, dimension_0, l_config);
  layer_0->InitAllWeightsTo1();
  // make sure it runs in a sane manner
  layer_0->ForwardPhase(test_graph->GetLocalFeatures());
  // pointer is to GPU memory: copy it over to a CPU source for verification
  const std::vector<galois::GNNFloat>& layer_0_forward_output =
      layer_0->CopyForwardOutputFromGPU();

  //////////////////////////////////////////////////////////////////////////////
  // sanity check output
  //////////////////////////////////////////////////////////////////////////////

  // check each row on each host: convert row into GID, and based on GID we
  // know what the ground truth is
  // row 0 = 3
  // row 1 = 6
  // row 2 = 12
  // row 3 = 18
  // row 4 = 24
  // row 5 = 30
  // row 6 = 15

  // row should correspond to LID
  for (size_t row = 0; row < test_graph->size(); row++) {
    // row -> GID
    size_t global_row = test_graph->GetGID(row);

    galois::GNNFloat ground_truth = 0.0;

    switch (global_row) {
    case 0:
      ground_truth = 3;
      break;
    case 1:
      ground_truth = 6;
      break;
    case 2:
      ground_truth = 12;
      break;
    case 3:
      ground_truth = 18;
      break;
    case 4:
      ground_truth = 24;
      break;
    case 5:
      ground_truth = 30;
      break;
    case 6:
      ground_truth = 15;
      break;
    default:
      GALOIS_LOG_FATAL("bad global row for test graph");
      break;
    }

    // size 2 columns
    for (size_t c = 0; c < 2; c++) {
      GALOIS_LOG_ASSERT(layer_0_forward_output[row * 2 + c] == ground_truth);
    }
  }

  //////////////////////////////////////////////////////////////////////////////

  std::vector<galois::GNNFloat> dummy_ones_v(test_graph->size() * 2, 1);
  galois::PointerWithSize<galois::GNNFloat> dummy_ones =
      layer_0->AllocateGPU(dummy_ones_v);
  // backward pass checking
  // layer 0 means that an empty weight matrix is returned since there is no
  // point passing back anything
  layer_0->BackwardPhase(test_graph->GetLocalFeatures(), &dummy_ones);
  const galois::PointerWithSize<galois::GNNFloat>& layer_0_backward_output =
      layer_0->CopyBackwardOutputFromGPU();

  //////////////////////////////////////////////////////////////////////////////
  // sanity check layer 0 backward output; all 0 because layer 0
  //////////////////////////////////////////////////////////////////////////////
  // since norm factors aren't invovled it is possible to do full assertions
  GALOIS_LOG_ASSERT(layer_0_backward_output.size() == test_graph->size() * 3);
  for (size_t i = 0; i < layer_0_backward_output.size(); i++) {
    GALOIS_LOG_ASSERT((layer_0_backward_output)[i] == 0);
  }

  //////////////////////////////////////////////////////////////////////////////
  // layer 1 to check backward output
  //////////////////////////////////////////////////////////////////////////////
  std::unique_ptr<galois::GraphConvolutionalLayer<char, void>> layer_1 =
      std::make_unique<galois::GraphConvolutionalLayer<char, void>>(
          1, *(test_graph.get()), &p_back, dimension_0, l_config);
  layer_1->InitAllWeightsTo1();
  layer_1->ForwardPhase(test_graph->GetLocalFeatures());
  const std::vector<galois::GNNFloat>& layer_1_forward_output =
      layer_1->CopyForwardOutputFromGPU();

  // same check for forward as before
  for (size_t row = 0; row < test_graph->size(); row++) {
    // row -> GID
    size_t global_row = test_graph->GetGID(row);

    galois::GNNFloat ground_truth = 0.0;

    switch (global_row) {
    case 0:
      ground_truth = 3;
      break;
    case 1:
      ground_truth = 6;
      break;
    case 2:
      ground_truth = 12;
      break;
    case 3:
      ground_truth = 18;
      break;
    case 4:
      ground_truth = 24;
      break;
    case 5:
      ground_truth = 30;
      break;
    case 6:
      ground_truth = 15;
      break;
    default:
      GALOIS_LOG_FATAL("bad global row for test graph");
      break;
    }

    // size 2 columns
    for (size_t c = 0; c < 2; c++) {
      GALOIS_LOG_ASSERT(layer_1_forward_output[row * 2 + c] == ground_truth);
    }
  }

  // since layer isn't 0 anymore, backward phase will actually return something
  dummy_ones_v.assign(test_graph->size() * 2, 1);
  layer_1->BackwardPhase(test_graph->GetLocalFeatures(), &dummy_ones);
  const galois::PointerWithSize<galois::GNNFloat>& layer_1_backward_output =
      layer_1->CopyBackwardOutputFromGPU();

  for (size_t row = 0; row < test_graph->size(); row++) {
    // row -> GID
    size_t global_row = test_graph->GetGID(row);

    galois::GNNFloat ground_truth = 0.0;

    switch (global_row) {
    case 0:
    case 6:
      ground_truth = 2;
      break;
    case 1:
    case 2:
    case 3:
    case 4:
    case 5:
      ground_truth = 4;
      break;
    default:
      GALOIS_LOG_FATAL("bad global row for test graph");
      break;
    }

    // size 3 columns
    for (size_t c = 0; c < 3; c++) {
      GALOIS_LOG_ASSERT((layer_1_backward_output)[row * 3 + c] == ground_truth);
    }
  }

  // TODO CVC
}
