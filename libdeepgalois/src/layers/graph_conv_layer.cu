#include "deepgalois/layers/graph_conv_layer.h"
#include "deepgalois/math_functions.hh"

namespace deepgalois {

void graph_conv_layer::malloc_and_init() {
  size_t x = input_dims[0];
  size_t y = input_dims[1];
  size_t z = output_dims[1];

  if (dropout_)
    CUDA_CHECK(cudaMalloc((void**)&dropout_mask, x * y * sizeof(mask_t)));
  float_malloc_device(x * y, in_temp);
  init_const_gpu(x * y, 0.0, in_temp);
  if (y <= z) {
    float_malloc_device(x * y, in_temp1);
    init_const_gpu(x * y, 0.0, in_temp1);
  }
  float_malloc_device(x * z, out_temp);
  init_const_gpu(x * z, 0.0, out_temp);
  float_malloc_device(y * z, d_W);
  auto init_range = sqrt(6.0 / (y + z));
  // Glorot & Bengio (AISTATS 2010)
  rng_uniform_gpu(y * z, -init_range, init_range, d_W);
  float_malloc_device(y * z, layer::d_weight_grad);
  init_const_gpu(y * z, 0.0, layer::d_weight_grad);
}

void graph_conv_layer::aggregate(size_t len, GraphGPU& g, const float_t* in,
                                 float_t* out) {
#ifdef USE_CUSPARSE
  deepgalois::update_all_csrmm(len, g, in, out, norm_, norm_consts);
#else
  deepgalois::update_all(len, g, in, out, norm_, norm_consts);
#endif
}

void graph_conv_layer::d_aggregate(size_t len, GraphGPU& g, const float_t* in,
                                   float_t* out) {
#ifdef USE_CUSPARSE
  deepgalois::update_all_csrmm(len, g, in, out, norm_, norm_consts);
#else
  deepgalois::update_all(len, g, in, out, norm_, norm_consts);
#endif
}

void graph_conv_layer::combine(size_t dim_x, size_t dim_y, const float_t* self,
                               const float_t* neighbors, float_t* out) {}

// GPU forward: compute output features
// NOTE: in_data will be used in back-prop, so it can not be modified
void graph_conv_layer::forward_propagation(const float_t* in_data,
                                           float_t* out_data) {
  size_t x = input_dims[0];
  size_t y = input_dims[1];
  size_t z = output_dims[1];

  // currently only support feature length <= 128
  if (z > MAX_NUM_CLASSES) {
    std::cout << "Currently support maximum hidden feature length of "
              << MAX_NUM_CLASSES << "\n";
    exit(0);
  }
  init_const_gpu(x * z, 0.0, out_temp);
  if (dropout_ && phase_ == net_phase::train)
    dropout_gpu(x * y, scale_, dropout_rate_, in_data, dropout_mask, in_temp);
  else
    copy_gpu(x * y, in_data, in_temp);
  if (y > z) {
    sgemm_gpu(CblasNoTrans, CblasNoTrans, x, z, y, 1.0, in_temp, d_W, 0.0,
              out_temp);
    graph_conv_layer::aggregate(z, *graph_gpu, out_temp, out_data);
  } else {
    graph_conv_layer::aggregate(y, *graph_gpu, in_temp, in_temp1);
    sgemm_gpu(CblasNoTrans, CblasNoTrans, x, z, y, 1.0, in_temp1, d_W, 0.0,
              out_data);
  }
  if (act_)
    relu_gpu(x * z, out_data, out_data);
}

// GPU backward: compute input gradients (in_grad) and weight gradients
// (d_weight_grad)
void graph_conv_layer::back_propagation(const float_t* in_data,
                                        const float_t* out_data,
                                        float_t* out_grad, float_t* in_grad) {
  size_t x = input_dims[0];
  size_t y = input_dims[1];
  size_t z = output_dims[1];

  if (act_)
    d_relu_gpu(x * z, out_grad, out_data, out_grad);
  if (y > z) {
    graph_conv_layer::d_aggregate(z, *graph_gpu, out_grad, out_temp);
    if (level_ != 0)
      sgemm_gpu(CblasNoTrans, CblasTrans, x, y, z, 1.0, out_temp, d_W, 0.0,
                in_grad);
    sgemm_gpu(CblasTrans, CblasNoTrans, y, z, x, 1.0, in_data, out_temp, 0.0,
              layer::d_weight_grad);
  } else {
    if (level_ != 0) {
      sgemm_gpu(CblasNoTrans, CblasTrans, x, y, z, 1.0, out_grad, d_W, 0.0,
                in_temp);
      graph_conv_layer::d_aggregate(y, *graph_gpu, in_temp, in_grad);
    }
    sgemm_gpu(CblasTrans, CblasNoTrans, y, z, x, 1.0, in_data, out_grad, 0.0,
              layer::d_weight_grad);
  }
  if (level_ != 0 && dropout_)
    d_dropout_gpu(x * y, scale_, dropout_rate_, in_grad, dropout_mask, in_grad);
}

acc_t graph_conv_layer::get_weight_decay_loss() {
  return l2_norm_gpu(input_dims[1] * output_dims[1], d_W);
}

} // namespace deepgalois
