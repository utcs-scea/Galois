#include "deepgalois/layers/sigmoid_loss_layer.h"
#include "deepgalois/math_functions.hh"
#include "gg.h"
#include "ggcuda.h"

namespace deepgalois {

sigmoid_loss_layer::sigmoid_loss_layer(unsigned level,
                                       std::vector<size_t> in_dims,
                                       std::vector<size_t> out_dims)
    : layer(level, in_dims, out_dims) {
  trainable_ = false;
  name_      = layer_type() + "_" + std::to_string(level);
}

sigmoid_loss_layer::~sigmoid_loss_layer() { float_free_device(loss); }

void sigmoid_loss_layer::malloc_and_init() {
  float_malloc_device(input_dims[0], loss);
}

void sigmoid_loss_layer::forward_propagation(const float_t* in_data,
                                             float_t* out_data) {
  init_const_gpu(input_dims[0], 0.0, loss);
  sigmoid_cross_entropy_gpu(input_dims[1], begin_, end_, in_data, d_masks_,
                            labels, loss, out_data);
}

void sigmoid_loss_layer::back_propagation(const float_t* in_data,
                                          const float_t* out_data,
                                          float_t* out_grad, float_t* in_grad) {
  d_sigmoid_cross_entropy_gpu(input_dims[1], begin_, end_, d_masks_, labels,
                              out_data, in_grad);
}

acc_t sigmoid_loss_layer::get_prediction_loss() {
  return masked_avg_loss_gpu(begin_, end_, count_, d_masks_, loss);
}

} // namespace deepgalois
