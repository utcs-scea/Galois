#pragma once
#include "layer.h"

namespace deepgalois {
// ReLU Layer
class relu_layer : public layer {
public:
  relu_layer(unsigned level, dims_t in_dims, dims_t out_dims)
      : layer(level, in_dims, out_dims) {
    trainable_ = false;
  }
  ~relu_layer() {}
  std::string layer_type() const override { return std::string("relu"); }
  virtual void forward_propagation(const float_t* in_data, float_t* out_data);
  virtual void back_propagation(const float_t* in_data, const float_t* out_data,
                                float_t* out_grad, float_t* in_grad);
};
} // namespace deepgalois
