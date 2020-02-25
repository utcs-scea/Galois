#include "layers/softmax_loss_layer.h"

softmax_loss_layer::softmax_loss_layer(unsigned level,
                                       std::vector<size_t> in_dims,
                                       std::vector<size_t> out_dims)
    : layer(level, in_dims, out_dims) {
  trainable_ = false;
  name_      = layer_type() + "_" + std::to_string(level);
  std::cout << name_ << ": allocating memory for intermediate data... ";
#ifdef CPU_ONLY
  loss = new float_t[in_dims[0]]; // error for each sample
#else
  loss_malloc_device(in_dims[0], loss);
#endif
  std::cout << "Done\n";
}
#ifdef CPU_ONLY
// TODO: need kernel fusion optimization
// 𝑦[i] = 𝑒^𝑥[i] / Σ 𝑒^𝑥[𝑘]
void softmax_loss_layer::forward_propagation(const float_t* in_data,
                                             float_t* out_data) {
  size_t len = input_dims[1];
  galois::do_all(galois::iterate(begin_, end_),
                 [&](const auto& i) {
                   if (masks_[i] == 1) { // masked
                     softmax(len, &in_data[len * i],
                             &out_data[len * i]); // normalize using softmax
                     // y is a one hot encoded vector for the labels
                     std::vector<acc_t> y(output_dims[1], 0.0); // ground truth
                     y[context->get_label(i)] = 1.0;            // one-hot
                     loss[i] = cross_entropy(len, &y[0], &out_data[len * i]);
                   }
                 },
                 galois::chunk_size<CHUNK_SIZE>(), galois::steal(),
                 galois::loopname("softmax-loss-fw"));
}

void softmax_loss_layer::back_propagation(const float_t* in_data,
                                          const float_t* out_data,
                                          float_t* out_grad, float_t* in_grad) {
  size_t len = input_dims[1];
  galois::do_all(galois::iterate(begin_, end_),
                 [&](const auto& i) {
                   if (masks_[i] == 1) { // masked
                     vec_t norm_grad(len);
                     std::vector<acc_t> y(len, 0.0); // ground truth
                     y[context->get_label(i)] = 1.0;
                     d_cross_entropy(len, &y[0], &out_data[len * i],
                                     &norm_grad[0]);
                     d_softmax(len, &in_data[len * i], &out_data[len * i],
                               &in_grad[len * i], &norm_grad[0]);
                   }
                 },
                 galois::chunk_size<CHUNK_SIZE>(), galois::steal(),
                 galois::loopname("softmax-loss-bw"));
}

acc_t softmax_loss_layer::get_masked_loss() {
  AccumF total_loss;
  AccumU valid_sample_count;
  total_loss.reset();
  valid_sample_count.reset();
  galois::do_all(galois::iterate(begin_, end_),
                 [&](const auto& i) {
                   if (masks_[i]) {
                     total_loss += loss[i];
                     valid_sample_count += 1;
                   }
                 },
                 galois::chunk_size<256>(), galois::steal(),
                 galois::loopname("getMaskedLoss"));
  assert(valid_sample_count.reduce() == count_);
  return total_loss.reduce() / (acc_t)count_;
}

#else // GPU implementation
void softmax_loss_layer::forward_propagation(const float_t* in_data,
                                             float_t* out_data) {
  softmax_cross_entropy_gpu(input_dims[0], input_dims[1], in_data, d_masks_,
                            context->d_labels, loss, out_data);
}

void softmax_loss_layer::back_propagation(const float_t* in_data,
                                          const float_t* out_data,
                                          float_t* out_grad, float_t* in_grad) {
  d_softmax_cross_entropy_gpu(input_dims[0], input_dims[1], in_data, d_masks_,
                              context->d_labels, out_data, in_grad);
}

acc_t softmax_loss_layer::get_masked_loss() {
  return masked_avg_loss(begin_, end_, count_, masks_, loss);
}
#endif
