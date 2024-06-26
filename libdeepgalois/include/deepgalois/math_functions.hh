/**
 * File inspired by similar one from TinyDNN
 * https://github.com/tiny-dnn/
 */
#ifndef _MATH_FUNCTIONS_
#define _MATH_FUNCTIONS_
#include <cmath>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include "deepgalois/types.h"

#ifdef USE_MKL
#include <mkl.h>
#else
extern "C" {
#include <cblas.h>
}
#endif

namespace deepgalois {

namespace math {

// single-precision dense matrix multiply
void sgemm_cpu(const CBLAS_TRANSPOSE TransA, const CBLAS_TRANSPOSE TransB,
               const int M, const int N, const int K, const float alpha,
               const float* A, const float* B, const float beta, float* C);

// single-precision sparse matrix dense matrix multiply, C = A * B, A is sparse
void csrmm_cpu(const int M, const int N, const int K, const int nnz,
               const float alpha, float* A_nonzeros, int* A_idx_ptr,
               int* A_nonzero_idx, const float* B, const float beta, float* C);

// matrix-vector multiply
void mvmul(const CBLAS_TRANSPOSE TransA, const int M, const int N,
           const float alpha, const float* A, const float* x, const float beta,
           float* y);

//! add 2 arrays for n elements
void vadd_cpu(size_t n, const float_t* a, const float_t* b, float_t* out);

//! multiply n elements of vector by scalar
void scal(size_t n, const float_t alpha, float_t* x);
void scale(size_t n, const float_t alpha, const float_t* x, float_t* y);
void mul_scalar(size_t n, const float_t alpha, const float_t* x, float_t* y);

//! do dot product of 2 vectors
float_t dot(size_t n, const float_t* x, const float_t* y);

// concatenation of two vectors into one
void concat(size_t n, const float_t* x, const float_t* y, float_t* z);

// SAXPY stands for “Single-precision A*X Plus Y"
void axpy(size_t n, const float_t a, float_t* x, float_t* y);

// Returns the index of the maximum value
int argmax(const size_t n, const float_t* x); // the arguments of the maxima

//! Computes half the L2 norm of a tensor without the sqrt: output = sum(t ** 2)
//! / 2
float_t l2_norm(size_t n, const float_t* a);

//! clear n elements of a vector
void clear_cpu(size_t n, float_t* in);

//! copy vector from in -> out; first len elements
void copy_cpu(size_t len, const float_t* in, float_t* out);

// dropout functions randomly remove weights
void dropout_cpu(size_t n, size_t m, float scale, float dropout_rate,
                 const float_t* in, mask_t* mask, float_t* out);

// dropout derivative: use existing dropouts in masks instead of generating
// them;
void d_dropout_cpu(size_t n, size_t m, float scale, const float_t* in,
                   mask_t* mask, float_t* out);

//! ReLU = keep if positive; and ReLU derivative: 1 if data > 0, 0 otherwise
void relu_cpu(size_t n, const float_t* in, float_t* out);
void d_relu_cpu(size_t n, const float_t* in, const float_t* data, float_t* out);

// Leaky ReLU
void leaky_relu(float_t epsilon, float_t in, float_t& out);
void d_leaky_relu(float_t epsilon, float_t in, float_t data, float_t& out);
void leaky_relu_cpu(size_t n, float_t epsilon, const float_t* in, float_t* out);
void d_leaky_relu_cpu(size_t n, float_t epsilon, const float_t* in,
                      const float_t* data, float_t* out);

// Loss function for single-class label (one-hot) data: softmax
void softmax(size_t n, const float_t* input, float_t* output);
void d_softmax(size_t n, const float_t* y, const float_t* p, float_t* dy,
               const float_t* dp);

// Cross entropy
float_t cross_entropy(size_t n, const float_t* y, const float_t* p);
void d_cross_entropy(size_t n, const float_t* y, const float_t* p, float_t* d);

// Loss function for multi-class label (one-hot) data: sigmoid
void sigmoid(size_t n, const float_t* input, float_t* output);
void d_sigmoid(size_t n, const float_t* y, const float_t* p, float_t* dy,
               const float_t* dp);

// dropout functions randomly remove weights
void dropout(float scale, float dropout_rate, const float_t* in, mask_t* mask,
             float_t* out);
void d_dropout(const float scale, const float_t* in, mask_t* mask,
               float_t* out);

//! transposes a matrix (malloc'd array)
void transpose(size_t x, size_t y, const float_t* in, float_t* out);

} // namespace math
} // namespace deepgalois

// GPU operators
bool isnan_gpu(int n,
               const float_t* array); // does array contain any 'nan' element
void init_const_gpu(int n, float_t value, float_t* array);
void copy_gpu(int len, const float_t* in, float_t* out);
void vadd_gpu(const int n, const float_t* a, const float_t* b,
              float_t* out); // vector add
void axpy_gpu(const int n, const float_t a, const float_t* x,
              float_t* y);                                   // axpy
void relu_gpu(const int n, const float_t* in, float_t* out); // ReLU
void d_relu_gpu(const int n, const float_t* in_diff, const float_t* data,
                float_t* out_diff); // ReLU derivative
void leaky_relu_gpu(const int n, const float_t epsilon, const float_t* in,
                    float_t* out); // Leaky ReLU
void d_leaky_relu_gpu(const int n, const float_t epsilon,
                      const float_t* in_diff, const float_t* data,
                      float_t* out_diff); // Leaky ReLU derivative
void dropout_gpu(int n, float scale, float dropout_rate, const float_t* in,
                 mask_t* masks, float_t* out); // dropout
void d_dropout_gpu(int n, float scale, float dropout_rate, const float_t* in,
                   const mask_t* masks, float_t* out); // dropout derivative
void sgemm_gpu(const CBLAS_TRANSPOSE TransA, const CBLAS_TRANSPOSE TransB,
               const int M, const int N, const int K, const float alpha,
               const float* A, const float* B, const float beta, float* C);
void matmul_gpu(const size_t x, const size_t y, const size_t z,
                const float_t* A, const float_t* B, float_t* C);
void matmul1D1D_gpu(const size_t dim_x, const size_t dim_y, const size_t dim_z,
                    const float_t* A, const float_t* B,
                    float_t* C); // matrix multiply
void csrmm_gpu(const int M, const int N, const int K, const int nnz,
               const float alpha, const float* A_nonzeros, const int* A_idx_ptr,
               const int* A_nonzero_idx, const float* B, const float beta,
               float* trans_C, float* C);
void softmax_cross_entropy_gpu(int len, int begin, int end,
                               const float_t* in_data, const mask_t* masks,
                               const label_t* labels, float_t* loss,
                               float_t* out_data);
void d_softmax_cross_entropy_gpu(int len, int bengin, int end,
                                 const mask_t* masks, const label_t* labels,
                                 const float_t* out_data, float_t* diff);
void sigmoid_cross_entropy_gpu(int len, int begin, int end,
                               const float_t* in_data, const mask_t* masks,
                               const label_t* labels, float_t* loss,
                               float_t* out_data);
void d_sigmoid_cross_entropy_gpu(int len, int bengin, int end,
                                 const mask_t* masks, const label_t* labels,
                                 const float_t* out_data, float_t* diff);
void scal_gpu(const int n, const float alpha, float* X);
void add_scalar_gpu(const int n, const float_t alpha, float_t* Y);
void rng_uniform_gpu(size_t n, const float_t a, const float_t b, float_t* r);
bool is_allocated_device(float_t* data);
void copy_masks_device(int n, mask_t* h_masks, mask_t*& d_masks);
void float_malloc_device(int n, float_t*& ptr);
void float_free_device(float_t*& ptr);
void float_copy_device(int n, float_t* h_ptr, float_t* d_ptr);
void uint8_malloc_device(int n, uint8_t*& ptr);
void uint8_free_device(uint8_t*& ptr);
void uint8_copy_device(int n, uint8_t* h_ptr, uint8_t* d_ptr);
acc_t masked_avg_loss_gpu(int begin, int end, int count, mask_t* masks,
                          float_t* loss);
acc_t l2_norm_gpu(int n, const float_t* in);
void l2_norm_gpu(size_t x, size_t y, const float_t* in, float_t* out);
void d_l2_norm_gpu(size_t x, size_t y, const float_t* in_data, float_t* in_diff,
                   float_t* out_diff);
#endif
