#include "math_functions.hh"
#include "context.h"

extern "C" {
#include <cblas.h>
//#include <clapack.h>
}

// flattern data into 1D before feed into the ReLU operater
__global__ void relu_kernel(const int n, const float_t* in, float_t* out) {
	CUDA_KERNEL_LOOP(index, n) {
		out[index] = in[index] > 0 ? in[index] : 0;
	}
}

void relu_gpu(const int n, const float_t *in_data, float_t* out_data) {
	relu_kernel<<<CUDA_GET_BLOCKS(n), CUDA_NUM_THREADS>>>(n, in_data, out_data);
}

__global__ void d_relu_kernel(const int n, const float_t* in_diff, const float_t* data, float_t* out_diff) {
	CUDA_KERNEL_LOOP(index, n) {
		out_diff[index] = data[index] > 0 ? in_diff[index] : 0;
	}
}

void d_relu_gpu(const int n, const float_t *in_diff, const float_t *data, float_t *out_diff) {
	d_relu_kernel<<<CUDA_GET_BLOCKS(n), CUDA_NUM_THREADS>>>(n, in_diff, data, out_diff);
}

void sgemm_gpu(const CBLAS_TRANSPOSE TransA, const CBLAS_TRANSPOSE TransB, 
	const int M, const int N, const int K, const float alpha, 
	const float* A, const float* B, const float beta, float* C) {
	// Note that cublas follows fortran order.
	int lda = (TransA == CblasNoTrans) ? K : M;
	int ldb = (TransB == CblasNoTrans) ? N : K;
	cublasOperation_t cuTransA = (TransA == CblasNoTrans) ? CUBLAS_OP_N : CUBLAS_OP_T;
	cublasOperation_t cuTransB = (TransB == CblasNoTrans) ? CUBLAS_OP_N : CUBLAS_OP_T;
	CUBLAS_CHECK(cublasSgemm(Context::cublas_handle(), cuTransB, cuTransA, N, M, K, &alpha, B, ldb, A, lda, &beta, C, N));
}

void matmul1D1D_gpu(const size_t dim_x, const size_t dim_y, const size_t dim_z, const float_t *A, const float_t *B, float_t *C) {
	const CBLAS_TRANSPOSE TransA = CblasNoTrans;
	const CBLAS_TRANSPOSE TransB = CblasNoTrans;
	sgemm_gpu(TransA, TransB, dim_x, dim_y, dim_z, 1.0, A, B, 0.0, C);
}

// the arguments of the maxima
int argmax_gpu(const size_t n, const float_t *x) {
	return 0;
}

void gemv_gpu(const CBLAS_TRANSPOSE TransA, const int M, const int N, 
	const float alpha, const float* A, const float* x, const float beta, float* y) {
	cublasOperation_t cuTransA = (TransA == CblasNoTrans) ? CUBLAS_OP_T : CUBLAS_OP_N;
	CUBLAS_CHECK(cublasSgemv(Context::cublas_handle(), cuTransA, N, M, &alpha, A, N, x, 1, &beta, y, 1));
}

void scal_gpu(const int N, const float alpha, float *X) {
	CUBLAS_CHECK(cublasSscal(Context::cublas_handle(), N, &alpha, X, 1));
}

void dot_gpu(const int n, const float* x, const float* y, float* out) {
	CUBLAS_CHECK(cublasSdot(Context::cublas_handle(), n, x, 1, y, 1, out));
}

void asum_gpu(const int n, const float* x, float* y) {
	CUBLAS_CHECK(cublasSasum(Context::cublas_handle(), n, x, 1, y));
}

void scale_gpu(const int n, const float alpha, const float *x, float* y) {
	CUBLAS_CHECK(cublasScopy(Context::cublas_handle(), n, x, 1, y, 1));
	CUBLAS_CHECK(cublasSscal(Context::cublas_handle(), n, &alpha, y, 1));
}

__global__ void set_kernel(const int n, const float_t alpha, float_t* y) {
	CUDA_KERNEL_LOOP(index, n) {
		y[index] = alpha;
	}
}

void set_gpu(const int N, const float_t alpha, float_t* Y) {
	if (alpha == 0) {
		CUDA_CHECK(cudaMemset(Y, 0, sizeof(float_t) * N));
		return;
	}
	set_kernel<<<CUDA_GET_BLOCKS(N), CUDA_NUM_THREADS>>>(N, alpha, Y);
}

__global__ void add_scalar_kernel(const int n, const float_t alpha, float_t* y) {
	CUDA_KERNEL_LOOP(index, n) {
		y[index] += alpha;
	}
}

void add_scalar_gpu(const int N, const float_t alpha, float_t* Y) {
	add_scalar_kernel<<<CUDA_GET_BLOCKS(N), CUDA_NUM_THREADS>>>(N, alpha, Y);
}

__global__ void vadd_kernel(const int n, const float_t* a, const float_t* b, float_t* y) {
	CUDA_KERNEL_LOOP(index, n) {
		y[index] = a[index] + b[index];
	}
}

void vadd_gpu(const int N, const float_t* a, const float_t* b, float_t* y) {
	vadd_kernel<<<CUDA_GET_BLOCKS(N), CUDA_NUM_THREADS>>>(N, a, b, y);
}

void softmax_cross_entropy_gpu(int x, int y, float_t *in_data, float_t *out_data) {
}
