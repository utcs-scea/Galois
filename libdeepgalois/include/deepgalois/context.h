#pragma once
/**
 * Based on common.hpp file of the Caffe deep learning library.
 */

#include <string>
#include <cassert>
#include "deepgalois/types.h"
#include "deepgalois/reader.h"
#include "deepgalois/GraphTypes.h"

#ifdef __GALOIS_HET_CUDA__
#include "deepgalois/cutils.h"
#endif

namespace deepgalois {

class Context {
public:
  Context();
  //! initializer for gpu; goes ahead and sets a few things
  Context(bool use_gpu)
      : is_device(use_gpu), n(0), num_classes(0), feat_len(0),
        is_single_class(true), is_selfloop_added(false), use_subgraph(false),
        h_labels(NULL), h_feats(NULL), d_labels(NULL), d_labels_subg(NULL),
        d_feats(NULL), d_feats_subg(NULL), norm_factors(NULL) {}
  ~Context();

  size_t read_graph(bool selfloop);
  size_t read_labels() {
    num_classes = reader.read_labels(is_single_class, h_labels);
    return num_classes;
  }
  size_t read_features() {
    feat_len = reader.read_features(h_feats);
    return feat_len;
  }
  size_t read_masks(std::string mask_type, size_t n, size_t& begin, size_t& end,
                    mask_t* masks) {
    return reader.read_masks(mask_type, n, begin, end, masks);
  }

  label_t get_label(size_t i) {
    return h_labels[i];
  } // single-class (one-hot) label
  // label_t get_label(size_t i, size_t j) { return labels[i*num_classes+j]; }
  // // multi-class label
  float_t* get_norm_factors_ptr() { return norm_factors; }
  float_t* get_norm_factors_subg_ptr() { return &norm_factors_subg[0]; }

  void set_dataset(std::string dataset_str) {
    dataset = dataset_str;
    reader.init(dataset);
  }
  void set_label_class(bool is_single = true) { is_single_class = is_single; }
  void set_use_subgraph(bool use_subg) { use_subgraph = use_subg; }
  void copy_data_to_device(); // copy labels and input features
  void norm_factor_computing(bool is_subgraph, int subg_id = 0);
  void gen_subgraph_labels(size_t m, const mask_t* masks);
  void gen_subgraph_feats(size_t m, const mask_t* masks);
  void createSubgraphs(int num_subgraphs);

#ifndef __GALOIS_HET_CUDA__
  Graph* graph_cpu; // the input graph, |V| = N
  std::vector<Graph*> subgraphs_cpu;
  void add_selfloop(Graph& og, Graph& g);
  //! returns pointer to the graph
  Graph* getGraphPointer() { return graph_cpu; }
  Graph* getSubgraphPointer(int id) { return subgraphs_cpu[id]; };
  float_t* get_feats_ptr() { return h_feats; }
  float_t* get_feats_subg_ptr() { return &h_feats_subg[0]; }
  label_t* get_labels_ptr() { return h_labels; }
  label_t* get_labels_subg_ptr() { return &h_labels_subg[0]; }
#else
  GraphGPU graph_gpu; // the input graph, |V| = N
  std::vector<GraphGPU*> subgraphs_gpu;
  GraphGPU* getGraphPointer() { return &graph_gpu; }
  GraphGPU* getSubgraphPointer(int id) { return subgraphs_gpu[id]; };
  float_t* get_feats_ptr() { return d_feats; }
  float_t* get_feats_subg_ptr() { return d_feats_subg; }
  label_t* get_labels_ptr() { return d_labels; }
  label_t* get_labels_subg_ptr() { return d_labels_subg; }
  inline static cublasHandle_t cublas_handle() { return cublas_handle_; }
  inline static cusparseHandle_t cusparse_handle() { return cusparse_handle_; }
  inline static cusparseMatDescr_t cusparse_matdescr() {
    return cusparse_matdescr_;
  }
  inline static curandGenerator_t curand_generator() {
    return curand_generator_;
  }
#endif

protected:
  std::string dataset;
  bool is_device;         // is this on device or host
  size_t n;               // number of samples: N
  size_t num_classes;     // number of classes: E
  size_t feat_len;        // input feature length: D
  bool is_single_class;   // single-class (one-hot) or multi-class label
  bool is_selfloop_added; // whether selfloop is added to the input graph
  bool use_subgraph;      // whether to use subgraph
  label_t* h_labels;      // labels for classification. Single-class label: Nx1,
                          // multi-class label: NxE
  float_t* h_feats;       // input features: N x D
  // label_t *h_labels_subg;      // labels for subgraph
  // float_t* h_feats_subg;       // input features for subgraph
  label_t* d_labels;      // labels on device
  label_t* d_labels_subg; // labels for subgraph on device
  float_t* d_feats;       // input features on device
  float_t* d_feats_subg;  // input features for subgraph on device
  float_t* norm_factors;  // normalization constant based on graph structure
  std::vector<label_t> h_labels_subg;     // labels for subgraph
  std::vector<float_t> h_feats_subg;      // input features for subgraph
  std::vector<float_t> norm_factors_subg; // normalization constant for subgraph
  // float_t* norm_factors_subg;  // normalization constant for subgraph
  Reader reader;

  void alloc_norm_factor();
  void alloc_subgraph_norm_factor(int subg_id);

#ifndef __GALOIS_HET_CUDA__
  void read_edgelist(const char* filename, bool symmetrize = false,
                     bool add_self_loop = false);
#else
  static cublasHandle_t cublas_handle_;         // used to call cuBLAS
  static cusparseHandle_t cusparse_handle_;     // used to call cuSPARSE
  static cusparseMatDescr_t cusparse_matdescr_; // used to call cuSPARSE
  static curandGenerator_t
      curand_generator_; // used to generate random numbers on GPU
#endif
};

} // namespace deepgalois
