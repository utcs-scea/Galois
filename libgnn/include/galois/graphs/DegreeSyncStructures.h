#include "galois/GNNTypes.h"

namespace galois {
namespace graphs {

extern uint32_t* gnn_degree_vec_1_;
extern uint32_t* gnn_degree_vec_2_;

extern galois::DynamicBitSet bitset_sampled_degrees_;
extern std::vector<galois::LargeArray<uint32_t>>* gnn_sampled_out_degrees_;

struct InitialDegreeSync {
  using ValTy = std::pair<uint32_t, uint32_t>;

  //! return a vector of floats to sync
  static ValTy extract(uint32_t lid, char&) {
    return std::make_pair(gnn_degree_vec_1_[lid], gnn_degree_vec_2_[lid]);
  }

  //! reduction is addition in this case; add received vector to
  //! own vector
  static bool reduce(uint32_t lid, char&, ValTy y) {
    gnn_degree_vec_1_[lid] += y.first;
    gnn_degree_vec_2_[lid] += y.second;
    if (y.first || y.second) {
      return true;
    } else {
      return false;
    }
  }

  //! No-op: readAny = overwritten anyways
  static void reset(uint32_t lid, char&) {
    gnn_degree_vec_1_[lid] = 0;
    gnn_degree_vec_2_[lid] = 0;
  }

  //! element wise set
  static void setVal(uint32_t lid, char&, ValTy y) {
    gnn_degree_vec_1_[lid] = y.first;
    gnn_degree_vec_2_[lid] = y.second;
  }

  // GPU options TODO for GPU
  static bool extract_batch(unsigned, uint8_t*, size_t*, DataCommMode*) {
    return false;
  }
  static bool extract_batch(unsigned, uint8_t*) { return false; }
  static bool reduce_batch(unsigned, uint8_t*, DataCommMode) { return false; }
  static bool reduce_mirror_batch(unsigned, uint8_t*, DataCommMode) {
    return false;
  }
  static bool setVal_batch(unsigned, uint8_t*, DataCommMode) { return false; }
  static bool extract_reset_batch(unsigned, uint8_t*, size_t*, DataCommMode*) {
    return false;
  }
  static bool extract_reset_batch(unsigned, uint8_t*) { return false; }
};

struct SubgraphDegreeSync {
  using ValTy = galois::gstl::Vector<uint32_t>;

  //! return a vector of floats to sync
  static ValTy extract(uint32_t lid, char&) {
    ValTy vec_to_send(gnn_sampled_out_degrees_->size());
    size_t count = 0;
    for (galois::LargeArray<uint32_t>& layer_degrees :
         *gnn_sampled_out_degrees_) {
      vec_to_send[count++] = layer_degrees[lid];
    }
    assert(count == vec_to_send.size());
    return vec_to_send;
  }

  static bool reduce(uint32_t lid, char&, ValTy y) {
    assert(y.size() == gnn_sampled_out_degrees_->size());
    for (size_t degree_index = 0; degree_index < y.size(); degree_index++) {
      (*gnn_sampled_out_degrees_)[degree_index][lid] += y[degree_index];
    }
    return true;
  }

  //! No-op: readAny = overwritten anyways; can probably get away with no-op
  static void reset(uint32_t lid, char&) {
    for (galois::LargeArray<uint32_t>& layer_degrees :
         *gnn_sampled_out_degrees_) {
      layer_degrees[lid] = 0;
    }
  }

  //! element wise set
  static void setVal(uint32_t lid, char&, ValTy y) {
    assert(y.size() == gnn_sampled_out_degrees_->size());
    for (size_t degree_index = 0; degree_index < y.size(); degree_index++) {
      (*gnn_sampled_out_degrees_)[degree_index][lid] = y[degree_index];
    }
  }

  // GPU options TODO for GPU
  static bool extract_batch(unsigned, uint8_t*, size_t*, DataCommMode*) {
    return false;
  }
  static bool extract_batch(unsigned, uint8_t*) { return false; }
  static bool reduce_batch(unsigned, uint8_t*, DataCommMode) { return false; }
  static bool reduce_mirror_batch(unsigned, uint8_t*, DataCommMode) {
    return false;
  }
  static bool setVal_batch(unsigned, uint8_t*, DataCommMode) { return false; }
  static bool extract_reset_batch(unsigned, uint8_t*, size_t*, DataCommMode*) {
    return false;
  }
  static bool extract_reset_batch(unsigned, uint8_t*) { return false; }
};

struct SubgraphDegreeBitset {
  static constexpr bool is_vector_bitset() { return false; }
  static constexpr bool is_valid() { return true; }
  static galois::DynamicBitSet& get() { return bitset_sampled_degrees_; }
  static void reset_range(size_t begin, size_t end) {
    bitset_sampled_degrees_.reset(begin, end);
  }
};

} // namespace graphs
} // namespace galois
