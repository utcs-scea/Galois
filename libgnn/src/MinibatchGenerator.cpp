#include "galois/MinibatchGenerator.h"
#include "galois/Galois.h"
#include <cassert>

void galois::MinibatchGenerator::OriginalGetNextMinibatch(
    std::vector<char>* batch_mask) {
  assert(current_position_ <= mask_to_minibatch_.size());
  assert(current_position_ <= master_bound_);
  assert(batch_mask->size() == mask_to_minibatch_.size());

  galois::ParallelSTL::fill(batch_mask->begin(), batch_mask->end(), 0);
  if (current_position_ >= master_bound_) {
    return;
  }

  size_t current_count = 0;
  // start from last positiion
  while (current_position_ < master_bound_) {
    if (mask_to_minibatch_[current_position_]) {
      // XXX and a master node; seed nodes only exist locally
      (*batch_mask)[current_position_] = 1;
      current_count++;
    }
    // break when minibatch is large enough
    current_position_++;
    if (current_count == minibatch_size_)
      break;
  }

  // advance current position to next set bit for next call (or to end to detect
  // no more minibatches
  while (!mask_to_minibatch_[current_position_] &&
         (current_position_ < master_bound_)) {
    current_position_++;
  }
}

void galois::MinibatchGenerator::ShuffleGetNextMinibatch(
    std::vector<char>* batch_mask) {
  size_t current_count = 0;
  galois::ParallelSTL::fill(batch_mask->begin(), batch_mask->end(), 0);
  // loops through a number of indices locally and sets
  while (current_position_ < all_indices_.size()) {
    (*batch_mask)[all_indices_[current_position_++]] = 1;
    current_count++;
    if (current_count == minibatch_size_)
      break;
  }
}

// used if all hosts have a global view of the same minibatch sequence
// (occurs if all hosts use same shuffle seed)
// Do not use unless you know what you are doing
void galois::MinibatchGenerator::DistributedShuffleGetNextMinibatch(
    std::vector<char>* batch_mask) {
  galois::ParallelSTL::fill(batch_mask->begin(), batch_mask->end(), 0);

  size_t current_count = 0;
  size_t global_minibatch_size =
      minibatch_size_ * galois::runtime::getSystemNetworkInterface().Num;
  while (current_position_ < all_indices_.size()) {
    size_t candidate_lid = all_indices_[current_position_++];
    if (candidate_lid < batch_mask->size() && candidate_lid < master_bound_) {
      (*batch_mask)[candidate_lid] = 1;
    }

    current_count++;
    if (current_count == global_minibatch_size)
      break;
  }
}

// used with distributed minibatch tracker which is deprecated; code not
// guaranteed to work
void galois::MinibatchGenerator::DistributedShuffleGetNextMinibatch(
    std::vector<char>* batch_mask, size_t num_to_get) {
  size_t current_count = 0;
  galois::ParallelSTL::fill(batch_mask->begin(), batch_mask->end(), 0);
  while (current_position_ < all_indices_.size()) {
    (*batch_mask)[all_indices_[current_position_++]] = 1;
    current_count++;
    if (current_count == num_to_get)
      break;
  }
}
