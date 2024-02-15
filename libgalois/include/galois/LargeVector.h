#ifndef GALOIS_LARGEVECTOR_H
#define GALOIS_LARGEVECTOR_H

#include <linux/memfd.h>
#include <linux/mman.h>
#include <sys/mman.h>
#include <unistd.h>
#include <atomic>
#include <cstddef>
#include <stdexcept>
#include <iterator>
#include <iostream>
#include <list>
#include <string.h>
#include <utility>

/*
 * A vector backed by huge pages. Guarantees addresss stability, so values do
 * not have to be moveable.
 *
 */
template <typename T>
class LargeVector {
protected:
  std::atomic_flag m_lock    = ATOMIC_FLAG_INIT;
  volatile size_t m_capacity = 1;
  volatile size_t m_size     = 0;
  T* volatile m_data         = nullptr;

  const int m_fd;
  std::list<std::pair<void*, size_t>> m_mappings;

  void lock() {
    while (m_lock.test_and_set(std::memory_order_acquire))
      ;
    ;
  }

  void unlock() { m_lock.clear(std::memory_order_release); }

  void grow() {
    // assumes lock is held or constructor
    m_capacity *= 2;

    const size_t file_size =
        (sizeof(T) * m_capacity + 2097152 - 1) & (~(2097152 - 1));

    if (ftruncate(m_fd, file_size) == -1)
      throw std::runtime_error(std::string("ftruncate: ") + strerror(errno));

    if (m_capacity * sizeof(T) > 0) {
      m_data = static_cast<T*>(
          mmap(nullptr, m_capacity * sizeof(T), PROT_READ | PROT_WRITE,
               MAP_SHARED | MAP_HUGETLB | MAP_HUGE_2MB, m_fd, 0));
      if (m_data == MAP_FAILED) {
        throw std::runtime_error(std::string("mmap failed: ") +
                                 strerror(errno));
      }

      m_mappings.push_back(std::make_pair(m_data, m_capacity * sizeof(T)));
    }
  }

public:
  LargeVector()
      : m_fd(memfd_create("LargeVector", MFD_HUGETLB | MFD_HUGE_2MB)) {
    if (m_fd == -1) {
      throw std::runtime_error(std::string("creating memfd: ") +
                               strerror(errno));
    }

    grow();
  }

  ~LargeVector() {
    for (; !m_mappings.empty(); m_mappings.pop_front())
      munmap(m_mappings.front().first, m_mappings.front().second);

    close(m_fd);
  }

  // no consistency guarantees with push_back, but size only increases
  uint64_t size() const noexcept { return m_size; }

  template <typename... Args>
  T& emplace_back(Args&&... args) {
    lock();
    if (m_size == m_capacity) {
      grow();
    }
    T& ret = *new (m_data + m_size++) T(std::forward<Args>(args)...);
    unlock();
    return ret;
  }

  T& push_back(const T& t) { return emplace_back(t); }
  T& push_back(T&& t) { return emplace_back(std::move(t)); }

  T& operator[](size_t index) const { return m_data[index]; }

  class iterator : public std::iterator<std::random_access_iterator_tag, T> {
  public:
    using iterator_category = std::random_access_iterator_tag;
    using value_type        = const T;
    using difference_type   = std::ptrdiff_t;
    using pointer           = const T*;
    using reference         = const T&;

  private:
    T* ptr;

  public:
    iterator(T* p) : ptr(p) {}

    iterator& operator++() {
      ++ptr;
      return *this;
    }

    iterator operator++(int) {
      iterator temp = *this;
      ++*this;
      return temp;
    }

    bool operator==(const iterator& other) const { return ptr == other.ptr; }

    bool operator!=(const iterator& other) const { return ptr != other.ptr; }

    T& operator*() const { return *ptr; }
  };

  iterator begin() { return iterator(m_data); }

  iterator end() { return iterator(m_data + m_size); }
};

namespace std {
template <typename T>
ostream& operator<<(std::ostream& os, const LargeVector<T>& vec) {
  for (uint64_t i = 0; i < vec.getSize(); i++) {
    os << vec[i];
    if (i < vec.getSize() - 1) {
      os << " ";
    }
  }
  return os;
}

template <typename T>
istream& operator>>(istream& is, LargeVector<T>& vec) {
  T value;
  while (is >> value) {
    vec.push_back(value);
  }
  return is;
}
} // namespace std

#endif