#include <sys/mman.h>
#include <cstddef>
#include <stdexcept>
#include <cstring>
#include <atomic>
#include <iterator>
#include <algorithm>

template <typename T>
class MyVector {
    public:
        T* data;
        uint64_t size;
        uint64_t capacity;
        std::atomic_flag lock = ATOMIC_FLAG_INIT;

        void lockSpin() {
            while (lock.test_and_set(std::memory_order_acquire)) {
                
            }
        }

        void unlockSpin() {
            lock.clear(std::memory_order_release);
        }

        void resize (uint64_t new_capacity) {
            lockSpin();
            if (new_capacity > capacity && new_capacity != capacity) {
                // std::cout << "Resizing" <<std::endl;
                // std::cout << "Current capacity " << capacity << std::endl;
                // std::cout << "New capacity " << 2*capacity << std::endl;
                // std::cout << "Current address " << std::hex << data << std::endl;
                void* newAddress = static_cast<char*>(static_cast<void*>(data)) + capacity*sizeof(T);
                // std::cout << "Additional block at " << std::hex << newAddress << std::endl;
                T* new_data = static_cast<T*>(mmap(newAddress, capacity*sizeof(T), PROT_READ | PROT_WRITE, MAP_FIXED | MAP_PRIVATE | MAP_ANONYMOUS | MAP_POPULATE | MAP_HUGETLB | MAP_HUGE_2MB, -1, 0));
                if (new_data == MAP_FAILED) {
                    throw std::runtime_error("Error resizing memory map: " + std::string(strerror(errno)));
                }
                capacity = 2*capacity;
            }
            unlockSpin();
        }

        MyVector() : data(nullptr), size(1), capacity(1) {
            // std::cout << "Constructor called" << std::endl;
        }

        void setSize (uint64_t _capacity) {
            uintptr_t hexAddress = 0x7F0000000000;
            void* desiredAddress = reinterpret_cast<void*>(hexAddress);
            capacity = 2*1024*1024/sizeof(T);
            size = 0;
            // std::cout << "Size of T is " << sizeof(T) << std::endl;
            data = static_cast<T*>(mmap(desiredAddress, capacity*sizeof(T), PROT_READ | PROT_WRITE, MAP_FIXED | MAP_PRIVATE | MAP_ANONYMOUS | MAP_POPULATE | MAP_HUGETLB | MAP_HUGE_2MB, -1, 0));
            if (data == MAP_FAILED) {
                throw std::runtime_error("Error mapping memory " + std::string(strerror(errno)));
            }
            // else {
            //     std::cout << "Memory allocated successfully" << std::endl;
            // }
            auto temp = data[0];
            // std::cout << "Can read data" << std::endl;
        }

        ~MyVector () {
            munmap (data, capacity*sizeof(T));
        }

        void push_back (const T& value) {
            if (size == capacity) {
                resize(capacity==0 ? 1 : capacity*2);
            }
            data[size++] = value;
        }

        T& operator[](uint64_t index) {
            if (index >= capacity) {
                // std::cout << "Coming here" << std::endl;
                resize(std::max(2*capacity, index+1));
                size = std::max(size, index);
            }
            return data[index];
        }

        uint64_t getSize() const {
            return size;
        }

        // Currently, the iterator is not thread safe? Does every thread have a separate iterator?
        class iterator : public std::iterator<std::forward_iterator_tag, T> {
            
            public:
                using iterator_category = std::forward_iterator_tag;
                using value_type = const T;
                using difference_type = std::ptrdiff_t;
                using pointer = const T*;
                using reference = const T&;

            private:
                T* ptr;

            public:
                iterator (T* p) : ptr(p) {}

                iterator& operator++() {
                    ++ptr;
                    return *this;
                }

                iterator operator++(int) {
                    iterator temp = *this;
                    ++*this;
                    return temp;
                }

                bool operator==(const iterator &other) const {
                    return ptr == other.ptr;
                }

                bool operator!=(const iterator &other) const {
                    return ptr != other.ptr;
                }

                T& operator*() const {
                    return *ptr;
            }
        };

        iterator begin() {
            return iterator(data);
        }

        iterator end() {
            return iterator(data + size);
        }
};

template <typename T>
std::ostream& operator<<(std::ostream& os, const MyVector<T>& vec) {
    for (uint64_t i=0; i<vec.getSize(); i++) {
        os << vec[i];
        if (i < vec.getSize() - 1) {
            os << " ";
        }
    }
    return os;
}

template <typename T>
std::istream& operator>>(std::istream& is, MyVector<T>& vec) {
    T value;
    while (is >> value) {
        vec.push_back(value);
    }
    return is;
}

template <typename T>
void swap (MyVector<T>& first, MyVector<T>& second) {
    if (first.getSize() == second.getSize()) {
        std::swap(first.data, second.data);
        std::swap(first.size, second.size);
        std::swap(first.capacity, second.capacity);
    }
    else {
        throw std::invalid_argument("Vectors must be of same size to swap");
    }
}