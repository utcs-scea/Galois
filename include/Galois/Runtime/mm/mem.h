// heap building blocks -*- C++ -*-
// strongly inspired by heap layers:
// http://www.heaplayers.org/
// FSB is modified from
// http://warp.povusers.org/FSBAllocator/

#include "Support/ThreadSafe/simple_lock.h"
#include "Galois/Runtime/PerCPU.h"
#include <iostream>
#include <string.h>

#ifndef _GALOISRUNTIME_MM_MEM_H
#define _GALOISRUNTIME_MM_MEM_H

namespace GaloisRuntime {
namespace MM {

//Base mmap wrapper
class mmapWrapper {
  static void* _alloc();
  static void _free(void*);
public:
  enum {AllocSize = 2*1024*1024,
	Alighment = 4*1024,
	AutoFree = 1};
  
  void* allocate(unsigned int size) {
    assert(size % AllocSize == 0);
    return _alloc();
  }

  void deallocate(void* ptr) {
    _free(ptr);
  }
};


//Per-thread heaps using Galois thread aware construct
template<class LocalHeap>
class ThreadAwarePrivateHeap {
  CPUSpaced<LocalHeap> heaps;

public:
  enum { AllocSize = LocalHeap::AllocSize };

  ThreadAwarePrivateHeap() : heaps(LocalHeap::merge) {}

  inline void* allocate(unsigned int size) {
    return heaps.get().allocate(size);
  }

  inline void deallocate(void* ptr) {
    heaps.get().deallocate(ptr);
  }
};

//Apply a lock to a heap
template<class RealHeap>
class LockedHeap : public RealHeap {
  threadsafe::simpleLock<int, true> lock;
public:
  enum { AllocSize = RealHeap::AllocSize };

  inline void* allocate(unsigned int size) {
    lock.lock();
    void* retval = RealHeap::allocate(size);
    lock.unlock();
    return retval;
  }
  
  inline void deallocate(void* ptr) {
    lock.lock();
    RealHeap::deallocate(ptr);
    lock.unlock();
  }
};

template<typename SourceHeap>
class ZeroOut : public SourceHeap {
public:
  enum { AllocSize = SourceHeap::AllocSize } ;
  inline void* allocate(unsigned int size) {
    void* retval = SourceHeap::allocate(size);
    memset(retval, 0, size);
    return retval;
  }

  inline void deallocate(void* ptr) {
    SourceHeap::deallocate(ptr);
  }
};

//Add a header to objects
template<typename Header, typename SourceHeap>
class AddHeader : public SourceHeap {
  enum { offset = (sizeof(Header) + (sizeof(double) - 1)) & ~(sizeof(double) - 1) };

public:
  inline void* allocate(unsigned int size) {
    //First increase the size of the header to be aligned to a double
    void* ptr = SourceHeap::allocate(size + offset);
    //Now return the offseted pointer
    return (char*)ptr + offset;
  }
  
  inline void deallocate(void* ptr) {
    SourceHeap::deallocate(getHeader(ptr));
  }

  inline static Header* getHeader(void* ptr) {
    return (Header*)((char*)ptr - offset);
  }
};

//Allow looking up parent heap pointers
template<class SourceHeap>
class OwnerTaggedHeap : public AddHeader<void*, SourceHeap> {
  typedef AddHeader<OwnerTaggedHeap*, SourceHeap> Src;
public:
  inline void* allocate(unsigned int size) {
    void* retval = Src::allocate(size);
    *(Src::getHeader(retval)) = this;
    return retval;
  }

  inline void deallocate(void* ptr) {
    assert(*(Src::getHeader(ptr)) == this);
    Src::deallocate(ptr);
  }

  inline static OwnerTaggedHeap* owner(void* ptr) {
    return *(OwnerTaggedHeap**)Src::getHeader(ptr);
  }
};

//Maintain a freelist
template<class SourceHeap>
class FreeListHeap : public SourceHeap {
  struct FreeNode {
    FreeNode* next;
  };
  FreeNode* head;

public:
  enum { AllocSize = SourceHeap::AllocSize };

  void clear() {
    while (head) {
      FreeNode* N = head;
      head = N->next;
      SourceHeap::deallocate(N);
    }
  }

  FreeListHeap() : head(0) {}
  ~FreeListHeap() {
    clear();
  }

  inline void* allocate(unsigned int size) {
    if (head) {
      void* ptr = head;
      head = head->next;
      return ptr;
    }
    return SourceHeap::allocate(size);
  }

  inline void deallocate(void* ptr) {
    if (!ptr) return;
    FreeNode* NH = (FreeNode*)ptr;
    NH->next = head;
    head = NH;
  }

};

//Maintain a freelist using a lock which doesn't cover SourceHeap
template<class SourceHeap>
class SelfLockFreeListHeap : public SourceHeap {
  struct FreeNode {
    FreeNode* next;
  };
  FreeNode* head;

public:
  enum { AllocSize = SourceHeap::AllocSize };

  void clear() {
    FreeNode* h = 0;
    do {
      h = head;
    } while (!__sync_bool_compare_and_swap(&head, h, 0));
    while (h) {
      FreeNode* N = h;
      h = N->next;
      SourceHeap::deallocate(N);
    }
  }

  SelfLockFreeListHeap() : head(0) {}
  ~SelfLockFreeListHeap() {
    clear();
  }

  inline void* allocate(unsigned int size) {
    static threadsafe::simpleLock<int, true> lock;

    lock.lock();
    FreeNode* OH = 0;
    FreeNode* NH = 0;
    do {
      OH = head;
      if (!OH) {
	lock.unlock();
	return SourceHeap::allocate(size);
      }
      NH = OH->next; //The lock protects this line
    } while (!__sync_bool_compare_and_swap(&head, OH, NH));
    lock.unlock();
    assert(OH);
    return (void*)OH;
  }

  inline void deallocate(void* ptr) {
    if (!ptr) return;
    FreeNode* OH;
    FreeNode* NH;
    do {
      OH = head;
      NH = (FreeNode*)ptr;
      NH->next = OH;
    } while (!__sync_bool_compare_and_swap(&head, OH, NH));
  }

};

template<unsigned ElemSize, typename SourceHeap>
class BlockAlloc : public SourceHeap {

  struct TyEq {
    double data[((ElemSize + sizeof(double) - 1) & ~(sizeof(double) - 1))/sizeof(double)];
  };

  struct Block_basic {
    union {
      Block_basic* next;
      double dummy;
    };
    TyEq data[1];
  };

  enum {BytesLeft = (SourceHeap::AllocSize - sizeof(Block_basic)),
	BytesLeftR = BytesLeft & ~(sizeof(double) - 1),
	FitLeft = BytesLeftR / sizeof(TyEq[1]),
	TotalFit = FitLeft + 1
  };

  struct Block {
    union {
      Block* next;
      double dummy;
    };
    TyEq data[TotalFit];
  };

  Block* head;
  int headIndex;

  void refill() {
    void* P = SourceHeap::allocate(SourceHeap::AllocSize);
    Block* BP = (Block*)P;
    BP->next = head;
    head = BP;
    headIndex = 0;
  }
public:
  enum { AllocSize = ElemSize };

  void clear() {
    while(head) {
      Block* B = head;
      head = B->next;
      SourceHeap::deallocate(B);
    }
  }

  BlockAlloc() :SourceHeap(), head(0), headIndex(0) {
    //    std::cerr << "BA " << TotalFit << " " << ElemSize << " " << sizeof(TyEq) << " " << sizeof(Block) << " " << SourceHeap::AllocSize << "\n";
    assert(sizeof(Block) <= SourceHeap::AllocSize);
  }

  ~BlockAlloc() {
    clear();
  }

  inline void* allocate(unsigned int size) {
    assert(size == ElemSize);
    if (!head || headIndex == TotalFit)
      refill();
    return &head->data[headIndex++];
  }

  inline void deallocate(void* ptr) {}

  static void merge(BlockAlloc& lhs, BlockAlloc& rhs) {
    if (!lhs.head) { //Steal
      lhs.head = rhs.head;
      lhs.headIndex = rhs.headIndex;
    } else { //Merge
      Block* B = lhs.head;
      while (B->next) { B = B->next; }
      B->next = rhs.head;
    }
    rhs.head = 0;
    rhs.headIndex = 0;
  }
};

//This implements a bump pointer though chunks of memory
template<typename SourceHeap>
class SimpleBumpPtr : public SourceHeap {

  struct Block {
    union {
      Block* next;
      double dummy; // for alignment
    };
  };

  Block* head;
  int offset;

  void refill() {
    void* P = SourceHeap::allocate(SourceHeap::AllocSize);
    Block* BP = (Block*)P;
    BP->next = head;
    head = BP;
    offset = sizeof(Block);
  }
public:
  enum { AllocSize = 0 };

  SimpleBumpPtr() :SourceHeap(), head(0), offset(0) {}
  ~SimpleBumpPtr() {
    clear();
  }

  void clear() {
    while(head) {
      Block* B = head;
      head = B->next;
      SourceHeap::deallocate(B);
    }
  }

  inline void* allocate(unsigned int size) {
    //increase to alignment
    size = (size + sizeof(double) - 1) & ~(sizeof(double) - 1);
    //Check current block
    if (!head || offset + size > SourceHeap::AllocSize)
      refill();
    //Make sure this will fit
    if (offset + size > SourceHeap::AllocSize) {
      assert(0 && "Too large");
      return 0;
    }
    char* retval = (char*)head;
    retval += offset;
    offset += size;
    return retval;
  }

  inline void deallocate(void* ptr) {}

  static void merge(SimpleBumpPtr& lhs, SimpleBumpPtr& rhs) {
    if (!lhs.head) { //Steal
      lhs.head = rhs.head;
      lhs.offset = rhs.offset;
    } else { //Merge
      Block* B = lhs.head;
      while (B->next) { B = B->next; }
      B->next = rhs.head;
    }
    rhs.head = 0;
    rhs.offset = 0;
  }
};

//This is the base source of memory for all allocators
//It maintains a freelist of hunks acquired from the system
template<typename RealBase>
class SystemBaseAllocator {
  static SelfLockFreeListHeap<RealBase> Source;
public:
  enum { AllocSize = SelfLockFreeListHeap<RealBase>::AllocSize };

  inline void* allocate(unsigned int size) {
    return Source.allocate(size);
  }

  inline void deallocate(void* ptr) {
    Source.deallocate(ptr);
  }
};

typedef SystemBaseAllocator<mmapWrapper> SystemBaseAlloc;





////////////////////////////////////////////////////////////////////////////////
// Now adapt to standard std allocators
////////////////////////////////////////////////////////////////////////////////

//A fixed size block allocator
template<typename Ty>
class FSBGaloisAllocator
{
  ThreadAwarePrivateHeap<BlockAlloc<sizeof(Ty), SystemBaseAlloc> > Alloc;

public:
  typedef size_t size_type;
  typedef ptrdiff_t difference_type;
  typedef Ty *pointer;
  typedef const Ty *const_pointer;
  typedef Ty& reference;
  typedef const Ty& const_reference;
  typedef Ty value_type;
  
  pointer address(reference val) const { return &val; }
  const_pointer address(const_reference val) const { return &val; }
  
  template<class Other>
  struct rebind
  {
    typedef FSBGaloisAllocator<Other> other;
  };

  template <class U>
  FSBGaloisAllocator ( const FSBGaloisAllocator<U>& ) throw() {}

  FSBGaloisAllocator() throw() {}
  
  pointer allocate(int x)
  {
    assert(x == 1);
    return static_cast<pointer>(Alloc.allocate(sizeof(Ty)));
  }
  
  void deallocate(pointer ptr, size_type)
  {
    Alloc.deallocate(ptr);
  }
  
  template<typename TyC>
  void construct(pointer ptr, const TyC& val)
  {
    new ((void *)ptr) Ty(val);
  }
  
  void destroy(pointer ptr)
  {
    ptr->Ty::~Ty();
  }
  
  size_type max_size() const throw() { return 1; }
};

//Keep a reference to an external allocator
template<typename Ty, typename AllocTy>
class ExternRefGaloisAllocator
{
public:
  AllocTy* Alloc; // Should be private except that makes copy hard

  typedef size_t size_type;
  typedef ptrdiff_t difference_type;
  typedef Ty *pointer;
  typedef const Ty *const_pointer;
  typedef Ty& reference;
  typedef const Ty& const_reference;
  typedef Ty value_type;
  
  pointer address(reference val) const { return &val; }
  const_pointer address(const_reference val) const { return &val; }
  
  template<class Other>
  struct rebind
  {
    typedef ExternRefGaloisAllocator<Other, AllocTy> other;
  };

  template <class U>
  ExternRefGaloisAllocator ( const ExternRefGaloisAllocator<U,AllocTy>& rhs)
    throw() {
    Alloc = rhs.Alloc;
  }

  explicit ExternRefGaloisAllocator(AllocTy* a) throw() 
  :Alloc(a)
  {}
  
  pointer allocate(int x)
  {
    return static_cast<pointer>(Alloc->allocate(x*sizeof(Ty)));
  }
  
  void deallocate(pointer ptr, size_type x)
  {
    Alloc->deallocate(ptr);
  }
  
  template<typename TyC>
  void construct(pointer ptr, const TyC& val)
  {
    new ((void *)ptr) Ty(val);
  }
  
  void destroy(pointer ptr)
  {
    ptr->Ty::~Ty();
  }
  
  size_type max_size() const throw() { return 1024*1024; }
};

}
}

#endif
