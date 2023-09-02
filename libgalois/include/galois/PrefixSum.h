#include <cstdint>
#include "galois/WaterFallLock.h"
#include "galois/Galois.h"
#include <variant>

using namespace galois;

inline void empty(std::monostate& a, uint64_t i){ (void) a; (void) i;}

template <typename T>
inline T equalizer(const T& t){return t;}

template<typename T>
inline void before (T l, uint64_t tid){l.template wait<1>(tid);}
template<typename T>
inline void after  (T l, uint64_t tid){l.template done<2>(tid);}

/*
 * @param src the source array
 * @param dst the destination array
 * @param lock a reference to a WaterFallLock, which should have length of the number of threads
 * @param paste a conduit assigned per thread in order to ensure cache_line padding for speed
 */
template <typename A, typename B,
          B (*transmute) (const A&), B (*scan_op) (const A& x, const B& y),
          B (*combiner) (const B& x, const B& y), template <typename C> typename Conduit>
class PrefixSum
{
  A* src;
  B* dst;
  using PArr = Conduit<B>;
  Conduit<B> paste;
  using WFLType = galois::WaterFallLock<Conduit<unsigned>>;
  WFLType lock;

  template<typename T>
  struct Arr
  {
    T* arr;

    Arr(T* arr) : arr(arr) {}

    template<typename i_type>
    T& operator[](i_type i) {return arr[i];}
  };


  template <typename A1, typename A2, A2 (*trans) (const A1&), A2 (*scan) (const A1& x, const A2& y),
           typename CTX, void (*before) (CTX&, uint64_t), void (*after) (CTX&, uint64_t),
            template<typename C> typename Holder, bool combine = false>
    inline void serial_pfxsum(Holder<A1> src, Holder<A2> dst, uint64_t ns, CTX ctx)
  {
    if(!combine) dst[0] = trans(src[0]);
    for(uint64_t i = 1; i < ns; i++)
    {
      before(ctx, i);
      dst[i] = scan(src[i], dst[i-1]);
      after(ctx, i);
    }
  }

  inline void parallel_pfxsum_phase_0(A* src, B* dst, uint64_t ns, B& paste_loc, uint64_t wfl_id)
  {
    serial_pfxsum<A, B, transmute, scan_op, std::monostate, empty, empty, Arr>(src, dst, ns, std::monostate());
    paste_loc = dst[ns - 1];
    lock.template done<1>(wfl_id);
  }

  inline void parallel_pfxsum_phase_1(uint64_t ns, uint64_t wfl_id)
  {
    if(!wfl_id)
    {
      lock.template done<2>(wfl_id);
      serial_pfxsum<B, B, equalizer, combiner, WFLType&, before<WFLType&>, after<WFLType&>, Conduit>(paste, paste, ns, lock);
    }
    else {
      lock.template wait<2>(wfl_id - 1);
    }
  }

  inline void parallel_pfxsum_phase_2(A* src, B* dst, uint64_t ns, const B& phase1_val, bool pfxsum)
  {
    if(pfxsum)
    {
      dst[0] = scan_op(src[0], phase1_val);
      serial_pfxsum<A, B, transmute, scan_op, std::monostate, empty, empty, Arr, true>(src, dst, ns, std::monostate());
    }
    else
    {
      for(uint64_t i = 0; i < ns; i++) dst[i] = combiner(phase1_val, dst[i]);
    }
  }

  inline void parallel_pfxsum_work(
      uint64_t phase0_ind, uint64_t phase0_sz,
      uint64_t phase2_ind, uint64_t phase2_sz,
      uint64_t wfl_id, uint64_t nt)
  {

    parallel_pfxsum_phase_0(&src[phase0_ind], &dst[phase0_ind], phase0_sz,
         paste[wfl_id], wfl_id);

    parallel_pfxsum_phase_1(nt, wfl_id);

    const B& paste_val = paste[wfl_id ? wfl_id - 1 : nt - 1];
    parallel_pfxsum_phase_2(&src[phase2_ind], &dst[phase2_ind], phase2_sz, paste_val, !wfl_id);
  }

  /*
   * @param ns the number of items to sum
   * @param wf_id this corresponds to the thread id
   * @param nt this is the number of threads
   */
  void parallel_pfxsum_op(uint64_t ns, uint64_t wf_id, uint64_t nt)
  {
    uint64_t div_sz = ns / (nt + 1);
    uint64_t bigs   = ns % (nt + 1);
    uint64_t mid    = nt >> 1;
    bool     is_mid = mid == wf_id;
    //Concentrate the big in the middle thread
    uint64_t phase0_sz  = is_mid ? div_sz + bigs : div_sz;
    uint64_t phase0_ind;
    if(wf_id <= mid) phase0_ind = div_sz * wf_id;
    else phase0_ind = bigs + (div_sz * wf_id);

    uint64_t phase2_sz  = phase0_sz;
    uint64_t phase2_ind = wf_id ? phase0_ind : ns - div_sz;
    parallel_pfxsum_work(phase0_ind, phase0_sz, phase2_ind, phase2_sz, wf_id, nt);
  }

public:

  PrefixSum(A* src, B* dst) : src(src), dst(dst), paste(B()), lock() {}

  void computePrefixSum(uint64_t ns)
  {
    galois::on_each([&](unsigned tid, unsigned numThreads)
        {
          this->parallel_pfxsum_op(ns, tid, numThreads);
        });
    this->lock.reset();
  }
  const char* name()
  {
    return typeid(PrefixSum<A,B, transmute, scan_op, combiner, Conduit>).name();
  }
};
