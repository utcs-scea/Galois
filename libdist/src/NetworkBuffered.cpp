/*
 * This file belongs to the Galois project, a C++ library for exploiting
 * parallelism. The code is being released under the terms of the 3-Clause BSD
 * License (a copy is located in LICENSE.txt at the top-level directory).
 *
 * Copyright (C) 2018, The University of Texas at Austin. All rights reserved.
 * UNIVERSITY EXPRESSLY DISCLAIMS ANY AND ALL WARRANTIES CONCERNING THIS
 * SOFTWARE AND DOCUMENTATION, INCLUDING ANY WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR ANY PARTICULAR PURPOSE, NON-INFRINGEMENT AND WARRANTIES OF
 * PERFORMANCE, AND ANY WARRANTY THAT MIGHT OTHERWISE ARISE FROM COURSE OF
 * DEALING OR USAGE OF TRADE.  NO WARRANTY IS EITHER EXPRESS OR IMPLIED WITH
 * RESPECT TO THE USE OF THE SOFTWARE OR DOCUMENTATION. Under no circumstances
 * shall University be liable for incidental, special, indirect, direct or
 * consequential damages or loss of profits, interruption of business, or
 * related expenses which may arise from use of Software or Documentation,
 * including but not limited to those resulting from defects in Software and/or
 * Documentation, or loss or inaccuracy of data of any kind.
 */

/**
 * @file NetworkBuffered.cpp
 *
 * Contains NetworkInterfaceBuffered, an implementation of a network interface
 * that buffers messages before sending them out.
 *
 * @todo document this file more
 */

#include "galois/runtime/Network.h"
#include "galois/runtime/Tracer.h"
#include "galois/runtime/ThreadTimer.h"
#include "galois/Threads.h"
#include "galois/concurrentqueue.h"

#include <thread>
#include <mutex>
#include <iostream>
#include <limits>
#include <unordered_map>
#include <tuple>
#include <functional>
#include <pthread.h>
#include <sched.h>
#include <unistd.h>

using namespace galois::runtime;
using namespace galois::substrate;

constexpr static const char* const NETWORK_NAME = "Network";

namespace {

/**
 * @class NetworkInterfaceBuffered
 *
 * Buffered network interface: messages are buffered before they are sent out.
 * A single worker thread is initialized to send/receive messages from/to
 * buffers.
 */
class NetworkInterfaceBuffered : public NetworkInterface {
  using NetworkInterface::ID;
  using NetworkInterface::Num;

  constexpr size_t COMM_MIN = 1400; //! bytes (sligtly smaller than an ethernet packet)

  unsigned long statSendNum;
  unsigned long statSendBytes;
  unsigned long statRecvNum;
  unsigned long statRecvBytes;
  bool anyReceivedMessages;

  unsigned numT;

  // using vTy = std::vector<uint8_t>;
  using vTy = galois::PODResizeableArray<uint8_t>;

  /**
   * Wrapper for dealing with MPI error codes. Program dies if the error code
   * isn't MPI_SUCCESS.
   *
   * @param rc Error code to check for success
   */
  static void handleError(int rc) {
      if (rc != MPI_SUCCESS) {
          MPI_Abort(MPI_COMM_WORLD, rc);
      }
  }

  /**
   * Get the host id of the caller.
   *
   * @returns host id of the caller with regard to the MPI setup
   */
  int getID() {
      int rank;
      handleError(MPI_Comm_rank(MPI_COMM_WORLD, &rank));
      return rank;
  }

  /**
   * Get the total number of hosts in the system.
   *
   * @returns number of hosts with regard to the MPI setup
   */
  int getNum() {
      int hostSize;
      handleError(MPI_Comm_size(MPI_COMM_WORLD, &hostSize));
      return hostSize;
  }

  struct recvMessage {
      uint32_t host; //!< destination of this message
      uint32_t tag;  //!< tag on message indicating distinct communication phases
      vTy data;      //!< data portion of message

      //! Default constructor initializes host and tag to large numbers.
      recvMessage() : host(~0), tag(~0) {}
      //! @param h Host to send message to
      //! @param t Tag to associate with message
      //! @param d Data to save in message
      recvMessage(uint32_t h, uint32_t t, vTy&& d) : host(h), tag(t), data(std::move(d)) {}
  };

  /**
   * Receive buffers for the buffered network interface
   */
  class recvBufferData {
      // single producer single consumer
      moodycamel::ConcurrentQueue<recvMessage> data;
      moodycamel::ProducerToken ptok;

      std::atomic<uint32_t> frontTag;
      recvMessage frontMsg;

      std::string pop_str;
      galois::StatTimer StatTimer_pop;
      std::string add_str;
      galois::StatTimer StatTimer_add;

  public:
      std::atomic<size_t> inflightRecvs = 0;

      recvBufferData() : ptok(data), frontTag(~0U), pop_str("RecvDataTimer_Pop"), StatTimer_pop(pop_str.c_str(), NETWORK_NAME), add_str("RecvDataTimer_Add"), StatTimer_add(add_str.c_str(), NETWORK_NAME) {}

      std::optional<RecvBuffer> tryPopMsg(uint32_t tag, uint32_t& src) {
          StatTimer_pop.start();

          if (frontTag == ~0U) { // no messages available
              StatTimer_pop.stop();
              return std::optional<RecvBuffer>();
          }
          else {
              if (frontTag != tag) {
                  StatTimer_pop.stop();
                  return std::optional<RecvBuffer>();
              }
              else {
                  src = frontMsg.host;
                  frontTag = ~0U;

                  --inflightRecvs;
                  StatTimer_pop.stop();
                  return std::optional<RecvBuffer>(RecvBuffer(std::move(frontMsg.data)));
              }
          }
      }

      // Worker thread interface
      void add(uint32_t host, uint32_t tag, vTy vec) {
          StatTimer_add.start();
          data.enqueue(ptok, recvMessage(host, tag, std::move(vec)));
          StatTimer_add.stop();
      }

      bool hasData(uint32_t tag) {
          if (frontTag == ~0U) {
              if (data.size_approx() != 0) {
                  bool success = data.try_dequeue_from_producer(ptok, frontMsg);
                  if (success) {
                      frontTag = frontMsg.tag;
                  }
              }
          }

          return frontTag == tag;
      }
  }; // end recv buffer class

  std::vector<recvBufferData> recvData;

  /**
   * Receive buffers for the buffered network interface
   */
  class recvBufferRemoteWork {
      // single producer single consumer
      moodycamel::ConcurrentQueue<vTy> data;
      moodycamel::ProducerToken ptok;

      std::string pop_str;
      galois::StatTimer StatTimer_pop;
      std::string add_str;
      galois::StatTimer StatTimer_add;

  public:
      std::atomic<size_t> inflightRecvs = 0;

      recvBufferRemoteWork() : ptok(data), pop_str("RecvWorkTimer_Pop"), StatTimer_pop(pop_str.c_str(), NETWORK_NAME), add_str("RecvWorkTimer_Add"), StatTimer_add(add_str.c_str(), NETWORK_NAME) {}

      std::optional<RecvBuffer> tryPopMsg() {
          StatTimer_pop.start();

          vTy vec;
          bool success = data.try_dequeue_from_producer(ptok, vec);
          if (success) {
              --inflightRecvs;
              StatTimer_pop.stop();
              return std::optional<RecvBuffer>(RecvBuffer(std::move(vec)));
          }
          else {
              return std::optional<RecvBuffer>();
          }
      }

      // Worker thread interface
      void add(vTy vec) {
          StatTimer_add.start();
          data.enqueue(ptok, std::move(vec));
          StatTimer_add.stop();
      }
  }; // end recv buffer class

  recvBufferRemoteWork recvRemoteWork;

  struct sendMessage {
      uint32_t tag;  //!< tag on message indicating distinct communication phases
      vTy data;      //!< data portion of message

      //! Default constructor initializes host and tag to large numbers.
      sendMessage() : tag(~0) {}
      //! @param t Tag to associate with message
      //! @param d Data to save in message
      sendMessage(uint32_t t) : tag(t) {}
      sendMessage(uint32_t t, vTy&& d) : tag(t), data(std::move(d)) {}
  };

  /**
   * Single producer single consumer with multiple tags
   */
  class sendBufferData {
      moodycamel::ConcurrentQueue<sendMessage> messages;
      moodycamel::ProducerToken ptok;

      std::atomic<size_t> flush;

      std::string pop_str;
      galois::StatTimer StatTimer_pop;
      std::string add_str;
      galois::StatTimer StatTimer_add;

  public:
      std::atomic<size_t> inflightSends = 0;

      sendBufferData() : ptok(messages), flush(0), pop_str("SendDataTimer_Pop"), StatTimer_pop(pop_str.c_str(), NETWORK_NAME), add_str("SendDataTimer_Add"), StatTimer_add(add_str.c_str(), NETWORK_NAME) {}

      void setFlush() {
          flush += 1;
      }

      bool checkFlush() {
          return flush > 0;
      }

      sendMessage pop() {
          StatTimer_pop.start();

          sendMessage m;
          bool success = messages.try_dequeue_from_producer(ptok, m);
          if (success) {
              flush -= 1;
              StatTimer_pop.stop();
              return m;
          }
          else {
              StatTimer_pop.stop();
              return sendMessage(~0U);
          }
      }

      void push(uint32_t tag, vTy b) {
          StatTimer_add.start();
          messages.enqueue(ptok, sendMessage(tag, std::move(b)));
          ++inflightSends;
          flush += 1;
          StatTimer_add.stop();
      }
  };

  std::vector<sendBufferData> sendData;

  class fastHeap {
    // Single producer multi-consumer
    moodycamel::ConcurrentQueue<uint8_t*> ptrQueue;
    std::vector<std::pair<void*, size_t>> mmaps;
    galois::substrate::SimpleLock allocationLock;
    moodyCamel::ProducerToken ptok;
    uint64_t blockSize;

    explicit fastHeap(uint64_t size) : ptrQueue(), mmaps(), allocationLock(), ptok(ptrQueue), blockSize(size) {
      static_assert(size % COMM_MIN == 0);
      void* ptr = mmap(nullptr, size, PROT_READ |PROT_WRITE, MAP_ANONYMOUS|MAP_PRIVATE|MAP_POPULATE, -1, 0);
      if(ptr == MAP_FAILED) throw "Failed to MMAP region";
      mmaps.emplace_back(ptr, size);
      for(std::uint64_t i = 0; i < size; i+=COMM_MIN) {
        ptrQueue.enqueue(ptok, ptr + i);
      }
    }
    ~fastHeap() {}

    // This will be in multiples of comm_min
    uint8_t* allocate() {
      uint8_t* ret = nullptr;
      ptrQueue.try_deque_from_producer(ptok, ret);
      return ret;
    }
    void deallocate(uint8_t* ptr) {
      ptrQueue.enqueue(ptok, ptr);
    }
  }

  /**
   * single producer single consumer with single tag
   */
  class sendBufferRemoteWork {
      moodycamel::ConcurrentQueue<FixedSizedSerializeBuffer<COMM_MIN,fastHeap>> messages;
      moodycamel::ProducerToken ptok;


      uint32_t len;
      FixedSizeSerializeBuffer<COMM_MIN, fastHeap> vec;

      std::atomic<size_t> flush;

      std::string pop_str;
      galois::StatTimer StatTimer_pop;
      std::string add_str;
      galois::StatTimer StatTimer_add;

  public:
      std::atomic<size_t> inflightSends = 0;

      sendBufferRemoteWork() : ptok(messages), len(0), vec(), flush(0), pop_str("SendWorkTimer_Pop"), StatTimer_pop(pop_str.c_str(), NETWORK_NAME), add_str("SendWorkTimer_Add"), StatTimer_add(add_str.c_str(), NETWORK_NAME) {}

      void setFlush() {
          if (len != 0) {
            messages.enqueue(ptok, std::move(vec));

            len = 0;
            vec.clear();

            ++inflightSends;
            flush += 1;
          }
      }

      bool checkFlush() {
          return flush > 0;
      }

      std::optional<vTy> pop() {
          StatTimer_pop.start();

          vTy m;
          bool success = messages.try_dequeue_from_producer(ptok, m);
          if (success) {
              flush -= 1;
              StatTimer_pop.stop();
              return std::optional<vTy>(std::move(m));
          }
          else {
              StatTimer_pop.stop();
              return std::optional<vTy>();
          }
      }

      void add(vTy b) {
          StatTimer_add.start();

          len += b.size();
          if (len > COMM_MIN) {
              if (len > static_cast<size_t>(std::numeric_limits<int>::max())) {
                  messages.enqueue(ptok, std::move(vec));

                  len = b.size();
                  vec.clear();
                  vec.insert(vec.end(), b.begin(), b.end());
              }
              else {
                  vec.insert(vec.end(), b.begin(), b.end());
                  messages.enqueue(ptok, std::move(vec));

                  len = 0;
                  vec.clear();
              }

              ++inflightSends;
              flush += 1;
          }
          else {
              vec.insert(vec.end(), b.begin(), b.end());
          }

          StatTimer_add.stop();
      }
  };

  std::vector<std::vector<sendBufferRemoteWork>> sendRemoteWork;

  /**
   * Message type to send/recv in this network IO layer.
   */
  struct mpiMessage {
      uint32_t host;
      unsigned tid;
      uint32_t tag;
      vTy data;
      MPI_Request req;

      mpiMessage(uint32_t host, uint32_t tag, vTy&& data) : host(host), tag(tag), data(std::move(data)) {}
      mpiMessage(uint32_t host, unsigned tid, uint32_t tag, vTy&& data) : host(host), tid(tid), tag(tag), data(std::move(data)) {}
  };

  std::deque<mpiMessage> sendInflight;

  void sendComplete() {
      if (!sendInflight.empty()) {
          int flag = 0;
          MPI_Status status;
          auto& f = sendInflight.front();
          int rv  = MPI_Test(&f.req, &flag, &status);
          handleError(rv);
          if (flag) {
              memUsageTracker.decrementMemUsage(f.data.size());
              if (f.tag == galois::runtime::terminationTag) {
                  --inflightTermination;
              }
              else if (f.tag == galois::runtime::remoteWorkTag) {
                --sendRemoteWork[f.host][f.tid].inflightSends;
              }
              else {
                --sendData[f.host].inflightSends;
              }
              sendInflight.pop_front();
          }
      }
  }

  std::string MPISend_str;
  galois::StatTimer StatTimer_MPISend;

  void send(uint32_t dest, unsigned tid, sendMessage m) {
      StatTimer_MPISend.start();
      sendInflight.emplace_back(dest, tid, m.tag, std::move(m.data));
      auto& f = sendInflight.back();
      int rv = MPI_Isend(f.data.data(), f.data.size(), MPI_BYTE, f.host, f.tag, MPI_COMM_WORLD, &f.req);
      handleError(rv);
      StatTimer_MPISend.stop();
  }

  std::deque<mpiMessage> recvInflight;

  std::string MPIRecv_str;
  galois::StatTimer StatTimer_MPIRecv;
  std::string disaggregate_str;
  galois::StatTimer StatTimer_disaggregate;

    // FIXME: Does synchronous recieves overly halt forward progress?
  void recvProbe() {
      int flag = 0;
      MPI_Status status;
      // check for new messages
      int rv = MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);
      handleError(rv);
      if (flag) {
#ifdef GALOIS_USE_BARE_MPI
          assert(status.MPI_TAG <= 32767);
          if (status.MPI_TAG != 32767) {
#endif
              StatTimer_MPIRecv.start();
              int nbytes;
              rv = MPI_Get_count(&status, MPI_BYTE, &nbytes);
              handleError(rv);
              recvInflight.emplace_back(status.MPI_SOURCE, status.MPI_TAG, nbytes);
              auto& m = recvInflight.back();
              memUsageTracker.incrementMemUsage(m.data.size());
              rv = MPI_Irecv(m.data.data(), nbytes, MPI_BYTE, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, &m.req);
              handleError(rv);
              StatTimer_MPIRecv.stop();
#ifdef GALOIS_USE_BARE_MPI
          }
#endif
      }

      // complete messages
      if (!recvInflight.empty()) {
          auto& m  = recvInflight.front();
          int flag = 0;
          rv       = MPI_Test(&m.req, &flag, MPI_STATUS_IGNORE);
          handleError(rv);
          if (flag) {
              StatTimer_disaggregate.start();
              if (m.tag == galois::runtime::terminationTag) {
                  hostTermination[m.host] = true;
                  //galois::gPrint("Host ", ID, " : received termination from Host ", m.host, "\n");
              }
              else if (m.tag == galois::runtime::remoteWorkTag) {
                  //galois::gPrint("Host ", ID, " : MPI_recv work from Host ", m.host, "\n");
                  ++recvRemoteWork.inflightRecvs;
                  recvRemoteWork.add(std::move(m.data));
              }
              else {
                  //galois::gPrint("Host ", ID, " : MPI_recv data from Host ", m.host, "\n");

                  ++recvData[m.host].inflightRecvs;
                  recvData[m.host].add(m.host, m.tag, std::move(m.data));
              }

              recvInflight.pop_front();
              StatTimer_disaggregate.stop();
          }
      }
  }

  void workerThread() {

      // Set thread affinity
      cpu_set_t cpuset;
      CPU_ZERO(&cpuset);           // Clear the CPU set
      CPU_SET(commCoreID, &cpuset);   // Set the specified core

      // Get the native handle of the std::thread
      pthread_t thread = pthread_self();

      // Set the CPU affinity of the thread
      if (pthread_setaffinity_np(thread, sizeof(cpu_set_t), &cpuset) != 0) {
          std::cerr << "Error setting thread affinity" << std::endl;
          return;
      }

      initializeMPI();
      galois::gDebug("[", NetworkInterface::ID, "] MPI initialized");
      ID = getID();
      Num = getNum();

      ready = 1;
      while (ready < 2) { /*fprintf(stderr, "[WaitOnReady-2]");*/
      };
      while (ready != 3) {
          for (unsigned i = 0; i < Num; ++i) {
              if (i != ID) {
                  // handle send queue
                  // 1. remote work
                  bool hostEmpty = true;
                  for (unsigned t=0; t<numT; t++) {
                      // push progress forward on the network IO
                      sendComplete();
                      recvProbe();

                      auto& srw = sendRemoteWork[i][t];
                      if (srw.checkFlush()) {
                          auto payload = srw.pop();
                          //galois::gPrint("Host ", ID, " : flush work to Host ", i, "\n");

                          if (payload.has_value()) {
                              memUsageTracker.incrementMemUsage(payload.value().size());
                              send(i, t, sendMessage(galois::runtime::remoteWorkTag, std::move(payload.value())));
                              //galois::gPrint("Host ", ID, " : MPI_Send work to Host ", i, "\n");

                              hostEmpty = false;
                          }
                      }
                  }

                  if(hostEmpty) { // wait until all works are sent
                      // 2. termination
                      if (sendTermination[i]) {
                          ++inflightTermination;
                          send(i, 0, sendMessage(galois::runtime::terminationTag));
                          //galois::gPrint("Host ", ID, " : MPI_Send termination to Host ", i, "\n");

                          sendTermination[i] = false;
                      }
                  }
                  // 3. data
                  auto& sd = sendData[i];
                  if (sd.checkFlush()) {
                      //galois::gPrint("Host ", ID, " : flush data to Host ", i, "\n");
                      sendMessage msg = sd.pop();

                      if (msg.tag != ~0U) {
                          memUsageTracker.incrementMemUsage(msg.data.size());
                          send(i, 0, std::move(msg));
                          //galois::gPrint("Host ", ID, " : MPI_Send data to Host ", i, "\n");
                      }
                  }
              }
          }
      }

      finalizeMPI();
  }

  std::thread worker;
  std::atomic<int> ready;

  std::atomic<size_t> inflightTermination;
  std::vector<std::atomic<bool>> sendTermination;
  std::vector<std::atomic<bool>> hostTermination;
  virtual void resetTermination() {
      for (unsigned i=0; i<Num; i++) {
          if (i == ID) {
              continue;
          }
          hostTermination[i] = false;
      }
  }

  bool checkTermination() {
      for (unsigned i=0; i<Num; i++) {
          if (i == ID) {
              continue;
          }
          if (hostTermination[i] == false) {
              return false;
          }
      }
      return true;
  }

public:
  NetworkInterfaceBuffered() : MPISend_str("MPISendTimer"), StatTimer_MPISend(MPISend_str.c_str(), NETWORK_NAME), MPIRecv_str("MPIRecvTimer"), StatTimer_MPIRecv(MPIRecv_str.c_str(), NETWORK_NAME), disaggregate_str("DisaggregateTimer"), StatTimer_disaggregate(disaggregate_str.c_str(), NETWORK_NAME) {
    ready               = 0;
    inflightTermination = 0;
    anyReceivedMessages = false;
    worker = std::thread(&NetworkInterfaceBuffered::workerThread, this);
    numT = galois::getActiveThreads();
    while (ready != 1) {};

    recvData = decltype(recvData)(Num);
    sendData = decltype(sendData)(Num);
    sendRemoteWork.resize(Num);
    for (auto& hostSendRemoteWork : sendRemoteWork) {
        std::vector<sendBufferRemoteWork> temp(numT);
        hostSendRemoteWork = std::move(temp);
    }
    sendTermination = decltype(sendTermination)(Num);
    hostTermination = decltype(hostTermination)(Num);
    for (unsigned i=0; i<Num; i++) {
        sendTermination[i] = false;
        if (i == ID) {
            hostTermination[i] = true;
        }
        else {
            hostTermination[i] = false;
        }
    }
    ready    = 2;
  }

  virtual ~NetworkInterfaceBuffered() {
    ready = 3;
    worker.join();
  }

  virtual void sendTagged(uint32_t dest, uint32_t tag, SendBuffer& buf,
                          int phase) {
    tag += phase;
    statSendNum += 1;
    statSendBytes += buf.size();

    auto& sd = sendData[dest];
    sd.push(tag, std::move(buf.getVec()));
  }

  template<typename T>
  virtual void sendWork(uint32_t dest, T t) {
  }

  virtual void sendWork(uint32_t dest, SendBuffer& buf) {
    statSendNum += 1;
    statSendBytes += buf.size();

    unsigned tid = galois::substrate::ThreadPool::getTID();
    auto& sd = sendRemoteWork[dest][tid];
    sd.add(std::move(buf.getVec()));
  }

  virtual std::optional<std::pair<uint32_t, RecvBuffer>>
  receiveTagged(uint32_t tag, int phase) {
      tag += phase;

      for (unsigned h=0; h<Num; h++) {
          if (h == ID) {
              continue;
          }

          auto& rq = recvData[h];
          if (rq.hasData(tag)) {
              uint32_t src;
              auto buf = rq.tryPopMsg(tag, src);
              if (buf) {
                  ++statRecvNum;
                  statRecvBytes += buf->size();
                  memUsageTracker.decrementMemUsage(buf->size());
                  anyReceivedMessages = true;
                  return std::optional<std::pair<uint32_t, RecvBuffer>>(std::make_pair(src, std::move(*buf)));
              }
          }
      }

      return std::optional<std::pair<uint32_t, RecvBuffer>>();
  }

  virtual std::optional<RecvBuffer>
  receiveRemoteWork() {
      auto buf = recvRemoteWork.tryPopMsg();
      if (buf) {
          ++statRecvNum;
          statRecvBytes += buf->size();
          memUsageTracker.decrementMemUsage(buf->size());
          anyReceivedMessages = true;
          return std::optional<RecvBuffer>(std::move(*buf));
      }
      else {
          return std::optional<RecvBuffer>();
      }
  }

  virtual std::optional<RecvBuffer>
  receiveRemoteWork(bool& terminateFlag) {
      terminateFlag = false;

      auto buf = recvRemoteWork.tryPopMsg();
      if (buf) {
          ++statRecvNum;
          statRecvBytes += buf->size();
          memUsageTracker.decrementMemUsage(buf->size());
          anyReceivedMessages = true;
          return std::optional<RecvBuffer>(std::move(*buf));
      }
      else {
          if (checkTermination()) {
              terminateFlag = true;
          }
          return std::optional<RecvBuffer>();
      }
  }

  virtual void flush() {
      flushData();
  }

  virtual void flushData() {
    for (auto& sd : sendData) {
        sd.setFlush();
    }
  }

  virtual void flushRemoteWork() {
    for (auto& hostSendRemoteWork : sendRemoteWork) {
        for (auto& threadSendRemoteWork : hostSendRemoteWork) {
            threadSendRemoteWork.setFlush();
        }
    }
  }

  virtual void broadcastTermination() {
      for (unsigned i=0; i<Num; i++) {
          if (i == ID) {
              continue;
          }
          else {
              sendTermination[i] = true;
          }
      }
  }

  virtual bool anyPendingSends() {
      if (inflightTermination > 0) {
          return true;
      }
      for (unsigned i=0; i<Num; i++) {
          if (sendData[i].inflightSends > 0) {
              return true;
          }
          for (unsigned t=0; t<numT; t++) {
              if (sendRemoteWork[i][t].inflightSends > 0) {
                  return true;
              }
          }
      }

      return false;
  }

  virtual bool anyPendingReceives() {
    if (anyReceivedMessages) { // might not be acted on by the computation yet
      anyReceivedMessages = false;
      // galois::gDebug("[", ID, "] receive out of buffer \n");
      return true;
    }

    if (recvRemoteWork.inflightRecvs > 0) {
        return true;
    }
    for (unsigned i=0; i<Num; i++) {
        if (recvData[i].inflightRecvs > 0) {
            return true;
        }
    }

    return false;
  }

  virtual unsigned long reportSendBytes() const { return statSendBytes; }
  virtual unsigned long reportSendMsgs() const { return statSendNum; }
  virtual unsigned long reportRecvBytes() const { return statRecvBytes; }
  virtual unsigned long reportRecvMsgs() const { return statRecvNum; }

  virtual std::vector<unsigned long> reportExtra() const {
    std::vector<unsigned long> retval(2);
    return retval;
  }

  virtual std::vector<std::pair<std::string, unsigned long>>
  reportExtraNamed() const {
    std::vector<std::pair<std::string, unsigned long>> retval(2);
    return retval;
  }
};

} // namespace

/**
 * Create a buffered network interface, or return one if already
 * created.
 */
NetworkInterface& galois::runtime::makeNetworkBuffered() {
  static std::atomic<NetworkInterfaceBuffered*> net;
  static substrate::SimpleLock m_mutex;

  // create the interface if it doesn't yet exist in the static variable
  auto* tmp = net.load();
  if (tmp == nullptr) {
    std::lock_guard<substrate::SimpleLock> lock(m_mutex);
    tmp = net.load();
    if (tmp == nullptr) {
      tmp = new NetworkInterfaceBuffered();
      net.store(tmp);
    }
  }

  return *tmp;
}
