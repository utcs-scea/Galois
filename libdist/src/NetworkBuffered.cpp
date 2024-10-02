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
#include "galois/Threads.h"
#include "galois/concurrentqueue.h"
#include "llvm/Support/CommandLine.h"

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

namespace cll = llvm::cl;

namespace {

cll::opt<int> commCoreID("commCoreID", cll::desc("The core id to pin the communication thread to"), cll::init(0));

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

  static const int COMM_MIN = 1400; //! bytes (sligtly smaller than an ethernet packet)

  unsigned long statSendNum;
  unsigned long statSendBytes;
  unsigned long statRecvNum;
  unsigned long statRecvBytes;
  bool anyReceivedMessages;

  unsigned int numT;

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

  struct message {
      uint32_t host; //!< destination of this message
      uint32_t tag;  //!< tag on message indicating distinct communication phases
      vTy data;      //!< data portion of message

      //! Default constructor initializes host and tag to large numbers.
      message() : host(~0), tag(~0) {}
      //! @param h Host to send message to
      //! @param t Tag to associate with message
      //! @param d Data to save in message
      message(uint32_t h, uint32_t t, vTy&& d)
          : host(h), tag(t), data(std::move(d)) {}

      //! A message is valid if there is data to be sent
      //! @returns true if data is non-empty
      bool valid() const { return !data.empty(); }
  };

  /**
   * Receive buffers for the buffered network interface
   */
  class recvBuffer {
    std::deque<message> data;
    SimpleLock rlock;
    // tag of head of queue
    std::atomic<uint32_t> dataPresent = ~0;

  public:
    std::optional<RecvBuffer> popMsg(uint32_t tag, std::atomic<size_t>& inflightRecvs, uint32_t& src) {
      if (data.empty() || data.front().tag != tag)
        return std::optional<RecvBuffer>();

      src = data.front().host;
      vTy vec(std::move(data.front().data));

      data.pop_front();
      --inflightRecvs;
      if (!data.empty()) {
        dataPresent = data.front().tag;
      } else {
        dataPresent = ~0;
      }

      return std::optional<RecvBuffer>(RecvBuffer(std::move(vec), sizeof(uint32_t)));
    }

    // Worker thread interface
    void add(message m) {
      std::lock_guard<SimpleLock> lg(rlock);
      if (data.empty()) {
        galois::runtime::trace("ADD LATEST ", m.tag);
        dataPresent = m.tag;
      }

      data.push_back(std::move(m));

      assert(data.back().data.size() !=
             (unsigned int)std::count(data.back().data.begin(),
                                      data.back().data.end(), 0));
    }

    bool hasData(uint32_t tag) { return dataPresent == tag; }

    size_t size() { return data.size(); }

    uint32_t getPresentTag() { return dataPresent; }

    bool try_lock() { return rlock.try_lock();}

    void unlock() { return rlock.unlock();}
  }; // end recv buffer class
  
  recvBuffer recvData;

  /**
   * Receive buffers for the buffered network interface
   */
  class concurrentRecvBuffer {
    // single producer multiple consumer
    moodycamel::ConcurrentQueue<message> data;
    unsigned int numT;
    moodycamel::ProducerToken ptok;

  public:
    concurrentRecvBuffer () : ptok(data) {
        numT = galois::getActiveThreads();
    }

    std::optional<RecvBuffer> tryPopMsg(std::atomic<size_t>& inflightRecvs, uint32_t& src) {
      message m;

      if (data.try_dequeue_from_producer(ptok, m)) {
          src = m.host;
          vTy vec(std::move(m.data));
          --inflightRecvs;

          return std::optional<RecvBuffer>(RecvBuffer(std::move(vec), sizeof(uint32_t)));
      }

      return std::optional<RecvBuffer>();
    }

    // Worker thread interface
    void add(message m) {
      data.enqueue(ptok, std::move(m));
    }

    size_t size() { return data.size_approx(); }
  }; // end recv buffer class

  concurrentRecvBuffer recvRemoteWork;
  
  /**
   * Base send buffer class for the buffered network interface
   */
  class sendBuffer {
  protected:
      struct msg {
          uint32_t tag;
          vTy data;
          msg() {}
          msg(uint32_t t, vTy& _data) : tag(t), data(std::move(_data)) {}
      };

      moodycamel::ConcurrentQueue<msg> messages;
      unsigned int numT;

      std::atomic<bool> flush;
      std::atomic<size_t> numBytes;

  public:
      sendBuffer () : flush(false), numBytes(0) {
          numT = galois::getActiveThreads();
      }

      size_t size() { return messages.size_approx(); }

      void setFlush() {
          flush = true;
      }
    
      bool checkFlush() {
          return flush;
      }
    
      virtual void assemble(std::vector<std::pair<uint32_t, vTy>>& payloads, std::atomic<size_t>& GALOIS_UNUSED(inflightSends)) = 0;

      virtual void add(uint32_t tag, vTy& b) = 0;
  }; // end send buffer class

  
  /**
   * Single producer single consumer with multiple tags
   */
  class sendBufferData : public sendBuffer {
      moodycamel::ProducerToken ptok;

  public:
      sendBufferData () : sendBuffer(), ptok(messages) {}
    
      virtual void assemble(std::vector<std::pair<uint32_t, vTy>>& payloads, std::atomic<size_t>& GALOIS_UNUSED(inflightSends)) {
          std::unordered_map<uint32_t, std::tuple<uint32_t, int, vTy>> tagMap;

          msg m;

          while (messages.try_dequeue_from_producer(ptok, m)) {
              union {
                  uint32_t a;
                  uint8_t b[sizeof(uint32_t)];
              } foo;

              if (tagMap.find(m.tag) != tagMap.end()) {
                  uint32_t& len = std::get<0>(tagMap[m.tag]);
                  int& num = std::get<1>(tagMap[m.tag]);
                  vTy& vec = std::get<2>(tagMap[m.tag]);
                  // do not let it go over the integer limit because MPI_Isend cannot deal with it
                  if ((len + num + m.data.size() + sizeof(uint32_t)) > static_cast<size_t>(std::numeric_limits<int>::max())) {
                      // first put onto payloads
                      payloads.push_back(std::make_pair(m.tag, std::move(vec)));
                      ++inflightSends;
                      // then reset values
                      len = m.data.size();
                      num = sizeof(uint32_t);
                      vTy vec_temp;
                    
                      foo.a = m.data.size();
                      vec_temp.insert(vec_temp.end(), &foo.b[0], &foo.b[sizeof(uint32_t)]);
                      vec_temp.insert(vec_temp.end(), m.data.begin(), m.data.end());

                      vec = std::move(vec_temp);
                  }
                  else {
                      len += m.data.size();
                      num += sizeof(uint32_t);

                      foo.a = m.data.size();
                      vec.insert(vec.end(), &foo.b[0], &foo.b[sizeof(uint32_t)]);
                      vec.insert(vec.end(), m.data.begin(), m.data.end());
                  }
              } else {
                  uint32_t len = m.data.size();
                  int num = sizeof(uint32_t);
                  vTy vec;
                
                  foo.a = m.data.size();
                  vec.insert(vec.end(), &foo.b[0], &foo.b[sizeof(uint32_t)]);
                  vec.insert(vec.end(), m.data.begin(), m.data.end());

                  tagMap[m.tag] = std::make_tuple(len, num, std::move(vec));
              }
            
              --inflightSends;
              numBytes -= m.data.size();
          }

          flush = false;

          // push all payloads in the map into the vector
          for (auto it=tagMap.begin(); it!=tagMap.end(); ++it) {
              payloads.push_back(std::make_pair(it->first, std::move(std::get<2>(it->second))));
              ++inflightSends;
          }

          // sort the payloads with respect to tag
          // lower tag value should go first
          std::sort(payloads.begin(), payloads.end(),
                    [](const auto& a, const auto& b) {
                        return a.first < b.first;
                    }
          );
      }

      virtual void add(uint32_t tag, vTy& b) {
          unsigned oldNumBytes = numBytes;
          numBytes += b.size();
          galois::runtime::trace("BufferedAdd", oldNumBytes, numBytes, tag, galois::runtime::printVec(b));
          messages.enqueue(ptok, msg(tag, b));

          if (numBytes >= COMM_MIN) {
              flush = true;
          }
      }
  };

  std::vector<sendBufferData> sendData;

  /**
   * multiple producer single consumer with single tag
   */
  class sendBufferRemoteWork : public sendBuffer {
      std::vector<moodycamel::ProducerToken> ptok;
      moodycamel::ConsumerToken ctok;

  public:
      sendBufferRemoteWork () : sendBuffer(), ctok(messages) {
          for (unsigned int t=0; t<numT; t++) {
              ptok.emplace_back(messages);
          }
      }
    
      virtual void assemble(std::vector<std::pair<uint32_t, vTy>>& payloads, std::atomic<size_t>& GALOIS_UNUSED(inflightSends)) {
          uint32_t tag = galois::runtime::remoteWorkTag;
          msg m;

          uint32_t len = 0;
          int num = 0;
          vTy vec;
          union {
              uint32_t a;
              uint8_t b[sizeof(uint32_t)];
          } foo;

          while (messages.try_dequeue(ctok, m)) {
              // do not let it go over the integer limit because MPI_Isend cannot deal with it
              if ((len + num + m.data.size() + sizeof(uint32_t)) > static_cast<size_t>(std::numeric_limits<int>::max())) {
                  // first put onto payloads
                  payloads.push_back(std::make_pair(tag, std::move(vec)));
                  ++inflightSends;
                  // then reset values
                  len = m.data.size();
                  num = sizeof(uint32_t);
                  vTy vec_temp;
                
                  foo.a = m.data.size();
                  vec_temp.insert(vec_temp.end(), &foo.b[0], &foo.b[sizeof(uint32_t)]);
                  vec_temp.insert(vec_temp.end(), m.data.begin(), m.data.end());

                  vec = std::move(vec_temp);
              }
              else {
                  len += m.data.size();
                  num += sizeof(uint32_t);

                  foo.a = m.data.size();
                  vec.insert(vec.end(), &foo.b[0], &foo.b[sizeof(uint32_t)]);
                  vec.insert(vec.end(), m.data.begin(), m.data.end());
              }
            
              --inflightSends;
              numBytes -= m.data.size();
          }

          flush = false;

          if (vec.size() != 0) {
              // push remaining payload into the vector
              payloads.push_back(std::make_pair(tag, std::move(vec)));
              ++inflightSends;
          }
      }

      void add(uint32_t tag, vTy& b) {
          unsigned oldNumBytes = numBytes;
          numBytes += b.size();
          galois::runtime::trace("BufferedAdd", oldNumBytes, numBytes, tag, galois::runtime::printVec(b));
          unsigned tid = galois::substrate::ThreadPool::getTID();
          messages.enqueue(ptok[tid], msg(tag, b));

          if (numBytes >= COMM_MIN) {
              flush = true;
          }
      }
  };

  std::vector<sendBufferRemoteWork> sendRemoteWork;

  /**
   * Message type to send/recv in this network IO layer.
   */
  struct mpiMessage {
      uint32_t host;
      uint32_t tag;
      vTy data;
      MPI_Request req;
        
      mpiMessage(uint32_t host, uint32_t tag, vTy&& data) : host(host), tag(tag), data(std::move(data)) {}
      mpiMessage(uint32_t host, uint32_t tag, size_t len) : host(host), tag(tag), data(len) {}
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
              sendInflight.pop_front();
              --inflightSends;
          }
      }
  }

  void send(message m) {
      sendInflight.emplace_back(m.host, m.tag, std::move(m.data));
      auto& f = sendInflight.back();
      galois::runtime::trace("MPI SEND", f.host, f.tag, f.data.size(), galois::runtime::printVec(f.data));
      int rv = MPI_Isend(f.data.data(), f.data.size(), MPI_BYTE, f.host, f.tag, MPI_COMM_WORLD, &f.req);
      handleError(rv);
  }
  
  std::deque<mpiMessage> recvInflight;

  uint32_t getSubMessageLen(vTy& data_array, size_t offset) {
      if ((data_array.size() - offset) > sizeof(uint32_t)) {
          union {
              uint8_t a[sizeof(uint32_t)];
              uint32_t b;
          } c;

          for (size_t i=0; i<sizeof(uint32_t); i++) {
              c.a[i] = data_array[offset + i];
          }

          return c.b;
      } else {
          return ~0;
      }
  }

  std::optional<vTy> getSubMessage(vTy& data_array, size_t& offset, uint32_t len) {
      if ((data_array.size() - offset) > len) {
          vTy vec;

          vec.insert(vec.end(), data_array.begin() + offset, data_array.begin() + offset + sizeof(uint32_t));
          offset += sizeof(uint32_t);
          vec.insert(vec.end(), data_array.begin() + offset, data_array.begin() + offset + len);
          offset += len;

          return std::optional<vTy>(std::move(vec));
      } else {
          return std::optional<vTy>();
      }
  }

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
              int nbytes;
              rv = MPI_Get_count(&status, MPI_BYTE, &nbytes);
              handleError(rv);
              recvInflight.emplace_back(status.MPI_SOURCE, status.MPI_TAG, nbytes);
              auto& m = recvInflight.back();
              memUsageTracker.incrementMemUsage(m.data.size());
              rv = MPI_Irecv(m.data.data(), nbytes, MPI_BYTE, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, &m.req);
              handleError(rv);
              galois::runtime::trace("MPI IRECV", status.MPI_SOURCE, status.MPI_TAG, m.data.size());
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
              if (m.tag == galois::runtime::terminationTag) {
                  hostTermination[m.host] = true;
                  //galois::gPrint("Host ", ID, " : received termination from Host ", m.host, "\n");
              }
              else {
                  //galois::gPrint("Host ", ID, " : MPI_recv from Host ", m.host, "\n");
                  assert(m.data.size() != (unsigned int)std::count(m.data.begin(), m.data.end(), 0));
                  galois::runtime::trace("BufferedRecieving", m.host, m.tag, galois::runtime::printVec(m.data));

                  // Disassemble the aggregated message
                  size_t offset = 0;
                  while (offset < m.data.size()) {
                      // Read the length of the next message
                      uint32_t len = getSubMessageLen(m.data, offset);
                      if (len == ~0U || len == 0) {
                          galois::gError("Cannot read the length of the received message!\n");
                          break; // Not enough data to read the message size
                      }

                      message sub_msg;
                      sub_msg.host = m.host;
                      sub_msg.tag = m.tag;
                      
                      // Read the message data
                      auto vec = getSubMessage(m.data, offset, len);
                      if (vec.has_value()) {
                          sub_msg.data = std::move(vec.value());
                      }
                      else {
                          galois::gError("Cannot read the entire received message (length mismatch)!\n");
                          break; // Not enough data for the complete message
                      }

                      // Add the disassembled message to the receive buffer
                      if (sub_msg.tag == galois::runtime::remoteWorkTag) {
                          recvRemoteWork.add(std::move(sub_msg));
                      }
                      else {
                          recvData.add(std::move(sub_msg));
                      }
                      
                      ++inflightRecvs;
                  }
              }
                
              recvInflight.pop_front();
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

      int core_id = sched_getcpu(); // Get the current CPU core
      galois::gPrint("Host ", ID, " : workerThread running on core ", core_id, "\n");

      ready = 1;
      while (ready < 2) { /*fprintf(stderr, "[WaitOnReady-2]");*/
      };
      while (ready != 3) {
          for (unsigned i = 0; i < Num; ++i) {
              // push progress forward on the network IO
              sendComplete();
              recvProbe();
          
              if (i != ID) {
                  // handle send queue
                  // 1. remote work
                  auto& srw = sendRemoteWork[i];
                  if (srw.checkFlush()) {
                      std::vector<std::pair<uint32_t, vTy>> payloads;
                      srw.assemble(payloads, inflightSends);
                      
                      for (size_t k=0; k<payloads.size(); k++) {
                          message msg;
                          msg.host = i;
                          msg.tag = payloads[k].first;
                          msg.data = std::move(payloads[k].second);

                          if (msg.tag != ~0U) {
                              galois::runtime::trace("BufferedSending", msg.host, msg.tag, galois::runtime::printVec(msg.data));
                              memUsageTracker.incrementMemUsage(msg.data.size());
                              send(std::move(msg));
                              //galois::gPrint("Host ", ID, " : MPI_Send work to Host ", i, "\n");
                          }
                      }

                      
                      // 2. termination
                      if (sendTermination[i]) {
                          ++inflightSends;
                          message msg;
                          msg.host = i;
                          msg.tag = galois::runtime::terminationTag;
                          
                          send(std::move(msg));
                          //galois::gPrint("Host ", ID, " : MPI_Send termination to Host ", i, "\n");

                          sendTermination[i] = false;
                      }
                  }
                  // 3. data
                  auto& sd = sendData[i];
                  if (sd.checkFlush()) {
                      std::vector<std::pair<uint32_t, vTy>> payloads;
                      sd.assemble(payloads, inflightSends);
                      
                      for (size_t k=0; k<payloads.size(); k++) {
                          message msg;
                          msg.host = i;
                          msg.tag = payloads[k].first;
                          msg.data = std::move(payloads[k].second);

                          if (msg.tag != ~0U) {
                              galois::runtime::trace("BufferedSending", msg.host, msg.tag, galois::runtime::printVec(msg.data));
                              memUsageTracker.incrementMemUsage(msg.data.size());
                              send(std::move(msg));
                              //galois::gPrint("Host ", ID, " : MPI_Send data to Host ", i, "\n");
                          }
                      }
                  }
              }
          }
      }
      
      finalizeMPI();
  }
  
  std::thread worker;
  std::atomic<int> ready;
  
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
  NetworkInterfaceBuffered() {
    inflightSends       = 0;
    inflightRecvs       = 0;
    ready               = 0;
    anyReceivedMessages = false;
    worker = std::thread(&NetworkInterfaceBuffered::workerThread, this);
    numT = galois::getActiveThreads();
    while (ready != 1) {};
    
    sendData = decltype(sendData)(Num);
    sendRemoteWork = decltype(sendRemoteWork)(Num);
    sendTermination = decltype(sendTermination)(Num);
    hostTermination = decltype(hostTermination)(Num);
    resetTermination();
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
    ++inflightSends;
    tag += phase;
    statSendNum += 1;
    statSendBytes += buf.size();
    galois::runtime::trace("sendTagged", dest, tag,
                           galois::runtime::printVec(buf.getVec()));
    
    auto& sd = sendData[dest];
    sd.add(tag, buf.getVec());
  }
  
  virtual void sendWork(uint32_t dest, SendBuffer& buf) {
    ++inflightSends;
    statSendNum += 1;
    statSendBytes += buf.size();
    galois::runtime::trace("sendRemoteWork", dest, galois::runtime::printVec(buf.getVec()));
    
    auto& sd = sendRemoteWork[dest];
    sd.add(galois::runtime::remoteWorkTag, buf.getVec());
  }

  virtual std::optional<std::pair<uint32_t, RecvBuffer>>
  receiveTagged(uint32_t tag, int phase) {
      tag += phase;

      if (recvData.hasData(tag)) {
          if (recvData.try_lock()) {
              uint32_t src;
              auto buf = recvData.popMsg(tag, inflightRecvs, src);
              recvData.unlock();
              if (buf) {
                  ++statRecvNum;
                  statRecvBytes += buf->size();
                  memUsageTracker.decrementMemUsage(buf->size());
                  galois::runtime::trace("recvTagged", src, tag, galois::runtime::printVec(buf->getVec()));
                  anyReceivedMessages = true;
                  return std::optional<std::pair<uint32_t, RecvBuffer>>(std::make_pair(src, std::move(*buf)));
              }
          }
      }

      galois::runtime::trace("recvTagged BLOCKED this by that", tag, recvData.getPresentTag());
      return std::optional<std::pair<uint32_t, RecvBuffer>>();
  }
  
  virtual std::optional<std::pair<uint32_t, RecvBuffer>>
  receiveRemoteWork(bool& terminateFlag) {
      terminateFlag = false;

      uint32_t src;
      auto buf = recvRemoteWork.tryPopMsg(inflightRecvs, src);
      if (buf) {
          ++statRecvNum;
          statRecvBytes += buf->size();
          memUsageTracker.decrementMemUsage(buf->size());
          galois::runtime::trace("recvRemoteWork", src, galois::runtime::printVec(buf->getVec()));
          anyReceivedMessages = true;
          return std::optional<std::pair<uint32_t, RecvBuffer>>(std::make_pair(src, std::move(*buf)));
      }
      else {
          if (checkTermination()) {
              terminateFlag = true;
              return std::optional<std::pair<uint32_t, RecvBuffer>>();
          }
          else {
              galois::runtime::trace("recvRemoteWork queue empty but not termination");
              return std::optional<std::pair<uint32_t, RecvBuffer>>();
          }
      }
  }
  
  virtual void flush() {
    for (auto& sd : sendData) {
        sd.setFlush();
    }
  }
  
  virtual void flushData() {
    for (auto& sd : sendData) {
        sd.setFlush();
    }
  }
  
  virtual void flushRemoteWork() {
    for (auto& sd : sendRemoteWork) {
        sd.setFlush();
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
      return (inflightSends > 0);
  }

  virtual bool anyPendingReceives() {
    if (anyReceivedMessages) { // might not be acted on by the computation yet
      anyReceivedMessages = false;
      // galois::gDebug("[", ID, "] receive out of buffer \n");
      return true;
    }
    // if (inflightRecvs > 0) {
    // galois::gDebug("[", ID, "] inflight receive: ", inflightRecvs, " \n");
    // }
    return (inflightRecvs > 0);
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
