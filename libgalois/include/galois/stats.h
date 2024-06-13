// SPDX-License-Identifier: BSD-2-Clause
// Copyright (c) 2023. University of Texas at Austin. All rights reserved.

#pragma once

#include <sys/resource.h>
#include <linux/perf_event.h>
#include <sys/ioctl.h>
#include <unistd.h>
#include <asm/unistd.h>
#include <cstdint>
#include <chrono>
#include <iostream>
#include <string>
#include <fstream>

class HighResTimer {
private:
  std::chrono::high_resolution_clock::time_point start_time;
  std::chrono::high_resolution_clock::time_point end_time;

public:
  HighResTimer() = default;

  void start() { start_time = std::chrono::high_resolution_clock::now(); }

  void stop() { end_time = std::chrono::high_resolution_clock::now(); }

  uint64_t getDurationNano() const {
    return std::chrono::duration_cast<std::chrono::nanoseconds>(end_time -
                                                                start_time)
        .count();
  }
};

class PerfEvent {
private:
  int fd;
  struct perf_event_attr pe;

  uint64_t perf_event_open(struct perf_event_attr* hw_event, pid_t pid, int cpu,
                           int group_fd, uint64_t flags) {
    return syscall(__NR_perf_event_open, hw_event, pid, cpu, group_fd, flags);
  }

public:
  PerfEvent(uint64_t type, uint64_t config) {
    memset(&pe, 0, sizeof(struct perf_event_attr));
    pe.type           = type;
    pe.size           = sizeof(struct perf_event_attr);
    pe.config         = config;
    pe.disabled       = 1;
    pe.exclude_kernel = 1; // Exclude kernel events
    pe.exclude_hv     = 1; // Exclude hypervisor events

    fd = perf_event_open(&pe, 0, -1, -1, 0);
    if (fd == -1) {
      std::cerr << "Error opening leader " << std::endl;
    }
  }

  void start() {
    if (ioctl(fd, PERF_EVENT_IOC_RESET, 0) == -1) {
      std::cerr << "Error in PERF_EVENT_IOC_RESET" << std::endl;
    }
    if (ioctl(fd, PERF_EVENT_IOC_ENABLE, 0) == -1) {
      std::cerr << "Error in PERF_EVENT_IOC_ENABLE" << std::endl;
    }
  }

  void stop() {
    if (ioctl(fd, PERF_EVENT_IOC_DISABLE, 0) == -1) {
      std::cerr << "Error in PERF_EVENT_IOC_DISABLE" << std::endl;
    }
  }

  uint64_t readValue() {
    uint64_t count;
    if (read(fd, &count, sizeof(uint64_t)) == -1) {
      std::cerr << "Error reading count" << std::endl;
    }
    return count;
  }

  ~PerfEvent() { close(fd); }
};

class ScopeBenchmarker {
private:
  std::string scopeName;
  std::string outputName;
  HighResTimer timer;
  PerfEvent cacheMissesEvent;
  PerfEvent cacheReferencesEvent;
  PerfEvent instructionsEvent;
  PerfEvent minorPageFaultsEvent;
  PerfEvent majorPageFaultsEvent;

  static uint64_t getMaxRSS() {
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
#if defined(__linux__)
    return usage.ru_maxrss;
#else
    return usage.ru_maxrss * (sysconf(_SC_PAGESIZE) / 1024);
#endif
  }

public:
  explicit ScopeBenchmarker(const std::string& _scopeName, const std::string& _outputName)
      : scopeName(_scopeName),
        outputName(_outputName),
        cacheMissesEvent(PERF_TYPE_HARDWARE, PERF_COUNT_HW_CACHE_MISSES),
        cacheReferencesEvent(PERF_TYPE_HARDWARE,
                             PERF_COUNT_HW_CACHE_REFERENCES),
        instructionsEvent(PERF_TYPE_HARDWARE, PERF_COUNT_HW_INSTRUCTIONS),
        minorPageFaultsEvent(PERF_TYPE_SOFTWARE, PERF_COUNT_SW_PAGE_FAULTS_MIN),
        majorPageFaultsEvent(PERF_TYPE_SOFTWARE,
                             PERF_COUNT_SW_PAGE_FAULTS_MAJ) {
    timer.start();
    cacheMissesEvent.start();
    cacheReferencesEvent.start();
    instructionsEvent.start();
    minorPageFaultsEvent.start();
    majorPageFaultsEvent.start();
  }

  ~ScopeBenchmarker() {
    timer.stop();
    cacheMissesEvent.stop();
    cacheReferencesEvent.stop();
    instructionsEvent.stop();
    minorPageFaultsEvent.stop();
    majorPageFaultsEvent.stop();

    uint64_t cacheMisses     = cacheMissesEvent.readValue();
    uint64_t cacheReferences = cacheReferencesEvent.readValue();
    uint64_t instructions    = instructionsEvent.readValue();
    uint64_t max_rss         = getMaxRSS();
    uint64_t minorPageFaults = minorPageFaultsEvent.readValue();
    uint64_t majorPageFaults = majorPageFaultsEvent.readValue();


    std::string output = "Benchmark results for " + scopeName + ":\n" +
                         "Duration: " + std::to_string(timer.getDurationNano()) +
                         " nanoseconds\n" +
                         "Max RSS: " + std::to_string(max_rss) + " KB\n" +
                         "Cache Misses: " + std::to_string(cacheMisses) + "\n" +
                         "Cache References: " + std::to_string(cacheReferences) + "\n" +
                         "Instructions: " + std::to_string(instructions) + "\n" +
                         "Minor Page Faults: " + std::to_string(minorPageFaults) + "\n" +
                         "Major Page Faults: " + std::to_string(majorPageFaults) + "\n";

    std::cout << output;

    std::ofstream fp(outputName, std::ios::out);
    fp << output;
    fp.close();
  }
};

#define BENCHMARK_SCOPE(name, fname) ScopeBenchmarker benchmarker##__LINE__(name, fname)
