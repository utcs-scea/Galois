#include "GNNBench/Start.h"

namespace cll = llvm::cl;

cll::opt<unsigned> num_threads("t", cll::desc("Number of threads (default 1)"),
                               cll::init(1));
cll::opt<unsigned> num_epochs("epochs",
                              cll::desc("Number of epochs (default 50)"),
                              cll::init(50));

#ifdef GALOIS_ENABLE_GPU
std::string personality_str(DevicePersonality p) {
  switch (p) {
  case DevicePersonality::CPU:
    return "CPU";
  case DevicePersonality::GPU_CUDA:
    return "GPU_CUDA";
  default:
    GALOIS_LOG_ASSERT(false && "Invalid personality");
    break;
  }
  return "";
}

cll::opt<int> num_nodes(
    "numNodes",
    cll::desc("Num of physical nodes with devices (default = num of hosts): "
              "detect GPU to use for each host automatically"),
    cll::init(-1));
cll::opt<std::string> personality_set(
    "pset",
    cll::desc("String specifying personality for hosts on each physical "
              "node. 'c'=CPU, 'g'=GPU (default 'c')"),
    cll::init("c"));
#endif

cll::opt<std::string>
    stat_file("statFile", cll::desc("Optional output file to print stats to"));

////////////////////////////////////////////////////////////////////////////////

static void PrintVersion(llvm::raw_ostream& out) {
  out << "D-Galois Benchmark Suite v" << galois::getVersion() << " ("
      << galois::getRevision() << ")\n";
  out.flush();
}

////////////////////////////////////////////////////////////////////////////////

void GNNBenchStart(int argc, char** argv, const char* app) {
  GNNBenchStart(argc, argv, app, nullptr, nullptr);
}

void GNNBenchStart(int argc, char** argv, const char* app, const char* desc,
                   const char* url) {
  llvm::cl::SetVersionPrinter(PrintVersion);
  llvm::cl::ParseCommandLineOptions(argc, argv);
  num_threads = galois::setActiveThreads(num_threads);
  galois::runtime::setStatFile(stat_file);

  auto& net = galois::runtime::getSystemNetworkInterface();

  if (net.ID == 0) {
    PrintVersion(llvm::outs());
    llvm::outs() << "Copyright (C) " << galois::getCopyrightYear()
                 << " The University of Texas at Austin\n";
    llvm::outs() << "http://iss.ices.utexas.edu/galois/\n\n";
    llvm::outs() << "application: " << (app ? app : "unspecified") << "\n";

    if (desc) {
      llvm::outs() << desc << "\n";
    }
    if (url) {
      llvm::outs()
          << "http://iss.ices.utexas.edu/?p=projects/galois/benchmarks/" << url
          << "\n";
    }
    llvm::outs() << "\n";
    llvm::outs().flush();

    std::ostringstream cmdout;

    for (int i = 0; i < argc; ++i) {
      cmdout << argv[i];
      if (i != argc - 1)
        cmdout << " ";
    }

    galois::runtime::reportParam("GNNBench", "CommandLine", cmdout.str());
    galois::runtime::reportParam("GNNBench", "Threads", num_threads);
    galois::runtime::reportParam("GNNBench", "Hosts", net.Num);
    galois::runtime::reportParam("GNNBench", "Run_UUID",
                                 galois::runtime::getRandUUID());
    galois::runtime::reportParam("GNNBench", "InputDirectory", input_directory);
    galois::runtime::reportParam("GNNBench", "Input", input_name);
    galois::runtime::reportParam("GNNBench", "PartitionScheme",
                                 GNNPartitionToString(partition_scheme));
    // XXX report the rest of the command line options
  }

  char name[256];
  gethostname(name, 256);
  galois::runtime::reportParam("GNNBench", "Hostname", name);

#ifdef GALOIS_ENABLE_GPU
  internal::heteroSetup();
#endif
}

#ifdef GALOIS_ENABLE_GPU
void internal::heteroSetup() {
  const unsigned my_host_id = galois::runtime::getHostID();

  auto& net = galois::runtime::getSystemNetworkInterface();

  if (num_nodes == -1) {
    num_nodes = net.Num;
  }

  GALOIS_LOG_ASSERT((net.Num % num_nodes) == 0);

  device_personality = DevicePersonality::CPU;
  if (personality_set.length() == (net.Num / num_nodes)) {
    switch (personality_set.c_str()[my_host_id % (net.Num / num_nodes)]) {
    case 'g':
      galois::gInfo(my_host_id, " chooses GPU");
      device_personality = DevicePersonality::GPU_CUDA;
      break;
    case 'c':
      galois::gInfo(my_host_id, " chooses CPU");
      device_personality = DevicePersonality::CPU;
      break;
    }

    if (device_personality == DevicePersonality::GPU_CUDA) {
      gpudevice = get_gpu_device_id(personality_set, num_nodes);
    } else {
      gpudevice = -1;
    }

    SetCUDADeviceId(gpudevice);
  } else {
    galois::gWarn(
        "Command line option -pset ignored because its string length is not "
        "equal to the number of processes/hosts on each physical node");
  }
}
#endif
