#pragma once

#include "galois/Galois.h"
#include "galois/Version.h"
#include "GNNBench/Input.h"

////////////////////////////////////////////////////////////////////////////////
// CLI
////////////////////////////////////////////////////////////////////////////////

extern llvm::cl::opt<unsigned> num_threads;
extern llvm::cl::opt<unsigned> num_runs;
extern llvm::cl::opt<unsigned> num_epochs;
extern llvm::cl::opt<std::string> stat_file;

////////////////////////////////////////////////////////////////////////////////
// Init functions
////////////////////////////////////////////////////////////////////////////////

//! Parses command line + setup some stats
void GNNBenchStart(int argc, char** argv, const char* app);
void GNNBenchStart(int argc, char** argv, const char* app, const char* desc,
                   const char* url);