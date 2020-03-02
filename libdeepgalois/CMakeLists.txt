cmake_minimum_required(VERSION 2.8)

SET(OPENBLAS_INC /net/ohm/export/cdgc/cxh/OpenBLAS-faraday/build/include)
SET(OPENBLAS_LIB /net/ohm/export/cdgc/cxh/OpenBLAS-faraday/build/lib)
set(CUB_ROOT "${CMAKE_SOURCE_DIR}/cub") # only required headers
set(MGPU_ROOT "${CMAKE_SOURCE_DIR}/moderngpu") # only required headers
SET(CUDA_INC /org/centers/cdgc/cuda/cuda-8.0/include)
SET(CUDA_LIB /org/centers/cdgc/cuda/cuda-8.0/lib64/)
include_directories(${OPENBLAS_INC})
include_directories(${CMAKE_SOURCE_DIR}/libgalois/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${CUDA_INC})
include_directories(${CMAKE_SOURCE_DIR}/libgpu/include)
include_directories("${CUB_ROOT}")
include_directories("${MGPU_ROOT}/src")
link_directories(${OPENBLAS_LIB})
link_directories(${CMAKE_SOURCE_DIR}/libgalois)

set(USE_CPU ON CACHE BOOL "Build DeepGalois without CUDA support")
if(USE_CPU)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DCPU_ONLY")
else()
  find_package(CUDA REQUIRED)
  set(CUDA_SEPARABLE_COMPILATION ON)
  set(CUDA_PROPAGATE_HOST_FLAGS OFF)
  set(CUDA_HOST_COMPILER g++)
  #list(APPEND CUDA_NVCC_FLAGS "-gencode arch=compute_60,code=sm_60")
  link_directories(${CUDA_LIB})
  link_directories(${CMAKE_SOURCE_DIR}/libgpu)
  set(CUDA_SOURCES
    src/math_functions.cu
	src/aggregator.cu
	src/optimizer.cu
	src/context.cu
	src/node.cu
  )
  cuda_add_library(dg_gpu ${CUDA_SOURCES})
  target_link_libraries(dg_gpu galois_gpu -lcudart -lcublas -lcurand)
  set_target_properties(dg_gpu PROPERTIES COMPILE_FLAGS "-DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_CUDA")
  set_target_properties(dg_gpu PROPERTIES CUDA_SEPERABLE_COMPILATION ON)
  #cuda_compile(MF_O src/math_functions.cu)
  #cuda_compile(AGG_O src/aggregator.cu)
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(sources
  src/layers/graph_conv_layer.cpp
  src/layers/softmax_loss_layer.cpp
  src/math_functions.cpp
  src/aggregator.cpp
  src/optimizer.cpp
  src/context.cpp
  src/node.cpp
  src/net.cpp
)
add_library(dg_cpu STATIC ${sources})

if(USE_CPU)
  target_link_libraries(dg_cpu galois_shmem gllvm)
else()
  target_link_libraries(dg_cpu galois_shmem gllvm galois_gpu)
endif()

target_link_libraries(dg_cpu ${MPI_CXX_LIBRARIES})
target_link_libraries(dg_cpu -lopenblas)
target_link_libraries(dg_cpu -lcudart -lcublas -lcurand)

target_include_directories(dg_cpu PUBLIC
  ${CMAKE_SOURCE_DIR}/libllvm/include
  ${CMAKE_SOURCE_DIR}/libgalois/include
  ${CMAKE_CURRENT_SOURCE_DIR}/include
)

set_target_properties(dg_cpu PROPERTIES
  INTERFACE_POSITION_INDEPENDENT_CODE On
  POSITION_INDEPENDENT_CODE On
)