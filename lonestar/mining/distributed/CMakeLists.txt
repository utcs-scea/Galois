include_directories(${CMAKE_SOURCE_DIR}/libgalois/include)
include_directories("${PROJECT_SOURCE_DIR}/libgluon/include")
include_directories("${PROJECT_SOURCE_DIR}/libcusp/include")
# include_directories("${PROJECT_SOURCE_DIR}/lonestar/mining/libminingbench/include")
include_directories("${PROJECT_SOURCE_DIR}/libpangolin/include")

add_subdirectory(triangle-counting)
