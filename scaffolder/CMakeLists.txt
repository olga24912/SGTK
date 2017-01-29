cmake_minimum_required(VERSION 3.3)
project(scaffolder)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3")

include_directories(src)

add_subdirectory(src)
add_subdirectory(tests)

add_executable(build mainBuildGraph.cpp)
add_executable(filter mainFilter.cpp)

target_link_libraries(build builder_lib contig_graph_lib)
target_link_libraries(filter filter_lib contig_graph_lib)