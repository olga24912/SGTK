cmake_minimum_required(VERSION 3.3)
project(scaffolder)

find_package (ZLIB)
find_package(Threads)

include_directories (${SEQAN_INCLUDE_DIRS})
add_definitions (${SEQAN_DEFINITIONS})
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS} -std=c++14 -O3 -DSEQAN_HAS_ZLIB=1")

include_directories(src)

add_subdirectory(src)
add_subdirectory(tests)

add_executable(build mainBuildGraph.cpp)
add_executable(filter mainFilter.cpp)
add_executable(merger mainMergeContigs.cpp)
add_executable(addInfoToGraph mainAddInfoToGraph.cpp)
add_executable(addBothPath mainAddBothPath.cpp)
add_executable(mergeGraph mainMergeGraph.cpp)
add_executable(splitFast mainSplitFastqFiles.cpp)


target_link_libraries(build builder_lib contig_graph_lib ${SEQAN_LIBRARIES} ${ZLIB_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(filter filter_lib contig_graph_lib)
target_link_libraries(merger contig_merger_lib contig_graph_lib)
target_link_libraries(addInfoToGraph contig_graph_lib)
target_link_libraries(addBothPath contig_graph_lib)