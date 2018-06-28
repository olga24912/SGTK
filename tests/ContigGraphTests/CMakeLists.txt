include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

add_executable(runContigGraphTests
        contigGraphTest.cpp)

target_link_libraries(runContigGraphTests gtest gtest_main)
target_link_libraries(runContigGraphTests contig_graph_lib)
target_link_libraries(runContigGraphTests logger_lib)