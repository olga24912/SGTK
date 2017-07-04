include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

set(SOURCE_FILES GenTests/SplitReadsTest.h GenTests/SplitReadsTest.cpp
        GenTests/RnaPairReadsTest.cpp GenTests/RnaPairReadsTest.h
        GenTests/DnaPairReadsTest.cpp GenTests/DnaPairReadsTest.h
        GenTests/PairReadsTest.cpp GenTests/PairReadsTest.h
        GenTests/GenTest.cpp GenTests/GenTest.h)

add_executable(runBuilderTests builderTest.cpp ${SOURCE_FILES})

target_link_libraries(runBuilderTests gtest gtest_main)
target_link_libraries(runBuilderTests builder_lib)
target_link_libraries(runBuilderTests contig_graph_lib)
target_link_libraries(runBuilderTests ${SEQAN_LIBRARIES} ${ZLIB_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})