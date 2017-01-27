include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

set(SOURCE_FILES GenTests/SplitReadsTest.h GenTests/SplitReadsTest.cpp)

add_executable(runBuilderTests builderTest.cpp ${SOURCE_FILES})

target_link_libraries(runBuilderTests gtest gtest_main)
target_link_libraries(runBuilderTests scaffolder_lib)