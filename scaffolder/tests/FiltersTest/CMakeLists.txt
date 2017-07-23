include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

add_executable(runFiltersTests contigGraphTest.cpp)

target_link_libraries(runFiltersTests gtest gtest_main)
target_link_libraries(runFiltersTests filter_lib)
target_link_libraries(runFiltersTests logger_lib)
