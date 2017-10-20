include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

add_executable(runFiltersTests contigGraphTest.cpp)

target_link_libraries(runFiltersTests gtest gtest_main)
target_link_libraries(runFiltersTests filter_lib logger_lib ${SEQAN_LIBRARIES} ${ZLIB_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} ${Boost_LIBRARIES})
