include_directories(.)

add_subdirectory(Filter)
add_subdirectory(Builder)
add_subdirectory(Logger)
add_subdirectory(ReadsSplitter)

add_executable(buildApp mainBuildGraph.cpp)
add_executable(filterApp mainFilter.cpp)
add_executable(mergeGraph mainMergeGraph.cpp)
add_executable(readSplitter mainSplitReads.cpp)

target_link_libraries(buildApp builder_lib logger_lib ${SEQAN_LIBRARIES} ${ZLIB_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} ${Boost_LIBRARIES})
target_link_libraries(filterApp filter_lib logger_lib ${SEQAN_LIBRARIES} ${ZLIB_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} ${Boost_LIBRARIES})
target_link_libraries(mergeGraph logger_lib ${CMAKE_THREAD_LIBS_INIT} ${Boost_LIBRARIES})
target_link_libraries(readSplitter logger_lib read_splitter_lib ${SEQAN_LIBRARIES} ${ZLIB_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} ${Boost_LIBRARIES})

install(TARGETS filterApp
        DESTINATION bin
        COMPONENT runtime)
install(TARGETS buildApp
        DESTINATION bin
        COMPONENT runtime)
install(TARGETS mergeGraph
        DESTINATION bin
        COMPONENT runtime)
install(TARGETS readSplitter
        DESTINATION bin
        COMPONENT runtime)
