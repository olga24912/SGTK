include_directories(.)

add_subdirectory(Filter)
add_subdirectory(Builder)
#add_subdirectory(ContigMerger)
add_subdirectory(Logger)
add_subdirectory(ReadsSplitter)

add_executable(build mainBuildGraph.cpp)
add_executable(filter mainFilter.cpp)
#add_executable(merger mainMergeContigs.cpp)
add_executable(addInfoToGraph mainAddInfoToGraph.cpp)
add_executable(addBothPath mainAddBothPath.cpp)
add_executable(mergeGraph mainMergeGraph.cpp)
#add_executable(splitFast mainSplitFastqFiles.cpp)
add_executable(readSplitter mainSplitReads.cpp)

target_link_libraries(build builder_lib logger_lib ${SEQAN_LIBRARIES} ${ZLIB_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} ${Boost_LIBRARIES})
target_link_libraries(filter filter_lib logger_lib ${SEQAN_LIBRARIES} ${ZLIB_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} ${Boost_LIBRARIES})
#target_link_libraries(merger contig_merger_lib contig_graph_lib logger_lib ${SEQAN_LIBRARIES} ${ZLIB_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} ${Boost_LIBRARIES})
target_link_libraries(addInfoToGraph logger_lib builder_lib ${CMAKE_THREAD_LIBS_INIT} ${Boost_LIBRARIES})
target_link_libraries(addBothPath logger_lib builder_lib)
target_link_libraries(mergeGraph logger_lib ${CMAKE_THREAD_LIBS_INIT} ${Boost_LIBRARIES})
target_link_libraries(readSplitter logger_lib read_splitter_lib)


install(TARGETS addBothPath
        DESTINATION bin
        COMPONENT runtime)
install(TARGETS addInfoToGraph
        DESTINATION bin
        COMPONENT runtime)
install(TARGETS filter
        DESTINATION bin
        COMPONENT runtime)
install(TARGETS build
        DESTINATION bin
        COMPONENT runtime)
install(TARGETS mergeGraph
        DESTINATION bin
        COMPONENT runtime)