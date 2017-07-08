include_directories(.)

add_subdirectory(Filter)
add_subdirectory(ContigGraph)
add_subdirectory(Builder)
add_subdirectory(ContigMerger)
add_subdirectory(Logger)

add_executable(build mainBuildGraph.cpp)
add_executable(filter mainFilter.cpp)
#add_executable(merger mainMergeContigs.cpp)
add_executable(addInfoToGraph mainAddInfoToGraph.cpp)
add_executable(addBothPath mainAddBothPath.cpp)
add_executable(mergeGraph mainMergeGraph.cpp)
add_executable(splitFast mainSplitFastqFiles.cpp)

target_link_libraries(build builder_lib contig_graph_lib logger_lib ${SEQAN_LIBRARIES} ${ZLIB_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} ${Boost_LIBRARIES})
target_link_libraries(filter filter_lib contig_graph_lib  logger_lib ${SEQAN_LIBRARIES} ${ZLIB_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} ${Boost_LIBRARIES})
#target_link_libraries(merger contig_merger_lib contig_graph_lib logger_lib ${SEQAN_LIBRARIES} ${ZLIB_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} ${Boost_LIBRARIES})
target_link_libraries(addInfoToGraph contig_graph_lib logger_lib ${CMAKE_THREAD_LIBS_INIT} ${Boost_LIBRARIES})
target_link_libraries(addBothPath contig_graph_lib logger_lib)
target_link_libraries(mergeGraph logger_lib ${CMAKE_THREAD_LIBS_INIT} ${Boost_LIBRARIES})

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