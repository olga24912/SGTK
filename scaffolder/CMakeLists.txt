cmake_minimum_required(VERSION 3.3)
project(scaffolder)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.cpp GraphControl.cpp GraphControl.h ReadsAlignment.cpp ReadsAlignment.h GraphBuilder/GraphBuilder.cpp GraphBuilder/GraphBuilder.h  GraphBuilder/RNASplitReadGraphBuilder.cpp GraphBuilder/RNASplitReadGraphBuilder.h GraphBuilder/PairReadGraphBuilder.cpp GraphBuilder/PairReadGraphBuilder.h GraphBuilder/DNAPairReadGraphBuilder.cpp GraphBuilder/DNAPairReadGraphBuilder.h GraphBuilder/RNAPairReadGraphBuilder.cpp GraphBuilder/RNAPairReadGraphBuilder.h Graph.cpp Graph.h)
add_executable(scaffolder ${SOURCE_FILES})