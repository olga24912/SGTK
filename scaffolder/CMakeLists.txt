cmake_minimum_required(VERSION 3.3)
project(scaffolder)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.cpp GraphControl.cpp GraphControl.h
        GraphBuilder/GraphBuilder.cpp GraphBuilder/GraphBuilder.h
        GraphBuilder/RNASplitReadGraphBuilder.cpp GraphBuilder/RNASplitReadGraphBuilder.h
        GraphBuilder/PairReadGraphBuilder.cpp GraphBuilder/PairReadGraphBuilder.h
        GraphBuilder/DNAPairReadGraphBuilder.cpp GraphBuilder/DNAPairReadGraphBuilder.h
        GraphBuilder/RNAPairReadGraphBuilder.cpp GraphBuilder/RNAPairReadGraphBuilder.h ConigGraph.cpp ConigGraph.h Tools/WorkWithOtherTools.cpp Tools/WorkWithOtherTools.h Tools/FastaToolsIn.cpp Tools/FastaToolsIn.h GenTests/SplitReadsTest.cpp GenTests/SplitReadsTest.h ReadsSplitter/ReadsSplitter.cpp ReadsSplitter/ReadsSplitter.h Tools/FastaToolsOut.cpp Tools/FastaToolsOut.h Tools/SeqanUtils.cpp Tools/SeqanUtils.h OutsideLib/optionparser.h ReadsSplitter/ReadsSplitter50.cpp ReadsSplitter/ReadsSplitter50.h)
add_executable(scaffolder ${SOURCE_FILES})