cmake_minimum_required(VERSION 3.3)
project(scaffolder)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.cpp GraphControl.cpp GraphControl.h
        GraphBuilder/GraphBuilder.cpp GraphBuilder/GraphBuilder.h
        GraphBuilder/RNASplitReadGraphBuilder.cpp GraphBuilder/RNASplitReadGraphBuilder.h
        GraphBuilder/PairReadGraphBuilder.cpp GraphBuilder/PairReadGraphBuilder.h
        GraphBuilder/DNAPairReadGraphBuilder.cpp GraphBuilder/DNAPairReadGraphBuilder.h
        GraphBuilder/RNAPairReadGraphBuilder.cpp GraphBuilder/RNAPairReadGraphBuilder.h
        ContigGraph/ContigGraph.cpp ContigGraph/ContigGraph.h
        Tools/WorkWithOtherTools.cpp Tools/WorkWithOtherTools.h
        GenTests/SplitReadsTest.cpp GenTests/SplitReadsTest.h
        ReadsSplitter/ReadsSplitter.cpp ReadsSplitter/ReadsSplitter.h
        Tools/SeqanUtils.cpp Tools/SeqanUtils.h OutsideLib/optionparser.h
        ReadsSplitter/ReadsSplitter50.cpp ReadsSplitter/ReadsSplitter50.h ReadsSplitter/SplitterByUnmappedEnd.cpp
        ReadsSplitter/SplitterByUnmappedEnd.h ContigGraph/ContigGraphPrinter.cpp ContigGraph/ContigGraphPrinter.h ContigGraph/Serialization.cpp ContigGraph/Serialization.h ContigGraph/InteractiveFilter.cpp ContigGraph/InteractiveFilter.h SamFileWriter/SamFileWriter.cpp SamFileWriter/SamFileWriter.h SamFileWriter/SamFileWriteEdge.cpp SamFileWriter/SamFileWriteEdge.h SamFileWriter/ContigMerger.cpp SamFileWriter/ContigMerger.h Scaffolder/ScafSimplePath.cpp Scaffolder/ScafSimplePath.h GraphBuilder/ReferenceGraphBuilder.cpp GraphBuilder/ReferenceGraphBuilder.h Tools/Utils.cpp Tools/Utils.h)
add_executable(scaffolder ${SOURCE_FILES} Tools/Utils.cpp Tools/Utils.h)