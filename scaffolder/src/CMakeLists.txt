include_directories(.)

set(SOURCE_FILES Builder/GraphControl.cpp Builder/GraphControl.h
        Builder/GraphBuilder/GraphBuilder.cpp Builder/GraphBuilder/GraphBuilder.h
        Builder/GraphBuilder/RNASplitReadGraphBuilder.cpp Builder/GraphBuilder/RNASplitReadGraphBuilder.h
        Builder/GraphBuilder/PairReadGraphBuilder.cpp Builder/GraphBuilder/PairReadGraphBuilder.h
        Builder/GraphBuilder/DNAPairReadGraphBuilder.cpp Builder/GraphBuilder/DNAPairReadGraphBuilder.h
        Builder/GraphBuilder/RNAPairReadGraphBuilder.cpp Builder/GraphBuilder/RNAPairReadGraphBuilder.h
        ContigGraph/ContigGraph.cpp ContigGraph/ContigGraph.h
        Tools/SystemTools.cpp Tools/SystemTools.h
        ../tests/GenTests/SplitReadsTest.cpp ../tests/GenTests/SplitReadsTest.h
        Builder/ReadsSplitter/ReadsSplitter.cpp Builder/ReadsSplitter/ReadsSplitter.h
        Tools/SeqanUtils.cpp Tools/SeqanUtils.h Builder/ExternalLib/optionparser.h
        Builder/ReadsSplitter/ReadsSplitter50.cpp Builder/ReadsSplitter/ReadsSplitter50.h Builder/ReadsSplitter/SplitterByUnmappedEnd.cpp
        Builder/ReadsSplitter/SplitterByUnmappedEnd.h Filter/ContigGraphPrinter.cpp Filter/ContigGraphPrinter.h
        ContigGraph/Serialization.cpp ContigGraph/Serialization.h Filter/InteractiveFilter.cpp
        Filter/InteractiveFilter.h Builder/SamFileWriter/SamFileWriter.cpp Builder/SamFileWriter/SamFileWriter.h
        Builder/SamFileWriter/SamFileWriteEdge.cpp Builder/SamFileWriter/SamFileWriteEdge.h Builder/SamFileWriter/ContigMerger.cpp
        Builder/SamFileWriter/ContigMerger.h Filter/Scaffolder/ScafSimplePath.cpp Filter/Scaffolder/ScafSimplePath.h
        Builder/GraphBuilder/ReferenceGraphBuilder.cpp Builder/GraphBuilder/ReferenceGraphBuilder.h Tools/Utils.cpp Tools/Utils.h)
add_library(scaffolder_lib ${SOURCE_FILES})