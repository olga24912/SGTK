include_directories(.)

set(SOURCE_FILES Builder/GraphControl.cpp Builder/GraphControl.h
        Builder/GraphBuilder/GraphBuilder.cpp Builder/GraphBuilder/GraphBuilder.h
        Builder/GraphBuilder/RNASplitReadGraphBuilder.cpp Builder/GraphBuilder/RNASplitReadGraphBuilder.h
        Builder/GraphBuilder/PairReadGraphBuilder.cpp Builder/GraphBuilder/PairReadGraphBuilder.h
        Builder/GraphBuilder/DNAPairReadGraphBuilder.cpp Builder/GraphBuilder/DNAPairReadGraphBuilder.h
        Builder/GraphBuilder/RNAPairReadGraphBuilder.cpp Builder/GraphBuilder/RNAPairReadGraphBuilder.h
        ContigGraph/ContigGraph.cpp ContigGraph/ContigGraph.h
        Builder/Tools/SystemAlignmentTools.cpp Builder/Tools/SystemAlignmentTools.h
        ../tests/GenTests/SplitReadsTest.cpp ../tests/GenTests/SplitReadsTest.h
        Builder/ReadsSplitter/ReadsSplitter.cpp Builder/ReadsSplitter/ReadsSplitter.h
        Builder/Tools/SeqanUtils.cpp Builder/Tools/SeqanUtils.h Builder/ExternalLib/optionparser.h
        Builder/ReadsSplitter/ReadsSplitter50.cpp Builder/ReadsSplitter/ReadsSplitter50.h
        Builder/ReadsSplitter/SplitterByUnmappedEnd.cpp Builder/ReadsSplitter/SplitterByUnmappedEnd.h
        Filter/Writers/ContigGraphPrinter.cpp Filter/Writers/ContigGraphPrinter.h
        Filter/CommanParsers/InteractiveFilter.cpp Filter/CommanParsers/InteractiveFilter.h
        Builder/SamFileWriter/SamFileWriter.cpp Builder/SamFileWriter/SamFileWriter.h
        Builder/SamFileWriter/SamFileWriteEdge.cpp Builder/SamFileWriter/SamFileWriteEdge.h
        Builder/SamFileWriter/ContigMerger.cpp Builder/SamFileWriter/ContigMerger.h
        Filter/Scaffolder/ScafSimplePath.cpp Filter/Scaffolder/ScafSimplePath.h
        Builder/GraphBuilder/ReferenceGraphBuilder.cpp Builder/GraphBuilder/ReferenceGraphBuilder.h
        Filter/Tools/Utils.cpp Filter/Tools/Utils.h
        Filter/Tools/SystemTools.cpp Filter/Tools/SystemTools.h Filter/Filters/Filter.h Filter/Filters/FilterAdapter.cpp Filter/Filters/FilterAdapter.h Filter/Filters/FilterMinWeight.cpp Filter/Filters/FilterMinWeight.h Filter/Filters/Filter.cpp Filter/Filters/FilterIgnore.cpp Filter/Filters/FilterIgnore.h Filter/Writers/DotWriter.cpp Filter/Writers/DotWriter.h Filter/Writers/Searcher.cpp Filter/Writers/Searcher.h Filter/Writers/Writer.h Filter/Writers/WriteFullGraph.cpp Filter/Writers/WriteFullGraph.h Filter/Writers/writeLocal.cpp Filter/Writers/writeLocal.h)
add_library(scaffolder_lib ${SOURCE_FILES})