find_package (ZLIB)
find_package (BZip2)

include_directories(.)

set(SOURCE_FILES GraphControl.cpp GraphControl.h
        GraphBuilder/GraphBuilder.cpp GraphBuilder/GraphBuilder.h
        GraphBuilder/RNASplitReadGraphBuilder.cpp GraphBuilder/RNASplitReadGraphBuilder.h
        GraphBuilder/PairReadGraphBuilder.cpp GraphBuilder/PairReadGraphBuilder.h
        GraphBuilder/DNAPairReadGraphBuilder.cpp GraphBuilder/DNAPairReadGraphBuilder.h
        GraphBuilder/RNAPairReadGraphBuilder.cpp GraphBuilder/RNAPairReadGraphBuilder.h
        Tools/SystemAlignmentTools.cpp Tools/SystemAlignmentTools.h
        ReadsSplitter/ReadsSplitter.cpp ReadsSplitter/ReadsSplitter.h
        ExternalLib/optionparser.h
        ReadsSplitter/ReadsSplitter50.cpp ReadsSplitter/ReadsSplitter50.h
        ReadsSplitter/SplitterByUnmappedEnd.cpp ReadsSplitter/SplitterByUnmappedEnd.h
        SamFileWriter/SamFileWriteEdge.cpp SamFileWriter/SamFileWriteEdge.h
        GraphBuilder/ReferenceGraphBuilder.cpp GraphBuilder/ReferenceGraphBuilder.h
        SamFileWriter/SamFileWriterEdgeFictive.cpp SamFileWriter/SamFileWriterEdgeFictive.h)

add_library(builder_lib ${SOURCE_FILES})