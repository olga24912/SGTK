find_package (ZLIB)
find_package (BZip2)

include_directories(.)

set(SOURCE_FILES GraphControl.cpp GraphControl.h
        GraphBuilder/GraphBuilder.cpp GraphBuilder/GraphBuilder.h
        GraphBuilder/PairReadGraphBuilder.cpp GraphBuilder/PairReadGraphBuilder.h
        GraphBuilder/RNASplitGraphBuilder.cpp GraphBuilder/RNASplitGraphBuilder.h
        GraphBuilder/ConnectionGraphBuilder.cpp GraphBuilder/ConnectionGraphBuilder.h
        GraphBuilder/DNAPairReadsGraphBuilder.cpp GraphBuilder/DNAPairReadsGraphBuilder.h
        GraphBuilder/LongGraphBuilder.cpp GraphBuilder/LongGraphBuilder.h
        GraphBuilder/FASTGGraphBuilder.cpp GraphBuilder/FASTGGraphBuilder.h
        GraphBuilder/GFAGraphBuilder.cpp GraphBuilder/GFAGraphBuilder.h
        GraphBuilder/GFA2GraphBuilder.cpp GraphBuilder/GFA2GraphBuilder.h
        ContigGraph/ContigGraph.cpp ContigGraph/ContigGraph.h)

add_library(builder_lib ${SOURCE_FILES})