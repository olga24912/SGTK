include_directories(.)

set(SOURCE_FILES ContigGraph.cpp ContigGraph.h
        SeqanUtils.cpp SeqanUtils.h)

add_library(contig_graph_lib ${SOURCE_FILES})