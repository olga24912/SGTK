#include "WriteFullGraph.h"

namespace filter {
    namespace writers {
        WriteFullGraph::WriteFullGraph(std::string fileName, ContigGraph *graph1, FileValidator *validator,
                                       DotWriterBuilder *builder) : Writer(graph1, validator, builder),
                                                                    fileName(fileName) {}

        void WriteFullGraph::write() {
            INFO("Write full graph");
            dotWriter->writeVertexSet(graph->getVertexList(), fileName);
        }
    }
}