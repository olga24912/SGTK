#include "WriteFullGraph.h"

WriteFullGraph::WriteFullGraph(std::string fileName, Filter *filter1, FileValidator *validator,
                               DotWriterBuilder* builder) : Writer(filter1, validator, builder), fileName(fileName) {}

void WriteFullGraph::write() {
    INFO("Write full graph");
    dotWriter->writeVertexSet(filter->getVertexList(), fileName);
}
