//
// Created by olga on 26.01.17.
//

#include "WriteFullGraph.h"

WriteFullGraph::WriteFullGraph(std::string fileName, Filter *filter1) : Writer(filter1), fileName(fileName) {}

void WriteFullGraph::write() {
    dotWriter.writeVertexSet(filter->getVertexList(), fileName);
}
