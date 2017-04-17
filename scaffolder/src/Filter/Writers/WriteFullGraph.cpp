//
// Created by olga on 26.01.17.
//

#include "WriteFullGraph.h"

WriteFullGraph::WriteFullGraph(std::string fileName, Filter *filter1, FileValidator *validator,
                               int maxVert, int maxEdge) : Writer(filter1, validator,  maxVert, maxEdge), fileName(fileName) {}

void WriteFullGraph::write() {
    dotWriter.writeVertexSet(filter->getVertexList(), fileName);
}
