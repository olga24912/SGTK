//
// Created by olga on 08.10.16.
//

#include "GraphBuilder.h"

void GraphBuilder::setMinContigLen(int minContigLen) {
    GraphBuilder::minContigLen = minContigLen;
}

void GraphBuilder::setMinEdgeWight(int minEdge) {
    GraphBuilder::minEdgeWight = minEdge;
}

void GraphBuilder::setGraph(ContigGraph *graph) {
    GraphBuilder::graph = graph;
    graph->setColor(getLibColor());
}

void GraphBuilder::setLibName(string libName) {
    this->libName = libName;
}

void GraphBuilder::setSamFileWriter(SamFileWriteEdge writer) {
    this->samFileWriter = writer;
}
