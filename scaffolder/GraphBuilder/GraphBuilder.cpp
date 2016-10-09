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

void GraphBuilder::setGraph(ConigGraph *graph) {
    GraphBuilder::graph = graph;
}
