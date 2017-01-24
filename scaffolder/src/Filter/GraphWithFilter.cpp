//
// Created by olga on 24.01.17.
//

#include "GraphWithFilter.h"

void ContigGraph::filterByContigLen(int minContigLen) {
    ContigGraph::minContigLen = minContigLen;
}

void ContigGraph::filterByEdgeWeight(int minEdgeWeight) {
    libMinEdgeWight[libMinEdgeWight.size() - 1] = minEdgeWeight;
}

void ContigGraph::setMinEdgeWeightForLib(int libNum, int minWeight) {
    if (libNum >= libMinEdgeWight.size()) return;
    libMinEdgeWight[libNum] = minWeight;
}

bool ContigGraph::isGoodEdge(int e) {
    return edgeWeight[e] >= libMinEdgeWight[edgeLib[e]];
}

bool ContigGraph::isGoodVertex(int v) {
    return targetLen[idByV[v]] >= minContigLen;
}
