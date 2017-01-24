//
// Created by olga on 08.10.16.
//

#include "GraphBuilder.h"

void GraphBuilder::setMinContigLen(int minContigLen) {
    GraphBuilder::minContigLen = minContigLen;
}

void GraphBuilder::setGraph(ContigGraph *graph) {
    GraphBuilder::graph = graph;
    graph->newLib(libName, getLibColor());
}

void GraphBuilder::setLibName(string libName) {
    this->libName = libName;
}

void GraphBuilder::setSamFileWriter(SamFileWriteEdge writer) {
    this->samFileWriter = writer;
}

string GraphBuilder::colorToString(int *color) {
    string res = "#";
    for (int i = 0; i < 3; ++i) {
        if (color[i] / 16 < 10) {
            res += (color[i] / 16) + '0';
        } else {
            res += (color[i] / 16) - 10 + 'a';
        }

        if (color[i] % 16 < 10) {
            res += (color[i] % 16) + '0';
        } else {
            res += (color[i] % 16) - 10 + 'a';
        }
    }
    return res;
}