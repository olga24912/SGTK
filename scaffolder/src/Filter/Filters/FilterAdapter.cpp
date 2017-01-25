#include "FilterAdapter.h"

int FilterAdapter::getVertexCount() {
    return graph.getVertexCount();
}

std::vector<int> FilterAdapter::getEdges(int v) {
    return graph.getEdges(v);
}

std::vector<int> FilterAdapter::getEdgesR(int v) {
    return graph.getEdgesR(v);
}

std::string FilterAdapter::getTargetName(int v) {
    return graph.getTargetName(v);
}

int FilterAdapter::getTargetLen(int v) {
    return graph.getTargetLength(v);
}

int FilterAdapter::getEdgeTo(int e) {
    return graph.getToVertex(e);
}

int FilterAdapter::getEdgeFrom(int e) {
    return graph.getFromVertex(e);
}

int FilterAdapter::getEdgeWieght(int e) {
    return graph.getEdgeWeight(e);
}

int FilterAdapter::getLibCount() {
    return graph.getLibNum();
}

int FilterAdapter::getEdgeLib(int e) {
    return graph.getEdgeLib(e);
}

std::string FilterAdapter::getLibName(int l) {
    return graph.getLibName(l);
}

std::string FilterAdapter::getLibColor(int l) {
    return graph.getLibColor(l);
}

std::vector<int> FilterAdapter::getVertexList() {
    int vertNum = getVertexCount();
    std::vector<int> res((unsigned)vertNum);
    for (int i = 0; i < vertNum; ++i) {
        res[i] = i;
    }
    return res;
}
