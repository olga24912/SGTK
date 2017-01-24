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

std::string FilterAdapter::getEdgeColor(int e) {
    return graph.getLibColor(graph.getEdgeLib(e));
}

std::string FilterAdapter::getEdgeLibName(int e) {
    return graph.getLibName(graph.getEdgeLib(e));
}