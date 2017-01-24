#include "Filter.h"

int Filter::getVertexCount() {
    assert(subfilter != nullptr);
    return subfilter->getVertexCount();
}

std::vector<int> Filter::getEdges(int v) {
    assert(subfilter != nullptr);
    return subfilter->getEdges();
}

std::vector<int> Filter::getEdgesR(int v) {
    assert(subfilter != nullptr);
    return subfilter->getEdgesR(v);
}

std::string Filter::getTargetName(int v);
int Filter::getTargetLen(int v);
int Filter::getEdgeTo(int e);
int Filter::getEdgeFrom(int e);
int Filter::getEdgeWieght(int e);
std::string Filter::getEdgeColor(int e);
std::string Filter::getEdgeLibName(int e);