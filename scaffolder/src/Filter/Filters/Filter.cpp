#include "Filter.h"

int Filter::getVertexCount() {
    assert(subfilter != nullptr);
    return subfilter->getVertexCount();
}

std::vector<int> Filter::getEdges(int v) {
    assert(subfilter != nullptr);
    return subfilter->getEdges(v);
}

std::vector<int> Filter::getEdgesR(int v) {
    assert(subfilter != nullptr);
    return subfilter->getEdgesR(v);
}

std::string Filter::getTargetName(int v) {
    assert(subfilter != nullptr);
    return subfilter->getTargetName(v);
}

int Filter::getTargetLen(int v) {
    assert(subfilter != nullptr);
    return subfilter->getTargetLen(v);
}

int Filter::getEdgeTo(int e) {
    assert(subfilter != nullptr);
    return subfilter->getEdgeTo(e);
}

int Filter::getEdgeFrom(int e) {
    assert(subfilter != nullptr);
    return subfilter->getEdgeFrom(e);
}

int Filter::getEdgeWieght(int e) {
    assert(subfilter != nullptr);
    return subfilter->getEdgeWieght(e);
}

int Filter::getLibCount() {
    assert(subfilter != nullptr);
    return subfilter->getLibCount();
}

int Filter::getEdgeLib(int e) {
    assert(subfilter != nullptr);
    return subfilter->getEdgeLib(e);
}

std::string Filter::getLibName(int l) {
    assert(subfilter != nullptr);
    return subfilter->getLibName(l);
}

std::string Filter::getLibColor(int l) {
    assert(subfilter != nullptr);
    return subfilter->getLibColor(l);
}

std::vector<int> Filter::getVertexList() {
    assert(subfilter != nullptr);
    return subfilter->getVertexList();
}
