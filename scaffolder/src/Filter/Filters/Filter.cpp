#include <fstream>
#include <iostream>
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

int Filter::getEdgeWeight(int e) {
    assert(subfilter != nullptr);
    return subfilter->getEdgeWeight(e);
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

void Filter::processQuery(Query query) {
    assert(subfilter != nullptr);
    return subfilter->processQuery(query);
}

void Filter::write(std::string fileName) {
    std::ofstream out(fileName);

    std::cerr << "Write -> gr" << " " << fileName << std::endl;

    out << getLibCount() << "\n";
    for (int i = 0; i < (int)getLibCount(); ++i) {
        out << "l " << i << " " << getLibColor(i) << " " << getLibName(i) << "\n";
    }

    int edgeCnt = 0;

    out << getVertexCount() << "\n";
    for (int i = 0; i < getVertexCount(); ++i) {
        out << "v " << i << " " << getTargetName(i) << " " << getTargetLen(i) << "\n";
        std::vector<int> edges = getEdges(i);
        for (int j = 0; j < (int)edges.size(); ++j) {
            edgeCnt = std::max(edgeCnt, edges[j]);
        }
    }

    out << edgeCnt << "\n";
    for (int i = 0; i < edgeCnt; ++i) {
        out << "e " << i << " " << getEdgeFrom(i) << " " << getEdgeTo(i) << " "
            << getEdgeLib(i) << " " << getEdgeWeight(i) << "\n";
    }

    out.close();
}
