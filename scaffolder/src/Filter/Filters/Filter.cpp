#include <fstream>
#include <iostream>
#include <algorithm>
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

std::vector<int> Filter::getLibList() {
    assert(subfilter != nullptr);
    return subfilter->getLibList();
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

    std::vector<int> libList = getLibList();
    out << libList.size() << "\n";
    std::vector<int> newColr(libList[libList.size() - 1] + 1, 0);
    for (int i = 0; i < (int)libList.size(); ++i) {
        out << "l " << i << " " << getLibColor(libList[i]) << " " << getLibName(libList[i]) << "\n";
        newColr[libList[i]] = i;
    }

    std::vector<int> edges;
    out << getVertexCount() << "\n";
    for (int i = 0; i < getVertexCount(); ++i) {
        out << "v " << i << " " << getTargetName(i) << " " << getTargetLen(i) << "\n";
        std::vector<int> curEd = getEdges(i);
        for (int j = 0; j < (int)curEd.size(); ++j) {
            edges.push_back(curEd[j]);
        }
    }

    std::sort(edges.begin(), edges.end());
    edges.resize(unique(edges.begin(), edges.end()) - edges.begin());
    out << edges.size() << "\n";
    for (int i = 0; i < edges.size(); ++i) {
        out << "e " << i << " " << getEdgeFrom(edges[i]) << " " << getEdgeTo(edges[i]) << " "
            << newColr[getEdgeLib(edges[i])] << " " << getEdgeWeight(edges[i]) << "\n";
    }

    out.close();
}

std::string Filter::getInfo(int e) {
    return subfilter->getInfo(e);
}

std::pair<int, int> Filter::getFirstCoord(int e) {
    return subfilter->getFirstCoord(e);
}

std::pair<int, int> Filter::getSecondCoord(int e) {
    return subfilter->getFirstCoord(e);
}

ContigGraph::Lib::Type Filter::getLibType(int l) {
    return subfilter->getLibType(l);
}

int Filter::getEdgeCoordB1(int e) {
    return subfilter->getEdgeCoordB1(e);
}

int Filter::getEdgeCoordE1(int e) {
    return subfilter->getEdgeCoordE1(e);
}

int Filter::getEdgeCoordB2(int e) {
    return subfilter->getEdgeCoordB2(e);
}

int Filter::getEdgeCoordE2(int e) {
    return subfilter->getEdgeCoordE2(e);
}