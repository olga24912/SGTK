#include "FilterAdapter.h"

int FilterAdapter::getVertexCount() {
    TRACE("getVertexCount");
    return graph.getVertexCount();
}

std::vector<int> FilterAdapter::getEdges(int v) {
    TRACE("getEdges v=" << v);
    return graph.getEdges(v);
}

std::vector<int> FilterAdapter::getEdgesR(int v) {
    TRACE("getEdgesR v=" << v);
    return graph.getEdgesR(v);
}

std::string FilterAdapter::getTargetName(int v) {
    TRACE("get target name v=" << v);
    return graph.getTargetName(v);
}

int FilterAdapter::getTargetLen(int v) {
    TRACE("get target len v=" << v);
    return graph.getTargetLength(v);
}

int FilterAdapter::getEdgeTo(int e) {
    TRACE("get vertex to for e=" << e);
    return graph.getToVertex(e);
}

int FilterAdapter::getEdgeFrom(int e) {
    TRACE("get edge from for e=" << e);
    return graph.getFromVertex(e);
}

int FilterAdapter::getEdgeWeight(int e) {
    TRACE("get edge weight e=" << e);
    return graph.getEdgeWeight(e);
}

std::vector<int> FilterAdapter::getLibList() {
    TRACE("get Lib list, lib num in graph=" << graph.getLibNum());
    std::vector<int> libList(graph.getLibNum(), 0);
    for (int i = 0; i < (int)libList.size(); ++i) {
        libList[i] = i;
    }
    return libList;
}

int FilterAdapter::getEdgeLib(int e) {
    TRACE("getEdgeLib e=" << e);
    return graph.getEdgeLib(e);
}

std::string FilterAdapter::getLibName(int l) {
    TRACE("getLibName l=" << l);
    return graph.getLibName(l);
}

std::string FilterAdapter::getLibColor(int l) {
    TRACE("getLibColor l=" << l);
    return graph.getLibColor(l);
}

std::vector<int> FilterAdapter::getVertexList() {
    TRACE("getVertexList");

    int vertNum = getVertexCount();
    std::vector<int> res((unsigned)vertNum);
    for (int i = 0; i < vertNum; ++i) {
        res[i] = i;
    }
    return res;
}

void FilterAdapter::processQuery(Query query) {
    TRACE("processQuery")
    if (query.type == query.UPLOAD_GRAPH) {
        std::string fileName = query.argv;
        TRACE("query upload graph fileName=" << fileName);
        graph = contig_graph::ContigGraph::read(fileName);
    }
}

std::string FilterAdapter::getInfo(int e) {
    TRACE("getInfo e=" << e);
    return graph.getEdgeInfo(e);
}

std::pair<int, int> FilterAdapter::getFirstCoord(int e) {
    TRACE("getFirstCoord e=" << e);
    return graph.getFirstCoord(e);
}

std::pair<int, int> FilterAdapter::getSecondCoord(int e) {
    TRACE("getSecondCoord e=" << e);
    return graph.getSecondCoord(e);
}

contig_graph::ContigGraph::Lib::Type FilterAdapter::getLibType(int l) {
    TRACE("getLibType l=" << l);
    return graph.getLibType(l);
}

int FilterAdapter::getEdgeCoordB1(int e) {
    TRACE("getEdgeCoordB1 e=" << e);
    return graph.getEdgeCoordB1(e);
}

int FilterAdapter::getEdgeCoordE1(int e) {
    TRACE("getEdgeCoordE1 e=" << e);
    return graph.getEdgeCoordE1(e);
}

int FilterAdapter::getEdgeCoordB2(int e) {
    TRACE("getEdgeCoordB2 e=" << e);
    return graph.getEdgeCoordB2(e);
}

int FilterAdapter::getEdgeCoordE2(int e) {
    TRACE("getEdgeCoordE2 e=" << e);
    return graph.getEdgeCoordE2(e);
}
