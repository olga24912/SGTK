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

int FilterAdapter::getEdgeWeight(int e) {
    return graph.getEdgeWeight(e);
}

std::vector<int> FilterAdapter::getLibList() {
    std::vector<int> libList(graph.getLibNum(), 0);
    std::cerr << graph.getLibNum() << std::endl;
    for (int i = 0; i < (int)libList.size(); ++i) {
        libList[i] = i;
    }
    return libList;
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

void FilterAdapter::processQuery(Query query) {
    if (query.type == query.UPLOAD_GRAPH) {
        std::string fileName = query.argv;
        graph = ContigGraph::read(fileName);
    }
}

std::string FilterAdapter::getInfo(int e) {
    return graph.getEdgeInfo(e);
}
