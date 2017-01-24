//
// Created by olga on 09.10.16.
//

#include "ContigGraph.h"

int ContigGraph::getLibNum() {
    return (int)libName.size();
}

int ContigGraph::getVertexCount() {
    return (int)graph.size();
}

int ContigGraph::getTargetLength(int v) const {
    return targetLen[v];
}

int ContigGraph::addVertex(int id, string name, int len) {
    assert(id == (int)graph.size());

    graph.push_back(vector<int>());
    graphR.push_back(vector<int>());
    targetId[name] = id;

    targetName.resize((size_t)id + 1);
    targetName[id] = name;

    targetLen.resize((size_t)id + 1);
    targetLen[id] = len;

    vrtsToEdge.push_back(unordered_map());
    return id;
}

int ContigGraph::incEdgeWeight(int v, int u) {
    int e = vrtsToEdge[v][u];
    if (!vrtsToEdge[v].count(u)) {
        e = (int)edgeWeight.size();
        edgeWeight.push_back(0);
        to.push_back(u);
        from.push_back(v);
        edgeLib.push_back((int)libName.size() - 1);

        vrtsToEdge[v][u] = e;

        graph[v].push_back(e);
        graphR[u].push_back(e);
    }

    edgeWeight[e] += 1;
    return e;
}

void ContigGraph::newLib(string name, string color) {
    libColor.push_back(color);
    libName.push_back(name);

    for (int i = 0; i < (int)vrtsToEdge.size(); ++i) {
        vrtsToEdge.clear();
    }
}

vector<int> ContigGraph::getEdges(int v) {
    return graph[v];
}

vector<int> ContigGraph::getEdgesR(int v) {
    return graphR[v];
}

int ContigGraph::getToVertex(int e) {
    return to[e];
}

int ContigGraph::getFromVertex(int e) {
    return from[e];
}