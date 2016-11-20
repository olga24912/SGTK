//
// Created by olga on 09.10.16.
//

#include "ContigGraph.h"

string ContigGraph::genRandomColor() {
    int color[3] = {rand() % 256, rand() % 256, rand() % 256};

    if (libColor.size() < 5) {
        string colors[5] = {"#FFFF00", "#ff0000", "#FF00FF", "#3366FF", "#33CC33"};
        return colors[libColor.size()];
    }

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

int ContigGraph::getLibNum() {
    return (int)libName.size();
}

int ContigGraph::getVertexCount() {
    return (int)graph.size();
}

int ContigGraph::getTargetLength(int id)const {
    return targetLen[id];
}

void ContigGraph::incTargetCover(int id, double x) {
    targetCoverage[id] += x;
}

int ContigGraph::addVertex(int id, string name, double cov, int len) {
    int v = (int)graph.size();
    graph.push_back(vector<int>());
    graphR.push_back(vector<int>());
    targetId[name] = id;
    vById.resize((unsigned long) max((int)vById.size(), id + 1), -1);
    vById[id] = v;

    idByV.push_back(id);

    targetName.resize(vById.size());
    targetName[id] = name;

    targetCoverage.resize(vById.size());
    targetCoverage[id] = cov;

    targetLen.resize(vById.size());
    targetLen[id] = len;

    start.push_back(0);

    edgeIdByVertexes.push_back(vector<int>((unsigned long) (v + 1), -1));

    for (int i = 0; i < v; ++i) {
        edgeIdByVertexes[i].push_back(-1);
    }
    return v;
}

int ContigGraph::incEdgeWeight(int vId, int uId) {
    int v = vById[vId], u = vById[uId];
    int e = edgeIdByVertexes[v][u];
    if (e == -1) {
        e = (int)edgeWeight.size();
        edgeWeight.push_back(0);
        to.push_back(u);
        from.push_back(v);
        edgeLib.push_back((int)libName.size() - 1);
        edgeIdByVertexes[v][u] = e;
        graph[v].push_back(e);
        graphR[u].push_back(e);
    }

    edgeWeight[e] += 1;
    return e;
}

vector<int> ContigGraph::getEdgesWeight(int v) {
    vector<int> w;
    for (int i = start[v]; i < (int)graph[v].size(); ++i) {
        w.push_back(edgeWeight[graph[v][i]]);
    }
    return w;
}

void ContigGraph::delEdges(int v, int k) {
    graph[v].resize(graph[v].size() - k);
}

void ContigGraph::sortEdgeByWeight(int v) {
    vector<pair<int, int> > edges;
    for (int i = start[v]; i < (int)graph[v].size(); ++i) {
        edges.push_back(make_pair(edgeWeight[graph[v][i]], graph[v][i]));
    }

    sort(edges.rbegin(), edges.rend());

    for (int i = start[v]; i < (int)graph[v].size(); ++i) {
        graph[v][i] = edges[i - start[v]].second;
    }
}

void ContigGraph::filterByContigLen(int minContigLen) {
    ContigGraph::minContigLen = minContigLen;
}

void ContigGraph::filterByEdgeWeight(int minEdgeWeight) {
    libMinEdgeWight[libMinEdgeWight.size() - 1] = minEdgeWeight;
}

void ContigGraph::newLib() {
    libColor.push_back(genRandomColor());
    libName.push_back("");
    libMinEdgeWight.push_back(0);
    minContigLen = 0;
    for (int i = 0; i < (int)graph.size(); ++i) {
        start[i] = (int)graph[i].size();
    }

    startEdgeNum = (int) edgeWeight.size();

    for (int i = 0; i < (int)edgeIdByVertexes.size(); ++i) {
        for (int j = 0; j < (int)edgeIdByVertexes[i].size(); ++j) {
            edgeIdByVertexes[i][j] = -1;
        }
    }

}

void ContigGraph::setLibName(string s) {
    libName[libName.size() - 1] = s;
}

void ContigGraph::setMinEdgeWeightForLib(int libNum, int minWeight) {
    if (libNum >= libMinEdgeWight.size()) return;
    libMinEdgeWight[libNum] = minWeight;
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

bool ContigGraph::isGoodEdge(int e) {
    return edgeWeight[e] >= libMinEdgeWight[edgeLib[e]];
}

bool ContigGraph::isGoodVertex(int v) {
    return targetLen[idByV[v]] >= minContigLen;
}

int ContigGraph::getFromVertex(int e) {
    return from[e];
}
