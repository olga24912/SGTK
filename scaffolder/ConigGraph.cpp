//
// Created by olga on 09.10.16.
//

#include "ConigGraph.h"

string ConigGraph::genRandomColor() {
    int color[3] = {rand() % 256, rand() % 256, rand() % 256};
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

int ConigGraph::getLibNum() {
    return libNum;
}

int ConigGraph::getVertexCount() {
    return (int)graph.size();
}

int ConigGraph::getTargetLength(int id)const {
    return targetLen[id];
}

void ConigGraph::incTargetCover(int id, double x) {
    targetCoverage[id] += x;
}

int ConigGraph::addVertex(int id, string name, double cov, int len) {
    int v = (int)graph.size();
    graph.push_back(vector<int>());
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

void ConigGraph::incEdgeWeight(int vId, int uId) {
    int v = vById[vId], u = vById[uId];
    int e = edgeIdByVertexes[v][u];
    if (e == -1) {
        e = edgeWeight.size();
        edgeWeight.push_back(0);
        to.push_back(u);
        edgeLib.push_back(libNum);
        edgeIdByVertexes[v][u] = e;
        graph[v].push_back(e);
    }

    edgeWeight[e] += 1;
}

vector<int> ConigGraph::getEdgesWeight(int v) {
    vector<int> w;
    for (int i = start[v]; i < (int)graph[v].size(); ++i) {
        w.push_back(edgeWeight[graph[v][i]]);
    }
    return w;
}

void ConigGraph::delEdges(int v, int k) {
    graph[v].resize(graph[v].size() - k);
}

void ConigGraph::sortEdgeByWeight(int v) {
    vector<pair<int, int> > edges;
    for (int i = start[v]; i < (int)graph[v].size(); ++i) {
        edges.push_back(make_pair(edgeWeight[graph[v][i]], graph[v][i]));
    }

    sort(edges.rbegin(), edges.rend());

    for (int i = start[v]; i < (int)graph[v].size(); ++i) {
        graph[v][i] = edges[i - start[v]].second;
    }
}

void ConigGraph::writeGraphDotFormat(string fileName) {
    ofstream out(fileName);

    out << "digraph {\n";

    for (int i = 0; i < targetName.size(); ++i) {
        out << "    \"" << targetName[i] << "\"[label=\" " << targetName[i] << "\nlen = " << targetLen[i]
        << ", cover = "<< targetCoverage[i] <<"\"];\n";
    }

    for (int v = 0; v < (int)graph.size(); ++v) {
        int vId = idByV[v];
        if (targetLen[vId] < minContigLen) continue;
        for (int j = 0; j < (int)graph[v].size(); ++j) {
            int e = graph[v][j];
            int u = to[e];
            int uId = idByV[u];
            out << "    \"" << targetName[vId] << "\" -> \"" << targetName[uId] << "\" [ ";
            out << "color = \"" << libColor[edgeLib[e]] << "\", ";
            out << "penwidth = "<< 1 + (int)log10(edgeWeight[e]) << ", ";
            out << "label = " << "\" weight = " << edgeWeight[e] << "\" ]\n";
        }
    }

    out << "}\n";

    out.close();
}

void ConigGraph::filterByContigLen(int minContigLen) {
    ConigGraph::minContigLen = minContigLen;
}

void ConigGraph::filterByEdgeWeight(int minEdgeWeight) {
    for (int v = 0; v < (int)graph.size(); ++v) {
        for (int i = start[v]; i < (int)graph[v].size(); ++i) {
            if (edgeWeight[graph[v][i]] < minEdgeWeight) {
                swap(graph[v][i], graph[v][graph[v].size() - 1]);
                graph[v].pop_back();
                --i;
            }
        }
    }
}

void ConigGraph::newLib() {
    ++libNum;
    libColor.push_back(genRandomColor());
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
