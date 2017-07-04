#include <queue>
#include <iostream>
#include "GraphSplitter.h"

std::vector<std::vector<int> > GraphSplitter::split(Filter *filter, std::vector<int> vert) {
    INFO("start split file on parts");
    this->filter = filter;
    this->vert = vert;

    clear();

    for (int v : vert) {
        while (used[v] == 0) {
            findNewComp(v);
        }
    }

    return res;
}

void GraphSplitter::clear() {
    DEBUG("clear");
    used.resize(0);
    res.resize(0);
    edgeCol.resize(0);

    used.resize((unsigned int)filter->getVertexCount(), 2);
    for (int v : vert) {
        used[v] = 0;
    }


    int edgeCnt = 0;
    for (int i = 0; i < (unsigned)used.size(); ++i) {
        for (int e : filter->getEdges(i)) {
            edgeCnt = std::max(edgeCnt, e + 1);
        }
    }

    edgeCol.resize(edgeCnt);
}

void GraphSplitter::findNewComp(int v) {
    TRACE("find new components for v=" << v);
    std::queue<int> que;
    int cntV = 0;
    int cntE = 0;

    int colNum = (int)res.size();
    res.push_back(std::vector<int>());

    que.push(v);
    while (que.size() > 0) {
        if (cntV == maxVert) break;
        int u = que.front();
        que.pop();

        int extraEdge = 0;
        for (int e : filter->getEdges(u)) {
            int w = filter->getEdgeTo(e);
            if (used[w] == 1) {
                ++extraEdge;
            }
        }

        for (int e : filter->getEdgesR(u)) {
            int w = filter->getEdgeFrom(e);
            if (used[w] == 1) {
                ++extraEdge;
            }
        }

        if (cntE + extraEdge <= maxEdge && cntV + 1 <= maxVert) {
            res[colNum].push_back(u);
            cntE += extraEdge;
            cntV += 1;
            used[u] = 1;
            for (int e : filter->getEdges(u)) {
                int w = filter->getEdgeTo(e);
                if (used[w] == 1) {
                    edgeCol[e] = colNum + 1;
                } else if (edgeCol[e] == 0 && used[w] == 0) {
                    used[w] = 3;
                    que.push(w);
                }
            }

            for (int e : filter->getEdgesR(u)) {
                int w = filter->getEdgeFrom(e);
                if (used[w] == 1) {
                    edgeCol[e] = colNum + 1;
                }  else if (edgeCol[e] == 0 && used[w] == 0) {
                    used[w] = 3;
                    que.push(w);
                }
            }

        } else {
            used[u] = 3;
        }
    }

    for (int i = 0; i < (int)used.size(); ++i) {
        if (used[i] == 3) used[i] = 0;
        if (used[i] == 1) {
            used[i] = 2;
            for (int e : filter->getEdges(i)) {
                if (edgeCol[e] == 0 && used[filter->getEdgeTo(e)] != 2) {
                    used[i] = 0;
                }
            }


            for (int e : filter->getEdgesR(i)) {
                if (edgeCol[e] == 0 && used[filter->getEdgeFrom(e)] != 2) {
                    used[i] = 0;
                }
            }
        }
    }
}
