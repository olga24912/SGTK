#include "Searcher.h"

std::vector<int> Searcher::findVertInLocalArea(int v, int dist) {
    std::vector<int> res;
    if (dist == 0) return res;
    res.push_back(v);
    for (int e : filter->getEdges(v)) {
        int u = filter->getEdgeTo(e);
        std::vector<int> add = findVertInLocalArea(u, dist - 1);
        for (int j = 0; j < (int)add.size(); ++j) {
            res.push_back(add[j]);
        }
    }

    for (int e : filter->getEdges(v)) {
        int u = filter->getEdgeFrom(e);
        std::vector<int> add = findVertInLocalArea(u, dist - 1);
        for (int j = 0; j < (int)add.size(); ++j) {
            res.push_back(add[j]);
        }
    }

    std::sort(res.begin(), res.end());
    res.resize((unsigned long) (std::unique(res.begin(), res.end()) - res.begin()));

    return res;
}

int Searcher::findComponent(int *col) {
    int n = filter->getVertexCount();
    std::vector<int> vert = filter->getVertexList();

    int cur = 1;
    for (int i = 0; i < n; ++i) {
        col[i] = 0;
    }

    for (int v : vert) {
        if (col[v] == 0) {
            dfsFindComponent(v, cur, col);
            ++cur;
        }
    }
    return cur;
}

void Searcher::dfsFindComponent(int v, int currentCol, int *color) {
    color[v] = currentCol;

    for (int e : filter->getEdges(v)) {
        int u = filter->getEdgeTo(e);
        if (color[u] == 0) {
            dfsFindComponent(u, currentCol, color);
        }
    }

    for (int e : filter->getEdgesR(v)) {
        int u = filter->getEdgeFrom(e);
        if (color[u] == 0) {
            dfsFindComponent(u, currentCol, color);
        }
    }
}
