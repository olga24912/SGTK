#include <queue>
#include <set>
#include "Searcher.h"

std::vector<int> Searcher::findVertInLocalArea(int v, int dist) {
    std::vector<int> res;
    std::queue<std::pair<int, int> > que;
    std::set<int> used;

    que.push(std::make_pair(v, 0));
    while (que.size() > 0) {
        int y = que.front().first, d = que.front().second;
        que.pop();

        res.push_back(y);
        if (d == dist) continue;
        for (int e : filter->getEdges(y)) {
            int u = filter->getEdgeTo(e);
            if (!used.count(u)) {
                used.insert(u);
                que.push(std::make_pair(u, d + 1));
            }
        }

        for (int e : filter->getEdgesR(y)) {
            int u = filter->getEdgeFrom(e);
            if (!used.count(u)) {
                used.insert(u);
                que.push(std::make_pair(u, d + 1));
            }
        }
    }

    std::sort(res.begin(), res.end());
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
