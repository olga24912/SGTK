#include <queue>
#include <set>
#include <Logger/logger.hpp>
#include "Searcher.h"

namespace filter {
    namespace writers {
        std::vector<int> Searcher::findVertInLocalArea(int v, int dist) {
            DEBUG("find vertex in local area v=" << v << "dist=" << dist);
            std::vector<int> res;
            std::queue<std::pair<int, int> > que;
            std::set<int> used;

            que.push(std::make_pair(v, 0));
            while (que.size() > 0) {
                int y = que.front().first, d = que.front().second;
                que.pop();

                res.push_back(y);
                if (d == dist) continue;
                for (int e : graph->getEdges(y)) {
                    int u = graph->getEdgeTo(e);
                    if (!used.count(u)) {
                        used.insert(u);
                        que.push(std::make_pair(u, d + 1));
                    }
                }

                for (int e : graph->getEdgesR(y)) {
                    int u = graph->getEdgeFrom(e);
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
            DEBUG("find components");
            int n = graph->getVertexCount();
            std::vector<int> vert = graph->getVertexList();

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
            TRACE("dfs for find components v=" << v << " currentCol=" << currentCol);
            color[v] = currentCol;

            for (int e : graph->getEdges(v)) {
                int u = graph->getEdgeTo(e);
                if (color[u] == 0) {
                    dfsFindComponent(u, currentCol, color);
                }
            }

            for (int e : graph->getEdgesR(v)) {
                int u = graph->getEdgeFrom(e);
                if (color[u] == 0) {
                    dfsFindComponent(u, currentCol, color);
                }
            }
        }
    }
}