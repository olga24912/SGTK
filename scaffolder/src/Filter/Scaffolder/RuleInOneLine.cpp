#include <queue>
#include "RuleInOneLine.h"

namespace filter {
    namespace scaffolder {
        void RuleInOneLine::simplifyGraph(filter::Filter *filter) {
            std::vector<int> vect = filter->getVertexList();
            for (int v : vect) {
                std::vector<int> edges = filter->getEdges(v);
                for (int i = 0; i < (int)edges.size(); ++i) {
                    for (int j = i + 1; j < (int)edges.size(); ++j) {
                        int u = filter->getEdgeTo(edges[i]);
                        int w = filter->getEdgeTo(edges[j]);

                        if (havePath(filter, u, w)) {
                            std::stringstream ss;
                            ss << edges[j];
                            filter->processQuery(Query(Query::SET_IGNORE_EDGE, ss.str()));
                        } else if (havePath(filter, w, u)) {
                            std::stringstream ss;
                            ss << edges[i];
                            filter->processQuery(Query(Query::SET_IGNORE_EDGE, ss.str()));
                        }
                    }
                }
            }
        }

        bool RuleInOneLine::havePath(Filter *filter, int u, int w) {
            std::vector<int> dist(filter->getVertexCount(), -1);
            std::queue<int> que;
            que.push(u);
            dist[u] = 0;
            while (que.size() != 0) {
                int v = que.front();
                que.pop();

                if (dist[v] == 20) return false;
                std::vector<int> edges = filter->getEdges(v);
                for (int edge : edges) {
                    int y = filter->getEdgeTo(edge);
                    if (dist[y] == -1) {
                        dist[y] = dist[v] + 1;
                        que.push(y);
                        if (y == w) return true;
                    }
                }

            }
            return false;
        }
    }
}
