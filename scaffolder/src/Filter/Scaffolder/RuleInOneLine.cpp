#include <queue>
#include "RuleInOneLine.h"

namespace filter {
    namespace scaffolder {
        void RuleInOneLine::simplifyGraph(filter::ContigGraph *filter) {
            INFO("start simplify graph");
            std::vector<int> vect = filter->getVertexList();
            for (int v : vect) {
                std::vector<int> edges = filter->getEdges(v);
                TRACE("vertex num = " << v << " edegs cnt=" << edges.size());
                if (edges.size() != 2) continue;
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
            INFO("finish simplify graph");
        }

        bool RuleInOneLine::havePath(ContigGraph *filter, int u, int w) {
            std::unordered_map<int, int> dist;
            std::queue<int> que;
            que.push(u);
            dist[u] = 0;
            while (que.size() != 0) {
                int v = que.front();
                que.pop();

                if (dist[v] == 10) return false;
                std::vector<int> edges = filter->getEdges(v);
                for (int edge : edges) {
                    int y = filter->getEdgeTo(edge);
                    if (dist.count(y) == 0) {
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
