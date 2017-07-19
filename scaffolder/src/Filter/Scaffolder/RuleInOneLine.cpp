#include <queue>
#include "RuleInOneLine.h"

namespace filter {
    namespace scaffolder {
        void RuleInOneLine::simplifyGraph(ContigGraph *graph) {
            INFO("start simplify graph");
            std::vector<int> vect = graph->getVertexList();
            for (int v : vect) {
                std::vector<int> edges = graph->getEdges(v);
                TRACE("vertex num = " << v << " edegs cnt=" << edges.size());
                if (edges.size() != 2) continue;
                for (int i = 0; i < (int)edges.size(); ++i) {
                    for (int j = i + 1; j < (int)edges.size(); ++j) {
                        int u = graph->getEdgeTo(edges[i]);
                        int w = graph->getEdgeTo(edges[j]);

                        if (havePath(graph, u, w)) {
                            graph->delEdge(edges[j]);
                        } else if (havePath(graph, w, u)) {
                            graph->delEdge(edges[i]);
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
