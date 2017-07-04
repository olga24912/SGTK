#include "ScaffoldStrategyUniqueConnection.h"

namespace filter {
    namespace scaffolder {
        void
        ScaffoldStrategyUniqueConnection::addConnection(Scaffolds *scaffolds, Filter *graph, std::vector<int> minW) {
            DEBUG("addConnection");
            topSort(graph);
            findCycle(graph);
            for (int i = 0; i < (int) topsort.size(); ++i) {
                int v = topsort[i];
                std::vector<int> edges = (graph->getEdges(v));
                if (edges.size() == 0) {
                    continue;
                }
                int u = graph->getEdgeTo(edges[0]);
                if (isUniquePair(v, u, graph) && color[v] != color[u] &&
                    graph->getEdgeWeight(edges[0]) >= minW[graph->getLibType(graph->getEdgeLib(edges[0]))]) {
                    scaffolds->addConnection(v, u);
                }
            }
        }

        bool ScaffoldStrategyUniqueConnection::isUniquePair(int v1, int v2, Filter *graph) {
            TRACE("isUniquePair v1=" << v1 << " v2=" << v2);
            std::vector<int> edges = graph->getEdges(v1);
            for (int e : edges) {
                int u = graph->getEdgeTo(e);
                if (u != v2) {
                    return false;
                }
            }

            edges = graph->getEdgesR(v2);
            for (int e : edges) {
                int u = graph->getEdgeFrom(e);
                if (u != v1) {
                    return false;
                }
            }

            return true;
        }
    }
}
