#include "ScaffoldStrategyUniqueConnection.h"

namespace filter {
    namespace scaffolder {
        void
        ScaffoldStrategyUniqueConnection::addConnection(Scaffolds *scaffolds, ContigGraph *graph, std::vector<int> minW) {
            DEBUG("addConnection");
            topSort(graph);
            findCycle(graph);
            DEBUG("color[9385] = " << color[9385] << " color[10182] = " << color[10182]);
            DEBUG("is unique pair 9385 10182 " << isUniquePair(9385, 10182, graph));
            for (int i = 0; i < (int) topsort.size(); ++i) {
                int v = topsort[i];
                std::vector<int> edges = (graph->getEdges(v));
                if (edges.size() == 0) {
                    continue;
                }
                int u = graph->getEdgeTo(edges[0]);
                if (v == 9385) DEBUG("9385 pair " << u);
                if (v == 10183) DEBUG("10183 pair " << u);
                if (v == 9385 && u == 10182) {
                    DEBUG("9385-10182 w=" << graph->getEdgeWeight(edges[0]));
                    DEBUG("color[9385] = " << color[v] << " color[10182] = " << color[u]);
                    DEBUG("is unique pair 9385 10182 " << isUniquePair(v, u, graph));
                    DEBUG(isUniquePair(v, u, graph));
                    DEBUG((color[v] != color[u]));
                    DEBUG((graph->getEdgeWeight(edges[0]) >= 3));

                }
                if (isUniquePair(v, u, graph) && color[v] != color[u] && graph->getEdgeWeight(edges[0]) >= 3) {
                    if ((v == 9385 && u == 10182) || (v == 10183 && u == 9384)) DEBUG("add connection 9385->10182");
                    scaffolds->addConnection(v, u);
                    scaffolds->addConnection(u^1, v^1);
                }
            }
        }

        bool ScaffoldStrategyUniqueConnection::isUniquePair(int v1, int v2, ContigGraph *graph) {
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
