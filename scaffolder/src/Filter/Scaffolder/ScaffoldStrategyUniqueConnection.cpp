#include "ScaffoldStrategyUniqueConnection.h"

void ScaffoldStrategyUniqueConnection::addConnection(Scaffolds *scaffolds, Filter *graph) {
    topSort(graph);
    findCycle(graph);
    for (int i = 0; i < (int)topsort.size(); ++i) {
        int v = topsort[i];
        std::vector<int> edges = (graph->getEdges(v));
        if (edges.size() == 0) {
            continue;
        }
        int u = graph->getEdgeTo(edges[0]);
        if (isUniquePair(v, u, graph) && color[v] != color[u]) {
            scaffolds->addConnection(v, u);
        }
    }
}

bool ScaffoldStrategyUniqueConnection::isUniquePair(int v1, int v2, Filter *graph) {
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
