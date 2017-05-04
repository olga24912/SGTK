#include "ScaffoldStrategyOneLine.h"

void ScaffoldStrategyOneLine::addConnection(Scaffolds *scaffolds, Filter *graph) {
    topSort(graph);
    findCycle(graph);
    addFirstConnection(scaffolds, graph);
    delEdgeFromDifPath(scaffolds, graph);
}

void ScaffoldStrategyOneLine::addFirstConnection(Scaffolds *scaffolds, Filter *graph) {
    for (int i = 0; i < (int)topsort.size(); ++i) {
        int v = topsort[i];

        std::vector<int> edges = graph->getEdges(v);
        if (edges.size() == 0) continue;
        int minu = graph->getEdgeTo(edges[0]);

        for (int e : edges) {
            int u = graph->getEdgeTo(e);

            if (topSortPos[u] < topSortPos[minu] && color[v] != color[u]) {
                minu = u;
            }
        }
        if (color[v] != color[minu] && scaffolds->lineId(v) != scaffolds->lineId(minu) &&
                scaffolds->isLast(v) && scaffolds->isFirst(minu)) {
            scaffolds->addConnection(v, minu);
        }
    }
}

void ScaffoldStrategyOneLine::delEdgeFromDifPath(Scaffolds *scaffolds, Filter *graph) {
    for (int i = (int)topsort.size() - 1; i >= 0; --i) {
        int v = topsort[i];
        std::vector<int> edges = graph->getEdges(v);
        if (edges.size() == 0) continue;

        for (int e : edges) {
            int u = graph->getEdgeTo(e);
            assert(graph->getEdgeFrom(e) == v);

            if (scaffolds->lineId(v) != scaffolds->lineId(u)) {
                scaffolds->brokeConnection(v);
                scaffolds->brokeConnectionTo(u);
            }

        }
    }
}