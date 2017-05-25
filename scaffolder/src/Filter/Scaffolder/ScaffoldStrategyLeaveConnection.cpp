#include "ScaffoldStrategyLeaveConnection.h"

void ScaffoldStrategyLeaveConnection::addConnection(Scaffolds *scaffolds, Filter *graph, std::vector<int> minW) {
    topSort(graph);
    findCycle(graph);
    addConnectionForLeaves(scaffolds, graph, 0, minW);
    addConnectionForLeaves(scaffolds, graph, 1, minW);
}

void ScaffoldStrategyLeaveConnection::addConnectionForLeaves(Scaffolds *scaffolds, Filter *graph, int dir, std::vector<int> minW) {
    int n = graph->getVertexCount();

    for (int i = 0; i < n; ++i) {
        std::vector<int> edges = dir == 0 ? graph->getEdges(i) : graph->getEdgesR(i);


        int u = -1, v = -1;
        std::vector<int> wu(ContigGraph::Lib::typeCnt, 0), wv(ContigGraph::Lib::typeCnt, 0);
        for (int e : edges) {
            int to = dir == 0 ? graph->getEdgeTo(e) : graph->getEdgeFrom(e);
            if (u == -1 || u == to) {
                u = to;
                wu[graph->getLibType(graph->getEdgeLib(e))] = std::max(wu[graph->getLibType(graph->getEdgeLib(e))], graph->getEdgeWeight(e));
            } else if (v == -1 || v == to) {
                v = to;
                wv[graph->getLibType(graph->getEdgeLib(e))] = std::max(wv[graph->getLibType(graph->getEdgeLib(e))], graph->getEdgeWeight(e));
            } else {
                u = -2;
                break;
            }
        }

        if (u < 0 || v < 0) {
            continue;
        }

        if (color[i] == color[u] || color[i] == color[v] || color[u] == color[v]) {
            continue;
        }

        if (deg(v, graph, dir) != 0) {
            std::swap(v, u);
            std::swap(wv, wu);
        }
        if (deg(i, graph, !dir) <= 1 &&
            deg(v, graph, !dir) == 1 && deg(u, graph, !dir) == 1 &&
            deg(v, graph, dir) == 0 && deg(u, graph, dir) <= 1) {
            int flagBigW = 0;
            for (int j = 0; j < ContigGraph::Lib::typeCnt; ++j) {
                if (wu[j] >= minW[j]) {
                    flagBigW = 1;
                }
            }
            if (dir == 0 && flagBigW) {
                scaffolds->addConnection(i, u);
            } else if (flagBigW) {
                scaffolds->addConnection(u, i);
            }
        }

    }
}
