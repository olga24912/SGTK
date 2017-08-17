#include "RuleCoord.h"


namespace filter {
    namespace scaffolder {
        using namespace contig_graph;
        void RuleCoord::simplifyGraph(ContigGraph *graph) {
            std::vector<int> vert = graph->getVertexList();
            for (int v : vert) {
                directSimpl(graph, v);
            }

            for (int v : vert) {
                revSimpl(graph, v);
            }

            for (int e : edgeForDel) {
                graph->delEdge(e);
            }
        }

        void RuleCoord::directSimpl(ContigGraph *graph, int v) {
            int cntD = cntEdge(graph, graph->getEdges(v));
            int cntR = cntEdge(graph, graph->getEdgesR(v));

            if (cntD > 1 && cntR == 1) {
                int u = graph->getEdgeFrom(graph->getEdgesR(v)[0]);
                if (cntEdge(graph, graph->getEdges(u)) != 1) return;

                int cb = graph->getEdgeCoordB2(graph->getEdgesR(v)[0]);

                std::vector<int> edges = graph->getEdges(v);
                for (int e : edges) {
                    int ce = graph->getEdgeCoordE1(e);
                    if (ce + 500 < cb) {
                        edgeForDel.insert(e);
                    }
                }
            }
        }

        void RuleCoord::revSimpl(ContigGraph *graph, int v) {
            int cntD = cntEdge(graph, graph->getEdgesR(v));
            int cntR = cntEdge(graph, graph->getEdges(v));

            if (cntD > 1 && cntR == 1) {
                int u = graph->getEdgeTo(graph->getEdges(v)[0]);
                if (cntEdge(graph, graph->getEdgesR(u)) != 1) return;

                int cb = graph->getEdgeCoordE1(graph->getEdges(v)[0]);

                std::vector<int> edges = graph->getEdgesR(v);
                for (int e : edges) {
                    int ce = graph->getEdgeCoordB2(e);
                    if (cb + 500 < ce) {
                        edgeForDel.insert(e);
                    }
                }
            }
        }

        int RuleCoord::cntEdge(ContigGraph *graph, std::vector<int> edges) {
            int cnt = 0;
            for (int i = 0; i < (int)edges.size(); ++i) {
                int was = 0;
                for (int j = i + 1; j < (int)edges.size(); ++j) {
                    if (graph->getEdgeFrom(edges[i]) == graph->getEdgeFrom(edges[j]) &&
                            graph->getEdgeTo(edges[i]) == graph->getEdgeTo(edges[j]) &&
                            sameCoord1(graph, edges[i], edges[j]) &&
                            sameCoord2(graph, edges[i], edges[j])) {
                        was = 1;
                        break;
                    }
                }
                if (was == 0) ++cnt;
            }

            return cnt;
        }
    }
}
