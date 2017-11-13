#include "RuleBigDeg.h"

namespace filter {
    namespace scaffolder {
        void RuleBigDeg::simplifyGraph(ContigGraph *graph) {
            std::vector<int> vert = graph->getVertexList();
            std::vector<int> wignore(graph->getVertexCount(), 0);
            std::vector<int> wignoreR(graph->getVertexCount(), 0);

            for (int v : vert) {
                std::vector<int> edges = graph->getEdges(v);
                if (edges.size() >= BIG_DEG) {
                    wignore[v] = graph->getEdgeWeight(edges[0]);
                    for (int e : edges) {
                        wignore[v] = std::max(wignore[v], graph->getEdgeWeight(e));
                    }
                }

                std::vector<int> edgesR = graph->getEdgesR(v);
                if (edgesR.size() >= BIG_DEG) {
                    wignoreR[v] = graph->getEdgeWeight(edgesR[0]);
                    for (int e : edgesR) {
                        wignoreR[v] = std::max(wignoreR[v], graph->getEdgeWeight(e));
                    }
                }
            }

            ignoreEdges(graph, wignore);
            ignoreEdgesR(graph, wignoreR);
        }

        void RuleBigDeg::ignoreEdges(ContigGraph *graph, std::vector<int> wig) {
            std::vector<int> vect = graph->getVertexList();
            for (int v : vect) {
                if (wig[v] != 0 && wig[v] <= MAX_WEIGHT) {
                    std::vector<int> edges = graph->getEdges(v);
                    for (int e : edges) {
                        graph->delEdge(e);
                    }
                }
            }
        }

        void RuleBigDeg::ignoreEdgesR(ContigGraph *graph, std::vector<int> wig) {
            std::vector<int> vect = graph->getVertexList();
            for (int v : vect) {
                if (wig[v] != 0 && wig[v] <= MAX_WEIGHT) {
                    std::vector<int> edges = graph->getEdgesR(v);
                    for (int e : edges) {
                        graph->delEdge(e);
                    }
                }
            }
        }
    }
}
