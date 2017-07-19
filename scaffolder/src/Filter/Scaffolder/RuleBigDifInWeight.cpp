#include "RuleBigDifInWeight.h"

namespace filter {
    namespace scaffolder {
        void RuleBigDifInWeight::simplifyGraph(ContigGraph *graph) {
            INFO("start simplify graph");
            std::vector<int> vert = graph->getVertexList();
            for (int v : vert) {
                delSmallEdges(graph, graph->getEdges(v));
                delSmallEdges(graph, graph->getEdgesR(v));
            }
            INFO("finish simplify graph");
        }

        void RuleBigDifInWeight::delSmallEdges(ContigGraph *graph, const std::vector<int> &edges) const {
            int maxW = 0;
            for (int e : edges) {
                if (graph->getEdgeWeight(e) >= maxW) {
                    maxW = graph->getEdgeWeight(e);
                }
            }

            for (int e : edges) {
                if (graph->getEdgeWeight(e) * maxDif <= maxW) {
                    graph->delEdge(e);
                }
            }
        }
    }
}
