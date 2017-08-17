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

            for (int e : edgeForDel) {
                graph->delEdge(e);
            }

            INFO("finish simplify graph");
        }

        void RuleBigDifInWeight::delSmallEdges(ContigGraph *graph, const std::vector<int> &edges) {
            int hasGoodConnection = 0;
            for (int e : edges) {
                if (!hasGoodConnection && isGoodConnection(graph, e, edges)) {
                    hasGoodConnection = 1;
                }
            }

            for (int e : edges) {
                if (hasGoodConnection && isBadConnection(graph, e, edges)) {
                    edgeForDel.insert(e);
                }
            }
        }

        bool RuleBigDifInWeight::isGoodConnection(ContigGraph *graph, int e, const std::vector<int>& edges) const {
            int wasPair = 0, wasS50 = 0;

            if (graph->getEdgeWeight(e) >= 5) {
                if (graph->getLibType(graph->getEdgeLib(e)) == ContigGraph::Lib::RNA_PAIR) {
                    wasPair = 1;
                } else if (graph->getLibType(graph->getEdgeLib(e)) == ContigGraph::Lib::RNA_SPLIT_50) {
                    wasS50 = 1;
                }
            }

            for (int e1 : edges) {
                if (graph->getEdgeWeight(e1) < 5) continue;
                if (graph->getEdgeFrom(e1) == graph->getEdgeFrom(e) &&
                        graph->getEdgeTo(e1) == graph->getEdgeTo(e) &&
                        sameCoord1(graph, e1, e) && sameCoord2(graph, e1, e)) {
                    if (graph->getLibType(graph->getEdgeLib(e1)) == ContigGraph::Lib::RNA_PAIR) {
                        wasPair = 1;
                    } else if (graph->getLibType(graph->getEdgeLib(e1)) == ContigGraph::Lib::RNA_SPLIT_50) {
                        wasS50 = 1;
                    }
                }
            }

            return wasPair && wasS50;
        }

        bool RuleBigDifInWeight::isBadConnection(ContigGraph *graph, int e, const std::vector<int>& edges) const {
            int wasPair = 0, wasS50 = 0;
            int w = 0;

            for (int e1 : edges) {
                if (graph->getEdgeFrom(e1) == graph->getEdgeFrom(e) &&
                    graph->getEdgeTo(e1) == graph->getEdgeTo(e) &&
                    sameCoord1(graph, e1, e) && sameCoord2(graph, e1, e)) {
                    if (graph->getLibType(graph->getEdgeLib(e1)) == ContigGraph::Lib::RNA_PAIR) {
                        if (wasS50) return false;
                        wasPair = 1;
                    } else if (graph->getLibType(graph->getEdgeLib(e1)) == ContigGraph::Lib::RNA_SPLIT_50) {
                        if (wasPair) return false;
                        wasS50 = 1;
                    }
                    w += graph->getEdgeWeight(e1);
                }
            }
            return w < 10;
        }
    }
}
