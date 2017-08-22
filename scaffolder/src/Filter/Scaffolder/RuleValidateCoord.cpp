#include "RuleValidateCoord.h"

namespace filter {
    namespace scaffolder {
        using namespace contig_graph;
        void RuleValidateCoord::simplifyGraph(ContigGraph *filter) {
            std::vector<int> vert = filter->getVertexList();

            for (int v : vert) {
                std::vector<int> edges = filter->getEdges(v);
                for (int e : edges) {
                    if (filter->getLibType(filter->getEdgeLib(e)) == ContigGraph::Lib::RNA_SPLIT_50) {
                        validateEdge(filter, e);
                    }
                }

            }

        }

        void RuleValidateCoord::validateEdge(ContigGraph *graph, int e) {
            int v = graph->getEdgeFrom(e);
            int u = graph->getEdgeTo(e);

            std::vector<int> edges = graph->getEdges(v);
            int needDel = 0;

            for (int e1 : edges) {
                if (graph->getEdgeTo(e1) == u &&
                        sameCoord1(graph, e1, e) &&
                        sameCoord2(graph, e1, e)) {
                    if (graph->getEdgeCoordE1(e1) - 15 > graph->getEdgeCoordE1(e) ||
                            graph->getEdgeCoordB2(e1) + 15 < graph->getEdgeCoordB2(e)) {
                        needDel = 1;
                    }
                }
            }

            if (needDel == 1) {
                for (int e1 : edges) {
                    if (graph->getEdgeTo(e1) == u &&
                        sameCoord1(graph, e1, e) &&
                        sameCoord2(graph, e1, e) && e1 != e) {
                        graph->delEdge(e1);
                    }
                }
                graph->delEdge(e);
            }
        }
    }
}
