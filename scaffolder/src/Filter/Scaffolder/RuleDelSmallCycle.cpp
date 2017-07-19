#include "RuleDelSmallCycle.h"

namespace filter {
    namespace scaffolder {
        void filter::scaffolder::RuleDelSmallCycle::simplifyGraph(ContigGraph *graph) {
            std::vector<int> vert = graph->getVertexList();
            for (int v : vert) {
                std::vector<int> edges = graph->getEdges(v);
                std::vector<int> egdesR = graph->getEdgesR(v);

                for (int e : edges) {
                    for (int er : egdesR) {
                        int u = graph->getEdgeTo(e);
                        int w = graph->getEdgeFrom(er);
                        if (u == w) {
                            if (MAX_DIF * graph->getEdgeWeight(e) <= graph->getEdgeWeight(er)) {
                                graph->delEdge(e);
                            } else if (MAX_DIF * graph->getEdgeWeight(er) <= graph->getEdgeWeight(e)) {
                                graph->delEdge(er);
                            } else {
                                graph->delEdge(e);
                                graph->delEdge(er);
                            }
                        }
                    }
                }
            }
        }
    }
}
