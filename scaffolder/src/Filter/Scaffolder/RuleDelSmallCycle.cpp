#include <set>
#include "RuleDelSmallCycle.h"

namespace filter {
    namespace scaffolder {
        void filter::scaffolder::RuleDelSmallCycle::simplifyGraph(ContigGraph *graph) {
            INFO("start del small cycle");
            std::vector<int> vert = graph->getVertexList();
            for (int v : vert) {
                std::vector<int> edges = graph->getEdges(v);
                std::vector<int> egdesR = graph->getEdgesR(v);

                std::set<int> edgeForDel;

                for (int e : edges) {
                    for (int er : egdesR) {
                        int u = graph->getEdgeTo(e);
                        int w = graph->getEdgeFrom(er);
                        if (u == w) {
                            edgeForDel.insert(e);
                            edgeForDel.insert(er);
                        }
                    }
                }

                for (int e : edgeForDel) {
                    graph->delEdge(e);
                }
            }
            INFO("finish del small cycle");
        }
    }
}
