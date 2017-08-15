#include "RuleDelSmallEdges.h"

using namespace filter::contig_graph;

namespace filter {
    namespace scaffolder {
        void RuleDelSmallEdges::simplifyGraph(ContigGraph *graph) {
            std::vector<int> vertList = graph->getVertexList();
            for (int v : vertList) {
                delEdges(graph, graph->getEdges(v));
                delEdges(graph, graph->getEdgesR(v));
            }

            for (int e : delEdge) {
                graph->delEdge(e);
            }
        }

        void RuleDelSmallEdges::delEdges(ContigGraph *graph, std::vector<int> edges) {
            std::vector<int> w;
            for (int e : edges) {
                if (isAlone(graph, e)) {
                    w.push_back(graph->getEdgeWeight(e));
                }
            }

            std::sort(w.begin(), w.end());
            int curW = startDel;
            for (int i = 0; i < (int)w.size(); ++i) {
                if (w[i] <= curW * difDel) {
                    curW = std::max(curW, w[i]);
                }
            }

            for (int e : edges) {
                if (isAlone(graph, e) && graph->getEdgeWeight(e) <= curW) {
                    delEdge.insert(e);
                }
            }
        }

        bool RuleDelSmallEdges::isAlone(ContigGraph *graph, int e) {
            int u = graph->getEdgeFrom(e), v = graph->getEdgeTo(e);

            std::vector<int> eds = graph->getEdges(u);

            for (int ed : eds) {
                if (ed != e) {
                    if (graph->getEdgeTo(ed) == v && sameCoord1(graph, ed, e) && sameCoord2(graph, ed, e)) {
                        return false;
                    }
                }
            }

            return true;
        }

    }
}
