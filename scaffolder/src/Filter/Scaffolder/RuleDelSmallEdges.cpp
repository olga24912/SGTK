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
                w.push_back(graph->getEdgeWeight(e));
            }

            std::sort(w.begin(), w.end());
            int curW = startDel;
            for (int i = 0; i < (int)w.size(); ++i) {
                if (w[i] <= curW * difDel) {
                    curW = std::max(curW, w[i]);
                }
            }

            for (int e : edges) {
                if (graph->getEdgeWeight(e) <= curW) {
                    delEdge.insert(e);
                }
            }
        }
    }
}
