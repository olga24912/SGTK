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

        void RuleDelSmallEdges::delEdges(ContigGraph *graph, const std::vector<int>& edges) {
            //if (edges.size() > 20) {
            //    INFO("edge cnt = " << edges.size());
            //}
            std::vector<int> w;
            std::vector<bool> notAlone(edges.size(), false);

            for (int i = 0; i < edges.size(); ++i) {
                int e = edges[i];
                for (int j = i + 1; notAlone[i] == false && j < edges.size() && graph->getEdgeFrom(edges[j]) ==
                        graph->getEdgeFrom(edges[i]) &&
                        graph->getEdgeTo(edges[i]) ==
                        graph->getEdgeTo(edges[j]); ++j) {
                    if (sameCoord1(graph, edges[i], edges[j]) && sameCoord2(graph, edges[i], edges[j])) {
                        notAlone[i] = true;
                        notAlone[j] = true;
                    }
                }

                if (notAlone[i] == false) {
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

            for (int i = 0; i < (int)edges.size(); ++i) {
                int e = edges[i];
                if (notAlone[i] == false && graph->getEdgeWeight(e) <= curW) {
                    delEdge.insert(e);
                }
            }
        }
    }
}
