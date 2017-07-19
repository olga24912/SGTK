#include "ValidatorError.h"

namespace filter {
    namespace writers {
        bool ValidatorError::isGoodVertexSet(std::vector<int> vert, ContigGraph *graph) {
            aligInfo.parseCoordFile(graph, coordFileName);

            int haveError = 0;
            for (int v : vert) {
                std::vector<int> edges = graph->getEdges(v);

                for (int e : edges) {
                    if (graph->getEdgeLib(e) != libError) continue;
                    int u = graph->getEdgeTo(e);
                    int was = 0;
                    for (int w : vert) {
                        if (w == u) {
                            was = 1;
                        }
                    }
                    if (was == 0) continue;

                    if (aligInfo.isCorrectEdge(graph, e) != InfoAboutContigsAlig::OK) {
                        haveError = 1;
                    }
                }
            }

            return haveError;
        }
    }
}