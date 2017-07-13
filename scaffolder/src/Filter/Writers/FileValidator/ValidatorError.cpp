#include "ValidatorError.h"

namespace filter {
    namespace writers {
        bool ValidatorError::isGoodVertexSet(std::vector<int> vert, filter::Filter *filter) {
            aligInfo.parseCoordFile(filter, coordFileName);

            int haveError = 0;
            for (int v : vert) {
                std::vector<int> edges = filter->getEdges(v);

                for (int e : edges) {
                    if (filter->getEdgeLib(e) != libError) continue;
                    int u = filter->getEdgeTo(e);
                    int was = 0;
                    for (int w : vert) {
                        if (w == u) {
                            was = 1;
                        }
                    }
                    if (was == 0) continue;

                    if (aligInfo.isCorrectEdge(filter, e) != InfoAboutContigsAlig::OK) {
                        haveError = 1;
                    }
                }
            }

            return haveError;
        }
    }
}