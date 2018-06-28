#include <set>
#include "ValidatorFork.h"

namespace filter {
    namespace writers {
        bool ValidatorFork::isGoodVertexSet(std::vector<int> vert, ContigGraph *graph) {
            DEBUG("start validation");
            for (int v : vert) {
                std::vector<int> edges = graph->getEdges(v);
                std::set<int> nb;

                for (int e : edges) {
                    int u = graph->getEdgeTo(e);
                    if (graph->getEdgeLib(e) != lib) {
                        continue;
                    }
                    nb.insert(e);
                }

                if (nb.size() > 1) {
                    return true;
                }
            }

            return false;
        }
    }
}