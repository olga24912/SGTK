
#include <set>
#include "ValidatorNotSimplePath.h"

namespace filter {
    namespace writers {
        bool ValidatorNotSimplePath::isGoodVertexSet(std::vector<int> vert,
                                                                      filter::contig_graph::ContigGraph *graph) {
            DEBUG("start validation");
            for (int v : vert) {
                std::vector<int> edges = graph->getEdges(v);
                std::set<int> nb;

                for (int e : edges) {
                    int u = graph->getEdgeTo(e);
                    nb.insert(u);
                }

                if (nb.size() > 1) {
                    return true;
                }
            }

            return false;
        }
    }
}