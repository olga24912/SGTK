#include <set>
#include "ValidatorFork.h"

bool ValidatorFork::isGoodVertexSet(std::vector<int> vert, Filter *filter) {
    DEBUG("start validation");
    for (int v : vert) {
        std::vector<int> edges = filter->getEdges(v);
        std::set<int> nb;

        for (int e : edges) {
            int u = filter->getEdgeTo(e);
            if (filter->getEdgeLib(e) != lib) {
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
